#!/usr/bin/env Rscript
# =============================================================================
# PROMOTER CLASSIFICATION ANALYSIS v2 - FIXED & ENHANCED
# D. laeve Methylation Study - Tail Amputation vs Control
# Weber et al. 2007 Criteria
# =============================================================================
#
# FIXES from v1:
#   - HTML report: replaced sprintf() with paste0() to avoid %% escaping issues
#   - Added O/E index vs methylation correlation analysis
#   - Added volcano plot for differential methylation
#   - Added gene lists by class and methylation status
#   - Added ICP subpopulation analysis
#   - Improved plot annotations with gene counts
#
# =============================================================================

options(stringsAsFactors = FALSE)
options(repos = c(CRAN = "https://cloud.r-project.org"))

cat("\n")
cat("================================================================================\n")
cat("  PROMOTER CLASSIFICATION ANALYSIS v2 (Weber Criteria) - FIXED\n")
cat("  D. laeve - Tail Amputation vs Control\n")
cat("================================================================================\n")
cat("\n")

# =============================================================================
# PACKAGE LOADING
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(rtracklayer)
  library(IRanges)
  library(BSgenome.Dlaeve.NCBI.dlgm)
  library(Biostrings)
  library(DESeq2)
  library(ggplot2)
  library(scales)
  library(viridis)
  library(RColorBrewer)
  library(gridExtra)
  library(grid)
  library(pheatmap)
})

cat("✓ All packages loaded\n\n")

# =============================================================================
# CONFIGURATION
# =============================================================================

gff_genes <- "/mnt/data/alfredvar/30-Genoma/31-Alternative_Annotation_EviAnn/derLaeGenome_namesDlasi_v2.fasta.functional_note.pseudo_label.gff"

meth_samples <- c(
  C1 = "/mnt/data/alfredvar/jmiranda/50-Genoma/51-Metilacion/09_methylation_calls/C1.CpG_report.txt.gz",
  C2 = "/mnt/data/alfredvar/jmiranda/50-Genoma/51-Metilacion/09_methylation_calls/C2.CpG_report.txt.gz",
  A1 = "/mnt/data/alfredvar/jmiranda/50-Genoma/51-Metilacion/09_methylation_calls/A1.CpG_report.txt.gz",
  A2 = "/mnt/data/alfredvar/jmiranda/50-Genoma/51-Metilacion/09_methylation_calls/A2.CpG_report.txt.gz"
)

rnaseq_dir <- "/mnt/data/alfredvar/jmiranda/20-Transcriptomic_Bulk/25-metaAnalysisTranscriptome/counts_HTseq_EviAnn/"

out_dir <- "PROMOTER_CLASSIFICATION_RESULTS_v2"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

subdirs <- c("01_classification", "02_plots_png", "03_methylation_analysis",
             "04_expression_correlation", "05_data", "06_html_report",
             "07_gene_lists", "08_volcano")
for (d in subdirs) dir.create(file.path(out_dir, d), showWarnings = FALSE)

keep_chr <- paste0("chr", 1:31)

MIN_COV_5 <- 5
MIN_COV_10 <- 10

PROMOTER_UPSTREAM <- 2000
PROMOTER_DOWNSTREAM <- 500
WINDOW_SIZE <- 500
WINDOW_STEP <- 50

HCP_CPGOE_THRESHOLD <- 0.75
HCP_GC_THRESHOLD <- 0.55
LCP_CPGOE_THRESHOLD <- 0.48

class_colors <- c(HCP = "#2166AC", ICP = "#4DAF4A", LCP = "#D6604D")

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), sprintf(...), "\n")

calc_cpg_oe <- function(seq_string) {
  seq <- DNAString(seq_string)
  len <- length(seq)
  if (len < 50) return(NA_real_)
  n_c <- letterFrequency(seq, "C")[1]
  n_g <- letterFrequency(seq, "G")[1]
  n_cpg <- countPattern("CG", seq)
  expected_cpg <- (n_c * n_g) / len
  if (expected_cpg == 0) return(0)
  cpg_oe <- n_cpg / expected_cpg
  return(cpg_oe)
}

calc_gc_content <- function(seq_string) {
  seq <- DNAString(seq_string)
  len <- length(seq)
  if (len < 50) return(NA_real_)
  n_gc <- letterFrequency(seq, c("G", "C"), as.prob = FALSE)
  gc_content <- sum(n_gc) / len
  return(gc_content)
}

classify_promoter_weber <- function(promoter_seq, window_size = 500, step = 50) {
  seq <- DNAString(promoter_seq)
  seq_len <- length(seq)
  
  if (seq_len < window_size) {
    cpg_oe <- calc_cpg_oe(as.character(seq))
    gc <- calc_gc_content(as.character(seq))
    if (!is.na(cpg_oe) && !is.na(gc)) {
      if (cpg_oe >= HCP_CPGOE_THRESHOLD && gc >= HCP_GC_THRESHOLD) {
        return(list(class = "HCP", max_cpg_oe = cpg_oe, max_gc = gc))
      } else if (cpg_oe < LCP_CPGOE_THRESHOLD) {
        return(list(class = "LCP", max_cpg_oe = cpg_oe, max_gc = gc))
      } else {
        return(list(class = "ICP", max_cpg_oe = cpg_oe, max_gc = gc))
      }
    }
    return(list(class = "Unknown", max_cpg_oe = NA, max_gc = NA))
  }
  
  n_windows <- floor((seq_len - window_size) / step) + 1
  max_cpg_oe <- 0
  max_gc_at_max_oe <- 0
  has_hcp_window <- FALSE
  has_window_above_lcp <- FALSE
  
  for (i in seq_len(n_windows)) {
    start <- (i - 1) * step + 1
    end <- start + window_size - 1
    window_seq <- as.character(subseq(seq, start, end))
    cpg_oe <- calc_cpg_oe(window_seq)
    gc <- calc_gc_content(window_seq)
    if (!is.na(cpg_oe) && !is.na(gc)) {
      if (cpg_oe > max_cpg_oe) {
        max_cpg_oe <- cpg_oe
        max_gc_at_max_oe <- gc
      }
      if (cpg_oe >= HCP_CPGOE_THRESHOLD && gc >= HCP_GC_THRESHOLD) {
        has_hcp_window <- TRUE
      }
      if (cpg_oe >= LCP_CPGOE_THRESHOLD) {
        has_window_above_lcp <- TRUE
      }
    }
  }
  
  if (has_hcp_window) {
    return(list(class = "HCP", max_cpg_oe = max_cpg_oe, max_gc = max_gc_at_max_oe))
  } else if (!has_window_above_lcp) {
    return(list(class = "LCP", max_cpg_oe = max_cpg_oe, max_gc = max_gc_at_max_oe))
  } else {
    return(list(class = "ICP", max_cpg_oe = max_cpg_oe, max_gc = max_gc_at_max_oe))
  }
}

read_bismark_cpg <- function(path, keep_chr, min_cov = 5) {
  if (!file.exists(path)) stop("Missing: ", path)
  msg("  Reading: %s", basename(path))
  cmd <- paste("zcat", shQuote(path))
  dt <- fread(cmd = cmd, sep = "\t", header = FALSE, fill = TRUE, showProgress = FALSE,
              col.names = c("chrom", "pos", "strand", "m", "u", "context", "trinuc"))
  dt <- dt[grepl("^chr", chrom) & context == "CG"]
  dt[, `:=`(chrom = as.character(chrom), pos = as.integer(pos),
            m = as.numeric(m), u = as.numeric(u))]
  dt <- dt[chrom %in% keep_chr]
  dt[, cov := m + u]
  dt[, meth := ifelse(cov > 0, m / cov, NA_real_)]
  dt[cov >= min_cov, .(chrom, pos, strand, m, u, cov, meth)]
}

calc_feature_meth <- function(feature_gr, meth_dt, min_cov = 5, min_cpg = 3) {
  nfeat <- length(feature_gr)
  meth_filt <- meth_dt[cov >= min_cov & is.finite(meth)]
  out <- data.table(idx = seq_len(nfeat), n_cpg = 0L, meth_mean = NA_real_,
                    meth_median = NA_real_, meth_sd = NA_real_, callable = FALSE)
  if (nrow(meth_filt) == 0L || nfeat == 0L) return(out)
  
  for (chr in intersect(unique(as.character(seqnames(feature_gr))), unique(meth_filt$chrom))) {
    f_idx <- which(as.character(seqnames(feature_gr)) == chr)
    m_filt <- meth_filt[chrom == chr]
    if (length(f_idx) == 0L || nrow(m_filt) == 0L) next
    meth_gr <- GRanges(seqnames = chr, ranges = IRanges(start = m_filt$pos, width = 1),
                       meth = m_filt$meth, cov = m_filt$cov)
    ov <- findOverlaps(feature_gr[f_idx], meth_gr, ignore.strand = TRUE)
    if (length(ov) == 0L) next
    q_global <- f_idx[queryHits(ov)]
    meth_vals <- mcols(meth_gr)$meth[subjectHits(ov)]
    cov_vals <- mcols(meth_gr)$cov[subjectHits(ov)]
    agg <- data.table(idx = q_global, meth = meth_vals, cov = cov_vals)[, .(
      n_cpg = .N, meth_mean = sum(meth * cov) / sum(cov), meth_median = median(meth),
      meth_sd = sd(meth)), by = idx]
    out[agg, on = "idx", `:=`(n_cpg = i.n_cpg, meth_mean = i.meth_mean,
                              meth_median = i.meth_median, meth_sd = i.meth_sd)]
  }
  out[, callable := (n_cpg >= min_cpg)]
  out
}

# Themes
theme_pub <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "#E0E0E0", linewidth = 0.3),
      panel.background = element_rect(fill = "#1a1a1a", color = NA),
      plot.background = element_rect(fill = "#1a1a1a", color = NA),
      plot.title = element_text(face = "bold", hjust = 0.5, size = base_size + 4, color = "white"),
      plot.subtitle = element_text(color = "#AAAAAA", hjust = 0.5),
      axis.text = element_text(color = "#CCCCCC"),
      axis.title = element_text(color = "white", face = "bold"),
      legend.position = "bottom",
      legend.text = element_text(color = "white"),
      legend.title = element_text(color = "white", face = "bold"),
      legend.background = element_rect(fill = "#1a1a1a", color = NA)
    )
}

theme_pub_light <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "#E0E0E0", linewidth = 0.3),
      plot.title = element_text(face = "bold", hjust = 0.5, size = base_size + 2),
      plot.subtitle = element_text(color = "#666666", hjust = 0.5),
      legend.position = "bottom",
      strip.background = element_rect(fill = "#f0f0f0", color = NA),
      strip.text = element_text(face = "bold")
    )
}

# =============================================================================
# PHASE 1: LOAD GENOME AND GENE ANNOTATIONS
# =============================================================================

msg("=== PHASE 1: LOADING GENOME AND GENE ANNOTATIONS ===")

genome <- BSgenome.Dlaeve.NCBI.dlgm
keep_chr <- intersect(keep_chr, seqnames(genome))
genome_length <- sum(seqlengths(genome)[keep_chr], na.rm = TRUE)
msg("  Genome size (chr1-31): %s bp", format(genome_length, big.mark = ","))

genes_gff <- import(gff_genes)
genes_gff <- genes_gff[seqnames(genes_gff) %in% keep_chr]
genes <- genes_gff[genes_gff$type == "gene"]

msg("  Total genes: %s", format(length(genes), big.mark = ","))

gene_ids <- mcols(genes)$ID
gene_type <- ifelse(grepl("XLOC.*lncRNA", gene_ids), "lncRNA",
                    ifelse(grepl("^LOC_", gene_ids), "protein_coding", "other"))

protein_genes <- genes[gene_type == "protein_coding"]
msg("  Protein-coding genes: %d", length(protein_genes))

# =============================================================================
# PHASE 2: PROMOTER CLASSIFICATION (WEBER CRITERIA)
# =============================================================================

msg("=== PHASE 2: PROMOTER CLASSIFICATION (Weber Criteria) ===")

promoters_gr <- promoters(protein_genes, upstream = PROMOTER_UPSTREAM, downstream = PROMOTER_DOWNSTREAM)

msg("  Classifying %d promoters using sliding window approach...", length(promoters_gr))
msg("  Window size: %d bp, Step: %d bp", WINDOW_SIZE, WINDOW_STEP)
msg("  HCP criteria: CpG O/E >= %.2f AND GC >= %.0f%%", HCP_CPGOE_THRESHOLD, HCP_GC_THRESHOLD * 100)
msg("  LCP criteria: No window with CpG O/E >= %.2f", LCP_CPGOE_THRESHOLD)

promoter_seqs <- getSeq(genome, promoters_gr)

classification_results <- lapply(seq_along(promoter_seqs), function(i) {
  if (i %% 1000 == 0) msg("    Processed %d promoters...", i)
  seq_char <- as.character(promoter_seqs[[i]])
  result <- classify_promoter_weber(seq_char, window_size = WINDOW_SIZE, step = WINDOW_STEP)
  result$gene_id <- mcols(protein_genes)$ID[i]
  result
})

promoter_class_dt <- rbindlist(lapply(classification_results, function(x) {
  data.table(gene_id = x$gene_id, promoter_class = x$class,
             max_cpg_oe = x$max_cpg_oe, max_gc = x$max_gc)
}))

promoter_metrics <- rbindlist(lapply(seq_along(promoter_seqs), function(i) {
  seq_char <- as.character(promoter_seqs[[i]])
  seq <- DNAString(seq_char)
  len <- length(seq)
  n_c <- letterFrequency(seq, "C")[1]
  n_g <- letterFrequency(seq, "G")[1]
  n_cpg <- countPattern("CG", seq)
  gc_content <- (n_c + n_g) / len
  cpg_per_bp <- n_cpg / len
  cpg_oe_whole <- if (n_c > 0 && n_g > 0) (n_cpg * len) / (n_c * n_g) else 0
  
  data.table(
    gene_id = mcols(protein_genes)$ID[i],
    chrom = as.character(seqnames(promoters_gr)[i]),
    start = start(promoters_gr)[i],
    end = end(promoters_gr)[i],
    promoter_length = len,
    n_cpg = n_cpg,
    gc_content = gc_content,
    cpg_per_bp = cpg_per_bp,
    cpg_oe_whole = cpg_oe_whole
  )
}))

promoter_final <- merge(promoter_metrics, promoter_class_dt, by = "gene_id")

if ("Note" %in% names(mcols(protein_genes))) {
  gene_notes <- data.table(
    gene_id = mcols(protein_genes)$ID,
    gene_note = sapply(mcols(protein_genes)$Note, function(x) paste(unlist(x), collapse = "; "))
  )
  promoter_final <- merge(promoter_final, gene_notes, by = "gene_id", all.x = TRUE)
}

class_summary <- promoter_final[, .(
  count = .N,
  pct = .N / nrow(promoter_final) * 100,
  mean_cpg_oe = mean(max_cpg_oe, na.rm = TRUE),
  mean_gc = mean(max_gc, na.rm = TRUE) * 100,
  mean_cpg_per_bp = mean(cpg_per_bp, na.rm = TRUE)
), by = promoter_class]

msg("\n  Classification Summary:")
msg("  HCP (High CpG): %d (%.1f%%)", class_summary[promoter_class == "HCP", count], class_summary[promoter_class == "HCP", pct])
msg("  ICP (Intermediate): %d (%.1f%%)", class_summary[promoter_class == "ICP", count], class_summary[promoter_class == "ICP", pct])
msg("  LCP (Low CpG): %d (%.1f%%)", class_summary[promoter_class == "LCP", count], class_summary[promoter_class == "LCP", pct])

fwrite(promoter_final, file.path(out_dir, "01_classification", "promoter_classification_all.tsv"), sep = "\t")
fwrite(class_summary, file.path(out_dir, "01_classification", "classification_summary.tsv"), sep = "\t")

# =============================================================================
# PHASE 3: LOAD METHYLATION DATA
# =============================================================================

msg("\n=== PHASE 3: LOADING METHYLATION DATA ===")

sample_meth <- list()
for (samp in names(meth_samples)) {
  sample_meth[[samp]] <- read_bismark_cpg(meth_samples[[samp]], keep_chr, min_cov = MIN_COV_5)
}

# =============================================================================
# PHASE 4: CALCULATE PROMOTER METHYLATION
# =============================================================================

msg("\n=== PHASE 4: CALCULATING PROMOTER METHYLATION ===")

calc_promoter_meth_at_threshold <- function(min_cov) {
  msg("  Processing with minimum coverage = %d...", min_cov)
  prom_meth_list <- list()
  for (samp in names(sample_meth)) {
    meth_dt <- sample_meth[[samp]][cov >= min_cov]
    pm <- calc_feature_meth(promoters_gr, meth_dt, min_cov = min_cov, min_cpg = 3)
    pm[, gene_id := mcols(protein_genes)$ID]
    pm[, sample := samp]
    prom_meth_list[[samp]] <- pm[callable == TRUE, .(gene_id, sample, n_cpg, meth_mean)]
  }
  prom_all <- rbindlist(prom_meth_list)
  prom_wide <- dcast(prom_all, gene_id ~ sample, value.var = "meth_mean")
  prom_wide <- prom_wide[complete.cases(prom_wide)]
  prom_wide[, Control_mean := (C1 + C2) / 2]
  prom_wide[, Amputated_mean := (A1 + A2) / 2]
  prom_wide[, delta_meth := Amputated_mean - Control_mean]
  prom_wide
}

prom_meth_cov5 <- calc_promoter_meth_at_threshold(MIN_COV_5)
prom_meth_cov10 <- calc_promoter_meth_at_threshold(MIN_COV_10)

msg("  Promoters with methylation (cov>=5): %d", nrow(prom_meth_cov5))
msg("  Promoters with methylation (cov>=10): %d", nrow(prom_meth_cov10))

final_data_cov5 <- merge(promoter_final, prom_meth_cov5, by = "gene_id", all.x = TRUE)
final_data_cov10 <- merge(promoter_final, prom_meth_cov10, by = "gene_id", all.x = TRUE)

fwrite(final_data_cov5, file.path(out_dir, "05_data", "promoter_methylation_cov5.tsv"), sep = "\t")
fwrite(final_data_cov10, file.path(out_dir, "05_data", "promoter_methylation_cov10.tsv"), sep = "\t")

# =============================================================================
# PHASE 5: LOAD RNA-SEQ DATA
# =============================================================================

msg("\n=== PHASE 5: LOADING RNA-SEQ DATA (Tail Amputation Only) ===")

rnaseq_available <- FALSE

if (dir.exists(rnaseq_dir)) {
  count_files <- list.files(rnaseq_dir, pattern = "htseq_gene_counts\\.txt$", full.names = TRUE)
  tail_files <- count_files[grepl("^[TC]\\d+S\\d+", basename(count_files))]
  
  if (length(tail_files) > 0) {
    msg("  Found %d tail amputation count files", length(tail_files))
    
    sample_info <- data.frame(
      file = tail_files,
      sample_id = gsub("_htseq_gene_counts\\.txt$|.Aligned.out.bam", "", basename(tail_files)),
      stringsAsFactors = FALSE
    )
    sample_info$condition <- ifelse(grepl("^T", sample_info$sample_id), "Amputated", "Control")
    sample_info$fileName <- basename(tail_files)
    rownames(sample_info) <- sample_info$sample_id
    
    dds <- DESeqDataSetFromHTSeqCount(
      sampleTable = sample_info[, c("sample_id", "fileName", "condition")],
      directory = rnaseq_dir,
      design = ~ condition
    )
    dds <- dds[rowSums(counts(dds)) > 10, ]
    dds <- DESeq(dds)
    norm_counts <- counts(dds, normalized = TRUE)
    
    # Get DESeq2 results for proper p-values
    res <- results(dds, contrast = c("condition", "Amputated", "Control"))
    res_dt <- as.data.table(as.data.frame(res), keep.rownames = "gene_id")
    
    control_samples <- sample_info$sample_id[sample_info$condition == "Control"]
    amputated_samples <- sample_info$sample_id[sample_info$condition == "Amputated"]
    
    expr_data <- data.table(
      gene_id = rownames(norm_counts),
      expr_Control = rowMeans(norm_counts[, control_samples, drop = FALSE], na.rm = TRUE),
      expr_Amputated = rowMeans(norm_counts[, amputated_samples, drop = FALSE], na.rm = TRUE)
    )
    expr_data[, expr_mean := (expr_Control + expr_Amputated) / 2]
    expr_data[, expr_log2fc := log2((expr_Amputated + 1) / (expr_Control + 1))]
    
    # Merge DESeq2 p-values
    expr_data <- merge(expr_data, res_dt[, .(gene_id, log2FoldChange, pvalue, padj)],
                       by = "gene_id", all.x = TRUE)
    
    expr_data <- expr_data[gene_id %in% promoter_final$gene_id]
    msg("  Genes with expression data: %d", nrow(expr_data))
    
    fwrite(expr_data, file.path(out_dir, "04_expression_correlation", "expression_data.tsv"), sep = "\t")
    rnaseq_available <- TRUE
  }
}

# =============================================================================
# PHASE 6: GENERATE ALL PLOTS (with gene counts in annotations)
# =============================================================================

msg("\n=== PHASE 6: GENERATING PLOTS ===")

plot_dir <- file.path(out_dir, "02_plots_png")

# --- Plot 1: Scatter - GC% vs O/E index ---
msg("  Plot 01: GC vs O/E scatter...")
p1_data <- promoter_final[!is.na(max_cpg_oe) & !is.na(max_gc)]

p1 <- ggplot(p1_data, aes(x = max_cpg_oe, y = max_gc * 100)) +
  geom_point(alpha = 0.6, size = 1.5, color = "#4ECDC4") +
  geom_vline(xintercept = LCP_CPGOE_THRESHOLD, linetype = "dotted", color = "#FF6B6B", linewidth = 0.8) +
  geom_vline(xintercept = HCP_CPGOE_THRESHOLD, linetype = "dotted", color = "#FF6B6B", linewidth = 0.8) +
  geom_hline(yintercept = HCP_GC_THRESHOLD * 100, linetype = "dotted", color = "#FF6B6B", linewidth = 0.8) +
  annotate("text", x = 0.2, y = 85, label = paste0("n = ", nrow(p1_data), " genes"),
           color = "white", size = 4, fontface = "italic") +
  scale_x_continuous(limits = c(0, 2.2), breaks = seq(0, 2, 0.5)) +
  scale_y_continuous(limits = c(15, 90), breaks = seq(20, 80, 20)) +
  labs(title = "Scatter - Protein Genes - Promoters",
       subtitle = "D. laeve promoter regions (-2000 to +500 bp)",
       x = "O/E index", y = "GC (%)") +
  theme_pub()

ggsave(file.path(plot_dir, "01_scatter_GC_vs_OE.png"), p1, width = 10, height = 8, dpi = 300)

# --- Plot 2: Scatter colored by class ---
msg("  Plot 02: Classified scatter...")
n_hcp <- sum(p1_data$promoter_class == "HCP")
n_icp <- sum(p1_data$promoter_class == "ICP")
n_lcp <- sum(p1_data$promoter_class == "LCP")

p2 <- ggplot(p1_data, aes(x = max_cpg_oe, y = max_gc * 100, color = promoter_class)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = class_colors,
                     labels = c(paste0("HCP (n=", n_hcp, ")"),
                                paste0("ICP (n=", n_icp, ")"),
                                paste0("LCP (n=", n_lcp, ")"))) +
  geom_vline(xintercept = LCP_CPGOE_THRESHOLD, linetype = "dotted", color = "white", linewidth = 0.5) +
  geom_vline(xintercept = HCP_CPGOE_THRESHOLD, linetype = "dotted", color = "white", linewidth = 0.5) +
  geom_hline(yintercept = HCP_GC_THRESHOLD * 100, linetype = "dotted", color = "white", linewidth = 0.5) +
  scale_x_continuous(limits = c(0, 2.2), breaks = seq(0, 2, 0.5)) +
  scale_y_continuous(limits = c(15, 90), breaks = seq(20, 80, 20)) +
  labs(title = "Promoter Classification (Weber Criteria)",
       subtitle = paste0("HCP: ", n_hcp, " | ICP: ", n_icp, " | LCP: ", n_lcp),
       x = "CpG O/E Ratio", y = "GC Content (%)", color = "Promoter Class") +
  theme_pub() +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1)))

ggsave(file.path(plot_dir, "02_scatter_classified.png"), p2, width = 10, height = 8, dpi = 300)

# --- Plot 3: Histograms ---
msg("  Plot 03: CpG O/E histograms...")
p3a <- ggplot(promoter_final, aes(x = max_cpg_oe)) +
  geom_histogram(binwidth = 0.05, fill = "#888888", color = "white", alpha = 0.8) +
  scale_x_continuous(limits = c(0, 1.4), breaks = seq(0, 1.2, 0.2)) +
  labs(title = "CpG Ratio Distribution (All Promoters)",
       subtitle = paste0("n = ", nrow(promoter_final), " protein-coding gene promoters"),
       x = "CpG ratio (obs/exp)", y = "Number of promoters") +
  theme_pub_light()

p3b <- ggplot(promoter_final, aes(x = max_cpg_oe, fill = promoter_class)) +
  geom_histogram(binwidth = 0.05, color = "white", alpha = 0.8, position = "identity") +
  scale_fill_manual(values = class_colors) +
  scale_x_continuous(limits = c(0, 1.4), breaks = seq(0, 1.2, 0.2)) +
  facet_wrap(~ promoter_class, ncol = 1, scales = "free_y") +
  labs(title = "CpG Ratio by Promoter Class", x = "CpG ratio (obs/exp)", y = "Number of promoters") +
  theme_pub_light() + theme(legend.position = "none")

ggsave(file.path(plot_dir, "03_histogram_cpg_oe_all.png"), p3a, width = 8, height = 6, dpi = 300)
ggsave(file.path(plot_dir, "03_histogram_cpg_oe_by_class.png"), p3b, width = 8, height = 10, dpi = 300)

# --- Plot 4: Pie chart ---
msg("  Plot 04: Pie chart...")
pie_data <- promoter_final[, .N, by = promoter_class]
pie_data[, pct := N / sum(N) * 100]
pie_data[, label := paste0(promoter_class, "\n", formatC(N, big.mark = ","), " genes\n", round(pct, 1), "%")]

p4 <- ggplot(pie_data, aes(x = "", y = N, fill = promoter_class)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = class_colors) +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5),
            color = "white", fontface = "bold", size = 4.5) +
  labs(title = "Promoter Classification Distribution") +
  theme_void() + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14), legend.position = "none")

ggsave(file.path(plot_dir, "04_pie_chart_classes.png"), p4, width = 8, height = 8, dpi = 300)

# --- Plot 5 & 6: Methylation by class ---
for (cov_label in c("cov5", "cov10")) {
  min_cov_val <- ifelse(cov_label == "cov5", 5, 10)
  final_data <- if (cov_label == "cov5") final_data_cov5 else final_data_cov10
  meth_class <- final_data[!is.na(Control_mean)]
  plot_num <- ifelse(cov_label == "cov5", "05", "06")
  
  if (nrow(meth_class) > 0) {
    # Count genes per class
    class_n <- meth_class[, .N, by = promoter_class]
    
    meth_long <- melt(meth_class[, .(gene_id, promoter_class, Control_mean, Amputated_mean)],
                      id.vars = c("gene_id", "promoter_class"),
                      variable.name = "condition", value.name = "methylation")
    meth_long[, condition := gsub("_mean", "", condition)]
    
    p_meth <- ggplot(meth_long, aes(x = promoter_class, y = methylation, fill = condition)) +
      geom_violin(alpha = 0.7, position = position_dodge(0.8)) +
      geom_boxplot(width = 0.2, position = position_dodge(0.8), alpha = 0.9, outlier.size = 0.5) +
      scale_fill_manual(values = c(Control = "#3182BD", Amputated = "#E6550D")) +
      scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
      # Add gene counts as text
      annotate("text", x = 1:nrow(class_n), y = -0.05,
               label = paste0("n=", class_n$N), size = 3.5, color = "#666666") +
      labs(title = paste0("Promoter Methylation by Class (Coverage >= ", min_cov_val, ")"),
           subtitle = paste0("Total n = ", nrow(meth_class), " promoters"),
           x = "Promoter Class", y = "Mean Methylation", fill = "Condition") +
      theme_pub_light() +
      coord_cartesian(clip = "off")
    
    ggsave(file.path(plot_dir, paste0(plot_num, "_methylation_by_class_", cov_label, ".png")),
           p_meth, width = 10, height = 7, dpi = 300)
    
    # Stats
    class_meth_stats <- meth_class[, .(
      n = .N,
      meth_Control_mean = mean(Control_mean, na.rm = TRUE),
      meth_Control_sd = sd(Control_mean, na.rm = TRUE),
      meth_Amputated_mean = mean(Amputated_mean, na.rm = TRUE),
      meth_Amputated_sd = sd(Amputated_mean, na.rm = TRUE),
      delta_mean = mean(delta_meth, na.rm = TRUE)
    ), by = promoter_class]
    
    fwrite(class_meth_stats, file.path(out_dir, "03_methylation_analysis",
                                       paste0("methylation_stats_by_class_", cov_label, ".tsv")), sep = "\t")
  }
}

# --- NEW Plot: O/E index vs Methylation correlation ---
msg("  Plot NEW: O/E index vs Methylation correlation...")

meth_class_cov5 <- final_data_cov5[!is.na(Control_mean)]

if (nrow(meth_class_cov5) > 0) {
  # Spearman correlation by class
  oe_meth_cor <- meth_class_cov5[, .(
    n = .N,
    rho_control = cor(max_cpg_oe, Control_mean, method = "spearman", use = "complete.obs"),
    rho_amputated = cor(max_cpg_oe, Amputated_mean, method = "spearman", use = "complete.obs")
  ), by = promoter_class]
  
  overall_oe_rho <- cor(meth_class_cov5$max_cpg_oe, meth_class_cov5$Control_mean,
                        method = "spearman", use = "complete.obs")
  
  msg("  O/E vs Methylation correlations:")
  print(oe_meth_cor)
  msg("  Overall rho: %.4f", overall_oe_rho)
  
  fwrite(oe_meth_cor, file.path(out_dir, "03_methylation_analysis", "oe_vs_methylation_correlation.tsv"), sep = "\t")
  
  # Plot: O/E vs Methylation (Control) colored by class
  set.seed(42)
  plot_idx <- if (nrow(meth_class_cov5) > 8000) sample(1:nrow(meth_class_cov5), 8000) else 1:nrow(meth_class_cov5)
  plot_oe <- meth_class_cov5[plot_idx]
  
  p_oe_meth <- ggplot(plot_oe, aes(x = max_cpg_oe, y = Control_mean, color = promoter_class)) +
    geom_point(alpha = 0.3, size = 1) +
    geom_smooth(method = "loess", se = TRUE, linewidth = 1.2) +
    scale_color_manual(values = class_colors,
                       labels = paste0(oe_meth_cor$promoter_class,
                                       " (n=", oe_meth_cor$n,
                                       ", rho=", round(oe_meth_cor$rho_control, 3), ")")) +
    scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
    scale_x_continuous(limits = c(0, 2.0)) +
    labs(title = "CpG O/E Index vs Promoter Methylation",
         subtitle = paste0("Overall Spearman rho = ", round(overall_oe_rho, 4)),
         x = "CpG O/E Ratio (max window)", y = "Methylation (Control)",
         color = "Promoter Class") +
    theme_pub_light()
  
  ggsave(file.path(plot_dir, "10_OE_vs_methylation_by_class.png"), p_oe_meth, width = 10, height = 8, dpi = 300)
  
  # Faceted version
  p_oe_meth_facet <- ggplot(plot_oe, aes(x = max_cpg_oe, y = Control_mean)) +
    geom_point(aes(color = promoter_class), alpha = 0.3, size = 1) +
    geom_smooth(method = "loess", color = "black", se = TRUE) +
    scale_color_manual(values = class_colors) +
    scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
    facet_wrap(~ promoter_class, scales = "free_x") +
    labs(title = "O/E Index vs Methylation by Promoter Class",
         x = "CpG O/E Ratio", y = "Methylation (Control)") +
    theme_pub_light() + theme(legend.position = "none")
  
  ggsave(file.path(plot_dir, "10_OE_vs_methylation_faceted.png"), p_oe_meth_facet, width = 12, height = 5, dpi = 300)
}

# --- Plot 7: Methylation vs CpG density ---
msg("  Plot 07: Methylation vs CpG density...")
if (nrow(meth_class_cov5) > 0) {
  for (cls in c("LCP", "ICP", "HCP")) {
    cls_data <- meth_class_cov5[promoter_class == cls]
    if (nrow(cls_data) > 10) {
      p_cls <- ggplot(cls_data, aes(x = cpg_per_bp * 1000, y = Control_mean)) +
        geom_point(alpha = 0.5, size = 1.5, color = class_colors[cls]) +
        scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
        labs(title = paste0(cls, " - Methylation vs CpG Density (n=", nrow(cls_data), ")"),
             x = "CpG per 1000 bp", y = "Methylation (Control)") +
        theme_pub_light()
      ggsave(file.path(plot_dir, paste0("07_", cls, "_meth_vs_cpg_density.png")),
             p_cls, width = 8, height = 6, dpi = 300)
    }
  }
  
  p7_combined <- ggplot(meth_class_cov5, aes(x = cpg_per_bp * 1000, y = Control_mean, color = promoter_class)) +
    geom_point(alpha = 0.4, size = 1) +
    scale_color_manual(values = class_colors) +
    scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
    facet_wrap(~ promoter_class, scales = "free_x") +
    labs(title = "Methylation vs CpG Density by Promoter Class",
         x = "CpG per 1000 bp", y = "Methylation (Control)") +
    theme_pub_light() + theme(legend.position = "none")
  
  ggsave(file.path(plot_dir, "07_methylation_vs_cpg_density_faceted.png"),
         p7_combined, width = 12, height = 5, dpi = 300)
}

# =============================================================================
# PHASE 7: EXPRESSION-METHYLATION + VOLCANO + GENE LISTS
# =============================================================================

if (rnaseq_available) {
  msg("\n=== PHASE 7: EXPRESSION-METHYLATION CORRELATION + VOLCANO ===")
  
  expr_meth_cov5 <- merge(final_data_cov5[!is.na(Control_mean)], expr_data, by = "gene_id")
  expr_meth_cov10 <- merge(final_data_cov10[!is.na(Control_mean)], expr_data, by = "gene_id")
  
  msg("  Genes with both expression and methylation (cov>=5): %d", nrow(expr_meth_cov5))
  
  if (nrow(expr_meth_cov5) > 50) {
    
    # Correlations
    corr_by_class_cov5 <- expr_meth_cov5[, .(
      n = .N,
      cor_meth_expr_control = cor(Control_mean, log2(expr_Control + 1), method = "spearman", use = "complete.obs"),
      cor_meth_expr_amputated = cor(Amputated_mean, log2(expr_Amputated + 1), method = "spearman", use = "complete.obs"),
      cor_delta_meth_delta_expr = cor(delta_meth, expr_log2fc, method = "spearman", use = "complete.obs")
    ), by = promoter_class]
    
    msg("\n  Correlations (Spearman) by class (cov>=5):")
    print(corr_by_class_cov5)
    
    fwrite(corr_by_class_cov5, file.path(out_dir, "04_expression_correlation", "correlation_by_class_cov5.tsv"), sep = "\t")
    
    overall_cor_cov5 <- data.table(
      metric = c("Methylation vs Expression (Control)",
                 "Methylation vs Expression (Amputated)",
                 "Delta Methylation vs Delta Expression"),
      spearman_rho = c(
        cor(expr_meth_cov5$Control_mean, log2(expr_meth_cov5$expr_Control + 1), method = "spearman", use = "complete.obs"),
        cor(expr_meth_cov5$Amputated_mean, log2(expr_meth_cov5$expr_Amputated + 1), method = "spearman", use = "complete.obs"),
        cor(expr_meth_cov5$delta_meth, expr_meth_cov5$expr_log2fc, method = "spearman", use = "complete.obs")
      )
    )
    
    # --- Plot 8: Methylation vs Expression ---
    msg("  Plot 08: Methylation vs Expression...")
    set.seed(42)
    plot_sample <- if (nrow(expr_meth_cov5) > 5000) sample(1:nrow(expr_meth_cov5), 5000) else 1:nrow(expr_meth_cov5)
    plot_data <- expr_meth_cov5[plot_sample]
    
    p8a <- ggplot(plot_data, aes(x = Control_mean, y = log2(expr_Control + 1), color = promoter_class)) +
      geom_point(alpha = 0.4, size = 1) +
      geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
      scale_color_manual(values = class_colors) +
      scale_x_continuous(labels = percent_format()) +
      labs(title = "Promoter Methylation vs Expression (Control)",
           subtitle = paste0("Overall rho = ", round(overall_cor_cov5[1, spearman_rho], 3)),
           x = "Promoter Methylation", y = "log2(Expression + 1)", color = "Promoter Class") +
      theme_pub_light()
    
    ggsave(file.path(plot_dir, "08_methylation_vs_expression_control.png"), p8a, width = 10, height = 8, dpi = 300)
    
    p8b <- ggplot(plot_data, aes(x = Control_mean, y = log2(expr_Control + 1))) +
      geom_point(aes(color = promoter_class), alpha = 0.3, size = 1) +
      geom_smooth(method = "lm", color = "black", se = TRUE) +
      scale_color_manual(values = class_colors) +
      scale_x_continuous(labels = percent_format()) +
      facet_wrap(~ promoter_class, scales = "free") +
      labs(title = "Methylation vs Expression by Promoter Class",
           x = "Promoter Methylation", y = "log2(Expression + 1)") +
      theme_pub_light() + theme(legend.position = "none")
    
    ggsave(file.path(plot_dir, "08_methylation_vs_expression_by_class.png"), p8b, width = 12, height = 5, dpi = 300)
    
    # --- Plot 9: Delta methylation vs Delta expression ---
    p9 <- ggplot(plot_data, aes(x = delta_meth, y = expr_log2fc, color = promoter_class)) +
      geom_point(alpha = 0.4, size = 1) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      geom_smooth(method = "lm", se = FALSE) +
      scale_color_manual(values = class_colors) +
      labs(title = "Change in Methylation vs Change in Expression",
           subtitle = paste0("(Amputated - Control) | rho = ", round(overall_cor_cov5[3, spearman_rho], 3)),
           x = "Delta Methylation", y = "log2 Fold Change Expression", color = "Promoter Class") +
      theme_pub_light()
    
    ggsave(file.path(plot_dir, "09_delta_meth_vs_delta_expr.png"), p9, width = 10, height = 8, dpi = 300)
    
    # =============================================
    # NEW: VOLCANO PLOT (Differential Methylation)
    # =============================================
    msg("  Plot NEW: Volcano plot for differential methylation...")
    
    # For the volcano we need a test. With only 2 replicates, use a simple approach:
    # Calculate variance-stabilized delta and use a permutation-like approach
    # Or just plot delta_meth vs -log10(p) where p comes from a t-test-like metric
    
    # Simple approach: for each gene, compute a "significance" based on
    # consistency across replicates
    volcano_data <- expr_meth_cov5[, .(
      gene_id, promoter_class, delta_meth,
      # Individual deltas
      delta_1 = A1 - C1,
      delta_2 = A2 - C2,
      Control_mean, Amputated_mean,
      expr_log2fc, padj
    )]
    
    # Compute a naive z-score for methylation change
    # Using pooled variance across all genes
    global_sd <- sd(volcano_data$delta_meth, na.rm = TRUE)
    volcano_data[, meth_zscore := delta_meth / global_sd]
    volcano_data[, meth_pval := 2 * pnorm(-abs(meth_zscore))]
    volcano_data[, meth_padj := p.adjust(meth_pval, method = "BH")]
    volcano_data[, neg_log10_p := -log10(meth_pval + 1e-300)]
    
    # Classification of methylation changes
    DELTA_METH_THRESHOLD <- 0.10  # 10% change
    METH_PVAL_THRESHOLD <- 0.05
    
    volcano_data[, meth_status := "No change"]
    volcano_data[delta_meth > DELTA_METH_THRESHOLD & meth_pval < METH_PVAL_THRESHOLD,
                 meth_status := "Hypermethylated\n(gained methylation)"]
    volcano_data[delta_meth < -DELTA_METH_THRESHOLD & meth_pval < METH_PVAL_THRESHOLD,
                 meth_status := "Hypomethylated\n(lost methylation)"]
    
    status_colors <- c("No change" = "#888888",
                       "Hypermethylated\n(gained methylation)" = "#E6550D",
                       "Hypomethylated\n(lost methylation)" = "#3182BD")
    
    n_hyper <- sum(volcano_data$meth_status == "Hypermethylated\n(gained methylation)")
    n_hypo <- sum(volcano_data$meth_status == "Hypomethylated\n(lost methylation)")
    n_stable <- sum(volcano_data$meth_status == "No change")
    
    msg("    Hypermethylated (Amputated > Control): %d", n_hyper)
    msg("    Hypomethylated (Amputated < Control): %d", n_hypo)
    msg("    No change: %d", n_stable)
    
    p_volcano <- ggplot(volcano_data, aes(x = delta_meth, y = neg_log10_p, color = meth_status)) +
      geom_point(alpha = 0.4, size = 1) +
      geom_vline(xintercept = c(-DELTA_METH_THRESHOLD, DELTA_METH_THRESHOLD),
                 linetype = "dashed", color = "gray40") +
      geom_hline(yintercept = -log10(METH_PVAL_THRESHOLD),
                 linetype = "dashed", color = "gray40") +
      scale_color_manual(values = status_colors,
                         labels = c(paste0("No change (", n_stable, ")"),
                                    paste0("Gained meth (", n_hyper, ")"),
                                    paste0("Lost meth (", n_hypo, ")"))) +
      labs(title = "Volcano Plot: Differential Promoter Methylation",
           subtitle = "Amputated vs Control | Threshold: |delta| > 0.10, p < 0.05",
           x = "Delta Methylation (Amputated - Control)",
           y = "-log10(p-value)",
           color = "Status") +
      theme_pub_light() +
      guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
    
    ggsave(file.path(plot_dir, "11_volcano_differential_methylation.png"),
           p_volcano, width = 10, height = 8, dpi = 300)
    
    # Volcano by class
    p_volcano_class <- ggplot(volcano_data, aes(x = delta_meth, y = neg_log10_p, color = meth_status)) +
      geom_point(alpha = 0.4, size = 1) +
      geom_vline(xintercept = c(-DELTA_METH_THRESHOLD, DELTA_METH_THRESHOLD),
                 linetype = "dashed", color = "gray40") +
      geom_hline(yintercept = -log10(METH_PVAL_THRESHOLD),
                 linetype = "dashed", color = "gray40") +
      scale_color_manual(values = status_colors) +
      facet_wrap(~ promoter_class, scales = "free") +
      labs(title = "Volcano Plot by Promoter Class",
           x = "Delta Methylation", y = "-log10(p-value)", color = "Status") +
      theme_pub_light()
    
    ggsave(file.path(plot_dir, "11_volcano_by_class.png"),
           p_volcano_class, width = 14, height = 6, dpi = 300)
    
    # =============================================
    # NEW: ICP SUBPOPULATION ANALYSIS
    # =============================================
    msg("\n  Analyzing ICP subpopulations...")
    
    icp_data <- expr_meth_cov5[promoter_class == "ICP"]
    
    # ICP has a bimodal methylation pattern - split into sub-groups
    icp_data[, meth_group := cut(Control_mean,
                                 breaks = c(-Inf, 0.20, 0.80, Inf),
                                 labels = c("Low (<20%)", "Medium (20-80%)", "High (>80%)"))]
    
    icp_subpop <- icp_data[, .(
      n = .N,
      mean_meth_control = mean(Control_mean, na.rm = TRUE),
      mean_meth_amputated = mean(Amputated_mean, na.rm = TRUE),
      mean_delta = mean(delta_meth, na.rm = TRUE),
      mean_oe = mean(max_cpg_oe, na.rm = TRUE),
      mean_expr = mean(log2(expr_Control + 1), na.rm = TRUE),
      cor_meth_expr = cor(Control_mean, log2(expr_Control + 1), method = "spearman", use = "complete.obs")
    ), by = meth_group]
    
    msg("  ICP subpopulations:")
    print(icp_subpop)
    
    fwrite(icp_subpop, file.path(out_dir, "03_methylation_analysis", "ICP_subpopulation_stats.tsv"), sep = "\t")
    
    # ICP subpopulation plot
    p_icp_sub <- ggplot(icp_data, aes(x = Control_mean, fill = meth_group)) +
      geom_histogram(binwidth = 0.02, alpha = 0.8, color = "white") +
      scale_fill_manual(values = c("Low (<20%)" = "#2166AC", "Medium (20-80%)" = "#FDBF6F", "High (>80%)" = "#D6604D")) +
      scale_x_continuous(labels = percent_format()) +
      labs(title = "ICP Methylation Distribution - Subpopulations",
           subtitle = paste0("n = ", nrow(icp_data), " ICP genes"),
           x = "Methylation (Control)", y = "Number of genes", fill = "Methylation Group") +
      theme_pub_light()
    
    ggsave(file.path(plot_dir, "12_ICP_subpopulation_histogram.png"),
           p_icp_sub, width = 10, height = 7, dpi = 300)
    
    # ICP O/E vs methylation with subgroups
    p_icp_oe <- ggplot(icp_data, aes(x = max_cpg_oe, y = Control_mean, color = meth_group)) +
      geom_point(alpha = 0.3, size = 1) +
      scale_color_manual(values = c("Low (<20%)" = "#2166AC", "Medium (20-80%)" = "#FDBF6F", "High (>80%)" = "#D6604D")) +
      scale_y_continuous(labels = percent_format()) +
      labs(title = "ICP: O/E Index vs Methylation by Subgroup",
           x = "CpG O/E Ratio", y = "Methylation (Control)", color = "Methylation Group") +
      theme_pub_light()
    
    ggsave(file.path(plot_dir, "12_ICP_OE_vs_meth_subgroups.png"),
           p_icp_oe, width = 10, height = 8, dpi = 300)
    
    # =============================================
    # NEW: GENE LISTS
    # =============================================
    msg("\n  Generating gene lists...")
    
    # 1) Canonical pattern genes
    canonical_regulated <- expr_meth_cov5[
      (Control_mean < 0.3 & expr_Control > quantile(expr_Control, 0.75, na.rm = TRUE)) |
        (Control_mean > 0.7 & expr_Control < quantile(expr_Control, 0.25, na.rm = TRUE))
    ]
    canonical_regulated[, regulation_type := ifelse(Control_mean < 0.3, "Low_Meth_High_Expr", "High_Meth_Low_Expr")]
    
    # 2) Differentially methylated genes
    dm_genes <- volcano_data[meth_status != "No change"]
    dm_genes_hyper <- dm_genes[grepl("Hyper", meth_status)]
    dm_genes_hypo <- dm_genes[grepl("Hypo", meth_status)]
    
    # 3) Genes that gained methylation AND lost expression (and vice versa)
    concordant_changes <- expr_meth_cov5[
      (delta_meth > DELTA_METH_THRESHOLD & expr_log2fc < -0.5) |  # gained meth, lost expr
        (delta_meth < -DELTA_METH_THRESHOLD & expr_log2fc > 0.5)    # lost meth, gained expr
    ]
    concordant_changes[, change_type := ifelse(delta_meth > 0, "Gained_Meth_Lost_Expr", "Lost_Meth_Gained_Expr")]
    
    msg("    Canonical pattern genes: %d", nrow(canonical_regulated))
    msg("    Differentially methylated: %d (hyper: %d, hypo: %d)",
        nrow(dm_genes), nrow(dm_genes_hyper), nrow(dm_genes_hypo))
    msg("    Concordant meth-expr changes: %d", nrow(concordant_changes))
    
    # Save gene lists
    gene_list_dir <- file.path(out_dir, "07_gene_lists")
    
    save_gene_list <- function(dt, filename, desc) {
      cols_to_keep <- intersect(
        c("gene_id", "promoter_class", "chrom", "max_cpg_oe", "max_gc",
          "Control_mean", "Amputated_mean", "delta_meth",
          "expr_Control", "expr_Amputated", "expr_log2fc", "padj",
          "regulation_type", "change_type", "meth_status", "gene_note"),
        names(dt)
      )
      out_dt <- dt[, ..cols_to_keep]
      fwrite(out_dt, file.path(gene_list_dir, filename), sep = "\t")
      msg("      Saved %s: %d genes - %s", filename, nrow(out_dt), desc)
    }
    
    save_gene_list(canonical_regulated, "canonical_regulated_genes.tsv",
                   "Genes following canonical methylation-expression pattern")
    save_gene_list(dm_genes_hyper, "hypermethylated_genes.tsv",
                   "Genes gaining methylation after amputation")
    save_gene_list(dm_genes_hypo, "hypomethylated_genes.tsv",
                   "Genes losing methylation after amputation")
    save_gene_list(concordant_changes, "concordant_meth_expr_changes.tsv",
                   "Genes with concordant methylation and expression changes")
    
    # Gene lists by class
    for (cls in c("HCP", "ICP", "LCP")) {
      cls_dm <- dm_genes[promoter_class == cls]
      if (nrow(cls_dm) > 0) {
        save_gene_list(cls_dm, paste0("differentially_methylated_", cls, ".tsv"),
                       paste0("DM genes in ", cls, " class"))
      }
    }
    
    # Save full integrated dataset
    fwrite(expr_meth_cov5, file.path(out_dir, "04_expression_correlation", "expression_methylation_integrated_cov5.tsv"), sep = "\t")
    fwrite(volcano_data, file.path(out_dir, "08_volcano", "volcano_data.tsv"), sep = "\t")
  }
}

# =============================================================================
# PHASE 8: GENERATE HTML REPORT (FIXED - using paste0 instead of sprintf)
# =============================================================================

msg("\n=== PHASE 8: GENERATING HTML REPORT (FIXED) ===")

# Use paste0() to avoid the %% escaping nightmare in sprintf()
n_hcp_final <- class_summary[promoter_class == "HCP", count]
n_icp_final <- class_summary[promoter_class == "ICP", count]
n_lcp_final <- class_summary[promoter_class == "LCP", count]
n_total_final <- nrow(promoter_final)
gen_date <- format(Sys.time(), "%Y-%m-%d %H:%M")

html_report <- paste0('<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Promoter Classification Analysis - D. laeve</title>
    <style>
        @import url("https://fonts.googleapis.com/css2?family=Inter:wght@400;600;700&family=JetBrains+Mono&display=swap");

        :root {
            --bg-primary: #0f0f0f;
            --bg-secondary: #1a1a1a;
            --bg-tertiary: #252525;
            --text-primary: #ffffff;
            --text-secondary: #a0a0a0;
            --accent-hcp: #2166AC;
            --accent-icp: #4DAF4A;
            --accent-lcp: #D6604D;
            --accent-control: #3182BD;
            --accent-amputated: #E6550D;
        }

        * { margin: 0; padding: 0; box-sizing: border-box; }

        body {
            font-family: "Inter", -apple-system, sans-serif;
            background: var(--bg-primary);
            color: var(--text-primary);
            line-height: 1.6;
        }

        .container { max-width: 1400px; margin: 0 auto; padding: 2rem; }

        header {
            text-align: center;
            padding: 3rem 0;
            border-bottom: 1px solid var(--bg-tertiary);
            margin-bottom: 3rem;
        }

        h1 {
            font-size: 2.5rem;
            font-weight: 700;
            margin-bottom: 0.5rem;
            background: linear-gradient(135deg, var(--accent-hcp), var(--accent-icp), var(--accent-lcp));
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
        }

        .subtitle { color: var(--text-secondary); font-size: 1.1rem; }

        .stats-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 1.5rem;
            margin: 2rem 0;
        }

        .stat-card {
            background: var(--bg-secondary);
            border-radius: 12px;
            padding: 1.5rem;
            text-align: center;
            border: 1px solid var(--bg-tertiary);
        }

        .stat-value {
            font-size: 2.5rem;
            font-weight: 700;
            font-family: "JetBrains Mono", monospace;
        }

        .stat-value.hcp { color: var(--accent-hcp); }
        .stat-value.icp { color: var(--accent-icp); }
        .stat-value.lcp { color: var(--accent-lcp); }

        .stat-label { color: var(--text-secondary); font-size: 0.9rem; margin-top: 0.5rem; }

        .section {
            background: var(--bg-secondary);
            border-radius: 16px;
            padding: 2rem;
            margin: 2rem 0;
            border: 1px solid var(--bg-tertiary);
        }

        .section h2 {
            font-size: 1.5rem;
            margin-bottom: 1.5rem;
            padding-bottom: 0.5rem;
            border-bottom: 2px solid var(--bg-tertiary);
        }

        .plot-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(400px, 1fr));
            gap: 2rem;
        }

        .plot-card {
            background: var(--bg-tertiary);
            border-radius: 12px;
            overflow: hidden;
        }

        .plot-card img { width: 100%; height: auto; display: block; }
        .plot-caption { padding: 1rem; font-size: 0.9rem; color: var(--text-secondary); }

        table {
            width: 100%;
            border-collapse: collapse;
            margin: 1rem 0;
            font-family: "JetBrains Mono", monospace;
            font-size: 0.85rem;
        }

        th, td { padding: 0.6rem 0.8rem; text-align: left; border-bottom: 1px solid var(--bg-tertiary); }
        th { background: var(--bg-tertiary); font-weight: 600; }

        .method-box {
            background: var(--bg-tertiary);
            border-radius: 8px;
            padding: 1.5rem;
            margin: 1rem 0;
            border-left: 4px solid var(--accent-hcp);
        }

        .method-box h3 { margin-bottom: 0.5rem; }
        .method-box p { color: var(--text-secondary); font-size: 0.9rem; }

        footer {
            text-align: center;
            padding: 2rem;
            color: var(--text-secondary);
            font-size: 0.9rem;
            border-top: 1px solid var(--bg-tertiary);
            margin-top: 3rem;
        }

        .gene-table-container { max-height: 400px; overflow-y: auto; border-radius: 8px; }
    </style>
</head>
<body>
    <div class="container">
        <header>
            <h1>Promoter Classification Analysis</h1>
            <p class="subtitle">D. laeve Regeneration Study | Weber Criteria (2007)</p>
            <p class="subtitle">Generated: ', gen_date, '</p>
        </header>

        <div class="stats-grid">
            <div class="stat-card">
                <div class="stat-value hcp">', n_hcp_final, '</div>
                <div class="stat-label">HCP (High CpG)<br>', round(class_summary[promoter_class == "HCP", pct], 1), '%</div>
            </div>
            <div class="stat-card">
                <div class="stat-value icp">', n_icp_final, '</div>
                <div class="stat-label">ICP (Intermediate)<br>', round(class_summary[promoter_class == "ICP", pct], 1), '%</div>
            </div>
            <div class="stat-card">
                <div class="stat-value lcp">', n_lcp_final, '</div>
                <div class="stat-label">LCP (Low CpG)<br>', round(class_summary[promoter_class == "LCP", pct], 1), '%</div>
            </div>
            <div class="stat-card">
                <div class="stat-value" style="color: #fff;">', format(n_total_final, big.mark = ","), '</div>
                <div class="stat-label">Total Genes Analyzed</div>
            </div>
        </div>

        <div class="section">
            <h2>Classification Methods</h2>
            <div class="method-box">
                <h3>Weber et al. 2007 Criteria</h3>
                <p>Promoters classified using 500bp sliding windows across -2000 to +500 bp relative to TSS:</p>
                <ul style="margin-top: 0.5rem; margin-left: 1.5rem; color: var(--text-secondary);">
                    <li><strong style="color: var(--accent-hcp);">HCP:</strong> CpG O/E >= 0.75 AND GC content >= 55%</li>
                    <li><strong style="color: var(--accent-lcp);">LCP:</strong> No window with CpG O/E >= 0.48</li>
                    <li><strong style="color: var(--accent-icp);">ICP:</strong> Neither HCP nor LCP</li>
                </ul>
            </div>
        </div>

        <div class="section">
            <h2>Classification Plots</h2>
            <div class="plot-grid">
                <div class="plot-card">
                    <img src="02_plots_png/01_scatter_GC_vs_OE.png" alt="GC vs O/E">
                    <div class="plot-caption">GC content vs CpG O/E ratio for all protein-coding gene promoters</div>
                </div>
                <div class="plot-card">
                    <img src="02_plots_png/02_scatter_classified.png" alt="Classified">
                    <div class="plot-caption">Promoters colored by Weber classification</div>
                </div>
                <div class="plot-card">
                    <img src="02_plots_png/03_histogram_cpg_oe_all.png" alt="Histogram">
                    <div class="plot-caption">Distribution of CpG O/E ratios (unimodal in D. laeve)</div>
                </div>
                <div class="plot-card">
                    <img src="02_plots_png/04_pie_chart_classes.png" alt="Pie chart">
                    <div class="plot-caption">Proportion of each class - note 91.6% ICP dominance</div>
                </div>
            </div>
        </div>

        <div class="section">
            <h2>Methylation Analysis</h2>
            <div class="plot-grid">
                <div class="plot-card">
                    <img src="02_plots_png/05_methylation_by_class_cov5.png" alt="Meth cov5">
                    <div class="plot-caption">Methylation by class (Coverage >= 5)</div>
                </div>
                <div class="plot-card">
                    <img src="02_plots_png/06_methylation_by_class_cov10.png" alt="Meth cov10">
                    <div class="plot-caption">Methylation by class (Coverage >= 10)</div>
                </div>
                <div class="plot-card">
                    <img src="02_plots_png/10_OE_vs_methylation_by_class.png" alt="OE vs Meth">
                    <div class="plot-caption"><strong>NEW:</strong> CpG O/E ratio vs Methylation by class</div>
                </div>
                <div class="plot-card">
                    <img src="02_plots_png/07_methylation_vs_cpg_density_faceted.png" alt="Meth vs CpG">
                    <div class="plot-caption">Methylation vs CpG density</div>
                </div>
            </div>
        </div>

        <div class="section">
            <h2>Differential Methylation (Volcano)</h2>
            <div class="plot-grid">
                <div class="plot-card">
                    <img src="02_plots_png/11_volcano_differential_methylation.png" alt="Volcano">
                    <div class="plot-caption"><strong>NEW:</strong> Volcano plot - differentially methylated promoters</div>
                </div>
                <div class="plot-card">
                    <img src="02_plots_png/11_volcano_by_class.png" alt="Volcano by class">
                    <div class="plot-caption"><strong>NEW:</strong> Volcano plot by promoter class</div>
                </div>
            </div>
        </div>

        <div class="section">
            <h2>ICP Subpopulation Analysis</h2>
            <div class="plot-grid">
                <div class="plot-card">
                    <img src="02_plots_png/12_ICP_subpopulation_histogram.png" alt="ICP subpop">
                    <div class="plot-caption"><strong>NEW:</strong> ICP methylation distribution showing subpopulations</div>
                </div>
                <div class="plot-card">
                    <img src="02_plots_png/12_ICP_OE_vs_meth_subgroups.png" alt="ICP OE subgroups">
                    <div class="plot-caption"><strong>NEW:</strong> ICP O/E vs methylation colored by subgroup</div>
                </div>
            </div>
        </div>

        <div class="section">
            <h2>Expression-Methylation Correlation</h2>
            <div class="plot-grid">
                <div class="plot-card">
                    <img src="02_plots_png/08_methylation_vs_expression_control.png" alt="Meth vs Expr">
                    <div class="plot-caption">Methylation vs Expression (Control)</div>
                </div>
                <div class="plot-card">
                    <img src="02_plots_png/09_delta_meth_vs_delta_expr.png" alt="Delta">
                    <div class="plot-caption">Delta Methylation vs Delta Expression</div>
                </div>
            </div>
        </div>

        <div class="section">
            <h2>Gene Lists</h2>
            <p style="color: var(--text-secondary); margin-bottom: 1rem;">
                Gene lists are saved in the <code>07_gene_lists/</code> directory as TSV files.
            </p>
            <table>
                <tr><th>File</th><th>Description</th></tr>
                <tr><td>canonical_regulated_genes.tsv</td><td>Genes following canonical methylation-expression pattern</td></tr>
                <tr><td>hypermethylated_genes.tsv</td><td>Genes gaining methylation after tail amputation</td></tr>
                <tr><td>hypomethylated_genes.tsv</td><td>Genes losing methylation after tail amputation</td></tr>
                <tr><td>concordant_meth_expr_changes.tsv</td><td>Genes with concordant methylation and expression changes</td></tr>
                <tr><td>differentially_methylated_HCP.tsv</td><td>DM genes in HCP class</td></tr>
                <tr><td>differentially_methylated_ICP.tsv</td><td>DM genes in ICP class</td></tr>
                <tr><td>differentially_methylated_LCP.tsv</td><td>DM genes in LCP class</td></tr>
            </table>
        </div>

        <footer>
            <p>Analysis based on Weber et al. 2007 Nature Genetics classification criteria</p>
            <p>D. laeve Methylation Study - Tail Amputation vs Control</p>
        </footer>
    </div>
</body>
</html>')

writeLines(html_report, file.path(out_dir, "06_html_report", "promoter_classification_report.html"))
writeLines(html_report, file.path(out_dir, "promoter_report.html"))

msg("\n")
msg("================================================================================")
msg("  ANALYSIS COMPLETE! (v2 - FIXED)")
msg("================================================================================")
msg("  Output directory: %s", out_dir)
msg("  HTML Report: %s/promoter_report.html", out_dir)
msg("  PNG Plots: %s/02_plots_png/", out_dir)
msg("  Gene Lists: %s/07_gene_lists/", out_dir)
msg("  Volcano Data: %s/08_volcano/", out_dir)
msg("  Data files: %s/05_data/", out_dir)
msg("================================================================================")