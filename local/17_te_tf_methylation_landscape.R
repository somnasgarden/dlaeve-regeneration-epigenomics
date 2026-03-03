#!/usr/bin/env Rscript
# =============================================================================
# TE Genomic Context, Methylation Landscape & TF Methylation Analysis
# =============================================================================
# Part 1: TE-overlapping DMPs/DMRs — genic vs intergenic breakdown
# Part 2: Genome methylation landscape with genomation (meta-gene profile)
# Part 3: DeepTFactor TF methylation — how many TFs are methylated & where
#
# Input:  results/01_methylation/tables/ (annotated DMPs/DMRs)
#         DATA/collapsed_te_age_data.tsv
#         DATA/derLaeGenome_namesDlasi_v2.fasta.functional_note.pseudo_label.gff
#         DATA/prediction_result.txt (DeepTFactor)
#         DATA/{C1,C2,A1,A2}.CpG_report.txt (Bismark)
# Output: results/20_te_tf_landscape/
# =============================================================================

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(RColorBrewer)
})

options(stringsAsFactors = FALSE)
options(scipen = 999)

theme_set(theme_minimal(base_size = 12) +
            theme(plot.title = element_text(size = 14, face = "bold")))

# ── Paths ────────────────────────────────────────────────────────────────────

# Waterfall: try local DATA path first, then cluster
DATA_DIR <- NULL
for (d in c("/mnt/c/Users/rafae/Projects/DATA",
            "C:/Users/rafae/Projects/DATA",
            "/mnt/data/alfredvar/rlopezt/Preliminary")) {
  if (dir.exists(d)) { DATA_DIR <- d; break }
}
if (is.null(DATA_DIR)) stop("Cannot find DATA directory")

RESULTS_BASE <- "results/01_methylation/tables"
OUTPUT_DIR   <- "results/20_te_tf_landscape"
dir.create(file.path(OUTPUT_DIR, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "plots"),  recursive = TRUE, showWarnings = FALSE)

GFF_FILE   <- file.path(DATA_DIR, "derLaeGenome_namesDlasi_v2.fasta.functional_note.pseudo_label.gff")
TE_FILE    <- file.path(DATA_DIR, "collapsed_te_age_data.tsv")
TF_FILE    <- file.path(DATA_DIR, "prediction_result.txt")
ANNOT_FILE <- file.path(DATA_DIR, "derLaeGenome_eviann_annotations.tsv")

save_both <- function(p, name, w = 10, h = 7) {
  ggsave(file.path(OUTPUT_DIR, "plots", paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(OUTPUT_DIR, "plots", paste0(name, ".png")), p, width = w, height = h, dpi = 300)
  cat("  Saved:", name, "\n")
}

cat("=== TE, TF & Methylation Landscape Analysis ===\n")
cat("DATA_DIR:", DATA_DIR, "\n")
cat("OUTPUT_DIR:", OUTPUT_DIR, "\n\n")

# Chromosome filter — chr1-31 only (exclude scaffolds)
keep_chr <- paste0("chr", 1:31)

# =============================================================================
# LOAD COMMON DATA
# =============================================================================

cat("--- Loading DMP/DMR annotations ---\n")

dmps <- fread(file.path(RESULTS_BASE, "dmps_annotated.txt"))
dmrs <- fread(file.path(RESULTS_BASE, "dmrs_annotated.txt"))

# Filter to chr1-31
dmps <- dmps[seqnames %in% keep_chr]
dmrs <- dmrs[seqnames %in% keep_chr]

cat("  DMPs:", nrow(dmps), "| DMRs:", nrow(dmrs), "(chr1-31 only)\n")
cat("  DMP regions:", paste(names(table(dmps$annotation)), collapse = ", "), "\n")
cat("  DMR regions:", paste(names(table(dmrs$annotation)), collapse = ", "), "\n\n")

# Build GRanges
dmp_gr <- GRanges(seqnames = dmps$seqnames,
                  ranges = IRanges(start = dmps$start, end = dmps$end),
                  annotation = dmps$annotation,
                  nearest_gene = dmps$nearest_gene,
                  methylation_diff = dmps$methylation_diff,
                  fdr = dmps$fdr)

dmr_gr <- GRanges(seqnames = dmrs$seqnames,
                  ranges = IRanges(start = dmrs$start, end = dmrs$end),
                  annotation = dmrs$annotation,
                  nearest_gene = dmrs$nearest_gene,
                  nCG = dmrs$nCG,
                  diff_methy = dmrs$diff_methy)


# =============================================================================
# PART 1: TE-OVERLAPPING DMPs/DMRs — GENIC vs INTERGENIC
# =============================================================================

cat("================================================================\n")
cat(" PART 1: TE-DMP/DMR Genomic Context                            \n")
cat("================================================================\n\n")

cat("Loading TE annotations...\n")
te_data <- fread(TE_FILE, sep = "\t", header = TRUE)
if ("chrom" %in% colnames(te_data)) setnames(te_data, "chrom", "chr")

# Filter to chr1-31
te_data <- te_data[chr %in% keep_chr]

# Parse class/family
te_data[, c("te_class", "te_family") := tstrsplit(class_family, "/", fixed = TRUE)]
te_data[is.na(te_family), te_family := te_class]

te_gr <- GRanges(seqnames = te_data$chr,
                 ranges = IRanges(start = te_data$start, end = te_data$end),
                 te_name = te_data$te_name,
                 te_class = te_data$te_class,
                 te_family = te_data$te_family)

cat("  TEs loaded (chr1-31):", length(te_gr), "\n\n")

# ── Find DMP-TE overlaps ──
cat("Finding DMP-TE overlaps...\n")
dmp_te_hits <- findOverlaps(dmp_gr, te_gr)
dmp_in_te_idx <- unique(queryHits(dmp_te_hits))
dmp_not_in_te_idx <- setdiff(seq_along(dmp_gr), dmp_in_te_idx)

cat("  DMPs overlapping TEs:", length(dmp_in_te_idx), "/", length(dmp_gr),
    sprintf("(%.1f%%)\n", 100 * length(dmp_in_te_idx) / length(dmp_gr)))
cat("  DMPs NOT in TEs:", length(dmp_not_in_te_idx), "/", length(dmp_gr),
    sprintf("(%.1f%%)\n\n", 100 * length(dmp_not_in_te_idx) / length(dmp_gr)))

# Genomic context of DMPs in TEs
dmp_te_context <- as.data.frame(table(mcols(dmp_gr[dmp_in_te_idx])$annotation))
colnames(dmp_te_context) <- c("Region", "Count")
dmp_te_context$Percentage <- round(100 * dmp_te_context$Count / sum(dmp_te_context$Count), 2)
dmp_te_context$Category <- "DMPs_in_TEs"

# Context of DMPs NOT in TEs (for comparison)
dmp_nonte_context <- as.data.frame(table(mcols(dmp_gr[dmp_not_in_te_idx])$annotation))
colnames(dmp_nonte_context) <- c("Region", "Count")
dmp_nonte_context$Percentage <- round(100 * dmp_nonte_context$Count / sum(dmp_nonte_context$Count), 2)
dmp_nonte_context$Category <- "DMPs_not_in_TEs"

# All DMPs for reference
dmp_all_context <- as.data.frame(table(mcols(dmp_gr)$annotation))
colnames(dmp_all_context) <- c("Region", "Count")
dmp_all_context$Percentage <- round(100 * dmp_all_context$Count / sum(dmp_all_context$Count), 2)
dmp_all_context$Category <- "All_DMPs"

dmp_context_combined <- rbind(dmp_te_context, dmp_nonte_context, dmp_all_context)

cat("DMP Genomic Context Breakdown:\n")
cat("\n  DMPs IN TEs:\n")
for (i in 1:nrow(dmp_te_context)) {
  cat(sprintf("    %-15s %5d (%5.1f%%)\n",
              dmp_te_context$Region[i], dmp_te_context$Count[i], dmp_te_context$Percentage[i]))
}
cat("\n  DMPs NOT in TEs:\n")
for (i in 1:nrow(dmp_nonte_context)) {
  cat(sprintf("    %-15s %5d (%5.1f%%)\n",
              dmp_nonte_context$Region[i], dmp_nonte_context$Count[i], dmp_nonte_context$Percentage[i]))
}
cat("\n  All DMPs:\n")
for (i in 1:nrow(dmp_all_context)) {
  cat(sprintf("    %-15s %5d (%5.1f%%)\n",
              dmp_all_context$Region[i], dmp_all_context$Count[i], dmp_all_context$Percentage[i]))
}

fwrite(dmp_context_combined, file.path(OUTPUT_DIR, "tables", "dmp_te_genomic_context.tsv"), sep = "\t")

# ── Find DMR-TE overlaps ──
cat("\nFinding DMR-TE overlaps...\n")
dmr_te_hits <- findOverlaps(dmr_gr, te_gr)
dmr_in_te_idx <- unique(queryHits(dmr_te_hits))
dmr_not_in_te_idx <- setdiff(seq_along(dmr_gr), dmr_in_te_idx)

cat("  DMRs overlapping TEs:", length(dmr_in_te_idx), "/", length(dmr_gr),
    sprintf("(%.1f%%)\n", 100 * length(dmr_in_te_idx) / length(dmr_gr)))

# DMR context
dmr_te_context <- as.data.frame(table(mcols(dmr_gr[dmr_in_te_idx])$annotation))
colnames(dmr_te_context) <- c("Region", "Count")
dmr_te_context$Percentage <- round(100 * dmr_te_context$Count / sum(dmr_te_context$Count), 2)
dmr_te_context$Category <- "DMRs_in_TEs"

dmr_nonte_context <- as.data.frame(table(mcols(dmr_gr[dmr_not_in_te_idx])$annotation))
colnames(dmr_nonte_context) <- c("Region", "Count")
dmr_nonte_context$Percentage <- round(100 * dmr_nonte_context$Count / sum(dmr_nonte_context$Count), 2)
dmr_nonte_context$Category <- "DMRs_not_in_TEs"

dmr_all_context <- as.data.frame(table(mcols(dmr_gr)$annotation))
colnames(dmr_all_context) <- c("Region", "Count")
dmr_all_context$Percentage <- round(100 * dmr_all_context$Count / sum(dmr_all_context$Count), 2)
dmr_all_context$Category <- "All_DMRs"

dmr_context_combined <- rbind(dmr_te_context, dmr_nonte_context, dmr_all_context)

cat("\n  DMRs IN TEs:\n")
for (i in 1:nrow(dmr_te_context)) {
  cat(sprintf("    %-15s %5d (%5.1f%%)\n",
              dmr_te_context$Region[i], dmr_te_context$Count[i], dmr_te_context$Percentage[i]))
}
cat("\n  DMRs NOT in TEs:\n")
for (i in 1:nrow(dmr_nonte_context)) {
  cat(sprintf("    %-15s %5d (%5.1f%%)\n",
              dmr_nonte_context$Region[i], dmr_nonte_context$Count[i], dmr_nonte_context$Percentage[i]))
}

fwrite(dmr_context_combined, file.path(OUTPUT_DIR, "tables", "dmr_te_genomic_context.tsv"), sep = "\t")

# ── TE class breakdown within genic vs intergenic ──
cat("\nTE class breakdown for DMP-overlapping TEs by genomic context...\n")

# Build a detailed table: for each DMP-TE overlap, record both annotation and TE class
dmp_te_detail <- data.table(
  dmp_idx = queryHits(dmp_te_hits),
  te_idx = subjectHits(dmp_te_hits)
)
dmp_te_detail[, annotation := mcols(dmp_gr)$annotation[dmp_idx]]
dmp_te_detail[, te_class := mcols(te_gr)$te_class[te_idx]]
dmp_te_detail[, te_family := mcols(te_gr)$te_family[te_idx]]
dmp_te_detail[, meth_diff := mcols(dmp_gr)$methylation_diff[dmp_idx]]
dmp_te_detail[, direction := ifelse(meth_diff > 0, "Hyper", "Hypo")]

# Cross-tabulation: annotation × TE class (unique DMPs)
te_class_by_region <- dmp_te_detail[, .N, by = .(annotation, te_class)]
te_class_by_region[, pct := round(100 * N / sum(N), 2), by = annotation]

cat("\n  TE class × Genomic region cross-table:\n")
te_cross <- dcast(te_class_by_region, te_class ~ annotation, value.var = "N", fill = 0)
print(te_cross)

fwrite(te_class_by_region, file.path(OUTPUT_DIR, "tables", "dmp_te_class_by_genomic_region.tsv"), sep = "\t")
fwrite(dmp_te_detail, file.path(OUTPUT_DIR, "tables", "dmp_te_overlap_detail.tsv"), sep = "\t")

# ── Fisher's exact test: are TE-DMPs enriched in intergenic vs gene body? ──
cat("\nFisher's exact test: TE-DMPs enriched in intergenic?\n")

n_te_intergenic <- sum(mcols(dmp_gr[dmp_in_te_idx])$annotation == "Intergenic")
n_te_genic <- sum(mcols(dmp_gr[dmp_in_te_idx])$annotation != "Intergenic")
n_nonte_intergenic <- sum(mcols(dmp_gr[dmp_not_in_te_idx])$annotation == "Intergenic")
n_nonte_genic <- sum(mcols(dmp_gr[dmp_not_in_te_idx])$annotation != "Intergenic")

fisher_mat <- matrix(c(n_te_intergenic, n_te_genic,
                        n_nonte_intergenic, n_nonte_genic),
                     nrow = 2, byrow = TRUE,
                     dimnames = list(c("In_TE", "Not_in_TE"),
                                     c("Intergenic", "Genic")))
fisher_res <- fisher.test(fisher_mat)

cat("  Contingency table:\n")
print(fisher_mat)
cat(sprintf("  OR = %.3f | p = %.2e\n", fisher_res$estimate, fisher_res$p.value))
cat(sprintf("  Interpretation: DMPs in TEs are %s for intergenic regions\n",
            ifelse(fisher_res$estimate > 1, "ENRICHED", "DEPLETED")))

# ── Visualization: DMP context in TE vs non-TE ──
p1 <- ggplot(dmp_context_combined, aes(x = Category, y = Percentage, fill = Region)) +
  geom_col(position = "stack", width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", Percentage)),
            position = position_stack(vjust = 0.5), size = 3.5) +
  scale_fill_manual(values = c("Gene Body" = "#2E86C1",
                                "Intergenic" = "#E74C3C",
                                "Promoter" = "#27AE60")) +
  labs(title = "DMP Genomic Context: TEs vs Non-TEs",
       subtitle = sprintf("DMPs in TEs: %d (%.1f%%) | Fisher OR = %.2f, p = %.2e",
                           length(dmp_in_te_idx),
                           100 * length(dmp_in_te_idx) / length(dmp_gr),
                           fisher_res$estimate, fisher_res$p.value),
       x = "", y = "Percentage (%)") +
  coord_flip()
save_both(p1, "01_dmp_genomic_context_te_vs_nonte", 10, 5)

# TE class stacked bar by region
p2 <- ggplot(te_class_by_region, aes(x = annotation, y = N, fill = te_class)) +
  geom_col(position = "stack") +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "TE Classes in DMP-Overlapping TEs by Genomic Region",
       x = "Genomic Annotation", y = "Number of DMP-TE Overlaps",
       fill = "TE Class") +
  coord_flip()
save_both(p2, "02_te_class_by_genomic_region", 10, 6)


# =============================================================================
# PART 2: GENOME METHYLATION LANDSCAPE (genomation)
# =============================================================================

cat("\n================================================================\n")
cat(" PART 2: Genome Methylation Landscape (genomation)              \n")
cat("================================================================\n\n")

# Check if genomation is available
has_genomation <- requireNamespace("genomation", quietly = TRUE)
if (!has_genomation) {
  cat("Installing genomation...\n")
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  BiocManager::install("genomation", update = FALSE, ask = FALSE)
  has_genomation <- requireNamespace("genomation", quietly = TRUE)
}

if (has_genomation) {
  library(genomation)

  cat("Loading GFF for gene coordinates...\n")
  gff <- import(GFF_FILE)
  genes <- gff[gff$type == "gene"]
  genes <- genes[as.character(seqnames(genes)) %in% keep_chr]
  cat("  Total genes (chr1-31):", length(genes), "\n")

  # Load CpG reports — construct methylation bedGraph
  # Bismark CpG_report format: chr, pos, strand, count_meth, count_unmeth, context, trinuc
  cat("\nLoading CpG methylation data from Bismark reports...\n")
  cat("  (This may take several minutes for large genomes)\n\n")

  cpg_report_files <- list(
    C1 = file.path(DATA_DIR, "C1.CpG_report.txt"),
    C2 = file.path(DATA_DIR, "C2.CpG_report.txt"),
    A1 = file.path(DATA_DIR, "A1.CpG_report.txt"),
    A2 = file.path(DATA_DIR, "A2.CpG_report.txt")
  )

  # Check which files exist
  available_files <- cpg_report_files[sapply(cpg_report_files, file.exists)]
  cat("  Available CpG reports:", paste(names(available_files), collapse = ", "), "\n")

  if (length(available_files) >= 2) {

    # Read and average methylation across samples
    # Read one file at a time to manage memory
    cat("  Reading CpG reports and computing average methylation...\n")

    all_meth <- NULL
    for (sample_name in names(available_files)) {
      cat("    Reading", sample_name, "... ")
      cpg <- fread(available_files[[sample_name]], header = FALSE,
                   col.names = c("chr", "pos", "strand", "count_meth",
                                 "count_unmeth", "context", "trinuc"),
                   select = c(1, 2, 4, 5))  # Only read needed columns

      # Filter to chr1-31 and CpGs with coverage
      cpg <- cpg[chr %in% keep_chr]
      cpg <- cpg[count_meth + count_unmeth >= 5]
      cpg[, beta := count_meth / (count_meth + count_unmeth)]

      if (is.null(all_meth)) {
        all_meth <- cpg[, .(chr, pos, beta)]
        setnames(all_meth, "beta", sample_name)
      } else {
        all_meth <- merge(all_meth, cpg[, .(chr, pos, beta)],
                          by = c("chr", "pos"), all = FALSE)
        setnames(all_meth, "beta", sample_name)
      }
      cat(nrow(cpg), "CpGs with coverage >= 5\n")
      rm(cpg); gc(verbose = FALSE)
    }

    # Compute mean methylation per condition
    control_cols <- intersect(c("C1", "C2"), names(all_meth))
    amputated_cols <- intersect(c("A1", "A2"), names(all_meth))

    all_meth[, mean_control := rowMeans(.SD, na.rm = TRUE), .SDcols = control_cols]
    all_meth[, mean_amputated := rowMeans(.SD, na.rm = TRUE), .SDcols = amputated_cols]
    all_meth[, mean_all := rowMeans(.SD, na.rm = TRUE), .SDcols = c(control_cols, amputated_cols)]

    cat("\n  CpGs with data in all samples:", nrow(all_meth), "\n")

    # Build GRanges for genomation
    meth_gr <- GRanges(seqnames = all_meth$chr,
                       ranges = IRanges(start = all_meth$pos, width = 1),
                       score = all_meth$mean_all,
                       control = all_meth$mean_control,
                       amputated = all_meth$mean_amputated)

    # ── Meta-gene profile ──
    cat("\n  Computing meta-gene methylation profile...\n")

    # Resize genes to include flanking regions
    # genomation expects consistent strand info
    gene_flanked <- genes
    strand(gene_flanked)[strand(gene_flanked) == "*"] <- "+"

    # ScoreMatrixBin: bin gene bodies + flanks into equal windows
    # 50 bins for upstream (2kb), 100 bins for gene body, 50 bins for downstream (2kb)
    sm_all <- ScoreMatrixBin(target = meth_gr,
                              windows = gene_flanked,
                              bin.num = 100,
                              strand.aware = TRUE)

    sm_ctrl <- ScoreMatrixBin(target = meth_gr[, "control"],
                               windows = gene_flanked,
                               bin.num = 100,
                               strand.aware = TRUE)

    sm_amp <- ScoreMatrixBin(target = meth_gr[, "amputated"],
                              windows = gene_flanked,
                              bin.num = 100,
                              strand.aware = TRUE)

    # Calculate mean across all genes per bin
    profile_all <- colMeans(sm_all, na.rm = TRUE)
    profile_ctrl <- colMeans(sm_ctrl, na.rm = TRUE)
    profile_amp <- colMeans(sm_amp, na.rm = TRUE)

    profile_df <- data.frame(
      bin = rep(1:100, 3),
      methylation = c(profile_all, profile_ctrl, profile_amp) * 100,
      condition = rep(c("All samples", "Control", "Amputated"), each = 100)
    )

    p3 <- ggplot(profile_df, aes(x = bin, y = methylation, color = condition)) +
      geom_line(linewidth = 1.2) +
      scale_color_manual(values = c("All samples" = "grey40",
                                     "Control" = "#2E86C1",
                                     "Amputated" = "#E74C3C")) +
      geom_vline(xintercept = c(1, 100), linetype = "dashed", color = "grey60") +
      annotate("text", x = 50, y = max(profile_df$methylation, na.rm = TRUE),
               label = "Gene Body", size = 4) +
      labs(title = "Meta-Gene Methylation Profile (D. laeve)",
           subtitle = sprintf("Averaged over %s genes | CpG coverage >= 5x",
                              format(length(genes), big.mark = ",")),
           x = "Gene Body Position (TSS -> TES, 100 bins)",
           y = "Mean CpG Methylation (%)",
           color = "Condition")
    save_both(p3, "03_metagene_methylation_landscape", 12, 6)

    # ── Heatmap of methylation across genes ──
    cat("  Generating methylation heatmap...\n")

    png(file.path(OUTPUT_DIR, "plots", "04_metagene_methylation_heatmap.png"),
        width = 1000, height = 800, res = 150)
    heatMatrix(sm_all,
               xcoords = c(1, 100),
               xlab = "Gene Body Position (TSS -> TES)",
               main = "CpG Methylation Across All Genes")
    dev.off()

    pdf(file.path(OUTPUT_DIR, "plots", "04_metagene_methylation_heatmap.pdf"),
        width = 10, height = 8)
    heatMatrix(sm_all,
               xcoords = c(1, 100),
               xlab = "Gene Body Position (TSS -> TES)",
               main = "CpG Methylation Across All Genes")
    dev.off()
    cat("  Saved: 04_metagene_methylation_heatmap\n")

    # ── Methylation landscape by gene size quartile ──
    gene_sizes <- width(genes)
    size_q <- quantile(gene_sizes, c(0.25, 0.5, 0.75))
    gene_size_class <- cut(gene_sizes,
                           breaks = c(0, size_q[1], size_q[2], size_q[3], Inf),
                           labels = c("Q1_small", "Q2_medium", "Q3_large", "Q4_xlarge"))

    sm_list <- list()
    for (q in levels(gene_size_class)) {
      idx <- which(gene_size_class == q)
      if (length(idx) > 50) {
        sm_q <- ScoreMatrixBin(target = meth_gr,
                                windows = gene_flanked[idx],
                                bin.num = 100,
                                strand.aware = TRUE)
        sm_list[[q]] <- colMeans(sm_q, na.rm = TRUE) * 100
      }
    }

    if (length(sm_list) > 0) {
      size_profile <- do.call(rbind, lapply(names(sm_list), function(q) {
        data.frame(bin = 1:100, methylation = sm_list[[q]], size_class = q)
      }))

      p4 <- ggplot(size_profile, aes(x = bin, y = methylation, color = size_class)) +
        geom_line(linewidth = 1) +
        scale_color_brewer(palette = "Set2") +
        labs(title = "Meta-Gene Methylation by Gene Size",
             x = "Gene Body Position (TSS -> TES, 100 bins)",
             y = "Mean CpG Methylation (%)",
             color = "Gene Size Quartile")
      save_both(p4, "05_metagene_by_gene_size", 12, 6)
    }

    # Save profile data
    fwrite(profile_df, file.path(OUTPUT_DIR, "tables", "metagene_methylation_profile.tsv"), sep = "\t")

    rm(all_meth, meth_gr, sm_all, sm_ctrl, sm_amp)
    gc(verbose = FALSE)

  } else {
    cat("  WARNING: CpG report files not found locally.\n")
    cat("  Run on cluster with access to Bismark CpG reports.\n")
  }
} else {
  cat("  WARNING: genomation not available. Install via BiocManager::install('genomation')\n")
}


# =============================================================================
# PART 3: DEEPTFACTOR TF METHYLATION ANALYSIS
# =============================================================================

cat("\n================================================================\n")
cat(" PART 3: TF Methylation (DeepTFactor)                          \n")
cat("================================================================\n\n")

cat("Loading DeepTFactor predictions...\n")
tf_pred <- fread(TF_FILE)
cat("  Total proteins:", nrow(tf_pred), "\n")
cat("  Predicted TFs:", sum(tf_pred$prediction == TRUE), "\n")
cat("  Non-TFs:", sum(tf_pred$prediction == FALSE), "\n\n")

# Extract gene ID from sequence_ID (e.g., LOC_00032486-mRNA-1 -> LOC_00032486)
tf_pred[, gene_id := sub("-mRNA-.*$", "", sequence_ID)]

# Keep unique genes (some have multiple mRNA isoforms)
# A gene is TF if ANY isoform is predicted TF
tf_genes <- tf_pred[, .(is_tf = any(prediction == TRUE),
                         max_score = max(score),
                         n_isoforms = .N), by = gene_id]
cat("  Unique genes:", nrow(tf_genes), "\n")
cat("  Unique TF genes:", sum(tf_genes$is_tf), "\n\n")

# ── Cross-reference TFs with DMPs ──
cat("Cross-referencing TFs with methylation data...\n")

# DMPs have nearest_gene column
dmp_genes <- unique(dmps$nearest_gene)
dmr_genes <- unique(dmrs$nearest_gene)

# Which TFs have DMPs?
tf_genes[, has_dmp := gene_id %in% dmp_genes]
tf_genes[, has_dmr := gene_id %in% dmr_genes]

n_tf_with_dmp <- sum(tf_genes$is_tf & tf_genes$has_dmp)
n_tf_total <- sum(tf_genes$is_tf)
n_nontf_with_dmp <- sum(!tf_genes$is_tf & tf_genes$has_dmp)
n_nontf_total <- sum(!tf_genes$is_tf)

cat(sprintf("  TFs with DMPs: %d / %d (%.1f%%)\n",
            n_tf_with_dmp, n_tf_total, 100 * n_tf_with_dmp / n_tf_total))
cat(sprintf("  Non-TFs with DMPs: %d / %d (%.1f%%)\n",
            n_nontf_with_dmp, n_nontf_total, 100 * n_nontf_with_dmp / n_nontf_total))

n_tf_with_dmr <- sum(tf_genes$is_tf & tf_genes$has_dmr)
n_nontf_with_dmr <- sum(!tf_genes$is_tf & tf_genes$has_dmr)
cat(sprintf("  TFs with DMRs: %d / %d (%.1f%%)\n",
            n_tf_with_dmr, n_tf_total, 100 * n_tf_with_dmr / n_tf_total))
cat(sprintf("  Non-TFs with DMRs: %d / %d (%.1f%%)\n",
            n_nontf_with_dmr, n_nontf_total, 100 * n_nontf_with_dmr / n_nontf_total))

# Fisher's test: TFs depleted for DMPs?
fisher_tf_mat <- matrix(c(n_tf_with_dmp, n_tf_total - n_tf_with_dmp,
                           n_nontf_with_dmp, n_nontf_total - n_nontf_with_dmp),
                        nrow = 2, byrow = TRUE,
                        dimnames = list(c("TF", "Non-TF"),
                                        c("Has_DMP", "No_DMP")))
fisher_tf <- fisher.test(fisher_tf_mat)

cat(sprintf("\n  Fisher's exact test (TF vs DMP):\n"))
cat(sprintf("    OR = %.3f | p = %.2e\n", fisher_tf$estimate, fisher_tf$p.value))
cat(sprintf("    TFs are %s for DMPs\n",
            ifelse(fisher_tf$estimate < 1, "DEPLETED", "ENRICHED")))

# Same for DMRs
fisher_dmr_mat <- matrix(c(n_tf_with_dmr, n_tf_total - n_tf_with_dmr,
                            n_nontf_with_dmr, n_nontf_total - n_nontf_with_dmr),
                         nrow = 2, byrow = TRUE,
                         dimnames = list(c("TF", "Non-TF"),
                                         c("Has_DMR", "No_DMR")))
fisher_tf_dmr <- fisher.test(fisher_dmr_mat)

cat(sprintf("  Fisher's exact test (TF vs DMR):\n"))
cat(sprintf("    OR = %.3f | p = %.2e\n", fisher_tf_dmr$estimate, fisher_tf_dmr$p.value))

# ── WHERE are TFs methylated? ──
cat("\nGenomic location of TF methylation:\n")

# Get DMPs that map to TF genes
tf_gene_ids <- tf_genes[is_tf == TRUE]$gene_id
dmps_in_tfs <- dmps[nearest_gene %in% tf_gene_ids]

if (nrow(dmps_in_tfs) > 0) {
  tf_dmp_regions <- as.data.frame(table(dmps_in_tfs$annotation))
  colnames(tf_dmp_regions) <- c("Region", "Count")
  tf_dmp_regions$Percentage <- round(100 * tf_dmp_regions$Count / sum(tf_dmp_regions$Count), 2)

  cat("  DMPs in TF genes by region:\n")
  for (i in 1:nrow(tf_dmp_regions)) {
    cat(sprintf("    %-15s %5d (%5.1f%%)\n",
                tf_dmp_regions$Region[i], tf_dmp_regions$Count[i], tf_dmp_regions$Percentage[i]))
  }

  # Compare with non-TF genes
  nontf_gene_ids <- tf_genes[is_tf == FALSE]$gene_id
  dmps_in_nontfs <- dmps[nearest_gene %in% nontf_gene_ids]
  nontf_dmp_regions <- as.data.frame(table(dmps_in_nontfs$annotation))
  colnames(nontf_dmp_regions) <- c("Region", "Count")
  nontf_dmp_regions$Percentage <- round(100 * nontf_dmp_regions$Count / sum(nontf_dmp_regions$Count), 2)

  cat("\n  DMPs in Non-TF genes by region:\n")
  for (i in 1:nrow(nontf_dmp_regions)) {
    cat(sprintf("    %-15s %5d (%5.1f%%)\n",
                nontf_dmp_regions$Region[i], nontf_dmp_regions$Count[i], nontf_dmp_regions$Percentage[i]))
  }

  # Direction of methylation in TFs
  tf_direction <- as.data.frame(table(ifelse(dmps_in_tfs$methylation_diff > 0, "Hyper", "Hypo")))
  colnames(tf_direction) <- c("Direction", "Count")
  tf_direction$Percentage <- round(100 * tf_direction$Count / sum(tf_direction$Count), 2)
  cat("\n  TF DMP direction:\n")
  for (i in 1:nrow(tf_direction)) {
    cat(sprintf("    %-8s %5d (%5.1f%%)\n",
                tf_direction$Direction[i], tf_direction$Count[i], tf_direction$Percentage[i]))
  }

  fwrite(tf_dmp_regions, file.path(OUTPUT_DIR, "tables", "tf_dmp_by_region.tsv"), sep = "\t")
}

# ── DMP count per TF gene ──
dmp_count_per_tf <- dmps[nearest_gene %in% tf_gene_ids, .N, by = nearest_gene]
setnames(dmp_count_per_tf, c("gene_id", "n_dmps"))
dmp_count_per_tf <- merge(dmp_count_per_tf, tf_genes[is_tf == TRUE, .(gene_id, max_score)],
                          by = "gene_id")

# Add gene names from EviAnn
if (file.exists(ANNOT_FILE)) {
  annot <- fread(ANNOT_FILE, header = FALSE, col.names = c("gene_id", "gene_name", "description"))
  annot_unique <- annot[!duplicated(gene_id)]
  dmp_count_per_tf <- merge(dmp_count_per_tf, annot_unique[, .(gene_id, gene_name, description)],
                            by = "gene_id", all.x = TRUE)
}

dmp_count_per_tf <- dmp_count_per_tf[order(-n_dmps)]

cat(sprintf("\n  Top 20 most methylated TFs:\n"))
top_tfs <- head(dmp_count_per_tf, 20)
for (i in 1:nrow(top_tfs)) {
  gene_label <- if ("gene_name" %in% colnames(top_tfs) && !is.na(top_tfs$gene_name[i]) && top_tfs$gene_name[i] != "") {
    paste0(top_tfs$gene_id[i], " (", top_tfs$gene_name[i], ")")
  } else {
    top_tfs$gene_id[i]
  }
  cat(sprintf("    %-35s %3d DMPs  (TF score: %.4f)\n",
              gene_label, top_tfs$n_dmps[i], top_tfs$max_score[i]))
}

fwrite(dmp_count_per_tf, file.path(OUTPUT_DIR, "tables", "tf_genes_dmp_counts.tsv"), sep = "\t")

# ── Visualization: TF vs Non-TF DMP rates ──
tf_summary <- data.frame(
  Category = c("TF", "Non-TF"),
  DMP_rate = c(100 * n_tf_with_dmp / n_tf_total,
               100 * n_nontf_with_dmp / n_nontf_total),
  DMR_rate = c(100 * n_tf_with_dmr / n_tf_total,
               100 * n_nontf_with_dmr / n_nontf_total),
  Total = c(n_tf_total, n_nontf_total)
)

tf_long <- tf_summary %>%
  pivot_longer(cols = c(DMP_rate, DMR_rate),
               names_to = "type", values_to = "rate") %>%
  mutate(type = gsub("_rate", "", type))

p5 <- ggplot(tf_long, aes(x = Category, y = rate, fill = type)) +
  geom_col(position = "dodge", width = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%", rate)),
            position = position_dodge(width = 0.6), vjust = -0.3, size = 4) +
  scale_fill_manual(values = c("DMP" = "#E74C3C", "DMR" = "#3498DB")) +
  labs(title = "TF vs Non-TF Methylation Rates",
       subtitle = sprintf("Fisher DMP OR = %.3f (p = %.2e) | TFs %s",
                           fisher_tf$estimate, fisher_tf$p.value,
                           ifelse(fisher_tf$estimate < 1, "DEPLETED", "ENRICHED")),
       x = "", y = "Genes with Methylation Change (%)",
       fill = "Change Type") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
save_both(p5, "05_tf_vs_nontf_methylation_rate", 8, 6)

# TF DMP region comparison
if (nrow(dmps_in_tfs) > 0) {
  region_compare <- rbind(
    data.frame(tf_dmp_regions, Category = "TF genes"),
    data.frame(nontf_dmp_regions, Category = "Non-TF genes")
  )

  p6 <- ggplot(region_compare, aes(x = Category, y = Percentage, fill = Region)) +
    geom_col(position = "stack", width = 0.6) +
    geom_text(aes(label = sprintf("%.1f%%", Percentage)),
              position = position_stack(vjust = 0.5), size = 3.5) +
    scale_fill_manual(values = c("Gene Body" = "#2E86C1",
                                  "Intergenic" = "#E74C3C",
                                  "Promoter" = "#27AE60")) +
    labs(title = "Where Are TF Genes Methylated?",
         subtitle = "Genomic region distribution of DMPs in TF vs Non-TF genes",
         x = "", y = "Percentage of DMPs (%)") +
    coord_flip()
  save_both(p6, "06_tf_dmp_region_distribution", 10, 5)
}

# Save TF summary
fwrite(tf_genes, file.path(OUTPUT_DIR, "tables", "tf_genes_methylation_status.tsv"), sep = "\t")
fwrite(tf_summary, file.path(OUTPUT_DIR, "tables", "tf_vs_nontf_summary.tsv"), sep = "\t")


# =============================================================================
# FINAL SUMMARY
# =============================================================================

cat("\n\n")
cat("================================================================\n")
cat("  ANALYSIS COMPLETE                                             \n")
cat("================================================================\n\n")

cat("PART 1 — TE Methylation Genomic Context:\n")
cat(sprintf("  DMPs in TEs: %d (%.1f%%)\n", length(dmp_in_te_idx),
            100 * length(dmp_in_te_idx) / length(dmp_gr)))
cat(sprintf("  DMRs in TEs: %d (%.1f%%)\n", length(dmr_in_te_idx),
            100 * length(dmr_in_te_idx) / length(dmr_gr)))
cat(sprintf("  TE-DMPs intergenic: %.1f%% vs non-TE DMPs intergenic: %.1f%%\n",
            dmp_te_context$Percentage[dmp_te_context$Region == "Intergenic"],
            dmp_nonte_context$Percentage[dmp_nonte_context$Region == "Intergenic"]))
cat(sprintf("  Fisher OR = %.3f, p = %.2e (TEs %s intergenic)\n",
            fisher_res$estimate, fisher_res$p.value,
            ifelse(fisher_res$estimate > 1, "enriched", "depleted")))

cat("\nPART 2 — Methylation Landscape:\n")
if (has_genomation && length(available_files) >= 2) {
  cat("  Meta-gene profiles computed for control and amputated\n")
  cat("  Heatmap generated across all annotated genes\n")
} else {
  cat("  Skipped (needs genomation + CpG reports)\n")
}

cat("\nPART 3 — TF Methylation:\n")
cat(sprintf("  Predicted TFs: %d genes\n", n_tf_total))
cat(sprintf("  TFs with DMPs: %d (%.1f%%)\n", n_tf_with_dmp, 100 * n_tf_with_dmp / n_tf_total))
cat(sprintf("  TFs with DMRs: %d (%.1f%%)\n", n_tf_with_dmr, 100 * n_tf_with_dmr / n_tf_total))
cat(sprintf("  Fisher OR (DMP): %.3f, p = %.2e — TFs %s\n",
            fisher_tf$estimate, fisher_tf$p.value,
            ifelse(fisher_tf$estimate < 1, "DEPLETED", "ENRICHED")))
cat(sprintf("  Most methylated TF: %s (%d DMPs)\n",
            dmp_count_per_tf$gene_id[1], dmp_count_per_tf$n_dmps[1]))

cat(sprintf("\nOutput: %s/ (%d tables, %d plots)\n",
            OUTPUT_DIR,
            length(list.files(file.path(OUTPUT_DIR, "tables"), pattern = "\\.tsv$")),
            length(list.files(file.path(OUTPUT_DIR, "plots"), pattern = "\\.(png|pdf)$"))))

cat("\n================================================================\n")
