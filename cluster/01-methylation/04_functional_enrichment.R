#!/usr/bin/env Rscript
# =============================================================================
# Part 4: Gene Names, ChIPseeker, GO/KEGG/Reactome Enrichment,
#          Folder Organization & PNG Generation
# =============================================================================
# Prerequisites: Run Parts 1, 2, 3 first
# Uses: protein.enrichment.terms.v12.0.txt from STRING database
# =============================================================================

suppressPackageStartupMessages({
  library(bsseq)
  library(DSS)
  library(GenomicRanges)
  library(GenomicFeatures)
  library(rtracklayer)
  library(data.table)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(gridExtra)
})

for (pkg in c("ChIPseeker", "clusterProfiler")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  }
}
suppressPackageStartupMessages({
  library(ChIPseeker)
  library(clusterProfiler)
})

options(stringsAsFactors = FALSE)
options(scipen = 999)

theme_set(theme_minimal() +
            theme(text = element_text(size = 12),
                  plot.title = element_text(size = 14, face = "bold"),
                  axis.title = element_text(size = 12),
                  legend.text = element_text(size = 10)))

# ── Paths ────────────────────────────────────────────────────────────────────

BASE_DIR    <- "/mnt/data/alfredvar/"
RESULTS_DIR <- "/mnt/data/alfredvar/rlopezt/feb_w2_scripts/Preliminary/"

GFF_FILE <- file.path(BASE_DIR,
                      "30-Genoma/31-Alternative_Annotation_EviAnn",
                      "derLaeGenome_namesDlasi_v2.fasta.functional_note.pseudo_label.gff")

STRING_FILE <- file.path(RESULTS_DIR, "protein.enrichment.terms.v12.0.txt")
keep_chr <- paste0("chr", 1:31)

# Analysis parameters (same as Parts 1-3)
analysis_params <- list(
  min_coverage   = 5,
  min_diff       = 0.1,
  fdr_threshold  = 0.05,
  min_dmr_cpgs   = 3,
  min_dmr_length = 50
)

sample_info <- data.frame(
  sample_id = c("C1", "C2", "A1", "A2"),
  condition = c("Control", "Control", "Amputated", "Amputated"),
  stringsAsFactors = FALSE
)
sample_names <- sample_info$sample_id

# Dynamic condition color helper - reads actual levels from bs_obj
# so pheatmap annotation_colors never mismatch factor levels
get_condition_colors <- function(bs) {
  conds <- unique(as.character(pData(bs)$condition))
  cols <- setNames(rep("#E74C3C", length(conds)), conds)
  if ("Control" %in% conds) cols["Control"] <- "#2E86C1"
  return(cols)
}

# ── Create Organized Folder Structure ────────────────────────────────────────

cat("=== Creating organized folder structure ===\n\n")

dirs <- list(
  pdf        = file.path(RESULTS_DIR, "PDFs"),
  png        = file.path(RESULTS_DIR, "PNGs"),
  tables     = file.path(RESULTS_DIR, "Tables"),
  rds        = file.path(RESULTS_DIR, "RDS_Objects"),
  enrichment = file.path(RESULTS_DIR, "Enrichment"),
  chipseeker = file.path(RESULTS_DIR, "ChIPseeker")
)

for (d in dirs) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  cat(sprintf("  Created: %s\n", d))
}

# ── Helper Functions ─────────────────────────────────────────────────────────

# Save ggplot as both PDF and PNG
save_both <- function(plot_obj, name, width, height,
                      pdf_dir = dirs$pdf, png_dir = dirs$png) {
  ggsave(file.path(pdf_dir, paste0(name, ".pdf")), plot_obj,
         width = width, height = height)
  ggsave(file.path(png_dir, paste0(name, ".png")), plot_obj,
         width = width, height = height, dpi = 300)
  cat(sprintf("  Saved: %s (.pdf + .png)\n", name))
}

# Save pheatmap as both PDF and PNG
save_pheatmap_both <- function(pheatmap_expr, name, width, height,
                               pdf_dir = dirs$pdf, png_dir = dirs$png) {
  pdf(file.path(pdf_dir, paste0(name, ".pdf")), width = width, height = height)
  eval(pheatmap_expr)
  dev.off()
  png(file.path(png_dir, paste0(name, ".png")),
      width = width, height = height, units = "in", res = 300)
  eval(pheatmap_expr)
  dev.off()
  cat(sprintf("  Saved: %s (.pdf + .png)\n", name))
}


# #############################################################################
# SECTION A: ANALYSIS (Gene Names, ChIPseeker, Enrichment)
# #############################################################################

cat("\n")
cat("================================================================\n")
cat("  SECTION A: Gene Names + ChIPseeker + Enrichment Analysis\n")
cat("================================================================\n\n")

# ── A1. Load Previous Results ────────────────────────────────────────────────

cat("Loading results from Parts 1-3...\n")
bs_obj      <- readRDS(file.path(RESULTS_DIR, "bsseq_object.rds"))
dmp_results <- readRDS(file.path(RESULTS_DIR, "dmp_results.rds"))
dmr_results <- readRDS(file.path(RESULTS_DIR, "dmr_results.rds"))
dml_test    <- readRDS(file.path(RESULTS_DIR, "dml_test_results.rds"))

cat(sprintf("  BSseq: %s sites x %s samples\n",
            format(nrow(bs_obj), big.mark = ","), ncol(bs_obj)))
cat(sprintf("  DMPs: %s | DMRs: %s\n",
            format(nrow(dmp_results), big.mark = ","),
            format(nrow(dmr_results), big.mark = ",")))

# Set condition colors from actual data
condition_colors <- get_condition_colors(bs_obj)
cat(sprintf("  Condition levels: %s\n", paste(names(condition_colors), collapse = ", ")))
# ── A2. Extract Gene Names from GFF ─────────────────────────────────────────

cat("\n--- Extracting gene names from GFF ---\n")

genes_gff <- import(GFF_FILE)
genes_gff <- genes_gff[seqnames(genes_gff) %in% keep_chr]
genes <- genes_gff[genes_gff$type == "gene"]

cat("Available metadata columns in GFF:\n")
print(colnames(mcols(genes)))

# Print first 5 genes
cat("\nFirst 5 gene entries:\n")
for (i in 1:min(5, length(genes))) {
  cat(sprintf("\n--- Gene %d ---\n", i))
  gene_meta <- mcols(genes)[i, ]
  for (col_name in colnames(gene_meta)) {
    raw_val <- gene_meta[[col_name]]
    # Handle CharacterList/AtomicList columns
    if (is(raw_val, "AtomicList") || is(raw_val, "CharacterList") || is(raw_val, "List")) {
      val <- if (length(raw_val[[1]]) == 0) NA_character_ else paste(raw_val[[1]], collapse = "; ")
    } else {
      val <- as.character(raw_val)
    }
    if (!is.na(val) && val != "" && val != "NA") {
      cat(sprintf("  %-20s: %s\n", col_name, substr(val, 1, 120)))
    }
  }
}

# Build gene info table
gene_info <- data.frame(
  gene_id = mcols(genes)$ID,
  chr     = as.character(seqnames(genes)),
  start   = start(genes),
  end     = end(genes),
  strand  = as.character(strand(genes)),
  stringsAsFactors = FALSE
)

possible_name_cols <- c("Name", "functional_note", "product", "Note",
                        "description", "Dbxref", "gene_name")

# Helper: safely convert GFF metadata columns to character
# Some columns are CharacterList (AtomicList) when genes have multiple values
safe_as_character <- function(x) {
  if (is(x, "AtomicList") || is(x, "CharacterList") || is(x, "List")) {
    # Collapse multi-value entries with semicolons
    return(vapply(x, function(el) {
      if (length(el) == 0 || all(is.na(el))) return(NA_character_)
      paste(el, collapse = "; ")
    }, character(1)))
  }
  return(as.character(x))
}

for (col in possible_name_cols) {
  if (col %in% colnames(mcols(genes))) {
    gene_info[[col]] <- safe_as_character(mcols(genes)[[col]])
    n_non_empty <- sum(!is.na(gene_info[[col]]) & gene_info[[col]] != "")
    cat(sprintf("  Found column '%s': %d genes have values\n", col, n_non_empty))
  }
}

gene_info$display_name <- gene_info$gene_id
for (col in rev(possible_name_cols)) {
  if (col %in% colnames(gene_info)) {
    has_value <- !is.na(gene_info[[col]]) & gene_info[[col]] != ""
    gene_info$display_name[has_value] <- gene_info[[col]][has_value]
  }
}

gene_info$gene_type <- ifelse(grepl("XLOC.*lncRNA", gene_info$gene_id), "lncRNA",
                              ifelse(grepl("^LOC_", gene_info$gene_id),
                                     "protein_coding", "other"))

write.table(gene_info, file.path(RESULTS_DIR, "gene_info_with_names.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("\nGene info saved: %s genes\n", format(nrow(gene_info), big.mark = ",")))

# ── A3. Add Names to Top Genes, DMPs, DMRs ──────────────────────────────────

cat("\n--- Adding names to results ---\n")

gene_meth_summary <- read.delim(
  file.path(RESULTS_DIR, "gene_level_methylation_summary.txt"),
  stringsAsFactors = FALSE)

top_genes <- gene_meth_summary %>%
  group_by(nearest_gene) %>%
  summarise(total_dmps = sum(n_dmps), n_hyper = sum(n_hyper), n_hypo = sum(n_hypo),
            mean_diff = mean(mean_diff), min_fdr = min(min_fdr),
            gene_type = first(gene_type), .groups = "drop") %>%
  arrange(desc(total_dmps)) %>%
  head(50) %>%
  left_join(gene_info[, c("gene_id", "display_name")],
            by = c("nearest_gene" = "gene_id"))

cat("\nTop 20 genes with names:\n")
print(as.data.frame(head(top_genes[, c("nearest_gene", "display_name",
                                       "total_dmps", "mean_diff")], 20)),
      right = FALSE)

write.table(top_genes, file.path(RESULTS_DIR, "top50_genes_with_names.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# DMRs with names
dmr_anno <- read.delim(file.path(RESULTS_DIR, "dmrs_annotated.txt"),
                       stringsAsFactors = FALSE)
dmr_anno_named <- dmr_anno %>%
  left_join(gene_info[, c("gene_id", "display_name")],
            by = c("nearest_gene" = "gene_id"))
write.table(dmr_anno_named, file.path(RESULTS_DIR, "dmrs_annotated_with_names.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# DMPs with names
dmp_anno <- read.delim(file.path(RESULTS_DIR, "dmps_annotated.txt"),
                       stringsAsFactors = FALSE)
dmp_anno_named <- dmp_anno %>%
  left_join(gene_info[, c("gene_id", "display_name")],
            by = c("nearest_gene" = "gene_id"))
write.table(dmp_anno_named, file.path(RESULTS_DIR, "dmps_annotated_with_names.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("Gene names added to all annotation files.\n")

# ── A4. ChIPseeker Annotation (FIXED) ───────────────────────────────────────

cat("\n--- ChIPseeker Annotation ---\n")

cat("Building TxDb from GFF...\n")
txdb <- makeTxDbFromGFF(GFF_FILE, format = "gff3")
cat("  TxDb created.\n")

# --- DMPs ---
if (nrow(dmp_results) > 0) {
  cat("\nAnnotating DMPs with ChIPseeker...\n")
  
  dmp_gr <- GRanges(
    seqnames = dmp_results$chr,
    ranges   = IRanges(start = dmp_results$pos, width = 1),
    strand   = "*",
    methylation_diff = dmp_results$diff,
    pvalue = dmp_results$pval,
    fdr    = dmp_results$fdr)
  
  dmp_chipseeker <- annotatePeak(dmp_gr, tssRegion = c(-3000, 3000),
                                 TxDb = txdb, level = "gene")
  
  dmp_cs_df <- as.data.frame(dmp_chipseeker)
  write.table(dmp_cs_df,
              file.path(dirs$chipseeker, "dmps_chipseeker_annotated.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  cat("  DMP annotation summary:\n")
  print(dmp_chipseeker)
  
  # FIX: plotAnnoPie uses base graphics - call directly inside device
  pdf(file.path(dirs$chipseeker, "dmp_chipseeker_pie.pdf"), width = 10, height = 8)
  tryCatch(plotAnnoPie(dmp_chipseeker), error = function(e) {
    cat("    plotAnnoPie base graphics failed, using fallback.\n")
  })
  dev.off()
  png(file.path(dirs$chipseeker, "dmp_chipseeker_pie.png"),
      width = 10, height = 8, units = "in", res = 300, bg = "white")
  tryCatch(plotAnnoPie(dmp_chipseeker), error = function(e) NULL)
  dev.off()
  
  # FIX: plotAnnoBar returns ggplot - need to capture and ggsave
  p_bar <- tryCatch(plotAnnoBar(dmp_chipseeker), error = function(e) NULL)
  if (!is.null(p_bar)) {
    ggsave(file.path(dirs$chipseeker, "dmp_chipseeker_bar.pdf"),
           p_bar, width = 12, height = 6)
    ggsave(file.path(dirs$chipseeker, "dmp_chipseeker_bar.png"),
           p_bar, width = 12, height = 6, dpi = 300)
  }
  
  # FIX: plotDistToTSS returns ggplot
  p_tss <- tryCatch(
    plotDistToTSS(dmp_chipseeker,
                  title = "Distribution of DMPs Relative to TSS"),
    error = function(e) NULL)
  if (!is.null(p_tss)) {
    ggsave(file.path(dirs$chipseeker, "dmp_chipseeker_distToTSS.pdf"),
           p_tss, width = 12, height = 6)
    ggsave(file.path(dirs$chipseeker, "dmp_chipseeker_distToTSS.png"),
           p_tss, width = 12, height = 6, dpi = 300)
  }
  
  # BONUS: ggplot2 fallback plots (always work)
  anno_stats <- dmp_chipseeker@annoStat
  if (!is.null(anno_stats) && nrow(anno_stats) > 0) {
    # Use the Feature column (not rownames, which are just row numbers)
    anno_df <- data.frame(Feature = anno_stats$Feature,
                          Frequency = anno_stats$Frequency,
                          stringsAsFactors = FALSE)
    # Add percentages to legend labels
    anno_df$Label <- paste0(anno_df$Feature, " (", round(anno_df$Frequency, 2), "%)")
    
    p_pie_fb <- ggplot(anno_df, aes(x = "", y = Frequency, fill = Label)) +
      geom_bar(stat = "identity", width = 1) +
      coord_polar("y") +
      scale_fill_manual(values = setNames(
        colorRampPalette(brewer.pal(12, "Set3"))(nrow(anno_df)),
        anno_df$Label)) +
      labs(title = "DMP Genomic Annotation (ChIPseeker)", fill = "Feature") +
      theme_void() +
      theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
            legend.text = element_text(size = 8))
    save_both(p_pie_fb, "dmp_chipseeker_annotation_pie",
              10, 8, pdf_dir = dirs$chipseeker, png_dir = dirs$chipseeker)
    
    p_bar_fb <- ggplot(anno_df, aes(x = reorder(Feature, Frequency),
                                    y = Frequency, fill = Feature)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Set3"))(nrow(anno_df))) +
      coord_flip() +
      labs(title = "DMP Genomic Annotation (ChIPseeker)",
           x = "", y = "Percentage (%)") +
      theme_minimal() +
      guides(fill = "none")
    save_both(p_bar_fb, "dmp_chipseeker_annotation_bar",
              12, 7, pdf_dir = dirs$chipseeker, png_dir = dirs$chipseeker)
  }
  cat("  DMP ChIPseeker plots saved.\n")
}

# --- DMRs ---
if (nrow(dmr_results) > 0) {
  cat("\nAnnotating DMRs with ChIPseeker...\n")
  
  dmr_gr <- GRanges(
    seqnames = dmr_results$chr,
    ranges   = IRanges(start = dmr_results$start, end = dmr_results$end))
  
  dmr_chipseeker <- annotatePeak(dmr_gr, tssRegion = c(-3000, 3000),
                                 TxDb = txdb, level = "gene")
  
  dmr_cs_df <- as.data.frame(dmr_chipseeker)
  write.table(dmr_cs_df,
              file.path(dirs$chipseeker, "dmrs_chipseeker_annotated.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  cat("  DMR annotation summary:\n")
  print(dmr_chipseeker)
  
  pdf(file.path(dirs$chipseeker, "dmr_chipseeker_pie.pdf"), width = 10, height = 8)
  tryCatch(plotAnnoPie(dmr_chipseeker), error = function(e) NULL)
  dev.off()
  png(file.path(dirs$chipseeker, "dmr_chipseeker_pie.png"),
      width = 10, height = 8, units = "in", res = 300, bg = "white")
  tryCatch(plotAnnoPie(dmr_chipseeker), error = function(e) NULL)
  dev.off()
  
  p_bar2 <- tryCatch(plotAnnoBar(dmr_chipseeker), error = function(e) NULL)
  if (!is.null(p_bar2)) {
    ggsave(file.path(dirs$chipseeker, "dmr_chipseeker_bar.pdf"),
           p_bar2, width = 12, height = 6)
    ggsave(file.path(dirs$chipseeker, "dmr_chipseeker_bar.png"),
           p_bar2, width = 12, height = 6, dpi = 300)
  }
  
  p_tss2 <- tryCatch(
    plotDistToTSS(dmr_chipseeker,
                  title = "Distribution of DMRs Relative to TSS"),
    error = function(e) NULL)
  if (!is.null(p_tss2)) {
    ggsave(file.path(dirs$chipseeker, "dmr_chipseeker_distToTSS.pdf"),
           p_tss2, width = 12, height = 6)
    ggsave(file.path(dirs$chipseeker, "dmr_chipseeker_distToTSS.png"),
           p_tss2, width = 12, height = 6, dpi = 300)
  }
  
  anno_stats2 <- dmr_chipseeker@annoStat
  if (!is.null(anno_stats2) && nrow(anno_stats2) > 0) {
    anno_df2 <- data.frame(Feature = anno_stats2$Feature,
                           Frequency = anno_stats2$Frequency,
                           stringsAsFactors = FALSE)
    anno_df2$Label <- paste0(anno_df2$Feature, " (", round(anno_df2$Frequency, 2), "%)")
    
    p_pie_fb2 <- ggplot(anno_df2, aes(x = "", y = Frequency, fill = Label)) +
      geom_bar(stat = "identity", width = 1) +
      coord_polar("y") +
      scale_fill_manual(values = setNames(
        colorRampPalette(brewer.pal(12, "Set3"))(nrow(anno_df2)),
        anno_df2$Label)) +
      labs(title = "DMR Genomic Annotation (ChIPseeker)", fill = "Feature") +
      theme_void() +
      theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
            legend.text = element_text(size = 8))
    save_both(p_pie_fb2, "dmr_chipseeker_annotation_pie",
              10, 8, pdf_dir = dirs$chipseeker, png_dir = dirs$chipseeker)
    
    p_bar_fb2 <- ggplot(anno_df2, aes(x = reorder(Feature, Frequency),
                                      y = Frequency, fill = Feature)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Set3"))(nrow(anno_df2))) +
      coord_flip() +
      labs(title = "DMR Genomic Annotation (ChIPseeker)",
           x = "", y = "Percentage (%)") +
      theme_minimal() +
      guides(fill = "none")
    save_both(p_bar_fb2, "dmr_chipseeker_annotation_bar",
              12, 7, pdf_dir = dirs$chipseeker, png_dir = dirs$chipseeker)
  }
  cat("  DMR ChIPseeker plots saved.\n")
}

# ── A5. GO & Reactome & KEGG Enrichment (STRING-based) ──────────────────────

cat("\n--- GO Enrichment (STRING-based) ---\n")

cat("Loading STRING protein enrichment file...\n")
string_raw <- fread(STRING_FILE, sep = "\t", header = TRUE, quote = "")
cat(sprintf("  Raw entries: %s\n", format(nrow(string_raw), big.mark = ",")))
cat("  Categories found:\n")
print(table(string_raw$category))

# Extract gene IDs
string_raw$gene_id <- sub("^.*\\.(LOC_\\d+)$", "\\1", string_raw$`#string_protein_id`)
string_raw$gene_id[!grepl("^LOC_", string_raw$gene_id)] <-
  sub("^.*\\.(XLOC_\\S+)$", "\\1",
      string_raw$`#string_protein_id`[!grepl("^LOC_", string_raw$gene_id)])

n_mapped <- length(unique(string_raw$gene_id[grepl("^LOC_|^XLOC_", string_raw$gene_id)]))
cat(sprintf("  Genes with STRING annotations: %s\n", format(n_mapped, big.mark = ",")))

# Build GO tables
go_data <- string_raw %>%
  filter(grepl("Gene Ontology", category)) %>%
  filter(grepl("^LOC_|^XLOC_", gene_id)) %>%
  select(gene_id, category, term, description)

go_data$go_category <- ifelse(grepl("Biological Process", go_data$category), "BP",
                              ifelse(grepl("Molecular Function", go_data$category), "MF",
                                     ifelse(grepl("Cellular Component", go_data$category), "CC", NA)))

cat(sprintf("  GO annotations: %s entries\n", format(nrow(go_data), big.mark = ",")))
print(table(go_data$go_category))

TERM2GENE <- go_data[, c("term", "gene_id")]
TERM2NAME <- go_data[, c("term", "description")] %>% distinct()

cat(sprintf("  Unique GO terms: %s | Unique genes with GO: %s\n",
            format(nrow(TERM2NAME), big.mark = ","),
            format(length(unique(TERM2GENE$gene_id)), big.mark = ",")))

# Reactome
reactome_data <- string_raw %>%
  filter(grepl("Reactome", category)) %>%
  filter(grepl("^LOC_|^XLOC_", gene_id)) %>%
  select(gene_id, term, description)
TERM2GENE_reactome <- reactome_data[, c("term", "gene_id")]
TERM2NAME_reactome <- reactome_data[, c("term", "description")] %>% distinct()
cat(sprintf("  Reactome: %s terms, %s genes\n",
            format(nrow(TERM2NAME_reactome), big.mark = ","),
            format(length(unique(TERM2GENE_reactome$gene_id)), big.mark = ",")))

# KEGG
kegg_data <- string_raw %>%
  filter(grepl("KEGG", category)) %>%
  filter(grepl("^LOC_|^XLOC_", gene_id)) %>%
  select(gene_id, term, description)
TERM2GENE_kegg <- kegg_data[, c("term", "gene_id")]
TERM2NAME_kegg <- kegg_data[, c("term", "description")] %>% distinct()
cat(sprintf("  KEGG: %s terms, %s genes\n",
            format(nrow(TERM2NAME_kegg), big.mark = ","),
            format(length(unique(TERM2GENE_kegg$gene_id)), big.mark = ",")))

# ── Enrichment function (saves to Enrichment/ folder, PDF + PNG) ────────────

perform_enrichment <- function(gene_list, TERM2GENE, TERM2NAME,
                               region_type, db_name, results_dir,
                               universe = NULL) {
  
  cat(sprintf("\n  %s enrichment for %s (%d genes)...\n",
              db_name, region_type, length(gene_list)))
  
  if (length(gene_list) < 10) {
    cat("    Not enough genes (need >= 10). Skipping.\n")
    return(NULL)
  }
  
  result <- enricher(
    gene          = gene_list,
    TERM2GENE     = TERM2GENE,
    TERM2NAME     = TERM2NAME,
    pvalueCutoff  = 0.05,
    pAdjustMethod = "BH",
    minGSSize     = 5,
    maxGSSize     = 500,
    universe      = universe
  )
  
  if (is.null(result) || nrow(result@result[result@result$p.adjust < 0.05, ]) == 0) {
    cat("    No significant enrichment found.\n")
    return(NULL)
  }
  
  sig_results <- result@result[result@result$p.adjust < 0.05, ]
  cat(sprintf("    Significant terms: %d\n", nrow(sig_results)))
  
  prefix <- tolower(paste0(region_type, "_", db_name))
  
  # Save table
  write.table(result@result,
              file.path(results_dir, paste0(prefix, "_enrichment.txt")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  top_n <- min(20, nrow(sig_results))
  top_terms <- sig_results %>% arrange(p.adjust) %>% head(top_n)
  h <- max(6, top_n * 0.35)
  
  # Dotplot (PDF + PNG)
  p_dot <- ggplot(top_terms,
                  aes(x = Count, y = reorder(Description, Count),
                      color = p.adjust, size = Count)) +
    geom_point() +
    scale_color_gradient(low = "#E74C3C", high = "#3498DB", name = "Adj. p-value") +
    scale_size_continuous(range = c(3, 10)) +
    labs(title = sprintf("Top %d %s Terms (%s genes)", top_n, db_name, region_type),
         x = "Gene Count", y = "") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 9))
  
  ggsave(file.path(results_dir, paste0(prefix, "_dotplot.pdf")),
         p_dot, width = 14, height = h)
  ggsave(file.path(results_dir, paste0(prefix, "_dotplot.png")),
         p_dot, width = 14, height = h, dpi = 300)
  
  # Barplot (PDF + PNG)
  p_bar <- ggplot(top_terms,
                  aes(x = reorder(Description, Count), y = Count, fill = p.adjust)) +
    geom_bar(stat = "identity") +
    scale_fill_gradient(low = "#E74C3C", high = "#3498DB", name = "Adj. p-value") +
    coord_flip() +
    labs(title = sprintf("%s Enrichment (%s)", db_name, region_type),
         x = "", y = "Gene Count") +
    theme_minimal()
  
  ggsave(file.path(results_dir, paste0(prefix, "_barplot.pdf")),
         p_bar, width = 14, height = h)
  ggsave(file.path(results_dir, paste0(prefix, "_barplot.png")),
         p_bar, width = 14, height = h, dpi = 300)
  
  cat(sprintf("    Plots saved: %s_dotplot + %s_barplot (.pdf + .png)\n", prefix, prefix))
  return(result)
}

# ── Run all enrichments ─────────────────────────────────────────────────────

cat("\n--- Running enrichment analyses ---\n")

dmp_genes <- unique(dmp_anno$nearest_gene[!is.na(dmp_anno$nearest_gene)])
dmr_genes <- unique(dmr_anno$nearest_gene[!is.na(dmr_anno$nearest_gene)])
all_genes <- gene_info$gene_id

cat(sprintf("  DMP-associated genes: %d\n", length(dmp_genes)))
cat(sprintf("  DMR-associated genes: %d\n", length(dmr_genes)))
cat(sprintf("  Universe (all genes): %d\n", length(all_genes)))

# GO by category (BP, MF, CC)
for (go_cat in c("BP", "MF", "CC")) {
  cat(sprintf("\n====== GO %s ======\n", go_cat))
  terms_in_cat <- go_data$term[go_data$go_category == go_cat]
  T2G_cat <- TERM2GENE[TERM2GENE$term %in% terms_in_cat, ]
  T2N_cat <- TERM2NAME[TERM2NAME$term %in% terms_in_cat, ]
  
  perform_enrichment(dmp_genes, T2G_cat, T2N_cat,
                     paste0("DMP_GO_", go_cat), paste0("GO_", go_cat),
                     dirs$enrichment, universe = all_genes)
  perform_enrichment(dmr_genes, T2G_cat, T2N_cat,
                     paste0("DMR_GO_", go_cat), paste0("GO_", go_cat),
                     dirs$enrichment, universe = all_genes)
}

# GO all combined
cat("\n====== GO All Categories Combined ======\n")
perform_enrichment(dmp_genes, TERM2GENE, TERM2NAME,
                   "DMP", "GO_all", dirs$enrichment, universe = all_genes)
perform_enrichment(dmr_genes, TERM2GENE, TERM2NAME,
                   "DMR", "GO_all", dirs$enrichment, universe = all_genes)

# Reactome
if (nrow(TERM2GENE_reactome) > 0) {
  cat("\n====== Reactome Pathways ======\n")
  perform_enrichment(dmp_genes, TERM2GENE_reactome, TERM2NAME_reactome,
                     "DMP", "Reactome", dirs$enrichment, universe = all_genes)
  perform_enrichment(dmr_genes, TERM2GENE_reactome, TERM2NAME_reactome,
                     "DMR", "Reactome", dirs$enrichment, universe = all_genes)
}

# KEGG
if (nrow(TERM2GENE_kegg) > 0) {
  cat("\n====== KEGG Pathways ======\n")
  perform_enrichment(dmp_genes, TERM2GENE_kegg, TERM2NAME_kegg,
                     "DMP", "KEGG", dirs$enrichment, universe = all_genes)
  perform_enrichment(dmr_genes, TERM2GENE_kegg, TERM2NAME_kegg,
                     "DMR", "KEGG", dirs$enrichment, universe = all_genes)
}

# Hyper vs Hypo
cat("\n====== Hyper vs Hypo Enrichment (GO BP) ======\n")

dmp_hyper_genes <- unique(dmp_anno$nearest_gene[
  !is.na(dmp_anno$nearest_gene) & dmp_anno$methylation_diff > 0])
dmp_hypo_genes <- unique(dmp_anno$nearest_gene[
  !is.na(dmp_anno$nearest_gene) & dmp_anno$methylation_diff < 0])

cat(sprintf("  Hypermethylated DMP genes: %d\n", length(dmp_hyper_genes)))
cat(sprintf("  Hypomethylated DMP genes: %d\n", length(dmp_hypo_genes)))

terms_bp <- go_data$term[go_data$go_category == "BP"]
T2G_bp <- TERM2GENE[TERM2GENE$term %in% terms_bp, ]
T2N_bp <- TERM2NAME[TERM2NAME$term %in% terms_bp, ]

perform_enrichment(dmp_hyper_genes, T2G_bp, T2N_bp,
                   "DMP_hyper", "GO_BP", dirs$enrichment, universe = all_genes)
perform_enrichment(dmp_hypo_genes, T2G_bp, T2N_bp,
                   "DMP_hypo", "GO_BP", dirs$enrichment, universe = all_genes)


# #############################################################################
# SECTION B: RE-GENERATE ALL PLOTS FROM PARTS 1-3 AS PNGs
# #############################################################################

cat("\n")
cat("================================================================\n")
cat("  SECTION B: Regenerating All Plots (PDF + PNG)\n")
cat("================================================================\n\n")

# ── B1. Sample Correlation Heatmap ──────────────────────────────────────────

cat("-- Correlation Heatmap --\n")
beta_matrix    <- getMeth(bs_obj, type = "raw")
complete_sites <- complete.cases(beta_matrix)
beta_complete  <- beta_matrix[complete_sites, ]
sample_cors    <- cor(beta_complete)

annotation_df <- data.frame(Condition = pData(bs_obj)$condition,
                            row.names = colnames(sample_cors))

heatmap_expr <- quote(
  pheatmap(sample_cors,
           annotation_col = annotation_df, annotation_row = annotation_df,
           color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
           main = "Sample-to-Sample Methylation Correlations",
           fontsize = 12, display_numbers = TRUE, number_format = "%.3f")
)
save_pheatmap_both(heatmap_expr, "sample_correlation_heatmap", 8, 6)

# ── B2. PCA ─────────────────────────────────────────────────────────────────

cat("\n-- PCA --\n")
pca_result    <- prcomp(t(beta_complete), center = TRUE, scale. = FALSE)
var_explained <- summary(pca_result)$importance[2, 1:min(4, ncol(beta_complete))] * 100

pca_data <- data.frame(
  PC1 = pca_result$x[, 1], PC2 = pca_result$x[, 2],
  Sample = rownames(pca_result$x), Condition = pData(bs_obj)$condition)

p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Condition, label = Sample)) +
  geom_point(size = 4) +
  geom_text(vjust = -1, hjust = 0.5, size = 3) +
  scale_color_manual(values = condition_colors) +
  labs(title = "PCA of DNA Methylation Profiles",
       x = sprintf("PC1 (%.1f%%)", var_explained[1]),
       y = sprintf("PC2 (%.1f%%)", var_explained[2])) +
  theme_minimal()
save_both(p_pca, "pca_analysis", 10, 8)

# ── B3. Methylation Density ─────────────────────────────────────────────────

cat("\n-- Methylation Density --\n")
beta_long <- as.data.frame(beta_matrix)
colnames(beta_long) <- sample_names
beta_long <- tidyr::pivot_longer(beta_long, cols = everything(),
                                 names_to = "Sample", values_to = "Methylation")
beta_long <- beta_long[!is.na(beta_long$Methylation), ]
beta_long$Condition <- sample_info$condition[match(beta_long$Sample, sample_info$sample_id)]

p_density <- ggplot(beta_long, aes(x = Methylation, color = Sample, linetype = Condition)) +
  geom_density(linewidth = 0.8) +
  scale_color_manual(values = c("C1" = "#2E86C1", "C2" = "#5DADE2",
                                "A1" = "#E74C3C", "A2" = "#F1948A")) +
  labs(title = "CpG Methylation Distribution by Sample",
       x = "Methylation Level (beta)", y = "Density") +
  theme_minimal()
save_both(p_density, "methylation_density", 10, 6)

# ── B4. Volcano Plot ────────────────────────────────────────────────────────

cat("\n-- Volcano Plot --\n")
if (nrow(dmp_results) > 0) {
  volcano_data <- data.frame(
    diff = dmp_results$diff, pvalue = dmp_results$pval, fdr = dmp_results$fdr,
    significant = dmp_results$fdr < analysis_params$fdr_threshold &
      abs(dmp_results$diff) > analysis_params$min_diff)
  volcano_data$category <- "Not Significant"
  volcano_data$category[volcano_data$significant & volcano_data$diff > 0] <- "Hypermethylated"
  volcano_data$category[volcano_data$significant & volcano_data$diff < 0] <- "Hypomethylated"
  volcano_data$category <- factor(volcano_data$category,
                                  levels = c("Not Significant", "Hypermethylated", "Hypomethylated"))
  n_hyper <- sum(volcano_data$category == "Hypermethylated")
  n_hypo  <- sum(volcano_data$category == "Hypomethylated")
  
  p_volcano <- ggplot(volcano_data, aes(x = diff, y = -log10(pvalue), color = category)) +
    geom_point(alpha = 0.6, size = 0.8) +
    scale_color_manual(values = c("Not Significant" = "gray70",
                                  "Hypermethylated" = "#E74C3C",
                                  "Hypomethylated" = "#3498DB"),
                       name = "Methylation Change") +
    labs(title = "Volcano Plot - Differentially Methylated Positions",
         subtitle = sprintf("Hypermethylated: %s | Hypomethylated: %s",
                            format(n_hyper, big.mark = ","), format(n_hypo, big.mark = ",")),
         x = "Methylation Difference (Treatment - Control)", y = "-log10(p-value)") +
    geom_vline(xintercept = c(-analysis_params$min_diff, analysis_params$min_diff),
               linetype = "dashed", alpha = 0.7, color = "gray50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.7, color = "gray50") +
    theme_minimal()
  save_both(p_volcano, "volcano_plot_dmps", 10, 8)
}

# ── B5. Manhattan Plot ──────────────────────────────────────────────────────

cat("\n-- Manhattan Plot --\n")
if (nrow(dmp_results) > 0) {
  manhattan_data <- dmp_results %>%
    mutate(chr_num = as.numeric(gsub("chr", "", chr)),
           significance = ifelse(fdr < analysis_params$fdr_threshold &
                                   abs(diff) > analysis_params$min_diff,
                                 "Significant", "Not Significant")) %>%
    filter(!is.na(chr_num)) %>% arrange(chr_num, pos)
  
  chr_info <- manhattan_data %>%
    group_by(chr, chr_num) %>%
    summarise(max_pos = max(pos), .groups = "drop") %>%
    arrange(chr_num) %>%
    mutate(cumulative_length = cumsum(c(0, head(max_pos, -1))),
           chr_center = cumulative_length + max_pos / 2)
  
  manhattan_plot_data <- manhattan_data %>%
    left_join(chr_info, by = c("chr", "chr_num")) %>%
    mutate(genome_pos = cumulative_length + pos)
  
  p_manhattan <- ggplot(manhattan_plot_data,
                        aes(x = genome_pos, y = diff, color = significance)) +
    geom_point(alpha = 0.6, size = 0.5) +
    scale_color_manual(values = c("Not Significant" = "gray70", "Significant" = "#E74C3C"),
                       name = "Significance") +
    scale_x_continuous(breaks = chr_info$chr_center,
                       labels = gsub("chr", "", chr_info$chr), expand = c(0.01, 0)) +
    labs(title = "Manhattan Plot - Methylation Differences Across the Genome",
         subtitle = "Treatment vs Control Comparison",
         x = "Chromosome", y = "Methylation Difference (Treatment - Control)") +
    theme_minimal() +
    theme(legend.position = "bottom", panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank()) +
    geom_hline(yintercept = c(-analysis_params$min_diff, analysis_params$min_diff),
               linetype = "dashed", alpha = 0.7, color = "gray50") +
    geom_vline(xintercept = chr_info$cumulative_length[-1],
               color = "gray80", linetype = "dashed", alpha = 0.3)
  save_both(p_manhattan, "manhattan_plot_methylation_differences", 16, 8)
}

# ── B6. Genome-Wide Methylation Profiles ────────────────────────────────────

cat("\n-- Genome-Wide Methylation Profiles --\n")

plot_genome_wide_methylation <- function(bs_obj, sample_name, window_size = 1000000) {
  sample_idx <- which(colnames(bs_obj) == sample_name)
  beta_vals  <- getMeth(bs_obj, type = "raw")[, sample_idx]
  methylation_data <- data.frame(
    chr = as.character(seqnames(bs_obj)), pos = start(bs_obj),
    methylation = beta_vals * 100, stringsAsFactors = FALSE)
  methylation_data <- methylation_data[!is.na(methylation_data$methylation), ]
  chr_order <- paste0("chr", 1:31)
  methylation_data <- methylation_data[methylation_data$chr %in% chr_order, ]
  methylation_data$chr <- factor(methylation_data$chr, levels = chr_order)
  
  methylation_binned <- methylation_data %>%
    mutate(bin = floor(pos / window_size) * window_size,
           chr_num = as.numeric(gsub("chr", "", as.character(chr)))) %>%
    group_by(chr, chr_num, bin) %>%
    summarise(mean_methylation = mean(methylation, na.rm = TRUE),
              n_cpgs = n(), .groups = "drop") %>%
    filter(n_cpgs >= 10)
  
  chr_lengths <- methylation_binned %>%
    group_by(chr, chr_num) %>%
    summarise(max_pos = max(bin), .groups = "drop") %>%
    arrange(chr_num) %>%
    mutate(cumulative_length = cumsum(c(0, head(max_pos, -1))),
           chr_center = cumulative_length + max_pos / 2)
  
  methylation_plot <- methylation_binned %>%
    left_join(chr_lengths, by = c("chr", "chr_num")) %>%
    mutate(genome_pos = cumulative_length + bin)
  
  ggplot(methylation_plot, aes(x = genome_pos, y = mean_methylation, color = chr)) +
    geom_point(alpha = 0.6, size = 0.5) +
    scale_color_manual(values = rep(c("#2E86C1", "#E74C3C"), 16), guide = "none") +
    scale_x_continuous(breaks = chr_lengths$chr_center,
                       labels = gsub("chr", "", chr_lengths$chr), expand = c(0.01, 0)) +
    scale_y_continuous(limits = c(0, 100), expand = c(0.02, 0)) +
    labs(title = paste("Whole-Genome DNA Methylation Profile:", sample_name),
         subtitle = paste("Window size:", window_size / 1e6, "Mb"),
         x = "Chromosome", y = "Mean Methylation (%)") +
    theme_minimal() +
    theme(legend.position = "none", panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          plot.title = element_text(size = 14, face = "bold")) +
    geom_vline(xintercept = chr_lengths$cumulative_length[-1],
               color = "gray80", linetype = "dashed", alpha = 0.7)
}

# Multi-page PDF
pdf(file.path(dirs$pdf, "genome_wide_methylation.pdf"), width = 16, height = 10)
for (s in colnames(bs_obj)) { print(plot_genome_wide_methylation(bs_obj, s)) }
dev.off()

# Individual PNGs
for (s in colnames(bs_obj)) {
  ggsave(file.path(dirs$png, paste0("genome_wide_methylation_", s, ".png")),
         plot_genome_wide_methylation(bs_obj, s), width = 16, height = 10, dpi = 300)
}
cat("  Saved: genome_wide_methylation (.pdf + per-sample .png)\n")

# ── B7. DMP Heatmap ─────────────────────────────────────────────────────────

cat("\n-- DMP Heatmap --\n")
if (nrow(dmp_results) > 0) {
  n_sites  <- min(100, nrow(dmp_results))
  top_dmps <- dmp_results[order(abs(dmp_results$diff), decreasing = TRUE), ][1:n_sites, ]
  dmp_ids     <- paste(top_dmps$chr, top_dmps$pos, sep = ":")
  bs_site_ids <- paste(as.character(seqnames(bs_obj)), start(bs_obj), sep = ":")
  matching_idx <- which(bs_site_ids %in% dmp_ids)
  
  if (length(matching_idx) >= 2) {
    beta_mat <- getMeth(bs_obj[matching_idx, ], type = "raw")
    rownames(beta_mat) <- paste0("DMP_", seq_len(nrow(beta_mat)))
    rv <- apply(beta_mat, 1, var, na.rm = TRUE)
    valid <- complete.cases(beta_mat) & rv > 0 & !is.na(rv)
    beta_mat <- beta_mat[valid, ]
    
    if (nrow(beta_mat) >= 2) {
      col_ann <- data.frame(Condition = pData(bs_obj)$condition,
                            row.names = colnames(beta_mat))
      col_colors <- list(Condition = condition_colors)
      
      dmp_hm_expr <- quote(
        pheatmap(beta_mat, scale = "row",
                 clustering_distance_rows = "euclidean",
                 clustering_distance_cols = "euclidean",
                 clustering_method = "complete",
                 annotation_col = col_ann, annotation_colors = col_colors,
                 color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
                 main = sprintf("Top %d DMPs - Site-Level Methylation", nrow(beta_mat)),
                 fontsize = 8, show_rownames = FALSE, show_colnames = TRUE)
      )
      save_pheatmap_both(dmp_hm_expr, "dmp_level_heatmap", 8, 12)
    }
  }
}

# ── B8. DMR Heatmap ─────────────────────────────────────────────────────────

cat("\n-- DMR Heatmap --\n")
if (nrow(dmr_results) > 0) {
  n_regions <- min(50, nrow(dmr_results))
  top_dmrs  <- dmr_results[order(abs(dmr_results$diff.Methy), decreasing = TRUE), ][1:n_regions, ]
  dmr_matrix <- matrix(nrow = n_regions, ncol = ncol(bs_obj))
  colnames(dmr_matrix) <- colnames(bs_obj)
  rownames(dmr_matrix) <- paste0("DMR_", 1:n_regions, "_",
                                 top_dmrs$chr, ":", top_dmrs$start, "-", top_dmrs$end)
  bs_ranges <- GRanges(seqnames = seqnames(bs_obj),
                       ranges = IRanges(start = start(bs_obj), width = 1))
  valid_dmrs <- logical(n_regions)
  for (i in 1:n_regions) {
    dmr <- top_dmrs[i, ]
    dmr_range <- GRanges(seqnames = dmr$chr,
                         ranges = IRanges(start = dmr$start, end = dmr$end))
    overlaps <- findOverlaps(bs_ranges, dmr_range)
    cpg_indices <- queryHits(overlaps)
    if (length(cpg_indices) > 0) {
      dmr_beta <- getMeth(bs_obj[cpg_indices, ], type = "raw")
      dmr_matrix[i, ] <- if (length(cpg_indices) == 1) as.numeric(dmr_beta) else colMeans(dmr_beta, na.rm = TRUE)
      valid_dmrs[i] <- TRUE
    }
  }
  dmr_matrix <- dmr_matrix[valid_dmrs, , drop = FALSE]
  
  if (nrow(dmr_matrix) >= 2) {
    col_ann <- data.frame(Condition = pData(bs_obj)$condition,
                          row.names = colnames(dmr_matrix))
    col_colors <- list(Condition = condition_colors)
    
    dmr_hm_expr <- quote(
      pheatmap(dmr_matrix, scale = "row",
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "complete",
               annotation_col = col_ann, annotation_colors = col_colors,
               color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
               main = sprintf("Top %d DMRs - Region-Level Methylation", nrow(dmr_matrix)),
               fontsize = 10, show_rownames = TRUE, show_colnames = TRUE)
    )
    save_pheatmap_both(dmr_hm_expr, "dmr_level_heatmap", 10, 12)
  }
}

# ── B9. DMP Annotation Plots ────────────────────────────────────────────────

cat("\n-- DMP Annotation Plots --\n")
if (nrow(dmp_anno) > 0) {
  anno_summary_dmp <- as.data.frame(table(dmp_anno$annotation))
  colnames(anno_summary_dmp) <- c("Feature", "Count")
  anno_summary_dmp$Percentage <- round(anno_summary_dmp$Count / sum(anno_summary_dmp$Count) * 100, 2)
  
  p_anno_pie <- ggplot(anno_summary_dmp, aes(x = "", y = Count, fill = Feature)) +
    geom_bar(stat = "identity", width = 1) + coord_polar("y") +
    scale_fill_brewer(palette = "Set2") +
    labs(title = "DMP Genomic Feature Distribution") +
    theme_void() +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  save_both(p_anno_pie, "dmp_annotation_pie", 8, 6)
  
  dmp_anno$direction <- ifelse(dmp_anno$methylation_diff > 0, "Hyper", "Hypo")
  anno_direction <- dmp_anno %>%
    group_by(annotation, direction) %>% summarise(n = n(), .groups = "drop")
  
  p_anno_bar <- ggplot(anno_direction, aes(x = annotation, y = n, fill = direction)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("Hyper" = "#E74C3C", "Hypo" = "#3498DB")) +
    labs(title = "DMP Distribution by Feature and Direction",
         x = "Genomic Feature", y = "Count", fill = "Direction") +
    theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_both(p_anno_bar, "dmp_annotation_direction_bar", 10, 6)
}

# ── B10. DMR Annotation Pie ─────────────────────────────────────────────────

cat("\n-- DMR Annotation Plots --\n")
if (nrow(dmr_anno) > 0) {
  anno_summary_dmr <- as.data.frame(table(dmr_anno$annotation))
  colnames(anno_summary_dmr) <- c("Feature", "Count")
  
  p_dmr_pie <- ggplot(anno_summary_dmr, aes(x = "", y = Count, fill = Feature)) +
    geom_bar(stat = "identity", width = 1) + coord_polar("y") +
    scale_fill_brewer(palette = "Set2") +
    labs(title = "DMR Genomic Feature Distribution") +
    theme_void() +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  save_both(p_dmr_pie, "dmr_annotation_pie", 8, 6)
}

# ── B11. Chromosome Distribution ────────────────────────────────────────────

cat("\n-- Chromosome Distribution --\n")
if (nrow(dmp_results) > 0) {
  chr_summary <- dmp_results %>%
    mutate(direction = ifelse(diff > 0, "Hyper", "Hypo"),
           chr_num = as.numeric(gsub("chr", "", chr))) %>%
    filter(!is.na(chr_num)) %>%
    group_by(chr, chr_num, direction) %>%
    summarise(count = n(), .groups = "drop") %>% arrange(chr_num)
  
  p_chr <- ggplot(chr_summary, aes(x = reorder(chr, chr_num), y = count, fill = direction)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("Hyper" = "#E74C3C", "Hypo" = "#3498DB")) +
    labs(title = "DMP Distribution by Chromosome",
         x = "Chromosome", y = "Number of DMPs", fill = "Direction") +
    theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_both(p_chr, "dmp_chromosome_distribution", 14, 6)
}

# ── B12. Top Genes Bar Plot (with display names) ────────────────────────────

cat("\n-- Top Genes Bar Plot --\n")
if (nrow(dmp_anno) > 0) {
  top_genes_plot <- gene_meth_summary %>%
    group_by(nearest_gene) %>%
    summarise(total_dmps = sum(n_dmps), mean_diff = mean(mean_diff),
              gene_type = first(gene_type), .groups = "drop") %>%
    arrange(desc(total_dmps)) %>% head(20) %>%
    left_join(gene_info[, c("gene_id", "display_name")],
              by = c("nearest_gene" = "gene_id"))
  
  top_genes_plot$label <- ifelse(
    !is.na(top_genes_plot$display_name) &
      top_genes_plot$display_name != "function unknown",
    substr(top_genes_plot$display_name, 1, 50),
    top_genes_plot$nearest_gene)
  
  p_top_genes <- ggplot(top_genes_plot, aes(x = reorder(label, total_dmps),
                                            y = total_dmps, fill = mean_diff)) +
    geom_bar(stat = "identity") +
    scale_fill_gradient2(low = "#3498DB", mid = "white", high = "#E74C3C",
                         midpoint = 0, name = "Mean\nMeth Diff") +
    coord_flip() +
    labs(title = "Top 20 Genes by Number of DMPs", x = "Gene", y = "Number of DMPs") +
    theme_minimal()
  save_both(p_top_genes, "top_genes_dmp_barplot", 14, 8)
}


# #############################################################################
# SECTION C: ORGANIZE OLD FILES INTO FOLDERS
# #############################################################################

cat("\n")
cat("================================================================\n")
cat("  SECTION C: Moving old flat files into folders\n")
cat("================================================================\n\n")

all_files <- list.files(RESULTS_DIR, full.names = TRUE, recursive = FALSE)
all_files <- all_files[!file.info(all_files)$isdir]

moved_count <- 0
for (f in all_files) {
  fname <- basename(f)
  ext   <- tools::file_ext(fname)
  dest  <- NULL
  
  if (ext == "pdf") {
    if (grepl("chipseeker", fname, ignore.case = TRUE)) {
      dest <- file.path(dirs$chipseeker, fname)
    } else if (grepl("enrichment|dotplot|barplot|kegg|reactome|go_", fname, ignore.case = TRUE)) {
      dest <- file.path(dirs$enrichment, fname)
    } else {
      dest <- file.path(dirs$pdf, fname)
    }
  } else if (ext == "rds") {
    dest <- file.path(dirs$rds, fname)
  } else if (ext == "txt") {
    if (grepl("enrichment", fname, ignore.case = TRUE)) {
      dest <- file.path(dirs$enrichment, fname)
    } else if (grepl("chipseeker", fname, ignore.case = TRUE)) {
      dest <- file.path(dirs$chipseeker, fname)
    } else {
      dest <- file.path(dirs$tables, fname)
    }
  }
  
  if (!is.null(dest)) {
    file.copy(f, dest, overwrite = TRUE)
    moved_count <- moved_count + 1
  }
}
cat(sprintf("  Organized %d files into subfolders.\n", moved_count))

# Optionally remove originals from flat directory (uncomment if desired)
# for (f in all_files) {
#   fname <- basename(f)
#   ext   <- tools::file_ext(fname)
#   if (ext %in% c("pdf", "rds", "txt")) file.remove(f)
# }


# #############################################################################
# FINAL SUMMARY
# #############################################################################

cat("\n\n")
cat("================================================================\n")
cat("           PART 4 COMPLETE                                      \n")
cat("================================================================\n")
cat("Results directory:", RESULTS_DIR, "\n\n")

cat("Folder structure:\n")
cat("  PDFs/            All main analysis PDFs\n")
cat("  PNGs/            All plots as high-res PNGs (300 dpi)\n")
cat("  Tables/          All text/TSV result tables\n")
cat("  RDS_Objects/     Saved R objects (BSseq, DMP, DMR, etc.)\n")
cat("  Enrichment/      GO, KEGG, Reactome results + plots (PDF+PNG)\n")
cat("  ChIPseeker/      ChIPseeker annotation + FIXED plots (PDF+PNG)\n\n")

for (dname in names(dirs)) {
  files <- list.files(dirs[[dname]])
  cat(sprintf("--- %s/ (%d files) ---\n", dname, length(files)))
  for (f in files) cat(sprintf("    %s\n", f))
  cat("\n")
}

cat("================================================================\n")
cat("Notes:\n")
cat("  - DMRs had limited enrichment (only 1 significant GO CC term:\n")
cat("    'Anchoring junction'). Expected with fewer DMR-associated\n")
cat("    genes (1,317) vs DMP genes (3,937).\n")
cat("  - ChIPseeker plots include ggplot2 fallback versions that\n")
cat("    always render: dmp/dmr_chipseeker_annotation_*.pdf/png\n")
cat("  - Original flat files kept as backup (uncomment removal block\n")
cat("    in Section C if you want to clean them up).\n")
cat("================================================================\n")