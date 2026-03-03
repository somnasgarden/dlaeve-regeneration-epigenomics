#!/usr/bin/env Rscript
# =============================================================================
# MXT Part 1: Methylation × Transcriptomics — DMP/DMR–Expression Correlation
# =============================================================================
# Goal: Correlate differential methylation (DMPs/DMRs) with gene expression
#        changes using ONLY tail amputation transcriptomes.
#
# Methylation samples:   C1, C2 (Control)  vs  A1, A2 (Amputated)
# Transcriptomics:       C1S1–C4S4 (Tail Ctrl) vs T2S6–T4S8 (Tail Amp)
#                        (T1S5 removed as outlier in WGCNA)
#
# Inputs:
#   - BSseq object, DMP/DMR results, ChIPseeker annotations (methylation pipeline)
#   - DESeq2 objects + normalized counts (transcriptomics pipeline, tail only)
#   - GFF annotation
#
# Outputs:
#   - Gene-level methylation–expression correlation tables
#   - Scatter plots by genomic region (Promoter, Exon, Intron, Intergenic, etc.)
#   - Heatmaps, barplots, and summary statistics
# =============================================================================

# ── Libraries ────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(bsseq)
  library(GenomicRanges)
  library(GenomicFeatures)
  library(rtracklayer)
  library(DESeq2)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(scales)
  library(ggrepel)
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
METH_DIR    <- "/mnt/data/alfredvar/rlopezt/Preliminary/"
TRANS_DIR   <- "/mnt/data/alfredvar/rlopezt/CorrelationMatrix/"
MXT_DIR     <- "/mnt/data/alfredvar/rlopezt/MXT/"

GFF_FILE    <- file.path(BASE_DIR,
                         "30-Genoma/31-Alternative_Annotation_EviAnn",
                         "derLaeGenome_namesDlasi_v2.fasta.functional_note.pseudo_label.gff")

# Create output directories
dirs <- list(
  main       = MXT_DIR,
  tables     = file.path(MXT_DIR, "Tables"),
  pdf        = file.path(MXT_DIR, "PDFs"),
  png        = file.path(MXT_DIR, "PNGs"),
  rds        = file.path(MXT_DIR, "RDS_Objects")
)
for (d in dirs) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# Helper: save ggplot as both PDF and PNG
save_both <- function(plot_obj, name, width, height) {
  ggsave(file.path(dirs$pdf, paste0(name, ".pdf")), plot_obj,
         width = width, height = height)
  ggsave(file.path(dirs$png, paste0(name, ".png")), plot_obj,
         width = width, height = height, dpi = 300)
  cat(sprintf("  Saved: %s (.pdf + .png)\n", name))
}

# Helper: save pheatmap as both
save_pheatmap_both <- function(pheatmap_expr, name, width, height) {
  pdf(file.path(dirs$pdf, paste0(name, ".pdf")), width = width, height = height)
  eval(pheatmap_expr)
  dev.off()
  png(file.path(dirs$png, paste0(name, ".png")),
      width = width, height = height, units = "in", res = 300)
  eval(pheatmap_expr)
  dev.off()
  cat(sprintf("  Saved: %s (.pdf + .png)\n", name))
}

cat("================================================================\n")
cat("  MXT Part 1: DMP/DMR × Expression Correlation\n")
cat("  Methylation × Transcriptomics Integration\n")
cat("================================================================\n\n")


# #############################################################################
# SECTION 1: LOAD METHYLATION DATA
# #############################################################################

cat("=== Loading Methylation Data ===\n\n")

bs_obj      <- readRDS(file.path(METH_DIR, "bsseq_object.rds"))
dmp_results <- readRDS(file.path(METH_DIR, "dmp_results.rds"))
dmr_results <- readRDS(file.path(METH_DIR, "dmr_results.rds"))

cat(sprintf("  BSseq: %s sites x %s samples\n",
            format(nrow(bs_obj), big.mark = ","), ncol(bs_obj)))
cat(sprintf("  DMPs: %s | DMRs: %s\n",
            format(nrow(dmp_results), big.mark = ","),
            format(nrow(dmr_results), big.mark = ",")))

# Load ChIPseeker annotations (produced by V1_4)
chipseeker_dir <- file.path(METH_DIR, "ChIPseeker")
dmp_cs_file <- file.path(chipseeker_dir, "dmps_chipseeker_annotated.txt")
dmr_cs_file <- file.path(chipseeker_dir, "dmrs_chipseeker_annotated.txt")

if (file.exists(dmp_cs_file)) {
  dmp_cs <- read.delim(dmp_cs_file, stringsAsFactors = FALSE)
  cat(sprintf("  ChIPseeker DMP annotations: %s rows\n", format(nrow(dmp_cs), big.mark = ",")))
} else {
  cat("  WARNING: ChIPseeker DMP file not found. Will annotate from scratch.\n")
  dmp_cs <- NULL
}

if (file.exists(dmr_cs_file)) {
  dmr_cs <- read.delim(dmr_cs_file, stringsAsFactors = FALSE)
  cat(sprintf("  ChIPseeker DMR annotations: %s rows\n", format(nrow(dmr_cs), big.mark = ",")))
} else {
  cat("  WARNING: ChIPseeker DMR file not found. Will annotate from scratch.\n")
  dmr_cs <- NULL
}

# Also load custom GFF-based annotations as fallback
dmp_anno_file <- file.path(METH_DIR, "dmps_annotated.txt")
dmr_anno_file <- file.path(METH_DIR, "dmrs_annotated.txt")

# Try tables subfolder if not in main
if (!file.exists(dmp_anno_file))
  dmp_anno_file <- file.path(METH_DIR, "Tables", "dmps_annotated.txt")
if (!file.exists(dmr_anno_file))
  dmr_anno_file <- file.path(METH_DIR, "Tables", "dmrs_annotated.txt")

dmp_anno <- if (file.exists(dmp_anno_file)) read.delim(dmp_anno_file, stringsAsFactors = FALSE) else NULL
dmr_anno <- if (file.exists(dmr_anno_file)) read.delim(dmr_anno_file, stringsAsFactors = FALSE) else NULL


# #############################################################################
# SECTION 2: LOAD TRANSCRIPTOMICS DATA — TAIL AMPUTATION ONLY
# #############################################################################

cat("\n=== Loading Transcriptomics Data (Tail Amputation Only) ===\n\n")

# Load DESeq2 objects from Part 2
trans_objects_file <- file.path(TRANS_DIR, "Part2_DEGs/data/Part2_objects.RData")
if (!file.exists(trans_objects_file)) {
  # Try Part1b
  trans_objects_file <- file.path(TRANS_DIR, "Part1b_BatchCorrection/data/Part1b_objects.RData")
}
cat(sprintf("  Loading: %s\n", trans_objects_file))
load(trans_objects_file)

cat(sprintf("  DESeq2 object: %s genes x %s samples\n",
            format(nrow(dds), big.mark = ","), ncol(dds)))

# Identify tail amputation samples ONLY
sample_metadata <- as.data.frame(colData(dds))
tail_samples <- rownames(sample_metadata)[sample_metadata$experiment == "tail_amputation"]
tail_ctrl    <- rownames(sample_metadata)[sample_metadata$experiment == "tail_amputation" &
                                            sample_metadata$condition == "control"]
tail_amp     <- rownames(sample_metadata)[sample_metadata$experiment == "tail_amputation" &
                                            sample_metadata$condition == "amputated"]

cat(sprintf("  Tail Control samples (n=%d): %s\n", length(tail_ctrl), paste(tail_ctrl, collapse = ", ")))
cat(sprintf("  Tail Amputated samples (n=%d): %s\n", length(tail_amp), paste(tail_amp, collapse = ", ")))

# Subset to tail amputation only
dds_tail <- dds[, c(tail_ctrl, tail_amp)]
colData(dds_tail)$condition <- droplevels(colData(dds_tail)$condition)

# Re-run DESeq2 on tail-only subset
cat("\n  Running DESeq2 on tail-only samples...\n")
design(dds_tail) <- ~ condition
dds_tail <- DESeq(dds_tail)

# Extract results: amputated vs control
res_tail <- results(dds_tail, contrast = c("condition", "amputated", "control"))
res_tail <- res_tail[!is.na(res_tail$padj), ]
res_tail_df <- as.data.frame(res_tail)
res_tail_df$gene_id <- rownames(res_tail_df)

cat(sprintf("  Tail DESeq2 results: %s genes tested\n", format(nrow(res_tail_df), big.mark = ",")))
cat(sprintf("  Tail DEGs (padj<0.05, |log2FC|>1): %d\n",
            sum(res_tail_df$padj < 0.05 & abs(res_tail_df$log2FoldChange) > 1, na.rm = TRUE)))

# Get normalized counts for tail samples
vsd_tail <- vst(dds_tail, blind = TRUE)
tail_counts <- assay(vsd_tail)

# Save tail-only results
write.table(res_tail_df, file.path(dirs$tables, "DESeq2_tail_amputated_vs_control.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(dds_tail, file.path(dirs$rds, "dds_tail_only.rds"))

cat("  Tail-only DESeq2 results saved.\n")


# #############################################################################
# SECTION 3: BUILD GENE-LEVEL METHYLATION–EXPRESSION TABLE
# #############################################################################

cat("\n=== Building Gene-Level Methylation–Expression Table ===\n\n")

# Step 3.1: Prepare DMP data with gene annotation and region
# Use ChIPseeker annotations (more detailed regions) if available

if (!is.null(dmp_cs)) {
  cat("  ChIPseeker DMP columns: ", paste(colnames(dmp_cs), collapse = ", "), "\n")
  
  # ChIPseeker file only has annotation columns — meth values come from dmp_results
  # Match by position (both are ordered the same from V1_4)
  if (nrow(dmp_cs) == nrow(dmp_results)) {
    dmp_gene_data <- data.frame(
      chr             = dmp_cs$seqnames,
      pos             = dmp_cs$start,
      meth_diff       = dmp_results$diff,
      pvalue          = dmp_results$pval,
      fdr             = dmp_results$fdr,
      gene_id         = dmp_cs$geneId,
      annotation_full = dmp_cs$annotation,
      stringsAsFactors = FALSE
    )
  } else {
    # Different lengths — merge by chr:pos
    cat("  WARNING: ChIPseeker and DMP results differ in length, merging by position...\n")
    dmp_cs$site_id <- paste(dmp_cs$seqnames, dmp_cs$start, sep = ":")
    dmp_results$site_id <- paste(dmp_results$chr, dmp_results$pos, sep = ":")
    dmp_merged <- merge(dmp_cs, dmp_results, by = "site_id")
    dmp_gene_data <- data.frame(
      chr             = dmp_merged$seqnames,
      pos             = dmp_merged$start,
      meth_diff       = dmp_merged$diff,
      pvalue          = dmp_merged$pval,
      fdr             = dmp_merged$fdr,
      gene_id         = dmp_merged$geneId,
      annotation_full = dmp_merged$annotation,
      stringsAsFactors = FALSE
    )
  }
  
  # Simplify ChIPseeker annotation to broad categories
  dmp_gene_data$region <- dmp_gene_data$annotation_full
  dmp_gene_data$region <- gsub(" \\(.*", "", dmp_gene_data$region)  # Remove parenthetical
  dmp_gene_data$region[grepl("Promoter", dmp_gene_data$region)] <- "Promoter"
  dmp_gene_data$region[grepl("Exon", dmp_gene_data$region)]     <- "Exon"
  dmp_gene_data$region[grepl("Intron", dmp_gene_data$region)]   <- "Intron"
  dmp_gene_data$region[grepl("UTR", dmp_gene_data$region)]      <- "UTR"
  dmp_gene_data$region[grepl("Downstream", dmp_gene_data$region)] <- "Downstream"
  dmp_gene_data$region[grepl("Intergenic|Distal", dmp_gene_data$region)] <- "Intergenic"
  
  cat("  DMP region distribution (ChIPseeker):\n")
  print(table(dmp_gene_data$region))
  
} else if (!is.null(dmp_anno)) {
  # Fallback to custom GFF annotation
  dmp_gene_data <- data.frame(
    chr        = dmp_anno$seqnames,
    pos        = dmp_anno$start,
    meth_diff  = dmp_anno$methylation_diff,
    pvalue     = if ("pvalue" %in% colnames(dmp_anno)) dmp_anno$pvalue else NA,
    fdr        = if ("fdr" %in% colnames(dmp_anno)) dmp_anno$fdr else NA,
    gene_id    = dmp_anno$nearest_gene,
    region     = dmp_anno$annotation,
    stringsAsFactors = FALSE
  )
  cat("  DMP region distribution (custom GFF):\n")
  print(table(dmp_gene_data$region))
  
} else {
  stop("No DMP annotation file found. Please run V1_3 or V1_4 first.")
}

# Remove entries with no gene assignment
dmp_gene_data <- dmp_gene_data[!is.na(dmp_gene_data$gene_id) & dmp_gene_data$gene_id != "", ]
cat(sprintf("\n  DMPs with gene assignment: %s\n", format(nrow(dmp_gene_data), big.mark = ",")))

# Step 3.2: Aggregate methylation by gene
gene_meth <- dmp_gene_data %>%
  group_by(gene_id) %>%
  summarise(
    n_dmps           = n(),
    mean_meth_diff   = mean(meth_diff, na.rm = TRUE),
    median_meth_diff = median(meth_diff, na.rm = TRUE),
    max_abs_diff     = max(abs(meth_diff), na.rm = TRUE),
    n_hyper          = sum(meth_diff > 0, na.rm = TRUE),
    n_hypo           = sum(meth_diff < 0, na.rm = TRUE),
    primary_region   = names(sort(table(region), decreasing = TRUE))[1],
    regions          = paste(unique(region), collapse = ";"),
    min_fdr          = min(fdr, na.rm = TRUE),
    .groups = "drop"
  )

cat(sprintf("  Unique genes with DMPs: %s\n", format(nrow(gene_meth), big.mark = ",")))

# Step 3.3: Aggregate methylation by gene AND region
gene_meth_by_region <- dmp_gene_data %>%
  group_by(gene_id, region) %>%
  summarise(
    n_dmps         = n(),
    mean_meth_diff = mean(meth_diff, na.rm = TRUE),
    .groups = "drop"
  )

# Step 3.4: Merge methylation with expression
mxt_data <- merge(gene_meth, res_tail_df,
                  by = "gene_id", all = FALSE)  # Inner join

cat(sprintf("  Genes with BOTH methylation and expression data: %s\n",
            format(nrow(mxt_data), big.mark = ",")))

# Classify expression changes
mxt_data$expr_sig <- mxt_data$padj < 0.05 & abs(mxt_data$log2FoldChange) > 1
mxt_data$meth_direction <- ifelse(mxt_data$mean_meth_diff > 0, "Hyper", "Hypo")
mxt_data$expr_direction <- ifelse(mxt_data$log2FoldChange > 0, "Up", "Down")

# Quadrant classification
mxt_data$quadrant <- paste(mxt_data$meth_direction, mxt_data$expr_direction, sep = "_")

cat("\n  Quadrant distribution:\n")
print(table(mxt_data$quadrant))

# By-region merged data
mxt_by_region <- merge(gene_meth_by_region, res_tail_df,
                       by = "gene_id", all = FALSE)

# Save tables
write.table(mxt_data, file.path(dirs$tables, "MXT_gene_level_meth_vs_expression.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(mxt_by_region, file.path(dirs$tables, "MXT_gene_region_level_meth_vs_expression.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)


# #############################################################################
# SECTION 4: DMR–EXPRESSION INTEGRATION
# #############################################################################

cat("\n=== DMR–Expression Integration ===\n\n")

if (!is.null(dmr_cs)) {
  cat("  ChIPseeker DMR columns: ", paste(colnames(dmr_cs), collapse = ", "), "\n")
  
  # ChIPseeker file only has annotation — meth values from dmr_results
  if (nrow(dmr_cs) == nrow(dmr_results)) {
    dmr_gene_data <- data.frame(
      chr            = dmr_cs$seqnames,
      start          = dmr_cs$start,
      end            = dmr_cs$end,
      gene_id        = dmr_cs$geneId,
      annotation     = dmr_cs$annotation,
      dmr_meth_diff  = dmr_results$diff.Methy,
      dmr_areaStat   = dmr_results$areaStat,
      dmr_nCG        = dmr_results$nCG,
      stringsAsFactors = FALSE
    )
  } else {
    cat("  WARNING: ChIPseeker and DMR results differ in length, merging by position...\n")
    dmr_cs$site_id <- paste(dmr_cs$seqnames, dmr_cs$start, dmr_cs$end, sep = ":")
    dmr_results$site_id <- paste(dmr_results$chr, dmr_results$start, dmr_results$end, sep = ":")
    dmr_merged <- merge(dmr_cs, dmr_results, by = "site_id")
    dmr_gene_data <- data.frame(
      chr            = dmr_merged$seqnames,
      start          = dmr_merged$start.x,
      end            = dmr_merged$end.x,
      gene_id        = dmr_merged$geneId,
      annotation     = dmr_merged$annotation,
      dmr_meth_diff  = dmr_merged$diff.Methy,
      dmr_areaStat   = dmr_merged$areaStat,
      dmr_nCG        = dmr_merged$nCG,
      stringsAsFactors = FALSE
    )
  }
} else if (!is.null(dmr_anno)) {
  dmr_gene_data <- data.frame(
    chr            = dmr_anno$seqnames,
    start          = dmr_anno$start,
    end            = dmr_anno$end,
    gene_id        = dmr_anno$nearest_gene,
    annotation     = dmr_anno$annotation,
    dmr_meth_diff  = if ("diff.Methy" %in% colnames(dmr_anno)) dmr_anno$diff.Methy else
      if ("methylation_diff" %in% colnames(dmr_anno)) dmr_anno$methylation_diff else NA,
    stringsAsFactors = FALSE
  )
} else {
  cat("  No DMR annotation available.\n")
  dmr_gene_data <- data.frame()
}

if (nrow(dmr_gene_data) > 0) {
  dmr_gene_data <- dmr_gene_data[!is.na(dmr_gene_data$gene_id) & dmr_gene_data$gene_id != "", ]
  
  # Simplify annotation
  dmr_gene_data$region <- dmr_gene_data$annotation
  dmr_gene_data$region <- gsub(" \\(.*", "", dmr_gene_data$region)
  dmr_gene_data$region[grepl("Promoter", dmr_gene_data$region)] <- "Promoter"
  dmr_gene_data$region[grepl("Exon", dmr_gene_data$region)]     <- "Exon"
  dmr_gene_data$region[grepl("Intron", dmr_gene_data$region)]   <- "Intron"
  dmr_gene_data$region[grepl("UTR", dmr_gene_data$region)]      <- "UTR"
  dmr_gene_data$region[grepl("Downstream", dmr_gene_data$region)] <- "Downstream"
  dmr_gene_data$region[grepl("Intergenic|Distal", dmr_gene_data$region)] <- "Intergenic"
  
  # Aggregate DMRs by gene
  gene_dmr <- dmr_gene_data %>%
    group_by(gene_id) %>%
    summarise(
      n_dmrs           = n(),
      mean_dmr_diff    = mean(dmr_meth_diff, na.rm = TRUE),
      total_areaStat   = sum(dmr_areaStat, na.rm = TRUE),
      primary_region   = names(sort(table(region), decreasing = TRUE))[1],
      .groups = "drop"
    )
  
  # Merge DMRs with expression
  mxt_dmr <- merge(gene_dmr, res_tail_df, by = "gene_id", all = FALSE)
  
  cat(sprintf("  DMRs with gene assignment: %s\n", nrow(dmr_gene_data)))
  cat(sprintf("  Genes with both DMR and expression data: %s\n", nrow(mxt_dmr)))
  
  write.table(mxt_dmr, file.path(dirs$tables, "MXT_DMR_gene_level_vs_expression.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}


# #############################################################################
# SECTION 5: CORRELATION ANALYSIS
# #############################################################################

cat("\n=== Correlation Analysis ===\n\n")

# 5.1 Overall methylation vs expression correlation
cor_overall <- cor.test(mxt_data$mean_meth_diff, mxt_data$log2FoldChange,
                        method = "spearman")
cat(sprintf("  Overall Spearman correlation: rho = %.4f, p = %.2e\n",
            cor_overall$estimate, cor_overall$p.value))

cor_pearson <- cor.test(mxt_data$mean_meth_diff, mxt_data$log2FoldChange,
                        method = "pearson")
cat(sprintf("  Overall Pearson correlation:  r = %.4f, p = %.2e\n",
            cor_pearson$estimate, cor_pearson$p.value))

# 5.2 Correlation by genomic region
cat("\n  Correlations by genomic region:\n")
region_cors <- mxt_by_region %>%
  group_by(region) %>%
  summarise(
    n_genes     = n(),
    spearman_r  = cor(mean_meth_diff, log2FoldChange, method = "spearman", use = "complete.obs"),
    spearman_p  = tryCatch(cor.test(mean_meth_diff, log2FoldChange, method = "spearman")$p.value,
                           error = function(e) NA),
    pearson_r   = cor(mean_meth_diff, log2FoldChange, method = "pearson", use = "complete.obs"),
    pearson_p   = tryCatch(cor.test(mean_meth_diff, log2FoldChange, method = "pearson")$p.value,
                           error = function(e) NA),
    .groups = "drop"
  ) %>%
  arrange(spearman_p)

print(as.data.frame(region_cors))

write.table(region_cors, file.path(dirs$tables, "MXT_correlation_by_region.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# 5.3 Summary statistics table
summary_stats <- data.frame(
  Metric = c(
    "Total genes with DMP + expression data",
    "Total genes with DMR + expression data",
    paste0("Overall Spearman (meth_diff vs log2FC): rho=",
           round(cor_overall$estimate, 4), ", p=", signif(cor_overall$p.value, 3)),
    paste0("Overall Pearson (meth_diff vs log2FC): r=",
           round(cor_pearson$estimate, 4), ", p=", signif(cor_pearson$p.value, 3)),
    "Genes: Hyper-methylated + Upregulated",
    "Genes: Hyper-methylated + Downregulated",
    "Genes: Hypo-methylated + Upregulated",
    "Genes: Hypo-methylated + Downregulated"
  ),
  Value = c(
    nrow(mxt_data),
    if (exists("mxt_dmr")) nrow(mxt_dmr) else 0,
    "",
    "",
    sum(mxt_data$quadrant == "Hyper_Up", na.rm = TRUE),
    sum(mxt_data$quadrant == "Hyper_Down", na.rm = TRUE),
    sum(mxt_data$quadrant == "Hypo_Up", na.rm = TRUE),
    sum(mxt_data$quadrant == "Hypo_Down", na.rm = TRUE)
  )
)

write.table(summary_stats, file.path(dirs$tables, "MXT_summary_statistics.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)


# #############################################################################
# SECTION 6: VISUALIZATIONS — DMP × EXPRESSION
# #############################################################################

cat("\n=== Generating Visualizations ===\n\n")

# 6.1 Main scatter: Methylation Difference vs Log2FC
p_main <- ggplot(mxt_data, aes(x = mean_meth_diff, y = log2FoldChange)) +
  geom_point(aes(color = primary_region), alpha = 0.4, size = 1.2) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8, linetype = "dashed") +
  geom_hline(yintercept = 0, color = "gray50", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "gray50", linetype = "dashed") +
  scale_color_brewer(palette = "Set2", name = "Genomic Region") +
  labs(
    title = "Methylation Difference vs Gene Expression Change — D. laeve",
    subtitle = sprintf("Tail amputation only | n=%s genes | Spearman ρ = %.3f (p = %.2e)",
                       format(nrow(mxt_data), big.mark = ","),
                       cor_overall$estimate, cor_overall$p.value),
    x = "Mean Methylation Difference (Amputated − Control)",
    y = "Log2 Fold Change (Expression)"
  )
save_both(p_main, "MXT_01_meth_diff_vs_log2FC_all", 12, 9)

# 6.2 Scatter faceted by genomic region
p_facet <- ggplot(mxt_by_region, aes(x = mean_meth_diff, y = log2FoldChange)) +
  geom_point(aes(color = region), alpha = 0.4, size = 1) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.7) +
  geom_hline(yintercept = 0, color = "gray50", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "gray50", linetype = "dashed") +
  facet_wrap(~ region, scales = "free") +
  scale_color_brewer(palette = "Set2") +
  labs(
    title = "Methylation vs Expression by Genomic Region — D. laeve",
    subtitle = "Tail amputation: Amputated vs Control",
    x = "Mean Methylation Difference at DMPs",
    y = "Log2 Fold Change (Expression)"
  ) +
  theme(legend.position = "none")
save_both(p_facet, "MXT_02_meth_vs_expr_by_region", 16, 12)

# 6.3 Correlation barplot by region
p_cor_bar <- ggplot(region_cors, aes(x = reorder(region, spearman_r),
                                     y = spearman_r, fill = spearman_r)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("r=%.3f\np=%.2e", spearman_r, spearman_p)),
            hjust = ifelse(region_cors$spearman_r > 0, -0.1, 1.1),
            size = 3) +
  scale_fill_gradient2(low = "#3498DB", mid = "white", high = "#E74C3C",
                       midpoint = 0, name = "Spearman ρ") +
  coord_flip() +
  labs(
    title = "Methylation–Expression Correlation by Genomic Region",
    subtitle = "Spearman correlation: DMPs in each region vs gene expression",
    x = "Genomic Region", y = "Spearman Correlation (ρ)"
  )
save_both(p_cor_bar, "MXT_03_correlation_by_region_barplot", 12, 7)

# 6.4 Quadrant plot with density
p_quad <- ggplot(mxt_data, aes(x = mean_meth_diff, y = log2FoldChange)) +
  geom_hex(bins = 50) +
  scale_fill_viridis_c(name = "Count", option = "inferno") +
  geom_hline(yintercept = 0, color = "white", linewidth = 0.8) +
  geom_vline(xintercept = 0, color = "white", linewidth = 0.8) +
  annotate("text", x = max(mxt_data$mean_meth_diff) * 0.7,
           y = max(mxt_data$log2FoldChange) * 0.9,
           label = paste("Hyper + Up\nn =", sum(mxt_data$quadrant == "Hyper_Up")),
           color = "white", size = 4, fontface = "bold") +
  annotate("text", x = max(mxt_data$mean_meth_diff) * 0.7,
           y = min(mxt_data$log2FoldChange) * 0.9,
           label = paste("Hyper + Down\nn =", sum(mxt_data$quadrant == "Hyper_Down")),
           color = "white", size = 4, fontface = "bold") +
  annotate("text", x = min(mxt_data$mean_meth_diff) * 0.7,
           y = max(mxt_data$log2FoldChange) * 0.9,
           label = paste("Hypo + Up\nn =", sum(mxt_data$quadrant == "Hypo_Up")),
           color = "white", size = 4, fontface = "bold") +
  annotate("text", x = min(mxt_data$mean_meth_diff) * 0.7,
           y = min(mxt_data$log2FoldChange) * 0.9,
           label = paste("Hypo + Down\nn =", sum(mxt_data$quadrant == "Hypo_Down")),
           color = "white", size = 4, fontface = "bold") +
  labs(
    title = "Methylation–Expression Quadrant Plot — D. laeve",
    subtitle = "Tail amputation: density of genes by methylation change × expression change",
    x = "Mean Methylation Difference (Amputated − Control)",
    y = "Log2 Fold Change (Expression)"
  )
save_both(p_quad, "MXT_04_quadrant_density", 11, 9)

# 6.5 Promoter-specific analysis (classic epigenomic expectation)
promoter_data <- mxt_by_region %>% filter(region == "Promoter")

if (nrow(promoter_data) > 10) {
  cor_promoter <- cor.test(promoter_data$mean_meth_diff, promoter_data$log2FoldChange,
                           method = "spearman")
  
  p_promoter <- ggplot(promoter_data, aes(x = mean_meth_diff, y = log2FoldChange)) +
    geom_point(alpha = 0.5, color = "#E74C3C", size = 1.5) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    labs(
      title = "Promoter Methylation vs Gene Expression — D. laeve",
      subtitle = sprintf("n=%d genes | Spearman ρ = %.3f (p = %.2e)\nExpectation: promoter hypermethylation → transcriptional silencing",
                         nrow(promoter_data), cor_promoter$estimate, cor_promoter$p.value),
      x = "Mean Methylation Difference at Promoter DMPs",
      y = "Log2 Fold Change (Expression)"
    )
  save_both(p_promoter, "MXT_05_promoter_meth_vs_expression", 10, 8)
}

# 6.6 Gene body analysis (inverse expectation: gene body meth → active transcription)
genebody_data <- mxt_by_region %>% filter(region %in% c("Exon", "Intron"))

if (nrow(genebody_data) > 10) {
  cor_genebody <- cor.test(genebody_data$mean_meth_diff, genebody_data$log2FoldChange,
                           method = "spearman")
  
  p_genebody <- ggplot(genebody_data, aes(x = mean_meth_diff, y = log2FoldChange, color = region)) +
    geom_point(alpha = 0.4, size = 1.5) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    scale_color_manual(values = c("Exon" = "#2E86C1", "Intron" = "#27AE60")) +
    labs(
      title = "Gene Body Methylation vs Gene Expression — D. laeve",
      subtitle = sprintf("n=%d gene-region pairs | Spearman ρ = %.3f (p = %.2e)\nExpectation: gene body methylation can correlate positively with expression",
                         nrow(genebody_data), cor_genebody$estimate, cor_genebody$p.value),
      x = "Mean Methylation Difference at Gene Body DMPs",
      y = "Log2 Fold Change (Expression)"
    )
  save_both(p_genebody, "MXT_06_genebody_meth_vs_expression", 10, 8)
}

# 6.7 Highlight significant genes (both methylation and expression)
sig_both <- mxt_data %>%
  filter(padj < 0.05, abs(log2FoldChange) > 1, n_dmps >= 3)

if (nrow(sig_both) > 0) {
  # Load gene names if available
  gene_info_file <- file.path(METH_DIR, "gene_info_with_names.txt")
  if (!file.exists(gene_info_file))
    gene_info_file <- file.path(METH_DIR, "Tables", "gene_info_with_names.txt")
  
  if (file.exists(gene_info_file)) {
    gene_names <- read.delim(gene_info_file, stringsAsFactors = FALSE)
    sig_both <- merge(sig_both, gene_names[, c("gene_id", "display_name")],
                      by = "gene_id", all.x = TRUE)
    sig_both$label <- ifelse(!is.na(sig_both$display_name) &
                               sig_both$display_name != "function unknown",
                             substr(sig_both$display_name, 1, 40),
                             sig_both$gene_id)
  } else {
    sig_both$label <- sig_both$gene_id
  }
  
  p_sig <- ggplot(mxt_data, aes(x = mean_meth_diff, y = log2FoldChange)) +
    geom_point(alpha = 0.2, color = "gray70", size = 0.8) +
    geom_point(data = sig_both, aes(color = primary_region), size = 2.5, alpha = 0.8) +
    geom_text_repel(data = head(sig_both[order(-abs(sig_both$log2FoldChange)), ], 20),
                    aes(label = label), size = 2.5, max.overlaps = 15) +
    scale_color_brewer(palette = "Set2", name = "Region") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    labs(
      title = "Significant Genes: Differential Methylation + Expression — D. laeve",
      subtitle = sprintf("padj<0.05, |log2FC|>1, ≥3 DMPs | n=%d genes highlighted", nrow(sig_both)),
      x = "Mean Methylation Difference",
      y = "Log2 Fold Change"
    )
  save_both(p_sig, "MXT_07_significant_genes_labeled", 14, 10)
  
  write.table(sig_both, file.path(dirs$tables, "MXT_significant_meth_AND_expr_genes.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

# 6.8 DMR correlation plot (if DMRs exist)
if (exists("mxt_dmr") && nrow(mxt_dmr) > 10) {
  cor_dmr <- cor.test(mxt_dmr$mean_dmr_diff, mxt_dmr$log2FoldChange,
                      method = "spearman")
  
  p_dmr <- ggplot(mxt_dmr, aes(x = mean_dmr_diff, y = log2FoldChange)) +
    geom_point(aes(size = n_dmrs, color = primary_region), alpha = 0.5) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
    scale_color_brewer(palette = "Set2", name = "Region") +
    scale_size_continuous(name = "# DMRs", range = c(1, 5)) +
    labs(
      title = "DMR Methylation vs Gene Expression — D. laeve",
      subtitle = sprintf("n=%d genes | Spearman ρ = %.3f (p = %.2e)",
                         nrow(mxt_dmr), cor_dmr$estimate, cor_dmr$p.value),
      x = "Mean DMR Methylation Difference",
      y = "Log2 Fold Change (Expression)"
    )
  save_both(p_dmr, "MXT_08_DMR_meth_vs_expression", 12, 9)
}

# 6.9 Heatmap: Top genes by combined signal
top_combined <- mxt_data %>%
  filter(!is.na(padj)) %>%
  mutate(combined_score = abs(mean_meth_diff) * abs(log2FoldChange) / padj) %>%
  arrange(desc(combined_score)) %>%
  head(50)

if (nrow(top_combined) > 5) {
  # Get expression for these genes across tail samples
  top_expr <- tail_counts[top_combined$gene_id[top_combined$gene_id %in% rownames(tail_counts)], ]
  
  if (nrow(top_expr) > 2) {
    col_ann <- data.frame(
      Condition = c(rep("Control", length(tail_ctrl)), rep("Amputated", length(tail_amp))),
      row.names = c(tail_ctrl, tail_amp)
    )
    col_colors <- list(Condition = c("Control" = "#2E86C1", "Amputated" = "#E74C3C"))
    
    hm_expr <- quote(
      pheatmap(top_expr[, c(tail_ctrl, tail_amp)],
               scale = "row",
               annotation_col = col_ann,
               annotation_colors = col_colors,
               show_rownames = FALSE,
               cluster_cols = FALSE,
               color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
               main = "Top 50 Genes: Highest Combined Meth+Expr Signal — Tail")
    )
    save_pheatmap_both(hm_expr, "MXT_09_heatmap_top_combined_genes", 10, 12)
  }
}

# 6.10 Region distribution among significant overlap genes
if (nrow(sig_both) > 0) {
  region_sig <- as.data.frame(table(sig_both$primary_region))
  colnames(region_sig) <- c("Region", "Count")
  
  p_region_sig <- ggplot(region_sig, aes(x = "", y = Count, fill = Region)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y") +
    scale_fill_brewer(palette = "Set2") +
    labs(title = "Genomic Region Distribution of Genes with\nSignificant Methylation AND Expression Changes") +
    theme_void() +
    theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5))
  save_both(p_region_sig, "MXT_10_region_pie_significant_genes", 9, 7)
}


# #############################################################################
# SECTION 7: SAVE ALL OBJECTS
# #############################################################################

cat("\n=== Saving R Objects ===\n\n")

save(mxt_data, mxt_by_region,
     res_tail_df, dmp_gene_data,
     gene_meth, region_cors,
     file = file.path(dirs$rds, "MXT_Part1_objects.RData"))

if (exists("mxt_dmr")) {
  save(mxt_dmr, dmr_gene_data, gene_dmr,
       file = file.path(dirs$rds, "MXT_Part1_DMR_objects.RData"))
}

cat("  Objects saved.\n")


# #############################################################################
# FINAL SUMMARY
# #############################################################################

cat("\n")
cat("================================================================\n")
cat("   MXT PART 1 COMPLETE — DMP/DMR × Expression Correlation\n")
cat("================================================================\n\n")

cat("KEY RESULTS:\n")
cat(sprintf("  Genes with DMP + expression: %s\n", format(nrow(mxt_data), big.mark = ",")))
if (exists("mxt_dmr")) {
  cat(sprintf("  Genes with DMR + expression: %s\n", format(nrow(mxt_dmr), big.mark = ",")))
}
cat(sprintf("  Overall Spearman rho: %.4f (p = %.2e)\n", cor_overall$estimate, cor_overall$p.value))
cat(sprintf("  Genes with sig. meth + sig. expr: %d\n", nrow(sig_both)))
cat("\n  Correlation by region:\n")
print(as.data.frame(region_cors[, c("region", "n_genes", "spearman_r", "spearman_p")]))

cat("\nOUTPUT FILES:\n")
cat(sprintf("  Tables:  %s\n", dirs$tables))
cat(sprintf("  PDFs:    %s\n", dirs$pdf))
cat(sprintf("  PNGs:    %s\n", dirs$png))
cat(sprintf("  RDS:     %s\n", dirs$rds))
cat("================================================================\n")