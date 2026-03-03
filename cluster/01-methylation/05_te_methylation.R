#!/usr/bin/env Rscript
# =============================================================================
# TE-METHYLATION ANALYSIS PIPELINE - OPTIMIZED VERSION
# =============================================================================
# FAST vectorized approach - should run in ~10-30 minutes, not 20+ hours
# =============================================================================

suppressPackageStartupMessages({
  library(bsseq)
  library(GenomicRanges)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggridges)
  library(scales)
  library(RColorBrewer)
})

options(stringsAsFactors = FALSE)
options(scipen = 999)

theme_set(theme_minimal() +
            theme(text = element_text(size = 12),
                  plot.title = element_text(size = 14, face = "bold")))

# ── Paths ────────────────────────────────────────────────────────────────────

RESULTS_DIR <- "/mnt/data/alfredvar/rlopezt/Preliminary/"
TE_FILE     <- "/mnt/data/alfredvar/30-Genoma/32-Repeats/age_of_transposons/collapsed_te_age_data.tsv"

OUTPUT_DIR  <- file.path(RESULTS_DIR, "TE_Analysis")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Analysis parameters
ANALYSIS_PARAMS <- list(
  min_diff       = 0.1,
  fdr_threshold  = 0.05
)

# =============================================================================
# SECTION 1: DATA LOADING
# =============================================================================

cat("=== Loading Data ===\n\n")

# Chromosome filter — chr1-31 only (exclude scaffolds)
keep_chr <- paste0("chr", 1:31)

bs_obj      <- readRDS(file.path(RESULTS_DIR, "bsseq_object.rds"))
dmp_results <- readRDS(file.path(RESULTS_DIR, "dmp_results.rds"))
dmr_results <- readRDS(file.path(RESULTS_DIR, "dmr_results.rds"))

# Filter to chr1-31
dmp_results <- dmp_results[dmp_results$chr %in% keep_chr, ]
dmr_results <- dmr_results[dmr_results$chr %in% keep_chr, ]
cat("  Filtered DMP/DMR data to chr1-31\n")

cat(sprintf("BSseq object: %s CpG sites x %s samples\n",
            format(nrow(bs_obj), big.mark = ","), ncol(bs_obj)))
cat(sprintf("DMPs: %s total\n", format(nrow(dmp_results), big.mark = ",")))
cat(sprintf("DMRs: %s total\n", format(nrow(dmr_results), big.mark = ",")))

# Load TE annotation
cat("\nLoading TE annotation...\n")
te_data <- fread(TE_FILE, sep = "\t", header = TRUE)

if ("chrom" %in% colnames(te_data)) {
  setnames(te_data, "chrom", "chr")
}

te_data[, c("te_class", "te_family") := tstrsplit(class_family, "/", fixed = TRUE)]
te_data[is.na(te_family), te_family := te_class]

# Filter to chr1-31
te_data <- te_data[chr %in% keep_chr]
te_data[, te_idx := .I]  # Add index column after filtering

cat(sprintf("TE annotations: %s TEs loaded (chr1-31)\n", format(nrow(te_data), big.mark = ",")))
cat("\nTE class distribution:\n")
print(table(te_data$te_class))

# =============================================================================
# SECTION 2: CREATE GENOMIC RANGES
# =============================================================================

cat("\n=== Creating GenomicRanges Objects ===\n\n")

cpg_gr <- granges(bs_obj)
cat(sprintf("CpG sites GRanges: %s positions\n", format(length(cpg_gr), big.mark = ",")))

# DMPs
dmp_gr <- GRanges(
  seqnames = dmp_results$chr,
  ranges   = IRanges(start = dmp_results$pos, width = 1),
  diff     = dmp_results$diff,
  pval     = dmp_results$pval,
  fdr      = dmp_results$fdr
)
mcols(dmp_gr)$direction <- ifelse(dmp_results$diff > 0, "Hyper", "Hypo")
mcols(dmp_gr)$significant <- dmp_results$fdr < ANALYSIS_PARAMS$fdr_threshold &
  abs(dmp_results$diff) > ANALYSIS_PARAMS$min_diff

cat(sprintf("DMP GRanges: %s positions\n", format(length(dmp_gr), big.mark = ",")))
cat(sprintf("  - Significant: %s\n", sum(mcols(dmp_gr)$significant)))

# DMRs
dmr_gr <- GRanges(
  seqnames  = dmr_results$chr,
  ranges    = IRanges(start = dmr_results$start, end = dmr_results$end),
  diff_methy = dmr_results$diff.Methy,
  areaStat  = dmr_results$areaStat,
  nCG       = dmr_results$nCG
)
mcols(dmr_gr)$direction <- ifelse(dmr_results$diff.Methy > 0, "Hyper", "Hypo")
cat(sprintf("DMR GRanges: %s regions\n", format(length(dmr_gr), big.mark = ",")))

# TEs
te_gr <- GRanges(
  seqnames = te_data$chr,
  ranges   = IRanges(start = te_data$start, end = te_data$end),
  te_idx   = te_data$te_idx,
  te_name  = te_data$te_name,
  te_class = te_data$te_class,
  te_family = te_data$te_family,
  kimura_div = te_data$kimura_div
)
cat(sprintf("TE GRanges: %s elements\n", format(length(te_gr), big.mark = ",")))

# =============================================================================
# SECTION 3: FAST VECTORIZED METHYLATION CALCULATION
# =============================================================================

cat("\n=== Calculating Methylation Per TE (FAST vectorized) ===\n\n")

# Get methylation matrix
beta_matrix <- getMeth(bs_obj, type = "raw")
colnames(beta_matrix) <- colnames(bs_obj)

control_samples <- c("C1", "C2")
treatment_samples <- c("A1", "A2")

# Find all CpG-TE overlaps
cat("Finding CpG-TE overlaps...\n")
overlaps <- findOverlaps(cpg_gr, te_gr)
cat(sprintf("  Total overlaps: %s\n", format(length(overlaps), big.mark = ",")))

# FAST: Extract data and use data.table for aggregation
cat("Extracting methylation values (vectorized)...\n")

cpg_idx <- queryHits(overlaps)
te_idx <- subjectHits(overlaps)

# Get methylation values for all overlapping CpGs at once
meth_control <- rowMeans(beta_matrix[cpg_idx, control_samples, drop = FALSE], na.rm = TRUE)
meth_treatment <- rowMeans(beta_matrix[cpg_idx, treatment_samples, drop = FALSE], na.rm = TRUE)

# Create data.table for fast aggregation
cat("Aggregating by TE (data.table)...\n")
overlap_dt <- data.table(
  te_idx = te_idx,
  meth_control = meth_control,
  meth_treatment = meth_treatment
)

# Aggregate by TE - MUCH faster than looping
te_meth_agg <- overlap_dt[, .(
  n_cpgs = .N,
  mean_meth_control = mean(meth_control, na.rm = TRUE),
  mean_meth_treatment = mean(meth_treatment, na.rm = TRUE)
), by = te_idx]

te_meth_agg[, meth_diff := mean_meth_treatment - mean_meth_control]
te_meth_agg[, mean_meth_all := (mean_meth_control + mean_meth_treatment) / 2]

# Merge with TE info
te_info_dt <- data.table(
  te_idx = te_data$te_idx,
  te_name = te_data$te_name,
  te_class = te_data$te_class,
  te_family = te_data$te_family,
  te_length = te_data$end - te_data$start,
  kimura_div = te_data$kimura_div
)

te_methylation <- merge(te_info_dt, te_meth_agg, by = "te_idx", all.x = TRUE)
te_methylation[is.na(n_cpgs), n_cpgs := 0]

cat(sprintf("\nTEs with CpG sites: %s / %s (%.1f%%)\n",
            sum(te_methylation$n_cpgs > 0),
            nrow(te_methylation),
            100 * sum(te_methylation$n_cpgs > 0) / nrow(te_methylation)))

# Save
fwrite(te_methylation, file.path(OUTPUT_DIR, "TE_methylation_per_element.txt"), sep = "\t")
cat("TE methylation data saved.\n")

# =============================================================================
# SECTION 4: METHYLATION BY CLASS AND FAMILY (Fast)
# =============================================================================

cat("\n=== Methylation by TE Class and Family ===\n\n")

# By CLASS
meth_by_class <- te_methylation[n_cpgs > 0, .(
  n_tes = .N,
  total_cpgs = sum(n_cpgs),
  mean_meth_control = mean(mean_meth_control, na.rm = TRUE),
  mean_meth_treatment = mean(mean_meth_treatment, na.rm = TRUE),
  mean_meth_diff = mean(meth_diff, na.rm = TRUE),
  sd_meth_control = sd(mean_meth_control, na.rm = TRUE)
), by = te_class][order(-n_tes)]

cat("Methylation by TE Class:\n")
print(meth_by_class)
fwrite(meth_by_class, file.path(OUTPUT_DIR, "methylation_by_TE_class.txt"), sep = "\t")

# By FAMILY
meth_by_family <- te_methylation[n_cpgs > 0, .(
  n_tes = .N,
  total_cpgs = sum(n_cpgs),
  mean_meth_control = mean(mean_meth_control, na.rm = TRUE),
  mean_meth_treatment = mean(mean_meth_treatment, na.rm = TRUE),
  mean_meth_diff = mean(meth_diff, na.rm = TRUE)
), by = .(te_class, te_family)][order(te_class, -n_tes)]

fwrite(meth_by_family, file.path(OUTPUT_DIR, "methylation_by_TE_family.txt"), sep = "\t")

# =============================================================================
# SECTION 5: DMPs IN TEs (Fast)
# =============================================================================

cat("\n=== Analyzing DMPs in TEs ===\n\n")

dmp_te_overlaps <- findOverlaps(dmp_gr, te_gr)

cat(sprintf("Total DMP-TE overlaps: %s\n", length(dmp_te_overlaps)))
cat(sprintf("Unique DMPs in TEs: %s / %s (%.1f%%)\n",
            length(unique(queryHits(dmp_te_overlaps))),
            length(dmp_gr),
            100 * length(unique(queryHits(dmp_te_overlaps))) / length(dmp_gr)))

# Create DMP-TE mapping
dmp_te_df <- data.table(
  dmp_idx = queryHits(dmp_te_overlaps),
  te_idx = subjectHits(dmp_te_overlaps)
)

# Add DMP info
dmp_te_df[, `:=`(
  chr = as.character(seqnames(dmp_gr))[dmp_idx],
  pos = start(dmp_gr)[dmp_idx],
  diff = mcols(dmp_gr)$diff[dmp_idx],
  fdr = mcols(dmp_gr)$fdr[dmp_idx],
  direction = mcols(dmp_gr)$direction[dmp_idx],
  significant = mcols(dmp_gr)$significant[dmp_idx]
)]

# Add TE info
dmp_te_df[, `:=`(
  te_name = mcols(te_gr)$te_name[te_idx],
  te_class = mcols(te_gr)$te_class[te_idx],
  te_family = mcols(te_gr)$te_family[te_idx]
)]

fwrite(dmp_te_df, file.path(OUTPUT_DIR, "DMPs_in_TEs_full.txt"), sep = "\t")

# Summary
te_methylation[, has_dmp := te_idx %in% unique(dmp_te_df$te_idx)]
te_methylation[, has_significant_dmp := te_idx %in% unique(dmp_te_df[significant == TRUE]$te_idx)]

cat(sprintf("\nTEs with any DMP: %s (%.1f%%)\n",
            sum(te_methylation$has_dmp),
            100 * sum(te_methylation$has_dmp) / nrow(te_methylation)))
cat(sprintf("TEs with significant DMP: %s (%.1f%%)\n",
            sum(te_methylation$has_significant_dmp),
            100 * sum(te_methylation$has_significant_dmp) / nrow(te_methylation)))

# =============================================================================
# SECTION 6: DMRs IN TEs (Fast)
# =============================================================================

cat("\n=== Analyzing DMRs in TEs ===\n\n")

dmr_te_overlaps <- findOverlaps(dmr_gr, te_gr)

cat(sprintf("Total DMR-TE overlaps: %s\n", length(dmr_te_overlaps)))
cat(sprintf("Unique DMRs overlapping TEs: %s / %s (%.1f%%)\n",
            length(unique(queryHits(dmr_te_overlaps))),
            length(dmr_gr),
            100 * length(unique(queryHits(dmr_te_overlaps))) / length(dmr_gr)))

dmr_te_df <- data.table(
  dmr_idx = queryHits(dmr_te_overlaps),
  te_idx = subjectHits(dmr_te_overlaps)
)

dmr_te_df[, `:=`(
  chr = as.character(seqnames(dmr_gr))[dmr_idx],
  dmr_start = start(dmr_gr)[dmr_idx],
  dmr_end = end(dmr_gr)[dmr_idx],
  diff_methy = mcols(dmr_gr)$diff_methy[dmr_idx],
  areaStat = mcols(dmr_gr)$areaStat[dmr_idx],
  direction = mcols(dmr_gr)$direction[dmr_idx]
)]

dmr_te_df[, `:=`(
  te_name = mcols(te_gr)$te_name[te_idx],
  te_class = mcols(te_gr)$te_class[te_idx],
  te_family = mcols(te_gr)$te_family[te_idx],
  te_start = start(te_gr)[te_idx],
  te_end = end(te_gr)[te_idx]
)]

dmr_te_df[, overlap_bp := pmin(dmr_end, te_end) - pmax(dmr_start, te_start) + 1]

fwrite(dmr_te_df, file.path(OUTPUT_DIR, "DMRs_in_TEs_full.txt"), sep = "\t")

te_methylation[, has_dmr := te_idx %in% unique(dmr_te_df$te_idx)]

cat(sprintf("TEs with DMR overlap: %s (%.1f%%)\n",
            sum(te_methylation$has_dmr),
            100 * sum(te_methylation$has_dmr) / nrow(te_methylation)))

# =============================================================================
# SECTION 7: SUMMARY STATISTICS
# =============================================================================

cat("\n=== Summary Statistics ===\n\n")

te_methylation[, category := fcase(
  has_significant_dmp | has_dmr, "With significant changes",
  has_dmp, "With non-significant DMPs only",
  default = "No methylation changes"
)]

cat("TE categories:\n")
print(table(te_methylation$category))

# Class summary
te_class_summary <- te_methylation[n_cpgs > 0, .(
  n_total = .N,
  n_with_dmp = sum(has_dmp),
  n_with_sig_dmp = sum(has_significant_dmp),
  n_with_dmr = sum(has_dmr),
  pct_with_dmp = 100 * sum(has_dmp) / .N,
  pct_with_sig_dmp = 100 * sum(has_significant_dmp) / .N,
  pct_with_dmr = 100 * sum(has_dmr) / .N,
  mean_meth = mean(mean_meth_all, na.rm = TRUE)
), by = te_class][order(-n_total)]

cat("\nTE class summary:\n")
print(te_class_summary)
fwrite(te_class_summary, file.path(OUTPUT_DIR, "TE_class_change_summary.txt"), sep = "\t")

# Family summary
te_family_summary <- te_methylation[n_cpgs > 0, .(
  n_total = .N,
  n_with_dmp = sum(has_dmp),
  n_with_sig_dmp = sum(has_significant_dmp),
  n_with_dmr = sum(has_dmr),
  pct_with_dmp = 100 * sum(has_dmp) / .N,
  mean_meth = mean(mean_meth_all, na.rm = TRUE),
  mean_meth_diff = mean(meth_diff, na.rm = TRUE)
), by = .(te_class, te_family)][n_total >= 10][order(te_class, -n_total)]

fwrite(te_family_summary, file.path(OUTPUT_DIR, "TE_family_change_summary.txt"), sep = "\t")

# Save te_methylation for plotting
saveRDS(te_methylation, file.path(OUTPUT_DIR, "te_methylation_summary.rds"))
saveRDS(dmp_te_df, file.path(OUTPUT_DIR, "dmp_te_df.rds"))
saveRDS(dmr_te_df, file.path(OUTPUT_DIR, "dmr_te_df.rds"))

# =============================================================================
# SECTION 8: VISUALIZATIONS
# =============================================================================

cat("\n=== Generating Visualizations ===\n\n")

save_plot <- function(plot_obj, name, width, height) {
  ggsave(file.path(OUTPUT_DIR, paste0(name, ".pdf")), plot_obj, 
         width = width, height = height)
  ggsave(file.path(OUTPUT_DIR, paste0(name, ".png")), plot_obj, 
         width = width, height = height, dpi = 300)
  cat(sprintf("  Saved: %s\n", name))
}

# Convert to data.frame for ggplot
te_meth_df <- as.data.frame(te_methylation)
dmp_te_plot <- as.data.frame(dmp_te_df)
dmr_te_plot <- as.data.frame(dmr_te_df)

# =========================================================================
# 8.1 RIDGE PLOTS - ALL TEs
# =========================================================================

cat("\n--- Ridge Plots: All TEs ---\n")

# Overall
te_meth_overall <- te_meth_df %>%
  filter(n_cpgs > 0, !is.na(mean_meth_control)) %>%
  select(mean_meth_control, mean_meth_treatment) %>%
  pivot_longer(cols = everything(),
               names_to = "condition", values_to = "methylation") %>%
  mutate(condition = ifelse(condition == "mean_meth_control", "Control", "Amputated"),
         category = "All TEs")

p_ridge_overall <- ggplot(te_meth_overall, 
                          aes(x = methylation * 100, y = category, fill = condition)) +
  geom_density_ridges(alpha = 0.6, scale = 0.9) +
  scale_fill_manual(values = c("Control" = "#2E86C1", "Amputated" = "#E74C3C")) +
  labs(title = "Overall TE Methylation Distribution",
       subtitle = sprintf("n = %s TEs with CpG sites", 
                          format(sum(te_meth_df$n_cpgs > 0), big.mark = ",")),
       x = "Mean Methylation (%)", y = "", fill = "Condition") +
  theme_minimal() +
  theme(axis.text.y = element_blank())
save_plot(p_ridge_overall, "ridge_methylation_ALL_TEs", 10, 4)

# By CLASS
te_meth_by_class <- te_meth_df %>%
  filter(n_cpgs > 0, !is.na(mean_meth_control)) %>%
  select(te_class, mean_meth_control, mean_meth_treatment) %>%
  pivot_longer(cols = c(mean_meth_control, mean_meth_treatment),
               names_to = "condition", values_to = "methylation") %>%
  mutate(condition = ifelse(condition == "mean_meth_control", "Control", "Amputated"))

p_ridge_class <- ggplot(te_meth_by_class, 
                        aes(x = methylation * 100, y = te_class, fill = condition)) +
  geom_density_ridges(alpha = 0.6, scale = 0.9) +
  scale_fill_manual(values = c("Control" = "#2E86C1", "Amputated" = "#E74C3C")) +
  labs(title = "TE Methylation Distribution by Class",
       x = "Mean Methylation (%)", y = "TE Class", fill = "Condition") +
  theme_minimal()
save_plot(p_ridge_class, "ridge_methylation_by_TE_class", 12, 8)

# By FAMILY (top 20)
top_families <- te_meth_df %>%
  filter(n_cpgs > 0) %>%
  count(te_class, te_family) %>%
  arrange(desc(n)) %>%
  head(20) %>%
  mutate(family_label = paste0(te_family, " (", te_class, ")"))

te_meth_by_family <- te_meth_df %>%
  filter(n_cpgs > 0, !is.na(mean_meth_control)) %>%
  inner_join(top_families, by = c("te_class", "te_family")) %>%
  select(family_label, mean_meth_control, mean_meth_treatment) %>%
  pivot_longer(cols = c(mean_meth_control, mean_meth_treatment),
               names_to = "condition", values_to = "methylation") %>%
  mutate(condition = ifelse(condition == "mean_meth_control", "Control", "Amputated"))

p_ridge_family <- ggplot(te_meth_by_family, 
                         aes(x = methylation * 100, y = family_label, fill = condition)) +
  geom_density_ridges(alpha = 0.6, scale = 0.9) +
  scale_fill_manual(values = c("Control" = "#2E86C1", "Amputated" = "#E74C3C")) +
  labs(title = "TE Methylation by Family (Top 20)",
       x = "Mean Methylation (%)", y = "TE Family", fill = "Condition") +
  theme_minimal()
save_plot(p_ridge_family, "ridge_methylation_by_TE_family", 14, 12)

# =========================================================================
# 8.2 RIDGE PLOTS - DMPs
# =========================================================================

cat("\n--- Ridge Plots: DMPs in TEs ---\n")

dmp_meth_data <- dmp_te_plot %>%
  mutate(meth_diff_pct = diff * 100)

# Overall
p_ridge_dmp_overall <- ggplot(dmp_meth_data, 
                              aes(x = meth_diff_pct, y = "All TEs", fill = direction)) +
  geom_density_ridges(alpha = 0.6, scale = 0.9) +
  geom_vline(xintercept = c(-10, 0, 10), linetype = c("dashed", "solid", "dashed"), 
             color = c("gray50", "black", "gray50")) +
  scale_fill_manual(values = c("Hyper" = "#E74C3C", "Hypo" = "#3498DB")) +
  labs(title = "DMP Methylation Difference in TEs",
       x = "Methylation Difference (%)", y = "", fill = "Direction") +
  theme_minimal()
save_plot(p_ridge_dmp_overall, "ridge_DMP_meth_diff_ALL_TEs", 10, 4)

# By class
p_ridge_dmp_class <- ggplot(dmp_meth_data, 
                            aes(x = meth_diff_pct, y = te_class, fill = direction)) +
  geom_density_ridges(alpha = 0.6, scale = 0.9) +
  geom_vline(xintercept = c(-10, 0, 10), linetype = c("dashed", "solid", "dashed"),
             color = c("gray50", "black", "gray50")) +
  scale_fill_manual(values = c("Hyper" = "#E74C3C", "Hypo" = "#3498DB")) +
  labs(title = "DMP Methylation Difference by TE Class",
       x = "Methylation Difference (%)", y = "TE Class", fill = "Direction") +
  theme_minimal()
save_plot(p_ridge_dmp_class, "ridge_DMP_meth_diff_by_TE_class", 12, 8)

# By family
dmp_family_counts <- dmp_meth_data %>%
  count(te_class, te_family) %>%
  arrange(desc(n)) %>%
  head(20) %>%
  mutate(family_label = paste0(te_family, " (", te_class, ")"))

dmp_meth_by_family <- dmp_meth_data %>%
  inner_join(dmp_family_counts, by = c("te_class", "te_family"))

p_ridge_dmp_family <- ggplot(dmp_meth_by_family, 
                             aes(x = meth_diff_pct, y = family_label, fill = direction)) +
  geom_density_ridges(alpha = 0.6, scale = 0.9) +
  geom_vline(xintercept = c(-10, 0, 10), linetype = c("dashed", "solid", "dashed"),
             color = c("gray50", "black", "gray50")) +
  scale_fill_manual(values = c("Hyper" = "#E74C3C", "Hypo" = "#3498DB")) +
  labs(title = "DMP Methylation Difference by TE Family (Top 20)",
       x = "Methylation Difference (%)", y = "TE Family", fill = "Direction") +
  theme_minimal()
save_plot(p_ridge_dmp_family, "ridge_DMP_meth_diff_by_TE_family", 14, 12)

# Significant only
dmp_sig_data <- dmp_meth_data %>% filter(significant)
if (nrow(dmp_sig_data) > 100) {
  p_ridge_dmp_sig <- ggplot(dmp_sig_data, 
                            aes(x = meth_diff_pct, y = te_class, fill = direction)) +
    geom_density_ridges(alpha = 0.6, scale = 0.9) +
    geom_vline(xintercept = c(-10, 0, 10), linetype = c("dashed", "solid", "dashed"),
               color = c("gray50", "black", "gray50")) +
    scale_fill_manual(values = c("Hyper" = "#E74C3C", "Hypo" = "#3498DB")) +
    labs(title = "Significant DMP Methylation Difference by TE Class",
         x = "Methylation Difference (%)", y = "TE Class", fill = "Direction") +
    theme_minimal()
  save_plot(p_ridge_dmp_sig, "ridge_DMP_significant_by_TE_class", 12, 8)
}

# =========================================================================
# 8.3 RIDGE PLOTS - DMRs
# =========================================================================

cat("\n--- Ridge Plots: DMRs in TEs ---\n")

dmr_meth_data <- dmr_te_plot %>%
  mutate(meth_diff_pct = diff_methy * 100)

# Overall
p_ridge_dmr_overall <- ggplot(dmr_meth_data, 
                              aes(x = meth_diff_pct, y = "All TEs", fill = direction)) +
  geom_density_ridges(alpha = 0.6, scale = 0.9) +
  geom_vline(xintercept = c(-10, 0, 10), linetype = c("dashed", "solid", "dashed"),
             color = c("gray50", "black", "gray50")) +
  scale_fill_manual(values = c("Hyper" = "#E74C3C", "Hypo" = "#3498DB")) +
  labs(title = "DMR Methylation Difference in TEs",
       x = "Methylation Difference (%)", y = "", fill = "Direction") +
  theme_minimal()
save_plot(p_ridge_dmr_overall, "ridge_DMR_meth_diff_ALL_TEs", 10, 4)

# By class
p_ridge_dmr_class <- ggplot(dmr_meth_data, 
                            aes(x = meth_diff_pct, y = te_class, fill = direction)) +
  geom_density_ridges(alpha = 0.6, scale = 0.9) +
  geom_vline(xintercept = c(-10, 0, 10), linetype = c("dashed", "solid", "dashed"),
             color = c("gray50", "black", "gray50")) +
  scale_fill_manual(values = c("Hyper" = "#E74C3C", "Hypo" = "#3498DB")) +
  labs(title = "DMR Methylation Difference by TE Class",
       x = "Methylation Difference (%)", y = "TE Class", fill = "Direction") +
  theme_minimal()
save_plot(p_ridge_dmr_class, "ridge_DMR_meth_diff_by_TE_class", 12, 8)

# By family
dmr_family_counts <- dmr_meth_data %>%
  count(te_class, te_family) %>%
  arrange(desc(n)) %>%
  head(20) %>%
  mutate(family_label = paste0(te_family, " (", te_class, ")"))

dmr_meth_by_family <- dmr_meth_data %>%
  inner_join(dmr_family_counts, by = c("te_class", "te_family"))

if (nrow(dmr_meth_by_family) > 50) {
  p_ridge_dmr_family <- ggplot(dmr_meth_by_family, 
                               aes(x = meth_diff_pct, y = family_label, fill = direction)) +
    geom_density_ridges(alpha = 0.6, scale = 0.9) +
    geom_vline(xintercept = c(-10, 0, 10), linetype = c("dashed", "solid", "dashed"),
               color = c("gray50", "black", "gray50")) +
    scale_fill_manual(values = c("Hyper" = "#E74C3C", "Hypo" = "#3498DB")) +
    labs(title = "DMR Methylation Difference by TE Family (Top 20)",
         x = "Methylation Difference (%)", y = "TE Family", fill = "Direction") +
    theme_minimal()
  save_plot(p_ridge_dmr_family, "ridge_DMR_meth_diff_by_TE_family", 14, 12)
}

# =========================================================================
# 8.4 BAR PLOTS
# =========================================================================

cat("\n--- Bar Plots ---\n")

# DMR count by class
dmr_class_counts <- dmr_meth_data %>%
  group_by(te_class, direction) %>%
  summarise(n = n(), .groups = "drop")

p_dmr_class <- ggplot(dmr_class_counts, 
                      aes(x = reorder(te_class, n), y = n, fill = direction)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Hyper" = "#E74C3C", "Hypo" = "#3498DB")) +
  coord_flip() +
  labs(title = "DMRs Overlapping TEs by Class",
       x = "TE Class", y = "Number of DMR overlaps", fill = "Direction") +
  theme_minimal()
save_plot(p_dmr_class, "barplot_DMR_count_by_TE_class", 10, 6)

# DMP count by class
dmp_class_counts <- dmp_meth_data %>%
  group_by(te_class, direction) %>%
  summarise(n = n(), .groups = "drop")

p_dmp_class <- ggplot(dmp_class_counts, 
                      aes(x = reorder(te_class, n), y = n, fill = direction)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Hyper" = "#E74C3C", "Hypo" = "#3498DB")) +
  coord_flip() +
  labs(title = "DMPs in TEs by Class",
       x = "TE Class", y = "Number of DMPs", fill = "Direction") +
  theme_minimal()
save_plot(p_dmp_class, "barplot_DMP_count_by_TE_class", 10, 6)

# TE change percentages
te_changes_long <- as.data.frame(te_class_summary) %>%
  select(te_class, pct_with_dmp, pct_with_sig_dmp, pct_with_dmr) %>%
  pivot_longer(cols = starts_with("pct"),
               names_to = "change_type", values_to = "percentage") %>%
  mutate(change_type = case_when(
    change_type == "pct_with_dmp" ~ "Any DMP",
    change_type == "pct_with_sig_dmp" ~ "Significant DMP",
    change_type == "pct_with_dmr" ~ "DMR overlap"
  ))

p_te_changes <- ggplot(te_changes_long, 
                       aes(x = reorder(te_class, percentage), y = percentage, fill = change_type)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Any DMP" = "#95A5A6", 
                               "Significant DMP" = "#E74C3C",
                               "DMR overlap" = "#3498DB")) +
  coord_flip() +
  labs(title = "Percentage of TEs with Methylation Changes",
       x = "TE Class", y = "Percentage (%)", fill = "Change Type") +
  theme_minimal()
save_plot(p_te_changes, "barplot_TE_with_changes_by_class", 12, 8)

# =========================================================================
# 8.5 VIOLIN PLOTS
# =========================================================================

cat("\n--- Violin Plots ---\n")

p_dmp_violin <- ggplot(dmp_meth_data, 
                       aes(x = te_class, y = meth_diff_pct, fill = te_class)) +
  geom_violin(alpha = 0.7) +
  geom_hline(yintercept = c(-10, 0, 10), linetype = c("dashed", "solid", "dashed"),
             color = c("gray50", "black", "gray50")) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Methylation Difference at DMPs in TEs",
       x = "TE Class", y = "Methylation Difference (%)") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p_dmp_violin, "violin_meth_diff_DMPs_in_TEs", 10, 6)

p_dmr_violin <- ggplot(dmr_meth_data, 
                       aes(x = te_class, y = meth_diff_pct, fill = te_class)) +
  geom_violin(alpha = 0.7) +
  geom_hline(yintercept = c(-10, 0, 10), linetype = c("dashed", "solid", "dashed"),
             color = c("gray50", "black", "gray50")) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Methylation Difference at DMRs in TEs",
       x = "TE Class", y = "Methylation Difference (%)") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p_dmr_violin, "violin_meth_diff_DMRs_in_TEs", 10, 6)

# =========================================================================
# 8.6 KIMURA (TE AGE) ANALYSIS
# =========================================================================

cat("\n--- Kimura (TE Age) Analysis ---\n")

te_age_meth <- te_meth_df %>%
  filter(n_cpgs > 0, !is.na(kimura_div), kimura_div < 60, !is.na(mean_meth_all))

p_kimura <- ggplot(te_age_meth, 
                   aes(x = kimura_div, y = mean_meth_all * 100, color = te_class)) +
  geom_point(alpha = 0.1, size = 0.3) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1) +
  facet_wrap(~ te_class) +
  labs(title = "TE Age (Kimura) vs. Methylation",
       x = "Kimura Divergence (%)", y = "Mean Methylation (%)") +
  theme_minimal() +
  theme(legend.position = "none")
save_plot(p_kimura, "scatter_kimura_vs_methylation", 14, 10)

# =========================================================================
# 8.7 WITH VS WITHOUT CHANGES COMPARISON
# =========================================================================

cat("\n--- Comparison Plots ---\n")

te_comparison <- te_meth_df %>%
  filter(n_cpgs > 0) %>%
  mutate(has_any_change = has_significant_dmp | has_dmr) %>%
  group_by(te_class, has_any_change) %>%
  summarise(
    n = n(),
    mean_meth = mean(mean_meth_all, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  mutate(status = ifelse(has_any_change, "With changes", "No changes"))

p_comparison <- ggplot(te_comparison, 
                       aes(x = te_class, y = mean_meth, fill = status)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("With changes" = "#E74C3C", "No changes" = "#95A5A6")) +
  labs(title = "Mean Methylation: TEs with vs without Changes",
       x = "TE Class", y = "Mean Methylation (%)", fill = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p_comparison, "barplot_methylation_with_vs_without_changes", 12, 6)

# =============================================================================
# FINAL SUMMARY
# =============================================================================

cat("\n\n")
cat("================================================================\n")
cat("           TE METHYLATION ANALYSIS COMPLETE                    \n")
cat("================================================================\n\n")

cat("KEY FINDINGS:\n\n")

cat("1. Overall TE Methylation:\n")
cat(sprintf("   - TEs with CpG sites: %s (%.1f%%)\n",
            format(sum(te_meth_df$n_cpgs > 0), big.mark = ","),
            100 * sum(te_meth_df$n_cpgs > 0) / nrow(te_meth_df)))
cat(sprintf("   - Mean methylation (Control): %.1f%%\n",
            mean(te_meth_df$mean_meth_control[te_meth_df$n_cpgs > 0], na.rm = TRUE) * 100))
cat(sprintf("   - Mean methylation (Amputated): %.1f%%\n",
            mean(te_meth_df$mean_meth_treatment[te_meth_df$n_cpgs > 0], na.rm = TRUE) * 100))

cat("\n2. DMPs in TEs:\n")
cat(sprintf("   - DMPs overlapping TEs: %s / %s (%.1f%%)\n",
            format(length(unique(dmp_te_df$dmp_idx)), big.mark = ","),
            format(length(dmp_gr), big.mark = ","),
            100 * length(unique(dmp_te_df$dmp_idx)) / length(dmp_gr)))
cat(sprintf("   - TEs with any DMP: %s\n", sum(te_meth_df$has_dmp)))
cat(sprintf("   - TEs with significant DMP: %s\n", sum(te_meth_df$has_significant_dmp)))

cat("\n3. DMRs in TEs:\n")
cat(sprintf("   - DMRs overlapping TEs: %s / %s (%.1f%%)\n",
            format(length(unique(dmr_te_df$dmr_idx)), big.mark = ","),
            format(length(dmr_gr), big.mark = ","),
            100 * length(unique(dmr_te_df$dmr_idx)) / length(dmr_gr)))
cat(sprintf("   - TEs with DMR overlap: %s\n", sum(te_meth_df$has_dmr)))

cat("\n\nOutput directory:", OUTPUT_DIR, "\n")
cat("================================================================\n")