#!/usr/bin/env Rscript
# =============================================================================
# Part 2: Differential Methylation Analysis & Advanced Visualization
# =============================================================================
# Prerequisites: Run Part 1 first (produces bsseq_object.rds)
# =============================================================================

suppressPackageStartupMessages({
  library(bsseq)
  library(DSS)
  library(GenomicRanges)
  library(data.table)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(dplyr)
  library(scales)
})

options(stringsAsFactors = FALSE)
options(scipen = 999)

theme_set(theme_minimal() +
            theme(text = element_text(size = 12),
                  plot.title = element_text(size = 14, face = "bold"),
                  axis.title = element_text(size = 12),
                  legend.text = element_text(size = 10)))

# >>> EDIT THIS: Same base directory as Part 1 <<<
BASE_DIR    <- "/mnt/data/alfredvar/"
RESULTS_DIR <- "/mnt/data/alfredvar/rlopezt/feb_w2_scripts/Preliminary/"


# Analysis parameters (same as Part 1)
analysis_params <- list(
  min_coverage   = 5,
  min_diff       = 0.1,
  fdr_threshold  = 0.05,
  min_dmr_cpgs   = 3,
  min_dmr_length = 50
)

# ── 1. Load BSseq Object ────────────────────────────────────────────────────

cat("Loading BSseq object from Part 1...\n")
bs_obj <- readRDS(file.path(RESULTS_DIR, "bsseq_object.rds"))
cat(sprintf("  Loaded: %s CpG sites x %s samples\n",
            format(nrow(bs_obj), big.mark = ","), ncol(bs_obj)))
print(pData(bs_obj)[, c("sample_id", "condition")])

# ── 2. Differential Methylation Analysis with DSS ────────────────────────────

cat("\n--- Differential Methylation Testing (DSS) ---\n")

group1_samples <- c("C1", "C2")     # Control
group2_samples <- c("A1", "A2")     # Treatment

cat("  Group 1 (Control):", paste(group1_samples, collapse = ", "), "\n")
cat("  Group 2 (Treatment):", paste(group2_samples, collapse = ", "), "\n")

# DML test with smoothing
cat("  Running DMLtest (smoothing = TRUE)... this may take a while.\n")
dml_test <- DMLtest(bs_obj,
                    group1 = group1_samples,
                    group2 = group2_samples,
                    smoothing = TRUE)

cat(sprintf("  DMLtest complete: %s CpG sites tested\n",
            format(nrow(dml_test), big.mark = ",")))

# Call DMPs (differentially methylated positions)
dmp_results <- callDML(dml_test,
                       p.threshold = analysis_params$fdr_threshold,
                       delta = analysis_params$min_diff)

cat(sprintf("  Significant DMPs: %s\n", format(nrow(dmp_results), big.mark = ",")))

# Call DMRs (differentially methylated regions)
dmr_results <- callDMR(dml_test,
                       p.threshold = analysis_params$fdr_threshold,
                       delta = analysis_params$min_diff,
                       minlen = analysis_params$min_dmr_length,
                       minCG = analysis_params$min_dmr_cpgs)

cat(sprintf("  Significant DMRs: %s\n", format(nrow(dmr_results), big.mark = ",")))

# ── 3. Save Results & Summary ────────────────────────────────────────────────

if (nrow(dmp_results) > 0) {
  write.table(dmp_results, file.path(RESULTS_DIR, "differentially_methylated_positions.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

if (nrow(dmr_results) > 0) {
  write.table(dmr_results, file.path(RESULTS_DIR, "differentially_methylated_regions.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

# Summary statistics
n_hyper_dmp <- ifelse(nrow(dmp_results) > 0,
                      sum(dmp_results$fdr < analysis_params$fdr_threshold &
                            dmp_results$diff > analysis_params$min_diff, na.rm = TRUE), 0)
n_hypo_dmp  <- ifelse(nrow(dmp_results) > 0,
                      sum(dmp_results$fdr < analysis_params$fdr_threshold &
                            dmp_results$diff < -analysis_params$min_diff, na.rm = TRUE), 0)

summary_stats <- data.frame(
  Metric = c("Total CpG sites analyzed",
             "Significant DMPs (FDR < 0.05)",
             "Significant DMRs",
             "Hypermethylated sites (Treatment > Control)",
             "Hypomethylated sites (Treatment < Control)",
             "Mean DMR length (bp)",
             "Mean CpGs per DMR"),
  Value = c(
    format(nrow(dml_test), big.mark = ","),
    format(nrow(dmp_results), big.mark = ","),
    format(nrow(dmr_results), big.mark = ","),
    format(n_hyper_dmp, big.mark = ","),
    format(n_hypo_dmp, big.mark = ","),
    ifelse(nrow(dmr_results) > 0, round(mean(dmr_results$length), 1), "NA"),
    ifelse(nrow(dmr_results) > 0, round(mean(dmr_results$nCG), 1), "NA")
  )
)

write.table(summary_stats, file.path(RESULTS_DIR, "analysis_summary.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nAnalysis Summary:\n")
print(summary_stats, row.names = FALSE)

# Save DML test for later use
saveRDS(dml_test, file.path(RESULTS_DIR, "dml_test_results.rds"))
saveRDS(dmp_results, file.path(RESULTS_DIR, "dmp_results.rds"))
saveRDS(dmr_results, file.path(RESULTS_DIR, "dmr_results.rds"))

# ── 4. Volcano Plot ──────────────────────────────────────────────────────────

cat("\n--- Generating Visualizations ---\n")

if (nrow(dmp_results) > 0) {
  volcano_data <- data.frame(
    diff   = dmp_results$diff,
    pvalue = dmp_results$pval,
    fdr    = dmp_results$fdr,
    significant = dmp_results$fdr < analysis_params$fdr_threshold &
      abs(dmp_results$diff) > analysis_params$min_diff
  )
  
  volcano_data$category <- "Not Significant"
  volcano_data$category[volcano_data$significant & volcano_data$diff > 0] <- "Hypermethylated"
  volcano_data$category[volcano_data$significant & volcano_data$diff < 0] <- "Hypomethylated"
  volcano_data$category <- factor(volcano_data$category,
                                  levels = c("Not Significant",
                                             "Hypermethylated",
                                             "Hypomethylated"))
  
  n_hyper <- sum(volcano_data$category == "Hypermethylated")
  n_hypo  <- sum(volcano_data$category == "Hypomethylated")
  
  p_volcano <- ggplot(volcano_data, aes(x = diff, y = -log10(pvalue), color = category)) +
    geom_point(alpha = 0.6, size = 0.8) +
    scale_color_manual(
      values = c("Not Significant" = "gray70",
                 "Hypermethylated" = "#E74C3C",
                 "Hypomethylated"  = "#3498DB"),
      name = "Methylation Change"
    ) +
    labs(
      title    = "Volcano Plot - Differentially Methylated Positions",
      subtitle = sprintf("Hypermethylated: %s | Hypomethylated: %s",
                         format(n_hyper, big.mark = ","),
                         format(n_hypo, big.mark = ",")),
      x = "Methylation Difference (Treatment - Control)",
      y = "-log10(p-value)"
    ) +
    geom_vline(xintercept = c(-analysis_params$min_diff, analysis_params$min_diff),
               linetype = "dashed", alpha = 0.7, color = "gray50") +
    geom_hline(yintercept = -log10(0.05),
               linetype = "dashed", alpha = 0.7, color = "gray50") +
    theme_minimal()
  
  ggsave(file.path(RESULTS_DIR, "volcano_plot_dmps.pdf"), p_volcano,
         width = 10, height = 8)
  cat("  Volcano plot saved.\n")
}

# ── 5. Manhattan Plot ────────────────────────────────────────────────────────

if (nrow(dmp_results) > 0) {
  # Extract chromosome numbers from your naming convention
  # Adjust the gsub pattern if your chroms are named differently
  manhattan_data <- dmp_results %>%
    mutate(
      chr_num = as.numeric(gsub("chr", "", chr)),
      significance = ifelse(
        fdr < analysis_params$fdr_threshold & abs(diff) > analysis_params$min_diff,
        "Significant", "Not Significant"
      )
    ) %>%
    filter(!is.na(chr_num)) %>%
    arrange(chr_num, pos)
  
  chr_info <- manhattan_data %>%
    group_by(chr, chr_num) %>%
    summarise(max_pos = max(pos), .groups = "drop") %>%
    arrange(chr_num) %>%
    mutate(
      cumulative_length = cumsum(c(0, head(max_pos, -1))),
      chr_center = cumulative_length + max_pos / 2
    )
  
  manhattan_plot_data <- manhattan_data %>%
    left_join(chr_info, by = c("chr", "chr_num")) %>%
    mutate(genome_pos = cumulative_length + pos)
  
  p_manhattan <- ggplot(manhattan_plot_data,
                        aes(x = genome_pos, y = diff, color = significance)) +
    geom_point(alpha = 0.6, size = 0.5) +
    scale_color_manual(
      values = c("Not Significant" = "gray70", "Significant" = "#E74C3C"),
      name = "Significance"
    ) +
    scale_x_continuous(
      breaks = chr_info$chr_center,
      labels = gsub("chr", "", chr_info$chr),
      expand = c(0.01, 0)
    ) +
    labs(
      title    = "Manhattan Plot - Methylation Differences Across the Genome",
      subtitle = "Treatment vs Control Comparison",
      x = "Chromosome",
      y = "Methylation Difference (Treatment - Control)"
    ) +
    theme_minimal() +
    theme(
      legend.position   = "bottom",
      panel.grid.minor  = element_blank(),
      panel.grid.major.x = element_blank()
    ) +
    geom_hline(yintercept = c(-analysis_params$min_diff, analysis_params$min_diff),
               linetype = "dashed", alpha = 0.7, color = "gray50") +
    geom_vline(xintercept = chr_info$cumulative_length[-1],
               color = "gray80", linetype = "dashed", alpha = 0.3)
  
  ggsave(file.path(RESULTS_DIR, "manhattan_plot_methylation_differences.pdf"),
         p_manhattan, width = 16, height = 8)
  cat("  Manhattan plot saved.\n")
}

# ── 6. Genome-Wide Methylation Profile ───────────────────────────────────────

plot_genome_wide_methylation <- function(bs_obj, sample_name, window_size = 1000000) {
  sample_idx <- which(colnames(bs_obj) == sample_name)
  beta_vals  <- getMeth(bs_obj, type = "raw")[, sample_idx]
  
  methylation_data <- data.frame(
    chr         = as.character(seqnames(bs_obj)),
    pos         = start(bs_obj),
    methylation = beta_vals * 100,
    stringsAsFactors = FALSE
  )
  
  methylation_data <- methylation_data[!is.na(methylation_data$methylation), ]
  
  # Order chromosomes numerically
  chr_order <- paste0("chr", 1:31)
  methylation_data <- methylation_data[methylation_data$chr %in% chr_order, ]
  methylation_data$chr <- factor(methylation_data$chr, levels = chr_order)
  
  methylation_binned <- methylation_data %>%
    mutate(
      bin     = floor(pos / window_size) * window_size,
      chr_num = as.numeric(gsub("chr", "", as.character(chr)))
    ) %>%
    group_by(chr, chr_num, bin) %>%
    summarise(
      mean_methylation = mean(methylation, na.rm = TRUE),
      n_cpgs = n(),
      .groups = "drop"
    ) %>%
    filter(n_cpgs >= 10)
  
  chr_lengths <- methylation_binned %>%
    group_by(chr, chr_num) %>%
    summarise(max_pos = max(bin), .groups = "drop") %>%
    arrange(chr_num) %>%
    mutate(
      cumulative_length = cumsum(c(0, head(max_pos, -1))),
      chr_center = cumulative_length + max_pos / 2
    )
  
  methylation_plot <- methylation_binned %>%
    left_join(chr_lengths, by = c("chr", "chr_num")) %>%
    mutate(genome_pos = cumulative_length + bin)
  
  ggplot(methylation_plot, aes(x = genome_pos, y = mean_methylation, color = chr)) +
    geom_point(alpha = 0.6, size = 0.5) +
    scale_color_manual(values = rep(c("#2E86C1", "#E74C3C"), 16), guide = "none") +
    scale_x_continuous(
      breaks = chr_lengths$chr_center,
      labels = gsub("chr", "", chr_lengths$chr),
      expand = c(0.01, 0)
    ) +
    scale_y_continuous(limits = c(0, 100), expand = c(0.02, 0)) +
    labs(
      title    = paste("Whole-Genome DNA Methylation Profile:", sample_name),
      subtitle = paste("Window size:", window_size / 1e6, "Mb"),
      x = "Chromosome",
      y = "Mean Methylation (%)"
    ) +
    theme_minimal() +
    theme(
      legend.position    = "none",
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      plot.title = element_text(size = 14, face = "bold")
    ) +
    geom_vline(xintercept = chr_lengths$cumulative_length[-1],
               color = "gray80", linetype = "dashed", alpha = 0.7)
}

pdf(file.path(RESULTS_DIR, "genome_wide_methylation.pdf"), width = 16, height = 10)
for (sample in colnames(bs_obj)) {
  p_genome <- plot_genome_wide_methylation(bs_obj, sample, window_size = 1000000)
  print(p_genome)
}
dev.off()
cat("  Genome-wide methylation profiles saved.\n")

# ── 7. DMP & DMR Heatmaps ───────────────────────────────────────────────────

# DMP-level heatmap (top 100 by effect size)
create_dmp_level_heatmap <- function(bs_obj, dmp_results, n_top = 100) {
  if (nrow(dmp_results) == 0) { cat("No DMPs - skipping heatmap.\n"); return() }
  
  n_sites  <- min(n_top, nrow(dmp_results))
  top_dmps <- dmp_results[order(abs(dmp_results$diff), decreasing = TRUE), ][1:n_sites, ]
  
  dmp_ids    <- paste(top_dmps$chr, top_dmps$pos, sep = ":")
  bs_site_ids <- paste(as.character(seqnames(bs_obj)), start(bs_obj), sep = ":")
  matching_idx <- which(bs_site_ids %in% dmp_ids)
  
  if (length(matching_idx) < 2) { cat("Not enough matching DMPs.\n"); return() }
  
  beta_mat <- getMeth(bs_obj[matching_idx, ], type = "raw")
  rownames(beta_mat) <- paste0("DMP_", seq_len(nrow(beta_mat)))
  
  # Filter rows with zero variance
  rv <- apply(beta_mat, 1, var, na.rm = TRUE)
  valid <- complete.cases(beta_mat) & rv > 0 & !is.na(rv)
  beta_mat <- beta_mat[valid, ]
  
  if (nrow(beta_mat) < 2) { cat("Not enough valid DMPs for heatmap.\n"); return() }
  
  col_ann <- data.frame(Condition = pData(bs_obj)$condition,
                        row.names = colnames(beta_mat))
  col_colors <- list(Condition = c("Control" = "#2E86C1", "Treatment" = "#E74C3C"))
  
  pdf(file.path(RESULTS_DIR, "dmp_level_heatmap.pdf"), width = 8, height = 12)
  pheatmap(beta_mat,
           scale = "row",
           clustering_distance_rows = "euclidean",
           clustering_distance_cols = "euclidean",
           clustering_method = "complete",
           annotation_col = col_ann,
           annotation_colors = col_colors,
           color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
           main = sprintf("Top %d DMPs - Site-Level Methylation", nrow(beta_mat)),
           fontsize = 8,
           show_rownames = FALSE,
           show_colnames = TRUE)
  dev.off()
  cat("  DMP heatmap saved.\n")
}

# DMR-level heatmap (top 50 by effect size)
create_dmr_level_heatmap <- function(bs_obj, dmr_results, n_top = 50) {
  if (nrow(dmr_results) == 0) { cat("No DMRs - skipping heatmap.\n"); return() }
  
  n_regions <- min(n_top, nrow(dmr_results))
  top_dmrs  <- dmr_results[order(abs(dmr_results$diff.Methy), decreasing = TRUE), ][1:n_regions, ]
  
  dmr_matrix <- matrix(nrow = n_regions, ncol = ncol(bs_obj))
  colnames(dmr_matrix) <- colnames(bs_obj)
  rownames(dmr_matrix) <- paste0("DMR_", 1:n_regions, "_",
                                 top_dmrs$chr, ":",
                                 top_dmrs$start, "-",
                                 top_dmrs$end)
  
  bs_ranges <- GRanges(seqnames = seqnames(bs_obj),
                       ranges = IRanges(start = start(bs_obj), width = 1))
  
  valid_dmrs <- logical(n_regions)
  for (i in 1:n_regions) {
    dmr <- top_dmrs[i, ]
    dmr_range <- GRanges(seqnames = dmr$chr,
                         ranges = IRanges(start = dmr$start, end = dmr$end))
    overlaps    <- findOverlaps(bs_ranges, dmr_range)
    cpg_indices <- queryHits(overlaps)
    
    if (length(cpg_indices) > 0) {
      dmr_beta <- getMeth(bs_obj[cpg_indices, ], type = "raw")
      if (length(cpg_indices) == 1) {
        dmr_matrix[i, ] <- as.numeric(dmr_beta)
      } else {
        dmr_matrix[i, ] <- colMeans(dmr_beta, na.rm = TRUE)
      }
      valid_dmrs[i] <- TRUE
    }
  }
  
  dmr_matrix <- dmr_matrix[valid_dmrs, , drop = FALSE]
  if (nrow(dmr_matrix) < 2) { cat("Not enough valid DMRs for heatmap.\n"); return() }
  
  col_ann <- data.frame(Condition = pData(bs_obj)$condition,
                        row.names = colnames(dmr_matrix))
  col_colors <- list(Condition = c("Control" = "#2E86C1", "Treatment" = "#E74C3C"))
  
  pdf(file.path(RESULTS_DIR, "dmr_level_heatmap.pdf"), width = 10, height = 12)
  pheatmap(dmr_matrix,
           scale = "row",
           clustering_distance_rows = "euclidean",
           clustering_distance_cols = "euclidean",
           clustering_method = "complete",
           annotation_col = col_ann,
           annotation_colors = col_colors,
           color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
           main = sprintf("Top %d DMRs - Region-Level Methylation", nrow(dmr_matrix)),
           fontsize = 10,
           show_rownames = TRUE,
           show_colnames = TRUE)
  dev.off()
  cat("  DMR heatmap saved.\n")
}

if (nrow(dmp_results) > 0) create_dmp_level_heatmap(bs_obj, dmp_results, n_top = 100)
if (nrow(dmr_results) > 0) create_dmr_level_heatmap(bs_obj, dmr_results, n_top = 50)

cat("\n========================================\n")
cat("Part 2 complete! Outputs in:", RESULTS_DIR, "\n")
cat("Files created:\n")
cat("  - differentially_methylated_positions.txt\n")
cat("  - differentially_methylated_regions.txt\n")
cat("  - analysis_summary.txt\n")
cat("  - dml_test_results.rds / dmp_results.rds / dmr_results.rds\n")
cat("  - volcano_plot_dmps.pdf\n")
cat("  - manhattan_plot_methylation_differences.pdf\n")
cat("  - genome_wide_methylation.pdf\n")
cat("  - dmp_level_heatmap.pdf\n")
cat("  - dmr_level_heatmap.pdf\n")
cat("========================================\n")
cat("\nProceed to Part 3 for annotation & enrichment analysis.\n")