#!/usr/bin/env Rscript
# =============================================================================
# WGCNA Parameter Sensitivity Analysis
# =============================================================================
# Purpose: Systematically test different variance filters and soft threshold
#          powers to find the parameter set that best achieves scale-free
#          topology while maintaining biologically meaningful modules.
#
# Tests:
#   - 5 variance filter levels (top 3K, 5K, 8K, 10K, MAD-based)
#   - Full power sweep (1-30) for each filter level
#   - Automated outlier detection (Z.k threshold)
#   - Module stability comparison with original WGCNA (ARI)
#
# Requires: 64GB RAM, 8 CPUs, ~3-5 hours
# Input:  results/02_rnaseq/Part5_WGCNA/inputs/normalized_counts.tsv
#         results/02_rnaseq/Part5_WGCNA/data/all_gene_module_assignments.tsv
# Output: results/wgcna_sensitivity/
# =============================================================================

library(WGCNA)
library(ggplot2)
library(dplyr)

enableWGCNAThreads(nThreads = 8)
options(stringsAsFactors = FALSE)

cat("=== WGCNA Parameter Sensitivity Analysis ===\n")
cat("CPUs:", parallel::detectCores(), "| Using 8 cores\n\n")

outdir <- "results/wgcna_sensitivity"
dir.create(file.path(outdir, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "data"), recursive = TRUE, showWarnings = FALSE)

# ── Load expression data ────────────────────────────────────────────────────

cat("Loading expression data...\n")

# Waterfall: cluster paths first, then local paths
expr_data <- NULL
meta_data <- NULL

# Try 1: Cluster CWD (Part4 RData — standard input on cluster)
rdata_file <- "Part4_Kmeans/data/Part4_objects.RData"
if (file.exists(rdata_file)) {
  load(rdata_file)
  expr_data <- normalized_counts
  meta_data <- df
  cat("  Loaded from Part4 objects (cluster CWD)\n")
}

# Try 2: Cluster WGCNA saved TSV inputs
if (is.null(expr_data)) {
  for (wpath in c("Part5_WGCNA/inputs/normalized_counts.tsv",
                   "results/02_rnaseq/Part5_WGCNA/inputs/normalized_counts.tsv")) {
    if (file.exists(wpath)) {
      expr_data <- read.delim(wpath, check.names = FALSE)
      mpath <- sub("normalized_counts.tsv", "sample_metadata.tsv", wpath)
      if (file.exists(mpath)) meta_data <- read.delim(mpath, check.names = FALSE)
      cat("  Loaded from:", wpath, "\n")
      break
    }
  }
}

if (is.null(expr_data)) stop("Cannot find expression data. Check paths.")

cat("  Dimensions:", nrow(expr_data), "genes x", ncol(expr_data), "samples\n")

# Remove NAs
expr_data <- expr_data[complete.cases(expr_data), ]
cat("  After NA removal:", nrow(expr_data), "genes\n")

# Load original module assignments for comparison (waterfall)
orig_modules <- NULL
for (mpath in c("Part5_WGCNA/data/all_gene_module_assignments.tsv",
                "results/02_rnaseq/Part5_WGCNA/data/all_gene_module_assignments.tsv")) {
  if (file.exists(mpath)) {
    orig_modules <- read.delim(mpath)
    cat("  Original modules from:", mpath, "\n")
    break
  }
}
if (is.null(orig_modules)) stop("Cannot find original module assignments.")
cat("  Original modules:", length(unique(orig_modules$Module_Color)), "colors,",
    nrow(orig_modules), "genes\n\n")


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 1: Automated Outlier Detection
# ══════════════════════════════════════════════════════════════════════════════

cat("--- Section 1: Automated Outlier Detection ---\n")

# Use standardized connectivity (Z.k) to find outliers
# WGCNA FAQ recommends Z.k < -2.5 as threshold
expr_t <- t(expr_data)

# Quick adjacency for connectivity
A <- adjacency(expr_t, power = 6, type = "signed")
k <- as.numeric(apply(A, 2, sum)) - 1
Z_k <- scale(k)

outlier_threshold <- -2.5
auto_outliers <- rownames(expr_t)[Z_k < outlier_threshold]

cat("  Standardized connectivity Z.k threshold:", outlier_threshold, "\n")
cat("  Automated outliers detected:", length(auto_outliers), "\n")
if (length(auto_outliers) > 0) {
  for (s in auto_outliers) {
    cat("    ", s, "Z.k =", round(Z_k[which(rownames(expr_t) == s)], 2), "\n")
  }
}

# Compare with original hard-coded outliers
orig_outliers <- c("T1S5", "dcrep4", "irrep7")
match_orig <- intersect(auto_outliers, orig_outliers)
new_outliers <- setdiff(auto_outliers, orig_outliers)
missed_orig <- setdiff(orig_outliers, auto_outliers)

cat("\n  Comparison with original outlier selection:\n")
cat("    Confirmed:", paste(match_orig, collapse = ", "), "\n")
cat("    New (auto only):", ifelse(length(new_outliers) > 0,
                                    paste(new_outliers, collapse = ", "), "none"), "\n")
cat("    Missed (orig only):", ifelse(length(missed_orig) > 0,
                                       paste(missed_orig, collapse = ", "), "none"), "\n")

# Use original outliers for consistency (but report the comparison)
outliers_to_remove <- orig_outliers[orig_outliers %in% rownames(expr_t)]
expr_t_clean <- expr_t[!rownames(expr_t) %in% outliers_to_remove, ]
cat("\n  Final: removed", length(outliers_to_remove), "outliers →",
    nrow(expr_t_clean), "samples\n\n")

# Save outlier comparison
outlier_table <- data.frame(
  sample = rownames(expr_t),
  Z_k = as.numeric(Z_k),
  is_auto_outlier = rownames(expr_t) %in% auto_outliers,
  is_orig_outlier = rownames(expr_t) %in% orig_outliers,
  removed = rownames(expr_t) %in% outliers_to_remove
)
write.table(outlier_table, file.path(outdir, "tables", "S01_outlier_detection.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 2: Variance Filter Sweep
# ══════════════════════════════════════════════════════════════════════════════

cat("--- Section 2: Variance Filter Sweep ---\n")

gene_variance <- apply(expr_data, 1, var)
gene_mad <- apply(expr_data, 1, mad)

# Define 5 filter strategies
filter_configs <- list(
  list(name = "top3000", n = 3000, desc = "Top 3,000 by variance"),
  list(name = "top5000", n = 5000, desc = "Top 5,000 by variance"),
  list(name = "top8000", n = 8000, desc = "Top 8,000 by variance"),
  list(name = "top10000", n = 10000, desc = "Top 10,000 by variance"),
  list(name = "MAD75", n = NA, desc = "Top 75th percentile by MAD")
)

filtered_matrices <- list()
for (cfg in filter_configs) {
  if (cfg$name == "MAD75") {
    # MAD-based: keep genes above 25th percentile of MAD
    mad_threshold <- quantile(gene_mad, 0.25)
    keep_genes <- names(gene_mad)[gene_mad > mad_threshold]
  } else {
    # Top N by variance
    keep_genes <- names(sort(gene_variance, decreasing = TRUE))[1:min(cfg$n, length(gene_variance))]
  }

  mat <- expr_t_clean[, keep_genes]
  filtered_matrices[[cfg$name]] <- mat
  cat(sprintf("  %-10s: %d genes (from %s)\n", cfg$name, ncol(mat), cfg$desc))
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 3: Soft Threshold Power Sweep
# ══════════════════════════════════════════════════════════════════════════════

cat("\n--- Section 3: Soft Threshold Power Sweep ---\n")

powers <- c(seq(1, 10, by = 1), seq(12, 30, by = 2))

sft_results <- data.frame()
best_params <- data.frame()

for (cfg in filter_configs) {
  cat(sprintf("\n  Testing %s (%d genes)...\n", cfg$name, ncol(filtered_matrices[[cfg$name]])))

  mat <- filtered_matrices[[cfg$name]]

  # Check for good genes/samples
  gsg <- goodSamplesGenes(mat, verbose = 0)
  if (!gsg$allOK) {
    mat <- mat[gsg$goodSamples, gsg$goodGenes]
    cat("    Removed bad genes/samples:", ncol(filtered_matrices[[cfg$name]]) - ncol(mat), "genes\n")
    filtered_matrices[[cfg$name]] <- mat
  }

  sft <- pickSoftThreshold(mat, powerVector = powers,
                           networkType = "signed", verbose = 0)

  fit_vals <- -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2]

  for (i in seq_along(powers)) {
    sft_results <- rbind(sft_results, data.frame(
      filter = cfg$name,
      n_genes = ncol(mat),
      power = powers[i],
      R2 = round(fit_vals[i], 4),
      mean_k = round(sft$fitIndices[i, 5], 2),
      median_k = round(sft$fitIndices[i, 6], 2),
      max_k = round(sft$fitIndices[i, 7], 2)
    ))
  }

  # Find best power for this filter
  max_r2 <- max(fit_vals)
  best_power_r2 <- powers[which.max(fit_vals)]

  # Power that first reaches R2 > 0.80
  above_80 <- which(fit_vals > 0.80)
  if (length(above_80) > 0) {
    best_power_80 <- powers[above_80[1]]
    best_r2_80 <- fit_vals[above_80[1]]
    mean_k_80 <- sft$fitIndices[above_80[1], 5]
  } else {
    best_power_80 <- best_power_r2
    best_r2_80 <- max_r2
    mean_k_80 <- sft$fitIndices[which.max(fit_vals), 5]
  }

  best_params <- rbind(best_params, data.frame(
    filter = cfg$name,
    n_genes = ncol(mat),
    max_R2 = round(max_r2, 4),
    power_at_max_R2 = best_power_r2,
    power_at_R2_80 = best_power_80,
    R2_at_best = round(best_r2_80, 4),
    mean_k_at_best = round(mean_k_80, 2),
    reaches_80 = any(fit_vals > 0.80)
  ))

  cat(sprintf("    Max R²=%.4f at power=%d | R²>0.80: %s (power=%d, mean_k=%.1f)\n",
              max_r2, best_power_r2,
              ifelse(any(fit_vals > 0.80), "YES", "NO"),
              best_power_80, mean_k_80))
}

write.table(sft_results, file.path(outdir, "tables", "S02_soft_threshold_sweep.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(best_params, file.path(outdir, "tables", "S03_best_parameters.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n  Parameter summary:\n")
print(best_params)


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 4: Quick Module Detection per Filter Level
# ══════════════════════════════════════════════════════════════════════════════

cat("\n--- Section 4: Module Detection per Filter Level ---\n")

module_results <- list()

for (i in seq_along(filter_configs)) {
  cfg <- filter_configs[[i]]
  mat <- filtered_matrices[[cfg$name]]
  bp <- best_params[i, ]

  # Use the power that achieves R2>0.80, or best available
  use_power <- bp$power_at_R2_80

  cat(sprintf("\n  Running blockwiseModules for %s (n=%d, power=%d)...\n",
              cfg$name, ncol(mat), use_power))

  net <- blockwiseModules(
    mat,
    power = use_power,
    networkType = "signed",
    TOMType = "signed",
    corType = "bicor",
    maxPOutliers = 0.05,
    minModuleSize = 30,
    reassignThreshold = 0,
    mergeCutHeight = 0.25,
    numericLabels = TRUE,
    pamRespectsDendro = FALSE,
    saveTOMs = FALSE,
    verbose = 1
  )

  module_colors <- labels2colors(net$colors)
  n_modules <- length(unique(module_colors)) - 1  # exclude grey
  grey_pct <- round(sum(module_colors == "grey") / length(module_colors) * 100, 1)

  module_results[[cfg$name]] <- data.frame(
    gene = colnames(mat),
    module = module_colors,
    stringsAsFactors = FALSE
  )

  cat(sprintf("    Modules: %d | Grey: %.1f%% | Largest: %d genes\n",
              n_modules, grey_pct, max(table(module_colors[module_colors != "grey"]))))
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 5: Module Stability Analysis (ARI)
# ══════════════════════════════════════════════════════════════════════════════

cat("\n--- Section 5: Module Stability Analysis ---\n")

# Adjusted Rand Index: measures agreement between clusterings
# ARI = 1 means perfect agreement, 0 means random
compute_ari <- function(labels1, labels2) {
  # Build contingency table
  ct <- table(labels1, labels2)
  n <- sum(ct)
  a <- sum(choose(ct, 2))
  b <- sum(choose(rowSums(ct), 2))
  c_val <- sum(choose(colSums(ct), 2))
  d <- choose(n, 2)
  expected <- b * c_val / d
  max_val <- (b + c_val) / 2
  if (max_val == expected) return(1)
  (a - expected) / (max_val - expected)
}

stability_results <- data.frame()

for (cfg in filter_configs) {
  # Compare with original: find shared genes
  new_mods <- module_results[[cfg$name]]
  shared_genes <- intersect(new_mods$gene, orig_modules$Gene)

  if (length(shared_genes) >= 100) {
    orig_labels <- orig_modules$Module_Color[match(shared_genes, orig_modules$Gene)]
    new_labels <- new_mods$module[match(shared_genes, new_mods$gene)]

    ari <- compute_ari(orig_labels, new_labels)
    n_shared <- length(shared_genes)
  } else {
    ari <- NA
    n_shared <- length(shared_genes)
  }

  n_mods <- length(unique(new_mods$module)) - 1
  grey_pct <- round(sum(new_mods$module == "grey") / nrow(new_mods) * 100, 1)

  stability_results <- rbind(stability_results, data.frame(
    filter = cfg$name,
    n_genes = nrow(new_mods),
    n_modules = n_mods,
    grey_pct = grey_pct,
    n_shared_with_orig = n_shared,
    ARI_vs_original = round(ari, 4)
  ))

  cat(sprintf("  %-10s: %d modules, %.1f%% grey, ARI=%.4f (shared=%d genes)\n",
              cfg$name, n_mods, grey_pct, ifelse(is.na(ari), 0, ari), n_shared))
}

write.table(stability_results, file.path(outdir, "tables", "S04_module_stability.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 6: Parameter Recommendation
# ══════════════════════════════════════════════════════════════════════════════

cat("\n--- Section 6: Parameter Recommendation ---\n")

# Merge all results
recommendation <- best_params %>%
  left_join(stability_results, by = "filter")

# Score: prioritize R2 >= 0.80, then ARI, then low grey
recommendation$score <- 0
recommendation$score <- recommendation$score +
  ifelse(recommendation$reaches_80, 3, 0) +         # R2 > 0.80 gets 3 points
  ifelse(!is.na(recommendation$ARI_vs_original),
         recommendation$ARI_vs_original * 2, 0) +     # ARI * 2
  (1 - recommendation$grey_pct / 100) +               # Low grey gets up to 1
  ifelse(recommendation$mean_k_at_best > 10, 1, 0)    # Mean k > 10 gets 1

recommendation <- recommendation %>% arrange(desc(score))

cat("\nFinal recommendation table:\n")
print(recommendation %>% select(filter, n_genes.x, max_R2, power_at_R2_80,
                                 R2_at_best, mean_k_at_best, n_modules,
                                 grey_pct, ARI_vs_original, score))

best_filter <- recommendation$filter[1]
best_power_final <- recommendation$power_at_R2_80[1]

cat(sprintf("\n  RECOMMENDED: filter=%s, power=%d (R²=%.4f, ARI=%.4f, %d modules)\n",
            best_filter, best_power_final,
            recommendation$R2_at_best[1],
            recommendation$ARI_vs_original[1],
            recommendation$n_modules[1]))

write.table(recommendation, file.path(outdir, "tables", "S05_parameter_recommendation.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Save recommendation for Part B
save(best_filter, best_power_final, filtered_matrices, outliers_to_remove,
     module_results, stability_results, recommendation,
     file = file.path(outdir, "data", "parameter_sensitivity_results.RData"))


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 7: Diagnostic Plots
# ══════════════════════════════════════════════════════════════════════════════

cat("\n--- Section 7: Diagnostic Plots ---\n")

# Plot 1: Scale-free topology R² vs power, faceted by filter
p1 <- ggplot(sft_results, aes(x = power, y = R2, color = filter)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  geom_hline(yintercept = 0.80, linetype = "dashed", color = "red", linewidth = 0.5) +
  facet_wrap(~ filter, ncol = 3) +
  labs(title = "Scale-Free Topology Fit vs Soft Threshold Power",
       subtitle = "Red dashed line = R² = 0.80 target",
       x = "Soft Threshold Power", y = "Scale-Free Topology R²",
       color = "Filter") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(outdir, "plots", "S01_R2_vs_power.png"), p1, width = 12, height = 8, dpi = 300)
ggsave(file.path(outdir, "plots", "S01_R2_vs_power.pdf"), p1, width = 12, height = 8)
cat("  Saved: S01_R2_vs_power\n")

# Plot 2: Mean connectivity vs power
p2 <- ggplot(sft_results, aes(x = power, y = mean_k, color = filter)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  facet_wrap(~ filter, ncol = 3, scales = "free_y") +
  labs(title = "Mean Connectivity vs Soft Threshold Power",
       x = "Soft Threshold Power", y = "Mean Connectivity",
       color = "Filter") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(outdir, "plots", "S02_connectivity_vs_power.png"), p2, width = 12, height = 8, dpi = 300)
ggsave(file.path(outdir, "plots", "S02_connectivity_vs_power.pdf"), p2, width = 12, height = 8)
cat("  Saved: S02_connectivity_vs_power\n")

# Plot 3: ARI and module count comparison
stab_plot <- stability_results %>%
  tidyr::pivot_longer(cols = c(n_modules, grey_pct, ARI_vs_original),
                      names_to = "metric", values_to = "value")

p3 <- ggplot(stab_plot, aes(x = filter, y = value, fill = filter)) +
  geom_col(width = 0.7) +
  facet_wrap(~ metric, scales = "free_y", ncol = 3) +
  labs(title = "Module Stability Across Filter Levels",
       x = "Variance Filter", y = "Value") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

ggsave(file.path(outdir, "plots", "S03_module_stability.png"), p3, width = 10, height = 5, dpi = 300)
ggsave(file.path(outdir, "plots", "S03_module_stability.pdf"), p3, width = 10, height = 5)
cat("  Saved: S03_module_stability\n")


# ══════════════════════════════════════════════════════════════════════════════
# SUMMARY
# ══════════════════════════════════════════════════════════════════════════════

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("PARAMETER SENSITIVITY SUMMARY\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat("Outlier detection:\n")
cat("  Automated (Z.k < -2.5):", length(auto_outliers), "outliers\n")
cat("  Original (hard-coded):", length(orig_outliers), "outliers\n")
cat("  Agreement:", length(match_orig), "/", length(orig_outliers), "\n\n")

cat("Scale-free topology:\n")
for (i in 1:nrow(best_params)) {
  bp <- best_params[i, ]
  cat(sprintf("  %-10s R²=%.4f at power=%d | R²>0.80: %s\n",
              bp$filter, bp$max_R2, bp$power_at_max_R2,
              ifelse(bp$reaches_80, "YES", "NO")))
}

cat(sprintf("\nRECOMMENDATION: %s filter, power %d\n", best_filter, best_power_final))
cat("  → Run 02_improved_wgcna.R next\n")

cat("\nOutput: results/wgcna_sensitivity/\n")
cat("  Tables:", length(list.files(file.path(outdir, "tables"), pattern = "\\.tsv$")), "\n")
cat("  Plots:", length(list.files(file.path(outdir, "plots"), pattern = "\\.png$")), "\n")

cat("\n=== Script 01 complete ===\n")
