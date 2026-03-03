#!/usr/bin/env Rscript
# =============================================================================
# Cluster Script 03: Bootstrap Confidence Intervals for Key Correlations
# =============================================================================
# Purpose: Generate bootstrap 95% CIs for all key correlations and effect
#          sizes reported in the paper. This strengthens the statistical
#          claims with non-parametric uncertainty estimates.
#
# Requires: 64GB RAM, 8 CPUs, ~2-4 hours
#
# What gets bootstrapped (10,000 resamples each):
#   1. Module × region Spearman correlations (48 tests)
#   2. TF DMP rate confidence intervals
#   3. Hub gene DMP rate CIs per module
#   4. DMP-DMR concordance CI
#   5. Methylation dosage-response correlation CI
#   6. DMP direction bias per module
#
# Input: Same as integration + novel findings scripts
# Output: results/bootstrap_ci/
# =============================================================================

library(dplyr)
library(parallel)

cat("=== Cluster 03: Bootstrap Confidence Intervals ===\n")
cat("CPUs:", detectCores(), "| Using 8 cores\n\n")

N_BOOT <- 10000
N_CORES <- 8

outdir <- "results/bootstrap_ci"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ── Load data ────────────────────────────────────────────────────────────────

gene_meth <- read.delim("results/03_integration/Tables/MXT_gene_level_meth_vs_expression.txt",
                         stringsAsFactors = FALSE)
modules <- read.delim("results/02_rnaseq/Part5_WGCNA/data/all_gene_module_assignments.tsv",
                       stringsAsFactors = FALSE)
hubs <- read.delim("results/02_rnaseq/Part5_WGCNA/data/hub_genes_all_modules.tsv",
                    stringsAsFactors = FALSE)
dmrs <- read.delim("results/01_methylation/dmrs_annotated.txt", stringsAsFactors = FALSE)

# Merge module info
gene_meth <- gene_meth %>%
  left_join(modules %>% select(Gene, Module_Color) %>% rename(gene_id = Gene), by = "gene_id")

# ── Helper: Bootstrap CI function ────────────────────────────────────────────

boot_ci <- function(data, stat_fn, n_boot = N_BOOT, n_cores = N_CORES, alpha = 0.05) {
  observed <- stat_fn(data)

  boot_stats <- mclapply(1:n_boot, function(i) {
    boot_data <- data[sample(nrow(data), replace = TRUE), ]
    stat_fn(boot_data)
  }, mc.cores = n_cores)

  boot_stats <- unlist(boot_stats)
  boot_stats <- boot_stats[is.finite(boot_stats)]

  ci <- quantile(boot_stats, probs = c(alpha/2, 1 - alpha/2))

  list(
    observed = observed,
    ci_lower = ci[1],
    ci_upper = ci[2],
    se = sd(boot_stats),
    n_boot = length(boot_stats)
  )
}

# ══════════════════════════════════════════════════════════════════════════════
# 1. Module × Region Spearman correlations with CIs
# ══════════════════════════════════════════════════════════════════════════════

cat("--- 1. Module × Region Spearman correlations ---\n")

mod_region_results <- data.frame()
module_colors <- unique(gene_meth$Module_Color)
module_colors <- module_colors[!is.na(module_colors) & module_colors != "grey"]

regions <- c("Promoter", "Intergenic", "Intron", "Gene Body")

for (mod in module_colors) {
  for (reg in regions) {
    sub <- gene_meth %>% filter(Module_Color == mod & primary_region == reg)
    if (nrow(sub) >= 15) {
      spearman_fn <- function(d) {
        ct <- suppressWarnings(cor.test(d$mean_meth_diff, d$log2FoldChange, method = "spearman"))
        ct$estimate
      }

      result <- boot_ci(sub, spearman_fn)

      mod_region_results <- rbind(mod_region_results, data.frame(
        module = mod,
        region = reg,
        n = nrow(sub),
        rho = round(result$observed, 4),
        ci_lower = round(result$ci_lower, 4),
        ci_upper = round(result$ci_upper, 4),
        se = round(result$se, 4),
        significant = !(result$ci_lower <= 0 & result$ci_upper >= 0),
        stringsAsFactors = FALSE
      ))

      cat(sprintf("  %-12s × %-12s n=%3d  rho=%.3f  95%%CI=[%.3f, %.3f] %s\n",
                  mod, reg, nrow(sub), result$observed,
                  result$ci_lower, result$ci_upper,
                  ifelse(!(result$ci_lower <= 0 & result$ci_upper >= 0), "*", "")))
    }
  }
}

write.table(mod_region_results, file.path(outdir, "boot_module_region_correlations.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("\n")

# ══════════════════════════════════════════════════════════════════════════════
# 2. Hub gene DMP rate CIs per module
# ══════════════════════════════════════════════════════════════════════════════

cat("--- 2. Hub gene DMP rate bootstrap ---\n")

genes_with_dmps <- unique(gene_meth$gene_id)

hub_ci_results <- data.frame()
for (mod in module_colors) {
  mod_genes <- modules$Gene[modules$Module_Color == mod]
  mod_hubs <- hubs$Gene[hubs$Module == mod & hubs$IsHub == TRUE]
  mod_nonhubs <- setdiff(mod_genes, mod_hubs)

  if (length(mod_hubs) >= 5) {
    # Hub DMP rate
    hub_dmp_rate <- mean(mod_hubs %in% genes_with_dmps)
    nonhub_dmp_rate <- mean(mod_nonhubs %in% genes_with_dmps)

    # Bootstrap the difference
    set.seed(42)
    boot_diffs <- mclapply(1:N_BOOT, function(i) {
      boot_hubs <- sample(mod_hubs, replace = TRUE)
      boot_nonhubs <- sample(mod_nonhubs, replace = TRUE)
      mean(boot_hubs %in% genes_with_dmps) - mean(boot_nonhubs %in% genes_with_dmps)
    }, mc.cores = N_CORES)

    boot_diffs <- unlist(boot_diffs)
    ci <- quantile(boot_diffs, probs = c(0.025, 0.975))
    obs_diff <- hub_dmp_rate - nonhub_dmp_rate

    hub_ci_results <- rbind(hub_ci_results, data.frame(
      module = mod,
      n_hubs = length(mod_hubs),
      hub_rate = round(hub_dmp_rate * 100, 1),
      nonhub_rate = round(nonhub_dmp_rate * 100, 1),
      diff = round(obs_diff * 100, 1),
      ci_lower = round(ci[1] * 100, 1),
      ci_upper = round(ci[2] * 100, 1),
      significant = !(ci[1] <= 0 & ci[2] >= 0),
      stringsAsFactors = FALSE
    ))

    cat(sprintf("  %-12s hubs=%.1f%% nonhubs=%.1f%%  diff=%.1f%%  95%%CI=[%.1f, %.1f] %s\n",
                mod, hub_dmp_rate*100, nonhub_dmp_rate*100, obs_diff*100,
                ci[1]*100, ci[2]*100,
                ifelse(!(ci[1] <= 0 & ci[2] >= 0), "*", "")))
  }
}

write.table(hub_ci_results, file.path(outdir, "boot_hub_dmp_rate.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("\n")

# ══════════════════════════════════════════════════════════════════════════════
# 3. DMP-DMR concordance CI
# ══════════════════════════════════════════════════════════════════════════════

cat("--- 3. DMP-DMR concordance bootstrap ---\n")

# Genes with both
dmr_genes <- unique(dmrs$nearest_gene)
both_genes <- intersect(genes_with_dmps, dmr_genes)

# Build concordance data
dmr_dir <- dmrs %>%
  group_by(nearest_gene) %>%
  summarise(dmr_dir = ifelse(mean(diff_methy) > 0, "Hyper", "Hypo"), .groups = "drop")

conc_data <- gene_meth %>%
  filter(gene_id %in% both_genes) %>%
  select(gene_id, meth_direction) %>%
  left_join(dmr_dir, by = c("gene_id" = "nearest_gene")) %>%
  filter(!is.na(dmr_dir)) %>%
  mutate(concordant = meth_direction == dmr_dir)

conc_fn <- function(d) mean(d$concordant)
conc_result <- boot_ci(conc_data, conc_fn)

cat(sprintf("  Concordance: %.1f%%  95%%CI=[%.1f%%, %.1f%%]\n",
            conc_result$observed * 100, conc_result$ci_lower * 100, conc_result$ci_upper * 100))

# ══════════════════════════════════════════════════════════════════════════════
# 4. Global methylation-expression correlation CI
# ══════════════════════════════════════════════════════════════════════════════

cat("\n--- 4. Global methylation-expression correlation ---\n")

global_fn <- function(d) {
  suppressWarnings(cor(d$mean_meth_diff, d$log2FoldChange, method = "spearman", use = "complete.obs"))
}
global_result <- boot_ci(gene_meth, global_fn)

cat(sprintf("  Global rho: %.4f  95%%CI=[%.4f, %.4f]\n",
            global_result$observed, global_result$ci_lower, global_result$ci_upper))

# ══════════════════════════════════════════════════════════════════════════════
# 5. Module DMP direction bias CIs
# ══════════════════════════════════════════════════════════════════════════════

cat("\n--- 5. Module DMP direction bias ---\n")

dir_ci_results <- data.frame()
for (mod in module_colors) {
  sub <- gene_meth %>% filter(Module_Color == mod & !is.na(meth_direction))
  if (nrow(sub) >= 20) {
    hypo_fn <- function(d) mean(d$meth_direction == "Hypo")
    result <- boot_ci(sub, hypo_fn)

    dir_ci_results <- rbind(dir_ci_results, data.frame(
      module = mod,
      n = nrow(sub),
      pct_hypo = round(result$observed * 100, 1),
      ci_lower = round(result$ci_lower * 100, 1),
      ci_upper = round(result$ci_upper * 100, 1),
      biased = !(result$ci_lower <= 0.5 & result$ci_upper >= 0.5),
      stringsAsFactors = FALSE
    ))

    cat(sprintf("  %-12s n=%3d  %%hypo=%.1f%%  95%%CI=[%.1f%%, %.1f%%] %s\n",
                mod, nrow(sub), result$observed*100,
                result$ci_lower*100, result$ci_upper*100,
                ifelse(!(result$ci_lower <= 0.5 & result$ci_upper >= 0.5), "BIASED", "")))
  }
}

write.table(dir_ci_results, file.path(outdir, "boot_module_direction_bias.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# SUMMARY
# ══════════════════════════════════════════════════════════════════════════════

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("BOOTSTRAP CI SUMMARY (", N_BOOT, " resamples each)\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat("Key CIs for the paper:\n")
cat("  Global meth-expr rho:", round(global_result$observed, 4),
    "[", round(global_result$ci_lower, 4), ",", round(global_result$ci_upper, 4), "]\n")
cat("  DMP-DMR concordance:", round(conc_result$observed*100, 1), "%",
    "[", round(conc_result$ci_lower*100, 1), ",", round(conc_result$ci_upper*100, 1), "]\n")

sig_mods <- mod_region_results %>% filter(significant)
cat("\n  Significant module × region correlations:", nrow(sig_mods), "\n")

sig_hubs <- hub_ci_results %>% filter(significant)
cat("  Significant hub DMP differences:", nrow(sig_hubs), "\n")

biased_mods <- dir_ci_results %>% filter(biased)
cat("  Modules with biased DMP direction:", nrow(biased_mods), "\n")

cat("\n=== Cluster script 03 complete ===\n")
