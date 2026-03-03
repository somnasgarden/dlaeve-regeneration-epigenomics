#!/usr/bin/env Rscript
# =============================================================================
# Cluster Script 02: Network Null Models & GENIE3 Coordination Analysis
# =============================================================================
# Purpose: Test whether methylation coordination in the GENIE3 regulatory
#          network is real by generating null models. This is the MOST
#          computationally intensive analysis (loads 500K edges, runs
#          10,000 permutations on network structure).
#
# Requires: 64GB RAM, 8 CPUs, ~4-8 hours
#
# Key questions:
#   1. Is 51.8% same-direction methylation in connected genes significant
#      beyond what we'd expect from gene proximity or module membership?
#   2. Do TF-target edges show stronger methylation coordination than
#      random gene pairs?
#   3. Are high-weight GENIE3 edges more coordinated than low-weight?
#   4. Module-stratified network coordination
#
# Input:
#   DATA/genie3_top500k.tsv
#   results/03_integration/Tables/MXT_gene_level_meth_vs_expression.txt
#   results/02_rnaseq/Part5_WGCNA/data/all_gene_module_assignments.tsv
#   results/05_tf_prediction/prediction_result.txt
#
# Output: results/network_null_models/
# =============================================================================

library(dplyr)
library(parallel)

cat("=== Cluster 02: Network Null Models & GENIE3 Coordination ===\n")
cat("CPUs:", detectCores(), "| Using 8 cores\n\n")

N_PERM <- 10000
N_CORES <- 8

outdir <- "results/network_null_models"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ── Load data ────────────────────────────────────────────────────────────────

cat("Loading GENIE3 (500K edges)... ")
genie3 <- read.delim("DATA/genie3_top500k.tsv", stringsAsFactors = FALSE)
cat(nrow(genie3), "edges loaded\n")

gene_meth <- read.delim("results/03_integration/Tables/MXT_gene_level_meth_vs_expression.txt",
                         stringsAsFactors = FALSE)
modules <- read.delim("results/02_rnaseq/Part5_WGCNA/data/all_gene_module_assignments.tsv",
                       stringsAsFactors = FALSE)

# Build methylation direction lookup
meth_dir <- setNames(gene_meth$meth_direction, gene_meth$gene_id)
meth_diff <- setNames(gene_meth$mean_meth_diff, gene_meth$gene_id)
mod_lookup <- setNames(modules$Module_Color, modules$Gene)

# ── Filter edges where both genes have DMPs ──────────────────────────────────

cat("Filtering edges where both genes have DMPs...\n")
dmp_genes <- unique(gene_meth$gene_id)

genie3$reg_has_dmp <- genie3$regulatoryGene %in% dmp_genes
genie3$tgt_has_dmp <- genie3$targetGene %in% dmp_genes
both_dmp <- genie3 %>% filter(reg_has_dmp & tgt_has_dmp)
cat("  Edges with both genes having DMPs:", nrow(both_dmp), "\n")

# Direction concordance
both_dmp$reg_dir <- meth_dir[both_dmp$regulatoryGene]
both_dmp$tgt_dir <- meth_dir[both_dmp$targetGene]
both_dmp$concordant <- both_dmp$reg_dir == both_dmp$tgt_dir

observed_concordance <- mean(both_dmp$concordant)
cat("  Observed concordance:", round(observed_concordance * 100, 1), "%\n")

# ══════════════════════════════════════════════════════════════════════════════
# NULL MODEL 1: Shuffle methylation labels (preserve network, randomize meth)
# ══════════════════════════════════════════════════════════════════════════════

cat("\n--- Null Model 1: Shuffle methylation labels ---\n")
cat("  ", N_PERM, "permutations...\n")

set.seed(42)
perm_concordance_1 <- mclapply(1:N_PERM, function(i) {
  # Shuffle methylation directions among DMP genes
  shuffled_dir <- sample(meth_dir)
  names(shuffled_dir) <- names(meth_dir)

  reg_d <- shuffled_dir[both_dmp$regulatoryGene]
  tgt_d <- shuffled_dir[both_dmp$targetGene]
  mean(reg_d == tgt_d, na.rm = TRUE)
}, mc.cores = N_CORES)

perm_concordance_1 <- unlist(perm_concordance_1)
p_null1 <- (sum(perm_concordance_1 >= observed_concordance) + 1) / (N_PERM + 1)
z_null1 <- (observed_concordance - mean(perm_concordance_1)) / sd(perm_concordance_1)

cat("  Null model 1 results:\n")
cat("    Observed concordance:", round(observed_concordance * 100, 2), "%\n")
cat("    Expected (null):", round(mean(perm_concordance_1) * 100, 2), "%\n")
cat("    Z-score:", round(z_null1, 2), "\n")
cat("    Permutation p-value:", format(p_null1, digits = 4), "\n")

# ══════════════════════════════════════════════════════════════════════════════
# NULL MODEL 2: Shuffle edges (preserve methylation, randomize network)
# ══════════════════════════════════════════════════════════════════════════════

cat("\n--- Null Model 2: Shuffle edges (keep meth, randomize connections) ---\n")

set.seed(123)
perm_concordance_2 <- mclapply(1:N_PERM, function(i) {
  # Randomly pair DMP genes (preserve number of edges)
  n_edges <- nrow(both_dmp)
  random_reg <- sample(dmp_genes, n_edges, replace = TRUE)
  random_tgt <- sample(dmp_genes, n_edges, replace = TRUE)

  reg_d <- meth_dir[random_reg]
  tgt_d <- meth_dir[random_tgt]
  mean(reg_d == tgt_d, na.rm = TRUE)
}, mc.cores = N_CORES)

perm_concordance_2 <- unlist(perm_concordance_2)
p_null2 <- (sum(perm_concordance_2 >= observed_concordance) + 1) / (N_PERM + 1)
z_null2 <- (observed_concordance - mean(perm_concordance_2)) / sd(perm_concordance_2)

cat("  Null model 2 results:\n")
cat("    Expected (random pairing):", round(mean(perm_concordance_2) * 100, 2), "%\n")
cat("    Z-score:", round(z_null2, 2), "\n")
cat("    Permutation p-value:", format(p_null2, digits = 4), "\n")

# ══════════════════════════════════════════════════════════════════════════════
# NULL MODEL 3: Module-aware null (control for shared module membership)
# ══════════════════════════════════════════════════════════════════════════════

cat("\n--- Null Model 3: Module-aware null ---\n")

both_dmp$reg_mod <- mod_lookup[both_dmp$regulatoryGene]
both_dmp$tgt_mod <- mod_lookup[both_dmp$targetGene]
both_dmp$same_module <- both_dmp$reg_mod == both_dmp$tgt_mod

# Concordance stratified by same/different module
same_mod_conc <- mean(both_dmp$concordant[both_dmp$same_module], na.rm = TRUE)
diff_mod_conc <- mean(both_dmp$concordant[!both_dmp$same_module], na.rm = TRUE)
cat("  Same-module edges:", sum(both_dmp$same_module, na.rm = TRUE),
    "concordance:", round(same_mod_conc * 100, 1), "%\n")
cat("  Cross-module edges:", sum(!both_dmp$same_module, na.rm = TRUE),
    "concordance:", round(diff_mod_conc * 100, 1), "%\n")

# Permute within modules to test if concordance exceeds module effect
set.seed(456)
perm_concordance_3 <- mclapply(1:N_PERM, function(i) {
  # Shuffle methylation within each module
  shuffled <- meth_dir
  for (mod in unique(modules$Module_Color)) {
    mod_genes_with_dmp <- intersect(modules$Gene[modules$Module_Color == mod], names(meth_dir))
    if (length(mod_genes_with_dmp) >= 2) {
      shuffled[mod_genes_with_dmp] <- sample(shuffled[mod_genes_with_dmp])
    }
  }

  reg_d <- shuffled[both_dmp$regulatoryGene]
  tgt_d <- shuffled[both_dmp$targetGene]
  mean(reg_d == tgt_d, na.rm = TRUE)
}, mc.cores = N_CORES)

perm_concordance_3 <- unlist(perm_concordance_3)
p_null3 <- (sum(perm_concordance_3 >= observed_concordance) + 1) / (N_PERM + 1)
z_null3 <- (observed_concordance - mean(perm_concordance_3)) / sd(perm_concordance_3)

cat("  Module-controlled null:\n")
cat("    Expected (within-module shuffle):", round(mean(perm_concordance_3) * 100, 2), "%\n")
cat("    Z-score:", round(z_null3, 2), "\n")
cat("    Permutation p-value:", format(p_null3, digits = 4), "\n")

# ══════════════════════════════════════════════════════════════════════════════
# ANALYSIS 4: Edge weight vs concordance
# ══════════════════════════════════════════════════════════════════════════════

cat("\n--- Analysis 4: Edge weight vs concordance ---\n")

both_dmp$weight_quartile <- cut(both_dmp$weight,
                                 breaks = quantile(both_dmp$weight, probs = c(0, 0.25, 0.5, 0.75, 1)),
                                 labels = c("Q1 (low)", "Q2", "Q3", "Q4 (high)"),
                                 include.lowest = TRUE)

weight_conc <- both_dmp %>%
  group_by(weight_quartile) %>%
  summarise(
    n_edges = n(),
    concordance = mean(concordant, na.rm = TRUE),
    .groups = "drop"
  )

cat("\nConcordance by edge weight quartile:\n")
for (i in 1:nrow(weight_conc)) {
  r <- weight_conc[i, ]
  cat(sprintf("  %-15s n=%5d  concordance=%.1f%%\n",
              r$weight_quartile, r$n_edges, r$concordance * 100))
}

# Trend test
ct <- cor.test(both_dmp$weight, as.numeric(both_dmp$concordant), method = "spearman")
cat("\nEdge weight vs concordance: rho=", round(ct$estimate, 4),
    ", p=", format(ct$p.value, digits = 4), "\n")

write.table(weight_conc, file.path(outdir, "edge_weight_concordance.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# ANALYSIS 5: Per-module network coordination
# ══════════════════════════════════════════════════════════════════════════════

cat("\n--- Analysis 5: Per-module network coordination ---\n")

module_net_results <- data.frame()
for (mod in module_colors) {
  mod_edges <- both_dmp %>%
    filter(reg_mod == mod | tgt_mod == mod)

  if (nrow(mod_edges) >= 20) {
    mod_conc <- mean(mod_edges$concordant, na.rm = TRUE)
    bt <- binom.test(sum(mod_edges$concordant, na.rm = TRUE),
                     sum(!is.na(mod_edges$concordant)), p = 0.5)

    module_net_results <- rbind(module_net_results, data.frame(
      module = mod,
      n_edges = nrow(mod_edges),
      concordance = round(mod_conc * 100, 1),
      p_value = bt$p.value,
      stringsAsFactors = FALSE
    ))

    cat(sprintf("  %-12s n=%5d  concordance=%.1f%%  p=%s\n",
                mod, nrow(mod_edges), mod_conc * 100, format(bt$p.value, digits = 3)))
  }
}

module_net_results$padj <- p.adjust(module_net_results$p_value, method = "BH")
write.table(module_net_results, file.path(outdir, "module_network_coordination.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# SUMMARY
# ══════════════════════════════════════════════════════════════════════════════

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("NETWORK NULL MODEL SUMMARY\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

summary_table <- data.frame(
  null_model = c("1. Shuffle methylation labels",
                 "2. Shuffle edges (random pairs)",
                 "3. Within-module shuffle"),
  observed_pct = round(observed_concordance * 100, 2),
  expected_pct = round(c(mean(perm_concordance_1), mean(perm_concordance_2), mean(perm_concordance_3)) * 100, 2),
  z_score = round(c(z_null1, z_null2, z_null3), 2),
  p_value = c(p_null1, p_null2, p_null3),
  stringsAsFactors = FALSE
)

print(summary_table)

write.table(summary_table, file.path(outdir, "null_model_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\nINTERPRETATION:\n")
cat("  If null model 1 is significant → methylation coordination is real\n")
cat("  If null model 2 is significant → network structure matters\n")
cat("  If null model 3 is significant → coordination exceeds module effect\n")
cat("  If null model 3 is NOT significant → module membership explains concordance\n\n")

cat("=== Cluster script 02 complete ===\n")
