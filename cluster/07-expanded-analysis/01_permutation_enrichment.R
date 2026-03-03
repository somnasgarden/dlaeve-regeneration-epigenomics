#!/usr/bin/env Rscript
# =============================================================================
# Cluster Script 01: Permutation-Based Enrichment Testing
# =============================================================================
# Purpose: Validate all Fisher's exact test enrichment results with
#          permutation-based null distributions. This provides non-parametric
#          p-values that don't assume independence and account for gene
#          structure (gene length, CpG density, etc.)
#
# Requires: 64GB RAM, 8 CPUs, ~2-4 hours
#
# What it tests (with 10,000 permutations each):
#   1. Module DMP enrichment (yellow, blue, etc.) — is it real?
#   2. Module DMR enrichment (blue) — is it real?
#   3. TF DMP depletion (DeepTFactor) — is it real?
#   4. Functional category enrichment (cytoskeleton, zinc fingers, etc.)
#   5. Hub gene DMP enrichment/depletion per module
#   6. Morphogenesis gene DMP depletion
#   7. Hox gene DMP protection (0/27)
#
# Input: Same as local/05-10 scripts (results/ + DATA/)
# Output: results/permutation_tests/
# =============================================================================

library(dplyr)
library(parallel)

cat("=== Cluster 01: Permutation-Based Enrichment Testing ===\n")
cat("CPUs available:", detectCores(), "\n")
cat("Using 8 cores for permutation testing\n\n")

N_PERM <- 10000
N_CORES <- 8

outdir <- "results/permutation_tests"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ── Load data ────────────────────────────────────────────────────────────────

cat("Loading data...\n")

# Gene-level methylation (which genes have DMPs)
gene_meth <- read.delim("results/03_integration/Tables/MXT_gene_level_meth_vs_expression.txt",
                         stringsAsFactors = FALSE)
genes_with_dmps <- unique(gene_meth$gene_id)

# Module assignments
modules <- read.delim("results/02_rnaseq/Part5_WGCNA/data/all_gene_module_assignments.tsv",
                       stringsAsFactors = FALSE)

# Hub genes
hubs <- read.delim("results/02_rnaseq/Part5_WGCNA/data/hub_genes_all_modules.tsv",
                    stringsAsFactors = FALSE)

# DEG results (universe of expressed genes)
deg <- read.csv("results/02_rnaseq/Part2_DEGs/data/DEresults_amputated_vs_control.csv",
                 stringsAsFactors = FALSE)
colnames(deg)[1] <- "gene_id"
universe <- unique(deg$gene_id)

# Annotations
annot <- read.delim("DATA/derLaeGenome_eviann_annotations.tsv",
                     header = FALSE, stringsAsFactors = FALSE)
colnames(annot) <- c("gene_id", "gene_symbol", "description")
annot <- annot %>% distinct(gene_id, .keep_all = TRUE)

# DMRs
dmrs <- read.delim("results/01_methylation/dmrs_annotated.txt", stringsAsFactors = FALSE)
dmr_genes <- unique(dmrs$nearest_gene)

# TF predictions
tf_pred <- read.delim("results/05_tf_prediction/prediction_result.txt", stringsAsFactors = FALSE)
tf_genes <- unique(gsub("-protein-.*", "", tf_pred$sequence_name[tf_pred$prediction == "True"]))

cat("  Universe:", length(universe), "expressed genes\n")
cat("  Genes with DMPs:", length(genes_with_dmps), "\n")
cat("  Genes with DMRs:", length(dmr_genes), "\n")
cat("  Predicted TFs:", length(tf_genes), "\n\n")

# ── Helper: Permutation test function ────────────────────────────────────────

permutation_enrichment <- function(target_genes, marker_genes, universe_genes,
                                    n_perm = N_PERM, n_cores = N_CORES) {
  # Observed overlap
  target_in_universe <- intersect(target_genes, universe_genes)
  marker_in_universe <- intersect(marker_genes, universe_genes)
  observed <- length(intersect(target_in_universe, marker_in_universe))

  n_target <- length(target_in_universe)
  n_universe <- length(universe_genes)

  # Permutation: randomly sample n_target genes and count overlaps
  set.seed(42)
  perm_counts <- mclapply(1:n_perm, function(i) {
    random_target <- sample(universe_genes, n_target, replace = FALSE)
    length(intersect(random_target, marker_in_universe))
  }, mc.cores = n_cores)

  perm_counts <- unlist(perm_counts)

  # Two-sided p-value
  p_enriched <- (sum(perm_counts >= observed) + 1) / (n_perm + 1)
  p_depleted <- (sum(perm_counts <= observed) + 1) / (n_perm + 1)
  p_twosided <- 2 * min(p_enriched, p_depleted)

  # Z-score
  z_score <- (observed - mean(perm_counts)) / sd(perm_counts)

  list(
    observed = observed,
    expected = mean(perm_counts),
    sd = sd(perm_counts),
    z_score = z_score,
    p_enriched = p_enriched,
    p_depleted = p_depleted,
    p_twosided = min(p_twosided, 1),
    direction = ifelse(observed > mean(perm_counts), "enriched", "depleted"),
    n_target = n_target,
    n_marker = length(marker_in_universe),
    perm_counts = perm_counts
  )
}

# ══════════════════════════════════════════════════════════════════════════════
# TEST 1: Module DMP enrichment
# ══════════════════════════════════════════════════════════════════════════════

cat("--- Test 1: Module DMP enrichment (", N_PERM, " permutations) ---\n")

module_colors <- unique(modules$Module_Color)
module_colors <- module_colors[module_colors != "grey"]

module_perm_results <- data.frame()
for (mod in module_colors) {
  mod_genes <- modules$Gene[modules$Module_Color == mod]
  result <- permutation_enrichment(mod_genes, genes_with_dmps, universe)

  cat(sprintf("  %-12s observed=%d, expected=%.1f, z=%.2f, p_enrich=%.4f, p_deplete=%.4f\n",
              mod, result$observed, result$expected, result$z_score,
              result$p_enriched, result$p_depleted))

  module_perm_results <- rbind(module_perm_results, data.frame(
    module = mod,
    n_module = result$n_target,
    observed_dmp = result$observed,
    expected_dmp = round(result$expected, 1),
    z_score = round(result$z_score, 3),
    p_enriched = result$p_enriched,
    p_depleted = result$p_depleted,
    p_twosided = result$p_twosided,
    direction = result$direction,
    stringsAsFactors = FALSE
  ))
}

module_perm_results$padj <- p.adjust(module_perm_results$p_twosided, method = "BH")
write.table(module_perm_results, file.path(outdir, "perm_module_dmp_enrichment.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("\n")

# ══════════════════════════════════════════════════════════════════════════════
# TEST 2: Module DMR enrichment
# ══════════════════════════════════════════════════════════════════════════════

cat("--- Test 2: Module DMR enrichment ---\n")

module_dmr_results <- data.frame()
for (mod in module_colors) {
  mod_genes <- modules$Gene[modules$Module_Color == mod]
  result <- permutation_enrichment(mod_genes, dmr_genes, universe)

  cat(sprintf("  %-12s observed=%d, expected=%.1f, z=%.2f, p_enrich=%.4f\n",
              mod, result$observed, result$expected, result$z_score, result$p_enriched))

  module_dmr_results <- rbind(module_dmr_results, data.frame(
    module = mod,
    n_module = result$n_target,
    observed_dmr = result$observed,
    expected_dmr = round(result$expected, 1),
    z_score = round(result$z_score, 3),
    p_enriched = result$p_enriched,
    p_depleted = result$p_depleted,
    direction = result$direction,
    stringsAsFactors = FALSE
  ))
}

module_dmr_results$padj <- p.adjust(module_dmr_results$p_enriched, method = "BH")
write.table(module_dmr_results, file.path(outdir, "perm_module_dmr_enrichment.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("\n")

# ══════════════════════════════════════════════════════════════════════════════
# TEST 3: TF DMP depletion (DeepTFactor)
# ══════════════════════════════════════════════════════════════════════════════

cat("--- Test 3: TF DMP depletion ---\n")

tf_in_universe <- intersect(tf_genes, universe)
cat("  TFs in expressed universe:", length(tf_in_universe), "\n")

tf_result <- permutation_enrichment(tf_in_universe, genes_with_dmps, universe)
cat(sprintf("  TF DMPs: observed=%d, expected=%.1f, z=%.2f, p_depleted=%.6f\n",
            tf_result$observed, tf_result$expected, tf_result$z_score, tf_result$p_depleted))

# ══════════════════════════════════════════════════════════════════════════════
# TEST 4: Hub gene DMP enrichment per module
# ══════════════════════════════════════════════════════════════════════════════

cat("\n--- Test 4: Hub gene DMP enrichment per module ---\n")

hub_perm_results <- data.frame()
for (mod in module_colors) {
  mod_hubs <- hubs$Gene[hubs$Module == mod & hubs$IsHub == TRUE]
  mod_genes <- modules$Gene[modules$Module_Color == mod]

  if (length(mod_hubs) >= 5) {
    # Test within module: are hubs more likely to have DMPs than non-hubs?
    mod_dmp_genes <- intersect(mod_genes, genes_with_dmps)
    result <- permutation_enrichment(mod_hubs, mod_dmp_genes, mod_genes)

    cat(sprintf("  %-12s hubs=%d, observed=%d, expected=%.1f, z=%.2f, p_enrich=%.4f, p_deplete=%.4f\n",
                mod, length(mod_hubs), result$observed, result$expected,
                result$z_score, result$p_enriched, result$p_depleted))

    hub_perm_results <- rbind(hub_perm_results, data.frame(
      module = mod,
      n_hubs = length(mod_hubs),
      observed = result$observed,
      expected = round(result$expected, 1),
      z_score = round(result$z_score, 3),
      p_enriched = result$p_enriched,
      p_depleted = result$p_depleted,
      direction = result$direction,
      stringsAsFactors = FALSE
    ))
  }
}

hub_perm_results$padj <- p.adjust(pmin(hub_perm_results$p_enriched, hub_perm_results$p_depleted) * 2, method = "BH")
write.table(hub_perm_results, file.path(outdir, "perm_hub_dmp_enrichment.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("\n")

# ══════════════════════════════════════════════════════════════════════════════
# TEST 5: Hox gene protection (exact test)
# ══════════════════════════════════════════════════════════════════════════════

cat("--- Test 5: Hox gene methylation protection ---\n")

hox_genes <- annot$gene_id[grepl("hox|homeobox", tolower(paste(annot$gene_symbol, annot$description)))]
hox_in_universe <- intersect(hox_genes, universe)
cat("  Hox genes in universe:", length(hox_in_universe), "\n")

hox_result <- permutation_enrichment(hox_in_universe, genes_with_dmps, universe)
cat(sprintf("  Hox DMPs: observed=%d, expected=%.1f, z=%.2f, p_depleted=%.6f\n",
            hox_result$observed, hox_result$expected, hox_result$z_score, hox_result$p_depleted))

# Probability of 0 DMPs by chance
p_zero <- sum(hox_result$perm_counts == 0) / N_PERM
cat(sprintf("  P(0 DMPs by chance): %.4f (in %d permutations)\n", p_zero, N_PERM))

# ══════════════════════════════════════════════════════════════════════════════
# TEST 6: Functional category enrichment
# ══════════════════════════════════════════════════════════════════════════════

cat("\n--- Test 6: Functional category DMP enrichment ---\n")

# Classify genes (same as script 08)
classify_gene <- function(symbol, description) {
  s <- tolower(paste(symbol, description, sep = " "))
  if (grepl("morphogen|wnt|notch|hedgehog|bmp|fgf|tgf|sox|pax|hox|homeobox|lhx|dlx|msx|otx|six|eya|fork.?head|foxp|tbx|gata|hand|twist|snail|slug|embryo|limb|regenerat", s)) return("Morphogenesis_Development")
  if (grepl("cell.?cycle|cyclin|cdk|mitot|mitosi|kinetochore|spindle|centromere|checkpoint|proliferat|aurora|polo|plk|cdc|condensin|cohesin", s)) return("Cell_Cycle_Proliferation")
  if (grepl("zinc.?finger|zf\\-|znf|c2h2|cchh", s)) return("Zinc_Finger")
  if (grepl("transcription.?factor|bhlh|basic.?helix|leucine.?zipper|bzip|hmg|high.?mobility|ets|myb|rel|nf.?kb|creb|nuclear.?receptor|homeodomain", s)) return("Core_TF")
  if (grepl("actin|myosin|tubulin|kinesin|dynein|microtubul|tropomyosin|cofilin|profilin|rho|rac|cdc42|focal.?adhesion", s)) return("Cytoskeleton_Motility")
  if (grepl("splicing|spliceosome|snrnp|rna.?helicase|rna.?polymerase|ribosom|translation|elongation.?factor|trna|rrna", s)) return("RNA_Processing")
  return(NA)
}

annot$category <- mapply(classify_gene, annot$gene_symbol, annot$description)
annot_classified <- annot %>% filter(!is.na(category))

categories <- unique(annot_classified$category)
cat_perm_results <- data.frame()

for (cat_name in categories) {
  cat_genes <- annot_classified$gene_id[annot_classified$category == cat_name]
  result <- permutation_enrichment(cat_genes, genes_with_dmps, universe)

  cat(sprintf("  %-25s observed=%d, expected=%.1f, z=%.2f, p_twosided=%.4f (%s)\n",
              cat_name, result$observed, result$expected, result$z_score,
              result$p_twosided, result$direction))

  cat_perm_results <- rbind(cat_perm_results, data.frame(
    category = cat_name,
    n_genes = result$n_target,
    observed = result$observed,
    expected = round(result$expected, 1),
    z_score = round(result$z_score, 3),
    p_enriched = result$p_enriched,
    p_depleted = result$p_depleted,
    p_twosided = result$p_twosided,
    direction = result$direction,
    stringsAsFactors = FALSE
  ))
}

cat_perm_results$padj <- p.adjust(cat_perm_results$p_twosided, method = "BH")
write.table(cat_perm_results, file.path(outdir, "perm_category_dmp_enrichment.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# SUMMARY
# ══════════════════════════════════════════════════════════════════════════════

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("PERMUTATION TEST SUMMARY (", N_PERM, " permutations each)\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat("Module DMP enrichment (padj < 0.05):\n")
sig_mod <- module_perm_results %>% filter(padj < 0.05)
if (nrow(sig_mod) > 0) {
  for (i in 1:nrow(sig_mod)) cat("  ", sig_mod$module[i], ":", sig_mod$direction[i], "padj=", sig_mod$padj[i], "\n")
} else { cat("  None significant\n") }

cat("\nModule DMR enrichment (padj < 0.05):\n")
sig_dmr <- module_dmr_results %>% filter(padj < 0.05)
if (nrow(sig_dmr) > 0) {
  for (i in 1:nrow(sig_dmr)) cat("  ", sig_dmr$module[i], ":", sig_dmr$direction[i], "padj=", sig_dmr$padj[i], "\n")
} else { cat("  None significant\n") }

cat("\nTF depletion: p_depleted=", tf_result$p_depleted, "\n")
cat("Hox protection: p_depleted=", hox_result$p_depleted, ", P(0 by chance)=", p_zero, "\n")

cat("\nHub DMP enrichment (padj < 0.1):\n")
sig_hub <- hub_perm_results %>% filter(padj < 0.1)
if (nrow(sig_hub) > 0) {
  for (i in 1:nrow(sig_hub)) cat("  ", sig_hub$module[i], ":", sig_hub$direction[i], "padj=", sig_hub$padj[i], "\n")
} else { cat("  None significant\n") }

cat("\nCategory enrichment (padj < 0.05):\n")
sig_cat <- cat_perm_results %>% filter(padj < 0.05)
if (nrow(sig_cat) > 0) {
  for (i in 1:nrow(sig_cat)) cat("  ", sig_cat$category[i], ":", sig_cat$direction[i], "padj=", sig_cat$padj[i], "\n")
} else { cat("  None significant\n") }

cat("\n=== Cluster script 01 complete ===\n")
cat("Output:", outdir, "\n")
cat("Tables:", length(list.files(outdir, pattern = "\\.tsv$")), "\n")
