#!/usr/bin/env Rscript
# =============================================================================
# WGCNA Parameter Exploration with Enrichment Tracking
# =============================================================================
# Purpose: Run WGCNA across a grid of variance filters and soft threshold
#          powers, then run enrichment on each configuration to track how
#          top terms change (especially sulfation/glycosaminoglycan terms
#          that emerge at high thresholds).
#
# Key questions:
#   1. At which thresholds does enrichment shift to sulfation/GAG terms?
#   2. How does connectivity change with threshold?
#   3. Which samples are true outliers across parameter sets?
#   4. How does module composition change?
#   5. Is the sulfation signal from a sparser network capturing ECM biology?
#
# Requires: 64GB RAM, 8 CPUs, ~6-10 hours
# Working dir on cluster: /mnt/data/alfredvar/rlopezt/CorrelationMatrix
#
# Input (waterfall: cluster path → local path):
#   Expression: Part4_Kmeans/data/Part4_objects.RData (cluster CWD)
#               OR results/02_rnaseq/Part5_WGCNA/inputs/normalized_counts.tsv (local)
#   Modules:    Part5_WGCNA/data/all_gene_module_assignments.tsv (cluster CWD)
#               OR results/02_rnaseq/Part5_WGCNA/data/all_gene_module_assignments.tsv (local)
#   STRING DB:  /mnt/data/alfredvar/rlopezt/Preliminary/protein.enrichment.terms.v12.0.txt
#               OR results/01_methylation/protein.enrichment.terms.v12.0.txt (local)
# Output: results/wgcna_exploration/ (relative to CWD)
# =============================================================================

library(WGCNA)
library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)

enableWGCNAThreads(nThreads = 8)
options(stringsAsFactors = FALSE)
options(scipen = 999)

cat("=== WGCNA Parameter Exploration with Enrichment Tracking ===\n")
cat("CPUs:", parallel::detectCores(), "| Using 8 cores\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n\n")

outdir <- "results/wgcna_exploration"
dir.create(file.path(outdir, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "data"), recursive = TRUE, showWarnings = FALSE)

save_both <- function(p, name, w = 10, h = 7) {
  ggsave(file.path(outdir, "plots", paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(outdir, "plots", paste0(name, ".png")), p, width = w, height = h, dpi = 300)
  cat("  Saved:", name, "\n")
}


# =============================================================================
# SECTION 1: DATA LOADING
# =============================================================================

cat("--- Section 1: Loading Data ---\n\n")

# Waterfall data loading — cluster paths first, then local paths
expr_data <- NULL
meta_data <- NULL

# Try 1: Cluster CWD (Part4 RData — standard WGCNA input on cluster)
rdata_file <- "Part4_Kmeans/data/Part4_objects.RData"
if (file.exists(rdata_file)) {
  load(rdata_file)
  expr_data <- normalized_counts
  meta_data <- df
  cat("  Loaded from Part4 objects (cluster CWD)\n")
}

# Try 2: Cluster WGCNA saved TSV inputs
if (is.null(expr_data)) {
  wgcna_input <- "Part5_WGCNA/inputs/normalized_counts.tsv"
  wgcna_meta  <- "Part5_WGCNA/inputs/sample_metadata.tsv"
  if (file.exists(wgcna_input)) {
    expr_data <- read.delim(wgcna_input, check.names = FALSE)
    if (file.exists(wgcna_meta)) meta_data <- read.delim(wgcna_meta, check.names = FALSE)
    cat("  Loaded from Part5_WGCNA inputs (cluster CWD)\n")
  }
}

# Try 3: Local repo paths (when running from repo root)
if (is.null(expr_data)) {
  local_input <- "results/02_rnaseq/Part5_WGCNA/inputs/normalized_counts.tsv"
  local_meta  <- "results/02_rnaseq/Part5_WGCNA/inputs/sample_metadata.tsv"
  if (file.exists(local_input)) {
    expr_data <- read.delim(local_input, check.names = FALSE)
    if (file.exists(local_meta)) meta_data <- read.delim(local_meta, check.names = FALSE)
    cat("  Loaded from local repo path\n")
  }
}

if (is.null(expr_data)) stop("Cannot find expression data. Check paths.")

cat("  Dimensions:", nrow(expr_data), "genes x", ncol(expr_data), "samples\n")
expr_data <- expr_data[complete.cases(expr_data), ]
cat("  After NA removal:", nrow(expr_data), "genes\n")

# Load original module assignments (waterfall: cluster → local)
orig_modules <- NULL
for (mpath in c("Part5_WGCNA/data/all_gene_module_assignments.tsv",
                "results/02_rnaseq/Part5_WGCNA/data/all_gene_module_assignments.tsv")) {
  if (file.exists(mpath)) {
    orig_modules <- read.delim(mpath)
    cat("  Original modules loaded from:", mpath, "\n")
    break
  }
}
if (is.null(orig_modules)) stop("Cannot find original module assignments.")

# Load STRING enrichment terms (waterfall: cluster → local)
enrichment_db <- NULL
for (epath in c("/mnt/data/alfredvar/rlopezt/Preliminary/protein.enrichment.terms.v12.0.txt",
                "Part5_WGCNA/inputs/protein.enrichment.terms.v12.0.txt",
                "results/01_methylation/protein.enrichment.terms.v12.0.txt",
                "results/01_methylation/tables/protein.enrichment.terms.v12.0.txt")) {
  if (file.exists(epath)) {
    enrichment_db <- read.delim(epath, header = TRUE, quote = "")
    # Standardize column names (file header starts with #)
    colnames(enrichment_db) <- c("string_protein_id", "category", "term", "description")
    # Strip species prefix: "12345.LOC_00001234" → "LOC_00001234"
    enrichment_db$gene_id <- sub("^[^.]+\\.", "", enrichment_db$string_protein_id)
    cat("  STRING enrichment DB loaded from:", epath, "\n")
    cat("  Categories:", paste(unique(enrichment_db$category), collapse = ", "), "\n")
    cat("  Genes in DB:", length(unique(enrichment_db$gene_id)), "\n")
    break
  }
}
if (is.null(enrichment_db)) {
  cat("  WARNING: STRING enrichment DB not found. Enrichment will be skipped.\n")
}

cat("\n")


# =============================================================================
# SECTION 2: OUTLIER DETECTION (Multiple Methods)
# =============================================================================

cat("--- Section 2: Outlier Detection ---\n\n")

expr_t <- t(expr_data)

# Method 1: Standardized connectivity (Z.k) — at multiple powers
outlier_methods <- data.frame()

for (test_power in c(6, 8, 12, 16)) {
  A <- adjacency(expr_t, power = test_power, type = "signed")
  k <- as.numeric(apply(A, 2, sum)) - 1
  Z_k <- scale(k)

  auto_out <- rownames(expr_t)[Z_k < -2.5]
  for (s in rownames(expr_t)) {
    outlier_methods <- rbind(outlier_methods, data.frame(
      sample = s,
      method = paste0("Zk_power", test_power),
      Z_k = as.numeric(Z_k[which(rownames(expr_t) == s)]),
      is_outlier = s %in% auto_out
    ))
  }
  cat(sprintf("  Z.k (power=%d): %d outliers\n", test_power, length(auto_out)))
  if (length(auto_out) > 0) cat("    ", paste(auto_out, collapse = ", "), "\n")
}

# Method 2: Hierarchical clustering height
sampleTree <- hclust(dist(expr_t), method = "average")
for (cut_h in c(100, 120, 150)) {
  clust <- cutreeStatic(sampleTree, cutHeight = cut_h, minSize = 3)
  outlier_names <- rownames(expr_t)[clust == 0]
  for (s in rownames(expr_t)) {
    outlier_methods <- rbind(outlier_methods, data.frame(
      sample = s,
      method = paste0("hclust_h", cut_h),
      Z_k = NA,
      is_outlier = s %in% outlier_names
    ))
  }
  cat(sprintf("  Hclust (h=%d): %d outliers\n", cut_h, length(outlier_names)))
  if (length(outlier_names) > 0) cat("    ", paste(outlier_names, collapse = ", "), "\n")
}

# Consensus: samples flagged by >= 3 methods
outlier_votes <- outlier_methods %>%
  group_by(sample) %>%
  summarise(n_flagged = sum(is_outlier), n_methods = n(), .groups = "drop") %>%
  mutate(consensus_outlier = n_flagged >= 3)

cat("\n  Consensus outliers (flagged by >= 3 methods):\n")
consensus_outliers <- outlier_votes %>% filter(consensus_outlier)
if (nrow(consensus_outliers) > 0) {
  for (i in 1:nrow(consensus_outliers)) {
    cat(sprintf("    %s: flagged %d/%d times\n",
                consensus_outliers$sample[i],
                consensus_outliers$n_flagged[i],
                consensus_outliers$n_methods[i]))
  }
} else {
  cat("    None — using original outliers (T1S5, dcrep4, irrep7)\n")
}

write.table(outlier_votes, file.path(outdir, "tables", "E01_outlier_consensus.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Use original outliers for consistency
outliers_to_remove <- c("T1S5", "dcrep4", "irrep7")
outliers_to_remove <- outliers_to_remove[outliers_to_remove %in% rownames(expr_t)]
expr_t_clean <- expr_t[!rownames(expr_t) %in% outliers_to_remove, ]
cat(sprintf("\n  Removed %d outliers → %d samples\n\n", length(outliers_to_remove), nrow(expr_t_clean)))


# =============================================================================
# SECTION 3: PARAMETER GRID — VARIANCE FILTERS × SOFT THRESHOLDS
# =============================================================================

cat("--- Section 3: Parameter Grid Sweep ---\n\n")

gene_variance <- apply(expr_data, 1, var)
gene_mad <- apply(expr_data, 1, mad)

# Variance filters — wider range to explore
filter_configs <- list(
  list(name = "top2000",  n = 2000,  desc = "Top 2,000 by variance"),
  list(name = "top3000",  n = 3000,  desc = "Top 3,000 by variance"),
  list(name = "top5000",  n = 5000,  desc = "Top 5,000 by variance"),
  list(name = "top8000",  n = 8000,  desc = "Top 8,000 by variance"),
  list(name = "top10000", n = 10000, desc = "Top 10,000 by variance"),
  list(name = "MAD75",    n = NA,    desc = "Top 75th percentile by MAD"),
  list(name = "MAD50",    n = NA,    desc = "Top 50th percentile by MAD")
)

# Build filtered matrices
filtered_matrices <- list()
for (cfg in filter_configs) {
  if (cfg$name == "MAD75") {
    keep <- names(gene_mad)[gene_mad > quantile(gene_mad, 0.25)]
  } else if (cfg$name == "MAD50") {
    keep <- names(gene_mad)[gene_mad > quantile(gene_mad, 0.50)]
  } else {
    keep <- names(sort(gene_variance, decreasing = TRUE))[1:min(cfg$n, length(gene_variance))]
  }
  mat <- expr_t_clean[, keep]
  gsg <- goodSamplesGenes(mat, verbose = 0)
  if (!gsg$allOK) mat <- mat[gsg$goodSamples, gsg$goodGenes]
  filtered_matrices[[cfg$name]] <- mat
  cat(sprintf("  %-10s: %d genes\n", cfg$name, ncol(mat)))
}

# Powers to test (focused range + extremes)
powers <- c(seq(4, 20, by = 2), 24, 28)

cat("\n  Testing", length(filter_configs), "filters x", length(powers), "powers =",
    length(filter_configs) * length(powers), "combinations\n\n")


# =============================================================================
# SECTION 4: SOFT THRESHOLD SWEEP + MODULE DETECTION
# =============================================================================

cat("--- Section 4: Module Detection per Parameter Set ---\n\n")

sft_results <- data.frame()
all_module_results <- list()
all_connectivity <- data.frame()

for (cfg in filter_configs) {
  mat <- filtered_matrices[[cfg$name]]
  cat(sprintf("  [%s] Soft threshold sweep...\n", cfg$name))

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

  # Run module detection at selected powers (low, medium, high)
  test_powers <- c(
    powers[which(fit_vals >= 0.8)[1]],  # First power at R2 >= 0.8
    8,   # Low threshold
    12,  # Medium (original)
    16,  # High
    20   # Very high
  )
  test_powers <- unique(na.omit(test_powers))
  test_powers <- test_powers[test_powers %in% powers]

  for (tp in test_powers) {
    cat(sprintf("    power=%d: blockwiseModules... ", tp))

    net <- blockwiseModules(
      mat,
      power = tp,
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
      verbose = 0
    )

    module_colors <- labels2colors(net$colors)
    names(module_colors) <- colnames(mat)
    n_modules <- length(unique(module_colors)) - 1
    grey_pct <- round(sum(module_colors == "grey") / length(module_colors) * 100, 1)

    # Connectivity metrics
    k_intra <- intramodularConnectivity(
      adjacency(mat, power = tp, type = "signed",
                corFnc = "bicor", corOptions = list(maxPOutliers = 0.05)),
      module_colors
    )

    mean_kWithin <- mean(k_intra$kWithin, na.rm = TRUE)
    mean_kTotal <- mean(k_intra$kTotal, na.rm = TRUE)

    all_connectivity <- rbind(all_connectivity, data.frame(
      filter = cfg$name,
      n_genes = ncol(mat),
      power = tp,
      n_modules = n_modules,
      grey_pct = grey_pct,
      mean_kTotal = round(mean_kTotal, 2),
      mean_kWithin = round(mean_kWithin, 2),
      max_module_size = max(table(module_colors[module_colors != "grey"])),
      min_module_size = min(table(module_colors[module_colors != "grey"]))
    ))

    # Store module assignments for enrichment
    key <- paste0(cfg$name, "_p", tp)
    all_module_results[[key]] <- data.frame(
      gene = colnames(mat),
      module = module_colors,
      stringsAsFactors = FALSE
    )

    cat(sprintf("%d modules, %.1f%% grey, mean_k=%.1f\n", n_modules, grey_pct, mean_kWithin))

    rm(net, k_intra)
    gc(verbose = FALSE)
  }
}

write.table(sft_results, file.path(outdir, "tables", "E02_soft_threshold_sweep.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(all_connectivity, file.path(outdir, "tables", "E03_connectivity_by_params.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)


# =============================================================================
# SECTION 5: ENRICHMENT ANALYSIS PER PARAMETER SET
# =============================================================================

cat("\n--- Section 5: Enrichment Analysis per Parameter Set ---\n\n")

# Fisher's exact enrichment function (same as cluster/03)
run_enrichment <- function(gene_list, background_genes, enrichment_db, category_filter = "Process") {
  if (is.null(enrichment_db) || length(gene_list) < 5) return(data.frame())

  db_sub <- enrichment_db[enrichment_db$category == category_filter, ]
  terms <- unique(db_sub$term)

  results <- data.frame()
  for (term in terms) {
    term_genes <- unique(db_sub$gene_id[db_sub$term == term])

    a <- length(intersect(gene_list, term_genes))
    b <- length(gene_list) - a
    c_val <- length(intersect(background_genes, term_genes)) - a
    d <- length(background_genes) - a - b - c_val

    if (a >= 2 && d > 0) {
      ft <- fisher.test(matrix(c(a, b, c_val, d), nrow = 2), alternative = "greater")
      desc <- db_sub$description[db_sub$term == term][1]
      results <- rbind(results, data.frame(
        term = term,
        description = desc,
        n_overlap = a,
        n_term = length(term_genes),
        n_module = length(gene_list),
        odds_ratio = round(ft$estimate, 3),
        p_value = ft$p.value,
        stringsAsFactors = FALSE
      ))
    }
  }

  if (nrow(results) > 0) {
    results$padj <- p.adjust(results$p_value, method = "BH")
    results <- results %>% arrange(p_value)
  }
  return(results)
}

# Sulfation/GAG-related search terms
sulfation_keywords <- c("sulfat", "sulphat", "glycosaminoglycan", "heparan",
                         "chondroitin", "keratan", "proteoglycan", "GAG",
                         "carbohydrate sulfotransferase", "extracellular matrix")

all_enrichment_top <- data.frame()
all_sulfation_hits <- data.frame()

if (!is.null(enrichment_db)) {

  background_genes <- unique(enrichment_db$gene_id)

  for (key in names(all_module_results)) {
    mod_data <- all_module_results[[key]]
    unique_mods <- unique(mod_data$module)
    unique_mods <- unique_mods[unique_mods != "grey"]

    cat(sprintf("  [%s] Enrichment for %d modules... ", key, length(unique_mods)))

    n_sulfation_total <- 0

    for (mod in unique_mods) {
      mod_genes <- mod_data$gene[mod_data$module == mod]

      for (cat_name in c("Process", "Function", "KEGG", "Reactome")) {
        enrich <- run_enrichment(mod_genes, background_genes, enrichment_db, cat_name)

        if (nrow(enrich) > 0) {
          sig_enrich <- enrich %>% filter(padj < 0.05)

          if (nrow(sig_enrich) > 0) {
            # Store top 3 terms per category per module
            top_terms <- head(sig_enrich, 3)
            for (j in 1:nrow(top_terms)) {
              all_enrichment_top <- rbind(all_enrichment_top, data.frame(
                param_key = key,
                module = mod,
                category = cat_name,
                rank = j,
                term = top_terms$term[j],
                description = top_terms$description[j],
                n_overlap = top_terms$n_overlap[j],
                odds_ratio = top_terms$odds_ratio[j],
                padj = top_terms$padj[j],
                stringsAsFactors = FALSE
              ))
            }

            # Check for sulfation terms
            sulfation_hits <- sig_enrich %>%
              filter(grepl(paste(sulfation_keywords, collapse = "|"),
                           description, ignore.case = TRUE))

            if (nrow(sulfation_hits) > 0) {
              n_sulfation_total <- n_sulfation_total + nrow(sulfation_hits)
              for (j in 1:nrow(sulfation_hits)) {
                all_sulfation_hits <- rbind(all_sulfation_hits, data.frame(
                  param_key = key,
                  module = mod,
                  category = cat_name,
                  description = sulfation_hits$description[j],
                  n_overlap = sulfation_hits$n_overlap[j],
                  odds_ratio = sulfation_hits$odds_ratio[j],
                  padj = sulfation_hits$padj[j],
                  stringsAsFactors = FALSE
                ))
              }
            }
          }
        }
      }
    }

    cat(sprintf("%d enriched terms, %d sulfation-related\n",
                sum(all_enrichment_top$param_key == key),
                n_sulfation_total))
  }

  write.table(all_enrichment_top, file.path(outdir, "tables", "E04_enrichment_top_terms.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(all_sulfation_hits, file.path(outdir, "tables", "E05_sulfation_enrichment.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
}


# =============================================================================
# SECTION 6: MODULE STABILITY (ARI)
# =============================================================================

cat("\n--- Section 6: Module Stability vs Original ---\n\n")

compute_ari <- function(labels1, labels2) {
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

for (key in names(all_module_results)) {
  new_mods <- all_module_results[[key]]
  shared <- intersect(new_mods$gene, orig_modules$Gene)

  if (length(shared) >= 100) {
    orig_labels <- orig_modules$Module_Color[match(shared, orig_modules$Gene)]
    new_labels <- new_mods$module[match(shared, new_mods$gene)]
    ari <- compute_ari(orig_labels, new_labels)
  } else {
    ari <- NA
  }

  parts <- strsplit(key, "_p")[[1]]
  filter_name <- parts[1]
  power_val <- as.numeric(parts[2])

  stability_results <- rbind(stability_results, data.frame(
    param_key = key,
    filter = filter_name,
    power = power_val,
    n_genes = nrow(new_mods),
    n_modules = length(unique(new_mods$module)) - 1,
    grey_pct = round(sum(new_mods$module == "grey") / nrow(new_mods) * 100, 1),
    n_shared = length(shared),
    ARI = round(ari, 4),
    stringsAsFactors = FALSE
  ))

  cat(sprintf("  %-25s ARI=%.4f  %d modules  %.1f%% grey\n",
              key, ifelse(is.na(ari), 0, ari),
              length(unique(new_mods$module)) - 1,
              sum(new_mods$module == "grey") / nrow(new_mods) * 100))
}

write.table(stability_results, file.path(outdir, "tables", "E06_module_stability.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)


# =============================================================================
# SECTION 7: ENRICHMENT SHIFT ANALYSIS
# =============================================================================

cat("\n--- Section 7: How Enrichment Changes with Parameters ---\n\n")

if (nrow(all_enrichment_top) > 0) {

  # Parse filter and power from param_key
  all_enrichment_top$filter <- sub("_p\\d+$", "", all_enrichment_top$param_key)
  all_enrichment_top$power <- as.numeric(sub(".*_p(\\d+)$", "\\1", all_enrichment_top$param_key))

  # Count unique enriched terms per parameter set
  terms_per_param <- all_enrichment_top %>%
    group_by(param_key, filter, power) %>%
    summarise(
      n_enriched_terms = n(),
      n_unique_descriptions = n_distinct(description),
      top_term = description[1],
      .groups = "drop"
    ) %>%
    arrange(filter, power)

  cat("  Enriched terms per parameter set:\n")
  print(as.data.frame(terms_per_param), row.names = FALSE)

  write.table(terms_per_param, file.path(outdir, "tables", "E07_enrichment_shift_summary.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)

  # Sulfation emergence tracking
  if (nrow(all_sulfation_hits) > 0) {
    all_sulfation_hits$filter <- sub("_p\\d+$", "", all_sulfation_hits$param_key)
    all_sulfation_hits$power <- as.numeric(sub(".*_p(\\d+)$", "\\1", all_sulfation_hits$param_key))

    sulfation_summary <- all_sulfation_hits %>%
      group_by(filter, power) %>%
      summarise(
        n_sulfation_terms = n(),
        best_padj = min(padj),
        terms = paste(unique(description), collapse = "; "),
        .groups = "drop"
      ) %>%
      arrange(filter, power)

    cat("\n  SULFATION ENRICHMENT EMERGENCE:\n")
    for (i in 1:nrow(sulfation_summary)) {
      r <- sulfation_summary[i, ]
      cat(sprintf("    %-10s power=%-2d  %d terms (best padj=%.2e)\n",
                  r$filter, r$power, r$n_sulfation_terms, r$best_padj))
      cat(sprintf("      %s\n", substr(r$terms, 1, 120)))
    }

    write.table(sulfation_summary,
                file.path(outdir, "tables", "E08_sulfation_emergence.tsv"),
                sep = "\t", row.names = FALSE, quote = FALSE)
  } else {
    cat("  No sulfation-related terms found in any parameter set.\n")
    cat("  This means either: (a) the terms use different vocabulary, or\n")
    cat("  (b) the effect only appears with specific enrichment databases.\n")
  }
}


# =============================================================================
# SECTION 8: CONNECTIVITY vs ENRICHMENT CORRELATION
# =============================================================================

cat("\n--- Section 8: Connectivity vs Enrichment ---\n\n")

if (nrow(all_enrichment_top) > 0 && nrow(all_connectivity) > 0) {

  # Merge connectivity with enrichment counts
  conn_enrich <- all_connectivity %>%
    mutate(param_key = paste0(filter, "_p", power)) %>%
    left_join(
      all_enrichment_top %>%
        group_by(param_key) %>%
        summarise(n_enriched = n(), .groups = "drop"),
      by = "param_key"
    )
  conn_enrich$n_enriched[is.na(conn_enrich$n_enriched)] <- 0

  # Also add sulfation count
  if (nrow(all_sulfation_hits) > 0) {
    conn_enrich <- conn_enrich %>%
      left_join(
        all_sulfation_hits %>%
          group_by(param_key) %>%
          summarise(n_sulfation = n(), .groups = "drop"),
        by = "param_key"
      )
    conn_enrich$n_sulfation[is.na(conn_enrich$n_sulfation)] <- 0
  }

  # Correlation: mean connectivity vs number of enriched terms
  if (nrow(conn_enrich) >= 5) {
    ct <- cor.test(conn_enrich$mean_kWithin, conn_enrich$n_enriched, method = "spearman")
    cat(sprintf("  Connectivity vs enrichment: rho = %.3f, p = %.4f\n",
                ct$estimate, ct$p.value))

    if ("n_sulfation" %in% colnames(conn_enrich)) {
      ct2 <- cor.test(conn_enrich$mean_kWithin,
                       conn_enrich$n_sulfation, method = "spearman")
      cat(sprintf("  Connectivity vs sulfation: rho = %.3f, p = %.4f\n",
                  ct2$estimate, ct2$p.value))
      cat(sprintf("  → %s connectivity correlates with %s sulfation enrichment\n",
                  ifelse(ct2$estimate < 0, "LOWER", "HIGHER"),
                  ifelse(ct2$estimate < 0, "MORE", "LESS")))
    }
  }

  write.table(conn_enrich, file.path(outdir, "tables", "E09_connectivity_vs_enrichment.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
}


# =============================================================================
# SECTION 9: VISUALIZATION
# =============================================================================

cat("\n--- Section 9: Visualization ---\n\n")

# Plot 1: R² vs power by filter
p1 <- ggplot(sft_results, aes(x = power, y = R2, color = filter)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  geom_hline(yintercept = 0.80, linetype = "dashed", color = "red") +
  facet_wrap(~ filter, ncol = 4) +
  labs(title = "Scale-Free Topology R² vs Soft Threshold Power",
       subtitle = "Red dashed = R² = 0.80 target",
       x = "Soft Threshold Power", y = "R²") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none")
save_both(p1, "E01_R2_vs_power", 14, 8)

# Plot 2: Mean connectivity by filter and power
p2 <- ggplot(sft_results, aes(x = power, y = log10(mean_k + 1), color = filter)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  facet_wrap(~ filter, ncol = 4, scales = "free_y") +
  labs(title = "Mean Connectivity vs Power (log scale)",
       subtitle = "Higher power → sparser network → lower connectivity",
       x = "Soft Threshold Power", y = "log10(Mean Connectivity + 1)") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none")
save_both(p2, "E02_connectivity_vs_power", 14, 8)

# Plot 3: Module count and grey% by parameter set
if (nrow(all_connectivity) > 0) {
  conn_long <- all_connectivity %>%
    pivot_longer(cols = c(n_modules, grey_pct, mean_kWithin),
                 names_to = "metric", values_to = "value")

  p3 <- ggplot(conn_long, aes(x = power, y = value, color = filter, group = filter)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    facet_wrap(~ metric, scales = "free_y", ncol = 3) +
    labs(title = "Network Properties Across Parameters",
         x = "Soft Threshold Power", y = "Value", color = "Filter") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")
  save_both(p3, "E03_network_properties", 14, 6)
}

# Plot 4: ARI stability heatmap
if (nrow(stability_results) > 1) {
  stab_wide <- stability_results %>%
    select(filter, power, ARI) %>%
    pivot_wider(names_from = power, values_from = ARI) %>%
    as.data.frame()
  rownames(stab_wide) <- stab_wide$filter
  stab_wide$filter <- NULL

  if (ncol(stab_wide) > 1 && nrow(stab_wide) > 1) {
    stab_mat <- as.matrix(stab_wide)
    stab_mat[is.na(stab_mat)] <- 0

    png(file.path(outdir, "plots", "E04_ARI_heatmap.png"), width = 1000, height = 600, res = 150)
    pheatmap(stab_mat,
             main = "Module Stability (ARI vs Original WGCNA)",
             color = colorRampPalette(c("white", "#3498DB", "#E74C3C"))(50),
             cluster_rows = FALSE, cluster_cols = FALSE,
             display_numbers = TRUE, number_format = "%.3f",
             fontsize = 10)
    dev.off()
    cat("  Saved: E04_ARI_heatmap\n")
  }
}

# Plot 5: Sulfation enrichment emergence
if (nrow(all_sulfation_hits) > 0) {
  sulfation_plot <- all_sulfation_hits %>%
    group_by(filter, power) %>%
    summarise(n_terms = n(), best_padj = min(padj), .groups = "drop")

  p5 <- ggplot(sulfation_plot, aes(x = power, y = n_terms, fill = filter)) +
    geom_col(position = "dodge", width = 1.5) +
    labs(title = "Sulfation/GAG Enrichment Emergence by Parameters",
         subtitle = "More sulfation terms appear at specific filter × power combinations",
         x = "Soft Threshold Power", y = "Number of Sulfation Terms",
         fill = "Variance Filter") +
    theme_minimal(base_size = 12)
  save_both(p5, "E05_sulfation_emergence", 12, 6)
}

# Plot 6: Connectivity vs enrichment scatter
if (exists("conn_enrich") && nrow(conn_enrich) > 3) {
  p6 <- ggplot(conn_enrich, aes(x = mean_kWithin, y = n_enriched,
                                 color = filter, size = n_modules)) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE, color = "grey40", linewidth = 0.5) +
    labs(title = "Mean Connectivity vs Number of Enriched Terms",
         subtitle = "Does sparser network reveal different biology?",
         x = "Mean kWithin", y = "Total Enriched Terms (padj < 0.05)",
         color = "Filter", size = "N Modules") +
    theme_minimal(base_size = 12)
  save_both(p6, "E06_connectivity_vs_enrichment", 10, 7)
}


# =============================================================================
# SECTION 10: INTERPRETATION SUMMARY
# =============================================================================

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("WGCNA PARAMETER EXPLORATION SUMMARY\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat("PARAMETER GRID:\n")
cat(sprintf("  Filters tested: %d\n", length(filter_configs)))
cat(sprintf("  Powers tested: %s\n", paste(powers, collapse = ", ")))
cat(sprintf("  Module configurations run: %d\n\n", length(all_module_results)))

cat("OUTLIER DETECTION:\n")
cat(sprintf("  Methods tested: %d\n", length(unique(outlier_methods$method))))
if (nrow(consensus_outliers) > 0) {
  cat(sprintf("  Consensus outliers: %s\n", paste(consensus_outliers$sample, collapse = ", ")))
} else {
  cat("  No consensus outliers beyond original (T1S5, dcrep4, irrep7)\n")
}

cat("\nCONNECTIVITY:\n")
if (nrow(all_connectivity) > 0) {
  for (f in unique(all_connectivity$filter)) {
    f_data <- all_connectivity[all_connectivity$filter == f, ]
    cat(sprintf("  %-10s: kWithin ranges from %.1f (power=%d) to %.1f (power=%d)\n",
                f, max(f_data$mean_kWithin), f_data$power[which.max(f_data$mean_kWithin)],
                min(f_data$mean_kWithin), f_data$power[which.min(f_data$mean_kWithin)]))
  }
}

cat("\nSULFATION ENRICHMENT:\n")
if (nrow(all_sulfation_hits) > 0) {
  cat(sprintf("  Total sulfation-related hits: %d across %d parameter sets\n",
              nrow(all_sulfation_hits), length(unique(all_sulfation_hits$param_key))))
  cat("  Parameter sets with sulfation:\n")
  for (pk in unique(all_sulfation_hits$param_key)) {
    n <- sum(all_sulfation_hits$param_key == pk)
    cat(sprintf("    %s: %d terms\n", pk, n))
  }
} else {
  cat("  No sulfation terms detected in any configuration.\n")
  cat("  Try: (1) check if STRING DB has sulfotransferase terms, or\n")
  cat("       (2) search for heparan/chondroitin sulfate manually in results\n")
}

cat("\nINTERPRETATION:\n")
cat("  High threshold → sparser network → fewer but tighter modules\n")
cat("  If sulfation appears at high thresholds:\n")
cat("    → These genes are VERY tightly co-expressed\n")
cat("    → They form a coherent ECM/proteoglycan module\n")
cat("    → At lower thresholds they're diluted into larger modules\n")
cat("    → Biologically: tight co-expression of ECM remodeling genes\n")
cat("       is consistent with coordinated wound healing response\n")
cat("  Less connectivity at high threshold:\n")
cat("    → Only the strongest correlations survive\n")
cat("    → Modules represent truly co-regulated gene programs\n")
cat("    → If sulfation genes cluster tightly, they are under\n")
cat("       shared regulatory control (possibly through methylation)\n")

cat(sprintf("\nOutput: %s/\n", outdir))
cat(sprintf("  Tables: %d files\n", length(list.files(file.path(outdir, "tables")))))
cat(sprintf("  Plots: %d files\n", length(list.files(file.path(outdir, "plots")))))

# Save all objects
save(all_module_results, all_connectivity, stability_results,
     all_enrichment_top, all_sulfation_hits, sft_results, outlier_votes,
     file = file.path(outdir, "data", "exploration_results.RData"))

cat("\n=== Script 04 complete ===\n")
