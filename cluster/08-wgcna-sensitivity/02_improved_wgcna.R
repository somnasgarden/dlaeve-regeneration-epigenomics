#!/usr/bin/env Rscript
# =============================================================================
# Improved WGCNA with Best Parameters + Bicor
# =============================================================================
# Purpose: Run WGCNA with the optimal parameters identified in script 01,
#          using biweight midcorrelation (bicor) for robustness. Save TOM
#          for methylation-network integration in script 03.
#
# Requires: 64GB RAM, 8 CPUs, ~1-2 hours
# Input:  results/wgcna_sensitivity/data/parameter_sensitivity_results.RData
# Output: results/wgcna_sensitivity/
# =============================================================================

library(WGCNA)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(RColorBrewer)

enableWGCNAThreads(nThreads = 8)
options(stringsAsFactors = FALSE)

cat("=== Improved WGCNA with Bicor ===\n")
cat("CPUs:", parallel::detectCores(), "| Using 8 cores\n\n")

outdir <- "results/wgcna_sensitivity"

# ── Load Part A results ─────────────────────────────────────────────────────

cat("Loading parameter sensitivity results...\n")
load(file.path(outdir, "data", "parameter_sensitivity_results.RData"))

cat("  Best filter:", best_filter, "\n")
cat("  Best power:", best_power_final, "\n")

expr_matrix <- filtered_matrices[[best_filter]]
cat("  Expression matrix:", nrow(expr_matrix), "samples x", ncol(expr_matrix), "genes\n\n")

# ── Load metadata ────────────────────────────────────────────────────────────

# Try loading metadata (waterfall: cluster → local)
meta_data <- NULL
for (mpath in c("Part5_WGCNA/inputs/sample_metadata.tsv",
                "results/02_rnaseq/Part5_WGCNA/inputs/sample_metadata.tsv")) {
  if (file.exists(mpath)) {
    meta_data <- read.delim(mpath, check.names = FALSE)
    cat("  Metadata from:", mpath, "\n")
    break
  }
}
if (is.null(meta_data)) {
  rdata_file <- "Part4_Kmeans/data/Part4_objects.RData"
  if (file.exists(rdata_file)) {
    e <- new.env()
    load(rdata_file, envir = e)
    if ("df" %in% ls(e)) meta_data <- e$df
    cat("  Metadata from Part4 objects\n")
  }
}

# Filter metadata to match expression samples
if (!is.null(meta_data)) {
  shared_samples <- intersect(rownames(expr_matrix), rownames(meta_data))
  meta_data <- meta_data[shared_samples, ]
  cat("  Metadata:", nrow(meta_data), "samples\n")
}

# Load original modules for comparison (waterfall: cluster → local)
orig_modules <- NULL
for (opath in c("Part5_WGCNA/data/all_gene_module_assignments.tsv",
                "results/02_rnaseq/Part5_WGCNA/data/all_gene_module_assignments.tsv")) {
  if (file.exists(opath)) {
    orig_modules <- read.delim(opath)
    cat("  Original modules from:", opath, "\n")
    break
  }
}
if (is.null(orig_modules)) stop("Cannot find original module assignments.")


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 1: Build Network with Bicor
# ══════════════════════════════════════════════════════════════════════════════

cat("\n--- Section 1: Building improved network ---\n")
cat("  Using bicor (biweight midcorrelation) for robustness\n")
cat("  Power:", best_power_final, "| networkType: signed | minModuleSize: 30\n\n")

net <- blockwiseModules(
  expr_matrix,
  power = best_power_final,
  networkType = "signed",
  TOMType = "signed",
  corType = "bicor",
  maxPOutliers = 0.05,
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = TRUE,
  saveTOMFileBase = file.path(outdir, "data", "improved_TOM"),
  verbose = 3
)

module_colors <- labels2colors(net$colors)
names(module_colors) <- colnames(expr_matrix)

n_modules <- length(unique(module_colors)) - 1
cat("\n  Modules detected:", n_modules, "+ grey\n")
cat("  Module sizes:\n")
mod_table <- sort(table(module_colors), decreasing = TRUE)
for (m in names(mod_table)) {
  cat(sprintf("    %-12s: %d genes\n", m, mod_table[m]))
}

# Save module assignments
mod_assign <- data.frame(
  Gene = colnames(expr_matrix),
  Module_Color = module_colors,
  Module_Number = net$colors,
  stringsAsFactors = FALSE
)
write.table(mod_assign, file.path(outdir, "data", "improved_module_assignments.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 2: Module Eigengenes
# ══════════════════════════════════════════════════════════════════════════════

cat("\n--- Section 2: Module Eigengenes ---\n")

MEs <- net$MEs
colnames(MEs) <- gsub("^ME", "", colnames(MEs))

# Convert numeric labels to color labels in ME columns
color_order <- labels2colors(sort(unique(net$colors)))
ME_colors <- data.frame(MEs)
colnames(ME_colors) <- paste0("ME", labels2colors(as.numeric(gsub("^ME", "", colnames(MEs)))))

write.table(ME_colors, file.path(outdir, "data", "improved_module_eigengenes.tsv"),
            sep = "\t", row.names = TRUE, quote = FALSE)


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 3: Module Comparison with Original
# ══════════════════════════════════════════════════════════════════════════════

cat("\n--- Section 3: Module Comparison with Original ---\n")

shared_genes <- intersect(names(module_colors), orig_modules$Gene)
cat("  Shared genes:", length(shared_genes), "\n")

orig_labels <- orig_modules$Module_Color[match(shared_genes, orig_modules$Gene)]
new_labels <- module_colors[shared_genes]

# Overlap table
overlap <- table(orig_labels, new_labels)
cat("  Overlap table dimensions:", nrow(overlap), "x", ncol(overlap), "\n")

write.table(as.data.frame.matrix(overlap),
            file.path(outdir, "tables", "S06_module_overlap_matrix.tsv"),
            sep = "\t", quote = FALSE)

# Overlap heatmap
if (nrow(overlap) > 1 && ncol(overlap) > 1) {
  png(file.path(outdir, "plots", "S04_module_overlap_heatmap.png"),
      width = 1200, height = 1000, res = 150)
  pheatmap(log2(overlap + 1),
           main = "Module Overlap: Original vs Improved WGCNA",
           color = colorRampPalette(c("white", "#2980B9", "#E74C3C"))(50),
           cluster_rows = FALSE, cluster_cols = FALSE,
           display_numbers = TRUE, number_format = "%d",
           fontsize_number = 7)
  dev.off()
  cat("  Saved: S04_module_overlap_heatmap\n")
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 4: Hub Genes and Module Membership
# ══════════════════════════════════════════════════════════════════════════════

cat("\n--- Section 4: Hub Genes and Module Membership ---\n")

# Calculate kME (module membership) using bicor
kME <- bicor(expr_matrix, ME_colors, use = "pairwise.complete.obs",
             maxPOutliers = 0.05)

# Build trait data if metadata available
if (!is.null(meta_data)) {
  trait_data <- data.frame(
    tail_control = as.numeric(meta_data$tissue == "tail" & meta_data$condition == "control"),
    tail_amputated = as.numeric(meta_data$tissue == "tail" & meta_data$condition == "amputated"),
    row.names = rownames(meta_data)
  )

  # Gene significance for tail amputation
  GS <- as.data.frame(bicor(expr_matrix,
                             trait_data[rownames(expr_matrix), "tail_amputated", drop = FALSE],
                             use = "pairwise.complete.obs", maxPOutliers = 0.05))
  colnames(GS) <- "GS_tail_amputated"
}

# Intramodular connectivity
adj <- adjacency(expr_matrix, power = best_power_final, type = "signed",
                 corFnc = "bicor", corOptions = list(maxPOutliers = 0.05))

all_hub_genes <- data.frame()
unique_mods <- unique(module_colors)
unique_mods <- unique_mods[unique_mods != "grey"]

for (mod in unique_mods) {
  mod_genes <- names(module_colors)[module_colors == mod]
  if (length(mod_genes) < 10) next

  # Intramodular connectivity
  mod_adj <- adj[mod_genes, mod_genes]
  kWithin <- rowSums(mod_adj) - 1

  # Module membership
  me_col <- paste0("ME", mod)
  if (me_col %in% colnames(kME)) {
    mm <- abs(kME[mod_genes, me_col])
  } else {
    mm <- rep(NA, length(mod_genes))
    names(mm) <- mod_genes
  }

  # Hub criteria: top 10% connectivity AND |kME| > 0.7
  k_thresh <- quantile(kWithin, 0.90)
  is_hub <- kWithin >= k_thresh & !is.na(mm) & mm > 0.7

  mod_hub_df <- data.frame(
    Gene = mod_genes,
    Module = mod,
    kWithin = round(kWithin, 4),
    kME = round(mm, 4),
    IsHub = is_hub,
    stringsAsFactors = FALSE
  )

  if (exists("GS")) {
    mod_hub_df$GS <- round(GS[mod_genes, 1], 4)
  }

  all_hub_genes <- rbind(all_hub_genes, mod_hub_df)
}

write.table(all_hub_genes, file.path(outdir, "data", "improved_hub_genes.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

n_hubs <- sum(all_hub_genes$IsHub)
cat("  Total hub genes:", n_hubs, "\n")
cat("  Hub genes per module:\n")
hub_counts <- all_hub_genes %>% filter(IsHub) %>% count(Module) %>% arrange(desc(n))
for (i in 1:nrow(hub_counts)) {
  cat(sprintf("    %-12s: %d hubs\n", hub_counts$Module[i], hub_counts$n[i]))
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 5: Wilcoxon ME Comparison (Tail Regeneration)
# ══════════════════════════════════════════════════════════════════════════════

cat("\n--- Section 5: Wilcoxon ME Comparison ---\n")

if (!is.null(meta_data)) {
  ctrl_samples <- rownames(meta_data)[meta_data$tissue == "tail" & meta_data$condition == "control"]
  amp_samples <- rownames(meta_data)[meta_data$tissue == "tail" & meta_data$condition == "amputated"]

  ctrl_samples <- intersect(ctrl_samples, rownames(ME_colors))
  amp_samples <- intersect(amp_samples, rownames(ME_colors))

  cat("  Control samples:", length(ctrl_samples), "\n")
  cat("  Amputated samples:", length(amp_samples), "\n")

  if (length(ctrl_samples) >= 2 && length(amp_samples) >= 2) {
    me_wilcox <- data.frame()
    for (me_col in colnames(ME_colors)) {
      mod <- gsub("^ME", "", me_col)
      if (mod == "grey") next

      ctrl_vals <- ME_colors[ctrl_samples, me_col]
      amp_vals <- ME_colors[amp_samples, me_col]

      wt <- wilcox.test(ctrl_vals, amp_vals, exact = FALSE)
      diff_mean <- mean(amp_vals) - mean(ctrl_vals)
      pooled_sd <- sqrt((var(ctrl_vals) * (length(ctrl_vals) - 1) +
                          var(amp_vals) * (length(amp_samples) - 1)) /
                         (length(ctrl_vals) + length(amp_samples) - 2))
      cohens_d <- ifelse(pooled_sd > 0, diff_mean / pooled_sd, 0)

      me_wilcox <- rbind(me_wilcox, data.frame(
        module = mod,
        ctrl_mean = round(mean(ctrl_vals), 4),
        amp_mean = round(mean(amp_vals), 4),
        diff = round(diff_mean, 4),
        cohens_d = round(cohens_d, 3),
        direction = ifelse(diff_mean > 0, "UP", "DOWN"),
        p_value = wt$p.value,
        stringsAsFactors = FALSE
      ))
    }

    me_wilcox$padj <- p.adjust(me_wilcox$p_value, method = "BH")
    me_wilcox <- me_wilcox %>% arrange(p_value)

    cat("\n  Module eigengene differences (tail ctrl vs amp):\n")
    for (i in 1:nrow(me_wilcox)) {
      r <- me_wilcox[i, ]
      sig <- ifelse(r$padj < 0.05, "***", ifelse(r$p_value < 0.1, "*", ""))
      cat(sprintf("    %-12s d=%.3f  direction=%-4s  p=%.4f  padj=%.4f %s\n",
                  r$module, r$cohens_d, r$direction, r$p_value, r$padj, sig))
    }

    write.table(me_wilcox, file.path(outdir, "tables", "S07_ME_wilcoxon_improved.tsv"),
                sep = "\t", row.names = FALSE, quote = FALSE)
  }
} else {
  cat("  No metadata available — skipping Wilcoxon comparison\n")
}


# ══════════════════════════════════════════════════════════════════════════════
# SAVE
# ══════════════════════════════════════════════════════════════════════════════

cat("\n--- Saving objects ---\n")

save(net, module_colors, ME_colors, kME, adj, all_hub_genes,
     best_filter, best_power_final,
     file = file.path(outdir, "data", "improved_WGCNA_objects.RData"))

cat("  Saved: improved_WGCNA_objects.RData\n")
cat("  TOM blocks saved by blockwiseModules\n")

cat("\n=== Script 02 complete ===\n")
cat("  → Run 03_methylation_network_integration.R next\n")
