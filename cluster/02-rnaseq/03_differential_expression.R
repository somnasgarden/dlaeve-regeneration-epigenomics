###############################################################################
# PART 2: Differential Expression & Scaling
###############################################################################

# UNCOMMENT the correct one:
# load("Part1_QC/data/Part1_objects.RData")          # if you SKIPPED batch correction
load("Part1b_BatchCorrection/data/Part1b_objects.RData")  # if you RAN batch correction

library(DESeq2)
library(dplyr)
library(ggplot2)
library(pheatmap)

# ── Output directories ───────────────────────────────────────────────────────
dir.create("Part2_DEGs/plots", recursive = TRUE, showWarnings = FALSE)
dir.create("Part2_DEGs/data",  recursive = TRUE, showWarnings = FALSE)

# ── 1. Run DESeq2 ────────────────────────────────────────────────────────────
cat("Running DESeq2...\n")
dds <- DESeq(dds)
cat("Available coefficients:\n")
print(resultsNames(dds))
cat("\n")

# ── 2. Extract DEGs per condition vs control ─────────────────────────────────
conditions_to_test <- levels(colData(dds)$condition)
conditions_to_test <- conditions_to_test[conditions_to_test != "control"]

all_degs <- list()
all_results <- list()

for (cond in conditions_to_test) {
  contrast_name <- paste0(cond, "_vs_control")
  res <- results(dds, contrast = c("condition", cond, "control"))
  res <- res[!is.na(res$padj), ]
  sig <- res[res$padj < 0.05 & abs(res$log2FoldChange) > 1, ]
  all_degs[[contrast_name]] <- rownames(sig)
  all_results[[contrast_name]] <- as.data.frame(res)
  cat(sprintf("  %s: %d DEGs\n", contrast_name, nrow(sig)))
}

# ── 3. Union of all DEGs ────────────────────────────────────────────────────
deg_union <- unique(unlist(all_degs))
cat(sprintf("\nTotal unique DEGs (union): %d\n\n", length(deg_union)))

# ── CHECKPOINT 2A ────────────────────────────────────────────────────────────
cat("══ CHECKPOINT 2A: DEG numbers ══\n")
if (length(deg_union) < 100) {
  cat("⚠ Very few DEGs. Consider relaxing: padj<0.1, |log2FC|>0.5\n\n")
} else if (length(deg_union) > 5000) {
  cat("⚠ Many DEGs. Consider tightening: |log2FC|>2\n\n")
} else {
  cat("✓ DEG count looks good for clustering.\n\n")
}

# ── 4. Scale ─────────────────────────────────────────────────────────────────
deg_expression <- normalized_counts[deg_union, ]
deg_scaled <- t(scale(t(deg_expression)))

nan_genes <- rowSums(is.nan(deg_scaled)) > 0
if (any(nan_genes)) {
  cat("Removing", sum(nan_genes), "zero-variance genes\n")
  deg_scaled <- deg_scaled[!nan_genes, ]
  deg_expression <- deg_expression[!nan_genes, ]
}

cat("Scaled matrix:", nrow(deg_scaled), "genes x", ncol(deg_scaled), "samples\n\n")

# ── 5. Heatmap of all DEGs ──────────────────────────────────────────────────
annotation_col <- data.frame(
  Experiment = colData(dds)$experiment,
  Condition  = colData(dds)$condition,
  Tissue     = colData(dds)$tissue,
  row.names  = colnames(dds)
)

pheatmap(deg_scaled,
         cluster_rows = TRUE, cluster_cols = TRUE,
         annotation_col = annotation_col,
         show_rownames = FALSE, show_colnames = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = paste0("All DEGs (n=", nrow(deg_scaled), ") - Unsorted"),
         fontsize = 8, fontsize_col = 6,
         filename = "Part2_DEGs/plots/all_DEGs_heatmap.png",
         width = 14, height = 10)

pheatmap(deg_scaled,
         cluster_rows = TRUE, cluster_cols = TRUE,
         annotation_col = annotation_col,
         show_rownames = FALSE, show_colnames = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = paste0("All DEGs (n=", nrow(deg_scaled), ") - Unsorted"),
         fontsize = 8, fontsize_col = 6,
         filename = "Part2_DEGs/plots/all_DEGs_heatmap.pdf",
         width = 14, height = 10)

# ── Save data ────────────────────────────────────────────────────────────────
for (name in names(all_results)) {
  write.csv(all_results[[name]], paste0("Part2_DEGs/data/DEresults_", name, ".csv"))
}
write.csv(data.frame(Gene = deg_union), "Part2_DEGs/data/DEG_union_list.csv", row.names = FALSE)
write.table(deg_scaled, "Part2_DEGs/data/DEG_scaled_expression.tsv", sep = "\t", quote = FALSE)

save(dds, vsd, normalized_counts, df,
     deg_union, deg_expression, deg_scaled, annotation_col, all_degs,
     file = "Part2_DEGs/data/Part2_objects.RData")

cat("══ CHECKPOINT 2B ══\n")
cat("Check Part2_DEGs/plots/all_DEGs_heatmap.png — do you see structure?\n")
cat("═══ Part 2 complete ═══\n")