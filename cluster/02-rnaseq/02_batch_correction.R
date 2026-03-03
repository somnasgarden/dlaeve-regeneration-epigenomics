###############################################################################
# PART 1b: Batch Correction (ONLY if PCA shows experiment dominates)
###############################################################################

load("Part1_QC/data/Part1_objects.RData")

library(DESeq2)
library(limma)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

# ── Output directories ───────────────────────────────────────────────────────
dir.create("Part1b_BatchCorrection/plots", recursive = TRUE, showWarnings = FALSE)
dir.create("Part1b_BatchCorrection/data",  recursive = TRUE, showWarnings = FALSE)

# ── 1. Batch correction ─────────────────────────────────────────────────────
cat("Removing: experiment effect | Preserving: condition effect\n\n")

design_preserve <- model.matrix(~ condition, data = colData(dds))

normalized_corrected <- removeBatchEffect(
  assay(vsd),
  batch  = colData(dds)$experiment,
  design = design_preserve
)

# ── 2. PCA on corrected data ────────────────────────────────────────────────
pca_result <- prcomp(t(normalized_corrected), center = TRUE, scale. = FALSE)
percent_var <- round(100 * (pca_result$sdev^2 / sum(pca_result$sdev^2)), 1)

pca_df <- data.frame(
  PC1        = pca_result$x[, 1],
  PC2        = pca_result$x[, 2],
  experiment = colData(dds)$experiment,
  condition  = colData(dds)$condition,
  tissue     = colData(dds)$tissue,
  sample     = colnames(normalized_corrected)
)

p1 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = experiment, shape = condition)) +
  geom_point(size = 4, alpha = 0.8) +
  labs(title = "PCA AFTER Batch Correction: by Experiment",
       subtitle = "Experiments should now be mixed",
       x = paste0("PC1: ", percent_var[1], "% variance"),
       y = paste0("PC2: ", percent_var[2], "% variance")) +
  theme_minimal(base_size = 12) + theme(plot.title = element_text(face = "bold"))
ggsave("Part1b_BatchCorrection/plots/PCA_corrected_by_experiment.png", p1, width = 10, height = 7, dpi = 300)
ggsave("Part1b_BatchCorrection/plots/PCA_corrected_by_experiment.pdf", p1, width = 10, height = 7)

p2 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = condition, shape = experiment)) +
  geom_point(size = 4, alpha = 0.8) +
  labs(title = "PCA AFTER Batch Correction: by Condition",
       subtitle = "Biological signal should be clearer",
       x = paste0("PC1: ", percent_var[1], "% variance"),
       y = paste0("PC2: ", percent_var[2], "% variance")) +
  theme_minimal(base_size = 12) + theme(plot.title = element_text(face = "bold"))
ggsave("Part1b_BatchCorrection/plots/PCA_corrected_by_condition.png", p2, width = 10, height = 7, dpi = 300)
ggsave("Part1b_BatchCorrection/plots/PCA_corrected_by_condition.pdf", p2, width = 10, height = 7)

# ── Save ─────────────────────────────────────────────────────────────────────
normalized_counts_uncorrected <- assay(vsd)
normalized_counts <- normalized_corrected

write.table(normalized_counts, "Part1b_BatchCorrection/data/normalized_counts_corrected.tsv",
            sep = "\t", quote = FALSE)

save(dds, vsd, normalized_counts, normalized_counts_uncorrected, df,
     file = "Part1b_BatchCorrection/data/Part1b_objects.RData")

cat("\nCompare Part1_QC/plots/ vs Part1b_BatchCorrection/plots/\n")
cat("═══ Part 1b complete ═══\n")