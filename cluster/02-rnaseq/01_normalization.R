###############################################################################
# PART 1: Setup, Normalization & Quality Control
###############################################################################

library(dplyr)
library(stringr)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

# ── Output directories ───────────────────────────────────────────────────────
dir.create("Part1_QC/plots", recursive = TRUE, showWarnings = FALSE)
dir.create("Part1_QC/data",  recursive = TRUE, showWarnings = FALSE)

# ── 1. Metadata + DESeq2 setup ───────────────────────────────────────────────
dir_counts <- "/mnt/data/alfredvar/jmiranda/20-Transcriptomic_Bulk/25-metaAnalysisTranscriptome/counts_HTseq_EviAnn/"
files <- list.files(dir_counts, pattern = "htseq_gene_counts\\.txt$", full.names = TRUE)
stopifnot(length(files) > 0)

bn <- basename(files)
df <- data.frame(
  file      = files,
  fileName  = bn,
  sample_id = tools::file_path_sans_ext(bn),
  stringsAsFactors = FALSE
)

infer_row <- function(b) {
  experiment <- NA_character_; tissue <- NA_character_; condition <- NA_character_
  if (str_starts(b, "T\\d+S\\d+")) {
    experiment <- "tail_amputation"; tissue <- "tail"; condition <- "amputated"
  } else if (str_starts(b, "C\\d+S\\d+")) {
    experiment <- "tail_amputation"; tissue <- "tail"; condition <- "control"
  } else if (str_starts(b, "R\\d+_")) {
    experiment <- "eye_amputation"; tissue <- "eye"; condition <- "amputated"
  } else if (str_starts(b, "C\\d+_")) {
    experiment <- "eye_amputation"; tissue <- "eye"; condition <- "control"
  } else if (str_starts(b, "irrep")) {
    experiment <- "bodywall_irradiation"; tissue <- "bodywall"; condition <- "irradiated"
  } else if (str_starts(b, "dcrep")) {
    experiment <- "bodywall_irradiation"; tissue <- "bodywall"; condition <- "control"
  } else if (str_starts(b, "fungicide_l30")) {
    experiment <- "bodywall_fungicide"; tissue <- "bodywall"; condition <- "fungicide"
  } else if (str_starts(b, "fungicide_l0")) {
    experiment <- "bodywall_fungicide"; tissue <- "bodywall"; condition <- "control"
  } else if (str_detect(b, "rmoverrep_fpDLHead")) {
    experiment <- "baseline_tissues"; tissue <- "head"; condition <- "control"
  } else if (str_detect(b, "rmoverrep_fpDLJuv")) {
    experiment <- "baseline_tissues"; tissue <- "juvenile"; condition <- "control"
  } else if (str_detect(b, "rmoverrep_fpDLOvo")) {
    experiment <- "baseline_tissues"; tissue <- "ovotestis"; condition <- "control"
  }
  list(experiment = experiment, tissue = tissue, condition = condition)
}

parsed <- lapply(bn, infer_row)
parsed <- do.call(rbind, lapply(parsed, as.data.frame, stringsAsFactors = FALSE))
df$experiment <- parsed$experiment
df$tissue     <- parsed$tissue
df$condition  <- parsed$condition

df$experiment <- factor(df$experiment,
                        levels = c("baseline_tissues","tail_amputation","eye_amputation",
                                   "bodywall_irradiation","bodywall_fungicide"))
df$condition <- factor(df$condition,
                       levels = c("control","irradiated","amputated","fungicide"))

df <- df %>% mutate(
  sample_id = str_remove(sample_id, "_htseq_gene_counts"),
  sample_id = str_remove(sample_id, ".Aligned.out.bam"),
  sample_id = str_remove(sample_id, "rmoverrep_fpDL")
)
rownames(df) <- df$sample_id

cat("\n══ CHECKPOINT 1A: Verify metadata ══\n")
print(table(df$experiment, df$condition))
print(table(df$experiment, df$tissue))
cat("Total samples:", nrow(df), "\n")
cat("NAs in experiment?", any(is.na(df$experiment)), "\n")
cat("NAs in condition?", any(is.na(df$condition)), "\n\n")

# ── 2. DESeq2 object ─────────────────────────────────────────────────────────
sampleTable <- df[, c("sample_id", "fileName", "condition", "tissue", "experiment")]

dds <- DESeqDataSetFromHTSeqCount(
  sampleTable = sampleTable,
  directory    = dir_counts,
  design       = ~ experiment + condition
)

cat("Genes before filtering:", nrow(dds), "\n")
dds <- dds[rowSums(counts(dds)) > 1, ]
cat("After rowSum > 1:", nrow(dds), "\n")
keep <- rowSums(counts(dds) > 5) >= 0.25 * ncol(dds)
dds <- dds[keep, ]
cat("After >5 in 25% samples:", nrow(dds), "\n\n")

# ── 3. VST normalization ─────────────────────────────────────────────────────
vsd <- vst(dds, blind = TRUE)
normalized_counts <- assay(vsd)

# ── 4. PCA plots ─────────────────────────────────────────────────────────────
pca_data <- plotPCA(vsd, intgroup = c("experiment", "condition"), returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

p1 <- ggplot(pca_data, aes(x = PC1, y = PC2, color = experiment, shape = condition)) +
  geom_point(size = 4, alpha = 0.8) +
  labs(title = "PCA: by Experiment (Batch)",
       x = paste0("PC1: ", percent_var[1], "% variance"),
       y = paste0("PC2: ", percent_var[2], "% variance")) +
  theme_minimal(base_size = 12) + theme(plot.title = element_text(face = "bold"))
ggsave("Part1_QC/plots/PCA_by_experiment.png", p1, width = 10, height = 7, dpi = 300)
ggsave("Part1_QC/plots/PCA_by_experiment.pdf", p1, width = 10, height = 7)

p2 <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, shape = experiment)) +
  geom_point(size = 4, alpha = 0.8) +
  labs(title = "PCA: by Condition",
       x = paste0("PC1: ", percent_var[1], "% variance"),
       y = paste0("PC2: ", percent_var[2], "% variance")) +
  theme_minimal(base_size = 12) + theme(plot.title = element_text(face = "bold"))
ggsave("Part1_QC/plots/PCA_by_condition.png", p2, width = 10, height = 7, dpi = 300)
ggsave("Part1_QC/plots/PCA_by_condition.pdf", p2, width = 10, height = 7)

p3 <- ggplot(pca_data, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 4, alpha = 0.8) +
  labs(title = "PCA: by Experiment:Condition",
       x = paste0("PC1: ", percent_var[1], "% variance"),
       y = paste0("PC2: ", percent_var[2], "% variance")) +
  theme_minimal(base_size = 12) + theme(plot.title = element_text(face = "bold"))
ggsave("Part1_QC/plots/PCA_by_group.png", p3, width = 12, height = 7, dpi = 300)
ggsave("Part1_QC/plots/PCA_by_group.pdf", p3, width = 12, height = 7)

# ── 5. Sample distance heatmap ───────────────────────────────────────────────
sample_dists <- dist(t(normalized_counts))
sample_dist_matrix <- as.matrix(sample_dists)
annotation_df <- data.frame(
  Experiment = colData(vsd)$experiment,
  Condition  = colData(vsd)$condition,
  Tissue     = colData(vsd)$tissue,
  row.names  = colnames(vsd)
)

pheatmap(sample_dist_matrix,
         clustering_distance_rows = sample_dists,
         clustering_distance_cols = sample_dists,
         annotation_row = annotation_df, annotation_col = annotation_df,
         color = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
         show_rownames = TRUE, show_colnames = FALSE,
         main = "Sample-to-Sample Distance (VST)", fontsize = 8,
         filename = "Part1_QC/plots/sample_distance_heatmap.png",
         width = 12, height = 10)

pheatmap(sample_dist_matrix,
         clustering_distance_rows = sample_dists,
         clustering_distance_cols = sample_dists,
         annotation_row = annotation_df, annotation_col = annotation_df,
         color = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
         show_rownames = TRUE, show_colnames = FALSE,
         main = "Sample-to-Sample Distance (VST)", fontsize = 8,
         filename = "Part1_QC/plots/sample_distance_heatmap.pdf",
         width = 12, height = 10)

# ── Save data ────────────────────────────────────────────────────────────────
write.csv(df, "Part1_QC/data/sample_metadata.csv", row.names = FALSE)
write.table(normalized_counts, "Part1_QC/data/normalized_counts_VST.tsv",
            sep = "\t", quote = FALSE)
save(dds, vsd, normalized_counts, df, file = "Part1_QC/data/Part1_objects.RData")

cat("\n══ CHECKPOINT 1B ══\n")
cat("Check Part1_QC/plots/ — if experiments dominate PCA → run Part1b\n")
cat("═══ Part 1 complete ═══\n")