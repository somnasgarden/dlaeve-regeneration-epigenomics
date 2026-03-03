#!/usr/bin/env Rscript
# =============================================================================
# Part 1: Setup, Data Loading, Quality Control & Exploratory Data Analysis
# =============================================================================
# Adapted for Deroceras laeve genome (BSgenome.Dlaeve.NCBI.dlgm)
# Samples: C1, C2 (Control) vs A1, A2 (Treatment)
# Input: Bismark CpG_report.txt.gz files
# =============================================================================

# ── 1. Install & Load Packages ───────────────────────────────────────────────

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

essential_packages <- c(
  "bsseq",            # Core BS-seq analysis
  "DSS",              # Differential methylation analysis
  "GenomicRanges",    # Genomic interval operations
  "rtracklayer",      # GFF/GTF import
  "data.table",       # Fast data manipulation
  "ggplot2",          # Plotting
  "pheatmap",         # Heatmaps
  "RColorBrewer",     # Color palettes
  "dplyr",            # Data manipulation
  "tidyr",            # Data reshaping
  "scales",           # Scale functions for ggplot2
  "gridExtra"         # Multiple plots arrangement
)

# Install missing packages
for (pkg in essential_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  }
}

# Load all libraries
suppressPackageStartupMessages({
  library(bsseq)
  library(DSS)
  library(GenomicRanges)
  library(rtracklayer)
  library(data.table)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(gridExtra)
  library(BSgenome.Dlaeve.NCBI.dlgm)
})

# Set global options
options(stringsAsFactors = FALSE)
options(scipen = 999)

# Configure ggplot2 theme
theme_set(theme_minimal() +
            theme(text = element_text(size = 12),
                  plot.title = element_text(size = 14, face = "bold"),
                  axis.title = element_text(size = 12),
                  legend.text = element_text(size = 10)))

cat("All packages loaded successfully.\n")

# ── 2. Define Paths & Parameters ─────────────────────────────────────────────

BASE_DIR <- "/mnt/data/alfredvar/"

GFF_FILE <- file.path(BASE_DIR,
                      "30-Genoma/31-Alternative_Annotation_EviAnn/derLaeGenome_namesDlasi_v2.fasta.functional_note.pseudo_label.gff")

METH_SAMPLES <- c(
  C1 = file.path(BASE_DIR, "jmiranda/50-Genoma/51-Metilacion/09_methylation_calls/C1.CpG_report.txt.gz"),
  C2 = file.path(BASE_DIR, "jmiranda/50-Genoma/51-Metilacion/09_methylation_calls/C2.CpG_report.txt.gz"),
  A1 = file.path(BASE_DIR, "jmiranda/50-Genoma/51-Metilacion/09_methylation_calls/A1.CpG_report.txt.gz"),
  A2 = file.path(BASE_DIR, "jmiranda/50-Genoma/51-Metilacion/09_methylation_calls/A2.CpG_report.txt.gz")
)

# *** FIX 1: Use an absolute path directly, NOT file.path(BASE_DIR, ...) ***
RESULTS_DIR <- "/mnt/data/alfredvar/rlopezt/feb_w2_scripts/Preliminary/"
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

# Analysis parameters
analysis_params <- list(
  min_coverage  = 5,     # Minimum read coverage per CpG site
  min_samples   = 2,     # Minimum samples with adequate coverage
  min_diff      = 0.1,   # Minimum methylation difference (10%)
  fdr_threshold = 0.05,  # False discovery rate threshold
  min_dmr_cpgs  = 3,     # Minimum CpGs required for DMR calling
  min_dmr_length = 50    # Minimum DMR length in base pairs
)

# Define sample metadata
sample_info <- data.frame(
  sample_id = c("C1", "C2", "A1", "A2"),
  condition = c("Control", "Control", "Amputated", "Amputated"),
  file_path = unname(METH_SAMPLES),
  stringsAsFactors = FALSE
)

cat("Parameters and paths configured.\n")
print(sample_info[, c("sample_id", "condition")])

# ── 3. Genome Info ───────────────────────────────────────────────────────────

genome <- BSgenome.Dlaeve.NCBI.dlgm

cat("\nFirst 40 sequence names in genome:\n")
print(head(seqnames(genome), 40))

keep_chr <- paste0("chr", 1:31)
keep_chr <- intersect(keep_chr, seqnames(genome))
genome_length <- sum(seqlengths(genome)[keep_chr], na.rm = TRUE)
cat(sprintf("  Genome size (filtered chroms): %s bp\n", format(genome_length, big.mark = ",")))
cat(sprintf("  Chromosomes kept: %d\n", length(keep_chr)))

# ── 4. Load GFF Annotation ──────────────────────────────────────────────────

cat("\nLoading GFF annotation...\n")
genes_gff <- import(GFF_FILE)
genes_gff <- genes_gff[seqnames(genes_gff) %in% keep_chr]

genes <- genes_gff[genes_gff$type == "gene"]
cat(sprintf("  Total genes: %s\n", format(length(genes), big.mark = ",")))

# Classify gene types
gene_ids  <- mcols(genes)$ID
gene_type <- ifelse(grepl("XLOC.*lncRNA", gene_ids), "lncRNA",
                    ifelse(grepl("^LOC_", gene_ids), "protein_coding", "other"))
cat("  Gene type breakdown:\n")
print(table(gene_type))

# ── 5. Load Bismark CpG Report Files ────────────────────────────────────────

# Bismark CpG_report.txt format (7 columns):
# chr | pos | strand | count_methylated | count_unmethylated | C-context | trinucleotide

read_cpg_report <- function(file_path, sample_name, min_cov = 5, chroms = NULL) {
  cat(sprintf("  Loading %s from: %s\n", sample_name, basename(file_path)))
  
  # *** FIX 2: Read all 7 columns, then subset afterwards ***
  cov_data <- fread(file_path,
                    col.names = c("chr", "pos", "strand",
                                  "meth_count", "unmeth_count",
                                  "context", "trinuc"))
  
  # Drop columns we don't need
  cov_data[, c("strand", "context", "trinuc") := NULL]
  
  # Calculate total coverage
  cov_data[, total_count := meth_count + unmeth_count]
  
  # Filter: minimum coverage + target chromosomes
  cov_data <- cov_data[total_count >= min_cov]
  if (!is.null(chroms)) {
    cov_data <- cov_data[chr %in% chroms]
  }
  
  # Calculate methylation percentage
  cov_data[, methylation := meth_count / total_count * 100]
  
  cov_data$sample <- sample_name
  
  cat(sprintf("    -> %s CpG sites after filtering (cov >= %d)\n",
              format(nrow(cov_data), big.mark = ","), min_cov))
  
  return(cov_data)
}

# Load all samples
cat("\nLoading methylation data...\n")
sample_data_list <- list()

for (i in 1:nrow(sample_info)) {
  sample_data_list[[sample_info$sample_id[i]]] <-
    read_cpg_report(sample_info$file_path[i],
                    sample_info$sample_id[i],
                    analysis_params$min_coverage,
                    keep_chr)
}

# ── 6. Find Common CpG Sites & Build BSseq Object ───────────────────────────

cat("\nFinding CpG sites common to all samples...\n")

get_site_id <- function(data) paste(data$chr, data$pos, sep = ":")
site_lists   <- lapply(sample_data_list, get_site_id)
common_sites <- Reduce(intersect, site_lists)
cat(sprintf("  Common CpG sites across all 4 samples: %s\n",
            format(length(common_sites), big.mark = ",")))

# Filter each sample to common sites
filter_to_common <- function(data, common) {
  sid <- get_site_id(data)
  data[sid %in% common]
}
sample_data_filtered <- lapply(sample_data_list, filter_to_common, common_sites)

# Sort identically
for (nm in names(sample_data_filtered)) {
  sample_data_filtered[[nm]] <- sample_data_filtered[[nm]][order(chr, pos)]
}

# Build BSseq object
ref_data <- sample_data_filtered[[1]]
gr <- GRanges(seqnames = ref_data$chr,
              ranges   = IRanges(start = ref_data$pos, width = 1))

sample_names <- sample_info$sample_id

M_matrix   <- matrix(nrow = nrow(ref_data), ncol = length(sample_names))
Cov_matrix <- matrix(nrow = nrow(ref_data), ncol = length(sample_names))

for (i in seq_along(sample_names)) {
  sd <- sample_data_filtered[[sample_names[i]]]
  M_matrix[, i]   <- sd$meth_count
  Cov_matrix[, i] <- sd$total_count
}

colnames(M_matrix) <- colnames(Cov_matrix) <- sample_names

bs_obj <- BSseq(gr = gr, M = M_matrix, Cov = Cov_matrix)
pData(bs_obj) <- sample_info[match(sample_names, sample_info$sample_id), ]
rownames(pData(bs_obj)) <- sample_names

cat(sprintf("  BSseq object created: %s sites x %s samples\n",
            format(nrow(bs_obj), big.mark = ","), ncol(bs_obj)))

# Save BSseq object
saveRDS(bs_obj, file.path(RESULTS_DIR, "bsseq_object.rds"))
cat("  BSseq object saved.\n")

# ── 7. Quality Control: Global Methylation Statistics ────────────────────────

cat("\n--- Quality Control ---\n")

calculate_methylation_stats <- function(bs_obj) {
  beta_matrix <- getMeth(bs_obj, type = "raw")
  stats_df <- data.frame()
  
  for (i in 1:ncol(bs_obj)) {
    beta_vals <- beta_matrix[, i]
    beta_vals <- beta_vals[!is.nan(beta_vals)]
    
    sample_stats <- data.frame(
      sample             = colnames(bs_obj)[i],
      condition          = pData(bs_obj)$condition[i],
      total_sites        = length(beta_vals),
      mean_methylation   = mean(beta_vals, na.rm = TRUE),
      median_methylation = median(beta_vals, na.rm = TRUE),
      sd_methylation     = sd(beta_vals, na.rm = TRUE),
      sites_high_meth    = sum(beta_vals > 0.8, na.rm = TRUE),
      sites_low_meth     = sum(beta_vals < 0.2, na.rm = TRUE),
      sites_intermediate = sum(beta_vals >= 0.2 & beta_vals <= 0.8, na.rm = TRUE)
    )
    stats_df <- rbind(stats_df, sample_stats)
  }
  return(stats_df)
}

methylation_stats <- calculate_methylation_stats(bs_obj)
write.table(methylation_stats,
            file.path(RESULTS_DIR, "global_methylation_stats.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("Global methylation statistics:\n")
print(methylation_stats)

# ── 8. Sample Correlation Heatmap ────────────────────────────────────────────

beta_matrix    <- getMeth(bs_obj, type = "raw")
complete_sites <- complete.cases(beta_matrix)
beta_complete  <- beta_matrix[complete_sites, ]

sample_cors <- cor(beta_complete)

annotation_df <- data.frame(
  Condition = pData(bs_obj)$condition,
  row.names = colnames(sample_cors)
)

pdf(file.path(RESULTS_DIR, "sample_correlation_heatmap.pdf"), width = 8, height = 6)
pheatmap(sample_cors,
         annotation_col = annotation_df,
         annotation_row = annotation_df,
         color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
         main = "Sample-to-Sample Methylation Correlations",
         fontsize = 12,
         display_numbers = TRUE,
         number_format = "%.3f")
dev.off()
cat("Correlation heatmap saved.\n")

# ── 9. PCA Analysis ─────────────────────────────────────────────────────────

pca_result   <- prcomp(t(beta_complete), center = TRUE, scale. = FALSE)
var_explained <- summary(pca_result)$importance[2, 1:min(4, ncol(beta_complete))] * 100

pca_data <- data.frame(
  PC1       = pca_result$x[, 1],
  PC2       = pca_result$x[, 2],
  Sample    = rownames(pca_result$x),
  Condition = pData(bs_obj)$condition
)

p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Condition, label = Sample)) +
  geom_point(size = 4) +
  geom_text(vjust = -1, hjust = 0.5, size = 3) +
  scale_color_manual(values = c("Control" = "#2E86C1", "Treatment" = "#E74C3C")) +
  labs(
    title = "PCA of DNA Methylation Profiles",
    x = sprintf("PC1 (%.1f%%)", var_explained[1]),
    y = sprintf("PC2 (%.1f%%)", var_explained[2])
  ) +
  theme_minimal()

ggsave(file.path(RESULTS_DIR, "pca_analysis.pdf"), p_pca, width = 10, height = 8)
cat("PCA plot saved.\n")

# ── 10. Methylation Distribution Plots ───────────────────────────────────────

beta_long <- as.data.frame(beta_matrix)
colnames(beta_long) <- sample_names
beta_long <- tidyr::pivot_longer(beta_long, cols = everything(),
                                 names_to = "Sample", values_to = "Methylation")
beta_long <- beta_long[!is.na(beta_long$Methylation), ]
beta_long$Condition <- sample_info$condition[match(beta_long$Sample, sample_info$sample_id)]

p_density <- ggplot(beta_long, aes(x = Methylation, color = Sample, linetype = Condition)) +
  geom_density(linewidth = 0.8) +
  scale_color_manual(values = c("C1" = "#2E86C1", "C2" = "#5DADE2",
                                "A1" = "#E74C3C", "A2" = "#F1948A")) +
  labs(title = "CpG Methylation Distribution by Sample",
       x = "Methylation Level (beta)", y = "Density") +
  theme_minimal()

ggsave(file.path(RESULTS_DIR, "methylation_density.pdf"), p_density, width = 10, height = 6)
cat("Methylation density plot saved.\n")

cat("\n========================================\n")
cat("Part 1 complete! Outputs in:", RESULTS_DIR, "\n")
cat("Files created:\n")
cat("  - bsseq_object.rds\n")
cat("  - global_methylation_stats.txt\n")
cat("  - sample_correlation_heatmap.pdf\n")
cat("  - pca_analysis.pdf\n")
cat("  - methylation_density.pdf\n")
cat("========================================\n")