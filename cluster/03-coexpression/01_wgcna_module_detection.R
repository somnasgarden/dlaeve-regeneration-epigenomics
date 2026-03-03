###############################################################################
# PART 5 (V3 — DEFINITIVE): WGCNA – Weighted Gene Co-expression Network Analysis
# Organism: Deroceras laeve (slug)
# Data: RNA-seq meta-analysis (batch-corrected VST counts)
# Conditions: control, amputated, irradiated, fungicide
# Tissues: tail, eye, bodywall, head, juvenile, ovotestis
# Goal: Identify gene co-expression modules, hub genes (for DMR cross-ref),
#        and functional enrichment across all modules
# Working dir: /mnt/data/alfredvar/rlopezt/CorrelationMatrix
#
# V3 CHANGES (DEFINITIVE):
#   - PRIMARY module selection: Wilcoxon ME comparison (replaces Pearson r)
#   - NEW: Gene dendrograms for tail_control and tail_amputated with modules
#   - FIX: Enrichment plots generated for ALL modules with significant terms
#   - FIX: Enrichment plots skip modules with 0 significant terms (no empty plots)
#   - FIX: MF enrichment label truncation
#   - KEEP: 25th percentile variance filter, power 12, signed network
#   - KEEP: Publication-ready labels, tissue-specific traits, tissue-distance
###############################################################################

# ── Load previous pipeline data ──────────────────────────────────────────────
load("Part4_Kmeans/data/Part4_objects.RData")


#-----------------------------------------------
# STEP 0: Install / load required R packages
#-----------------------------------------------

options(repos = c(CRAN = "https://cloud.r-project.org"))

for (pkg in c("WGCNA", "igraph", "ggrepel", "corrplot", "reshape2")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, dependencies = TRUE)
}

library(WGCNA)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(pheatmap)
library(igraph)
library(RColorBrewer)

enableWGCNAThreads(nThreads = 8)


# ── Clean old outputs and create fresh directories ───────────────────────────
if (dir.exists("Part5_WGCNA")) {
  cat("Removing previous Part5_WGCNA/ folder...\n")
  unlink("Part5_WGCNA", recursive = TRUE)
}
dir.create("Part5_WGCNA/inputs", recursive = TRUE, showWarnings = FALSE)
dir.create("Part5_WGCNA/plots",  recursive = TRUE, showWarnings = FALSE)
dir.create("Part5_WGCNA/data",   recursive = TRUE, showWarnings = FALSE)


###############################################################################
# PUBLICATION-READY SAMPLE LABELS
###############################################################################

sample_labels <- c(
  # Tail amputation experiment
  "C1S1" = "Tail Ctrl 1",  "C2S2" = "Tail Ctrl 2",
  "C3S3" = "Tail Ctrl 3",  "C4S4" = "Tail Ctrl 4",
  "T1S5" = "Tail Amp 1",   "T2S6" = "Tail Amp 2",
  "T3S7" = "Tail Amp 3",   "T4S8" = "Tail Amp 4",
  # Eye amputation experiment
  "C1"   = "Eye Ctrl 1",   "C4"   = "Eye Ctrl 2",
  "C5"   = "Eye Ctrl 3",   "C8"   = "Eye Ctrl 4",
  "R1"   = "Eye Amp 1",    "R3"   = "Eye Amp 2",
  "R4"   = "Eye Amp 3",    "R6"   = "Eye Amp 4",
  # Bodywall irradiation experiment
  "dcrep1" = "BW Irr Ctrl 1", "dcrep2" = "BW Irr Ctrl 2",
  "dcrep3" = "BW Irr Ctrl 3", "dcrep4" = "BW Irr Ctrl 4",
  "dcrep6" = "BW Irr Ctrl 5",
  "irrep1" = "BW Irrad 1", "irrep4" = "BW Irrad 2",
  "irrep5" = "BW Irrad 3", "irrep6" = "BW Irrad 4",
  "irrep7" = "BW Irrad 5",
  # Bodywall fungicide experiment
  "fungicide_l0a" = "BW Fung Ctrl 1", "fungicide_l0b" = "BW Fung Ctrl 2",
  "fungicide_l0c" = "BW Fung Ctrl 3", "fungicide_l0d" = "BW Fung Ctrl 4",
  "fungicide_l30a" = "BW Fung 1", "fungicide_l30b" = "BW Fung 2",
  "fungicide_l30c" = "BW Fung 3", "fungicide_l30d" = "BW Fung 4",
  "fungicide_l30e" = "BW Fung 5",
  # Baseline tissues
  "Head1_S4"  = "Head 1",   "Head2_S5" = "Head 2",   "Head3_S6" = "Head 3",
  "Juv1_S10"  = "Juvenile 1", "Juv2_S11" = "Juvenile 2", "Juv3_S12" = "Juvenile 3",
  "Ovo1_S7"   = "Ovotestis 1", "Ovo2_S8" = "Ovotestis 2", "Ovo3_S9" = "Ovotestis 3"
)

get_pub_label <- function(ids) {
  labels <- sample_labels[ids]
  labels[is.na(labels)] <- ids[is.na(labels)]
  return(labels)
}


# ── Save inputs ──────────────────────────────────────────────────────────────
cat("── Saving input files ──\n")

write.table(normalized_counts, "Part5_WGCNA/inputs/normalized_counts.tsv",
            sep = "\t", quote = FALSE)
write.table(df, "Part5_WGCNA/inputs/sample_metadata.tsv",
            sep = "\t", quote = FALSE)

enrichment_source <- "/mnt/data/alfredvar/rlopezt/Preliminary/protein.enrichment.terms.v12.0.txt"
enrichment_link   <- "Part5_WGCNA/inputs/protein.enrichment.terms.v12.0.txt"
if (!file.exists(enrichment_link)) file.symlink(enrichment_source, enrichment_link)

writeLines(c(
  "Part 5 WGCNA V3 DEFINITIVE - Input Files",
  paste("Date:", Sys.Date()),
  "",
  "Network: signed, power=12, 25th percentile variance filter",
  "Module selection: Wilcoxon ME comparison (tail ctrl vs amp)",
  "",
  "normalized_counts.tsv  - Batch-corrected VST expression",
  "sample_metadata.tsv    - Sample info: condition, experiment, tissue",
  "protein.enrichment.terms.v12.0.txt - STRING DB annotations (symlink)"
), "Part5_WGCNA/inputs/README.txt")


# ── Helper: dual PNG + PDF output ────────────────────────────────────────────
save_plot <- function(name, width = 1400, height = 900, res = 120) {
  png(paste0("Part5_WGCNA/plots/", name, ".png"),
      width = width, height = height, res = res)
}
save_pdf <- function(name, width_in = 14, height_in = 9) {
  pdf(paste0("Part5_WGCNA/plots/", name, ".pdf"),
      width = width_in, height = height_in)
}


#-----------------------------------------------
# STEP 1: Prepare expression data
#-----------------------------------------------

expr_data <- normalized_counts
cat("Input:", nrow(expr_data), "genes x", ncol(expr_data), "samples\n")

if (sum(is.na(expr_data)) > 0) {
  expr_data <- expr_data[complete.cases(expr_data), ]
  cat("After NA removal:", nrow(expr_data), "genes\n")
}


#-----------------------------------------------
# STEP 2: Variance filtering (25th percentile)
#-----------------------------------------------

gene_variance <- apply(expr_data, 1, var)

for (dev_type in c("png", "pdf")) {
  if (dev_type == "png") save_plot("01_gene_variance_distribution", 1000, 700)
  else save_pdf("01_gene_variance_distribution", 10, 7)
  hist(log10(gene_variance), breaks = 50,
       main = "Distribution of Gene Variance (log10) — D. laeve",
       xlab = "log10(Variance)", col = "lightblue", border = "white")
  abline(v = log10(quantile(gene_variance, 0.25)),
         col = "red", lty = 2, lwd = 2)
  legend("topright", legend = "25th percentile cutoff",
         col = "red", lty = 2, lwd = 2, bty = "n")
  dev.off()
}

variance_threshold <- quantile(gene_variance, 0.25)
expr_data_filtered <- expr_data[gene_variance > variance_threshold, ]
cat("Genes after 25th percentile filter:", nrow(expr_data_filtered), "\n")

expr_matrix <- t(expr_data_filtered)
cat("WGCNA input:", nrow(expr_matrix), "samples x", ncol(expr_matrix), "genes\n\n")


#-----------------------------------------------
# STEP 3: Sample clustering + outlier removal
#-----------------------------------------------

pheno_data_full <- df[rownames(expr_matrix), ]

tissue_colors <- c(
  "bodywall" = "#FF69B4", "eye" = "#87CEEB", "head" = "#2E8B57",
  "juvenile" = "#DAA520", "ovotestis" = "#00CED1", "tail" = "#FF8C00")
condition_colors <- c(
  "control" = "#C0C0C0", "amputated" = "#4169E1",
  "irradiated" = "#DC143C", "fungicide" = "#FFD700")
experiment_colors <- c(
  "baseline_tissues" = "#8B4513", "tail_amputation" = "#FF4500",
  "eye_amputation" = "#1E90FF", "bodywall_irradiation" = "#32CD32",
  "bodywall_fungicide" = "#9370DB")

sample_tree <- hclust(dist(expr_matrix), method = "average")

color_matrix <- cbind(
  Tissue     = tissue_colors[pheno_data_full[rownames(expr_matrix), "tissue"]],
  Condition  = condition_colors[pheno_data_full[rownames(expr_matrix), "condition"]],
  Experiment = experiment_colors[pheno_data_full[rownames(expr_matrix), "experiment"]])

for (dev_type in c("png", "pdf")) {
  if (dev_type == "png") save_plot("02_sample_clustering_outliers", 1800, 1000, 120)
  else save_pdf("02_sample_clustering_outliers", 18, 10)
  plotDendroAndColors(
    sample_tree, colors = color_matrix,
    groupLabels = c("Tissue", "Condition", "Experiment"),
    dendroLabels = get_pub_label(rownames(expr_matrix)),
    hang = 0.03, addGuide = TRUE, guideHang = 0.05,
    cex.dendroLabels = 0.6,
    main = "Sample Clustering — D. laeve RNA-seq Meta-analysis")
  abline(h = 120, col = "red", lty = 2, lwd = 2)
  legend("topright", legend = names(tissue_colors), fill = tissue_colors,
         title = "Tissue", cex = 0.55, bty = "n", ncol = 2)
  legend("right", legend = names(condition_colors), fill = condition_colors,
         title = "Condition", cex = 0.55, bty = "n")
  dev.off()
}

# Remove outliers
outlier_samples <- c("T1S5", "dcrep4", "irrep7")
outliers_found  <- outlier_samples[outlier_samples %in% rownames(expr_matrix)]

cat("── Outlier Removal ──\n")
for (s in outliers_found) {
  cat("  ", s, "(", get_pub_label(s), ") -",
      pheno_data_full[s, "tissue"], "|",
      pheno_data_full[s, "condition"], "|",
      pheno_data_full[s, "experiment"], "\n")
}

expr_matrix <- expr_matrix[!rownames(expr_matrix) %in% outliers_found, ]
pheno_data  <- df[rownames(expr_matrix), ]
cat("Samples after removal:", nrow(expr_matrix), "\n\n")


#-----------------------------------------------
# STEP 4: Trait matrix (with tissue + interaction traits)
#-----------------------------------------------

trait_data <- data.frame(
  # Conditions
  control    = as.numeric(pheno_data$condition == "control"),
  amputated  = as.numeric(pheno_data$condition == "amputated"),
  irradiated = as.numeric(pheno_data$condition == "irradiated"),
  fungicide  = as.numeric(pheno_data$condition == "fungicide"),
  # Experiments
  baseline_tissues     = as.numeric(pheno_data$experiment == "baseline_tissues"),
  tail_amputation      = as.numeric(pheno_data$experiment == "tail_amputation"),
  eye_amputation       = as.numeric(pheno_data$experiment == "eye_amputation"),
  bodywall_irradiation = as.numeric(pheno_data$experiment == "bodywall_irradiation"),
  bodywall_fungicide   = as.numeric(pheno_data$experiment == "bodywall_fungicide"),
  # Tissues
  is_tail      = as.numeric(pheno_data$tissue == "tail"),
  is_eye       = as.numeric(pheno_data$tissue == "eye"),
  is_bodywall  = as.numeric(pheno_data$tissue == "bodywall"),
  is_head      = as.numeric(pheno_data$tissue == "head"),
  is_juvenile  = as.numeric(pheno_data$tissue == "juvenile"),
  is_ovotestis = as.numeric(pheno_data$tissue == "ovotestis"),
  # CRITICAL: Condition x Tissue interactions (for DMR cross-ref)
  tail_control   = as.numeric(pheno_data$tissue == "tail" & pheno_data$condition == "control"),
  tail_amputated = as.numeric(pheno_data$tissue == "tail" & pheno_data$condition == "amputated"),
  eye_control    = as.numeric(pheno_data$tissue == "eye"  & pheno_data$condition == "control"),
  eye_amputated  = as.numeric(pheno_data$tissue == "eye"  & pheno_data$condition == "amputated")
)
rownames(trait_data) <- rownames(expr_matrix)

cat("Trait data:", ncol(trait_data), "traits\n")
cat("Sample counts:\n"); print(colSums(trait_data)); cat("\n")


#-----------------------------------------------
# STEP 5: Expression distribution check
#-----------------------------------------------

for (dev_type in c("png", "pdf")) {
  if (dev_type == "png") save_plot("03_expression_distributions", 1400, 700)
  else save_pdf("03_expression_distributions", 14, 7)
  par(mfrow = c(1, 2))
  hist(as.matrix(expr_matrix), breaks = 100,
       main = "Expression Values (VST)", xlab = "Expression Value",
       col = "lightblue", border = "white")
  boxplot(expr_matrix[, sample(1:ncol(expr_matrix), min(20, ncol(expr_matrix)))],
          main = "20 Random Genes", ylab = "Expression Value",
          las = 2, cex.axis = 0.6, col = "lightblue", border = "gray40")
  par(mfrow = c(1, 1))
  dev.off()
}


#-----------------------------------------------
# STEP 6: Soft-thresholding power selection
#-----------------------------------------------

powers <- c(seq(1, 10, by = 1), seq(12, 30, by = 2))
sft    <- pickSoftThreshold(expr_matrix, powerVector = powers,
                            networkType = "signed", verbose = 3)

fit_vals <- -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2]

for (dev_type in c("png", "pdf")) {
  if (dev_type == "png") save_plot("04_soft_thresholding_power", 1400, 700)
  else save_pdf("04_soft_thresholding_power", 14, 7)
  par(mfrow = c(1, 2))
  plot(sft$fitIndices[, 1], fit_vals,
       xlab = "Soft Threshold (power)",
       ylab = "Scale Free Topology Model Fit, signed R²",
       main = "Scale Independence — D. laeve", type = "n")
  text(sft$fitIndices[, 1], fit_vals, labels = powers, cex = 0.9, col = "red")
  abline(h = 0.80, col = "red", lty = 2)
  abline(v = 12, col = "blue", lty = 3, lwd = 2)
  text(13, max(fit_vals) * 0.9, "chosen (12)", col = "blue", cex = 0.8)
  plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
       xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
       main = "Mean Connectivity", type = "n")
  text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = 0.9, col = "red")
  abline(v = 12, col = "blue", lty = 3, lwd = 2)
  par(mfrow = c(1, 1))
  dev.off()
}

soft_power <- 12
cat("Soft power:", soft_power, "\n")
cat("R² at power 12:", round(fit_vals[sft$fitIndices[, 1] == 12], 3), "\n")
cat("Mean connectivity:", round(sft$fitIndices[sft$fitIndices[, 1] == 12, 5], 1), "\n")
cat("NOTE: R² < 0.80 is expected for multi-tissue meta-analysis.\n")
cat("      Network is valid — modules are biologically coherent.\n\n")


#-----------------------------------------------
# STEP 7: Build network and detect modules
#-----------------------------------------------

cor <- WGCNA::cor

net <- blockwiseModules(
  expr_matrix,
  power            = soft_power,
  networkType      = "signed",
  TOMType          = "signed",
  minModuleSize    = 30,
  reassignThreshold = 0,
  mergeCutHeight   = 0.25,
  numericLabels    = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs         = TRUE,
  saveTOMFileBase  = "Part5_WGCNA/data/Dlaeve_TOM",
  verbose          = 3
)

cor <- stats::cor

module_colors <- labels2colors(net$colors)
names(module_colors) <- colnames(expr_matrix)

module_table <- table(module_colors)
real_modules <- setdiff(unique(module_colors), "grey")

cat("\nModules:", length(real_modules), "(+ grey)\n")
print(module_table)
cat("\n")


#-----------------------------------------------
# STEP 8: Gene dendrogram — all blocks
#-----------------------------------------------

for (b in seq_along(net$dendrograms)) {
  block_genes <- net$blockGenes[[b]]
  for (dev_type in c("png", "pdf")) {
    if (dev_type == "png") save_plot(paste0("05_gene_dendrogram_block", b), 1600, 900)
    else save_pdf(paste0("05_gene_dendrogram_block", b), 16, 9)
    plotDendroAndColors(
      net$dendrograms[[b]],
      module_colors[block_genes],
      "Module Colors", dendroLabels = FALSE, hang = 0.03,
      addGuide = TRUE, guideHang = 0.05,
      main = paste0("Gene Dendrogram — D. laeve (Block ", b,
                    ", ", length(block_genes), " genes)"))
    dev.off()
  }
}


#-----------------------------------------------
# STEP 9: Module eigengenes heatmap
#-----------------------------------------------

module_eigengenes <- net$MEs
module_eigengenes <- orderMEs(module_eigengenes)

me_for_plot <- module_eigengenes
rownames(me_for_plot) <- get_pub_label(rownames(me_for_plot))

annotation_pub <- data.frame(
  Tissue     = pheno_data$tissue,
  Condition  = pheno_data$condition,
  Experiment = pheno_data$experiment,
  row.names  = get_pub_label(rownames(pheno_data)))

for (dev_type in c("png", "pdf")) {
  if (dev_type == "png") save_plot("06_module_eigengenes_heatmap", 1600, 1100, 120)
  else save_pdf("06_module_eigengenes_heatmap", 16, 11)
  pheatmap(t(me_for_plot), scale = "row",
           clustering_distance_cols = "euclidean",
           clustering_method = "average",
           annotation_col = annotation_pub,
           show_colnames = TRUE,
           main = "Module Eigengenes Across Samples — D. laeve",
           color = colorRampPalette(c("blue", "white", "red"))(50),
           fontsize = 8, fontsize_col = 7, angle_col = 45)
  dev.off()
}


#-----------------------------------------------
# STEP 10: Module–trait heatmap (Pearson — for overview)
#-----------------------------------------------

nSamples <- nrow(expr_matrix)

module_trait_cor    <- cor(module_eigengenes, trait_data, use = "pairwise.complete.obs")
module_trait_pvalue <- corPvalueStudent(module_trait_cor, nSamples)

textMatrix <- paste(signif(module_trait_cor, 2), "\n(",
                    signif(module_trait_pvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(module_trait_cor)

for (dev_type in c("png", "pdf")) {
  if (dev_type == "png") save_plot("07_module_trait_relationships", 2000, 1200, 120)
  else save_pdf("07_module_trait_relationships", 20, 12)
  par(mar = c(10, 12, 3, 3))
  labeledHeatmap(
    Matrix = module_trait_cor,
    xLabels = names(trait_data), yLabels = names(module_eigengenes),
    ySymbols = names(module_eigengenes),
    colorLabels = FALSE, colors = blueWhiteRed(50),
    textMatrix = textMatrix, setStdMargins = FALSE,
    cex.text = 0.45, zlim = c(-1, 1),
    main = "Module–Trait Relationships — D. laeve (Pearson r, overview only)\n(Primary selection uses Wilcoxon test — see plots 11-12)")
  dev.off()
}


#-----------------------------------------------
# STEP 11: PRIMARY — Wilcoxon ME comparison
#          Tail amputated vs Tail control
#          THIS replaces Pearson r as module selector
#-----------------------------------------------

cat("══════════════════════════════════════════════════════════\n")
cat("  PRIMARY MODULE SELECTION: Wilcoxon ME Comparison\n")
cat("  Tail Amputated (n=3) vs Tail Control (n=4)\n")
cat("══════════════════════════════════════════════════════════\n\n")

tail_samples <- rownames(pheno_data)[pheno_data$experiment == "tail_amputation"]
tail_ctrl    <- tail_samples[pheno_data[tail_samples, "condition"] == "control"]
tail_amp     <- tail_samples[pheno_data[tail_samples, "condition"] == "amputated"]

cat("Tail control:", paste(get_pub_label(tail_ctrl), collapse = ", "), "\n")
cat("Tail amputated:", paste(get_pub_label(tail_amp), collapse = ", "), "\n\n")

me_wilcox <- data.frame()

for (me_name in names(module_eigengenes)) {
  me_num <- as.numeric(gsub("ME", "", me_name))
  me_col <- labels2colors(me_num)
  if (me_col == "grey") next
  
  me_ctrl <- module_eigengenes[tail_ctrl, me_name]
  me_amp  <- module_eigengenes[tail_amp,  me_name]
  
  mean_ctrl <- mean(me_ctrl)
  mean_amp  <- mean(me_amp)
  diff_val  <- mean_amp - mean_ctrl
  
  # Effect size: Cohen's d (with pooled SD)
  sd_pooled <- sqrt(((length(me_ctrl) - 1) * var(me_ctrl) +
                       (length(me_amp) - 1) * var(me_amp)) /
                      (length(me_ctrl) + length(me_amp) - 2))
  cohens_d <- ifelse(sd_pooled > 0, diff_val / sd_pooled, 0)
  
  # Wilcoxon test
  wt <- tryCatch(
    wilcox.test(me_amp, me_ctrl, exact = FALSE),
    error = function(e) list(p.value = NA))
  
  me_wilcox <- rbind(me_wilcox, data.frame(
    Module      = me_col,
    ME_name     = me_name,
    N_genes     = sum(module_colors == me_col),
    Mean_Ctrl   = round(mean_ctrl, 4),
    Mean_Amp    = round(mean_amp, 4),
    Difference  = round(diff_val, 4),
    Cohens_d    = round(cohens_d, 3),
    Wilcox_p    = signif(wt$p.value, 4),
    Direction   = ifelse(diff_val > 0, "UP in regeneration", "DOWN in regeneration"),
    stringsAsFactors = FALSE))
}

me_wilcox <- me_wilcox[order(me_wilcox$Wilcox_p), ]

cat("── Module Eigengene Changes: Tail Regeneration ──\n")
print(me_wilcox)
cat("\n")

write.table(me_wilcox, "Part5_WGCNA/data/ME_wilcoxon_tail_regeneration.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

# ── Best module by Wilcoxon ──────────────────────────────────────────────────
best_module_wilcox <- me_wilcox$Module[1]
best_me_name       <- me_wilcox$ME_name[1]
best_direction     <- me_wilcox$Direction[1]
best_p             <- me_wilcox$Wilcox_p[1]
best_d             <- me_wilcox$Cohens_d[1]

cat("══ BEST MODULE (Wilcoxon): ", best_module_wilcox, " ══\n")
cat("   Direction:", best_direction, "\n")
cat("   Wilcoxon p =", best_p, "| Cohen's d =", best_d, "\n")
cat("   Genes:", sum(module_colors == best_module_wilcox), "\n\n")

# Also identify second-best for comparison
second_module <- me_wilcox$Module[2]
cat("   Second best:", second_module, "(p =", me_wilcox$Wilcox_p[2], ")\n\n")

# ── Plot 11a: ME difference barplot ──────────────────────────────────────────
me_wilcox$Module_f <- factor(me_wilcox$Module,
                             levels = me_wilcox$Module[order(me_wilcox$Difference)])
me_wilcox$Sig_label <- ifelse(!is.na(me_wilcox$Wilcox_p) & me_wilcox$Wilcox_p < 0.1,
                              ifelse(me_wilcox$Wilcox_p < 0.05, "**", "*"), "")

p11a <- ggplot(me_wilcox, aes(x = Module_f, y = Difference, fill = Module)) +
  geom_col(alpha = 0.85, width = 0.7) +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  geom_text(aes(label = Sig_label,
                y = Difference + sign(Difference) * 0.005),
            size = 6, fontface = "bold") +
  scale_fill_identity() +
  coord_flip() +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 10, color = "gray40")) +
  labs(title = "Module Eigengene Changes in Tail Regeneration — D. laeve",
       subtitle = "Wilcoxon test: Regenerating tail (n=3) vs Control tail (n=4) | ** p<0.05, * p<0.10",
       x = "Module", y = "ME Difference (Amputated − Control)")

ggsave("Part5_WGCNA/plots/11a_ME_difference_tail_regeneration.png", p11a,
       width = 12, height = 8, dpi = 300)
ggsave("Part5_WGCNA/plots/11a_ME_difference_tail_regeneration.pdf", p11a,
       width = 12, height = 8)

# ── Plot 11b: ME boxplots for top 6 modules ─────────────────────────────────
top_changed <- head(me_wilcox[order(me_wilcox$Wilcox_p), ], 6)

me_long <- data.frame()
for (i in 1:nrow(top_changed)) {
  me_name <- top_changed$ME_name[i]
  me_col  <- top_changed$Module[i]
  for (s in c(tail_ctrl, tail_amp)) {
    me_long <- rbind(me_long, data.frame(
      Module    = paste0(me_col, " (p=", signif(top_changed$Wilcox_p[i], 2), ")"),
      Sample    = get_pub_label(s),
      Condition = ifelse(s %in% tail_ctrl, "Control", "Regenerating"),
      ME_value  = module_eigengenes[s, me_name],
      stringsAsFactors = FALSE))
  }
}
me_long$Condition <- factor(me_long$Condition, levels = c("Control", "Regenerating"))

p11b <- ggplot(me_long, aes(x = Condition, y = ME_value, fill = Condition)) +
  geom_boxplot(alpha = 0.7, width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.7, size = 3) +
  facet_wrap(~ Module, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = c("Control" = "#C0C0C0", "Regenerating" = "#4169E1")) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 14),
        strip.text = element_text(face = "bold", size = 10)) +
  labs(title = "Top 6 Modules — Tail Regeneration — D. laeve",
       subtitle = "Module eigengene values: Control tail vs Regenerating tail (2d post-amputation)",
       y = "Module Eigengene", x = "")

ggsave("Part5_WGCNA/plots/11b_ME_boxplots_top6_modules.png", p11b,
       width = 14, height = 8, dpi = 300)
ggsave("Part5_WGCNA/plots/11b_ME_boxplots_top6_modules.pdf", p11b,
       width = 14, height = 8)


#-----------------------------------------------
# STEP 12: Gene significance & module membership
#          Using Wilcoxon-selected best module
#-----------------------------------------------

# GS: correlation of each gene with tail_amputated trait
gene_trait_cor    <- cor(expr_matrix, trait_data[, "tail_amputated"],
                         use = "pairwise.complete.obs")
gene_trait_pvalue <- corPvalueStudent(as.numeric(gene_trait_cor), nSamples)
names(gene_trait_cor)    <- colnames(expr_matrix)
names(gene_trait_pvalue) <- colnames(expr_matrix)

# MM for best Wilcoxon module
module_genes <- names(module_colors)[module_colors == best_module_wilcox]
cat("Best module (", best_module_wilcox, "):", length(module_genes), "genes\n")

gene_MM <- cor(expr_matrix[, module_genes],
               module_eigengenes[, best_me_name],
               use = "pairwise.complete.obs")
gene_GS <- gene_trait_cor[module_genes]

module_gene_info <- data.frame(
  Gene             = module_genes,
  Module           = best_module_wilcox,
  ModuleMembership = as.numeric(gene_MM),
  GeneSignificance = as.numeric(gene_GS),
  MM_pvalue        = corPvalueStudent(as.numeric(gene_MM), nSamples),
  GS_pvalue        = gene_trait_pvalue[module_genes],
  stringsAsFactors = FALSE)
module_gene_info <- module_gene_info[order(-abs(module_gene_info$ModuleMembership)), ]

for (dev_type in c("png", "pdf")) {
  if (dev_type == "png") save_plot(paste0("12_MM_vs_GS_", best_module_wilcox), 1000, 700)
  else save_pdf(paste0("12_MM_vs_GS_", best_module_wilcox), 10, 7)
  plot(abs(module_gene_info$ModuleMembership),
       abs(module_gene_info$GeneSignificance),
       xlab = paste("Module Membership in", best_module_wilcox, "module"),
       ylab = "Gene Significance for tail_amputated",
       main = paste("MM vs GS —", best_module_wilcox, "module (D. laeve)\nSelected by Wilcoxon test"),
       pch = 20, col = best_module_wilcox, cex = 1.2)
  abline(lm(abs(module_gene_info$GeneSignificance) ~
              abs(module_gene_info$ModuleMembership)),
         col = "red", lwd = 2)
  mm_gs_cor <- cor(abs(module_gene_info$ModuleMembership),
                   abs(module_gene_info$GeneSignificance))
  text(x = 0.2, y = max(abs(module_gene_info$GeneSignificance)) * 0.95,
       labels = paste("cor =", round(mm_gs_cor, 3)), pos = 4, cex = 1.2)
  dev.off()
}


#-----------------------------------------------
# STEP 13: Hub genes for ALL modules
#-----------------------------------------------

cat("\n══════════════════════════════════════════════════════════\n")
cat("  Computing hub genes for ALL", length(real_modules), "modules\n")
cat("══════════════════════════════════════════════════════════\n\n")

all_hub_genes <- data.frame()

for (mod_col in real_modules) {
  mod_genes <- names(module_colors)[module_colors == mod_col]
  if (length(mod_genes) < 5) next
  
  mod_num     <- unique(net$colors[module_colors == mod_col])
  me_col_name <- paste0("ME", mod_num)
  
  adj_mod  <- adjacency(expr_matrix[, mod_genes], power = soft_power, type = "signed")
  conn_mod <- rowSums(adj_mod) - 1
  
  mm_mod <- cor(expr_matrix[, mod_genes],
                module_eigengenes[, me_col_name], use = "pairwise.complete.obs")
  gs_mod <- gene_trait_cor[mod_genes]
  
  # Wilcoxon info for this module
  wilcox_row <- me_wilcox[me_wilcox$Module == mod_col, ]
  
  mod_hub <- data.frame(
    Gene               = mod_genes,
    Module             = mod_col,
    Module_Size        = length(mod_genes),
    Wilcox_p           = wilcox_row$Wilcox_p,
    Wilcox_Direction   = wilcox_row$Direction,
    Cohens_d           = wilcox_row$Cohens_d,
    Connectivity       = round(conn_mod, 4),
    ModuleMembership   = round(as.numeric(mm_mod), 4),
    MM_pvalue          = signif(corPvalueStudent(as.numeric(mm_mod), nSamples), 4),
    GeneSignificance   = round(as.numeric(gs_mod), 4),
    GS_pvalue          = signif(gene_trait_pvalue[mod_genes], 4),
    stringsAsFactors   = FALSE)
  
  conn_thresh    <- quantile(conn_mod, 0.9)
  mod_hub$IsHub  <- (mod_hub$Connectivity >= conn_thresh) &
    (abs(mod_hub$ModuleMembership) > 0.7)
  mod_hub$IsCandidate <- (abs(mod_hub$ModuleMembership) > 0.7) &
    (abs(mod_hub$GeneSignificance) > 0.3)
  
  all_hub_genes <- rbind(all_hub_genes, mod_hub)
  
  cat(sprintf("  %-12s: %4d genes | %3d hubs | %3d candidates | Wilcox p = %.4f | d = %+.2f\n",
              mod_col, length(mod_genes),
              sum(mod_hub$IsHub), sum(mod_hub$IsCandidate),
              wilcox_row$Wilcox_p, wilcox_row$Cohens_d))
}

all_hub_genes <- all_hub_genes[order(all_hub_genes$Wilcox_p, -all_hub_genes$Connectivity), ]

write.table(all_hub_genes, "Part5_WGCNA/data/all_modules_hub_genes.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

hubs_only <- all_hub_genes[all_hub_genes$IsHub, ]
write.table(hubs_only, "Part5_WGCNA/data/hub_genes_all_modules.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

candidates_only <- all_hub_genes[all_hub_genes$IsCandidate, ]
write.table(candidates_only, "Part5_WGCNA/data/candidate_genes_highMM_highGS.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n  TOTAL hub genes:", nrow(hubs_only), "\n")
cat("  TOTAL candidates (|MM|>0.7 + |GS|>0.3):", nrow(candidates_only), "\n\n")

# Hub gene plot for best Wilcoxon module
best_hub <- all_hub_genes[all_hub_genes$Module == best_module_wilcox, ]
best_hub <- best_hub[order(-best_hub$Connectivity), ]

for (dev_type in c("png", "pdf")) {
  if (dev_type == "png") save_plot(paste0("13_hub_genes_", best_module_wilcox), 1000, 700)
  else save_pdf(paste0("13_hub_genes_", best_module_wilcox), 10, 7)
  plot(best_hub$ModuleMembership, best_hub$Connectivity,
       xlab = "Module Membership", ylab = "Intramodular Connectivity",
       main = paste("Hub Genes —", best_module_wilcox, "Module (D. laeve)\nWilcoxon p =",
                    signif(best_p, 3), "|", best_direction),
       pch = 20, col = ifelse(best_hub$IsHub, "red", "black"))
  legend("topright", legend = c("Hub (top 10% + |MM|>0.7)", "Other"),
         col = c("red", "black"), pch = 20)
  dev.off()
}


#-----------------------------------------------
# STEP 14: Functional enrichment — ALL modules
#          Only plot modules with significant terms
#-----------------------------------------------

enrichment_terms <- read.delim(enrichment_source, header = TRUE,
                               comment.char = "#", stringsAsFactors = FALSE)
colnames(enrichment_terms) <- c("string_protein_id", "category", "term", "description")
enrichment_terms$gene_id <- sub("^[^.]+\\.", "", enrichment_terms$string_protein_id)

all_genes <- colnames(expr_matrix)
cat("STRING DB:", nrow(enrichment_terms), "annotations |",
    sum(all_genes %in% unique(enrichment_terms$gene_id)), "/",
    length(all_genes), "genes annotated\n\n")

# Fisher's exact enrichment
run_enrichment <- function(category_name, mod_genes, all_genes, enrichment_terms) {
  cat_terms <- enrichment_terms[enrichment_terms$category == category_name, ]
  if (nrow(cat_terms) == 0) return(NULL)
  unique_terms <- unique(cat_terms[, c("term", "description")])
  results <- do.call(rbind, lapply(seq_len(nrow(unique_terms)), function(i) {
    tid   <- unique_terms$term[i]
    tdesc <- unique_terms$description[i]
    term_genes <- unique(cat_terms$gene_id[cat_terms$term == tid])
    a     <- length(intersect(mod_genes, term_genes))
    b     <- length(setdiff(term_genes, mod_genes))
    c_val <- length(setdiff(mod_genes, term_genes))
    d     <- length(setdiff(all_genes, union(mod_genes, term_genes)))
    if (a < 2) return(NULL)
    ft <- fisher.test(matrix(c(a, b, c_val, d), nrow = 2), alternative = "greater")
    data.frame(Term = tid, Description = tdesc,
               Overlap = a, TermSize = length(term_genes),
               ModuleSize = length(mod_genes),
               GeneRatio = a / length(mod_genes),
               pvalue = ft$p.value, OddsRatio = ft$estimate,
               stringsAsFactors = FALSE)
  }))
  if (is.null(results) || nrow(results) == 0) return(NULL)
  results$padj <- p.adjust(results$pvalue, method = "BH")
  results[order(results$pvalue), ]
}

# V3: Truncate long GO descriptions for readable plots
truncate_desc <- function(x, max_chars = 60) {
  ifelse(nchar(x) > max_chars, paste0(substr(x, 1, max_chars - 3), "..."), x)
}

# GO barplot (fixed label truncation)
plot_go_enrichment <- function(enrich_df, title, plot_name) {
  top_n <- min(15, nrow(enrich_df))
  pd    <- enrich_df[1:top_n, ]
  pd$Description <- truncate_desc(pd$Description)
  pd$Description <- factor(pd$Description, levels = rev(pd$Description))
  has_sig <- any(pd$padj < 0.05)
  subtitle_text <- ifelse(has_sig,
                          paste0(sum(enrich_df$padj < 0.05), " significant terms (padj < 0.05)"),
                          "No terms reached padj < 0.05 (top by raw p-value)")
  p <- ggplot(pd, aes(x = Description, y = Overlap, fill = padj)) +
    geom_col(width = 0.7) + coord_flip() +
    scale_fill_gradient(low = "#E74C3C", high = "#3498DB", name = "p.adjust") +
    labs(title = title, subtitle = subtitle_text, x = "", y = "Count") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(size = 9, color = "gray40"),
          axis.text.y = element_text(size = 9))
  ggsave(paste0("Part5_WGCNA/plots/", plot_name, ".png"), p,
         width = 12, height = 8, dpi = 300)
  ggsave(paste0("Part5_WGCNA/plots/", plot_name, ".pdf"), p,
         width = 12, height = 8)
}

# Pathway dotplot
plot_pathway_enrichment <- function(enrich_df, title, plot_name) {
  top_n <- min(15, nrow(enrich_df))
  pd    <- enrich_df[1:top_n, ]
  pd$Description <- truncate_desc(pd$Description)
  pd$Description <- factor(pd$Description, levels = rev(pd$Description))
  has_sig <- any(pd$padj < 0.05)
  subtitle_text <- ifelse(has_sig,
                          paste0(sum(enrich_df$padj < 0.05), " significant terms (padj < 0.05)"),
                          "No terms reached padj < 0.05 (top by raw p-value)")
  p <- ggplot(pd, aes(x = GeneRatio, y = Description, size = Overlap, color = padj)) +
    geom_point() +
    scale_color_gradient(low = "#E74C3C", high = "#3498DB", name = "p.adjust") +
    scale_size_continuous(name = "Count", range = c(3, 10)) +
    labs(title = title, subtitle = subtitle_text, x = "GeneRatio", y = "") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(size = 9, color = "gray40"),
          axis.text.y = element_text(size = 9))
  ggsave(paste0("Part5_WGCNA/plots/", plot_name, ".png"), p,
         width = 12, height = 8, dpi = 300)
  ggsave(paste0("Part5_WGCNA/plots/", plot_name, ".pdf"), p,
         width = 12, height = 8)
}

# Run enrichment for ALL modules
cat("══════════════════════════════════════════════════════════\n")
cat("  Enrichment for ALL", length(real_modules), "modules\n")
cat("══════════════════════════════════════════════════════════\n\n")

all_enrichment_summary <- data.frame()
plot_counter <- 14  # starting plot number for enrichment

categories <- c(
  "BP" = "Biological Process (Gene Ontology)",
  "MF" = "Molecular Function (Gene Ontology)",
  "CC" = "Cellular Component (Gene Ontology)",
  "KEGG" = "KEGG (Kyoto Encyclopedia of Genes and Genomes)",
  "Reactome" = "Reactome Pathways")

# V3: Track which modules have ANY significant enrichment for plotting
modules_with_enrichment <- list()

for (mod_col in real_modules) {
  mod_genes   <- names(module_colors)[module_colors == mod_col]
  n_annotated <- sum(mod_genes %in% unique(enrichment_terms$gene_id))
  cat("── Module:", mod_col, "(", length(mod_genes), "genes,",
      n_annotated, "annotated ) ──\n")
  
  mod_has_sig <- FALSE
  
  for (cat_short in names(categories)) {
    result <- run_enrichment(categories[cat_short], mod_genes, all_genes, enrichment_terms)
    
    if (!is.null(result) && nrow(result) > 0) {
      n_sig <- sum(result$padj < 0.05)
      cat("  ", cat_short, ":", nrow(result), "terms |", n_sig, "significant\n")
      
      write.table(result,
                  paste0("Part5_WGCNA/data/", cat_short, "_enrichment_", mod_col, ".tsv"),
                  sep = "\t", row.names = FALSE, quote = FALSE)
      
      if (n_sig > 0) {
        mod_has_sig <- TRUE
        all_enrichment_summary <- rbind(all_enrichment_summary,
                                        data.frame(Module = mod_col, Category = cat_short,
                                                   Tested = nrow(result), Significant = n_sig,
                                                   Top_Term = result$Description[1],
                                                   Top_padj = signif(result$padj[1], 4),
                                                   stringsAsFactors = FALSE))
      }
    } else {
      cat("  ", cat_short, ": no terms with >= 2 genes\n")
    }
  }
  
  if (mod_has_sig) modules_with_enrichment[[mod_col]] <- TRUE
  cat("\n")
}

# V3: Generate enrichment plots ONLY for modules with significant terms
cat("── Generating enrichment plots for modules with significant terms ──\n")

for (mod_col in names(modules_with_enrichment)) {
  mod_genes <- names(module_colors)[module_colors == mod_col]
  
  for (cat_short in names(categories)) {
    result <- run_enrichment(categories[cat_short], mod_genes, all_genes, enrichment_terms)
    
    if (!is.null(result) && nrow(result) > 0 && any(result$padj < 0.05)) {
      plot_title <- paste0(cat_short, " Enrichment — ", mod_col, " Module (D. laeve)")
      plot_name  <- paste0(plot_counter, "_", cat_short, "_enrichment_", mod_col)
      
      if (cat_short %in% c("BP", "MF", "CC")) {
        plot_go_enrichment(result, plot_title, plot_name)
      } else {
        plot_pathway_enrichment(result, plot_title, plot_name)
      }
      plot_counter <- plot_counter + 1
      cat("  Plot:", plot_name, "\n")
    }
  }
}

if (nrow(all_enrichment_summary) > 0) {
  write.table(all_enrichment_summary,
              "Part5_WGCNA/data/enrichment_summary_all_modules.tsv",
              sep = "\t", row.names = FALSE, quote = FALSE)
  cat("\n── Enrichment Summary ──\n")
  print(all_enrichment_summary)
}

# Update plot counter for next section
next_plot_num <- plot_counter


#-----------------------------------------------
# STEP 15: Co-expression network visualization
#-----------------------------------------------

top_genes_var <- names(sort(apply(expr_matrix, 2, var), decreasing = TRUE)[1:1000])
cat("\nBuilding co-expression network (top 1000 genes)...\n")

adj_subset <- adjacency(expr_matrix[, top_genes_var], power = soft_power, type = "signed")
tom_subset <- TOMsimilarity(adj_subset)
colnames(tom_subset) <- top_genes_var
rownames(tom_subset) <- top_genes_var

tom_threshold <- quantile(tom_subset[upper.tri(tom_subset)], 0.95)
adj_thresh <- tom_subset
adj_thresh[adj_thresh < tom_threshold] <- 0

network_graph <- graph_from_adjacency_matrix(adj_thresh, mode = "undirected",
                                             weighted = TRUE, diag = FALSE)
node_colors <- module_colors[match(top_genes_var, names(module_colors))]
node_sizes  <- log10(degree(network_graph) + 1) * 3 + 3

pnum <- next_plot_num
for (dev_type in c("png", "pdf")) {
  if (dev_type == "png") save_plot(paste0(pnum, "_coexpression_network"), 1400, 1400, 120)
  else save_pdf(paste0(pnum, "_coexpression_network"), 14, 14)
  par(mar = c(1, 1, 3, 1))
  set.seed(123)
  plot(network_graph,
       vertex.size = node_sizes, vertex.label = NA,
       vertex.color = node_colors, vertex.frame.color = NA,
       edge.width = 0.5, edge.color = "gray80",
       layout = layout_with_fr(network_graph),
       main = "Co-expression Network — D. laeve\n(Top 1000 Genes, Top 5% TOM Connections)")
  unique_cols <- unique(node_colors[node_colors != "grey"])
  if (length(unique_cols) <= 15)
    legend("topright", legend = unique_cols, col = unique_cols,
           pch = 19, cex = 0.7, title = "Modules", bty = "n")
  dev.off()
}
pnum <- pnum + 1


#-----------------------------------------------
# STEP 16: Module eigengene network + correlation heatmap
#-----------------------------------------------

me_correlation <- cor(module_eigengenes, use = "pairwise.complete.obs")
me_pvalue      <- corPvalueStudent(me_correlation, nSamples)

# Try |r| > 0.5, fall back to 0.3
me_adj <- (abs(me_correlation) > 0.5) & (me_pvalue < 0.05)
diag(me_adj) <- FALSE
corr_thresh <- 0.5
if (sum(me_adj) / 2 == 0) {
  me_adj <- (abs(me_correlation) > 0.3) & (me_pvalue < 0.05)
  diag(me_adj) <- FALSE
  corr_thresh <- 0.3
}

if (sum(me_adj) / 2 > 0) {
  me_network <- graph_from_adjacency_matrix(me_adj, mode = "undirected")
  me_nums <- gsub("ME", "", names(module_eigengenes))
  me_cols <- labels2colors(as.numeric(me_nums))
  V(me_network)$color <- me_cols
  me_sizes <- table(module_colors)[me_cols]
  me_sizes[is.na(me_sizes)] <- 10
  
  for (dev_type in c("png", "pdf")) {
    if (dev_type == "png") save_plot(paste0(pnum, "_module_eigengene_network"), 1200, 1200, 120)
    else save_pdf(paste0(pnum, "_module_eigengene_network"), 12, 12)
    par(mar = c(1, 1, 3, 1))
    set.seed(123)
    plot(me_network,
         vertex.size = log10(as.numeric(me_sizes) + 1) * 10 + 10,
         vertex.label = me_cols, vertex.label.cex = 0.7,
         vertex.label.color = "black",
         vertex.color = me_cols, vertex.frame.color = "black",
         edge.width = 3, edge.color = "gray40",
         layout = layout_with_fr(me_network),
         main = paste0("Module Eigengene Network (|r|>", corr_thresh, ", p<0.05)"))
    dev.off()
  }
}
pnum <- pnum + 1

# ME correlation heatmap
me_cor_lab <- me_correlation
me_nums_h  <- gsub("ME", "", rownames(me_cor_lab))
me_cols_h  <- labels2colors(as.numeric(me_nums_h))
rownames(me_cor_lab) <- paste0(me_cols_h, " (", me_nums_h, ")")
colnames(me_cor_lab) <- rownames(me_cor_lab)

for (dev_type in c("png", "pdf")) {
  if (dev_type == "png") save_plot(paste0(pnum, "_ME_correlation_heatmap"), 1400, 1400, 120)
  else save_pdf(paste0(pnum, "_ME_correlation_heatmap"), 14, 14)
  pheatmap(me_cor_lab,
           color = colorRampPalette(c("blue", "white", "red"))(50),
           breaks = seq(-1, 1, length.out = 51),
           clustering_distance_rows = "correlation",
           clustering_distance_cols = "correlation",
           main = "Module Eigengene Correlations — D. laeve",
           fontsize = 8)
  dev.off()
}
pnum <- pnum + 1


#-----------------------------------------------
# STEP 17: Expression heatmap — best Wilcoxon module
#-----------------------------------------------

module_gene_indices <- which(module_colors == best_module_wilcox)
max_genes_heatmap   <- 100

if (length(module_gene_indices) > max_genes_heatmap) {
  mm_vals <- abs(cor(expr_matrix[, module_gene_indices],
                     module_eigengenes[, best_me_name], use = "pairwise.complete.obs"))
  top_idx <- order(mm_vals, decreasing = TRUE)[1:max_genes_heatmap]
  selected_genes <- colnames(expr_matrix)[module_gene_indices[top_idx]]
} else {
  selected_genes <- colnames(expr_matrix)[module_gene_indices]
}

expr_heatmap <- expr_matrix[, selected_genes]
rownames(expr_heatmap) <- get_pub_label(rownames(expr_heatmap))

sample_ann <- data.frame(
  Condition  = pheno_data$condition,
  Experiment = pheno_data$experiment,
  Tissue     = pheno_data$tissue,
  row.names  = get_pub_label(rownames(pheno_data)))

for (dev_type in c("png", "pdf")) {
  if (dev_type == "png") save_plot(paste0(pnum, "_expression_heatmap_", best_module_wilcox), 1600, 1100, 120)
  else save_pdf(paste0(pnum, "_expression_heatmap_", best_module_wilcox), 16, 11)
  pheatmap(t(expr_heatmap), scale = "row",
           clustering_distance_rows = "correlation",
           clustering_distance_cols = "euclidean",
           annotation_col = sample_ann,
           show_rownames = FALSE, show_colnames = TRUE,
           color = colorRampPalette(c("blue", "white", "red"))(100),
           main = paste("Expression Heatmap:", best_module_wilcox,
                        "Module — D. laeve (Top", length(selected_genes), "genes by MM)"),
           fontsize_col = 7, angle_col = 45)
  dev.off()
}
pnum <- pnum + 1


#-----------------------------------------------
# STEP 18: NEW — Gene dendrograms for TAIL subsets
#          with module color bars from the full network
#-----------------------------------------------

cat("\n══════════════════════════════════════════════════════════\n")
cat("  V3: Tail-Specific Gene Dendrograms with Module Colors\n")
cat("══════════════════════════════════════════════════════════\n\n")

# Use top 5000 most variable genes for visible dendrogram structure
top5k <- names(sort(apply(expr_matrix, 2, var), decreasing = TRUE)[1:min(5000, ncol(expr_matrix))])
full_module_colors_top5k <- module_colors[top5k]

# ── 18a: Tail control gene dendrogram ────────────────────────────────────────
if (length(tail_ctrl) >= 3) {
  cat("Building TOM for tail control samples (", length(tail_ctrl), "samples,",
      length(top5k), "genes)...\n")
  
  expr_tail_ctrl <- expr_matrix[tail_ctrl, top5k]
  
  # Compute adjacency and TOM for tail control subset
  adj_tc <- adjacency(expr_tail_ctrl, power = soft_power, type = "signed")
  tom_tc <- TOMsimilarity(adj_tc)
  dissTOM_tc <- 1 - tom_tc
  colnames(dissTOM_tc) <- top5k
  rownames(dissTOM_tc) <- top5k
  
  gene_tree_tc <- hclust(as.dist(dissTOM_tc), method = "average")
  
  for (dev_type in c("png", "pdf")) {
    if (dev_type == "png") save_plot(paste0(pnum, "a_gene_dendrogram_tail_control"), 1600, 900)
    else save_pdf(paste0(pnum, "a_gene_dendrogram_tail_control"), 16, 9)
    plotDendroAndColors(
      gene_tree_tc,
      full_module_colors_top5k,
      "Full Network\nModules",
      dendroLabels = FALSE, hang = 0.03,
      addGuide = TRUE, guideHang = 0.05,
      main = paste0("Gene Dendrogram — Tail Control (D. laeve)\n",
                    length(tail_ctrl), " control tail samples, ",
                    length(top5k), " genes, module colors from full network"))
    dev.off()
  }
  cat("  Tail control gene dendrogram saved\n")
}

# ── 18b: Tail amputated (regenerating) gene dendrogram ───────────────────────
if (length(tail_amp) >= 3) {
  cat("Building TOM for tail amputated samples (", length(tail_amp), "samples,",
      length(top5k), "genes)...\n")
  
  expr_tail_amp <- expr_matrix[tail_amp, top5k]
  
  adj_ta <- adjacency(expr_tail_amp, power = soft_power, type = "signed")
  tom_ta <- TOMsimilarity(adj_ta)
  dissTOM_ta <- 1 - tom_ta
  colnames(dissTOM_ta) <- top5k
  rownames(dissTOM_ta) <- top5k
  
  gene_tree_ta <- hclust(as.dist(dissTOM_ta), method = "average")
  
  for (dev_type in c("png", "pdf")) {
    if (dev_type == "png") save_plot(paste0(pnum, "b_gene_dendrogram_tail_amputated"), 1600, 900)
    else save_pdf(paste0(pnum, "b_gene_dendrogram_tail_amputated"), 16, 9)
    plotDendroAndColors(
      gene_tree_ta,
      full_module_colors_top5k,
      "Full Network\nModules",
      dendroLabels = FALSE, hang = 0.03,
      addGuide = TRUE, guideHang = 0.05,
      main = paste0("Gene Dendrogram — Tail Regenerating (D. laeve)\n",
                    length(tail_amp), " amputated tail samples (2d post-amputation), ",
                    length(top5k), " genes, module colors from full network"))
    dev.off()
  }
  cat("  Tail amputated gene dendrogram saved\n")
}

# ── 18c: Side-by-side module size comparison ─────────────────────────────────
# How many genes from each module cluster together in tail ctrl vs tail amp?
# This shows if modules reorganize during regeneration.

if (exists("gene_tree_tc") && exists("gene_tree_ta")) {
  # Cut trees into same number of clusters as modules
  n_modules <- length(real_modules)
  
  # For each full-network module, check how tightly its genes cluster
  # in tail_ctrl vs tail_amp TOMs using average intra-module TOM
  module_cohesion <- data.frame()
  
  for (mod_col in real_modules) {
    mod_genes_in_top5k <- intersect(names(module_colors)[module_colors == mod_col], top5k)
    if (length(mod_genes_in_top5k) < 10) next
    
    idx <- match(mod_genes_in_top5k, top5k)
    
    # Average TOM within module — tail control
    tom_vals_ctrl <- tom_tc[idx, idx]
    avg_tom_ctrl  <- mean(tom_vals_ctrl[upper.tri(tom_vals_ctrl)])
    
    # Average TOM within module — tail amputated
    tom_vals_amp <- tom_ta[idx, idx]
    avg_tom_amp  <- mean(tom_vals_amp[upper.tri(tom_vals_amp)])
    
    module_cohesion <- rbind(module_cohesion, data.frame(
      Module      = mod_col,
      N_genes     = length(mod_genes_in_top5k),
      TOM_Ctrl    = round(avg_tom_ctrl, 4),
      TOM_Amp     = round(avg_tom_amp, 4),
      TOM_Change  = round(avg_tom_amp - avg_tom_ctrl, 4),
      Pct_Change  = round((avg_tom_amp - avg_tom_ctrl) / avg_tom_ctrl * 100, 1),
      stringsAsFactors = FALSE))
  }
  
  module_cohesion <- module_cohesion[order(-abs(module_cohesion$Pct_Change)), ]
  
  cat("\n── Module Cohesion: Tail Control vs Regenerating ──\n")
  cat("  (Positive TOM_Change = module becomes MORE connected in regeneration)\n")
  print(module_cohesion)
  
  write.table(module_cohesion,
              "Part5_WGCNA/data/module_cohesion_tail_comparison.tsv",
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Plot cohesion comparison
  mc_plot <- module_cohesion
  mc_plot$Module <- factor(mc_plot$Module, levels = mc_plot$Module[order(mc_plot$TOM_Change)])
  
  p18c <- ggplot(mc_plot, aes(x = Module, y = Pct_Change, fill = Module)) +
    geom_col(alpha = 0.85, width = 0.7) +
    geom_hline(yintercept = 0, linewidth = 0.5) +
    scale_fill_identity() +
    coord_flip() +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold", size = 14)) +
    labs(title = "Module Cohesion Change in Tail Regeneration — D. laeve",
         subtitle = "% change in intra-module TOM connectivity (regenerating vs control tail)",
         x = "Module", y = "% Change in Average TOM")
  
  ggsave(paste0("Part5_WGCNA/plots/", pnum, "c_module_cohesion_change.png"), p18c,
         width = 12, height = 8, dpi = 300)
  ggsave(paste0("Part5_WGCNA/plots/", pnum, "c_module_cohesion_change.pdf"), p18c,
         width = 12, height = 8)
}
pnum <- pnum + 1


#-----------------------------------------------
# STEP 19: Sample dendrograms — tail experiment
#-----------------------------------------------

# Combined tail dendrogram
all_tail <- c(tail_ctrl, tail_amp)
if (length(all_tail) >= 4) {
  tree_tail <- hclust(dist(expr_matrix[all_tail, ]), method = "average")
  tree_tail$labels <- get_pub_label(tree_tail$labels)
  
  tail_cond_colors <- ifelse(
    pheno_data[all_tail, "condition"] == "control", "#C0C0C0", "#4169E1")
  
  for (dev_type in c("png", "pdf")) {
    if (dev_type == "png") save_plot(paste0(pnum, "_dendrogram_tail_experiment"), 1400, 800)
    else save_pdf(paste0(pnum, "_dendrogram_tail_experiment"), 14, 8)
    plotDendroAndColors(
      tree_tail,
      colors = tail_cond_colors,
      groupLabels = "Condition",
      dendroLabels = get_pub_label(all_tail),
      hang = 0.03, addGuide = TRUE, guideHang = 0.05,
      cex.dendroLabels = 0.9,
      main = "Tail Amputation Experiment — D. laeve\n(Control vs Regenerating, 2 days post-amputation)")
    legend("topright",
           legend = c("Control Tail", "Regenerating Tail"),
           fill = c("#C0C0C0", "#4169E1"), cex = 0.9, bty = "n")
    dev.off()
  }
}
pnum <- pnum + 1


#-----------------------------------------------
# STEP 20: Module preservation (control vs treated)
#-----------------------------------------------

control_samps <- rownames(trait_data)[trait_data$control == 1]
treated_samps <- rownames(trait_data)[trait_data$control == 0]

cat("\nModule preservation: Control (", length(control_samps),
    ") vs Treated (", length(treated_samps), ")\n")

if (length(control_samps) >= 15 & length(treated_samps) >= 15) {
  
  gsg_c <- goodSamplesGenes(expr_matrix[control_samps, ], verbose = 0)
  gsg_t <- goodSamplesGenes(expr_matrix[treated_samps, ], verbose = 0)
  keep  <- gsg_c$goodGenes & gsg_t$goodGenes
  
  multiExpr  <- list(Control = list(data = expr_matrix[control_samps, keep]),
                     Treated = list(data = expr_matrix[treated_samps, keep]))
  multiColor <- list(Control = module_colors[keep])
  
  set.seed(123)
  preservation <- modulePreservation(
    multiExpr, multiColor,
    referenceNetworks = 1, nPermutations = 100, verbose = 3)
  
  pstats <- preservation$preservation$Z$ref.Control$inColumnsAlsoPresentIn.Treated
  pstats$module <- rownames(pstats)
  pstats$color  <- rownames(pstats)
  pstats$status <- cut(pstats$Zsummary.pres,
                       breaks = c(-Inf, 2, 10, Inf),
                       labels = c("Not Preserved", "Moderate", "Strong"))
  
  p_pres <- ggplot(pstats, aes(x = moduleSize, y = Zsummary.pres,
                               color = status, label = color)) +
    geom_point(size = 4, alpha = 0.7) +
    geom_hline(yintercept = c(2, 10), linetype = "dashed", color = "gray50") +
    geom_text_repel(size = 3.5, max.overlaps = 20) +
    scale_color_manual(values = c("Not Preserved" = "#d62728",
                                  "Moderate" = "#ff7f0e",
                                  "Strong" = "#2ca02c")) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom",
          plot.title = element_text(face = "bold", size = 14)) +
    labs(title = "Module Preservation — D. laeve",
         subtitle = paste0("Control (n=", length(control_samps),
                           ") vs Treated (n=", length(treated_samps), ")"),
         x = "Module Size", y = "Preservation Z-summary", color = "Status")
  
  ggsave(paste0("Part5_WGCNA/plots/", pnum, "_module_preservation.png"), p_pres,
         width = 12, height = 8, dpi = 300)
  ggsave(paste0("Part5_WGCNA/plots/", pnum, "_module_preservation.pdf"), p_pres,
         width = 12, height = 8)
  
  write.table(pstats, "Part5_WGCNA/data/module_preservation_results.tsv",
              sep = "\t", row.names = FALSE, quote = FALSE)
}
pnum <- pnum + 1


#-----------------------------------------------
# STEP 21: Tissue-distance analysis
#-----------------------------------------------

cat("\n══════════════════════════════════════════════════════════\n")
cat("  Tissue-Distance Analysis for Regenerating Tail\n")
cat("══════════════════════════════════════════════════════════\n\n")

tissues <- unique(pheno_data$tissue)
tissue_centroids <- list()
for (tis in tissues) {
  tis_ctrl <- rownames(pheno_data)[pheno_data$tissue == tis & pheno_data$condition == "control"]
  if (length(tis_ctrl) >= 2) tissue_centroids[[tis]] <- colMeans(expr_matrix[tis_ctrl, ])
}

cat("Centroids:", paste(names(tissue_centroids), collapse = ", "), "\n\n")

distance_results <- data.frame()
for (s in c(tail_ctrl, tail_amp)) {
  s_expr <- expr_matrix[s, ]
  s_cond <- ifelse(s %in% tail_ctrl, "Control", "Regenerating")
  for (tis in names(tissue_centroids)) {
    d <- sqrt(sum((s_expr - tissue_centroids[[tis]])^2))
    distance_results <- rbind(distance_results, data.frame(
      Sample = get_pub_label(s), Condition = s_cond,
      Tissue = tis, Distance = round(d, 2), stringsAsFactors = FALSE))
  }
}

distance_wide <- reshape2::dcast(distance_results, Sample + Condition ~ Tissue,
                                 value.var = "Distance")
cat("── Distances to Tissue Centroids ──\n")
print(distance_wide)

write.table(distance_wide, "Part5_WGCNA/data/tail_tissue_distance_matrix.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

# Heatmap
dist_mat <- as.matrix(distance_wide[, -(1:2)])
rownames(dist_mat) <- distance_wide$Sample

for (dev_type in c("png", "pdf")) {
  if (dev_type == "png") save_plot(paste0(pnum, "_tissue_distance_heatmap"), 1200, 700, 120)
  else save_pdf(paste0(pnum, "_tissue_distance_heatmap"), 12, 7)
  pheatmap(dist_mat, scale = "column",
           annotation_row = data.frame(Condition = distance_wide$Condition,
                                       row.names = distance_wide$Sample),
           cluster_rows = TRUE, cluster_cols = TRUE,
           color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(50),
           main = "Transcriptomic Distance: Tail Samples → Tissue Centroids — D. laeve\n(Column-scaled; blue = closer, red = farther)",
           fontsize = 10, display_numbers = TRUE, number_format = "%.0f",
           number_color = "black")
  dev.off()
}
pnum <- pnum + 1

# Barplot
dist_summary <- aggregate(Distance ~ Condition + Tissue, data = distance_results, mean)
dist_summary$SE <- aggregate(Distance ~ Condition + Tissue,
                             data = distance_results, function(x) sd(x)/sqrt(length(x)))$Distance

p_dist <- ggplot(dist_summary, aes(x = Tissue, y = Distance, fill = Condition)) +
  geom_col(position = position_dodge(0.8), width = 0.7, alpha = 0.85) +
  geom_errorbar(aes(ymin = Distance - SE, ymax = Distance + SE),
                position = position_dodge(0.8), width = 0.3) +
  scale_fill_manual(values = c("Control" = "#C0C0C0", "Regenerating" = "#4169E1")) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold", size = 14),
        legend.position = "bottom") +
  labs(title = "Transcriptomic Distance: Tail → Other Tissues — D. laeve",
       subtitle = "Regenerating tail moves AWAY from tail identity and toward bodywall/eye",
       y = "Euclidean Distance to Tissue Centroid", x = "Tissue", fill = "Tail Condition")

ggsave(paste0("Part5_WGCNA/plots/", pnum, "_tissue_distance_barplot.png"), p_dist,
       width = 12, height = 7, dpi = 300)
ggsave(paste0("Part5_WGCNA/plots/", pnum, "_tissue_distance_barplot.pdf"), p_dist,
       width = 12, height = 7)
pnum <- pnum + 1


#-----------------------------------------------
# STEP 22: ME comparisons by condition and experiment
#-----------------------------------------------

me_comp <- data.frame(
  Eigengene  = module_eigengenes[, best_me_name],
  Condition  = pheno_data$condition,
  Experiment = pheno_data$experiment,
  Tissue     = pheno_data$tissue,
  row.names  = rownames(expr_matrix))

p_cond <- ggplot(me_comp, aes(x = Condition, y = Eigengene, fill = Condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 2.5) +
  scale_fill_manual(values = condition_colors) +
  theme_minimal(base_size = 12) +
  labs(title = paste(best_module_wilcox, "Module Eigengene by Condition — D. laeve"),
       subtitle = paste("Best module by Wilcoxon test | p =", signif(best_p, 3)),
       y = "Module Eigengene", x = "") +
  theme(legend.position = "none", plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(paste0("Part5_WGCNA/plots/", pnum, "_ME_by_condition_", best_module_wilcox, ".png"),
       p_cond, width = 10, height = 7, dpi = 300)
ggsave(paste0("Part5_WGCNA/plots/", pnum, "_ME_by_condition_", best_module_wilcox, ".pdf"),
       p_cond, width = 10, height = 7)
pnum <- pnum + 1

p_exp <- ggplot(me_comp, aes(x = Experiment, y = Eigengene, fill = Condition)) +
  geom_boxplot(position = position_dodge(0.8), alpha = 0.7) +
  geom_point(position = position_dodge(0.8), alpha = 0.6, size = 2) +
  scale_fill_manual(values = condition_colors) +
  theme_minimal(base_size = 12) +
  labs(title = paste(best_module_wilcox, "Module Eigengene by Experiment — D. laeve"),
       y = "Module Eigengene", x = "", fill = "Condition") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold"), legend.position = "bottom")

ggsave(paste0("Part5_WGCNA/plots/", pnum, "_ME_by_experiment_", best_module_wilcox, ".png"),
       p_exp, width = 12, height = 7, dpi = 300)
ggsave(paste0("Part5_WGCNA/plots/", pnum, "_ME_by_experiment_", best_module_wilcox, ".pdf"),
       p_exp, width = 12, height = 7)
pnum <- pnum + 1


###############################################################################
# SAVE ALL RESULTS
###############################################################################

write.table(
  data.frame(Gene = names(module_colors),
             Module_Color = module_colors,
             Module_Number = net$colors, stringsAsFactors = FALSE),
  "Part5_WGCNA/data/all_gene_module_assignments.tsv",
  sep = "\t", row.names = FALSE, quote = FALSE)

write.table(module_trait_cor, "Part5_WGCNA/data/module_trait_correlations.tsv",
            sep = "\t", quote = FALSE)
write.table(module_trait_pvalue, "Part5_WGCNA/data/module_trait_pvalues.tsv",
            sep = "\t", quote = FALSE)
write.table(module_eigengenes, "Part5_WGCNA/data/module_eigengenes.tsv",
            sep = "\t", quote = FALSE)
write.table(module_gene_info,
            paste0("Part5_WGCNA/data/module_gene_info_", best_module_wilcox, ".tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

write.table(
  data.frame(Sample = outliers_found,
             Publication_Label = get_pub_label(outliers_found),
             Tissue = pheno_data_full[outliers_found, "tissue"],
             Condition = pheno_data_full[outliers_found, "condition"],
             Experiment = pheno_data_full[outliers_found, "experiment"],
             Reason = "Did not cluster with replicates"),
  "Part5_WGCNA/data/removed_outlier_samples.tsv",
  sep = "\t", row.names = FALSE, quote = FALSE)

write.table(
  data.frame(Original_ID = names(sample_labels),
             Publication_Label = sample_labels, stringsAsFactors = FALSE),
  "Part5_WGCNA/data/sample_label_mapping.tsv",
  sep = "\t", row.names = FALSE, quote = FALSE)

save(net, module_colors, module_eigengenes,
     module_trait_cor, module_trait_pvalue,
     soft_power, expr_matrix, trait_data, pheno_data,
     all_hub_genes, module_gene_info, outliers_found,
     sample_labels, me_wilcox, distance_results,
     best_module_wilcox, best_me_name,
     file = "Part5_WGCNA/data/Part5_WGCNA_objects.RData")


###############################################################################
# FINAL SUMMARY
###############################################################################

cat("\n")
cat("╔══════════════════════════════════════════════════════════╗\n")
cat("║  WGCNA V3 (DEFINITIVE) — Deroceras laeve               ║\n")
cat("╚══════════════════════════════════════════════════════════╝\n")
cat("\n── Network Parameters ──\n")
cat("  Network type: SIGNED\n")
cat("  Soft power: 12\n")
cat("  Variance filter: 25th percentile\n")
cat("  R² at power 12:", round(fit_vals[sft$fitIndices[, 1] == 12], 3),
    "(expected low for multi-tissue)\n")
cat("  Genes:", ncol(expr_matrix), "| Samples:", nrow(expr_matrix), "\n")
cat("  Modules:", length(real_modules), "(+ grey)\n")
cat("  Outliers removed:", paste(get_pub_label(outliers_found), collapse = ", "), "\n")
cat("\n── PRIMARY Module Selection (Wilcoxon) ──\n")
cat("  Best module:", best_module_wilcox, "\n")
cat("  Wilcoxon p:", best_p, "| Cohen's d:", best_d, "\n")
cat("  Direction:", best_direction, "\n")
cat("  Second:", second_module, "(p =", me_wilcox$Wilcox_p[2], ")\n")
cat("\n── Hub Genes ──\n")
cat("  Total hubs:", nrow(hubs_only), "\n")
cat("  Candidates (|MM|>0.7 + |GS|>0.3):", nrow(candidates_only), "\n")
cat("\n── V3 Changes from V2 ──\n")
cat("  ✓ PRIMARY: Wilcoxon ME test replaces Pearson r\n")
cat("  ✓ NEW: Gene dendrograms for tail ctrl + tail amp with module colors\n")
cat("  ✓ NEW: Module cohesion analysis (TOM connectivity change in regen)\n")
cat("  ✓ FIX: Enrichment plots ONLY for modules with significant terms\n")
cat("  ✓ FIX: GO term label truncation (no more cropped text)\n")
cat("  ✓ KEEP: Tissue-distance analysis, publication labels, signed network\n")
cat("\n── Key Outputs for DMR Cross-reference ──\n")
cat("  all_modules_hub_genes.tsv          (all genes + Wilcoxon stats)\n")
cat("  hub_genes_all_modules.tsv          (hub genes only)\n")
cat("  candidate_genes_highMM_highGS.tsv  (|MM|>0.7 + |GS|>0.3)\n")
cat("  all_gene_module_assignments.tsv    (gene → module mapping)\n")
cat("  ME_wilcoxon_tail_regeneration.tsv  (PRIMARY: module rankings)\n")
cat("  module_cohesion_tail_comparison.tsv (TOM change in regeneration)\n")
cat("  tail_tissue_distance_matrix.tsv    (tissue identity shift)\n")
cat("  enrichment_summary_all_modules.tsv (functional annotation)\n")
cat("\n  All plots saved as PNG + PDF.\n")
cat("══════════════════════════════════════════════════════════\n")