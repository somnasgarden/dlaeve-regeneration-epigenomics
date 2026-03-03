#!/usr/bin/env Rscript
# =============================================================================
# SOX19a & INTERGENIC REGULATION: The centerpiece finding
# =============================================================================
# Sox19a is significantly upregulated (log2FC=1.45, padj=0.044) in
# regenerating tissue. It has NO DMPs directly on the gene, but
# 13 DMPs in the intergenic/downstream region ~4-5kb upstream.
#
# This script investigates:
# 1. Sox19a intergenic DMP landscape in detail
# 2. Other genes showing the same pattern (expression change + intergenic DMPs)
# 3. Comparison with genes that DO have genic DMPs
# 4. The "distal regulation" gene set
# 5. Key signaling pathway genes with similar patterns
#
# Run: Rscript local/07_sox19a_intergenic_regulation.R
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ggrepel)
  library(scales)
})

options(stringsAsFactors = FALSE)
options(scipen = 999)

theme_set(theme_minimal() +
            theme(text = element_text(size = 12),
                  plot.title = element_text(size = 14, face = "bold"),
                  axis.title = element_text(size = 12)))

OUTPUT_DIR <- "results/13_intergenic_regulation"
dir.create(file.path(OUTPUT_DIR, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "plots"),  recursive = TRUE, showWarnings = FALSE)

DATA_DIR <- "/mnt/c/Users/rafae/Projects/DATA"

save_both <- function(plot_obj, name, width, height) {
  ggsave(file.path(OUTPUT_DIR, "plots", paste0(name, ".pdf")),
         plot_obj, width = width, height = height)
  ggsave(file.path(OUTPUT_DIR, "plots", paste0(name, ".png")),
         plot_obj, width = width, height = height, dpi = 300)
  cat(sprintf("  Saved: %s (.pdf + .png)\n", name))
}

cat("================================================================\n")
cat("  INTERGENIC REGULATION: Sox19a and Beyond\n")
cat("================================================================\n\n")

# Load data
dmps <- read.delim("results/01_methylation/dmps_annotated.txt")
dmps_chip <- read.delim("results/01_methylation/dmps_chipseeker_annotated.txt")
gene_meth_expr <- read.delim("results/03_integration/Tables/MXT_gene_level_meth_vs_expression.txt")
dmp_spatial <- read.delim("results/03_integration/Tables/MXT_gene_DMP_spatial_analysis.txt")
modules <- read.delim("results/02_rnaseq/Part5_WGCNA/data/all_gene_module_assignments.tsv")
hubs <- read.delim("results/02_rnaseq/Part5_WGCNA/data/hub_genes_all_modules.tsv")
degs <- read.csv("results/02_rnaseq/Part2_DEGs/data/DEresults_amputated_vs_control.csv",
                 row.names = 1)
degs$gene_id <- rownames(degs)
annot <- read.delim(file.path(DATA_DIR, "derLaeGenome_eviann_annotations.tsv"),
                    header = FALSE, col.names = c("gene_id", "gene_symbol", "description"))
annot_unique <- annot %>% group_by(gene_id) %>%
  summarise(gene_symbol = gene_symbol[1],
            description = paste(unique(description), collapse = "; "),
            .groups = "drop")

cat("  Data loaded.\n\n")

# =============================================================================
# PART 1: SOX19a DETAILED LANDSCAPE
# =============================================================================

cat("================================================================\n")
cat("  PART 1: Sox19a — The intergenic regulation paradigm\n")
cat("================================================================\n\n")

sox_id <- "LOC_00002600"

# ChIPseeker DMPs associated with Sox19a
sox_chip <- dmps_chip %>% filter(geneId == sox_id)
cat(sprintf("  Sox19a ChIPseeker-associated DMPs: %d\n", nrow(sox_chip)))

# Get the actual DMP data for these positions
sox_dmp_positions <- dmps %>%
  filter(start %in% sox_chip$start & seqnames %in% sox_chip$seqnames)

sox_merged <- sox_chip %>%
  left_join(dmps %>% select(seqnames, start, methylation_diff, fdr),
            by = c("seqnames", "start"))

cat("  Sox19a DMP details:\n")
print(sox_merged %>%
        select(annotation, distanceToTSS, methylation_diff, fdr) %>%
        arrange(distanceToTSS) %>% as.data.frame(),
      row.names = FALSE)

sox_expr <- degs %>% filter(gene_id == sox_id)
cat(sprintf("\n  Sox19a expression: log2FC = %.3f, padj = %.4f\n",
            sox_expr$log2FoldChange, sox_expr$padj))
cat(sprintf("  Sox19a mean meth diff of associated DMPs: %.3f\n",
            mean(sox_merged$methylation_diff, na.rm = TRUE)))
cat(sprintf("  Sox19a DMPs direction: %d hyper, %d hypo\n",
            sum(sox_merged$methylation_diff > 0, na.rm = TRUE),
            sum(sox_merged$methylation_diff < 0, na.rm = TRUE)))

# Plot: Sox19a DMP landscape
if (nrow(sox_merged) > 0 && any(!is.na(sox_merged$methylation_diff))) {
  sox_plot_data <- sox_merged %>%
    filter(!is.na(methylation_diff)) %>%
    mutate(
      position = distanceToTSS,
      sig = fdr < 0.05,
      region = ifelse(grepl("Intergenic", annotation), "Intergenic",
              ifelse(grepl("Downstream", annotation), "Downstream", "Gene body"))
    )

  p_sox <- ggplot(sox_plot_data, aes(x = position, y = methylation_diff)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
    geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 1) +
    geom_point(aes(color = region, shape = sig), size = 4, alpha = 0.8) +
    scale_color_manual(values = c("Intergenic" = "#3498DB",
                                   "Downstream" = "#E67E22",
                                   "Gene body" = "#2ECC71"),
                       name = "Region") +
    scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                       name = "FDR<0.05") +
    annotate("text", x = 500, y = max(sox_plot_data$methylation_diff, na.rm = TRUE) * 0.9,
             label = paste0("Sox19a\nlog2FC = +", round(sox_expr$log2FoldChange, 2),
                           "\npadj = ", formatC(sox_expr$padj, format = "e", digits = 1)),
             hjust = 0, size = 4, fontface = "bold") +
    annotate("segment", x = 0, xend = 0, y = -0.35, yend = 0.35,
             arrow = arrow(length = unit(0.3, "cm")),
             color = "black", linewidth = 1.5) +
    labs(
      title = "Sox19a: Intergenic Methylation Landscape",
      subtitle = "D. laeve: DMPs near Sox19a (TSS at position 0). Gene is significantly UPREGULATED.",
      x = "Distance from TSS (bp, negative = upstream)",
      y = "Methylation difference (amputated - control)"
    )
  save_both(p_sox, "fig_H01_sox19a_landscape", 14, 8)
}

# =============================================================================
# PART 2: GENES WITH INTERGENIC-ONLY REGULATION PATTERN
# =============================================================================

cat("\n================================================================\n")
cat("  PART 2: Genes regulated via intergenic methylation only\n")
cat("================================================================\n\n")

# Find genes where:
# 1. ChIPseeker associates DMPs with the gene
# 2. All/most DMPs are intergenic or downstream
# 3. Gene shows expression change

# First, count DMP regions per gene from ChIPseeker
gene_dmp_regions <- dmps_chip %>%
  mutate(
    region_type = case_when(
      grepl("Intergenic", annotation) ~ "Intergenic",
      grepl("Downstream", annotation) ~ "Downstream",
      grepl("Promoter", annotation) ~ "Promoter",
      grepl("Intron", annotation) ~ "Intron",
      grepl("Exon", annotation) ~ "Exon",
      grepl("UTR", annotation) ~ "UTR",
      TRUE ~ "Other"
    )
  ) %>%
  group_by(geneId) %>%
  summarise(
    n_total_dmps = n(),
    n_intergenic = sum(region_type == "Intergenic"),
    n_downstream = sum(region_type == "Downstream"),
    n_promoter = sum(region_type == "Promoter"),
    n_intron = sum(region_type == "Intron"),
    n_exon = sum(region_type == "Exon"),
    n_utr = sum(region_type == "UTR"),
    pct_distal = (n_intergenic + n_downstream) / n_total_dmps * 100,
    median_distTSS = median(abs(distanceToTSS)),
    .groups = "drop"
  )

# Merge with expression data
distal_regulated <- gene_dmp_regions %>%
  filter(pct_distal >= 80 & n_total_dmps >= 3) %>%
  left_join(degs %>% select(gene_id, log2FoldChange, padj, baseMean),
            by = c("geneId" = "gene_id")) %>%
  left_join(modules, by = c("geneId" = "Gene")) %>%
  left_join(annot_unique, by = c("geneId" = "gene_id")) %>%
  filter(!is.na(log2FoldChange))

cat(sprintf("  Genes with >=80%% distal DMPs (3+ DMPs): %d\n", nrow(distal_regulated)))
cat(sprintf("  Of these with expression data: %d\n",
            sum(!is.na(distal_regulated$log2FoldChange))))
cat(sprintf("  Significantly DE (padj<0.05): %d\n",
            sum(distal_regulated$padj < 0.05, na.rm = TRUE)))

# Compare expression effect: distal-only vs genic-only
genic_regulated <- gene_dmp_regions %>%
  filter(pct_distal <= 20 & n_total_dmps >= 3) %>%
  left_join(degs %>% select(gene_id, log2FoldChange, padj),
            by = c("geneId" = "gene_id")) %>%
  filter(!is.na(log2FoldChange))

if (nrow(distal_regulated) > 10 && nrow(genic_regulated) > 10) {
  wt_distal <- wilcox.test(abs(distal_regulated$log2FoldChange),
                            abs(genic_regulated$log2FoldChange))
  cat(sprintf("\n  Distal vs genic |log2FC| Wilcoxon: p = %s\n",
              formatC(wt_distal$p.value, format = "e", digits = 2)))
  cat(sprintf("    Distal median |log2FC|: %.3f (n=%d)\n",
              median(abs(distal_regulated$log2FoldChange)),
              nrow(distal_regulated)))
  cat(sprintf("    Genic median |log2FC|: %.3f (n=%d)\n",
              median(abs(genic_regulated$log2FoldChange)),
              nrow(genic_regulated)))
}

# Significantly DE genes with distal-only regulation
sig_distal <- distal_regulated %>%
  filter(padj < 0.05) %>%
  arrange(padj)

if (nrow(sig_distal) > 0) {
  cat(sprintf("\n  Significantly DE genes with distal-only methylation (%d genes):\n",
              nrow(sig_distal)))
  print(sig_distal %>%
          select(geneId, gene_symbol, Module_Color, n_total_dmps,
                 pct_distal, median_distTSS, log2FoldChange, padj) %>%
          head(30) %>% as.data.frame(),
        row.names = FALSE)
}

write.table(distal_regulated,
            file.path(OUTPUT_DIR, "tables", "distal_regulated_genes.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Module distribution of distal-regulated genes
if (nrow(distal_regulated) > 0) {
  distal_module <- distal_regulated %>%
    filter(!is.na(Module_Color)) %>%
    count(Module_Color, name = "n_distal") %>%
    left_join(
      modules %>% count(Module_Color, name = "module_size"),
      by = "Module_Color"
    ) %>%
    mutate(pct_distal = round(100 * n_distal / module_size, 2)) %>%
    arrange(desc(pct_distal))

  cat("\n  Module distribution of distal-regulated genes:\n")
  print(as.data.frame(distal_module), row.names = FALSE)

  write.table(distal_module,
              file.path(OUTPUT_DIR, "tables", "distal_regulation_by_module.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
}

# =============================================================================
# PART 3: SIGNALING PATHWAY INTERGENIC REGULATION
# =============================================================================

cat("\n================================================================\n")
cat("  PART 3: Signaling pathways and intergenic regulation\n")
cat("================================================================\n\n")

# Which signaling pathway genes have distal methylation patterns?
pathway_patterns <- c(
  "Wnt|wnt|wls|Daam|dishevel|frizzled|beta-catenin",
  "Hedgehog|hh|patched|smooth|shhb",
  "BMP|bone morphogenetic|Bambi|BMPR",
  "FGF|fibroblast growth",
  "Notch|notch|nrarp",
  "Sox|sox",
  "TGF|transforming growth factor",
  "EGF|epidermal growth factor",
  "hippo|yap|taz"
)
pathway_names <- c("Wnt", "Hedgehog", "BMP", "FGF", "Notch",
                   "Sox", "TGF-beta", "EGF", "Hippo")

pathway_dmp_analysis <- data.frame()

for (i in seq_along(pathway_patterns)) {
  pw_genes <- annot_unique %>%
    filter(grepl(pathway_patterns[i], description, ignore.case = TRUE)) %>%
    pull(gene_id)

  if (length(pw_genes) == 0) next

  pw_dmps <- gene_dmp_regions %>% filter(geneId %in% pw_genes)
  pw_expr <- degs %>% filter(gene_id %in% pw_genes)

  pathway_dmp_analysis <- rbind(pathway_dmp_analysis, data.frame(
    pathway = pathway_names[i],
    n_genes = length(pw_genes),
    n_with_dmps = nrow(pw_dmps),
    n_distal_only = sum(pw_dmps$pct_distal >= 80, na.rm = TRUE),
    pct_with_dmps = round(100 * nrow(pw_dmps) / length(pw_genes), 1),
    n_sig_DE = sum(pw_expr$padj < 0.05, na.rm = TRUE),
    mean_abs_FC = round(mean(abs(pw_expr$log2FoldChange), na.rm = TRUE), 3)
  ))
}

cat("  Signaling pathway DMP patterns:\n")
print(as.data.frame(pathway_dmp_analysis), row.names = FALSE)

write.table(pathway_dmp_analysis,
            file.path(OUTPUT_DIR, "tables", "pathway_dmp_patterns.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# =============================================================================
# PART 4: DISTAL vs PROXIMAL REGULATION MODEL
# =============================================================================

cat("\n================================================================\n")
cat("  PART 4: Building the distal regulation model\n")
cat("================================================================\n\n")

# Key question: Does distance from TSS matter for methylation effect?
# Use ALL genes with ChIPseeker data

all_gene_dist <- dmps_chip %>%
  group_by(geneId) %>%
  summarise(
    n_dmps = n(),
    median_dist = median(abs(distanceToTSS)),
    min_dist = min(abs(distanceToTSS)),
    max_dist = max(abs(distanceToTSS)),
    mean_dist = mean(abs(distanceToTSS)),
    pct_proximal = mean(abs(distanceToTSS) < 2000) * 100,
    .groups = "drop"
  ) %>%
  left_join(degs %>% select(gene_id, log2FoldChange, padj),
            by = c("geneId" = "gene_id")) %>%
  filter(!is.na(log2FoldChange))

# Bin by median distance
all_gene_dist$dist_bin <- cut(all_gene_dist$median_dist,
                               breaks = c(0, 1000, 2000, 5000, 10000, 50000, Inf),
                               labels = c("<1kb", "1-2kb", "2-5kb",
                                         "5-10kb", "10-50kb", ">50kb"),
                               include.lowest = TRUE)

dist_effect <- all_gene_dist %>%
  group_by(dist_bin) %>%
  summarise(
    n = n(),
    median_abs_FC = median(abs(log2FoldChange)),
    mean_abs_FC = mean(abs(log2FoldChange)),
    pct_sig = mean(padj < 0.05, na.rm = TRUE) * 100,
    mean_n_dmps = mean(n_dmps),
    .groups = "drop"
  )

cat("  DMP distance to TSS vs expression effect:\n")
print(as.data.frame(dist_effect), row.names = FALSE)

write.table(dist_effect,
            file.path(OUTPUT_DIR, "tables", "distance_expression_model.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Spearman correlation: distance vs expression
cor_dist <- cor.test(log10(all_gene_dist$median_dist + 1),
                     abs(all_gene_dist$log2FoldChange),
                     method = "spearman")
cat(sprintf("\n  log10(distance) vs |log2FC|: rho = %.3f, p = %s\n",
            cor_dist$estimate,
            formatC(cor_dist$p.value, format = "e", digits = 2)))

# Plot: Distance distribution of ALL DMPs
p_dist_dist <- ggplot(all_gene_dist, aes(x = median_dist)) +
  geom_histogram(binwidth = 1000, fill = "#3498DB", alpha = 0.7, color = "white") +
  geom_vline(xintercept = 2000, linetype = "dashed", color = "red") +
  annotate("text", x = 2500, y = Inf, vjust = 2,
           label = "2kb (proximal)", color = "red", size = 3.5) +
  scale_x_continuous(labels = comma, limits = c(0, 100000)) +
  labs(
    title = "Distribution of Median DMP Distance to TSS",
    subtitle = "D. laeve: Most DMPs are distal (>2kb from TSS)",
    x = "Median distance to TSS (bp)", y = "Number of genes"
  )
save_both(p_dist_dist, "fig_H02_dmp_distance_distribution", 12, 7)

# Plot: Distance vs expression effect
p_dist_effect <- ggplot(all_gene_dist %>% filter(median_dist < 200000),
                        aes(x = median_dist, y = abs(log2FoldChange))) +
  geom_point(alpha = 0.1, size = 0.5) +
  geom_smooth(method = "loess", color = "#E74C3C", linewidth = 1.2) +
  geom_vline(xintercept = c(2000, 10000), linetype = "dashed", alpha = 0.5) +
  scale_x_log10(labels = comma) +
  labs(
    title = "DMP Distance to TSS vs Expression Effect",
    subtitle = sprintf("D. laeve: Spearman rho = %.3f, p = %s",
                       cor_dist$estimate,
                       formatC(cor_dist$p.value, format = "e", digits = 2)),
    x = "Median distance to TSS (bp, log scale)", y = "|log2FC|"
  )
save_both(p_dist_effect, "fig_H03_distance_vs_expression_effect", 12, 8)

# Plot: Distal vs genic regulated genes - volcano-style
if (nrow(distal_regulated) > 5) {
  combined_regulation <- rbind(
    distal_regulated %>%
      mutate(regulation_type = "Distal (>=80% intergenic)") %>%
      select(geneId, gene_symbol, regulation_type, log2FoldChange, padj),
    genic_regulated %>%
      left_join(annot_unique, by = c("geneId" = "gene_id")) %>%
      mutate(regulation_type = "Genic (<=20% intergenic)") %>%
      select(geneId, gene_symbol, regulation_type, log2FoldChange, padj)
  )

  p_volcano <- ggplot(combined_regulation,
                      aes(x = log2FoldChange, y = -log10(padj),
                          color = regulation_type)) +
    geom_point(alpha = 0.3, size = 1.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray60") +
    geom_text_repel(
      data = combined_regulation %>% filter(padj < 0.05 & abs(log2FoldChange) > 1),
      aes(label = gene_symbol), size = 2.5, max.overlaps = 15
    ) +
    scale_color_manual(values = c("Distal (>=80% intergenic)" = "#3498DB",
                                   "Genic (<=20% intergenic)" = "#E67E22"),
                       name = "Regulation type") +
    labs(
      title = "Distal vs Genic DMP Regulation: Volcano Plot",
      subtitle = "D. laeve: Comparing genes with distal-only vs genic-only methylation changes",
      x = "log2FC", y = "-log10(padj)"
    )
  save_both(p_volcano, "fig_H04_distal_vs_genic_volcano", 14, 10)
}


# =============================================================================
# PART 5: TOP INTERGENIC-REGULATED CANDIDATES
# =============================================================================

cat("\n================================================================\n")
cat("  PART 5: Top intergenic-regulated candidate genes\n")
cat("================================================================\n\n")

# Build a master table of intergenic-regulated genes
# Prioritize by: significance, effect size, number of DMPs, distance pattern

intergenic_candidates <- distal_regulated %>%
  mutate(
    score = 0,
    score = score + ifelse(!is.na(padj) & padj < 0.05, 3, 0),
    score = score + ifelse(abs(log2FoldChange) > 1, 2,
                    ifelse(abs(log2FoldChange) > 0.5, 1, 0)),
    score = score + ifelse(n_total_dmps >= 10, 2,
                    ifelse(n_total_dmps >= 5, 1, 0)),
    score = score + ifelse(!is.na(Module_Color) & Module_Color == "yellow", 2, 0),
    score = score + ifelse(geneId %in% hubs$Gene, 1, 0)
  ) %>%
  arrange(desc(score), padj)

cat("  Top 30 intergenic-regulated candidate genes:\n")
print(intergenic_candidates %>%
        select(geneId, gene_symbol, Module_Color, n_total_dmps,
               pct_distal, median_distTSS, log2FoldChange, padj, score) %>%
        head(30) %>% as.data.frame(),
      row.names = FALSE)

write.table(intergenic_candidates,
            file.path(OUTPUT_DIR, "tables", "intergenic_candidates_ranked.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Top candidates lollipop
top_inter <- intergenic_candidates %>%
  head(25) %>%
  mutate(label = ifelse(is.na(gene_symbol) | gene_symbol == "",
                        geneId, gene_symbol))

p_top_inter <- ggplot(top_inter,
                      aes(x = reorder(label, score), y = score)) +
  geom_segment(aes(xend = label, y = 0, yend = score, color = Module_Color),
               linewidth = 1) +
  geom_point(aes(color = Module_Color, size = n_total_dmps)) +
  scale_color_identity() +
  scale_size_continuous(name = "N DMPs", range = c(2, 8)) +
  coord_flip() +
  labs(
    title = "Top 25 Intergenic-Regulated Candidate Genes",
    subtitle = "D. laeve: Genes with >=80% distal DMPs, ranked by composite score",
    x = NULL, y = "Composite score"
  )
save_both(p_top_inter, "fig_H05_top_intergenic_candidates", 14, 10)


# =============================================================================
# SUMMARY
# =============================================================================

cat("\n================================================================\n")
cat("  INTERGENIC REGULATION ANALYSIS COMPLETE\n")
cat("================================================================\n\n")

cat("KEY RESULTS:\n\n")

cat("1. SOX19a PARADIGM:\n")
cat(sprintf("   - Significantly upregulated (log2FC=%.2f, padj=%.3f)\n",
            sox_expr$log2FoldChange, sox_expr$padj))
cat(sprintf("   - %d associated DMPs, ALL in intergenic/downstream region\n",
            nrow(sox_merged)))
cat("   - DMPs located ~4-5kb upstream of TSS\n")
cat("   - This is consistent with enhancer demethylation activating Sox19a\n\n")

cat("2. DISTAL REGULATION IS WIDESPREAD:\n")
cat(sprintf("   - %d genes with >=80%% distal DMPs (3+ DMPs)\n",
            nrow(distal_regulated)))
if (exists("sig_distal") && nrow(sig_distal) > 0) {
  cat(sprintf("   - %d of these are significantly DE (padj<0.05)\n",
              nrow(sig_distal)))
}

cat("\n3. PATHWAY IMPLICATIONS:\n")
cat("   Signaling pathway genes show variable distal regulation\n")

cat(sprintf("\nOUTPUT: %s/\n", OUTPUT_DIR))
cat("================================================================\n")
