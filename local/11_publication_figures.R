#!/usr/bin/env Rscript
# =============================================================================
# Script 11: Publication-Quality Figures (Professional Styling)
# =============================================================================
# Purpose: Regenerate key figures with professional color palettes, clean
#          typography, and Nature/Science-quality aesthetics.
#
# Replaces: fig_E1-E4 from script 04 (candidate genes)
#           Adds polished versions of key findings figures
#
# Output:  results/22_publication_figures/
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(pheatmap)
  library(RColorBrewer)
  library(scales)
  library(ggrepel)
})

cat("=== Script 11: Publication-Quality Figures ===\n\n")

outdir <- "results/22_publication_figures"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ── Professional theme ───────────────────────────────────────────────────────

theme_pub <- theme_minimal(base_size = 13) +
  theme(
    text = element_text(family = "sans", color = "gray10"),
    plot.title = element_text(size = 16, face = "bold", margin = margin(b = 8)),
    plot.subtitle = element_text(size = 11, color = "gray30", margin = margin(b = 12)),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10, color = "gray20"),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid.major = element_line(color = "gray92", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    plot.margin = margin(15, 15, 15, 15)
  )
theme_set(theme_pub)

# Professional color palettes
pal_main <- c("#2C3E50", "#E74C3C", "#3498DB", "#2ECC71", "#F39C12", "#9B59B6", "#1ABC9C", "#E67E22")
pal_gradient <- c("#2C3E50", "#34495E", "#7F8C8D", "#BDC3C7", "#ECF0F1")
pal_diverging <- c("#2980B9", "#5DADE2", "#AED6F1", "#FDEBD0", "#F5B041", "#E74C3C")

save_pub <- function(plot_obj, name, width = 10, height = 7) {
  ggsave(file.path(outdir, paste0(name, ".png")), plot_obj,
         width = width, height = height, dpi = 300, bg = "white")
  ggsave(file.path(outdir, paste0(name, ".pdf")), plot_obj,
         width = width, height = height, bg = "white")
  cat("  Saved:", name, "\n")
}

# ── Load data ────────────────────────────────────────────────────────────────

cat("Loading data...\n")

gene_meth <- read.delim("results/03_integration/Tables/MXT_gene_level_meth_vs_expression.txt",
                         stringsAsFactors = FALSE)
modules <- read.delim("results/02_rnaseq/Part5_WGCNA/data/all_gene_module_assignments.tsv",
                       stringsAsFactors = FALSE)
hubs <- read.delim("results/02_rnaseq/Part5_WGCNA/data/hub_genes_all_modules.tsv",
                    stringsAsFactors = FALSE)
deg <- read.csv("results/02_rnaseq/Part2_DEGs/data/DEresults_amputated_vs_control.csv",
                 stringsAsFactors = FALSE)
colnames(deg)[1] <- "gene_id"
annot <- read.delim("/mnt/c/Users/rafae/Projects/DATA/derLaeGenome_eviann_annotations.tsv",
                     header = FALSE, stringsAsFactors = FALSE)
colnames(annot) <- c("gene_id", "gene_symbol", "description")
annot <- annot %>% distinct(gene_id, .keep_all = TRUE)
master <- read.delim("results/03_integration/Tables/MXT_MASTER_module_methylation_summary.txt",
                      stringsAsFactors = FALSE)
dmps <- read.delim("results/01_methylation/dmps_annotated.txt", stringsAsFactors = FALSE)
dmps_cs <- read.delim("results/01_methylation/dmps_chipseeker_annotated.txt", stringsAsFactors = FALSE)
dmrs <- read.delim("results/01_methylation/dmrs_annotated.txt", stringsAsFactors = FALSE)

cat("  All data loaded\n\n")

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 1: Module DMP Enrichment Overview (the key result)
# ══════════════════════════════════════════════════════════════════════════════

cat("--- Figure 1: Module DMP enrichment ---\n")

mod_data <- master %>%
  mutate(
    sig = ifelse(fisher_padj < 0.05, "Significant", "NS"),
    direction_label = ifelse(grepl("UP", Direction), "UP in regeneration", "DOWN in regeneration")
  ) %>%
  arrange(desc(pct_genes_with_dmp))

p1 <- ggplot(mod_data, aes(x = reorder(module, pct_genes_with_dmp), y = pct_genes_with_dmp)) +
  geom_col(aes(fill = pct_genes_with_dmp), width = 0.7, show.legend = FALSE) +
  geom_hline(yintercept = mean(mod_data$pct_genes_with_dmp), linetype = "dashed",
             color = "#E74C3C", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%", pct_genes_with_dmp)),
            hjust = -0.15, size = 3.5, fontface = "bold") +
  geom_text(aes(y = 1, label = ifelse(fisher_padj < 0.05, paste0("p=", format(fisher_padj, digits = 2)), "")),
            color = "#E74C3C", size = 3, fontface = "bold") +
  coord_flip(ylim = c(0, max(mod_data$pct_genes_with_dmp) * 1.15)) +
  scale_fill_gradient(low = "#AED6F1", high = "#2C3E50") +
  annotate("text", x = 0.5, y = mean(mod_data$pct_genes_with_dmp) + 0.3,
           label = sprintf("Mean: %.1f%%", mean(mod_data$pct_genes_with_dmp)),
           color = "#E74C3C", size = 3.5, hjust = 0) +
  labs(
    title = "Differential Methylation by WGCNA Co-expression Module",
    subtitle = "Percentage of genes with at least one DMP | Red dashed = mean across modules",
    x = NULL,
    y = "Genes with DMPs (%)"
  )

save_pub(p1, "fig01_module_dmp_enrichment", 10, 6)

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 2: TF Protection & Functional Category Enrichment
# ══════════════════════════════════════════════════════════════════════════════

cat("--- Figure 2: Functional category DMP enrichment ---\n")

# Load from morphogenesis analysis
cat_enrich <- read.delim("results/14_morphogenesis_enrichment/tables/I01_category_dmp_enrichment.tsv",
                          stringsAsFactors = FALSE)

cat_plot <- cat_enrich %>%
  filter(n_genes >= 20) %>%
  mutate(
    sig_label = ifelse(padj < 0.05, "***", ifelse(padj < 0.1, "*", "")),
    fill_color = case_when(
      padj < 0.05 & direction == "Enriched" ~ "#E74C3C",
      padj < 0.05 & direction == "Depleted" ~ "#3498DB",
      padj < 0.1 & direction == "Enriched" ~ "#F5B7B1",
      padj < 0.1 & direction == "Depleted" ~ "#AED6F1",
      TRUE ~ "#BDC3C7"
    ),
    category = gsub("_", " ", category)
  )

bg_pct <- cat_plot$background_pct[1]

p2 <- ggplot(cat_plot, aes(x = reorder(category, pct_with_dmp), y = pct_with_dmp)) +
  geom_hline(yintercept = bg_pct, linetype = "dashed", color = "#7F8C8D", linewidth = 0.5) +
  geom_col(aes(fill = fill_color), width = 0.7, show.legend = FALSE) +
  scale_fill_identity() +
  geom_text(aes(label = sprintf("%.1f%% %s", pct_with_dmp, sig_label)),
            hjust = -0.1, size = 3.5) +
  coord_flip(ylim = c(0, max(cat_plot$pct_with_dmp) * 1.2)) +
  annotate("text", x = 0.5, y = bg_pct + 0.5,
           label = sprintf("Background: %.1f%%", bg_pct),
           color = "#7F8C8D", size = 3.5, hjust = 0) +
  labs(
    title = "DMP Enrichment by Functional Gene Category",
    subtitle = "Red = enriched | Blue = depleted | *** padj<0.05, * padj<0.1",
    x = NULL,
    y = "Genes with DMPs (%)"
  )

save_pub(p2, "fig02_category_dmp_enrichment", 11, 7)

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 3: Hub Gene Methylation by Module
# ══════════════════════════════════════════════════════════════════════════════

cat("--- Figure 3: Hub gene methylation ---\n")

# Build hub vs non-hub DMP rates
genes_with_dmps <- unique(gene_meth$gene_id)
module_colors <- unique(modules$Module_Color)
module_colors <- module_colors[module_colors != "grey"]

hub_data <- data.frame()
for (mod in module_colors) {
  mod_genes <- modules$Gene[modules$Module_Color == mod]
  mod_hubs <- hubs$Gene[hubs$Module == mod & hubs$IsHub == TRUE]
  mod_nonhubs <- setdiff(mod_genes, mod_hubs)

  if (length(mod_hubs) >= 5) {
    hub_rate <- mean(mod_hubs %in% genes_with_dmps) * 100
    nonhub_rate <- mean(mod_nonhubs %in% genes_with_dmps) * 100

    hub_data <- rbind(hub_data, data.frame(
      module = mod,
      Hub = hub_rate,
      `Non-hub` = nonhub_rate,
      n_hubs = length(mod_hubs),
      check.names = FALSE, stringsAsFactors = FALSE
    ))
  }
}

hub_long <- hub_data %>%
  pivot_longer(cols = c("Hub", "Non-hub"), names_to = "type", values_to = "pct_dmp")

p3 <- ggplot(hub_long, aes(x = reorder(module, -pct_dmp), y = pct_dmp, fill = type)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_manual(values = c("Hub" = "#E74C3C", "Non-hub" = "#3498DB"),
                    name = "Gene Type") +
  labs(
    title = "Hub vs Non-Hub Gene DMP Rates by Module",
    subtitle = "Blue module hubs: 26.3% DMP rate (padj=0.005) | Black hubs: 0% (padj=0.020)",
    x = NULL,
    y = "Genes with DMPs (%)"
  ) +
  theme(legend.position = c(0.85, 0.85),
        legend.background = element_rect(fill = "white", color = "gray80"))

save_pub(p3, "fig03_hub_methylation_by_module", 11, 7)

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 4: DMP Genomic Region Distribution
# ══════════════════════════════════════════════════════════════════════════════

cat("--- Figure 4: DMP genomic regions ---\n")

dmps_cs$simple_region <- case_when(
  grepl("Promoter", dmps_cs$annotation) ~ "Promoter",
  grepl("Intron", dmps_cs$annotation) ~ "Intron",
  grepl("Exon", dmps_cs$annotation) ~ "Exon",
  grepl("UTR", dmps_cs$annotation) ~ "UTR",
  grepl("Downstream", dmps_cs$annotation) ~ "Downstream",
  grepl("Intergenic", dmps_cs$annotation) ~ "Intergenic",
  TRUE ~ "Other"
)

# Merge with methylation diff for direction
dmps_merged <- cbind(dmps_cs, methylation_diff = dmps$methylation_diff)
dmps_merged$direction <- ifelse(dmps_merged$methylation_diff > 0, "Hypermethylated", "Hypomethylated")

region_dir <- dmps_merged %>%
  group_by(simple_region, direction) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(simple_region) %>%
  mutate(pct = n / sum(n) * 100, total = sum(n)) %>%
  ungroup() %>%
  mutate(simple_region = reorder(simple_region, -total))

p4 <- ggplot(region_dir, aes(x = simple_region, y = n, fill = direction)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c("Hypermethylated" = "#E74C3C", "Hypomethylated" = "#3498DB"),
                    name = "Direction") +
  geom_text(data = region_dir %>% group_by(simple_region) %>% summarise(total = sum(n), .groups = "drop"),
            aes(x = simple_region, y = total, label = format(total, big.mark = ",")),
            inherit.aes = FALSE, vjust = -0.3, size = 3.5, fontface = "bold") +
  labs(
    title = "DMP Distribution Across Genomic Regions",
    subtitle = sprintf("18,754 DMPs | Intergenic + Intronic = 63%% of all DMPs"),
    x = NULL,
    y = "Number of DMPs"
  ) +
  theme(legend.position = c(0.85, 0.85),
        legend.background = element_rect(fill = "white", color = "gray80"))

save_pub(p4, "fig04_dmp_genomic_regions", 10, 7)

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 5: Candidate Genes Lollipop (IMPROVED)
# ══════════════════════════════════════════════════════════════════════════════

cat("--- Figure 5: Candidate genes (improved) ---\n")

# Build candidate scoring from scratch with correct DEG path
candidate_data <- modules %>%
  select(gene_id = Gene, module = Module_Color) %>%
  left_join(hubs %>% select(Gene, Module, IsHub, IsCandidate, Connectivity) %>%
              distinct(Gene, .keep_all = TRUE) %>% rename(gene_id = Gene), by = "gene_id") %>%
  left_join(gene_meth %>% select(gene_id, n_dmps, mean_meth_diff, max_abs_diff, primary_region),
            by = "gene_id") %>%
  left_join(deg %>% select(gene_id, log2FoldChange, padj) %>% rename(deg_padj = padj),
            by = "gene_id") %>%
  left_join(annot %>% select(gene_id, gene_symbol), by = "gene_id") %>%
  mutate(
    IsHub = ifelse(is.na(IsHub), FALSE, IsHub),
    IsCandidate = ifelse(is.na(IsCandidate), FALSE, IsCandidate),
    has_dmp = !is.na(n_dmps),
    has_sig_dmp = has_dmp & !is.na(max_abs_diff) & max_abs_diff > 0.15,
    is_de = !is.na(deg_padj) & deg_padj < 0.05,
    large_effect = !is.na(max_abs_diff) & max_abs_diff > 0.25
  )

# Score
candidate_data$score <- 0
candidate_data$score <- candidate_data$score + ifelse(candidate_data$IsHub, 3, 0)
candidate_data$score <- candidate_data$score + ifelse(candidate_data$IsCandidate, 2, 0)
candidate_data$score <- candidate_data$score + ifelse(candidate_data$has_sig_dmp, 3, 0)
candidate_data$score <- candidate_data$score + ifelse(candidate_data$has_dmp & !candidate_data$has_sig_dmp, 1, 0)
candidate_data$score <- candidate_data$score + ifelse(candidate_data$large_effect, 2, 0)
candidate_data$score <- candidate_data$score + ifelse(candidate_data$is_de, 3, 0)
candidate_data$score <- candidate_data$score + ifelse(candidate_data$module == "yellow", 2, 0)
candidate_data$score <- candidate_data$score + ifelse(candidate_data$module == "blue" & candidate_data$IsHub & candidate_data$has_dmp, 2, 0)

candidate_data <- candidate_data %>% arrange(desc(score))

# Top 30 for lollipop
top30 <- head(candidate_data, 30)
top30$label <- ifelse(!is.na(top30$gene_symbol) & top30$gene_symbol != "",
                      paste0(top30$gene_symbol, " (", top30$gene_id, ")"),
                      top30$gene_id)

# Module color mapping for a nice palette
module_pal <- c(
  "yellow" = "#F1C40F", "blue" = "#2980B9", "red" = "#E74C3C",
  "green" = "#27AE60", "brown" = "#8B4513", "turquoise" = "#1ABC9C",
  "black" = "#2C3E50", "pink" = "#FF69B4", "magenta" = "#C0392B",
  "purple" = "#8E44AD", "tan" = "#D2B48C", "greenyellow" = "#ADFF2F", "grey" = "#95A5A6"
)

p5 <- ggplot(top30, aes(x = reorder(label, score), y = score)) +
  geom_segment(aes(xend = reorder(label, score), y = 0, yend = score),
               linewidth = 1.2, color = "gray70") +
  geom_point(aes(fill = module), shape = 21, size = 5, stroke = 0.5, color = "white") +
  scale_fill_manual(values = module_pal, name = "Module") +
  coord_flip() +
  labs(
    title = "Top 30 Candidate Genes for Experimental Validation",
    subtitle = "Score: Hub(+3) + DE(+3) + Sig DMP(+3) + Large effect(+2) + Candidate(+2) + Yellow(+2) + Blue hub+DMP(+2) + DMP(+1)",
    x = NULL,
    y = "Composite Evidence Score"
  ) +
  theme(
    axis.text.y = element_text(size = 8),
    legend.position = "right"
  )

save_pub(p5, "fig05_top30_candidates_lollipop", 14, 10)

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 6: Evidence Matrix Heatmap (FIXED)
# ══════════════════════════════════════════════════════════════════════════════

cat("--- Figure 6: Evidence matrix (fixed) ---\n")

top40 <- head(candidate_data, 40)
top40$label <- ifelse(!is.na(top40$gene_symbol) & top40$gene_symbol != "",
                      paste0(top40$gene_symbol, " (", top40$gene_id, ")"),
                      top40$gene_id)

feature_mat <- data.frame(
  Hub = as.integer(top40$IsHub),
  `Sig DMP` = as.integer(top40$has_sig_dmp),
  `Large effect` = as.integer(top40$large_effect),
  DE = as.integer(top40$is_de),
  `Yellow mod` = as.integer(top40$module == "yellow"),
  `Blue hub+DMP` = as.integer(top40$module == "blue" & top40$IsHub & top40$has_dmp),
  Candidate = as.integer(top40$IsCandidate),
  check.names = FALSE
)
rownames(feature_mat) <- top40$label
feature_mat <- as.matrix(feature_mat)

ann_row <- data.frame(Module = top40$module, Score = top40$score, row.names = top40$label)

# Build clean module color annotation
unique_mods <- unique(top40$module)
ann_colors <- list(
  Module = module_pal[unique_mods],
  Score = colorRampPalette(c("#AED6F1", "#2C3E50"))(5)
)
names(ann_colors$Score) <- seq(min(top40$score), max(top40$score), length.out = 5)

png(file.path(outdir, "fig06_evidence_matrix.png"), width = 11, height = 14, units = "in", res = 300)
pheatmap(feature_mat,
         cluster_rows = FALSE, cluster_cols = FALSE,
         color = c("#ECF0F1", "#E74C3C"),
         legend = FALSE,
         annotation_row = ann_row,
         annotation_colors = ann_colors,
         main = "Top 40 Candidate Genes: Evidence Matrix",
         fontsize = 10, fontsize_row = 8,
         border_color = "gray85",
         cellwidth = 28, cellheight = 14,
         gaps_col = c(1, 3, 4, 5, 6))
dev.off()

pdf(file.path(outdir, "fig06_evidence_matrix.pdf"), width = 11, height = 14)
pheatmap(feature_mat,
         cluster_rows = FALSE, cluster_cols = FALSE,
         color = c("#ECF0F1", "#E74C3C"),
         legend = FALSE,
         annotation_row = ann_row,
         annotation_colors = ann_colors,
         main = "Top 40 Candidate Genes: Evidence Matrix",
         fontsize = 10, fontsize_row = 8,
         border_color = "gray85",
         cellwidth = 28, cellheight = 14,
         gaps_col = c(1, 3, 4, 5, 6))
dev.off()
cat("  Saved: fig06_evidence_matrix\n")

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 7: Four-Layer Model Diagram
# ══════════════════════════════════════════════════════════════════════════════

cat("--- Figure 7: Four-layer model ---\n")

model_data <- data.frame(
  layer = factor(rep(c("Layer 1:\nProtected Regulators",
                        "Layer 2:\nTargeted Modules",
                        "Layer 3:\nDistal Regulation",
                        "Layer 4:\nEffector Targets"), c(3, 4, 3, 3)),
                 levels = c("Layer 1:\nProtected Regulators",
                            "Layer 2:\nTargeted Modules",
                            "Layer 3:\nDistal Regulation",
                            "Layer 4:\nEffector Targets")),
  feature = c("Hox genes\n(0/27 DMPs)", "Homeodomain\nTFs (1.7%)", "Core TFs\n(OR=0.77)",
              "Yellow\n(cell cycle)", "Blue\n(dual DMP+DMR)", "Brown\n(hypomethylated)", "Black hubs\n(0% DMPs)",
              "Intergenic\n(31%)", "Intronic\n(32%)", "DMP clusters\n(2,109)",
              "Zinc fingers\n(22.4%)", "Cytoskeleton\n(OR=1.45)", "RNA processing\n(depleted)"),
  value = c(0, 1.7, 7.7, 19.5, 17.6, 14.8, 0, 31, 32, 14, 22.4, 21.2, 12.0),
  type = factor(c(rep("Protected", 3), rep("Targeted", 4),
                   rep("Distal", 3), rep("Effector", 3)),
                levels = c("Protected", "Targeted", "Distal", "Effector")),
  stringsAsFactors = FALSE
)

p7 <- ggplot(model_data, aes(x = feature, y = value, fill = type)) +
  geom_col(width = 0.7, alpha = 0.9) +
  facet_wrap(~layer, scales = "free_x", nrow = 1) +
  scale_fill_manual(
    values = c("Protected" = "#27AE60", "Targeted" = "#E74C3C",
               "Distal" = "#2980B9", "Effector" = "#F39C12"),
    name = "Role"
  ) +
  geom_text(aes(label = ifelse(value > 0, paste0(value, "%"), "0%")),
            vjust = -0.3, size = 3.5, fontface = "bold") +
  labs(
    title = "Four-Layer Model of Epigenetic Regulation During Regeneration",
    subtitle = "D. laeve protects core regulators while targeting functional modules through distal methylation",
    x = NULL,
    y = "DMP rate or percentage"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    legend.position = "bottom",
    strip.background = element_rect(fill = "gray95", color = "gray80")
  )

save_pub(p7, "fig07_four_layer_model", 16, 8)

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 8: DMP-DMR Concordance
# ══════════════════════════════════════════════════════════════════════════════

cat("--- Figure 8: DMP-DMR concordance ---\n")

dmr_dir <- dmrs %>%
  group_by(nearest_gene) %>%
  summarise(mean_dmr_meth = mean(diff_methy), .groups = "drop")

concordance <- gene_meth %>%
  inner_join(dmr_dir, by = c("gene_id" = "nearest_gene")) %>%
  mutate(
    concordant = (mean_meth_diff > 0 & mean_dmr_meth > 0) | (mean_meth_diff < 0 & mean_dmr_meth < 0)
  )

p8 <- ggplot(concordance, aes(x = mean_meth_diff, y = mean_dmr_meth)) +
  geom_hline(yintercept = 0, color = "gray80", linewidth = 0.4) +
  geom_vline(xintercept = 0, color = "gray80", linewidth = 0.4) +
  geom_point(aes(color = concordant), alpha = 0.4, size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_manual(
    values = c("TRUE" = "#27AE60", "FALSE" = "#E74C3C"),
    labels = c("Discordant", "Concordant"),
    name = "Direction"
  ) +
  annotate("text", x = -0.6, y = 0.7,
           label = sprintf("Concordance: %.1f%%\np = 1.4 x 10^-160",
                           mean(concordance$concordant) * 100),
           size = 4.5, fontface = "bold", color = "#2C3E50") +
  labs(
    title = "DMP vs DMR Methylation Direction Per Gene",
    subtitle = sprintf("%d genes with both DMPs and DMRs | Near-perfect directional agreement", nrow(concordance)),
    x = "Mean DMP methylation difference",
    y = "Mean DMR methylation difference"
  )

save_pub(p8, "fig08_dmp_dmr_concordance", 9, 8)

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 9: Methylation Threshold (no dose-response)
# ══════════════════════════════════════════════════════════════════════════════

cat("--- Figure 9: Methylation threshold ---\n")

gene_meth$abs_meth <- abs(gene_meth$mean_meth_diff)
gene_meth$meth_bin <- cut(gene_meth$abs_meth,
                          breaks = c(0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 1),
                          labels = c("<10%", "10-15%", "15-20%", "20-25%",
                                     "25-30%", "30-40%", "40-50%", ">50%"))

threshold_data <- gene_meth %>%
  filter(!is.na(meth_bin)) %>%
  left_join(deg %>% select(gene_id, padj) %>% rename(deg_padj = padj), by = "gene_id") %>%
  group_by(meth_bin) %>%
  summarise(
    n = n(),
    n_de = sum(!is.na(deg_padj) & deg_padj < 0.05),
    pct_de = n_de / n * 100,
    .groups = "drop"
  )

p9 <- ggplot(threshold_data, aes(x = meth_bin, y = pct_de)) +
  geom_col(fill = "#3498DB", alpha = 0.8, width = 0.7) +
  geom_text(aes(label = sprintf("n=%d", n)), vjust = -0.3, size = 3.2, color = "gray40") +
  geom_hline(yintercept = mean(threshold_data$pct_de), linetype = "dashed",
             color = "#E74C3C", linewidth = 0.5) +
  annotate("text", x = 7.5, y = max(threshold_data$pct_de) * 0.9,
           label = "Cochran-Armitage\ntrend p = 0.178\n(NOT significant)",
           size = 4, fontface = "bold.italic", color = "#7F8C8D") +
  labs(
    title = "Methylation Does NOT Act as a Rheostat",
    subtitle = "DE rate does not increase with methylation effect size | Methylation is a switch, not a dial",
    x = "Absolute methylation difference",
    y = "Differentially expressed (%)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_pub(p9, "fig09_methylation_threshold", 10, 7)

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 10: Module Evidence Summary Heatmap
# ══════════════════════════════════════════════════════════════════════════════

cat("--- Figure 10: Module evidence heatmap ---\n")

# Read DMR module data
dmr_mod <- read.delim("results/08_dmr_spatial/tables/J01_dmr_module_enrichment.tsv",
                       stringsAsFactors = FALSE)

evidence <- master %>%
  select(module, pct_genes_with_dmp, fisher_padj, Cohens_d, Direction) %>%
  left_join(dmr_mod %>% select(module, pct_with_dmr, odds_ratio) %>%
              rename(dmr_pct = pct_with_dmr, dmr_or = odds_ratio), by = "module")

heat_mat <- evidence %>%
  mutate(
    `DMP %` = pct_genes_with_dmp,
    `DMP -log10(p)` = -log10(pmax(fisher_padj, 1e-10)),
    `DMR OR` = ifelse(is.na(dmr_or), 1, dmr_or),
    `Expression |d|` = abs(Cohens_d)
  ) %>%
  select(module, `DMP %`, `DMP -log10(p)`, `DMR OR`, `Expression |d|`)

mat <- as.matrix(heat_mat[, -1])
rownames(mat) <- heat_mat$module
mat_scaled <- scale(mat)

png(file.path(outdir, "fig10_module_evidence_heatmap.png"), width = 10, height = 7, units = "in", res = 300)
pheatmap(mat_scaled,
         color = colorRampPalette(c("#2980B9", "#ECF0F1", "#E74C3C"))(100),
         main = "Module-Level Evidence Matrix (Z-scored)",
         fontsize = 12,
         fontsize_row = 11,
         fontsize_col = 11,
         cluster_cols = FALSE,
         border_color = "gray90",
         display_numbers = TRUE,
         number_format = "%.1f",
         number_color = "gray20",
         cellwidth = 50, cellheight = 28)
dev.off()

pdf(file.path(outdir, "fig10_module_evidence_heatmap.pdf"), width = 10, height = 7)
pheatmap(mat_scaled,
         color = colorRampPalette(c("#2980B9", "#ECF0F1", "#E74C3C"))(100),
         main = "Module-Level Evidence Matrix (Z-scored)",
         fontsize = 12,
         fontsize_row = 11,
         fontsize_col = 11,
         cluster_cols = FALSE,
         border_color = "gray90",
         display_numbers = TRUE,
         number_format = "%.1f",
         number_color = "gray20",
         cellwidth = 50, cellheight = 28)
dev.off()
cat("  Saved: fig10_module_evidence_heatmap\n")

# ══════════════════════════════════════════════════════════════════════════════

cat("\n=== Script 11 complete ===\n")
cat("Output:", outdir, "\n")
cat("Figures:", length(list.files(outdir, pattern = "\\.png$")), "PNG +",
    length(list.files(outdir, pattern = "\\.pdf$")), "PDF\n")
