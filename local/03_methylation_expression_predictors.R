#!/usr/bin/env Rscript
# =============================================================================
# Which genomic regions best predict expression from methylation?
# =============================================================================
# Analyzes module × region Spearman correlations between DMP methylation
# and gene expression log2FC to identify the strongest predictive relationships.
#
# Key question: Where does methylation actually matter for expression?
#
# Inputs:
#   - results/03_integration/Tables/MXT_module_region_correlation.txt
#   - results/03_integration/Tables/MXT_correlation_by_region.txt
#   - results/03_integration/Tables/MXT_intergenic_distance_vs_expression.txt
#
# Run: Rscript local/03_methylation_expression_predictors.R
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(pheatmap)
  library(RColorBrewer)
  library(scales)
})

options(stringsAsFactors = FALSE)
options(scipen = 999)

theme_set(theme_minimal() +
            theme(text = element_text(size = 12),
                  plot.title = element_text(size = 14, face = "bold"),
                  axis.title = element_text(size = 12),
                  legend.text = element_text(size = 10)))

# ── Paths ────────────────────────────────────────────────────────────────────

INT_DIR    <- "results/03_integration/Tables"
OUTPUT_DIR <- "results/10_meth_expr_predictors"

dir.create(file.path(OUTPUT_DIR, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "plots"),  recursive = TRUE, showWarnings = FALSE)

save_both <- function(plot_obj, name, width, height) {
  ggsave(file.path(OUTPUT_DIR, "plots", paste0(name, ".pdf")),
         plot_obj, width = width, height = height)
  ggsave(file.path(OUTPUT_DIR, "plots", paste0(name, ".png")),
         plot_obj, width = width, height = height, dpi = 300)
  cat(sprintf("  Saved: %s (.pdf + .png)\n", name))
}

save_pheatmap_both <- function(pheatmap_expr, name, width, height) {
  pdf(file.path(OUTPUT_DIR, "plots", paste0(name, ".pdf")),
      width = width, height = height)
  eval(pheatmap_expr)
  dev.off()
  png(file.path(OUTPUT_DIR, "plots", paste0(name, ".png")),
      width = width, height = height, units = "in", res = 300)
  eval(pheatmap_expr)
  dev.off()
  cat(sprintf("  Saved: %s (.pdf + .png)\n", name))
}

cat("================================================================\n")
cat("  Methylation-Expression Predictors by Region & Module\n")
cat("================================================================\n\n")


# ── Load data ────────────────────────────────────────────────────────────────

mod_region <- read.delim(file.path(INT_DIR, "MXT_module_region_correlation.txt"))
global_region <- read.delim(file.path(INT_DIR, "MXT_correlation_by_region.txt"))

cat(sprintf("  Module-region correlations: %d entries\n", nrow(mod_region)))
cat(sprintf("  Global region correlations: %d regions\n", nrow(global_region)))

# Mark significance
mod_region$sig <- mod_region$spearman_p < 0.05
mod_region$sig_label <- ifelse(mod_region$spearman_p < 0.001, "***",
                        ifelse(mod_region$spearman_p < 0.01, "**",
                        ifelse(mod_region$spearman_p < 0.05, "*", "")))
mod_region$neg_log10_p <- -log10(mod_region$spearman_p)
mod_region$direction <- ifelse(mod_region$spearman_r > 0, "Positive", "Negative")


# #############################################################################
# FIGURE D1: Heatmap — Module x Region Spearman r
# #############################################################################

cat("\n=== Generating methylation-expression predictor plots ===\n\n")

# Pivot to matrix
cor_wide <- mod_region %>%
  select(module, region, spearman_r) %>%
  pivot_wider(names_from = region, values_from = spearman_r)

cor_mat <- as.matrix(cor_wide[, -1])
rownames(cor_mat) <- cor_wide$module
cor_mat[is.na(cor_mat)] <- 0

# Significance stars for annotation
sig_wide <- mod_region %>%
  select(module, region, sig_label) %>%
  pivot_wider(names_from = region, values_from = sig_label, values_fill = "")

sig_mat <- as.matrix(sig_wide[, -1])
rownames(sig_mat) <- sig_wide$module
sig_mat[is.na(sig_mat)] <- ""

# Reorder columns by region importance
region_order <- c("Promoter", "Intron", "Intergenic", "UTR", "Exon", "Downstream")
region_order <- intersect(region_order, colnames(cor_mat))
cor_mat <- cor_mat[, region_order]
sig_mat <- sig_mat[, region_order]

hm_expr <- quote(
  pheatmap(cor_mat,
           cluster_rows = TRUE, cluster_cols = FALSE,
           display_numbers = sig_mat,
           fontsize_number = 14,
           number_color = "black",
           color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
           breaks = seq(-0.85, 0.85, length.out = 101),
           main = "Methylation-Expression Correlation: Module x Region (Spearman r)",
           fontsize = 11,
           angle_col = 0)
)
save_pheatmap_both(hm_expr, "fig_D1_meth_expr_correlation_heatmap", 10, 8)

# #############################################################################
# FIGURE D2: Dot plot — Significant correlations ranked by effect size
# #############################################################################

sig_only <- mod_region %>%
  filter(sig) %>%
  mutate(
    combo = paste(module, region, sep = " : "),
    abs_r = abs(spearman_r)
  ) %>%
  arrange(desc(abs_r))

p_dotplot <- ggplot(sig_only,
                    aes(x = spearman_r,
                        y = reorder(combo, abs_r),
                        size = n_obs,
                        color = region)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
  geom_point(alpha = 0.8) +
  scale_color_brewer(palette = "Set2", name = "Genomic Region") +
  scale_size_continuous(name = "N genes", range = c(2, 8)) +
  labs(
    title = "Significant Methylation-Expression Correlations (p<0.05)",
    subtitle = "D. laeve: Spearman r between DMP methylation and log2FC, by module and region",
    x = "Spearman r (methylation vs expression)",
    y = NULL
  ) +
  theme(axis.text.y = element_text(size = 9))
save_both(p_dotplot, "fig_D2_significant_meth_expr_correlations", 13, 10)

# #############################################################################
# FIGURE D3: Faceted by region — which modules show strongest effect?
# #############################################################################

# Only regions with >= 3 module-level entries
region_counts <- mod_region %>% count(region) %>% filter(n >= 3)

p_facet <- ggplot(mod_region %>% filter(region %in% region_counts$region),
                  aes(x = reorder(module, spearman_r),
                      y = spearman_r,
                      fill = module)) +
  geom_col(alpha = 0.85) +
  geom_hline(yintercept = 0, color = "gray40") +
  geom_text(aes(label = sig_label), vjust = ifelse(mod_region$spearman_r[mod_region$region %in% region_counts$region] >= 0, -0.3, 1.3), size = 4) +
  scale_fill_identity() +
  facet_wrap(~region, scales = "free_x", ncol = 3) +
  coord_flip() +
  labs(
    title = "Methylation-Expression Correlation by Region Across Modules",
    subtitle = "D. laeve: Spearman r. * p<0.05, ** p<0.01, *** p<0.001",
    x = NULL, y = "Spearman r"
  ) +
  theme(strip.text = element_text(face = "bold", size = 11))
save_both(p_facet, "fig_D3_meth_expr_by_region_faceted", 16, 12)

# #############################################################################
# FIGURE D4: Global vs module-specific — why global washes out
# #############################################################################

global_long <- global_region %>%
  select(region, spearman_r, spearman_p) %>%
  mutate(module = "GLOBAL", n_obs = global_region$n_genes, sig = spearman_p < 0.05)

combined <- bind_rows(
  mod_region %>%
    filter(sig) %>%
    select(module, region, spearman_r, spearman_p, n_obs, sig),
  global_long
) %>%
  mutate(is_global = module == "GLOBAL")

p_global <- ggplot(combined %>% filter(region %in% c("Promoter", "Intron", "Intergenic", "UTR")),
                   aes(x = module, y = spearman_r,
                       fill = is_global, alpha = is_global)) +
  geom_col() +
  geom_hline(yintercept = 0, color = "gray40") +
  scale_fill_manual(values = c("FALSE" = "#4DBEEE", "TRUE" = "#E74C3C"), guide = "none") +
  scale_alpha_manual(values = c("FALSE" = 0.7, "TRUE" = 1), guide = "none") +
  facet_wrap(~region, scales = "free_x", ncol = 2) +
  coord_flip() +
  labs(
    title = "Why Global Correlations Wash Out: Opposing Module Effects",
    subtitle = "D. laeve: Red = global (non-significant), blue = significant module-level correlations",
    x = NULL, y = "Spearman r"
  ) +
  theme(strip.text = element_text(face = "bold", size = 11))
save_both(p_global, "fig_D4_global_vs_module_correlations", 14, 10)

# #############################################################################
# FIGURE D5: Top predictive module-region combos — bubble plot
# #############################################################################

top_predictors <- sig_only %>%
  head(15)

p_bubble <- ggplot(top_predictors,
                   aes(x = spearman_r, y = neg_log10_p,
                       size = n_obs, color = region, shape = direction)) +
  geom_point(alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "red", alpha = 0.5) +
  ggrepel::geom_text_repel(
    aes(label = paste(module, region, sep = "\n")),
    size = 3, max.overlaps = 20, box.padding = 0.5
  ) +
  scale_color_brewer(palette = "Set2", name = "Region") +
  scale_size_continuous(name = "N genes", range = c(3, 10)) +
  scale_shape_manual(values = c("Positive" = 17, "Negative" = 16), name = "Direction") +
  labs(
    title = "Top 15 Methylation-Expression Predictors by Module & Region",
    subtitle = "D. laeve: Strongest Spearman correlations between DMP methylation and log2FC",
    x = "Spearman r", y = expression(-log[10](p))
  )
save_both(p_bubble, "fig_D5_top_predictors_bubble", 14, 10)


# #############################################################################
# SUMMARY TABLE
# #############################################################################

cat("\n=== Summary: Best Methylation-Expression Predictors ===\n\n")

summary_table <- sig_only %>%
  select(module, region, n_obs, spearman_r, spearman_p, direction) %>%
  mutate(
    spearman_r = round(spearman_r, 3),
    spearman_p = signif(spearman_p, 3),
    interpretation = case_when(
      direction == "Negative" & region == "Promoter" ~ "Classical silencing",
      direction == "Positive" & region == "Promoter" ~ "Paradoxical activation",
      direction == "Negative" & region == "Intron" ~ "Intronic repression",
      direction == "Positive" & region == "Intron" ~ "Intronic activation",
      direction == "Negative" & region == "Intergenic" ~ "Enhancer/silencer repression",
      direction == "Positive" & region == "Intergenic" ~ "Enhancer/silencer activation",
      direction == "Negative" & region %in% c("UTR", "Downstream") ~ "Distal repression",
      direction == "Positive" & region %in% c("UTR", "Exon") ~ "Genic activation",
      TRUE ~ paste(direction, "effect")
    )
  )

cat("  All significant module-region correlations (p<0.05):\n\n")
print(as.data.frame(summary_table), row.names = FALSE)

write.table(summary_table,
            file.path(OUTPUT_DIR, "tables", "meth_expr_predictors_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Region summary
cat("\n\n  Region-level summary of significant correlations:\n")
region_summary <- sig_only %>%
  group_by(region) %>%
  summarise(
    n_significant = n(),
    mean_abs_r = round(mean(abs(spearman_r)), 3),
    max_abs_r = round(max(abs(spearman_r)), 3),
    n_positive = sum(direction == "Positive"),
    n_negative = sum(direction == "Negative"),
    modules = paste(module, collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_abs_r))

print(as.data.frame(region_summary))

write.table(region_summary,
            file.path(OUTPUT_DIR, "tables", "region_predictive_power_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)


# #############################################################################
# FINAL REPORT
# #############################################################################

cat("\n================================================================\n")
cat("  METHYLATION-EXPRESSION PREDICTOR ANALYSIS COMPLETE\n")
cat("================================================================\n\n")

cat("KEY FINDINGS:\n\n")
cat("1. GLOBAL correlations are non-significant in ALL regions.\n")
cat("   This is because different modules have OPPOSITE effects that cancel out.\n\n")
cat("2. Top predictive combinations (by |Spearman r|):\n")
for (i in seq_len(min(10, nrow(summary_table)))) {
  row <- summary_table[i, ]
  cat(sprintf("   %2d. %s + %s: r = %+.3f (p = %s) -- %s\n",
              i, row$module, row$region, row$spearman_r,
              formatC(row$spearman_p, format = "e", digits = 1),
              row$interpretation))
}
cat("\n3. The relationship between methylation and expression is\n")
cat("   ENTIRELY module-specific and region-specific in D. laeve.\n")
cat("================================================================\n")
