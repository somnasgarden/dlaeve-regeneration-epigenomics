#!/usr/bin/env Rscript
# =============================================================================
# Threshold Sensitivity Analysis: Module & Gene Enrichment
# =============================================================================
# How do enrichment results change when we vary:
#   Part A — Enrichment padj cutoff (which GO/KEGG/Reactome terms are significant?)
#   Part B — DMP stringency (FDR + |meth_diff|) on module DMP enrichment (Fisher's)
#   Part C — Combined view
#
# Inputs (all from results/):
#   - results/02_rnaseq/Part5_WGCNA/data/*_enrichment_*.tsv  (60 files)
#   - results/02_rnaseq/Part5_WGCNA/data/all_gene_module_assignments.tsv
#   - results/01_methylation/dmps_annotated.txt  (18,754 DMPs with FDR + diff)
#   - results/01_methylation/dmps_chipseeker_annotated.txt  (ChIPseeker regions)
#
# Run locally: Rscript local/02_threshold_sensitivity_analysis.R
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

WGCNA_DIR  <- "results/02_rnaseq/Part5_WGCNA/data"
METH_DIR   <- "results/01_methylation"
OUTPUT_DIR <- "results/06_threshold_sensitivity"

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
cat("  Threshold Sensitivity Analysis\n")
cat("================================================================\n\n")


# #############################################################################
# PART A: ENRICHMENT padj THRESHOLD SENSITIVITY
# #############################################################################

cat("=== Part A: Enrichment padj Sensitivity ===\n\n")

# Load all enrichment TSVs
categories <- c("BP", "MF", "CC", "KEGG", "Reactome")
modules <- c("yellow", "turquoise", "red", "brown", "blue", "green",
             "black", "pink", "purple", "magenta", "greenyellow", "tan")

padj_thresholds <- c(0.001, 0.005, 0.01, 0.05, 0.1, 0.2)

# Build master enrichment table
all_enrichment <- list()
for (cat_name in categories) {
  for (mod in modules) {
    fname <- file.path(WGCNA_DIR, paste0(cat_name, "_enrichment_", mod, ".tsv"))
    if (!file.exists(fname)) next
    df <- read.delim(fname)
    if (nrow(df) == 0) next
    df$module   <- mod
    df$category <- cat_name
    all_enrichment[[paste(cat_name, mod)]] <- df
  }
}
enrich_all <- do.call(rbind, all_enrichment)
cat(sprintf("  Loaded %s enrichment terms across %d modules × %d categories\n",
            format(nrow(enrich_all), big.mark = ","),
            length(unique(enrich_all$module)),
            length(unique(enrich_all$category))))

# A.1: Count significant terms per module at each padj threshold
padj_sensitivity <- expand.grid(module = modules, threshold = padj_thresholds,
                                stringsAsFactors = FALSE)
padj_sensitivity$n_sig_terms     <- 0
padj_sensitivity$n_sig_BP        <- 0
padj_sensitivity$n_sig_MF        <- 0
padj_sensitivity$n_sig_CC        <- 0
padj_sensitivity$n_sig_KEGG      <- 0
padj_sensitivity$n_sig_Reactome  <- 0

for (i in seq_len(nrow(padj_sensitivity))) {
  mod   <- padj_sensitivity$module[i]
  thr   <- padj_sensitivity$threshold[i]
  sub   <- enrich_all[enrich_all$module == mod & enrich_all$padj < thr, ]
  padj_sensitivity$n_sig_terms[i] <- nrow(sub)
  for (cat_name in categories) {
    col_name <- paste0("n_sig_", cat_name)
    padj_sensitivity[[col_name]][i] <- sum(sub$category == cat_name)
  }
}

write.table(padj_sensitivity,
            file.path(OUTPUT_DIR, "tables", "enrichment_padj_sensitivity.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("  Enrichment padj sensitivity computed.\n\n")

# A.2: Heatmap — modules × padj thresholds → total significant terms
hm_data <- padj_sensitivity %>%
  select(module, threshold, n_sig_terms) %>%
  pivot_wider(names_from = threshold, values_from = n_sig_terms)

hm_mat <- as.matrix(hm_data[, -1])
rownames(hm_mat) <- hm_data$module
colnames(hm_mat) <- paste0("padj<", colnames(hm_mat))

# Order rows by total terms at padj < 0.05
order_idx <- order(-hm_mat[, which(colnames(hm_mat) == "padj<0.05")])
hm_mat <- hm_mat[order_idx, ]

hm_expr <- quote(
  pheatmap(hm_mat,
           cluster_rows = FALSE, cluster_cols = FALSE,
           display_numbers = TRUE, number_format = "%d",
           color = colorRampPalette(c("white", "#FEE08B", "#E74C3C"))(50),
           main = "Significant Enrichment Terms per Module at Different padj Thresholds",
           fontsize = 11, fontsize_number = 9,
           angle_col = 0)
)
save_pheatmap_both(hm_expr, "fig_A1_enrichment_padj_heatmap", 10, 8)

# A.3: Line plot — how terms accumulate with relaxed threshold
padj_agg <- padj_sensitivity %>%
  group_by(module, threshold) %>%
  summarise(n_terms = sum(n_sig_terms), .groups = "drop")

# Only plot modules that have any enrichment
active_modules <- padj_agg %>% filter(n_terms > 0) %>% pull(module) %>% unique()

p_padj_line <- ggplot(padj_agg %>% filter(module %in% active_modules),
                      aes(x = threshold, y = n_terms, color = module, group = module)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  scale_color_identity() +
  scale_x_continuous(breaks = padj_thresholds, trans = "log10",
                     labels = padj_thresholds) +
  labs(
    title = "Enrichment Terms Gained at Relaxed padj Thresholds — D. laeve",
    subtitle = "Modules colored by WGCNA color assignment",
    x = "padj Threshold (log scale)", y = "Number of Significant Terms"
  )
save_both(p_padj_line, "fig_A2_enrichment_padj_lineplot", 12, 8)

# A.4: Per-category breakdown at padj < 0.05 vs padj < 0.1
compare_thresholds <- padj_sensitivity %>%
  filter(threshold %in% c(0.05, 0.1)) %>%
  select(module, threshold, n_sig_BP, n_sig_MF, n_sig_CC, n_sig_KEGG, n_sig_Reactome) %>%
  pivot_longer(cols = starts_with("n_sig_"),
               names_to = "category", values_to = "n_terms") %>%
  mutate(category = gsub("n_sig_", "", category),
         threshold_label = paste0("padj < ", threshold))

p_cat_compare <- ggplot(compare_thresholds %>% filter(n_terms > 0),
                        aes(x = module, y = n_terms, fill = threshold_label)) +
  geom_col(position = "dodge", alpha = 0.85) +
  facet_wrap(~category, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c("padj < 0.05" = "#2166AC", "padj < 0.1" = "#92C5DE"),
                    name = "Threshold") +
  coord_flip() +
  labs(
    title = "Enrichment Terms: padj < 0.05 vs padj < 0.1 — D. laeve",
    subtitle = "Terms gained by relaxing threshold from 0.05 to 0.1",
    x = NULL, y = "Number of Significant Terms"
  ) +
  theme(axis.text.y = element_text(size = 9))
save_both(p_cat_compare, "fig_A3_enrichment_005_vs_01_by_category", 16, 10)


# #############################################################################
# PART B: DMP STRINGENCY SENSITIVITY ON MODULE ENRICHMENT
# #############################################################################

cat("\n=== Part B: DMP Threshold Sensitivity on Module Enrichment ===\n\n")

# B.1: Load DMP data + module assignments
dmps <- read.delim(file.path(METH_DIR, "dmps_annotated.txt"))
cat(sprintf("  DMPs loaded: %s\n", format(nrow(dmps), big.mark = ",")))

# ChIPseeker annotations for genomic regions
cs_file <- file.path(METH_DIR, "dmps_chipseeker_annotated.txt")
if (file.exists(cs_file)) {
  dmp_cs <- read.delim(cs_file)
  cat(sprintf("  ChIPseeker annotations: %s\n", format(nrow(dmp_cs), big.mark = ",")))
  # Use ChIPseeker gene assignment if available
  if (nrow(dmp_cs) == nrow(dmps)) {
    dmps$cs_gene_id   <- dmp_cs$geneId
    dmps$cs_annotation <- dmp_cs$annotation
  }
}

# Module assignments
mod_assign <- read.delim(file.path(WGCNA_DIR, "all_gene_module_assignments.tsv"))
cat(sprintf("  Module assignments: %s genes\n", format(nrow(mod_assign), big.mark = ",")))

# Use best available gene mapping
dmps$gene_id <- if ("cs_gene_id" %in% colnames(dmps)) dmps$cs_gene_id else dmps$nearest_gene

# Simplify annotations to broad regions
if ("cs_annotation" %in% colnames(dmps)) {
  dmps$region <- dmps$cs_annotation
} else {
  dmps$region <- dmps$annotation
}
dmps$region <- gsub(" \\(.*", "", dmps$region)
dmps$region[grepl("Promoter", dmps$region)] <- "Promoter"
dmps$region[grepl("Exon", dmps$region)]     <- "Exon"
dmps$region[grepl("Intron", dmps$region)]   <- "Intron"
dmps$region[grepl("UTR", dmps$region)]      <- "UTR"
dmps$region[grepl("Downstream", dmps$region)] <- "Downstream"
dmps$region[grepl("Intergenic|Distal", dmps$region)] <- "Intergenic"

# Clean unmapped
dmps <- dmps[!is.na(dmps$gene_id) & dmps$gene_id != "", ]
dmps$abs_diff <- abs(dmps$methylation_diff)

cat(sprintf("  DMPs with gene assignment: %s\n", format(nrow(dmps), big.mark = ",")))

# B.2: Define threshold grid
fdr_thresholds  <- c(1, 0.1, 0.05, 0.01, 0.001)
diff_thresholds <- c(0.1, 0.15, 0.2, 0.25, 0.3)

real_modules <- setdiff(unique(mod_assign$Module_Color), "grey")
module_sizes <- mod_assign %>%
  filter(Module_Color != "grey") %>%
  count(Module_Color, name = "module_size")

total_genes_network <- sum(module_sizes$module_size)

# B.3: Run Fisher's exact test for each threshold combination
cat("\n  Running Fisher's tests across threshold grid...\n")

fisher_grid <- list()
dmp_count_grid <- list()

for (fdr_thr in fdr_thresholds) {
  for (diff_thr in diff_thresholds) {
    # Filter DMPs
    dmps_filt <- dmps[dmps$fdr <= fdr_thr & dmps$abs_diff >= diff_thr, ]
    n_dmps_filt <- nrow(dmps_filt)

    # Map to genes in network
    dmp_gene_module <- merge(
      dmps_filt[, c("gene_id", "methylation_diff", "region")],
      mod_assign[, c("Gene", "Module_Color")],
      by.x = "gene_id", by.y = "Gene"
    )
    dmp_gene_module <- dmp_gene_module[dmp_gene_module$Module_Color != "grey", ]

    # Unique DMP genes in network
    dmp_genes_all <- unique(dmp_gene_module$gene_id)
    total_dmp_genes <- length(dmp_genes_all)

    # Store DMP counts
    dmp_count_grid[[paste(fdr_thr, diff_thr)]] <- data.frame(
      fdr_threshold  = fdr_thr,
      diff_threshold = diff_thr,
      n_dmps         = n_dmps_filt,
      n_genes_with_dmp = total_dmp_genes,
      stringsAsFactors = FALSE
    )

    # Fisher's test per module
    for (mod in real_modules) {
      mod_size <- module_sizes$module_size[module_sizes$Module_Color == mod]
      mod_dmp_genes <- unique(dmp_gene_module$gene_id[dmp_gene_module$Module_Color == mod])
      mod_dmp <- length(mod_dmp_genes)
      n_dmps_mod <- sum(dmp_gene_module$Module_Color == mod)

      # 2x2: DMP_gene/not × in_module/not
      a     <- mod_dmp
      b     <- mod_size - mod_dmp
      c_val <- total_dmp_genes - mod_dmp
      d     <- total_genes_network - mod_size - c_val

      ft <- tryCatch(
        fisher.test(matrix(c(a, c_val, b, d), nrow = 2), alternative = "greater"),
        error = function(e) list(p.value = NA, estimate = NA)
      )

      fisher_grid[[paste(fdr_thr, diff_thr, mod)]] <- data.frame(
        fdr_threshold     = fdr_thr,
        diff_threshold    = diff_thr,
        module            = mod,
        module_size       = mod_size,
        n_dmps_in_module  = n_dmps_mod,
        n_genes_with_dmp  = mod_dmp,
        pct_genes_with_dmp = 100 * mod_dmp / mod_size,
        fisher_p          = ft$p.value,
        odds_ratio        = as.numeric(ft$estimate),
        stringsAsFactors  = FALSE
      )
    }
  }
}

fisher_results <- do.call(rbind, fisher_grid)
dmp_counts     <- do.call(rbind, dmp_count_grid)

# Adjust p-values within each threshold combination
fisher_results <- fisher_results %>%
  group_by(fdr_threshold, diff_threshold) %>%
  mutate(fisher_padj = p.adjust(fisher_p, method = "BH")) %>%
  ungroup()

write.table(fisher_results,
            file.path(OUTPUT_DIR, "tables", "DMP_threshold_fisher_results.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

write.table(dmp_counts,
            file.path(OUTPUT_DIR, "tables", "DMP_threshold_counts.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat(sprintf("  Fisher's tests completed: %d combinations × %d modules = %d tests\n",
            nrow(dmp_counts), length(real_modules), nrow(fisher_results)))

# B.4: Tile plot — FDR × |diff| → Fisher padj for key modules
key_modules <- c("yellow", "blue", "green", "red", "brown", "turquoise")
key_modules <- intersect(key_modules, real_modules)

fisher_key <- fisher_results %>%
  filter(module %in% key_modules) %>%
  mutate(
    neg_log10_padj = -log10(fisher_padj),
    fdr_label      = ifelse(fdr_threshold == 1, "No FDR filter",
                            paste0("FDR<", fdr_threshold)),
    diff_label     = paste0("|diff|>=", diff_threshold),
    sig_label      = ifelse(fisher_padj < 0.05, "*", "")
  )

p_tile <- ggplot(fisher_key,
                 aes(x = diff_label, y = fdr_label, fill = neg_log10_padj)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sig_label), size = 6, color = "white") +
  facet_wrap(~module, ncol = 3) +
  scale_fill_gradient2(low = "#3498DB", mid = "#FEE08B", high = "#E74C3C",
                       midpoint = -log10(0.05),
                       name = expression(-log[10](padj))) +
  labs(
    title = "Module DMP Enrichment Sensitivity to DMP Thresholds — D. laeve",
    subtitle = "Fisher's exact test for DMP enrichment per module. * = padj < 0.05",
    x = "Methylation Difference Threshold",
    y = "FDR Threshold"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        strip.text = element_text(face = "bold", size = 11))
save_both(p_tile, "fig_B1_module_DMP_enrichment_threshold_tile", 14, 10)

# B.5: Line plot — how number of DMP genes per module changes with |diff| threshold
# (at FDR < 0.05)
gene_counts_by_diff <- fisher_results %>%
  filter(fdr_threshold == 0.05) %>%
  select(diff_threshold, module, n_genes_with_dmp, pct_genes_with_dmp)

p_gene_line <- ggplot(gene_counts_by_diff,
                      aes(x = diff_threshold, y = n_genes_with_dmp,
                          color = module, group = module)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  scale_color_identity() +
  labs(
    title = "Genes with DMPs per Module at Varying |meth_diff| Thresholds (FDR<0.05)",
    subtitle = "D. laeve — How many module genes retain DMPs as we increase stringency",
    x = "Minimum |Methylation Difference|", y = "Number of Genes with DMPs"
  )
save_both(p_gene_line, "fig_B2_genes_with_DMP_vs_diff_threshold", 12, 8)

# B.6: Percentage version
p_pct_line <- ggplot(gene_counts_by_diff,
                     aes(x = diff_threshold, y = pct_genes_with_dmp,
                         color = module, group = module)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  scale_color_identity() +
  labs(
    title = "% Module Genes with DMPs at Varying |meth_diff| Thresholds (FDR<0.05)",
    subtitle = "D. laeve — Which modules retain proportionally more DMPs at high stringency?",
    x = "Minimum |Methylation Difference|", y = "% of Module Genes with DMPs"
  )
save_both(p_pct_line, "fig_B3_pct_genes_with_DMP_vs_diff_threshold", 12, 8)

# B.7: DMP count summary barplot across thresholds
p_dmp_counts <- ggplot(dmp_counts,
                       aes(x = factor(diff_threshold),
                           y = n_dmps,
                           fill = factor(fdr_threshold))) +
  geom_col(position = "dodge", alpha = 0.85) +
  scale_fill_brewer(palette = "Blues", name = "FDR Threshold",
                    labels = function(x) ifelse(x == "1", "No filter", paste0("<", x))) +
  labs(
    title = "Total DMPs Retained at Each Threshold Combination — D. laeve",
    x = "Minimum |Methylation Difference|", y = "Number of DMPs"
  )
save_both(p_dmp_counts, "fig_B4_total_DMPs_by_threshold", 12, 7)

# B.8: Region breakdown at different |diff| thresholds (FDR < 0.05)
region_by_threshold <- list()
for (diff_thr in diff_thresholds) {
  dmps_filt <- dmps[dmps$fdr <= 0.05 & dmps$abs_diff >= diff_thr, ]
  region_counts <- as.data.frame(table(dmps_filt$region))
  colnames(region_counts) <- c("region", "count")
  region_counts$diff_threshold <- diff_thr
  region_counts$pct <- 100 * region_counts$count / nrow(dmps_filt)
  region_by_threshold[[as.character(diff_thr)]] <- region_counts
}
region_df <- do.call(rbind, region_by_threshold)

p_region <- ggplot(region_df,
                   aes(x = factor(diff_threshold), y = pct, fill = region)) +
  geom_col(position = "stack") +
  scale_fill_brewer(palette = "Set2", name = "Genomic Region") +
  labs(
    title = "DMP Genomic Region Distribution at Different |diff| Thresholds (FDR<0.05)",
    subtitle = "D. laeve — Does the regional profile shift with stringency?",
    x = "Minimum |Methylation Difference|", y = "% of DMPs"
  )
save_both(p_region, "fig_B5_DMP_region_by_diff_threshold", 11, 7)


# #############################################################################
# PART C: COMBINED SUMMARY
# #############################################################################

cat("\n=== Part C: Combined Summary ===\n\n")

# C.1: Which modules are robust to stringent thresholds?
robust_summary <- fisher_results %>%
  mutate(is_sig = fisher_padj < 0.05) %>%
  group_by(module) %>%
  summarise(
    n_combos_tested = n(),
    n_combos_sig    = sum(is_sig),
    pct_robust      = round(100 * mean(is_sig), 1),
    sig_at_default  = fisher_padj[fdr_threshold == 0.05 & diff_threshold == 0.1] < 0.05,
    sig_at_strict   = fisher_padj[fdr_threshold == 0.01 & diff_threshold == 0.25] < 0.05,
    .groups = "drop"
  ) %>%
  arrange(desc(pct_robust))

cat("  Module robustness to DMP threshold variation:\n")
print(as.data.frame(robust_summary))

write.table(robust_summary,
            file.path(OUTPUT_DIR, "tables", "module_robustness_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# C.2: Enrichment terms that appear at padj<0.05 vs those only at padj<0.1
new_at_01 <- enrich_all %>%
  filter(padj >= 0.05 & padj < 0.1) %>%
  arrange(padj) %>%
  select(module, category, Term, Description, Overlap, padj)

if (nrow(new_at_01) > 0) {
  cat(sprintf("\n  Terms gained by relaxing to padj < 0.1: %d\n", nrow(new_at_01)))
  cat("  Top 20:\n")
  print(head(as.data.frame(new_at_01), 20))

  write.table(new_at_01,
              file.path(OUTPUT_DIR, "tables", "enrichment_terms_gained_at_padj_01.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
}

# C.3: Summary of current vs alternative thresholds for key modules
current_vs_alt <- fisher_results %>%
  filter(module %in% key_modules) %>%
  mutate(
    scenario = case_when(
      fdr_threshold == 1   & diff_threshold == 0.1  ~ "Lenient (no FDR, |diff|>0.1)",
      fdr_threshold == 0.05 & diff_threshold == 0.1  ~ "Default (FDR<0.05, |diff|>0.1)",
      fdr_threshold == 0.05 & diff_threshold == 0.2  ~ "Moderate (FDR<0.05, |diff|>0.2)",
      fdr_threshold == 0.01 & diff_threshold == 0.25 ~ "Strict (FDR<0.01, |diff|>0.25)",
      fdr_threshold == 0.001 & diff_threshold == 0.3 ~ "Very strict (FDR<0.001, |diff|>0.3)",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(scenario)) %>%
  select(scenario, module, n_dmps_in_module, n_genes_with_dmp,
         pct_genes_with_dmp, fisher_padj) %>%
  arrange(module, scenario)

cat("\n  Key modules across named scenarios:\n")
print(as.data.frame(current_vs_alt))

write.table(current_vs_alt,
            file.path(OUTPUT_DIR, "tables", "key_modules_scenario_comparison.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# C.4: Scenario comparison plot for key modules
current_vs_alt$scenario_f <- factor(current_vs_alt$scenario,
                                     levels = c("Lenient (no FDR, |diff|>0.1)",
                                                "Default (FDR<0.05, |diff|>0.1)",
                                                "Moderate (FDR<0.05, |diff|>0.2)",
                                                "Strict (FDR<0.01, |diff|>0.25)",
                                                "Very strict (FDR<0.001, |diff|>0.3)"))

p_scenario <- ggplot(current_vs_alt,
                     aes(x = scenario_f, y = -log10(fisher_padj),
                         color = module, group = module)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", alpha = 0.5) +
  scale_color_identity() +
  annotate("text", x = 0.7, y = -log10(0.05) + 0.15,
           label = "padj = 0.05", size = 3, color = "red") +
  labs(
    title = "Module DMP Enrichment Across Stringency Scenarios — D. laeve",
    subtitle = "How does each module's Fisher enrichment respond to tighter DMP thresholds?",
    x = NULL, y = expression(-log[10](Fisher~padj))
  ) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 9))
save_both(p_scenario, "fig_C1_module_enrichment_scenarios", 14, 8)


# #############################################################################
# FINAL REPORT
# #############################################################################

cat("\n================================================================\n")
cat("  THRESHOLD SENSITIVITY ANALYSIS COMPLETE\n")
cat("================================================================\n\n")

cat("OUTPUT FILES:\n")
cat(sprintf("  Tables: %s/tables/\n", OUTPUT_DIR))
cat(sprintf("  Plots:  %s/plots/\n", OUTPUT_DIR))
cat("\nTables generated:\n")
cat("  - enrichment_padj_sensitivity.tsv          (Part A)\n")
cat("  - DMP_threshold_fisher_results.tsv          (Part B)\n")
cat("  - DMP_threshold_counts.tsv                  (Part B)\n")
cat("  - module_robustness_summary.tsv              (Part C)\n")
cat("  - enrichment_terms_gained_at_padj_01.tsv    (Part C)\n")
cat("  - key_modules_scenario_comparison.tsv        (Part C)\n")
cat("\nPlots generated:\n")
cat("  Part A: Enrichment padj sensitivity\n")
cat("    fig_A1  Heatmap: modules × padj thresholds\n")
cat("    fig_A2  Line plot: terms gained at relaxed thresholds\n")
cat("    fig_A3  Category comparison: padj<0.05 vs padj<0.1\n")
cat("  Part B: DMP threshold sensitivity\n")
cat("    fig_B1  Tile: FDR × |diff| → Fisher padj per module\n")
cat("    fig_B2  Line: genes with DMPs vs |diff| threshold\n")
cat("    fig_B3  Line: % genes with DMPs vs |diff| threshold\n")
cat("    fig_B4  Bar: total DMPs by threshold combination\n")
cat("    fig_B5  Stacked: region distribution by |diff| threshold\n")
cat("  Part C: Combined\n")
cat("    fig_C1  Scenario comparison for key modules\n")
cat("================================================================\n")
