#!/usr/bin/env Rscript
# =============================================================================
# Candidate Genes for Cloning: Sulfation + Key Module Hub Genes
# =============================================================================
# Goal: Prioritize genes for experimental validation (cloning) based on:
#   1. Sulfation pathway membership (appears at relaxed WGCNA thresholds)
#   2. DMP/DMR status (differentially methylated)
#   3. Hub gene status (highly connected in WGCNA modules)
#   4. Differential expression (significant log2FC)
#   5. Module membership (yellow = DMP-enriched, red = sulfation-enriched)
#
# Outputs a ranked candidate list with biological annotations.
#
# Run: Rscript local/04_candidate_genes_for_cloning.R
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
INT_DIR    <- "results/03_integration/Tables"
OUTPUT_DIR <- "results/11_candidate_genes"

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
cat("  Candidate Genes for Cloning\n")
cat("================================================================\n\n")


# #############################################################################
# 1. LOAD ALL DATA
# #############################################################################

cat("=== Loading data ===\n\n")

# Module assignments
mod_assign <- read.delim(file.path(WGCNA_DIR, "all_gene_module_assignments.tsv"))
cat(sprintf("  Module assignments: %s genes\n", format(nrow(mod_assign), big.mark = ",")))

# Hub genes
hub_genes <- read.delim(file.path(WGCNA_DIR, "all_modules_hub_genes.tsv"))
cat(sprintf("  Hub genes: %d hubs, %d candidates\n",
            sum(hub_genes$IsHub), sum(hub_genes$IsCandidate)))

# DMP data
dmps <- read.delim(file.path(METH_DIR, "dmps_annotated.txt"))
dmps_cs <- read.delim(file.path(METH_DIR, "dmps_chipseeker_annotated.txt"))
cat(sprintf("  DMPs: %s\n", format(nrow(dmps), big.mark = ",")))

# Gene names
gene_names_file <- file.path(METH_DIR, "gene_info_with_names.txt")
gene_names <- if (file.exists(gene_names_file)) {
  read.delim(gene_names_file)
} else NULL

# DEGs (tail amputated vs control)
deg_file <- file.path(INT_DIR, "DESeq2_tail_amputated_vs_control.txt")
degs <- if (file.exists(deg_file)) read.delim(deg_file) else NULL
if (!is.null(degs)) cat(sprintf("  DEGs: %s genes\n", format(nrow(degs), big.mark = ",")))

# Sulfation gene list
sulf_file <- "results/06_threshold_sensitivity/tables/sulfation_genes_from_string.tsv"
sulf_genes <- if (file.exists(sulf_file)) read.delim(sulf_file) else NULL
if (!is.null(sulf_genes)) cat(sprintf("  Sulfation genes: %d\n", nrow(sulf_genes)))

# MXT master table for module-level DMP enrichment
master_file <- file.path(INT_DIR, "MXT_MASTER_module_methylation_summary.txt")
master <- if (file.exists(master_file)) read.delim(master_file) else NULL

# Hub methylation status
hub_meth_file <- file.path(INT_DIR, "MXT_hub_genes_methylation_status.txt")
hub_meth <- if (file.exists(hub_meth_file)) read.delim(hub_meth_file) else NULL


# #############################################################################
# 2. BUILD DMP SUMMARY PER GENE
# #############################################################################

cat("\n=== Building gene-level DMP summary ===\n\n")

# Use ChIPseeker gene IDs
dmps$gene_id <- dmps_cs$geneId
dmps$region  <- dmps_cs$annotation
dmps$region  <- gsub(" \\(.*", "", dmps$region)
dmps$region[grepl("Promoter", dmps$region)] <- "Promoter"
dmps$region[grepl("Exon", dmps$region)]     <- "Exon"
dmps$region[grepl("Intron", dmps$region)]   <- "Intron"
dmps$region[grepl("UTR", dmps$region)]      <- "UTR"
dmps$region[grepl("Downstream", dmps$region)] <- "Downstream"
dmps$region[grepl("Intergenic|Distal", dmps$region)] <- "Intergenic"

dmps_clean <- dmps[!is.na(dmps$gene_id) & dmps$gene_id != "", ]

gene_dmp <- dmps_clean %>%
  group_by(gene_id) %>%
  summarise(
    n_dmps         = n(),
    n_sig_dmps     = sum(fdr < 0.05),
    mean_meth_diff = mean(methylation_diff, na.rm = TRUE),
    max_abs_diff   = max(abs(methylation_diff), na.rm = TRUE),
    n_hyper        = sum(methylation_diff > 0),
    n_hypo         = sum(methylation_diff < 0),
    regions        = paste(unique(region), collapse = ";"),
    primary_region = names(sort(table(region), decreasing = TRUE))[1],
    .groups = "drop"
  )

cat(sprintf("  Genes with DMPs: %d\n", nrow(gene_dmp)))


# #############################################################################
# 3. BUILD MASTER CANDIDATE TABLE
# #############################################################################

cat("\n=== Building master candidate table ===\n\n")

# Start with all genes in WGCNA network
candidates <- mod_assign %>%
  select(gene_id = Gene, module = Module_Color)

# Add hub gene info
if (!is.null(hub_genes)) {
  hub_info <- hub_genes %>%
    select(gene_id = Gene, Module, IsHub, IsCandidate, Connectivity, ModuleMembership) %>%
    distinct(gene_id, .keep_all = TRUE)
  candidates <- merge(candidates, hub_info, by = "gene_id", all.x = TRUE)
  candidates$IsHub[is.na(candidates$IsHub)] <- FALSE
  candidates$IsCandidate[is.na(candidates$IsCandidate)] <- FALSE
}

# Add DMP info
candidates <- merge(candidates, gene_dmp, by = "gene_id", all.x = TRUE)
candidates$has_dmp <- !is.na(candidates$n_dmps) & candidates$n_dmps > 0
candidates$has_sig_dmp <- !is.na(candidates$n_sig_dmps) & candidates$n_sig_dmps > 0

# Add DEG info
if (!is.null(degs)) {
  deg_cols <- intersect(c("gene_id", "Gene", "log2FoldChange", "padj", "pvalue"), colnames(degs))
  if ("Gene" %in% colnames(degs) && !"gene_id" %in% colnames(degs)) {
    degs$gene_id <- degs$Gene
  }
  if ("gene_id" %in% colnames(degs)) {
    deg_info <- degs %>%
      select(any_of(c("gene_id", "log2FoldChange", "padj"))) %>%
      distinct(gene_id, .keep_all = TRUE)
    candidates <- merge(candidates, deg_info, by = "gene_id", all.x = TRUE)
    candidates$is_deg <- !is.na(candidates$padj) & candidates$padj < 0.05
  }
}

# Add sulfation annotation
if (!is.null(sulf_genes)) {
  candidates$is_sulfation <- candidates$gene_id %in% sulf_genes$gene_id
  # Add sulfation term count
  sulf_info <- sulf_genes %>% select(gene_id, n_sulfation_terms = n_terms)
  candidates <- merge(candidates, sulf_info, by = "gene_id", all.x = TRUE)
  candidates$n_sulfation_terms[is.na(candidates$n_sulfation_terms)] <- 0
} else {
  candidates$is_sulfation <- FALSE
  candidates$n_sulfation_terms <- 0
}

# Add gene names
if (!is.null(gene_names)) {
  name_cols <- intersect(c("gene_id", "display_name", "gene_name", "Name"), colnames(gene_names))
  if (length(name_cols) >= 2) {
    name_info <- gene_names[, name_cols]
    colnames(name_info)[1] <- "gene_id"
    if (ncol(name_info) >= 2) colnames(name_info)[2] <- "gene_name"
    name_info <- distinct(name_info, gene_id, .keep_all = TRUE)
    candidates <- merge(candidates, name_info, by = "gene_id", all.x = TRUE)
  }
}

cat(sprintf("  Total candidates: %s genes\n", format(nrow(candidates), big.mark = ",")))
cat(sprintf("  Sulfation genes in network: %d\n", sum(candidates$is_sulfation)))
cat(sprintf("  Hub genes: %d\n", sum(candidates$IsHub)))
cat(sprintf("  Genes with DMPs: %d\n", sum(candidates$has_dmp)))


# #############################################################################
# 4. SCORING SYSTEM FOR CANDIDATE PRIORITIZATION
# #############################################################################

cat("\n=== Scoring candidates ===\n\n")

candidates$score <- 0

# Hub gene: +3
candidates$score <- candidates$score + ifelse(candidates$IsHub, 3, 0)

# Candidate gene (high MM + high GS): +2
candidates$score <- candidates$score + ifelse(candidates$IsCandidate, 2, 0)

# Has significant DMPs (FDR<0.05): +3
candidates$score <- candidates$score + ifelse(candidates$has_sig_dmp, 3, 0)

# Has any DMPs: +1
candidates$score <- candidates$score + ifelse(candidates$has_dmp & !candidates$has_sig_dmp, 1, 0)

# Large methylation effect (|diff| > 0.25): +2
candidates$score <- candidates$score +
  ifelse(!is.na(candidates$max_abs_diff) & candidates$max_abs_diff > 0.25, 2, 0)

# DEG (significant expression change): +2
if ("is_deg" %in% colnames(candidates)) {
  candidates$score <- candidates$score + ifelse(!is.na(candidates$is_deg) & candidates$is_deg, 2, 0)
}

# Sulfation pathway: +4 (key biological interest)
candidates$score <- candidates$score + ifelse(candidates$is_sulfation, 4, 0)

# Yellow module (DMP-enriched): +2
candidates$score <- candidates$score + ifelse(candidates$module == "yellow", 2, 0)

# Red module (sulfation-enriched): +1
candidates$score <- candidates$score + ifelse(candidates$module == "red", 1, 0)

# High connectivity: +1 (top 20%)
if ("Connectivity" %in% colnames(candidates)) {
  conn_threshold <- quantile(candidates$Connectivity, 0.8, na.rm = TRUE)
  candidates$score <- candidates$score +
    ifelse(!is.na(candidates$Connectivity) & candidates$Connectivity > conn_threshold, 1, 0)
}

# Sort by score
candidates <- candidates[order(-candidates$score), ]

cat("  Scoring criteria:\n")
cat("    +4  Sulfation pathway gene\n")
cat("    +3  Hub gene (top 10% connectivity + |MM|>0.7)\n")
cat("    +3  Significant DMP (FDR<0.05)\n")
cat("    +2  WGCNA candidate (|MM|>0.7 + |GS|>0.3)\n")
cat("    +2  Large methylation effect (|diff|>0.25)\n")
cat("    +2  Differentially expressed (DEG padj<0.05)\n")
cat("    +2  Yellow module (DMP-enriched)\n")
cat("    +1  Red module (sulfation-enriched)\n")
cat("    +1  High connectivity (top 20%)\n")
cat("    +1  Has DMP (not FDR-significant)\n")

cat(sprintf("\n  Score distribution:\n"))
print(table(candidates$score))


# #############################################################################
# 5. TOP CANDIDATES
# #############################################################################

cat("\n=== Top Candidate Genes for Cloning ===\n\n")

# Top 50 candidates
top50 <- head(candidates, 50)

# Select columns for display
display_cols <- intersect(c("gene_id", "gene_name", "module", "score",
                            "IsHub", "is_sulfation", "n_dmps", "n_sig_dmps",
                            "max_abs_diff", "mean_meth_diff", "primary_region",
                            "log2FoldChange", "padj", "Connectivity",
                            "ModuleMembership", "n_sulfation_terms"), colnames(top50))

cat("  Top 30 candidates:\n")
print(as.data.frame(head(top50[, display_cols], 30)), row.names = FALSE)

# Save full list
write.table(candidates[candidates$score >= 3, display_cols],
            file.path(OUTPUT_DIR, "tables", "candidate_genes_scored.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Save top 50 specifically
write.table(top50[, display_cols],
            file.path(OUTPUT_DIR, "tables", "top50_candidates_for_cloning.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)


# #############################################################################
# 6. SULFATION-SPECIFIC CANDIDATES
# #############################################################################

cat("\n=== Sulfation Pathway Candidates ===\n\n")

sulf_candidates <- candidates %>%
  filter(is_sulfation) %>%
  arrange(desc(score))

cat(sprintf("  Sulfation genes in WGCNA network: %d\n", nrow(sulf_candidates)))
cat(sprintf("  With DMPs: %d\n", sum(sulf_candidates$has_dmp)))
cat(sprintf("  Hub genes: %d\n", sum(sulf_candidates$IsHub)))

cat("\n  Top sulfation candidates:\n")
print(as.data.frame(head(sulf_candidates[, display_cols], 20)), row.names = FALSE)

write.table(sulf_candidates[, display_cols],
            file.path(OUTPUT_DIR, "tables", "sulfation_candidates.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Sulfation genes by module
sulf_by_module <- sulf_candidates %>%
  group_by(module) %>%
  summarise(
    n_genes = n(),
    n_hub = sum(IsHub),
    n_dmp = sum(has_dmp),
    n_sig_dmp = sum(has_sig_dmp),
    mean_score = round(mean(score), 1),
    .groups = "drop"
  ) %>%
  arrange(desc(n_genes))

cat("\n  Sulfation genes by module:\n")
print(as.data.frame(sulf_by_module))

write.table(sulf_by_module,
            file.path(OUTPUT_DIR, "tables", "sulfation_by_module.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)


# #############################################################################
# 7. MODULE-SPECIFIC TOP CANDIDATES
# #############################################################################

cat("\n=== Module-Specific Top Candidates ===\n\n")

key_modules <- c("yellow", "red", "blue", "green", "turquoise", "brown")

for (mod in key_modules) {
  mod_top <- candidates %>%
    filter(module == mod) %>%
    head(10)

  if (nrow(mod_top) > 0) {
    cat(sprintf("  %s module (top 10):\n", toupper(mod)))
    print(as.data.frame(mod_top[, intersect(display_cols, colnames(mod_top))]), row.names = FALSE)
    cat("\n")

    write.table(
      candidates %>% filter(module == mod) %>% head(20) %>% select(any_of(display_cols)),
      file.path(OUTPUT_DIR, "tables", paste0("top20_", mod, "_module.tsv")),
      sep = "\t", row.names = FALSE, quote = FALSE
    )
  }
}


# #############################################################################
# 8. VISUALIZATIONS
# #############################################################################

cat("\n=== Generating plots ===\n\n")

# 8.1 Score distribution by module (top candidates)
top_for_plot <- candidates %>% filter(score >= 5)

if (nrow(top_for_plot) > 0) {
  p_score_mod <- ggplot(top_for_plot, aes(x = module, fill = module)) +
    geom_bar(alpha = 0.85) +
    scale_fill_identity() +
    labs(
      title = "High-Scoring Candidate Genes by WGCNA Module",
      subtitle = "D. laeve: Genes with score >= 5 (multiple lines of evidence)",
      x = "Module", y = "Number of Candidates"
    )
  save_both(p_score_mod, "fig_E1_candidates_by_module", 12, 7)
}

# 8.2 Top 30 candidates lollipop chart
top30 <- head(candidates, 30)
if (!is.null(top30$gene_name)) {
  top30$label <- ifelse(!is.na(top30$gene_name) & top30$gene_name != "",
                        paste0(top30$gene_id, " (", top30$gene_name, ")"),
                        top30$gene_id)
} else {
  top30$label <- top30$gene_id
}

p_lollipop <- ggplot(top30, aes(x = reorder(label, score), y = score, color = module)) +
  geom_segment(aes(xend = reorder(label, score), y = 0, yend = score), linewidth = 1) +
  geom_point(size = 3) +
  scale_color_identity() +
  coord_flip() +
  labs(
    title = "Top 30 Candidate Genes for Cloning",
    subtitle = "D. laeve: Ranked by composite score (methylation + expression + network + sulfation)",
    x = NULL, y = "Composite Score"
  ) +
  theme(axis.text.y = element_text(size = 8))
save_both(p_lollipop, "fig_E2_top30_candidates_lollipop", 14, 10)

# 8.3 Sulfation genes: methylation vs expression
if (sum(sulf_candidates$has_dmp) > 5 && "log2FoldChange" %in% colnames(sulf_candidates)) {
  sulf_with_data <- sulf_candidates %>%
    filter(has_dmp, !is.na(log2FoldChange))

  if (nrow(sulf_with_data) > 3) {
    p_sulf_scatter <- ggplot(sulf_with_data,
                             aes(x = mean_meth_diff, y = log2FoldChange,
                                 color = module, size = n_dmps)) +
      geom_point(alpha = 0.7) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
      scale_color_identity() +
      scale_size_continuous(name = "# DMPs", range = c(2, 8)) +
      labs(
        title = "Sulfation Genes: Methylation vs Expression",
        subtitle = "D. laeve: Genes in sulfation pathways with DMPs",
        x = "Mean Methylation Difference (Amputated - Control)",
        y = "log2 Fold Change (Expression)"
      )

    if (requireNamespace("ggrepel", quietly = TRUE)) {
      top_sulf <- sulf_with_data %>%
        arrange(desc(score)) %>%
        head(15)
      label_col <- if ("gene_name" %in% colnames(top_sulf)) "gene_name" else "gene_id"
      p_sulf_scatter <- p_sulf_scatter +
        ggrepel::geom_text_repel(
          data = top_sulf,
          aes_string(label = label_col),
          size = 3, max.overlaps = 15
        )
    }
    save_both(p_sulf_scatter, "fig_E3_sulfation_meth_vs_expr", 12, 9)
  }
}

# 8.4 Heatmap of top candidates: features matrix
top_hm <- head(candidates, 40)
feature_mat <- data.frame(
  Hub = as.integer(top_hm$IsHub),
  Sulfation = as.integer(top_hm$is_sulfation),
  Sig_DMP = as.integer(top_hm$has_sig_dmp),
  Any_DMP = as.integer(top_hm$has_dmp),
  Yellow_mod = as.integer(top_hm$module == "yellow"),
  Red_mod = as.integer(top_hm$module == "red")
)
if ("is_deg" %in% colnames(top_hm)) {
  feature_mat$DEG <- as.integer(!is.na(top_hm$is_deg) & top_hm$is_deg)
}
if ("max_abs_diff" %in% colnames(top_hm)) {
  feature_mat$Large_effect <- as.integer(!is.na(top_hm$max_abs_diff) & top_hm$max_abs_diff > 0.25)
}

rownames(feature_mat) <- if ("gene_name" %in% colnames(top_hm)) {
  ifelse(!is.na(top_hm$gene_name) & top_hm$gene_name != "",
         paste0(top_hm$gene_id, " (", top_hm$gene_name, ")"),
         top_hm$gene_id)
} else {
  top_hm$gene_id
}

feature_mat <- as.matrix(feature_mat)

# Module color annotation
ann_row <- data.frame(Module = top_hm$module, row.names = rownames(feature_mat))
ann_colors <- list(Module = setNames(top_hm$module, top_hm$module))

hm_expr <- quote(
  pheatmap(feature_mat,
           cluster_rows = FALSE, cluster_cols = FALSE,
           color = c("white", "#E74C3C"),
           legend = FALSE,
           annotation_row = ann_row,
           annotation_colors = ann_colors,
           main = "Top 40 Candidate Genes: Evidence Matrix",
           fontsize = 9, fontsize_row = 7)
)
save_pheatmap_both(hm_expr, "fig_E4_candidates_evidence_matrix", 10, 14)


# #############################################################################
# FINAL REPORT
# #############################################################################

cat("\n================================================================\n")
cat("  CANDIDATE GENE ANALYSIS COMPLETE\n")
cat("================================================================\n\n")

cat("KEY RESULTS:\n\n")

cat(sprintf("Total genes scored: %s\n", format(nrow(candidates), big.mark = ",")))
cat(sprintf("Genes with score >= 5: %d\n", sum(candidates$score >= 5)))
cat(sprintf("Genes with score >= 8: %d\n", sum(candidates$score >= 8)))
cat(sprintf("Genes with score >= 10: %d\n\n", sum(candidates$score >= 10)))

cat("Top 10 candidates for cloning:\n")
for (i in 1:min(10, nrow(candidates))) {
  row <- candidates[i, ]
  name <- if ("gene_name" %in% colnames(row) && !is.na(row$gene_name) && row$gene_name != "")
    paste0(" (", row$gene_name, ")") else ""
  sulf_tag <- if (row$is_sulfation) " [SULFATION]" else ""
  hub_tag <- if (row$IsHub) " [HUB]" else ""
  dmp_tag <- if (row$has_sig_dmp) paste0(" [DMP:", row$n_sig_dmps, "]") else ""

  cat(sprintf("  %2d. %s%s — %s module, score=%d%s%s%s\n",
              i, row$gene_id, name, row$module, row$score,
              sulf_tag, hub_tag, dmp_tag))
}

cat(sprintf("\nSulfation pathway summary:\n"))
cat(sprintf("  Sulfation genes in network: %d\n", sum(candidates$is_sulfation)))
cat(sprintf("  With DMPs: %d\n", sum(candidates$is_sulfation & candidates$has_dmp)))
cat(sprintf("  Hub + sulfation: %d\n", sum(candidates$is_sulfation & candidates$IsHub)))

cat("\n\nOUTPUT FILES:\n")
cat(sprintf("  Tables: %s/tables/\n", OUTPUT_DIR))
cat(sprintf("  Plots:  %s/plots/\n", OUTPUT_DIR))
cat("\n  Tables:\n")
cat("    candidate_genes_scored.tsv          — All genes with score >= 3\n")
cat("    top50_candidates_for_cloning.tsv    — Top 50 by composite score\n")
cat("    sulfation_candidates.tsv            — Sulfation pathway genes\n")
cat("    sulfation_by_module.tsv             — Sulfation gene module distribution\n")
cat("    top20_<module>_module.tsv           — Top 20 per key module\n")
cat("\n  Plots:\n")
cat("    fig_E1  Candidates by module\n")
cat("    fig_E2  Top 30 lollipop chart\n")
cat("    fig_E3  Sulfation: methylation vs expression\n")
cat("    fig_E4  Evidence matrix heatmap\n")
cat("================================================================\n")
