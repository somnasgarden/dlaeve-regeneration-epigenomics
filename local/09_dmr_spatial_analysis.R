#!/usr/bin/env Rscript
# =============================================================================
# Script 09: DMR Spatial Analysis & Regulatory Element Identification
# =============================================================================
# Purpose: Analyze Differentially Methylated Regions (DMRs) as stronger
#          regulatory evidence than individual DMPs. DMRs represent coordinated
#          multi-CpG methylation changes that are more likely to reflect
#          functional regulatory elements (enhancers, silencers, insulators).
#
# Key questions:
#   1. How do DMRs relate to DMPs — are they in the same loci?
#   2. What is the DMR size distribution and CpG density?
#   3. Do DMR genes overlap with DMP-enriched modules?
#   4. Are DMRs closer to TSSs or more distal?
#   5. Which genes have both DMPs and DMRs (strongest evidence)?
#   6. Are there DMR "hotspots" — genomic regions with multiple DMRs?
#   7. DMR-gene expression relationship
#
# Input:
#   results/01_methylation/dmrs_annotated.txt
#   results/01_methylation/dmrs_chipseeker_annotated.txt
#   results/01_methylation/dmps_annotated.txt
#   results/03_integration/Tables/MXT_gene_level_meth_vs_expression.txt
#   results/02_rnaseq/Part5_WGCNA/data/all_gene_module_assignments.tsv
#   results/02_rnaseq/Part2_DEGs/data/DEresults_amputated_vs_control.csv
#   DATA/derLaeGenome_eviann_annotations.tsv
#
# Output:  results/08_dmr_spatial/
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

cat("=== Script 09: DMR Spatial Analysis & Regulatory Element Identification ===\n\n")

outdir <- "results/08_dmr_spatial"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
fig_dir <- file.path(outdir, "figures")
tab_dir <- file.path(outdir, "tables")
dir.create(fig_dir, showWarnings = FALSE)
dir.create(tab_dir, showWarnings = FALSE)

# ── Load data ────────────────────────────────────────────────────────────────

cat("Loading data...\n")

dmrs <- read.delim("results/01_methylation/dmrs_annotated.txt", stringsAsFactors = FALSE)
cat("  DMRs:", nrow(dmrs), "\n")

dmrs_cs <- read.delim("results/01_methylation/dmrs_chipseeker_annotated.txt", stringsAsFactors = FALSE)
cat("  DMRs (ChIPseeker):", nrow(dmrs_cs), "\n")

dmps <- read.delim("results/01_methylation/dmps_annotated.txt", stringsAsFactors = FALSE)
cat("  DMPs:", nrow(dmps), "\n")

gene_meth <- read.delim("results/03_integration/Tables/MXT_gene_level_meth_vs_expression.txt",
                         stringsAsFactors = FALSE)
cat("  Gene-level meth:", nrow(gene_meth), "genes\n")

modules <- read.delim("results/02_rnaseq/Part5_WGCNA/data/all_gene_module_assignments.tsv",
                       stringsAsFactors = FALSE)

deg <- read.csv("results/02_rnaseq/Part2_DEGs/data/DEresults_amputated_vs_control.csv",
                 stringsAsFactors = FALSE)
colnames(deg)[1] <- "gene_id"

annot <- read.delim("/mnt/c/Users/rafae/Projects/DATA/derLaeGenome_eviann_annotations.tsv",
                     header = FALSE, stringsAsFactors = FALSE)
colnames(annot) <- c("gene_id", "gene_symbol", "description")
annot <- annot %>% distinct(gene_id, .keep_all = TRUE)

# ── Part 1: DMR basic statistics ─────────────────────────────────────────────

cat("\n--- Part 1: DMR basic statistics ---\n")

cat("Total DMRs:", nrow(dmrs), "\n")
cat("Total CpGs in DMRs:", sum(dmrs$nCG), "\n")
cat("Mean CpGs per DMR:", round(mean(dmrs$nCG), 1), "\n")
cat("Median CpGs per DMR:", median(dmrs$nCG), "\n")
cat("Mean DMR length:", round(mean(dmrs$width), 0), "bp\n")
cat("Median DMR length:", median(dmrs$width), "bp\n")
cat("Max DMR length:", max(dmrs$width), "bp\n")

# CpG density
dmrs$cpg_density <- dmrs$nCG / dmrs$width * 1000  # CpGs per kb
cat("Mean CpG density:", round(mean(dmrs$cpg_density), 1), "CpGs/kb\n")

# Direction
n_hyper <- sum(dmrs$diff_methy > 0)
n_hypo <- sum(dmrs$diff_methy < 0)
cat("\nDirection: ", n_hyper, " hyper (", round(n_hyper/nrow(dmrs)*100,1), "%), ",
    n_hypo, " hypo (", round(n_hypo/nrow(dmrs)*100,1), "%)\n")

bt <- binom.test(n_hypo, nrow(dmrs), p = 0.5)
cat("Binomial test (hypo bias): p=", format(bt$p.value, digits = 4), "\n")

# Genomic annotation
cat("\nGenomic region distribution:\n")
region_tab <- sort(table(dmrs$annotation), decreasing = TRUE)
for (r in names(region_tab)) {
  cat(sprintf("  %-15s %4d (%.1f%%)\n", r, region_tab[r], region_tab[r]/nrow(dmrs)*100))
}

# ── Part 2: DMR size classes ────────────────────────────────────────────────

cat("\n--- Part 2: DMR size classes ---\n")

dmrs$size_class <- cut(dmrs$width,
                       breaks = c(0, 200, 500, 1000, 2000, 5000, Inf),
                       labels = c("0-200bp", "200-500bp", "500bp-1kb", "1-2kb", "2-5kb", ">5kb"))

dmrs$cpg_class <- cut(dmrs$nCG,
                      breaks = c(0, 5, 10, 20, 50, Inf),
                      labels = c("3-5", "6-10", "11-20", "21-50", ">50"))

size_summary <- dmrs %>%
  group_by(size_class) %>%
  summarise(
    n = n(),
    mean_ncg = round(mean(nCG), 1),
    mean_meth = round(mean(abs(diff_methy)), 3),
    pct_hypo = round(sum(diff_methy < 0) / n() * 100, 1),
    .groups = "drop"
  )

cat("\nDMR size distribution:\n")
for (i in 1:nrow(size_summary)) {
  r <- size_summary[i, ]
  cat(sprintf("  %-12s n=%4d  mean_CpGs=%.1f  |meth|=%.3f  %%hypo=%.1f\n",
              r$size_class, r$n, r$mean_ncg, r$mean_meth, r$pct_hypo))
}

# ── Part 3: DMR-gene overlaps ───────────────────────────────────────────────

cat("\n--- Part 3: DMR-gene relationships ---\n")

# Genes with DMRs
dmr_genes <- unique(dmrs$nearest_gene)
cat("Unique genes near DMRs:", length(dmr_genes), "\n")

# Genes with DMPs
dmp_genes <- unique(gene_meth$gene_id)
cat("Unique genes with DMPs:", length(dmp_genes), "\n")

# Overlap
both <- intersect(dmr_genes, dmp_genes)
only_dmr <- setdiff(dmr_genes, dmp_genes)
only_dmp <- setdiff(dmp_genes, dmr_genes)
cat("Genes with BOTH DMRs + DMPs:", length(both), "\n")
cat("Genes with DMRs only:", length(only_dmr), "\n")
cat("Genes with DMPs only:", length(only_dmp), "\n")

# ── Part 4: DMR module enrichment ───────────────────────────────────────────

cat("\n--- Part 4: DMR enrichment across WGCNA modules ---\n")

# Merge DMR genes with modules
dmr_gene_data <- data.frame(gene_id = dmr_genes, has_dmr = TRUE, stringsAsFactors = FALSE) %>%
  right_join(modules %>% select(Gene, Module_Color) %>% rename(gene_id = Gene), by = "gene_id") %>%
  mutate(has_dmr = ifelse(is.na(has_dmr), FALSE, has_dmr))

# Module-level Fisher's test
module_colors <- unique(dmr_gene_data$Module_Color)
module_colors <- module_colors[module_colors != "grey"]

dmr_module_results <- data.frame()
for (mod in module_colors) {
  in_mod <- dmr_gene_data$Module_Color == mod
  a <- sum(in_mod & dmr_gene_data$has_dmr)
  b <- sum(in_mod & !dmr_gene_data$has_dmr)
  c <- sum(!in_mod & dmr_gene_data$has_dmr)
  d <- sum(!in_mod & !dmr_gene_data$has_dmr)

  mat <- matrix(c(a, c, b, d), nrow = 2)
  ft <- fisher.test(mat)

  dmr_module_results <- rbind(dmr_module_results, data.frame(
    module = mod,
    total_genes = a + b,
    n_with_dmr = a,
    pct_with_dmr = round(a / (a + b) * 100, 1),
    odds_ratio = round(ft$estimate, 3),
    pvalue = ft$p.value,
    stringsAsFactors = FALSE
  ))
}

dmr_module_results$padj <- p.adjust(dmr_module_results$pvalue, method = "BH")
dmr_module_results <- dmr_module_results %>% arrange(pvalue)

bg_pct <- round(sum(dmr_gene_data$has_dmr) / nrow(dmr_gene_data) * 100, 1)
cat("Background DMR rate:", bg_pct, "%\n\n")

for (i in 1:nrow(dmr_module_results)) {
  r <- dmr_module_results[i, ]
  sig <- ifelse(r$padj < 0.05, "***", ifelse(r$padj < 0.1, "*", ""))
  cat(sprintf("  %-12s %4d/%4d (%.1f%%)  OR=%.3f  padj=%.4f %s\n",
              r$module, r$n_with_dmr, r$total_genes, r$pct_with_dmr,
              r$odds_ratio, r$padj, sig))
}

write.table(dmr_module_results, file.path(tab_dir, "J01_dmr_module_enrichment.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ── Part 5: DMR distance to TSS ─────────────────────────────────────────────

cat("\n--- Part 5: DMR distance to TSS analysis ---\n")

cat("Median DMR distance to TSS:", median(abs(dmrs_cs$distanceToTSS)), "bp\n")
cat("Mean DMR distance to TSS:", round(mean(abs(dmrs_cs$distanceToTSS))), "bp\n")

# Compare with DMP distances
dmps_cs <- read.delim("results/01_methylation/dmps_chipseeker_annotated.txt", stringsAsFactors = FALSE)
cat("Median DMP distance to TSS:", median(abs(dmps_cs$distanceToTSS)), "bp\n")

wt <- wilcox.test(abs(dmrs_cs$distanceToTSS), abs(dmps_cs$distanceToTSS))
cat("Wilcoxon test (DMR vs DMP distance): p=", format(wt$p.value, digits = 4), "\n")

# Categorize DMR distance
dmrs_cs$dist_class <- cut(abs(dmrs_cs$distanceToTSS),
                          breaks = c(0, 1000, 3000, 10000, 50000, Inf),
                          labels = c("<1kb", "1-3kb", "3-10kb", "10-50kb", ">50kb"))

cat("\nDMR distance-to-TSS distribution:\n")
dist_tab <- table(dmrs_cs$dist_class)
for (d in names(dist_tab)) {
  cat(sprintf("  %-10s %4d (%.1f%%)\n", d, dist_tab[d], dist_tab[d]/nrow(dmrs_cs)*100))
}

# ── Part 6: Strongest evidence genes (DMP + DMR + DE) ──────────────────────

cat("\n--- Part 6: Triple-evidence genes (DMP + DMR + DE) ---\n")

# Genes with DMPs
dmp_info <- gene_meth %>% select(gene_id, n_dmps, mean_meth_diff, primary_region, quadrant)

# Genes with DMRs
dmr_info <- dmrs %>%
  group_by(nearest_gene) %>%
  summarise(
    n_dmrs = n(),
    total_dmr_cpgs = sum(nCG),
    mean_dmr_meth = mean(diff_methy),
    max_dmr_length = max(width),
    dmr_region = annotation[which.max(abs(diff_methy))],
    .groups = "drop"
  ) %>%
  rename(gene_id = nearest_gene)

# DEG info
deg_info <- deg %>%
  filter(!is.na(padj)) %>%
  select(gene_id, log2FoldChange, padj) %>%
  rename(deg_padj = padj, deg_log2FC = log2FoldChange)

# Module info
mod_info <- modules %>% select(Gene, Module_Color) %>% rename(gene_id = Gene)

# Merge everything
triple <- dmp_info %>%
  inner_join(dmr_info, by = "gene_id") %>%
  left_join(deg_info, by = "gene_id") %>%
  left_join(mod_info, by = "gene_id") %>%
  left_join(annot %>% select(gene_id, gene_symbol, description), by = "gene_id") %>%
  mutate(is_de = !is.na(deg_padj) & deg_padj < 0.05)

cat("Genes with both DMPs and DMRs:", nrow(triple), "\n")
cat("Of these, also differentially expressed:", sum(triple$is_de), "\n")

# Show top triple-evidence genes
triple_de <- triple %>%
  filter(is_de) %>%
  mutate(evidence_score = n_dmps + total_dmr_cpgs * 2 + abs(deg_log2FC) * 5) %>%
  arrange(desc(evidence_score))

cat("\nTop 30 triple-evidence genes (DMP + DMR + DE):\n")
show_n <- min(30, nrow(triple_de))
if (show_n > 0) {
  for (i in 1:show_n) {
    r <- triple_de[i, ]
    cat(sprintf("  %s (%s): %d DMPs, %d DMRs (%d CpGs), log2FC=%.2f (padj=%.4f), module=%s, region=%s\n",
                r$gene_id, ifelse(is.na(r$gene_symbol), "?", r$gene_symbol),
                r$n_dmps, r$n_dmrs, r$total_dmr_cpgs, r$deg_log2FC, r$deg_padj,
                ifelse(is.na(r$Module_Color), "none", r$Module_Color), r$primary_region))
  }
}

write.table(triple_de, file.path(tab_dir, "J02_triple_evidence_genes.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ── Part 7: DMR hotspots ────────────────────────────────────────────────────

cat("\n--- Part 7: DMR hotspots (genomic windows with multiple DMRs) ---\n")

# Create 100kb windows and count DMRs
dmrs$window_100kb <- paste(dmrs$seqnames, floor(dmrs$start / 100000) * 100000, sep = ":")

hotspots <- dmrs %>%
  group_by(window_100kb) %>%
  summarise(
    n_dmrs = n(),
    total_cpgs = sum(nCG),
    chr = seqnames[1],
    window_start = floor(start[1] / 100000) * 100000,
    mean_meth = mean(diff_methy),
    genes = paste(unique(nearest_gene), collapse = ","),
    n_genes = length(unique(nearest_gene)),
    .groups = "drop"
  ) %>%
  filter(n_dmrs >= 3) %>%
  arrange(desc(n_dmrs))

cat("100kb windows with 3+ DMRs:", nrow(hotspots), "\n")
cat("\nTop DMR hotspots:\n")
show_n <- min(20, nrow(hotspots))
for (i in 1:show_n) {
  r <- hotspots[i, ]
  # Get gene names
  gene_ids <- strsplit(r$genes, ",")[[1]]
  gene_names <- annot %>% filter(gene_id %in% gene_ids) %>% pull(gene_symbol)
  name_str <- paste(unique(gene_names[!is.na(gene_names)]), collapse = ", ")
  if (nchar(name_str) > 50) name_str <- paste0(substr(name_str, 1, 47), "...")

  cat(sprintf("  %s:%s  %d DMRs, %d CpGs, %d genes: %s\n",
              r$chr, format(r$window_start, big.mark = ","),
              r$n_dmrs, r$total_cpgs, r$n_genes, name_str))
}

write.table(hotspots, file.path(tab_dir, "J03_dmr_hotspots_100kb.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ── Part 8: DMR vs DMP concordance ──────────────────────────────────────────

cat("\n--- Part 8: DMR-DMP methylation direction concordance ---\n")

# For genes with both DMPs and DMRs, are directions concordant?
concordance <- triple %>%
  mutate(
    dmp_direction = ifelse(mean_meth_diff > 0, "hyper", "hypo"),
    dmr_direction = ifelse(mean_dmr_meth > 0, "hyper", "hypo"),
    concordant = dmp_direction == dmr_direction
  )

cat("Genes with both DMPs and DMRs:", nrow(concordance), "\n")
cat("Concordant direction:", sum(concordance$concordant), "(", round(mean(concordance$concordant)*100,1), "%)\n")
cat("Discordant direction:", sum(!concordance$concordant), "(", round(mean(!concordance$concordant)*100,1), "%)\n")

bt2 <- binom.test(sum(concordance$concordant), nrow(concordance), p = 0.5)
cat("Binomial test vs 50/50: p=", format(bt2$p.value, digits = 4), "\n")

# ── Part 9: DMR effect size vs gene expression ──────────────────────────────

cat("\n--- Part 9: DMR effect size vs gene expression ---\n")

# Merge DMRs with expression
dmr_expr <- dmr_info %>%
  left_join(deg_info, by = "gene_id") %>%
  filter(!is.na(deg_log2FC))

if (nrow(dmr_expr) > 0) {
  cor_test <- cor.test(dmr_expr$mean_dmr_meth, dmr_expr$deg_log2FC, method = "spearman")
  cat("DMR methylation vs expression (all genes):\n")
  cat("  Spearman rho=", round(cor_test$estimate, 4), ", p=", format(cor_test$p.value, digits = 4), "\n")
  cat("  N =", nrow(dmr_expr), "\n")

  # By region
  dmr_by_region <- dmrs %>%
    left_join(deg %>% select(gene_id, log2FoldChange) %>% rename(nearest_gene = gene_id), by = "nearest_gene") %>%
    filter(!is.na(log2FoldChange))

  regions <- c("Gene Body", "Intergenic", "Promoter")
  for (reg in regions) {
    sub <- dmr_by_region %>% filter(annotation == reg)
    if (nrow(sub) >= 10) {
      ct <- cor.test(sub$diff_methy, sub$log2FoldChange, method = "spearman")
      cat(sprintf("  %-15s n=%4d  rho=%.4f  p=%s\n",
                  reg, nrow(sub), ct$estimate, format(ct$p.value, digits = 4)))
    }
  }
}

# ── Part 10: Large DMRs — putative regulatory elements ──────────────────────

cat("\n--- Part 10: Large DMRs as putative regulatory elements ---\n")

large_dmrs <- dmrs %>%
  filter(nCG >= 15 | width >= 1000) %>%
  left_join(annot %>% select(gene_id, gene_symbol, description),
            by = c("nearest_gene" = "gene_id")) %>%
  left_join(deg_info, by = c("nearest_gene" = "gene_id")) %>%
  left_join(mod_info, by = c("nearest_gene" = "gene_id")) %>%
  arrange(desc(nCG))

cat("Large DMRs (>=15 CpGs OR >=1kb):", nrow(large_dmrs), "\n\n")

show_n <- min(30, nrow(large_dmrs))
for (i in 1:show_n) {
  r <- large_dmrs[i, ]
  de_str <- ifelse(!is.na(r$deg_padj) & r$deg_padj < 0.05,
                   sprintf("DE log2FC=%.2f", r$deg_log2FC), "not DE")
  cat(sprintf("  %s:%s-%s  %d CpGs, %dbp, meth=%.3f, %s  gene=%s (%s), %s, module=%s\n",
              r$seqnames, format(r$start, big.mark = ","), format(r$end, big.mark = ","),
              r$nCG, r$width, r$diff_methy, r$annotation,
              r$nearest_gene, ifelse(is.na(r$gene_symbol), "?", r$gene_symbol),
              de_str, ifelse(is.na(r$Module_Color), "none", r$Module_Color)))
}

write.table(large_dmrs, file.path(tab_dir, "J04_large_dmrs.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# FIGURES
# ══════════════════════════════════════════════════════════════════════════════

cat("\n--- Generating figures ---\n")

# ── Figure J01: DMR size and CpG density ─────────────────────────────────────

p1 <- ggplot(dmrs, aes(x = width, y = nCG)) +
  geom_point(aes(color = diff_methy > 0), alpha = 0.4, size = 1.5) +
  geom_smooth(method = "lm", color = "gray30", se = TRUE) +
  scale_color_manual(values = c("TRUE" = "#E74C3C", "FALSE" = "#3498DB"),
                     labels = c("Hypomethylated", "Hypermethylated"),
                     name = "Direction") +
  scale_x_log10(labels = comma) +
  labs(
    title = "DMR Size vs CpG Count",
    subtitle = sprintf("%d DMRs total (%d hyper, %d hypo)", nrow(dmrs), n_hyper, n_hypo),
    x = "DMR length (bp, log10)",
    y = "Number of CpGs in DMR"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(fig_dir, "J01_dmr_size_vs_cpg.png"), p1, width = 8, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "J01_dmr_size_vs_cpg.pdf"), p1, width = 8, height = 6)
cat("  Saved J01\n")

# ── Figure J02: DMR genomic region distribution ─────────────────────────────

# Parse ChIPseeker annotation
dmrs_cs$simple_region <- case_when(
  grepl("Promoter", dmrs_cs$annotation) ~ "Promoter",
  grepl("Intron", dmrs_cs$annotation) ~ "Intron",
  grepl("Exon", dmrs_cs$annotation) ~ "Exon",
  grepl("UTR", dmrs_cs$annotation) ~ "UTR",
  grepl("Downstream", dmrs_cs$annotation) ~ "Downstream",
  grepl("Intergenic", dmrs_cs$annotation) ~ "Intergenic",
  TRUE ~ "Other"
)

region_counts <- dmrs_cs %>%
  group_by(simple_region) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(pct = n / sum(n) * 100)

p2 <- ggplot(region_counts, aes(x = reorder(simple_region, -pct), y = pct, fill = simple_region)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = sprintf("%d\n(%.1f%%)", n, pct)), vjust = -0.1, size = 3.5) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "DMR Genomic Region Distribution (ChIPseeker)",
    subtitle = sprintf("%d DMRs annotated", nrow(dmrs_cs)),
    x = "", y = "Percentage of DMRs"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(fig_dir, "J02_dmr_region_distribution.png"), p2, width = 8, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "J02_dmr_region_distribution.pdf"), p2, width = 8, height = 6)
cat("  Saved J02\n")

# ── Figure J03: DMR module enrichment ────────────────────────────────────────

mod_plot <- dmr_module_results %>%
  mutate(
    sig = ifelse(padj < 0.05, "padj < 0.05", ifelse(padj < 0.1, "padj < 0.1", "NS")),
    module = factor(module, levels = module[order(pct_with_dmr)])
  )

p3 <- ggplot(mod_plot, aes(x = pct_with_dmr, y = module)) +
  geom_vline(xintercept = bg_pct, linetype = "dashed", color = "gray50") +
  geom_segment(aes(x = bg_pct, xend = pct_with_dmr, y = module, yend = module),
               color = "gray70") +
  geom_point(aes(size = total_genes, color = sig)) +
  scale_color_manual(values = c("padj < 0.05" = "firebrick", "padj < 0.1" = "darkorange", "NS" = "gray40")) +
  labs(
    title = "DMR Enrichment by WGCNA Module",
    subtitle = sprintf("Background: %.1f%% genes have nearby DMR", bg_pct),
    x = "% genes with DMR",
    y = "",
    color = "Significance",
    size = "Module size"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(fig_dir, "J03_dmr_module_enrichment.png"), p3, width = 9, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "J03_dmr_module_enrichment.pdf"), p3, width = 9, height = 6)
cat("  Saved J03\n")

# ── Figure J04: DMR vs DMP effect concordance ───────────────────────────────

if (nrow(concordance) > 0) {
  p4 <- ggplot(concordance, aes(x = mean_meth_diff, y = mean_dmr_meth)) +
    geom_hline(yintercept = 0, color = "gray70") +
    geom_vline(xintercept = 0, color = "gray70") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(aes(color = concordant), alpha = 0.5, size = 2) +
    scale_color_manual(values = c("TRUE" = "#2ECC71", "FALSE" = "#E74C3C"),
                       labels = c("Discordant", "Concordant"),
                       name = "Direction") +
    labs(
      title = "DMP vs DMR Methylation Direction per Gene",
      subtitle = sprintf("%.1f%% concordant (p=%s, binomial test)",
                         mean(concordance$concordant)*100, format(bt2$p.value, digits = 3)),
      x = "Mean DMP methylation difference",
      y = "Mean DMR methylation difference"
    ) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))

  ggsave(file.path(fig_dir, "J04_dmp_dmr_concordance.png"), p4, width = 8, height = 7, dpi = 300)
  ggsave(file.path(fig_dir, "J04_dmp_dmr_concordance.pdf"), p4, width = 8, height = 7)
  cat("  Saved J04\n")
}

# ── Figure J05: DMR distance to TSS histogram ───────────────────────────────

p5 <- ggplot(dmrs_cs, aes(x = abs(distanceToTSS) / 1000)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "white", alpha = 0.8) +
  geom_vline(xintercept = median(abs(dmrs_cs$distanceToTSS)) / 1000,
             linetype = "dashed", color = "firebrick") +
  scale_x_log10(labels = comma) +
  annotate("text", x = median(abs(dmrs_cs$distanceToTSS)) / 1000 * 2, y = Inf,
           label = sprintf("Median = %.1f kb", median(abs(dmrs_cs$distanceToTSS)) / 1000),
           vjust = 2, color = "firebrick", fontface = "bold") +
  labs(
    title = "DMR Distance to Nearest TSS",
    subtitle = sprintf("%d DMRs | Most DMRs are distal (>3kb from TSS)", nrow(dmrs_cs)),
    x = "Distance to TSS (kb, log10)",
    y = "Count"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(fig_dir, "J05_dmr_distance_to_tss.png"), p5, width = 8, height = 5, dpi = 300)
ggsave(file.path(fig_dir, "J05_dmr_distance_to_tss.pdf"), p5, width = 8, height = 5)
cat("  Saved J05\n")

# ── Figure J06: DMR effect size vs expression ────────────────────────────────

if (nrow(dmr_expr) > 0) {
  p6 <- ggplot(dmr_expr, aes(x = mean_dmr_meth, y = deg_log2FC)) +
    geom_hline(yintercept = 0, color = "gray70") +
    geom_vline(xintercept = 0, color = "gray70") +
    geom_point(alpha = 0.3, size = 1.5, color = "steelblue") +
    geom_smooth(method = "lm", color = "firebrick", se = TRUE) +
    labs(
      title = "DMR Methylation Change vs Gene Expression Change",
      subtitle = sprintf("Spearman rho=%.4f, p=%s, n=%d",
                         cor_test$estimate, format(cor_test$p.value, digits = 3), nrow(dmr_expr)),
      x = "Mean DMR methylation difference",
      y = "log2 Fold Change (expression)"
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))

  ggsave(file.path(fig_dir, "J06_dmr_meth_vs_expression.png"), p6, width = 8, height = 6, dpi = 300)
  ggsave(file.path(fig_dir, "J06_dmr_meth_vs_expression.pdf"), p6, width = 8, height = 6)
  cat("  Saved J06\n")
}

# ══════════════════════════════════════════════════════════════════════════════
# SUMMARY
# ══════════════════════════════════════════════════════════════════════════════

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("SUMMARY: DMR Spatial Analysis\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat("1. DMR landscape:\n")
cat("   - ", nrow(dmrs), " DMRs covering ", sum(dmrs$nCG), " CpGs\n")
cat("   - Median length:", median(dmrs$width), "bp | Median CpGs:", median(dmrs$nCG), "\n")
cat("   - Direction:", n_hyper, "hyper,", n_hypo, "hypo (binomial p=",
    format(bt$p.value, digits = 3), ")\n\n")

cat("2. DMR-DMP overlap:\n")
cat("   -", length(both), "genes have both DMPs and DMRs (strongest evidence)\n")
cat("   - DMP-DMR concordance:", round(mean(concordance$concordant)*100,1),
    "% (p=", format(bt2$p.value, digits = 3), ")\n\n")

cat("3. Triple-evidence genes (DMP + DMR + DE):", sum(triple$is_de), "\n\n")

cat("4. DMRs are predominantly distal:\n")
cat("   - Median distance to TSS:", median(abs(dmrs_cs$distanceToTSS)), "bp\n")
cat("   - DMRs are further from TSSs than DMPs (p=",
    format(wt$p.value, digits = 3), ")\n\n")

cat("5. DMR-expression correlation:\n")
cat("   - Overall rho=", round(cor_test$estimate, 4), ", p=",
    format(cor_test$p.value, digits = 3), "\n\n")

cat("=== Script 09 complete ===\n")
cat("Output:", outdir, "\n")
cat("Figures:", length(list.files(fig_dir, pattern = "\\.png$")), "PNG +",
    length(list.files(fig_dir, pattern = "\\.pdf$")), "PDF\n")
cat("Tables:", length(list.files(tab_dir, pattern = "\\.tsv$")), "TSV\n")
