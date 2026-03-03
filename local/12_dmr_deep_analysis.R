#!/usr/bin/env Rscript
# =============================================================================
# Script 12: DMR Deep Analysis — All DMR-Based Questions
# =============================================================================
# Addresses three user questions:
#   1. "deep analysis only uses dmps, what if we use dmrs?"
#      → DMR-based functional category enrichment, module enrichment, TF analysis
#   2. "turquoise tf more methylated, does it pass the fisher tests?"
#      → Per-module TF Fisher's exact tests (DMP and DMR)
#   3. "DMRs with many CpGs — which genes, regions, CpG%, O/E, expression?"
#      → Case-by-case analysis of large DMRs (nCG >= 20)
#   4. "TFs with DMRs?"
#      → Full TF-DMR cross-analysis
#
# Input: results/ directory
# Output: results/09_dmr_deep_analysis/
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

cat("=== Script 12: DMR Deep Analysis ===\n\n")

# ── Output dirs ──────────────────────────────────────────────────────────────

outdir <- "results/09_dmr_deep_analysis"
fig_dir <- file.path(outdir, "figures")
tab_dir <- file.path(outdir, "tables")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

# ── Professional theme ───────────────────────────────────────────────────────

theme_pub <- theme_minimal(base_size = 13) +
  theme(
    text = element_text(color = "gray10"),
    plot.title = element_text(size = 15, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 11, color = "gray30"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    legend.text = element_text(size = 10),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray92"),
    strip.text = element_text(size = 11, face = "bold")
  )

save_both <- function(p, name, w = 8, h = 6) {
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
  cat(sprintf("  Saved: %s\n", name))
}

mod_palette <- c(
  yellow = "#F1C40F", blue = "#2980B9", red = "#E74C3C",
  green = "#27AE60", brown = "#8B4513", turquoise = "#1ABC9C",
  black = "#2C3E50", pink = "#E91E63", magenta = "#9C27B0",
  purple = "#7B1FA2", greenyellow = "#CDDC39", tan = "#D2B48C", grey = "#999999"
)

# ── Load data ────────────────────────────────────────────────────────────────

cat("Loading data...\n")

# DMRs
dmrs <- read.delim("results/01_methylation/dmrs_annotated.txt", stringsAsFactors = FALSE)
cat("  DMRs:", nrow(dmrs), "\n")

# Gene-level methylation (DMP-based)
gene_meth <- read.delim("results/03_integration/Tables/MXT_gene_level_meth_vs_expression.txt",
                         stringsAsFactors = FALSE)
cat("  Genes with DMPs:", nrow(gene_meth), "\n")

# Module assignments
modules <- read.delim("results/02_rnaseq/Part5_WGCNA/data/all_gene_module_assignments.tsv",
                       stringsAsFactors = FALSE)
cat("  Module assignments:", nrow(modules), "\n")

# DEGs
deg <- read.csv("results/02_rnaseq/Part2_DEGs/data/DEresults_amputated_vs_control.csv",
                 stringsAsFactors = FALSE)
colnames(deg)[1] <- "gene_id"
cat("  DEG results:", nrow(deg), "\n")

# EviAnn annotations
annot <- read.delim("/mnt/c/Users/rafae/Projects/DATA/derLaeGenome_eviann_annotations.tsv",
                     header = FALSE, stringsAsFactors = FALSE)
colnames(annot) <- c("gene_id", "gene_symbol", "description")
annot <- annot %>% distinct(gene_id, .keep_all = TRUE)
cat("  Annotations:", nrow(annot), "\n")

# Promoter CpG classification (for O/E index)
prom <- read.delim("results/04_promoter/01_classification/promoter_classification_all.tsv",
                    stringsAsFactors = FALSE)
cat("  Promoter classifications:", nrow(prom), "\n")

# Hub genes
hubs <- read.delim("results/02_rnaseq/Part5_WGCNA/data/hub_genes_all_modules.tsv",
                    stringsAsFactors = FALSE)

# TF predictions
tf_file <- "/mnt/c/Users/rafae/Projects/DATA/prediction_result.txt"
if (!file.exists(tf_file)) tf_file <- "results/05_tf_prediction/prediction_result.txt"
tf_pred <- read.delim(tf_file, stringsAsFactors = FALSE)

# ── Build DMR-gene mapping ──────────────────────────────────────────────────

cat("\nBuilding DMR-gene mapping...\n")

# Annotate DMRs with gene symbols, expression, modules
dmrs <- dmrs %>%
  left_join(annot, by = c("nearest_gene" = "gene_id")) %>%
  left_join(deg %>% select(gene_id, log2FoldChange, padj) %>%
              rename(deg_log2FC = log2FoldChange, deg_padj = padj),
            by = c("nearest_gene" = "gene_id")) %>%
  left_join(modules %>% select(Gene, Module_Color) %>%
              rename(nearest_gene = Gene), by = "nearest_gene") %>%
  left_join(prom %>% select(gene_id, cpg_oe_whole, gc_content, cpg_per_bp, promoter_class) %>%
              rename(nearest_gene = gene_id), by = "nearest_gene")

# Classify DMR direction
dmrs$dmr_direction <- ifelse(dmrs$diff_methy > 0, "Hyper", "Hypo")

# Gene-level DMR summary
gene_dmr <- dmrs %>%
  group_by(nearest_gene) %>%
  summarise(
    n_dmrs = n(),
    total_cpgs = sum(nCG),
    mean_dmr_meth = mean(diff_methy),
    max_nCG = max(nCG),
    max_width = max(width),
    primary_region = names(which.max(table(annotation))),
    dmr_direction = ifelse(mean(diff_methy) > 0, "Hyper", "Hypo"),
    .groups = "drop"
  ) %>%
  left_join(annot, by = c("nearest_gene" = "gene_id")) %>%
  left_join(modules %>% select(Gene, Module_Color) %>%
              rename(nearest_gene = Gene), by = "nearest_gene") %>%
  left_join(deg %>% select(gene_id, log2FoldChange, padj) %>%
              rename(deg_log2FC = log2FoldChange, deg_padj = padj),
            by = c("nearest_gene" = "gene_id")) %>%
  left_join(prom %>% select(gene_id, cpg_oe_whole, gc_content, promoter_class) %>%
              rename(nearest_gene = gene_id), by = "nearest_gene")

cat("  Genes with DMRs:", nrow(gene_dmr), "\n")

# ── Functional category classification ───────────────────────────────────────

cat("\nClassifying genes into functional categories...\n")

classify_gene <- function(symbol, desc) {
  txt <- tolower(paste(symbol, desc))
  if (grepl("transcription factor|zinc finger|homeobox|homeodomain|bhlh|hmg|forkhead|ets domain|pax|sox|gata|maf|myb|nuclear receptor|t-box|tbx|hox|nf-kappa|nfkb|smad|stat[0-9]|runx|irf|fox[a-z]", txt))
    return("Transcription_Factor")
  if (grepl("morphogen|homeobox|hox|wnt|hedgehog|notch|bmp|pattern|axis|anterior|posterior|dorsal|ventral|segment|gastrulat|limb|organizer|embryo[^n]|polarity|developmental|differentiat|regenerat|blastema|stem cell|pluripoten", txt))
    return("Morphogenesis_Development")
  if (grepl("kinase|phosphatas|receptor|signal|mapk|ras|rho|rac|cdc42|gtpase|camp|cgmp|calcium|calmodulin|notch|wnt|tgf|smad|jak|stat|erk|jnk|p38|pi3k|akt|mtor|growth factor|cytokine", txt))
    return("Signaling")
  if (grepl("cyclin|cdk|mitosis|meiosis|cell cycle|cell division|spindle|centrosome|kinetochore|checkpoint|proliferat|centriole", txt))
    return("Cell_Cycle_Proliferation")
  if (grepl("histone|chromatin|methyltrans|demethylas|acetyltrans|deacetylas|bromo|chromo|polycomb|trithorax|swi/snf|nucleosome|epigenet|dna methyl|hdac|hat |hmt |sirt|ezh|set domain|jmj|tet[0-9]|dnmt", txt))
    return("Chromatin_Epigenetic")
  if (grepl("cytoskelet|actin|myosin|tubulin|kinesin|dynein|intermediate filament|microtubul|motor protein|tropomyosin|spectrin|lamin|cell motil|migrat|chemotax|adhesion|cadherin|integrin|focal adhesion", txt))
    return("Cytoskeleton_Motility")
  if (grepl("apoptosis|programmed cell death|caspase|bcl-|bax|death domain|tnf|fas ligand|trail|autophagy|necroptosis", txt))
    return("Apoptosis")
  if (grepl("ribosom|spliceosom|mrna|rrna|trna|snrna|translation|splicing|polyadenyl|deadenyl|rna helicase|rna polymer|rna bind|rna process|exosome complex", txt))
    return("RNA_Processing")
  if (grepl("dna repair|dna damage|base excision|nucleotide excision|mismatch repair|double.strand break|homologous recomb|non-homologous|rad[0-9]|brca|atm |atr |checkpoint kinase|poly.adp", txt))
    return("DNA_Repair")
  if (grepl("immune|immun|toll|interleukin|interferon|mhc|antigen|antibod|complement|innate|adaptive|inflamm|cytokine|chemokine|nfkb|defensin|lectin|phagocyt|macrophage", txt))
    return("Immune_Defense")
  if (grepl("metabol|enzyme|oxidoreduct|transferase|hydrolase|lyase|isomerase|ligase|cytochrome|oxidase|reductase|synthase|dehydrogenase|carboxylase|kinase(?!.*receptor)|atp(?!.*transport)|glycolysis|tca|krebs|pentose|fatty acid|lipid biosynth|amino acid biosynth", txt, perl = TRUE))
    return("Metabolism")
  if (grepl("transport|channel|pump|carrier|symport|antiport|abc transport|solute carrier|ion channel|aquaporin|vesicl|endocyt|exocyt|secretion|golgi|endoplasmic", txt))
    return("Transport_Channel")
  if (grepl("collagen|extracellular matrix|ecm|laminin|fibronectin|proteoglycan|glycosaminoglycan|elastin|basement membrane|matrix metalloprot|structural", txt))
    return("ECM_Structural")
  return("Other_Unknown")
}

# Classify all genes (not just those with DMPs/DMRs)
all_genes <- annot %>%
  mutate(category = mapply(classify_gene, gene_symbol, description))
cat("  Classified", nrow(all_genes), "genes\n")

# Mark which genes have DMRs
all_genes$has_dmr <- all_genes$gene_id %in% gene_dmr$nearest_gene
all_genes$has_dmp <- all_genes$gene_id %in% gene_meth$gene_id

# Add module
all_genes <- all_genes %>%
  left_join(modules %>% select(Gene, Module_Color), by = c("gene_id" = "Gene"))

cat("  Genes with DMRs:", sum(all_genes$has_dmr), "\n")
cat("  Genes with DMPs:", sum(all_genes$has_dmp), "\n\n")


# ══════════════════════════════════════════════════════════════════════════════
# PART 1: DMR-based functional category enrichment
# ══════════════════════════════════════════════════════════════════════════════

cat("=== PART 1: DMR-Based Functional Category Enrichment ===\n")

categories <- unique(all_genes$category)
bg_rate <- mean(all_genes$has_dmr)
cat("  Background DMR rate:", round(bg_rate * 100, 1), "%\n\n")

dmr_enrichment <- data.frame()
for (cat_name in categories) {
  in_cat <- all_genes$category == cat_name
  ct <- table(in_cat, all_genes$has_dmr)
  ft <- fisher.test(ct)

  dmr_enrichment <- rbind(dmr_enrichment, data.frame(
    category = cat_name,
    n_genes = sum(in_cat),
    n_with_dmr = sum(in_cat & all_genes$has_dmr),
    pct_with_dmr = round(mean(all_genes$has_dmr[in_cat]) * 100, 1),
    background_pct = round(bg_rate * 100, 1),
    odds_ratio = round(ft$estimate, 3),
    pvalue = ft$p.value,
    direction = ifelse(ft$estimate > 1, "Enriched", "Depleted"),
    stringsAsFactors = FALSE
  ))
}
dmr_enrichment$padj <- p.adjust(dmr_enrichment$pvalue, method = "BH")
dmr_enrichment <- dmr_enrichment %>% arrange(pvalue)

cat("DMR enrichment by functional category:\n")
for (i in 1:nrow(dmr_enrichment)) {
  r <- dmr_enrichment[i, ]
  sig <- ifelse(r$padj < 0.05, "***", ifelse(r$padj < 0.1, "*", ""))
  cat(sprintf("  %-25s n=%4d  DMR%%=%.1f%%  OR=%.3f  padj=%.4f %s %s\n",
              r$category, r$n_genes, r$pct_with_dmr, r$odds_ratio,
              r$padj, r$direction, sig))
}

write.table(dmr_enrichment, file.path(tab_dir, "L01_dmr_category_enrichment.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ── Compare DMP vs DMR enrichment side-by-side ───────────────────────────────

cat("\n--- DMP vs DMR enrichment comparison ---\n")

# DMP enrichment
dmp_bg_rate <- mean(all_genes$has_dmp)
dmp_enrichment <- data.frame()
for (cat_name in categories) {
  in_cat <- all_genes$category == cat_name
  ct <- table(in_cat, all_genes$has_dmp)
  ft <- fisher.test(ct)
  dmp_enrichment <- rbind(dmp_enrichment, data.frame(
    category = cat_name,
    dmp_or = round(ft$estimate, 3),
    dmp_padj = p.adjust(ft$p.value, method = "BH")[1],
    stringsAsFactors = FALSE
  ))
}
dmp_enrichment$dmp_padj <- p.adjust(sapply(categories, function(cn) {
  ct <- table(all_genes$category == cn, all_genes$has_dmp)
  fisher.test(ct)$p.value
}), method = "BH")

comparison <- dmr_enrichment %>%
  select(category, dmr_or = odds_ratio, dmr_padj = padj, dmr_dir = direction) %>%
  left_join(dmp_enrichment, by = "category") %>%
  mutate(dmp_dir = ifelse(dmp_or > 1, "Enriched", "Depleted"))

cat("\n  Category                  DMP_OR  DMR_OR  Match?\n")
for (i in 1:nrow(comparison)) {
  r <- comparison[i, ]
  match <- r$dmp_dir == r$dmr_dir
  cat(sprintf("  %-25s  %.3f   %.3f   %s\n",
              r$category, r$dmp_or, r$dmr_or, ifelse(match, "YES", "DIFFERENT")))
}

write.table(comparison, file.path(tab_dir, "L02_dmp_vs_dmr_enrichment.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Plot DMR enrichment
p1 <- ggplot(dmr_enrichment %>% filter(category != "Other_Unknown"),
             aes(x = reorder(category, odds_ratio), y = odds_ratio,
                 fill = direction)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  geom_text(aes(label = sprintf("%.2f", odds_ratio)),
            hjust = ifelse(dmr_enrichment$odds_ratio[dmr_enrichment$category != "Other_Unknown"] > 1, -0.1, 1.1),
            size = 3.5) +
  coord_flip() +
  scale_fill_manual(values = c("Enriched" = "#E74C3C", "Depleted" = "#3498DB")) +
  labs(title = "DMR Enrichment by Functional Category",
       subtitle = "Fisher's exact test | Dashed line = no enrichment (OR=1)",
       x = NULL, y = "Odds Ratio", fill = "Direction") +
  theme_pub
save_both(p1, "L01_dmr_category_enrichment", w = 10, h = 7)


# ══════════════════════════════════════════════════════════════════════════════
# PART 2: Per-Module TF Fisher's Tests (DMP and DMR)
# ══════════════════════════════════════════════════════════════════════════════

cat("\n\n=== PART 2: Per-Module TF Fisher's Tests ===\n")

# TF genes
tf_genes <- all_genes$gene_id[all_genes$category == "Transcription_Factor"]
cat("  Total TF genes:", length(tf_genes), "\n\n")

module_colors <- unique(modules$Module_Color)
module_colors <- module_colors[!is.na(module_colors) & module_colors != "grey"]

tf_fisher_results <- data.frame()
for (mod in module_colors) {
  mod_genes <- modules$Gene[modules$Module_Color == mod]
  mod_tfs <- intersect(mod_genes, tf_genes)
  mod_nontfs <- setdiff(mod_genes, tf_genes)

  if (length(mod_tfs) >= 5) {
    # DMP Fisher's test: TFs vs non-TFs within this module
    tf_dmp <- sum(mod_tfs %in% gene_meth$gene_id)
    tf_nodmp <- length(mod_tfs) - tf_dmp
    nontf_dmp <- sum(mod_nontfs %in% gene_meth$gene_id)
    nontf_nodmp <- length(mod_nontfs) - nontf_dmp

    ct_dmp <- matrix(c(tf_dmp, tf_nodmp, nontf_dmp, nontf_nodmp), nrow = 2)
    ft_dmp <- fisher.test(ct_dmp)

    # DMR Fisher's test: TFs vs non-TFs within this module
    dmr_gene_ids <- unique(gene_dmr$nearest_gene)
    tf_dmr <- sum(mod_tfs %in% dmr_gene_ids)
    tf_nodmr <- length(mod_tfs) - tf_dmr
    nontf_dmr <- sum(mod_nontfs %in% dmr_gene_ids)
    nontf_nodmr <- length(mod_nontfs) - nontf_dmr

    ct_dmr <- matrix(c(tf_dmr, tf_nodmr, nontf_dmr, nontf_nodmr), nrow = 2)
    ft_dmr <- fisher.test(ct_dmr)

    tf_fisher_results <- rbind(tf_fisher_results, data.frame(
      module = mod,
      n_tfs = length(mod_tfs),
      n_nontfs = length(mod_nontfs),
      tf_dmp_pct = round(tf_dmp / length(mod_tfs) * 100, 1),
      nontf_dmp_pct = round(nontf_dmp / length(mod_nontfs) * 100, 1),
      dmp_or = round(ft_dmp$estimate, 3),
      dmp_pvalue = ft_dmp$p.value,
      tf_dmr_pct = round(tf_dmr / length(mod_tfs) * 100, 1),
      nontf_dmr_pct = round(nontf_dmr / length(mod_nontfs) * 100, 1),
      dmr_or = round(ft_dmr$estimate, 3),
      dmr_pvalue = ft_dmr$p.value,
      stringsAsFactors = FALSE
    ))
  }
}

tf_fisher_results$dmp_padj <- p.adjust(tf_fisher_results$dmp_pvalue, method = "BH")
tf_fisher_results$dmr_padj <- p.adjust(tf_fisher_results$dmr_pvalue, method = "BH")
tf_fisher_results <- tf_fisher_results %>% arrange(dmp_pvalue)

cat("\nPer-module TF Fisher's test (DMP):\n")
cat("  Module       TFs   TF_DMP%%  NonTF_DMP%%  OR      padj\n")
for (i in 1:nrow(tf_fisher_results)) {
  r <- tf_fisher_results[i, ]
  sig <- ifelse(r$dmp_padj < 0.05, "***", ifelse(r$dmp_padj < 0.1, "*", ""))
  cat(sprintf("  %-12s %3d   %5.1f%%    %5.1f%%      %.3f   %.4f %s\n",
              r$module, r$n_tfs, r$tf_dmp_pct, r$nontf_dmp_pct, r$dmp_or, r$dmp_padj, sig))
}

cat("\nPer-module TF Fisher's test (DMR):\n")
cat("  Module       TFs   TF_DMR%%  NonTF_DMR%%  OR      padj\n")
for (i in 1:nrow(tf_fisher_results)) {
  r <- tf_fisher_results[i, ]
  sig <- ifelse(r$dmr_padj < 0.05, "***", ifelse(r$dmr_padj < 0.1, "*", ""))
  cat(sprintf("  %-12s %3d   %5.1f%%    %5.1f%%      %.3f   %.4f %s\n",
              r$module, r$n_tfs, r$tf_dmr_pct, r$nontf_dmr_pct, r$dmr_or, r$dmr_padj, sig))
}

write.table(tf_fisher_results, file.path(tab_dir, "L03_module_tf_fisher_tests.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Plot module TF DMP/DMR rates
tf_plot_data <- tf_fisher_results %>%
  select(module, tf_dmp_pct, nontf_dmp_pct, tf_dmr_pct, nontf_dmr_pct) %>%
  pivot_longer(-module, names_to = "metric", values_to = "pct") %>%
  mutate(
    gene_type = ifelse(grepl("^tf_", metric), "TF", "Non-TF"),
    meth_type = ifelse(grepl("dmp", metric), "DMP", "DMR")
  )

p2 <- ggplot(tf_plot_data, aes(x = module, y = pct, fill = gene_type)) +
  geom_col(position = "dodge", width = 0.7) +
  facet_wrap(~ meth_type, scales = "free_y") +
  scale_fill_manual(values = c("TF" = "#E74C3C", "Non-TF" = "#BDC3C7")) +
  labs(title = "TF vs Non-TF Methylation by Module",
       subtitle = "Fisher's exact test per module | DMP and DMR rates",
       x = NULL, y = "% Genes with Methylation Change", fill = "Gene Type") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_both(p2, "L02_module_tf_dmp_dmr_rates", w = 12, h = 6)


# ══════════════════════════════════════════════════════════════════════════════
# PART 3: TFs with DMRs — Full Analysis
# ══════════════════════════════════════════════════════════════════════════════

cat("\n\n=== PART 3: TFs with DMRs ===\n")

# Identify TF genes with DMRs
tf_with_dmr <- gene_dmr %>%
  filter(nearest_gene %in% tf_genes) %>%
  arrange(desc(total_cpgs))

cat("  TFs with DMRs:", nrow(tf_with_dmr), "\n")
cat("  TFs with DMPs:", sum(tf_genes %in% gene_meth$gene_id), "\n")
cat("  TFs with both:", sum(tf_with_dmr$nearest_gene %in% gene_meth$gene_id), "\n\n")

# Detailed TF-DMR table
tf_dmr_detail <- tf_with_dmr %>%
  mutate(
    is_de = !is.na(deg_padj) & deg_padj < 0.05,
    is_hub = nearest_gene %in% hubs$Gene[hubs$IsHub == TRUE],
    has_dmp = nearest_gene %in% gene_meth$gene_id
  ) %>%
  select(nearest_gene, gene_symbol, description, Module_Color,
         n_dmrs, total_cpgs, max_nCG, max_width, primary_region,
         mean_dmr_meth, dmr_direction,
         deg_log2FC, deg_padj, is_de, is_hub, has_dmp,
         cpg_oe_whole, gc_content, promoter_class) %>%
  arrange(desc(total_cpgs))

# Print top TFs with DMRs
cat("Top 20 TFs with most DMR CpGs:\n")
cat(sprintf("  %-15s %-12s %-8s %5s %5s %5s %-10s %6s %6s %-5s\n",
            "Gene_ID", "Symbol", "Module", "nDMR", "CpGs", "MaxCG", "Region",
            "MethD", "log2FC", "DE"))
for (i in 1:min(20, nrow(tf_dmr_detail))) {
  r <- tf_dmr_detail[i, ]
  cat(sprintf("  %-15s %-12s %-8s %5d %5d %5d %-10s %+.3f %+.3f %-5s\n",
              r$nearest_gene,
              ifelse(is.na(r$gene_symbol), "-", r$gene_symbol),
              ifelse(is.na(r$Module_Color), "-", r$Module_Color),
              r$n_dmrs, r$total_cpgs, r$max_nCG, r$primary_region,
              r$mean_dmr_meth,
              ifelse(is.na(r$deg_log2FC), NA, r$deg_log2FC),
              ifelse(r$is_de, "YES", "no")))
}

write.table(tf_dmr_detail, file.path(tab_dir, "L04_tfs_with_dmrs.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Per-module TF DMR summary
tf_dmr_by_module <- tf_dmr_detail %>%
  filter(!is.na(Module_Color)) %>%
  group_by(Module_Color) %>%
  summarise(
    n_tfs_with_dmr = n(),
    total_dmr_cpgs = sum(total_cpgs),
    mean_cpgs = round(mean(total_cpgs), 1),
    n_hyper = sum(dmr_direction == "Hyper"),
    n_hypo = sum(dmr_direction == "Hypo"),
    pct_hypo = round(sum(dmr_direction == "Hypo") / n() * 100, 1),
    n_de = sum(is_de),
    n_hub = sum(is_hub),
    .groups = "drop"
  ) %>%
  arrange(desc(n_tfs_with_dmr))

cat("\nTF DMRs by module:\n")
cat(sprintf("  %-12s %4s %5s %6s %5s %5s %6s %3s %3s\n",
            "Module", "nTF", "CpGs", "Mean", "Hyper", "Hypo", "%Hypo", "DE", "Hub"))
for (i in 1:nrow(tf_dmr_by_module)) {
  r <- tf_dmr_by_module[i, ]
  cat(sprintf("  %-12s %4d %5d %6.1f %5d %5d %5.1f%% %3d %3d\n",
              r$Module_Color, r$n_tfs_with_dmr, r$total_dmr_cpgs, r$mean_cpgs,
              r$n_hyper, r$n_hypo, r$pct_hypo, r$n_de, r$n_hub))
}

write.table(tf_dmr_by_module, file.path(tab_dir, "L05_tf_dmr_by_module.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)


# ══════════════════════════════════════════════════════════════════════════════
# PART 4: Large DMR Case Studies (nCG >= 20)
# ══════════════════════════════════════════════════════════════════════════════

cat("\n\n=== PART 4: Large DMR Case Studies (nCG >= 20) ===\n")

large_dmrs <- dmrs %>%
  filter(nCG >= 20) %>%
  arrange(desc(nCG)) %>%
  mutate(
    is_de = !is.na(deg_padj) & deg_padj < 0.05,
    cpg_density = nCG / width * 1000  # CpGs per kb
  )

cat("  DMRs with >= 20 CpGs:", nrow(large_dmrs), "\n\n")

# Separate hyper and hypo
large_hyper <- large_dmrs %>% filter(diff_methy > 0)
large_hypo <- large_dmrs %>% filter(diff_methy < 0)
cat("  Hypermethylated:", nrow(large_hyper), "\n")
cat("  Hypomethylated:", nrow(large_hypo), "\n\n")

# Case-by-case table
cat("CASE-BY-CASE: All large DMRs (nCG >= 20)\n")
cat(paste(rep("=", 120), collapse = ""), "\n")
cat(sprintf("%-5s %-15s %-12s %-10s %4s %+7s %6s %6s %+7s %-5s %-8s %-10s\n",
            "Rank", "Gene_ID", "Symbol", "Region", "nCG", "MethD", "CpG/kb", "CpG_OE",
            "log2FC", "DE", "Module", "PromClass"))
cat(paste(rep("-", 120), collapse = ""), "\n")

for (i in 1:nrow(large_dmrs)) {
  r <- large_dmrs[i, ]
  cat(sprintf("%-5d %-15s %-12s %-10s %4d %+.4f %6.1f %6.3f %+.4f %-5s %-8s %-10s\n",
              i, r$nearest_gene,
              ifelse(is.na(r$gene_symbol), "-", substr(r$gene_symbol, 1, 12)),
              ifelse(is.na(r$annotation), "-", substr(r$annotation, 1, 10)),
              r$nCG,
              r$diff_methy,
              r$cpg_density,
              ifelse(is.na(r$cpg_oe_whole), NA, r$cpg_oe_whole),
              ifelse(is.na(r$deg_log2FC), 0, r$deg_log2FC),
              ifelse(r$is_de, "YES", "no"),
              ifelse(is.na(r$Module_Color), "-", r$Module_Color),
              ifelse(is.na(r$promoter_class), "-", r$promoter_class)))
}

write.table(large_dmrs %>%
              select(seqnames, start, end, width, annotation, nearest_gene, gene_symbol,
                     description, nCG, diff_methy, dmr_direction, cpg_density,
                     cpg_oe_whole, gc_content, promoter_class,
                     deg_log2FC, deg_padj, is_de, Module_Color),
            file.path(tab_dir, "L06_large_dmrs_case_studies.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ── Summarize the hypermethylated large DMRs specifically ────────────────────

cat("\n\n--- HYPERMETHYLATED large DMRs (diff_methy > 0, nCG >= 20) ---\n")
if (nrow(large_hyper) > 0) {
  cat(sprintf("  N = %d DMRs\n", nrow(large_hyper)))
  cat(sprintf("  Mean nCG = %.1f, Mean |meth diff| = %.3f\n",
              mean(large_hyper$nCG), mean(abs(large_hyper$diff_methy))))
  cat(sprintf("  Regions: %s\n", paste(names(table(large_hyper$annotation)), "=",
                                        table(large_hyper$annotation), collapse = ", ")))
  cat("\n  Genes near hypermethylated large DMRs:\n")
  for (i in 1:nrow(large_hyper)) {
    r <- large_hyper[i, ]
    cat(sprintf("    %2d. %-15s %-15s nCG=%d meth=%+.3f log2FC=%+.3f %-10s OE=%.3f %s\n",
                i, r$nearest_gene,
                ifelse(is.na(r$gene_symbol), "-", r$gene_symbol),
                r$nCG, r$diff_methy,
                ifelse(is.na(r$deg_log2FC), 0, r$deg_log2FC),
                ifelse(is.na(r$annotation), "-", r$annotation),
                ifelse(is.na(r$cpg_oe_whole), NA, r$cpg_oe_whole),
                ifelse(is.na(r$description), "", substr(r$description, 1, 40))))
  }
}

cat("\n--- HYPOMETHYLATED large DMRs (diff_methy < 0, nCG >= 20) ---\n")
if (nrow(large_hypo) > 0) {
  cat(sprintf("  N = %d DMRs\n", nrow(large_hypo)))
  cat(sprintf("  Mean nCG = %.1f, Mean |meth diff| = %.3f\n",
              mean(large_hypo$nCG), mean(abs(large_hypo$diff_methy))))
}

# ── CpG O/E and density analysis of large DMRs ──────────────────────────────

cat("\n\n--- CpG O/E and density analysis of large DMRs ---\n")
large_with_oe <- large_dmrs %>% filter(!is.na(cpg_oe_whole))
if (nrow(large_with_oe) > 0) {
  cat(sprintf("  DMRs with promoter O/E data: %d / %d\n", nrow(large_with_oe), nrow(large_dmrs)))
  cat(sprintf("  Mean CpG O/E: %.3f (genome-wide mean ~%.3f)\n",
              mean(large_with_oe$cpg_oe_whole), mean(prom$cpg_oe_whole, na.rm = TRUE)))
  cat(sprintf("  Mean GC%%: %.1f%%\n", mean(large_with_oe$gc_content, na.rm = TRUE) * 100))

  # Compare hyper vs hypo O/E
  hyper_oe <- large_with_oe %>% filter(diff_methy > 0)
  hypo_oe <- large_with_oe %>% filter(diff_methy < 0)
  if (nrow(hyper_oe) >= 3 & nrow(hypo_oe) >= 3) {
    wt <- wilcox.test(hyper_oe$cpg_oe_whole, hypo_oe$cpg_oe_whole)
    cat(sprintf("  Hyper O/E: %.3f vs Hypo O/E: %.3f  Wilcoxon p=%.4f\n",
                mean(hyper_oe$cpg_oe_whole), mean(hypo_oe$cpg_oe_whole), wt$p.value))
  }
}

# Plot: large DMR CpG count vs methylation change
p3 <- ggplot(large_dmrs, aes(x = nCG, y = diff_methy,
                               color = dmr_direction, size = width)) +
  geom_point(alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  geom_text(data = large_dmrs %>% filter(nCG >= 27),
            aes(label = ifelse(is.na(gene_symbol), nearest_gene, gene_symbol)),
            size = 3, hjust = -0.1, vjust = 0.5, show.legend = FALSE) +
  scale_color_manual(values = c("Hyper" = "#E74C3C", "Hypo" = "#3498DB")) +
  scale_size_continuous(range = c(2, 8), name = "DMR Width (bp)") +
  labs(title = "Large DMRs: CpG Count vs Methylation Change",
       subtitle = "DMRs with >= 20 CpGs | Labels for nCG >= 27",
       x = "Number of CpGs in DMR", y = "Mean Methylation Difference",
       color = "Direction") +
  theme_pub
save_both(p3, "L03_large_dmr_cpg_vs_meth", w = 10, h = 7)

# Plot: large DMR region distribution
p4 <- ggplot(large_dmrs, aes(x = annotation, fill = dmr_direction)) +
  geom_bar(position = "dodge", width = 0.7) +
  scale_fill_manual(values = c("Hyper" = "#E74C3C", "Hypo" = "#3498DB")) +
  labs(title = "Large DMR Genomic Regions (nCG >= 20)",
       subtitle = "Distribution of hypermethylated vs hypomethylated large DMRs",
       x = "Genomic Region", y = "Count", fill = "Direction") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_both(p4, "L04_large_dmr_regions", w = 8, h = 6)


# ══════════════════════════════════════════════════════════════════════════════
# PART 5: DMR Module Enrichment with Direction
# ══════════════════════════════════════════════════════════════════════════════

cat("\n\n=== PART 5: DMR Module Enrichment with Direction ===\n")

dmr_module_detail <- data.frame()
for (mod in module_colors) {
  mod_genes <- modules$Gene[modules$Module_Color == mod]
  mod_with_dmr <- intersect(mod_genes, gene_dmr$nearest_gene)

  if (length(mod_genes) >= 20) {
    # Fisher's test: DMR enrichment
    in_mod <- all_genes$gene_id %in% mod_genes
    ct <- table(in_mod, all_genes$has_dmr)
    ft <- fisher.test(ct)

    # Direction
    mod_dmr_data <- gene_dmr %>% filter(nearest_gene %in% mod_with_dmr)
    n_hyper <- sum(mod_dmr_data$dmr_direction == "Hyper")
    n_hypo <- sum(mod_dmr_data$dmr_direction == "Hypo")

    # Binomial test for direction bias
    if ((n_hyper + n_hypo) >= 5) {
      bt <- binom.test(n_hypo, n_hyper + n_hypo, p = 0.5)
      dir_p <- bt$p.value
    } else {
      dir_p <- NA
    }

    dmr_module_detail <- rbind(dmr_module_detail, data.frame(
      module = mod,
      n_genes = length(mod_genes),
      n_with_dmr = length(mod_with_dmr),
      pct_dmr = round(length(mod_with_dmr) / length(mod_genes) * 100, 1),
      odds_ratio = round(ft$estimate, 3),
      pvalue = ft$p.value,
      n_hyper = n_hyper,
      n_hypo = n_hypo,
      pct_hypo = round(n_hypo / (n_hyper + n_hypo) * 100, 1),
      direction_p = dir_p,
      stringsAsFactors = FALSE
    ))
  }
}

dmr_module_detail$padj <- p.adjust(dmr_module_detail$pvalue, method = "BH")
dmr_module_detail$dir_padj <- p.adjust(dmr_module_detail$direction_p, method = "BH")
dmr_module_detail <- dmr_module_detail %>% arrange(pvalue)

cat("\nDMR enrichment per module:\n")
for (i in 1:nrow(dmr_module_detail)) {
  r <- dmr_module_detail[i, ]
  sig <- ifelse(r$padj < 0.05, "***", ifelse(r$padj < 0.1, "*", ""))
  dir_sig <- ifelse(!is.na(r$dir_padj) & r$dir_padj < 0.05, " DIR_BIASED", "")
  cat(sprintf("  %-12s n=%4d  DMR%%=%.1f%%  OR=%.3f  padj=%.4f %s  Hyper=%d Hypo=%d (%.1f%% hypo)%s\n",
              r$module, r$n_genes, r$pct_dmr, r$odds_ratio, r$padj, sig,
              r$n_hyper, r$n_hypo, r$pct_hypo, dir_sig))
}

write.table(dmr_module_detail, file.path(tab_dir, "L07_dmr_module_enrichment_detail.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)


# ══════════════════════════════════════════════════════════════════════════════
# PART 6: Hub Genes with DMRs
# ══════════════════════════════════════════════════════════════════════════════

cat("\n\n=== PART 6: Hub Gene DMR Analysis ===\n")

hub_dmr_results <- data.frame()
for (mod in module_colors) {
  mod_genes <- modules$Gene[modules$Module_Color == mod]
  mod_hubs <- hubs$Gene[hubs$Module == mod & hubs$IsHub == TRUE]
  mod_nonhubs <- setdiff(mod_genes, mod_hubs)

  if (length(mod_hubs) >= 5) {
    hub_dmr_rate <- mean(mod_hubs %in% gene_dmr$nearest_gene)
    nonhub_dmr_rate <- mean(mod_nonhubs %in% gene_dmr$nearest_gene)

    ct <- matrix(c(
      sum(mod_hubs %in% gene_dmr$nearest_gene),
      sum(!(mod_hubs %in% gene_dmr$nearest_gene)),
      sum(mod_nonhubs %in% gene_dmr$nearest_gene),
      sum(!(mod_nonhubs %in% gene_dmr$nearest_gene))
    ), nrow = 2)
    ft <- fisher.test(ct)

    hub_dmr_results <- rbind(hub_dmr_results, data.frame(
      module = mod,
      n_hubs = length(mod_hubs),
      hub_dmr_pct = round(hub_dmr_rate * 100, 1),
      nonhub_dmr_pct = round(nonhub_dmr_rate * 100, 1),
      odds_ratio = round(ft$estimate, 3),
      pvalue = ft$p.value,
      stringsAsFactors = FALSE
    ))
  }
}

hub_dmr_results$padj <- p.adjust(hub_dmr_results$pvalue, method = "BH")
hub_dmr_results <- hub_dmr_results %>% arrange(pvalue)

cat("\nHub vs non-hub DMR rates per module:\n")
for (i in 1:nrow(hub_dmr_results)) {
  r <- hub_dmr_results[i, ]
  sig <- ifelse(r$padj < 0.05, "***", ifelse(r$padj < 0.1, "*", ""))
  cat(sprintf("  %-12s hubs=%5.1f%%  nonhubs=%5.1f%%  OR=%.3f  padj=%.4f %s\n",
              r$module, r$hub_dmr_pct, r$nonhub_dmr_pct, r$odds_ratio, r$padj, sig))
}

write.table(hub_dmr_results, file.path(tab_dir, "L08_hub_dmr_enrichment.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)


# ══════════════════════════════════════════════════════════════════════════════
# PART 7: Global TF Fisher's test (DMP and DMR)
# ══════════════════════════════════════════════════════════════════════════════

cat("\n\n=== PART 7: Global TF Methylation Tests ===\n")

# DMP: TFs vs non-TFs globally
tf_has_dmp <- sum(tf_genes %in% gene_meth$gene_id)
tf_no_dmp <- length(tf_genes) - tf_has_dmp
nontf_has_dmp <- sum(!(all_genes$gene_id %in% tf_genes) & all_genes$has_dmp)
nontf_no_dmp <- sum(!(all_genes$gene_id %in% tf_genes) & !all_genes$has_dmp)
ct_global_dmp <- matrix(c(tf_has_dmp, tf_no_dmp, nontf_has_dmp, nontf_no_dmp), nrow = 2)
ft_global_dmp <- fisher.test(ct_global_dmp)

cat(sprintf("  Global TF DMP: OR=%.3f, p=%.2e\n", ft_global_dmp$estimate, ft_global_dmp$p.value))
cat(sprintf("    TFs: %d/%d = %.1f%%\n", tf_has_dmp, length(tf_genes),
            tf_has_dmp / length(tf_genes) * 100))

# DMR: TFs vs non-TFs globally
tf_has_dmr <- sum(tf_genes %in% gene_dmr$nearest_gene)
tf_no_dmr <- length(tf_genes) - tf_has_dmr
nontf_has_dmr <- sum(!(all_genes$gene_id %in% tf_genes) & all_genes$has_dmr)
nontf_no_dmr <- sum(!(all_genes$gene_id %in% tf_genes) & !all_genes$has_dmr)
ct_global_dmr <- matrix(c(tf_has_dmr, tf_no_dmr, nontf_has_dmr, nontf_no_dmr), nrow = 2)
ft_global_dmr <- fisher.test(ct_global_dmr)

cat(sprintf("  Global TF DMR: OR=%.3f, p=%.2e\n", ft_global_dmr$estimate, ft_global_dmr$p.value))
cat(sprintf("    TFs: %d/%d = %.1f%%\n", tf_has_dmr, length(tf_genes),
            tf_has_dmr / length(tf_genes) * 100))

# Test by TF subcategory
cat("\n  TF subcategory DMR enrichment:\n")
tf_subcats <- all_genes %>%
  filter(category == "Transcription_Factor") %>%
  mutate(subcat = case_when(
    grepl("zinc finger", tolower(paste(gene_symbol, description))) ~ "Zinc_Finger",
    grepl("homeobox|homeodomain|hox", tolower(paste(gene_symbol, description))) ~ "Homeodomain",
    grepl("bhlh", tolower(paste(gene_symbol, description))) ~ "bHLH",
    grepl("forkhead|fox", tolower(paste(gene_symbol, description))) ~ "Forkhead",
    grepl("hmg|sox", tolower(paste(gene_symbol, description))) ~ "HMG_Sox",
    TRUE ~ "Other_TF"
  ))

for (sc in unique(tf_subcats$subcat)) {
  sc_genes <- tf_subcats$gene_id[tf_subcats$subcat == sc]
  n_with <- sum(sc_genes %in% gene_dmr$nearest_gene)
  n_total <- length(sc_genes)
  if (n_total >= 5) {
    cat(sprintf("    %-15s: %d/%d = %.1f%% with DMRs\n",
                sc, n_with, n_total, n_with / n_total * 100))
  }
}


# ══════════════════════════════════════════════════════════════════════════════
# PART 8: DMR-Expression Integration
# ══════════════════════════════════════════════════════════════════════════════

cat("\n\n=== PART 8: DMR-Expression Integration ===\n")

# Genes with DMRs that are also DE
dmr_de <- gene_dmr %>%
  filter(!is.na(deg_padj) & deg_padj < 0.05)
cat("  Genes with DMR + DE:", nrow(dmr_de), "\n")

# DMR direction vs expression direction
if (nrow(dmr_de) > 0) {
  dmr_de$expr_dir <- ifelse(dmr_de$deg_log2FC > 0, "Up", "Down")
  concordance <- dmr_de %>%
    mutate(
      expected = case_when(
        dmr_direction == "Hypo" & expr_dir == "Up" ~ "Classical",
        dmr_direction == "Hyper" & expr_dir == "Down" ~ "Classical",
        TRUE ~ "Non-classical"
      )
    )
  cat(sprintf("  Classical (hypo→up, hyper→down): %d/%d = %.1f%%\n",
              sum(concordance$expected == "Classical"), nrow(concordance),
              mean(concordance$expected == "Classical") * 100))

  if (nrow(concordance) >= 5) {
    bt <- binom.test(sum(concordance$expected == "Classical"), nrow(concordance), p = 0.5)
    cat(sprintf("  Binomial test p = %.4f\n", bt$p.value))
  }
}

# Region-stratified DMR-expression correlation
cat("\n  DMR-expression correlation by region:\n")
gene_dmr_expr <- gene_dmr %>% filter(!is.na(deg_log2FC))
for (reg in c("Promoter", "Gene Body", "Intergenic", "Intron")) {
  sub <- gene_dmr_expr %>% filter(primary_region == reg)
  if (nrow(sub) >= 10) {
    ct <- cor.test(sub$mean_dmr_meth, sub$deg_log2FC, method = "spearman")
    cat(sprintf("    %-12s n=%3d  rho=%.4f  p=%.4f\n",
                reg, nrow(sub), ct$estimate, ct$p.value))
  }
}


# ══════════════════════════════════════════════════════════════════════════════
# SUMMARY
# ══════════════════════════════════════════════════════════════════════════════

cat("\n\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("DMR DEEP ANALYSIS SUMMARY\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat("KEY FINDINGS:\n\n")

cat("1. DMR Category Enrichment (vs DMP):\n")
sig_dmr <- dmr_enrichment %>% filter(padj < 0.05)
if (nrow(sig_dmr) > 0) {
  for (i in 1:nrow(sig_dmr)) {
    cat(sprintf("   - %s: OR=%.3f, padj=%.4f (%s)\n",
                sig_dmr$category[i], sig_dmr$odds_ratio[i], sig_dmr$padj[i], sig_dmr$direction[i]))
  }
} else {
  cat("   No categories reach significance after correction\n")
}

cat("\n2. Per-Module TF Fisher's Tests:\n")
cat(sprintf("   Global TF DMP: OR=%.3f (p=%.2e)\n", ft_global_dmp$estimate, ft_global_dmp$p.value))
cat(sprintf("   Global TF DMR: OR=%.3f (p=%.2e)\n", ft_global_dmr$estimate, ft_global_dmr$p.value))
turq_row <- tf_fisher_results %>% filter(module == "turquoise")
if (nrow(turq_row) > 0) {
  cat(sprintf("   Turquoise TF DMP: OR=%.3f, padj=%.4f (TF%%=%.1f%% vs NonTF%%=%.1f%%)\n",
              turq_row$dmp_or, turq_row$dmp_padj, turq_row$tf_dmp_pct, turq_row$nontf_dmp_pct))
  cat(sprintf("   Turquoise TF DMR: OR=%.3f, padj=%.4f (TF%%=%.1f%% vs NonTF%%=%.1f%%)\n",
              turq_row$dmr_or, turq_row$dmr_padj, turq_row$tf_dmr_pct, turq_row$nontf_dmr_pct))
}

cat("\n3. TFs with DMRs:", nrow(tf_with_dmr), "\n")

cat("\n4. Large DMRs (nCG >= 20):", nrow(large_dmrs), "\n")
cat("   Hyper:", nrow(large_hyper), " Hypo:", nrow(large_hypo), "\n")
cat("   Largest:", large_dmrs$nCG[1], "CpGs near", large_dmrs$gene_symbol[1], "\n")

cat("\n5. DMR Module Enrichment:\n")
sig_mod <- dmr_module_detail %>% filter(padj < 0.05)
if (nrow(sig_mod) > 0) {
  for (i in 1:nrow(sig_mod)) {
    cat(sprintf("   - %s: OR=%.3f, padj=%.4f\n",
                sig_mod$module[i], sig_mod$odds_ratio[i], sig_mod$padj[i]))
  }
}

cat("\n\nOutput files:\n")
cat("  Figures:", length(list.files(fig_dir, pattern = "\\.png$")), "PNG +",
    length(list.files(fig_dir, pattern = "\\.pdf$")), "PDF\n")
cat("  Tables:", length(list.files(tab_dir, pattern = "\\.tsv$")), "TSV\n")

cat("\n=== Script 12 complete ===\n")
