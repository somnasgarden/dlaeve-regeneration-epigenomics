#!/usr/bin/env Rscript
# =============================================================================
# Script 08: Morphogenesis & Developmental Gene DMP Enrichment
# =============================================================================
# Purpose: Test whether DMPs are non-randomly enriched near morphogenesis/
#          developmental genes. The user observed that DMPs cluster near
#          morphogenesis-related genes — is this statistically significant?
#
# Methods:
#   - Map gene annotations to functional categories (morphogenesis, signaling,
#     TFs, metabolism, structural, etc.)
#   - Fisher's exact tests for DMP enrichment per category
#   - Distance-to-TSS analysis for developmental vs non-developmental DMPs
#   - Module × functional category cross-tabulation
#   - Pathway-level DMP burden analysis
#
# Input:
#   results/03_integration/Tables/MXT_gene_level_meth_vs_expression.txt
#   results/02_rnaseq/Part5_WGCNA/data/all_gene_module_assignments.tsv
#   results/02_rnaseq/Part5_WGCNA/data/enrichment_summary_all_modules.tsv
#   results/02_rnaseq/Part2_DEGs/data/DEresults_amputated_vs_control.csv
#   DATA/derLaeGenome_eviann_annotations.tsv
#   results/01_methylation/dmps_chipseeker_annotated.txt
#
# Output:  results/14_morphogenesis_enrichment/
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

cat("=== Script 08: Morphogenesis & Developmental Gene DMP Enrichment ===\n\n")

# Output directory
outdir <- "results/14_morphogenesis_enrichment"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
fig_dir <- file.path(outdir, "figures")
tab_dir <- file.path(outdir, "tables")
dir.create(fig_dir, showWarnings = FALSE)
dir.create(tab_dir, showWarnings = FALSE)

# ── Load data ────────────────────────────────────────────────────────────────

cat("Loading data...\n")

# Gene-level methylation
gene_meth <- read.delim("results/03_integration/Tables/MXT_gene_level_meth_vs_expression.txt",
                         stringsAsFactors = FALSE)
cat("  Gene-level methylation:", nrow(gene_meth), "genes with DMPs\n")

# Module assignments
modules <- read.delim("results/02_rnaseq/Part5_WGCNA/data/all_gene_module_assignments.tsv",
                       stringsAsFactors = FALSE)
cat("  Module assignments:", nrow(modules), "genes\n")

# DEG results
deg <- read.csv("results/02_rnaseq/Part2_DEGs/data/DEresults_amputated_vs_control.csv",
                 stringsAsFactors = FALSE)
colnames(deg)[1] <- "gene_id"
cat("  DEG results:", nrow(deg), "genes\n")

# EviAnn annotations
annot <- read.delim("/mnt/c/Users/rafae/Projects/DATA/derLaeGenome_eviann_annotations.tsv",
                     header = FALSE, stringsAsFactors = FALSE)
colnames(annot) <- c("gene_id", "gene_symbol", "description")
# Keep unique gene-level annotations
annot <- annot %>% distinct(gene_id, .keep_all = TRUE)
cat("  Annotations:", nrow(annot), "unique genes\n")

# ChIPseeker DMPs
dmps <- read.delim("results/01_methylation/dmps_chipseeker_annotated.txt",
                    stringsAsFactors = FALSE)
cat("  ChIPseeker DMPs:", nrow(dmps), "positions\n")

# DMP-level annotation (for methylation_diff)
dmps_annot <- read.delim("results/01_methylation/dmps_annotated.txt",
                          stringsAsFactors = FALSE)
cat("  Annotated DMPs:", nrow(dmps_annot), "positions\n")

# ── Part 1: Classify genes into functional categories ────────────────────────

cat("\n--- Part 1: Functional gene classification ---\n")

# Define keyword-based classification (case-insensitive)
classify_gene <- function(symbol, description) {
  s <- tolower(paste(symbol, description, sep = " "))

  # Morphogenesis & development
  if (grepl("morphogen|wnt|notch|hedgehog|bmp|fgf|tgf|nodal|sox|pax|hox|homeobox|lhx|dlx|msx|otx|six|eya|dach|fork.?head|foxp|tbx|gata|hand|twist|snail|slug|dorsal|ventral|anterior|posterior|axis|pattern|segment|gastrulat|neurulat|blastula|embryo|limb|appendage|regenerat", s)) {
    return("Morphogenesis_Development")
  }

  # Cell cycle & proliferation
  if (grepl("cell.?cycle|cyclin|cdk|mitot|mitosi|kinetochore|spindle|centromere|checkpoint|proliferat|aurora|polo|plk|cdc|condensin|cohesin|separase|securin|anaphase|metaphase|telophase|cytokinesis", s)) {
    return("Cell_Cycle_Proliferation")
  }

  # DNA repair & damage response
  if (grepl("dna.?repair|dna.?damage|recombina|nucleotide.?excision|base.?excision|mismatch.?repair|double.?strand|rad5|brca|atm|atr|fanconi|homologous.?recomb|non.?homologous|xeroderma|exonuclease|endonuclease", s)) {
    return("DNA_Repair")
  }

  # Chromatin & epigenetic
  if (grepl("histone|chromatin|methyltransferase|demethylase|acetyltransferase|deacetylase|hdac|hat|polycomb|trithorax|swi.?snf|nucleosome|heterochromatin|euchromatin|hp1|cbx|ezh|suz|eed|ring1|bmi1|dnmt|tet|uhrf", s)) {
    return("Chromatin_Epigenetic")
  }

  # Signaling pathways
  if (grepl("kinase|phosphatase|receptor|signal|mapk|erk|jnk|jak|stat|pi3k|akt|mtor|ras|raf|mek|src|abl|egfr|vegf|pdgf|igf|insulin|calcium|calmodulin|camp|cgmp|adenylyl|guanylyl|plc|pkc|diacylglycerol|inositol", s)) {
    return("Signaling")
  }

  # Transcription factors (broader)
  if (grepl("transcription.?factor|zinc.?finger|zf\\-|znf|c2h2|cchh|bhlh|basic.?helix|leucine.?zipper|bzip|hmg|high.?mobility|ets|mads|myb|rel|nf.?kb|ap.?1|creb|nuclear.?receptor|hormone.?receptor", s)) {
    return("Transcription_Factor")
  }

  # ECM & structural
  if (grepl("collagen|fibronectin|laminin|integrin|cadherin|protocadherin|extracellular.?matrix|ecm|matrix.?metallo|mmp|timp|elastin|fibrillin|tenascin|vitronectin|proteoglycan|heparan|chondroit|hyaluron|glycosaminoglycan", s)) {
    return("ECM_Structural")
  }

  # Immune & defense
  if (grepl("immun|defensin|toll|tlr|interleukin|cytokine|chemokine|complement|lectin|antimicrobial|pathogen|innate|interferon|nod|nlr|rig|cgas|sting|inflamm", s)) {
    return("Immune_Defense")
  }

  # Metabolism
  if (grepl("dehydrogenase|oxidoreductase|transferase|hydrolase|lyase|isomerase|ligase|synthase|synthetase|reductase|oxidase|metabolism|metabol|glycolys|krebs|tca|electron.?transport|atp.?synth|cytochrome|nadh|fadh|sult|sulfotransferase|sulfat", s)) {
    return("Metabolism")
  }

  # Transport
  if (grepl("transport|channel|pump|slc|abc.?transporter|aquaporin|ion.?channel|potassium|sodium|chloride|calcium.?channel|trp|voltage.?gated|glutamate.?receptor|gaba.?receptor|acetylcholine.?receptor|serotonin.?receptor", s)) {
    return("Transport_Channel")
  }

  # RNA processing
  if (grepl("splicing|spliceosome|snrnp|rna.?helicase|rna.?polymerase|ribosom|translation|elongation.?factor|initiation.?factor|trna|rrna|mrna|poly.?a|deadenylase|decapping", s)) {
    return("RNA_Processing")
  }

  # Cytoskeleton & motility
  if (grepl("actin|myosin|tubulin|kinesin|dynein|microtubul|microfilament|intermediate.?filament|tropomyosin|cofilin|profilin|arp2|wasp|wave|rho|rac|cdc42|lamellipodia|filopodia|focal.?adhesion", s)) {
    return("Cytoskeleton_Motility")
  }

  # Apoptosis
  if (grepl("apoptos|caspase|bcl|bax|bak|bid|bim|mcl|xiap|diablo|smac|cytochrome.?c|death.?domain|fadd|tradd|trail|fas.?ligand|programmed.?cell.?death", s)) {
    return("Apoptosis")
  }

  return("Other_Unknown")
}

# Classify all annotated genes
annot$category <- mapply(classify_gene, annot$gene_symbol, annot$description)

cat("Gene functional categories:\n")
cat_counts <- sort(table(annot$category), decreasing = TRUE)
for (cat_name in names(cat_counts)) {
  cat(sprintf("  %-30s %5d genes\n", cat_name, cat_counts[cat_name]))
}

# ── Part 2: DMP enrichment by functional category ───────────────────────────

cat("\n--- Part 2: Fisher's exact test for DMP enrichment per category ---\n")

# Build gene universe: all genes in DEG results (expressed)
universe <- deg$gene_id

# Genes with DMPs
genes_with_dmps <- unique(gene_meth$gene_id)

# Merge category info
gene_data <- data.frame(gene_id = universe, stringsAsFactors = FALSE) %>%
  left_join(annot %>% select(gene_id, gene_symbol, category), by = "gene_id") %>%
  mutate(
    has_dmp = gene_id %in% genes_with_dmps,
    category = ifelse(is.na(category), "Unannotated", category)
  )

cat("Universe:", length(universe), "expressed genes\n")
cat("With DMPs:", sum(gene_data$has_dmp), "genes\n")

# Fisher's exact test per category
categories <- unique(gene_data$category)
categories <- categories[categories != "Unannotated"]

fisher_results <- data.frame()
for (cat_name in categories) {
  in_cat <- gene_data$category == cat_name
  # 2x2 table: category vs DMP
  a <- sum(in_cat & gene_data$has_dmp)   # in category, has DMP
  b <- sum(in_cat & !gene_data$has_dmp)  # in category, no DMP
  c <- sum(!in_cat & gene_data$has_dmp)  # not in category, has DMP
  d <- sum(!in_cat & !gene_data$has_dmp) # not in category, no DMP

  mat <- matrix(c(a, c, b, d), nrow = 2)
  ft <- fisher.test(mat)

  pct_in <- ifelse(a + b > 0, a / (a + b) * 100, 0)
  pct_out <- ifelse(c + d > 0, c / (c + d) * 100, 0)

  fisher_results <- rbind(fisher_results, data.frame(
    category = cat_name,
    n_genes = a + b,
    n_with_dmp = a,
    pct_with_dmp = round(pct_in, 1),
    background_pct = round(pct_out, 1),
    odds_ratio = round(ft$estimate, 3),
    pvalue = ft$p.value,
    direction = ifelse(ft$estimate > 1, "Enriched", "Depleted"),
    stringsAsFactors = FALSE
  ))
}

fisher_results$padj <- p.adjust(fisher_results$pvalue, method = "BH")
fisher_results <- fisher_results %>% arrange(pvalue)

cat("\nFisher's exact test results (DMP enrichment by functional category):\n")
for (i in 1:nrow(fisher_results)) {
  r <- fisher_results[i, ]
  sig <- ifelse(r$padj < 0.05, "***", ifelse(r$padj < 0.1, "*", ""))
  cat(sprintf("  %-30s OR=%.3f  padj=%.4f  %s/%s (%.1f%% vs %.1f%%)  %s %s\n",
              r$category, r$odds_ratio, r$padj,
              r$n_with_dmp, r$n_genes, r$pct_with_dmp, r$background_pct,
              r$direction, sig))
}

write.table(fisher_results, file.path(tab_dir, "I01_category_dmp_enrichment.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ── Part 3: Focus on Morphogenesis_Development category ──────────────────────

cat("\n--- Part 3: Morphogenesis/Development deep dive ---\n")

morph_genes <- gene_data %>% filter(category == "Morphogenesis_Development")
cat("Morphogenesis/Development genes:", nrow(morph_genes), "\n")
cat("  With DMPs:", sum(morph_genes$has_dmp), "(", round(mean(morph_genes$has_dmp) * 100, 1), "%)\n")

# Get the actual DMP details for morphogenesis genes
morph_with_dmps <- morph_genes %>%
  filter(has_dmp) %>%
  left_join(gene_meth, by = "gene_id")

if (nrow(morph_with_dmps) > 0) {
  cat("  Mean DMPs per morph gene:", round(mean(morph_with_dmps$n_dmps, na.rm = TRUE), 1), "\n")
  cat("  Mean |meth diff|:", round(mean(abs(morph_with_dmps$mean_meth_diff), na.rm = TRUE), 3), "\n")

  # Add DEG info
  morph_with_dmps <- morph_with_dmps %>%
    left_join(deg %>% select(gene_id, log2FoldChange, padj) %>% rename(deg_padj = padj),
              by = "gene_id")

  # DE significance
  morph_de <- morph_with_dmps %>% filter(!is.na(deg_padj) & deg_padj < 0.05)
  cat("  Morph genes with DMPs AND differentially expressed:", nrow(morph_de), "\n")

  if (nrow(morph_de) > 0) {
    cat("\n  Key morphogenesis genes (DMP + DE):\n")
    morph_de_show <- morph_de %>%
      arrange(deg_padj) %>%
      head(20)
    for (i in 1:nrow(morph_de_show)) {
      r <- morph_de_show[i, ]
      cat(sprintf("    %s (%s): %d DMPs, meth=%.3f, log2FC=%.2f, padj=%.4f, region=%s\n",
                  r$gene_id, r$gene_symbol, r$n_dmps, r$mean_meth_diff,
                  r$log2FoldChange.x, r$deg_padj, r$primary_region))
    }
  }
}

# ── Part 4: DMP burden comparison across categories ──────────────────────────

cat("\n--- Part 4: DMP burden per gene by category ---\n")

# For genes with DMPs, compare DMP count and effect size by category
burden <- gene_meth %>%
  left_join(annot %>% select(gene_id, category), by = "gene_id") %>%
  mutate(category = ifelse(is.na(category), "Unannotated", category)) %>%
  filter(category != "Unannotated")

burden_summary <- burden %>%
  group_by(category) %>%
  summarise(
    n_genes = n(),
    mean_dmps = mean(n_dmps),
    median_dmps = median(n_dmps),
    mean_abs_meth = mean(abs(mean_meth_diff)),
    max_dmps = max(n_dmps),
    pct_hyper = mean(meth_direction == "Hyper", na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  arrange(desc(mean_dmps))

cat("\nDMP burden by category (genes with DMPs only):\n")
for (i in 1:nrow(burden_summary)) {
  r <- burden_summary[i, ]
  cat(sprintf("  %-30s n=%4d  mean_dmps=%.1f  median=%.0f  max=%.0f  |meth|=%.3f  %%hyper=%.0f\n",
              r$category, r$n_genes, r$mean_dmps, r$median_dmps, r$max_dmps,
              r$mean_abs_meth, r$pct_hyper))
}

# Kruskal-Wallis test for DMP count differences
kw_dmps <- kruskal.test(n_dmps ~ category, data = burden)
cat("\nKruskal-Wallis test (DMP count ~ category): chi^2=",
    round(kw_dmps$statistic, 2), ", p=", format(kw_dmps$p.value, digits = 4), "\n")

# Kruskal-Wallis for methylation effect size
kw_meth <- kruskal.test(abs(mean_meth_diff) ~ category, data = burden)
cat("Kruskal-Wallis test (|meth diff| ~ category): chi^2=",
    round(kw_meth$statistic, 2), ", p=", format(kw_meth$p.value, digits = 4), "\n")

write.table(burden_summary, file.path(tab_dir, "I02_dmp_burden_by_category.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ── Part 5: Morphogenesis genes in WGCNA modules ────────────────────────────

cat("\n--- Part 5: Morphogenesis gene distribution across WGCNA modules ---\n")

morph_modules <- gene_data %>%
  filter(category == "Morphogenesis_Development") %>%
  left_join(modules %>% select(Gene, Module_Color) %>% rename(gene_id = Gene), by = "gene_id")

morph_mod_table <- morph_modules %>%
  filter(!is.na(Module_Color)) %>%
  group_by(Module_Color) %>%
  summarise(
    n_morph = n(),
    n_morph_dmp = sum(has_dmp),
    pct_with_dmp = round(n_morph_dmp / n_morph * 100, 1),
    .groups = "drop"
  ) %>%
  arrange(desc(n_morph))

cat("\nMorphogenesis genes by WGCNA module:\n")
for (i in 1:nrow(morph_mod_table)) {
  r <- morph_mod_table[i, ]
  cat(sprintf("  %-12s %4d morph genes, %3d with DMPs (%.1f%%)\n",
              r$Module_Color, r$n_morph, r$n_morph_dmp, r$pct_with_dmp))
}

# Which modules have more morphogenesis genes than expected?
total_genes_per_module <- modules %>%
  group_by(Module_Color) %>%
  summarise(total = n(), .groups = "drop")

morph_mod_enrich <- morph_mod_table %>%
  left_join(total_genes_per_module, by = "Module_Color") %>%
  mutate(
    pct_morph = n_morph / total * 100,
    expected_pct = nrow(morph_modules) / nrow(modules) * 100
  )

cat("\nModule enrichment for morphogenesis genes (expected:",
    round(nrow(morph_modules) / nrow(modules) * 100, 1), "%):\n")
for (i in 1:nrow(morph_mod_enrich)) {
  r <- morph_mod_enrich[i, ]
  enriched <- ifelse(r$pct_morph > r$expected_pct * 1.5, "ENRICHED",
                     ifelse(r$pct_morph < r$expected_pct * 0.5, "DEPLETED", ""))
  cat(sprintf("  %-12s %.1f%% morph genes (%d/%d) %s\n",
              r$Module_Color, r$pct_morph, r$n_morph, r$total, enriched))
}

write.table(morph_mod_enrich, file.path(tab_dir, "I03_morphogenesis_module_distribution.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ── Part 6: Distance to TSS for developmental vs other genes ────────────────

cat("\n--- Part 6: DMP distance-to-TSS by gene category ---\n")

# Merge DMPs with gene categories
dmp_with_cat <- dmps %>%
  left_join(annot %>% select(gene_id, category), by = c("geneId" = "gene_id")) %>%
  mutate(category = ifelse(is.na(category), "Unannotated", category))

# Compare distanceToTSS for morph vs others
morph_dmps <- dmp_with_cat %>% filter(category == "Morphogenesis_Development")
other_dmps <- dmp_with_cat %>% filter(category != "Morphogenesis_Development" & category != "Unannotated")

cat("DMPs near morphogenesis genes:", nrow(morph_dmps), "\n")
cat("DMPs near other annotated genes:", nrow(other_dmps), "\n")

if (nrow(morph_dmps) > 0 & nrow(other_dmps) > 0) {
  cat("Morphogenesis DMPs - median distance to TSS:", median(abs(morph_dmps$distanceToTSS)), "\n")
  cat("Other DMPs - median distance to TSS:", median(abs(other_dmps$distanceToTSS)), "\n")

  wt <- wilcox.test(abs(morph_dmps$distanceToTSS), abs(other_dmps$distanceToTSS))
  cat("Wilcoxon rank-sum test: W=", wt$statistic, ", p=", format(wt$p.value, digits = 4), "\n")
}

# ── Part 7: Specific developmental pathway analysis ─────────────────────────

cat("\n--- Part 7: Specific developmental pathway DMP analysis ---\n")

# Define key developmental pathways
pathways <- list(
  Wnt = c("wnt", "frizzled", "dishevelled", "beta-catenin", "ctnnb", "axin", "apc", "gsk3", "lef", "tcf", "sfrp", "dkk", "rspo"),
  Notch = c("notch", "delta", "jagged", "serrate", "hes", "hey", "rbpj", "maml", "numb", "mib"),
  Hedgehog = c("hedgehog", "patched", "smoothened", "gli", "sufu", "ptch", "smo"),
  BMP_TGFb = c("bmp", "tgf", "smad", "noggin", "chordin", "follistatin", "activin", "inhibin", "gdf"),
  FGF = c("fgf", "fgfr"),
  Sox = c("sox", "sry"),
  Hox = c("hox"),
  Pax = c("pax"),
  EGF_Ras = c("egf", "egfr", "erbb", "ras", "raf", "mek", "erk", "mapk")
)

pathway_results <- data.frame()
for (pw_name in names(pathways)) {
  kw <- pathways[[pw_name]]
  pattern <- paste(kw, collapse = "|")

  # Find genes in this pathway
  pw_genes <- annot %>%
    filter(grepl(pattern, tolower(paste(gene_symbol, description))))

  n_total <- nrow(pw_genes)
  n_with_dmp <- sum(pw_genes$gene_id %in% genes_with_dmps)

  # Get DMP details
  pw_dmps <- gene_meth %>% filter(gene_id %in% pw_genes$gene_id)

  # DE info
  pw_de <- deg %>% filter(gene_id %in% pw_genes$gene_id & !is.na(padj) & padj < 0.05)

  pathway_results <- rbind(pathway_results, data.frame(
    pathway = pw_name,
    n_genes = n_total,
    n_with_dmp = n_with_dmp,
    pct_with_dmp = ifelse(n_total > 0, round(n_with_dmp / n_total * 100, 1), 0),
    mean_dmps = ifelse(nrow(pw_dmps) > 0, round(mean(pw_dmps$n_dmps), 1), 0),
    n_de = nrow(pw_de),
    n_both = sum(pw_de$gene_id %in% genes_with_dmps),
    stringsAsFactors = FALSE
  ))
}

cat("\nDevelopmental pathway DMP summary:\n")
bg_pct <- round(length(genes_with_dmps) / length(universe) * 100, 1)
cat("(Background rate: ", bg_pct, "% of expressed genes have DMPs)\n\n")

for (i in 1:nrow(pathway_results)) {
  r <- pathway_results[i, ]
  status <- ifelse(r$pct_with_dmp > bg_pct * 1.5, "ENRICHED",
                   ifelse(r$pct_with_dmp < bg_pct * 0.5, "DEPLETED", "similar"))
  cat(sprintf("  %-15s %3d genes, %3d with DMP (%.1f%%), %d DE, %d both  [%s]\n",
              r$pathway, r$n_genes, r$n_with_dmp, r$pct_with_dmp,
              r$n_de, r$n_both, status))
}

write.table(pathway_results, file.path(tab_dir, "I04_developmental_pathway_dmps.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ── Part 8: Detailed gene lists for key pathways ────────────────────────────

cat("\n--- Part 8: Key pathway gene details ---\n")

for (pw_name in c("Wnt", "Notch", "Sox", "Hox", "BMP_TGFb")) {
  kw <- pathways[[pw_name]]
  pattern <- paste(kw, collapse = "|")

  pw_genes <- annot %>%
    filter(grepl(pattern, tolower(paste(gene_symbol, description)))) %>%
    left_join(gene_meth %>% select(gene_id, n_dmps, mean_meth_diff, primary_region), by = "gene_id") %>%
    left_join(deg %>% select(gene_id, log2FoldChange, padj), by = "gene_id") %>%
    mutate(
      has_dmp = !is.na(n_dmps),
      is_de = !is.na(padj) & padj < 0.05
    ) %>%
    arrange(desc(has_dmp), desc(is_de))

  cat(sprintf("\n  %s pathway (%d genes):\n", pw_name, nrow(pw_genes)))
  show_n <- min(15, nrow(pw_genes))
  for (i in 1:show_n) {
    r <- pw_genes[i, ]
    dmp_info <- ifelse(r$has_dmp, sprintf("%d DMPs (%.3f)", r$n_dmps, r$mean_meth_diff), "no DMPs")
    de_info <- ifelse(r$is_de, sprintf("DE log2FC=%.2f", r$log2FoldChange), "not DE")
    cat(sprintf("    %s (%s): %s, %s\n", r$gene_id, r$gene_symbol, dmp_info, de_info))
  }
}

# ── Part 9: Category × Module × DMP cross-tabulation ────────────────────────

cat("\n--- Part 9: Category x Module x DMP cross-tabulation ---\n")

cross_tab <- gene_data %>%
  filter(category != "Unannotated") %>%
  left_join(modules %>% select(Gene, Module_Color) %>% rename(gene_id = Gene), by = "gene_id") %>%
  filter(!is.na(Module_Color)) %>%
  group_by(category, Module_Color) %>%
  summarise(
    n_genes = n(),
    n_dmp = sum(has_dmp),
    pct_dmp = round(n_dmp / n_genes * 100, 1),
    .groups = "drop"
  ) %>%
  filter(n_genes >= 10)  # Only categories with enough genes

# Find the most interesting combinations
cross_tab_sig <- cross_tab %>%
  filter(pct_dmp >= 25 | (n_dmp >= 5 & pct_dmp >= 20)) %>%
  arrange(desc(pct_dmp))

cat("\nTop category × module combinations (high DMP rate, >=10 genes):\n")
for (i in 1:min(20, nrow(cross_tab_sig))) {
  r <- cross_tab_sig[i, ]
  cat(sprintf("  %-30s × %-12s %d/%d (%.1f%%)\n",
              r$category, r$Module_Color, r$n_dmp, r$n_genes, r$pct_dmp))
}

write.table(cross_tab, file.path(tab_dir, "I05_category_module_crosstab.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ── Part 10: Morphogenesis DMP genomic region preference ────────────────────

cat("\n--- Part 10: Morphogenesis DMP genomic region distribution ---\n")

if (nrow(morph_dmps) > 0) {
  # Parse ChIPseeker annotation into simple categories
  morph_dmps$region <- case_when(
    grepl("Promoter", morph_dmps$annotation) ~ "Promoter",
    grepl("Intron", morph_dmps$annotation) ~ "Intron",
    grepl("Exon", morph_dmps$annotation) ~ "Exon",
    grepl("UTR", morph_dmps$annotation) ~ "UTR",
    grepl("Downstream", morph_dmps$annotation) ~ "Downstream",
    grepl("Intergenic", morph_dmps$annotation) ~ "Intergenic",
    TRUE ~ "Other"
  )

  morph_region <- table(morph_dmps$region)
  other_region_dmps <- dmp_with_cat %>%
    filter(category != "Morphogenesis_Development" & category != "Unannotated") %>%
    mutate(region = case_when(
      grepl("Promoter", annotation) ~ "Promoter",
      grepl("Intron", annotation) ~ "Intron",
      grepl("Exon", annotation) ~ "Exon",
      grepl("UTR", annotation) ~ "UTR",
      grepl("Downstream", annotation) ~ "Downstream",
      grepl("Intergenic", annotation) ~ "Intergenic",
      TRUE ~ "Other"
    ))
  other_region <- table(other_region_dmps$region)

  region_compare <- data.frame(
    region = names(morph_region),
    morph_n = as.integer(morph_region),
    morph_pct = round(as.integer(morph_region) / sum(morph_region) * 100, 1),
    stringsAsFactors = FALSE
  ) %>%
    left_join(
      data.frame(
        region = names(other_region),
        other_n = as.integer(other_region),
        other_pct = round(as.integer(other_region) / sum(other_region) * 100, 1),
        stringsAsFactors = FALSE
      ),
      by = "region"
    ) %>%
    arrange(desc(morph_pct))

  cat("\nGenomic region distribution - Morphogenesis vs Other:\n")
  for (i in 1:nrow(region_compare)) {
    r <- region_compare[i, ]
    diff <- r$morph_pct - r$other_pct
    cat(sprintf("  %-12s Morph: %.1f%% (%d)   Other: %.1f%% (%d)   diff: %+.1f%%\n",
                r$region, r$morph_pct, r$morph_n, r$other_pct, r$other_n, diff))
  }

  # DMP direction for morphogenesis genes
  morph_dmps_dir <- dmps_annot %>%
    filter(nearest_gene %in% morph_genes$gene_id[morph_genes$has_dmp])

  if (nrow(morph_dmps_dir) > 0) {
    n_hyper <- sum(morph_dmps_dir$methylation_diff > 0)
    n_hypo <- sum(morph_dmps_dir$methylation_diff < 0)
    cat("\nMorphogenesis DMPs direction: ", n_hyper, " hyper (", round(n_hyper/(n_hyper+n_hypo)*100,1),
        "%), ", n_hypo, " hypo (", round(n_hypo/(n_hyper+n_hypo)*100,1), "%)\n")

    bt <- binom.test(n_hyper, n_hyper + n_hypo, p = 0.5)
    cat("Binomial test vs 50/50: p=", format(bt$p.value, digits = 4), "\n")
  }
}

# ══════════════════════════════════════════════════════════════════════════════
# FIGURES
# ══════════════════════════════════════════════════════════════════════════════

cat("\n--- Generating figures ---\n")

# ── Figure I01: Category DMP enrichment forest plot ──────────────────────────

fig_data <- fisher_results %>%
  filter(n_genes >= 20) %>%
  mutate(
    category = reorder(category, odds_ratio),
    sig = ifelse(padj < 0.05, "padj < 0.05", ifelse(padj < 0.1, "padj < 0.1", "NS"))
  )

p1 <- ggplot(fig_data, aes(x = odds_ratio, y = category)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  geom_point(aes(size = n_genes, color = sig)) +
  geom_text(aes(label = sprintf("%.1f%%", pct_with_dmp)), hjust = -0.3, size = 3) +
  scale_color_manual(values = c("padj < 0.05" = "firebrick", "padj < 0.1" = "darkorange", "NS" = "gray40")) +
  scale_size_continuous(range = c(2, 8), name = "N genes") +
  labs(
    title = "DMP Enrichment by Functional Gene Category",
    subtitle = sprintf("Background rate: %.1f%% of expressed genes have DMPs", bg_pct),
    x = "Odds Ratio (Fisher's exact test)",
    y = "",
    color = "Significance"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(fig_dir, "I01_category_dmp_enrichment.png"), p1,
       width = 10, height = 7, dpi = 300)
ggsave(file.path(fig_dir, "I01_category_dmp_enrichment.pdf"), p1,
       width = 10, height = 7)
cat("  Saved I01_category_dmp_enrichment\n")

# ── Figure I02: DMP burden boxplot by category ───────────────────────────────

burden_plot <- burden %>%
  filter(category %in% c("Morphogenesis_Development", "Cell_Cycle_Proliferation",
                          "Signaling", "Transcription_Factor", "ECM_Structural",
                          "Metabolism", "Chromatin_Epigenetic", "DNA_Repair",
                          "Cytoskeleton_Motility", "Transport_Channel"))

p2 <- ggplot(burden_plot, aes(x = reorder(category, n_dmps, FUN = median), y = n_dmps)) +
  geom_boxplot(aes(fill = category), show.legend = FALSE, outlier.alpha = 0.3) +
  coord_flip() +
  scale_y_log10() +
  labs(
    title = "DMP Burden per Gene by Functional Category",
    subtitle = "Among genes that have at least one DMP",
    x = "",
    y = "Number of DMPs per gene (log10 scale)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(fig_dir, "I02_dmp_burden_by_category.png"), p2,
       width = 10, height = 7, dpi = 300)
ggsave(file.path(fig_dir, "I02_dmp_burden_by_category.pdf"), p2,
       width = 10, height = 7)
cat("  Saved I02_dmp_burden_by_category\n")

# ── Figure I03: Developmental pathway heatmap ────────────────────────────────

pw_heat_data <- pathway_results %>%
  mutate(
    enrichment = pct_with_dmp / bg_pct,
    label = sprintf("%.0f%%\n(%d/%d)", pct_with_dmp, n_with_dmp, n_genes)
  )

p3 <- ggplot(pw_heat_data, aes(x = 1, y = reorder(pathway, pct_with_dmp))) +
  geom_tile(aes(fill = pct_with_dmp), width = 0.8) +
  geom_text(aes(label = label), size = 3.5) +
  geom_hline(yintercept = 0.5, color = "gray") +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                       midpoint = bg_pct, name = "% with DMP") +
  labs(
    title = "DMP Rate in Developmental Signaling Pathways",
    subtitle = sprintf("Background: %.1f%% | Red = enriched, Blue = depleted", bg_pct),
    x = "", y = ""
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank()
  )

ggsave(file.path(fig_dir, "I03_pathway_dmp_heatmap.png"), p3,
       width = 6, height = 7, dpi = 300)
ggsave(file.path(fig_dir, "I03_pathway_dmp_heatmap.pdf"), p3,
       width = 6, height = 7)
cat("  Saved I03_pathway_dmp_heatmap\n")

# ── Figure I04: Category × Module heatmap ────────────────────────────────────

heat_data <- cross_tab %>%
  filter(n_genes >= 10) %>%
  select(category, Module_Color, pct_dmp) %>%
  pivot_wider(names_from = Module_Color, values_from = pct_dmp, values_fill = 0)

# Convert to matrix for pheatmap
heat_mat <- as.matrix(heat_data[, -1])
rownames(heat_mat) <- heat_data$category

# Only keep modules with enough data
col_sums <- colSums(heat_mat > 0)
heat_mat <- heat_mat[, col_sums >= 3, drop = FALSE]

if (ncol(heat_mat) >= 2 & nrow(heat_mat) >= 3) {
  library(pheatmap)

  # Module colors for annotation
  mod_colors <- colnames(heat_mat)
  names(mod_colors) <- colnames(heat_mat)

  png(file.path(fig_dir, "I04_category_module_heatmap.png"), width = 12, height = 8, units = "in", res = 300)
  pheatmap(heat_mat,
           color = colorRampPalette(c("white", "orange", "firebrick"))(50),
           main = "DMP Rate (%) by Functional Category × WGCNA Module",
           fontsize = 10,
           fontsize_row = 9,
           fontsize_col = 10,
           cluster_cols = TRUE,
           cluster_rows = TRUE,
           display_numbers = TRUE,
           number_format = "%.0f",
           number_color = "gray20",
           border_color = "gray90")
  dev.off()

  pdf(file.path(fig_dir, "I04_category_module_heatmap.pdf"), width = 12, height = 8)
  pheatmap(heat_mat,
           color = colorRampPalette(c("white", "orange", "firebrick"))(50),
           main = "DMP Rate (%) by Functional Category × WGCNA Module",
           fontsize = 10,
           fontsize_row = 9,
           fontsize_col = 10,
           cluster_cols = TRUE,
           cluster_rows = TRUE,
           display_numbers = TRUE,
           number_format = "%.0f",
           number_color = "gray20",
           border_color = "gray90")
  dev.off()
  cat("  Saved I04_category_module_heatmap\n")
}

# ── Figure I05: Morphogenesis region distribution ────────────────────────────

if (exists("region_compare") && nrow(region_compare) > 0) {
  region_long <- region_compare %>%
    select(region, morph_pct, other_pct) %>%
    pivot_longer(cols = c(morph_pct, other_pct),
                 names_to = "group", values_to = "percentage") %>%
    mutate(group = ifelse(group == "morph_pct", "Morphogenesis", "Other genes"))

  p5 <- ggplot(region_long, aes(x = reorder(region, -percentage), y = percentage, fill = group)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    scale_fill_manual(values = c("Morphogenesis" = "#E74C3C", "Other genes" = "#3498DB")) +
    labs(
      title = "Genomic Region Distribution of DMPs",
      subtitle = "Morphogenesis/Development genes vs Other annotated genes",
      x = "Genomic Region",
      y = "Percentage of DMPs (%)",
      fill = ""
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "top",
      panel.grid.minor = element_blank()
    )

  ggsave(file.path(fig_dir, "I05_morphogenesis_region_distribution.png"), p5,
         width = 9, height = 6, dpi = 300)
  ggsave(file.path(fig_dir, "I05_morphogenesis_region_distribution.pdf"), p5,
         width = 9, height = 6)
  cat("  Saved I05_morphogenesis_region_distribution\n")
}

# ══════════════════════════════════════════════════════════════════════════════
# SUMMARY
# ══════════════════════════════════════════════════════════════════════════════

cat("\n" , paste(rep("=", 70), collapse = ""), "\n")
cat("SUMMARY: Morphogenesis & Developmental Gene DMP Enrichment\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# Find the morphogenesis result
morph_row <- fisher_results %>% filter(category == "Morphogenesis_Development")
if (nrow(morph_row) > 0) {
  cat("KEY RESULT: Morphogenesis/Development genes\n")
  cat("  DMP rate:", morph_row$pct_with_dmp, "% (vs", morph_row$background_pct, "% background)\n")
  cat("  Odds ratio:", morph_row$odds_ratio, "\n")
  cat("  Fisher's p-value:", format(morph_row$pvalue, digits = 4), "\n")
  cat("  Adjusted p-value:", format(morph_row$padj, digits = 4), "\n")
  cat("  Direction:", morph_row$direction, "\n\n")
}

# Cell cycle
cc_row <- fisher_results %>% filter(category == "Cell_Cycle_Proliferation")
if (nrow(cc_row) > 0) {
  cat("Cell Cycle/Proliferation genes:\n")
  cat("  DMP rate:", cc_row$pct_with_dmp, "% vs", cc_row$background_pct, "% background\n")
  cat("  OR=", cc_row$odds_ratio, ", padj=", format(cc_row$padj, digits = 4), "\n\n")
}

# TFs
tf_row <- fisher_results %>% filter(category == "Transcription_Factor")
if (nrow(tf_row) > 0) {
  cat("Transcription Factor genes:\n")
  cat("  DMP rate:", tf_row$pct_with_dmp, "% vs", tf_row$background_pct, "% background\n")
  cat("  OR=", tf_row$odds_ratio, ", padj=", format(tf_row$padj, digits = 4), "\n\n")
}

cat("\nSignificantly enriched categories (padj < 0.05):\n")
sig_enriched <- fisher_results %>% filter(padj < 0.05 & direction == "Enriched")
if (nrow(sig_enriched) > 0) {
  for (i in 1:nrow(sig_enriched)) {
    cat("  +", sig_enriched$category[i], "OR=", sig_enriched$odds_ratio[i], "\n")
  }
} else {
  cat("  None\n")
}

cat("\nSignificantly depleted categories (padj < 0.05):\n")
sig_depleted <- fisher_results %>% filter(padj < 0.05 & direction == "Depleted")
if (nrow(sig_depleted) > 0) {
  for (i in 1:nrow(sig_depleted)) {
    cat("  -", sig_depleted$category[i], "OR=", sig_depleted$odds_ratio[i], "\n")
  }
} else {
  cat("  None\n")
}

cat("\n=== Script 08 complete ===\n")
cat("Output:", outdir, "\n")
cat("Figures:", length(list.files(fig_dir, pattern = "\\.png$")), "PNG +",
    length(list.files(fig_dir, pattern = "\\.pdf$")), "PDF\n")
cat("Tables:", length(list.files(tab_dir, pattern = "\\.tsv$")), "TSV\n")
