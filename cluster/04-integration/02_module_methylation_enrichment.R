#!/usr/bin/env Rscript
# =============================================================================
# MXT Part 2: WGCNA Modules × Methylation — Module & Hub Gene Overlap
# =============================================================================
# Goal: Determine which WGCNA co-expression modules are most affected by
#        differential methylation (DMPs/DMRs), whether hub genes are impacted,
#        and in which genomic regions the methylation changes occur.
#
# Inputs:
#   - WGCNA objects: module assignments, hub genes, eigengenes, Wilcoxon results
#   - Methylation: DMP/DMR results + ChIPseeker annotations
#   - MXT Part 1 outputs (gene-level methylation–expression table)
#
# Outputs:
#   - Module-level DMP/DMR overlap tables
#   - Hub gene methylation status
#   - Region breakdown per module
#   - Enrichment of DMPs in specific modules (Fisher's exact)
#   - Visualizations: barplots, heatmaps, network-style plots
# =============================================================================

# ── Libraries ────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(scales)
  library(ggrepel)
})

# Try loading WGCNA (needed for labels2colors)
if (requireNamespace("WGCNA", quietly = TRUE)) {
  library(WGCNA)
}

options(stringsAsFactors = FALSE)
options(scipen = 999)

theme_set(theme_minimal() +
            theme(text = element_text(size = 12),
                  plot.title = element_text(size = 14, face = "bold"),
                  axis.title = element_text(size = 12),
                  legend.text = element_text(size = 10)))

# ── Paths ────────────────────────────────────────────────────────────────────

METH_DIR    <- "/mnt/data/alfredvar/rlopezt/Preliminary/"
TRANS_DIR   <- "/mnt/data/alfredvar/rlopezt/CorrelationMatrix/"
MXT_DIR     <- "/mnt/data/alfredvar/rlopezt/MXT/"

# Output directories (reuse MXT structure)
dirs <- list(
  tables     = file.path(MXT_DIR, "Tables"),
  pdf        = file.path(MXT_DIR, "PDFs"),
  png        = file.path(MXT_DIR, "PNGs"),
  rds        = file.path(MXT_DIR, "RDS_Objects")
)
for (d in dirs) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# Helpers
save_both <- function(plot_obj, name, width, height) {
  ggsave(file.path(dirs$pdf, paste0(name, ".pdf")), plot_obj,
         width = width, height = height)
  ggsave(file.path(dirs$png, paste0(name, ".png")), plot_obj,
         width = width, height = height, dpi = 300)
  cat(sprintf("  Saved: %s (.pdf + .png)\n", name))
}

save_pheatmap_both <- function(pheatmap_expr, name, width, height) {
  pdf(file.path(dirs$pdf, paste0(name, ".pdf")), width = width, height = height)
  eval(pheatmap_expr)
  dev.off()
  png(file.path(dirs$png, paste0(name, ".png")),
      width = width, height = height, units = "in", res = 300)
  eval(pheatmap_expr)
  dev.off()
  cat(sprintf("  Saved: %s (.pdf + .png)\n", name))
}

cat("================================================================\n")
cat("  MXT Part 2: WGCNA Modules × Methylation\n")
cat("================================================================\n\n")


# #############################################################################
# SECTION 1: LOAD DATA
# #############################################################################

cat("=== Loading Data ===\n\n")

# 1.1 WGCNA objects
wgcna_file <- file.path(TRANS_DIR, "Part5_WGCNA/data/Part5_WGCNA_objects.RData")
cat(sprintf("  Loading WGCNA: %s\n", wgcna_file))
load(wgcna_file)

cat(sprintf("  WGCNA modules: %d (+ grey)\n", length(setdiff(unique(module_colors), "grey"))))
cat(sprintf("  Total genes in network: %s\n", format(length(module_colors), big.mark = ",")))
cat(sprintf("  Hub genes: %d | Candidate genes: %d\n",
            sum(all_hub_genes$IsHub), sum(all_hub_genes$IsCandidate)))
cat(sprintf("  Best Wilcoxon module: %s (p = %s)\n",
            best_module_wilcox, signif(me_wilcox$Wilcox_p[1], 3)))

# Module assignments
module_assignments <- data.frame(
  gene_id = names(module_colors),
  module  = module_colors,
  stringsAsFactors = FALSE
)

# Also load TSV version for backup
module_tsv <- file.path(TRANS_DIR, "Part5_WGCNA/data/all_gene_module_assignments.tsv")
if (file.exists(module_tsv)) {
  module_assign_tsv <- read.delim(module_tsv, stringsAsFactors = FALSE)
  cat(sprintf("  Module TSV: %s genes\n", nrow(module_assign_tsv)))
}

# Hub genes table
hub_tsv <- file.path(TRANS_DIR, "Part5_WGCNA/data/all_modules_hub_genes.tsv")
hub_genes_df <- if (file.exists(hub_tsv)) {
  read.delim(hub_tsv, stringsAsFactors = FALSE)
} else {
  all_hub_genes
}

# Wilcoxon ME results
wilcox_tsv <- file.path(TRANS_DIR, "Part5_WGCNA/data/ME_wilcoxon_tail_regeneration.tsv")
wilcox_results <- if (file.exists(wilcox_tsv)) {
  read.delim(wilcox_tsv, stringsAsFactors = FALSE)
} else {
  me_wilcox
}

cat(sprintf("  Hub genes data: %s entries\n", format(nrow(hub_genes_df), big.mark = ",")))

# 1.2 Methylation data
dmp_results <- readRDS(file.path(METH_DIR, "dmp_results.rds"))
dmr_results <- readRDS(file.path(METH_DIR, "dmr_results.rds"))

cat(sprintf("  DMPs: %s | DMRs: %s\n",
            format(nrow(dmp_results), big.mark = ","),
            format(nrow(dmr_results), big.mark = ",")))

# 1.3 ChIPseeker annotations
chipseeker_dir <- file.path(METH_DIR, "ChIPseeker")
dmp_cs_file <- file.path(chipseeker_dir, "dmps_chipseeker_annotated.txt")
dmr_cs_file <- file.path(chipseeker_dir, "dmrs_chipseeker_annotated.txt")

dmp_cs <- if (file.exists(dmp_cs_file)) read.delim(dmp_cs_file, stringsAsFactors = FALSE) else NULL
dmr_cs <- if (file.exists(dmr_cs_file)) read.delim(dmr_cs_file, stringsAsFactors = FALSE) else NULL

# Fallback custom annotations
dmp_anno_file <- file.path(METH_DIR, "dmps_annotated.txt")
if (!file.exists(dmp_anno_file)) dmp_anno_file <- file.path(METH_DIR, "Tables", "dmps_annotated.txt")
dmr_anno_file <- file.path(METH_DIR, "dmrs_annotated.txt")
if (!file.exists(dmr_anno_file)) dmr_anno_file <- file.path(METH_DIR, "Tables", "dmrs_annotated.txt")

dmp_anno <- if (file.exists(dmp_anno_file)) read.delim(dmp_anno_file, stringsAsFactors = FALSE) else NULL
dmr_anno <- if (file.exists(dmr_anno_file)) read.delim(dmr_anno_file, stringsAsFactors = FALSE) else NULL

# 1.4 MXT Part 1 objects (if available)
mxt1_file <- file.path(dirs$rds, "MXT_Part1_objects.RData")
has_mxt1 <- file.exists(mxt1_file)
if (has_mxt1) {
  load(mxt1_file)
  cat(sprintf("  MXT Part 1 data loaded: %s genes\n", nrow(mxt_data)))
}

# 1.5 Gene names
gene_info_file <- file.path(METH_DIR, "gene_info_with_names.txt")
if (!file.exists(gene_info_file))
  gene_info_file <- file.path(METH_DIR, "Tables", "gene_info_with_names.txt")
gene_names <- if (file.exists(gene_info_file)) {
  read.delim(gene_info_file, stringsAsFactors = FALSE)
} else NULL


# #############################################################################
# SECTION 2: BUILD DMP/DMR → GENE → MODULE MAPPING
# #############################################################################

cat("\n=== Building DMP/DMR → Gene → Module Mapping ===\n\n")

# 2.1 DMP → gene mapping (with region)
if (!is.null(dmp_cs)) {
  cat("  ChIPseeker DMP columns: ", paste(colnames(dmp_cs), collapse = ", "), "\n")
  
  # ChIPseeker file only has annotation — meth values from dmp_results
  if (nrow(dmp_cs) == nrow(dmp_results)) {
    dmp_gene_map <- data.frame(
      gene_id    = dmp_cs$geneId,
      meth_diff  = dmp_results$diff,
      fdr        = dmp_results$fdr,
      annotation = dmp_cs$annotation,
      stringsAsFactors = FALSE
    )
  } else {
    dmp_cs$site_id <- paste(dmp_cs$seqnames, dmp_cs$start, sep = ":")
    dmp_results$site_id <- paste(dmp_results$chr, dmp_results$pos, sep = ":")
    dmp_merged <- merge(dmp_cs, dmp_results, by = "site_id")
    dmp_gene_map <- data.frame(
      gene_id    = dmp_merged$geneId,
      meth_diff  = dmp_merged$diff,
      fdr        = dmp_merged$fdr,
      annotation = dmp_merged$annotation,
      stringsAsFactors = FALSE
    )
  }
} else if (!is.null(dmp_anno)) {
  dmp_gene_map <- data.frame(
    gene_id    = dmp_anno$nearest_gene,
    meth_diff  = dmp_anno$methylation_diff,
    fdr        = if ("fdr" %in% colnames(dmp_anno)) dmp_anno$fdr else NA,
    annotation = dmp_anno$annotation,
    stringsAsFactors = FALSE
  )
} else {
  stop("No DMP annotation found.")
}

# Simplify annotation to broad regions
dmp_gene_map$region <- dmp_gene_map$annotation
dmp_gene_map$region <- gsub(" \\(.*", "", dmp_gene_map$region)
dmp_gene_map$region[grepl("Promoter", dmp_gene_map$region)] <- "Promoter"
dmp_gene_map$region[grepl("Exon", dmp_gene_map$region)]     <- "Exon"
dmp_gene_map$region[grepl("Intron", dmp_gene_map$region)]   <- "Intron"
dmp_gene_map$region[grepl("UTR", dmp_gene_map$region)]      <- "UTR"
dmp_gene_map$region[grepl("Downstream", dmp_gene_map$region)] <- "Downstream"
dmp_gene_map$region[grepl("Intergenic|Distal", dmp_gene_map$region)] <- "Intergenic"

# Remove unmapped
dmp_gene_map <- dmp_gene_map[!is.na(dmp_gene_map$gene_id) & dmp_gene_map$gene_id != "", ]

# Add module assignment
dmp_gene_map <- merge(dmp_gene_map, module_assignments, by = "gene_id", all.x = TRUE)
dmp_gene_map$in_network <- !is.na(dmp_gene_map$module)

cat(sprintf("  DMPs with gene assignment: %s\n", format(nrow(dmp_gene_map), big.mark = ",")))
cat(sprintf("  DMPs in WGCNA network: %s (%.1f%%)\n",
            format(sum(dmp_gene_map$in_network), big.mark = ","),
            100 * mean(dmp_gene_map$in_network)))

# 2.2 DMR → gene mapping (with region)
if (!is.null(dmr_cs)) {
  if (nrow(dmr_cs) == nrow(dmr_results)) {
    dmr_gene_map <- data.frame(
      gene_id     = dmr_cs$geneId,
      dmr_diff    = dmr_results$diff.Methy,
      annotation  = dmr_cs$annotation,
      stringsAsFactors = FALSE
    )
  } else {
    dmr_cs$site_id <- paste(dmr_cs$seqnames, dmr_cs$start, dmr_cs$end, sep = ":")
    dmr_results$site_id <- paste(dmr_results$chr, dmr_results$start, dmr_results$end, sep = ":")
    dmr_merged <- merge(dmr_cs, dmr_results, by = "site_id")
    dmr_gene_map <- data.frame(
      gene_id     = dmr_merged$geneId,
      dmr_diff    = dmr_merged$diff.Methy,
      annotation  = dmr_merged$annotation,
      stringsAsFactors = FALSE
    )
  }
} else if (!is.null(dmr_anno)) {
  dmr_gene_map <- data.frame(
    gene_id    = dmr_anno$nearest_gene,
    dmr_diff   = if ("diff.Methy" %in% colnames(dmr_anno)) dmr_anno$diff.Methy else
      if ("methylation_diff" %in% colnames(dmr_anno)) dmr_anno$methylation_diff else NA,
    annotation = dmr_anno$annotation,
    stringsAsFactors = FALSE
  )
} else {
  dmr_gene_map <- data.frame()
}

if (nrow(dmr_gene_map) > 0) {
  dmr_gene_map$region <- dmr_gene_map$annotation
  dmr_gene_map$region <- gsub(" \\(.*", "", dmr_gene_map$region)
  dmr_gene_map$region[grepl("Promoter", dmr_gene_map$region)] <- "Promoter"
  dmr_gene_map$region[grepl("Exon", dmr_gene_map$region)]     <- "Exon"
  dmr_gene_map$region[grepl("Intron", dmr_gene_map$region)]   <- "Intron"
  dmr_gene_map$region[grepl("UTR", dmr_gene_map$region)]      <- "UTR"
  dmr_gene_map$region[grepl("Downstream", dmr_gene_map$region)] <- "Downstream"
  dmr_gene_map$region[grepl("Intergenic|Distal", dmr_gene_map$region)] <- "Intergenic"
  
  dmr_gene_map <- dmr_gene_map[!is.na(dmr_gene_map$gene_id) & dmr_gene_map$gene_id != "", ]
  dmr_gene_map <- merge(dmr_gene_map, module_assignments, by = "gene_id", all.x = TRUE)
  dmr_gene_map$in_network <- !is.na(dmr_gene_map$module)
  
  cat(sprintf("  DMRs with gene assignment: %s\n", nrow(dmr_gene_map)))
  cat(sprintf("  DMRs in WGCNA network: %s (%.1f%%)\n",
              sum(dmr_gene_map$in_network),
              100 * mean(dmr_gene_map$in_network)))
}


# #############################################################################
# SECTION 3: MODULE-LEVEL DMP/DMR ENRICHMENT
# #############################################################################

cat("\n=== Module-Level DMP/DMR Analysis ===\n\n")

real_modules <- setdiff(unique(module_assignments$module), "grey")

# 3.1 Count DMPs per module
module_dmp_counts <- dmp_gene_map %>%
  filter(in_network, module != "grey") %>%
  group_by(module) %>%
  summarise(
    n_dmps            = n(),
    n_genes_with_dmp  = n_distinct(gene_id),
    n_hyper           = sum(meth_diff > 0, na.rm = TRUE),
    n_hypo            = sum(meth_diff < 0, na.rm = TRUE),
    mean_meth_diff    = mean(meth_diff, na.rm = TRUE),
    .groups = "drop"
  )

# Add module size
module_sizes <- module_assignments %>%
  filter(module != "grey") %>%
  count(module, name = "module_size")

module_dmp_counts <- merge(module_dmp_counts, module_sizes, by = "module", all.y = TRUE)
module_dmp_counts$n_dmps[is.na(module_dmp_counts$n_dmps)] <- 0
module_dmp_counts$n_genes_with_dmp[is.na(module_dmp_counts$n_genes_with_dmp)] <- 0
module_dmp_counts$pct_genes_with_dmp <- 100 * module_dmp_counts$n_genes_with_dmp / module_dmp_counts$module_size

# 3.2 Count DMRs per module
if (nrow(dmr_gene_map) > 0) {
  module_dmr_counts <- dmr_gene_map %>%
    filter(in_network, module != "grey") %>%
    group_by(module) %>%
    summarise(
      n_dmrs            = n(),
      n_genes_with_dmr  = n_distinct(gene_id),
      mean_dmr_diff     = mean(dmr_diff, na.rm = TRUE),
      .groups = "drop"
    )
  
  module_dmp_counts <- merge(module_dmp_counts, module_dmr_counts, by = "module", all.x = TRUE)
  module_dmp_counts$n_dmrs[is.na(module_dmp_counts$n_dmrs)] <- 0
  module_dmp_counts$n_genes_with_dmr[is.na(module_dmp_counts$n_genes_with_dmr)] <- 0
}

# 3.3 Fisher's exact test: DMP enrichment per module
total_genes_network <- nrow(module_assignments[module_assignments$module != "grey", ])
dmp_genes_in_network <- unique(dmp_gene_map$gene_id[dmp_gene_map$in_network & dmp_gene_map$module != "grey"])
total_dmp_genes <- length(dmp_genes_in_network)

fisher_results <- list()
for (mod in real_modules) {
  mod_size <- module_sizes$module_size[module_sizes$module == mod]
  mod_dmp  <- module_dmp_counts$n_genes_with_dmp[module_dmp_counts$module == mod]
  if (is.na(mod_dmp)) mod_dmp <- 0
  
  # 2×2 table: DMP_gene / not_DMP_gene × in_module / not_in_module
  a <- mod_dmp                             # DMP genes in module
  b <- mod_size - mod_dmp                  # Non-DMP genes in module
  c_val <- total_dmp_genes - mod_dmp       # DMP genes outside module
  d <- total_genes_network - mod_size - c_val  # Non-DMP genes outside module
  
  ft <- tryCatch(
    fisher.test(matrix(c(a, c_val, b, d), nrow = 2), alternative = "greater"),
    error = function(e) list(p.value = NA, estimate = NA))
  
  fisher_results[[mod]] <- data.frame(
    module = mod,
    fisher_p = ft$p.value,
    odds_ratio = as.numeric(ft$estimate),
    stringsAsFactors = FALSE
  )
}
fisher_df <- do.call(rbind, fisher_results)
fisher_df$fisher_padj <- p.adjust(fisher_df$fisher_p, method = "BH")

module_dmp_counts <- merge(module_dmp_counts, fisher_df, by = "module")

# 3.4 Add Wilcoxon results from WGCNA
module_dmp_counts <- merge(module_dmp_counts,
                           wilcox_results[, c("Module", "Wilcox_p", "Cohens_d", "Direction")],
                           by.x = "module", by.y = "Module", all.x = TRUE)

# Sort by number of DMPs
module_dmp_counts <- module_dmp_counts[order(-module_dmp_counts$n_dmps), ]

cat("  Module-level DMP summary:\n")
print(as.data.frame(module_dmp_counts[, c("module", "module_size", "n_dmps",
                                          "n_genes_with_dmp", "pct_genes_with_dmp",
                                          "fisher_padj", "Wilcox_p")]))

write.table(module_dmp_counts,
            file.path(dirs$tables, "MXT_module_DMP_DMR_summary.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)


# #############################################################################
# SECTION 4: HUB GENE METHYLATION STATUS
# #############################################################################

cat("\n=== Hub Gene Methylation Analysis ===\n\n")

# 4.1 Which hub genes have DMPs?
hubs_only <- hub_genes_df[hub_genes_df$IsHub == TRUE, ]
candidates_only <- hub_genes_df[hub_genes_df$IsCandidate == TRUE, ]

# DMP info per gene (aggregated)
gene_dmp_summary <- dmp_gene_map %>%
  group_by(gene_id) %>%
  summarise(
    n_dmps         = n(),
    mean_meth_diff = mean(meth_diff, na.rm = TRUE),
    n_hyper        = sum(meth_diff > 0, na.rm = TRUE),
    n_hypo         = sum(meth_diff < 0, na.rm = TRUE),
    regions        = paste(unique(region), collapse = ";"),
    primary_region = names(sort(table(region), decreasing = TRUE))[1],
    .groups = "drop"
  )

# DMR info per gene
if (nrow(dmr_gene_map) > 0) {
  gene_dmr_summary <- dmr_gene_map %>%
    group_by(gene_id) %>%
    summarise(
      n_dmrs        = n(),
      mean_dmr_diff = mean(dmr_diff, na.rm = TRUE),
      dmr_regions   = paste(unique(region), collapse = ";"),
      .groups = "drop"
    )
} else {
  gene_dmr_summary <- data.frame(gene_id = character(), n_dmrs = integer(),
                                 mean_dmr_diff = numeric(), dmr_regions = character())
}

# Annotate hub genes
hubs_meth <- merge(hubs_only, gene_dmp_summary, by.x = "Gene", by.y = "gene_id", all.x = TRUE)
hubs_meth <- merge(hubs_meth, gene_dmr_summary, by.x = "Gene", by.y = "gene_id", all.x = TRUE)
hubs_meth$has_dmp <- !is.na(hubs_meth$n_dmps) & hubs_meth$n_dmps > 0
hubs_meth$has_dmr <- !is.na(hubs_meth$n_dmrs) & hubs_meth$n_dmrs > 0

# Add gene names
if (!is.null(gene_names)) {
  hubs_meth <- merge(hubs_meth, gene_names[, c("gene_id", "display_name")],
                     by.x = "Gene", by.y = "gene_id", all.x = TRUE)
}

cat(sprintf("  Total hub genes: %d\n", nrow(hubs_meth)))
cat(sprintf("  Hub genes with DMPs: %d (%.1f%%)\n",
            sum(hubs_meth$has_dmp), 100 * mean(hubs_meth$has_dmp)))
cat(sprintf("  Hub genes with DMRs: %d (%.1f%%)\n",
            sum(hubs_meth$has_dmr), 100 * mean(hubs_meth$has_dmr)))

# Same for candidates
cand_meth <- merge(candidates_only, gene_dmp_summary, by.x = "Gene", by.y = "gene_id", all.x = TRUE)
cand_meth <- merge(cand_meth, gene_dmr_summary, by.x = "Gene", by.y = "gene_id", all.x = TRUE)
cand_meth$has_dmp <- !is.na(cand_meth$n_dmps) & cand_meth$n_dmps > 0
cand_meth$has_dmr <- !is.na(cand_meth$n_dmrs) & cand_meth$n_dmrs > 0

if (!is.null(gene_names)) {
  cand_meth <- merge(cand_meth, gene_names[, c("gene_id", "display_name")],
                     by.x = "Gene", by.y = "gene_id", all.x = TRUE)
}

cat(sprintf("  Candidate genes with DMPs: %d / %d (%.1f%%)\n",
            sum(cand_meth$has_dmp), nrow(cand_meth), 100 * mean(cand_meth$has_dmp)))

# Hub genes per module affected by DMPs
hub_module_summary <- hubs_meth %>%
  group_by(Module) %>%
  summarise(
    total_hubs     = n(),
    hubs_with_dmp  = sum(has_dmp),
    hubs_with_dmr  = sum(has_dmr),
    pct_dmp        = round(100 * mean(has_dmp), 1),
    .groups = "drop"
  ) %>%
  arrange(desc(hubs_with_dmp))

cat("\n  Hub gene DMP overlap by module:\n")
print(as.data.frame(hub_module_summary))

# Save tables
write.table(hubs_meth, file.path(dirs$tables, "MXT_hub_genes_methylation_status.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(cand_meth, file.path(dirs$tables, "MXT_candidate_genes_methylation_status.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(hub_module_summary, file.path(dirs$tables, "MXT_hub_module_DMP_summary.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Highlight: affected hubs in best Wilcoxon module
best_hubs_affected <- hubs_meth %>%
  filter(Module == best_module_wilcox, has_dmp | has_dmr) %>%
  arrange(desc(Connectivity))

if (nrow(best_hubs_affected) > 0) {
  cat(sprintf("\n  Hub genes in BEST module (%s) with methylation changes: %d\n",
              best_module_wilcox, nrow(best_hubs_affected)))
  cols_to_show <- intersect(c("Gene", "display_name", "Connectivity", "ModuleMembership",
                              "n_dmps", "mean_meth_diff", "primary_region",
                              "n_dmrs", "mean_dmr_diff"), colnames(best_hubs_affected))
  print(as.data.frame(head(best_hubs_affected[, cols_to_show], 15)))
  
  write.table(best_hubs_affected,
              file.path(dirs$tables, paste0("MXT_hub_genes_", best_module_wilcox, "_affected.txt")),
              sep = "\t", quote = FALSE, row.names = FALSE)
}


# #############################################################################
# SECTION 5: REGION BREAKDOWN PER MODULE
# #############################################################################

cat("\n=== DMP Genomic Region Distribution per Module ===\n\n")

module_region <- dmp_gene_map %>%
  filter(in_network, module != "grey") %>%
  group_by(module, region) %>%
  summarise(n_dmps = n(), .groups = "drop")

# Wide format
module_region_wide <- module_region %>%
  pivot_wider(names_from = region, values_from = n_dmps, values_fill = 0)

write.table(module_region_wide, file.path(dirs$tables, "MXT_module_DMP_region_breakdown.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)


# #############################################################################
# SECTION 6: VISUALIZATIONS
# #############################################################################

cat("\n=== Generating Visualizations ===\n\n")

# 6.1 Module DMP count barplot (colored by module)
module_dmp_counts$module_f <- factor(module_dmp_counts$module,
                                     levels = module_dmp_counts$module[order(module_dmp_counts$n_dmps)])

p_mod_dmp <- ggplot(module_dmp_counts, aes(x = module_f, y = n_dmps, fill = module)) +
  geom_col(alpha = 0.85) +
  geom_text(aes(label = n_dmps), hjust = -0.2, size = 3) +
  scale_fill_identity() +
  coord_flip() +
  labs(
    title = "Number of DMPs per WGCNA Module — D. laeve",
    subtitle = "Methylation × Transcriptomics integration",
    x = "Module", y = "Number of DMPs"
  )
save_both(p_mod_dmp, "MXT_11_DMPs_per_module", 12, 8)

# 6.2 Percentage of module genes with DMPs
p_mod_pct <- ggplot(module_dmp_counts, aes(x = reorder(module, pct_genes_with_dmp),
                                           y = pct_genes_with_dmp, fill = module)) +
  geom_col(alpha = 0.85) +
  geom_text(aes(label = sprintf("%.1f%%", pct_genes_with_dmp)), hjust = -0.2, size = 3) +
  scale_fill_identity() +
  coord_flip() +
  labs(
    title = "Percentage of Module Genes Affected by DMPs — D. laeve",
    subtitle = "Proportion of genes in each module that overlap with differentially methylated positions",
    x = "Module", y = "% Genes with DMPs"
  )
save_both(p_mod_pct, "MXT_12_pct_genes_with_DMP_per_module", 12, 8)

# 6.3 Fisher enrichment dot plot
sig_modules <- module_dmp_counts %>% filter(fisher_padj < 0.2)

if (nrow(sig_modules) > 0) {
  p_fisher <- ggplot(module_dmp_counts, aes(x = reorder(module, -fisher_padj),
                                            y = -log10(fisher_padj),
                                            size = n_genes_with_dmp,
                                            color = module)) +
    geom_point(alpha = 0.8) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    scale_color_identity() +
    scale_size_continuous(name = "Genes with DMP", range = c(2, 10)) +
    coord_flip() +
    labs(
      title = "DMP Enrichment per WGCNA Module (Fisher's Exact Test)",
      subtitle = "Red dashed line: padj = 0.05",
      x = "Module", y = "-log10(Adjusted p-value)"
    )
  save_both(p_fisher, "MXT_13_DMP_enrichment_fisher_modules", 12, 8)
}

# 6.4 Stacked bar: DMP regions per module (top 10 modules)
top_modules <- head(module_dmp_counts$module[order(-module_dmp_counts$n_dmps)], 10)
region_top <- module_region %>% filter(module %in% top_modules)

if (nrow(region_top) > 0) {
  p_region_stack <- ggplot(region_top, aes(x = reorder(module, n_dmps),
                                           y = n_dmps, fill = region)) +
    geom_col(position = "stack") +
    scale_fill_brewer(palette = "Set2", name = "Genomic Region") +
    coord_flip() +
    labs(
      title = "DMP Genomic Region Distribution per Module (Top 10)",
      x = "Module", y = "Number of DMPs"
    )
  save_both(p_region_stack, "MXT_14_DMP_region_per_module_stacked", 14, 8)
}

# 6.5 Hub gene methylation status barplot by module
if (nrow(hub_module_summary) > 0) {
  hub_long <- hub_module_summary %>%
    select(Module, hubs_with_dmp, hubs_with_dmr,
           no_meth = total_hubs) %>%
    mutate(no_meth = no_meth - hubs_with_dmp) %>%
    pivot_longer(cols = c(hubs_with_dmp, hubs_with_dmr, no_meth),
                 names_to = "status", values_to = "count")
  
  hub_long$status <- factor(hub_long$status,
                            levels = c("no_meth", "hubs_with_dmp", "hubs_with_dmr"),
                            labels = c("No methylation change", "Has DMPs", "Has DMRs"))
  
  p_hub <- ggplot(hub_long, aes(x = reorder(Module, count), y = count, fill = status)) +
    geom_col(position = "stack") +
    scale_fill_manual(values = c("No methylation change" = "#95A5A6",
                                 "Has DMPs" = "#E74C3C",
                                 "Has DMRs" = "#3498DB"),
                      name = "Methylation Status") +
    coord_flip() +
    labs(
      title = "Hub Gene Methylation Status per WGCNA Module — D. laeve",
      subtitle = "Hub genes: top 10% connectivity + |MM| > 0.7",
      x = "Module", y = "Number of Hub Genes"
    )
  save_both(p_hub, "MXT_15_hub_genes_meth_status_per_module", 12, 8)
}

# 6.6 Bubble plot: Module Wilcoxon significance × DMP burden
p_bubble <- ggplot(module_dmp_counts,
                   aes(x = -log10(Wilcox_p),
                       y = pct_genes_with_dmp,
                       size = n_dmps,
                       color = module)) +
  geom_point(alpha = 0.8) +
  scale_color_identity() +
  scale_size_continuous(name = "# DMPs", range = c(2, 12)) +
  geom_text_repel(aes(label = module), size = 3.5, max.overlaps = 15) +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "red", alpha = 0.7) +
  labs(
    title = "Module Significance (WGCNA Wilcoxon) vs DMP Burden — D. laeve",
    subtitle = "Modules with high expression change AND high methylation change are top targets",
    x = "-log10(Wilcoxon p-value)\n(Expression change in tail regeneration)",
    y = "% Module Genes with DMPs\n(Methylation burden)"
  )
save_both(p_bubble, "MXT_16_module_wilcoxon_vs_DMP_burden", 13, 10)

# 6.7 Heatmap: Module × Region DMP counts
if (nrow(module_region_wide) > 2) {
  hm_data <- as.matrix(module_region_wide[, -1])
  rownames(hm_data) <- module_region_wide$module
  
  # Only show modules with DMPs
  hm_data <- hm_data[rowSums(hm_data) > 0, , drop = FALSE]
  
  if (nrow(hm_data) > 1) {
    hm_expr <- quote(
      pheatmap(hm_data,
               scale = "row",
               cluster_rows = TRUE, cluster_cols = TRUE,
               color = colorRampPalette(c("white", "#FEE08B", "#E74C3C"))(50),
               main = "DMP Distribution: Module × Genomic Region — D. laeve",
               fontsize = 10)
    )
    save_pheatmap_both(hm_expr, "MXT_17_module_region_DMP_heatmap", 10, 10)
  }
}

# 6.8 Hyper vs Hypo DMPs per module
if (all(c("n_hyper", "n_hypo") %in% colnames(module_dmp_counts))) {
  hypo_hyper <- module_dmp_counts %>%
    filter(n_dmps > 0) %>%
    select(module, n_hyper, n_hypo) %>%
    pivot_longer(cols = c(n_hyper, n_hypo), names_to = "direction", values_to = "count")
  
  hypo_hyper$direction <- ifelse(hypo_hyper$direction == "n_hyper", "Hypermethylated", "Hypomethylated")
  
  p_dir <- ggplot(hypo_hyper, aes(x = reorder(module, count), y = count, fill = direction)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Hypermethylated" = "#E74C3C", "Hypomethylated" = "#3498DB")) +
    coord_flip() +
    labs(
      title = "Hyper- vs Hypo-methylated DMPs per Module — D. laeve",
      x = "Module", y = "Number of DMPs", fill = "Direction"
    )
  save_both(p_dir, "MXT_18_hyper_hypo_DMPs_per_module", 12, 8)
}

# 6.9 Affected hub genes in best module: Connectivity vs Methylation
if (nrow(best_hubs_affected) > 0) {
  p_hub_best <- ggplot(best_hubs_affected,
                       aes(x = Connectivity, y = mean_meth_diff,
                           color = primary_region, size = n_dmps)) +
    geom_point(alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    scale_color_brewer(palette = "Set2", name = "Primary Region") +
    scale_size_continuous(name = "# DMPs", range = c(2, 8)) +
    labs(
      title = sprintf("Hub Genes in %s Module: Connectivity vs Methylation", best_module_wilcox),
      subtitle = sprintf("Hub genes affected by DMPs/DMRs in the best Wilcoxon module (p = %s)",
                         signif(me_wilcox$Wilcox_p[1], 3)),
      x = "Intramodular Connectivity",
      y = "Mean Methylation Difference"
    )
  save_both(p_hub_best, paste0("MXT_19_hub_genes_", best_module_wilcox, "_connectivity_vs_meth"), 11, 8)
}


# #############################################################################
# SECTION 7: COMBINED SUMMARY TABLE
# #############################################################################

cat("\n=== Generating Combined Summary ===\n\n")

# Merge module DMP/DMR info with Wilcoxon info for final master table
master_table <- module_dmp_counts %>%
  arrange(desc(n_dmps)) %>%
  select(module, module_size, n_dmps, n_genes_with_dmp, pct_genes_with_dmp,
         n_hyper, n_hypo, mean_meth_diff,
         fisher_padj, Wilcox_p, Cohens_d, Direction,
         any_of(c("n_dmrs", "n_genes_with_dmr")))

# Add hub info
master_table <- merge(master_table, hub_module_summary[, c("Module", "total_hubs",
                                                           "hubs_with_dmp", "hubs_with_dmr")],
                      by.x = "module", by.y = "Module", all.x = TRUE)
master_table[is.na(master_table)] <- 0
master_table <- master_table[order(-master_table$n_dmps), ]

write.table(master_table,
            file.path(dirs$tables, "MXT_MASTER_module_methylation_summary.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("  Master table saved.\n")


# #############################################################################
# SECTION 8: SAVE OBJECTS
# #############################################################################

save(module_dmp_counts, module_region, hub_module_summary,
     hubs_meth, cand_meth, best_hubs_affected,
     dmp_gene_map, dmr_gene_map, gene_dmp_summary,
     fisher_df, master_table,
     file = file.path(dirs$rds, "MXT_Part2_objects.RData"))

cat("  Objects saved.\n")


# #############################################################################
# FINAL SUMMARY
# #############################################################################

cat("\n")
cat("================================================================\n")
cat("   MXT PART 2 COMPLETE — Modules × Methylation\n")
cat("================================================================\n\n")

cat("KEY RESULTS:\n\n")

cat("Top 5 modules by DMP count:\n")
print(head(master_table[, c("module", "module_size", "n_dmps", "n_genes_with_dmp",
                            "pct_genes_with_dmp", "fisher_padj", "Wilcox_p")], 5))

cat(sprintf("\nBest Wilcoxon module (%s):\n", best_module_wilcox))
best_row <- master_table[master_table$module == best_module_wilcox, ]
if (nrow(best_row) > 0) {
  cat(sprintf("  Module size: %d genes\n", best_row$module_size))
  cat(sprintf("  DMPs: %d | Genes with DMPs: %d (%.1f%%)\n",
              best_row$n_dmps, best_row$n_genes_with_dmp, best_row$pct_genes_with_dmp))
  cat(sprintf("  Hub genes with DMPs: %d / %d\n", best_row$hubs_with_dmp, best_row$total_hubs))
  cat(sprintf("  Wilcoxon p: %s | Fisher DMP enrichment padj: %s\n",
              signif(best_row$Wilcox_p, 3), signif(best_row$fisher_padj, 3)))
}

cat(sprintf("\nTotal hub genes affected by DMPs: %d / %d (%.1f%%)\n",
            sum(hubs_meth$has_dmp), nrow(hubs_meth), 100 * mean(hubs_meth$has_dmp)))
cat(sprintf("Total hub genes affected by DMRs: %d / %d (%.1f%%)\n",
            sum(hubs_meth$has_dmr), nrow(hubs_meth), 100 * mean(hubs_meth$has_dmr)))

cat("\nOUTPUT FILES:\n")
cat(sprintf("  Tables:  %s\n", dirs$tables))
cat(sprintf("  PDFs:    %s\n", dirs$pdf))
cat(sprintf("  PNGs:    %s\n", dirs$png))
cat(sprintf("  RDS:     %s\n", dirs$rds))
cat("================================================================\n")