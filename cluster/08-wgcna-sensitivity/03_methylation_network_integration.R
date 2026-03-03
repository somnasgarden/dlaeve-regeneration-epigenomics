#!/usr/bin/env Rscript
# =============================================================================
# Methylation-Network Integration: The Honest Network Parameter
# =============================================================================
# Purpose: Calculate 6 metrics linking WGCNA network topology to methylation
#          status. These metrics answer: does network structure predict
#          epigenetic targeting?
#
# Metrics:
#   1. kWithin vs DMP rate (connectivity predicts methylation?)
#   2. kME vs methylation status (module membership predicts targeting?)
#   3. TOM vs methylation concordance (co-expressed genes co-methylated?)
#   4. Weighted Methylation Coordination Score (WMCS)
#   5. Hub-periphery methylation gradient
#   6. Eigengene-methylation burden correlation
#
# Requires: 64GB RAM, 8 CPUs, ~2-4 hours
# Input:  results/wgcna_sensitivity/data/improved_WGCNA_objects.RData
#         results/03_integration/Tables/MXT_gene_level_meth_vs_expression.txt
#         results/01_methylation/dmrs_annotated.txt
# Output: results/wgcna_sensitivity/
# =============================================================================

library(WGCNA)
library(ggplot2)
library(dplyr)
library(parallel)

enableWGCNAThreads(nThreads = 8)
options(stringsAsFactors = FALSE)

cat("=== Methylation-Network Integration ===\n")
cat("CPUs:", detectCores(), "| Using 8 cores\n\n")

N_CORES <- 8
outdir <- "results/wgcna_sensitivity"

# ── Load data ────────────────────────────────────────────────────────────────

cat("Loading WGCNA objects...\n")
load(file.path(outdir, "data", "improved_WGCNA_objects.RData"))
cat("  Modules:", length(unique(module_colors)) - 1, "+ grey\n")
cat("  Hub genes:", sum(all_hub_genes$IsHub), "\n")

cat("Loading methylation data...\n")
# Waterfall: cluster MXT paths → local paths
gene_meth <- NULL
for (gpath in c("/mnt/data/alfredvar/rlopezt/MXT/Tables/MXT_gene_level_meth_vs_expression.txt",
                "results/03_integration/Tables/MXT_gene_level_meth_vs_expression.txt")) {
  if (file.exists(gpath)) {
    gene_meth <- read.delim(gpath)
    cat("  Gene methylation from:", gpath, "\n")
    break
  }
}
if (is.null(gene_meth)) stop("Cannot find MXT gene-level table.")
cat("  Genes with DMPs:", nrow(gene_meth), "\n")

# Waterfall for DMRs: cluster → local
dmrs <- NULL
for (dpath in c("/mnt/data/alfredvar/rlopezt/Preliminary/dmrs_annotated.txt",
                "results/01_methylation/dmrs_annotated.txt",
                "results/01_methylation/tables/dmrs_annotated.txt")) {
  if (file.exists(dpath)) {
    dmrs <- read.delim(dpath)
    cat("  DMRs from:", dpath, "\n")
    break
  }
}
if (is.null(dmrs)) stop("Cannot find DMR annotated table.")
dmr_genes <- unique(dmrs$nearest_gene)
cat("  DMR genes:", length(dmr_genes), "\n")

# Build lookups
meth_lookup <- gene_meth %>%
  select(gene_id, n_dmps, mean_meth_diff, primary_region, meth_direction) %>%
  distinct(gene_id, .keep_all = TRUE)

has_dmp <- setNames(rep(FALSE, length(module_colors)), names(module_colors))
has_dmp[names(module_colors) %in% meth_lookup$gene_id] <- TRUE

n_dmps <- setNames(rep(0, length(module_colors)), names(module_colors))
n_dmps[meth_lookup$gene_id[meth_lookup$gene_id %in% names(module_colors)]] <-
  meth_lookup$n_dmps[meth_lookup$gene_id %in% names(module_colors)]

meth_direction <- setNames(rep(NA_character_, length(module_colors)), names(module_colors))
matched <- meth_lookup$gene_id[meth_lookup$gene_id %in% names(module_colors)]
meth_direction[matched] <- meth_lookup$meth_direction[match(matched, meth_lookup$gene_id)]

has_dmr <- setNames(names(module_colors) %in% dmr_genes, names(module_colors))

# Region lookup
region_lookup <- setNames(rep(NA_character_, length(module_colors)), names(module_colors))
region_lookup[matched] <- meth_lookup$primary_region[match(matched, meth_lookup$gene_id)]

cat("  Genes in network with DMPs:", sum(has_dmp), "\n")
cat("  Genes in network with DMRs:", sum(has_dmr), "\n\n")


# ── Get module list ──────────────────────────────────────────────────────────

unique_mods <- unique(module_colors)
unique_mods <- unique_mods[unique_mods != "grey"]
cat("Analyzing", length(unique_mods), "modules\n\n")


# ══════════════════════════════════════════════════════════════════════════════
# METRIC 1: kWithin vs DMP Rate
# ══════════════════════════════════════════════════════════════════════════════

cat("--- Metric 1: Intramodular Connectivity vs DMP Rate ---\n")

metric1_results <- data.frame()

for (mod in unique_mods) {
  mod_genes <- names(module_colors)[module_colors == mod]
  if (length(mod_genes) < 30) next

  # Get kWithin from hub genes table
  mod_hubs <- all_hub_genes %>% filter(Module == mod)
  if (nrow(mod_hubs) == 0) next

  # Spearman: kWithin vs n_dmps
  mod_hubs$n_dmps_val <- n_dmps[mod_hubs$Gene]
  mod_hubs$has_dmp_val <- has_dmp[mod_hubs$Gene]

  ct_count <- suppressWarnings(
    cor.test(mod_hubs$kWithin, mod_hubs$n_dmps_val, method = "spearman")
  )
  ct_binary <- suppressWarnings(
    cor.test(mod_hubs$kWithin, as.numeric(mod_hubs$has_dmp_val), method = "spearman")
  )

  metric1_results <- rbind(metric1_results, data.frame(
    module = mod,
    n_genes = nrow(mod_hubs),
    n_with_dmp = sum(mod_hubs$has_dmp_val),
    pct_dmp = round(mean(mod_hubs$has_dmp_val) * 100, 1),
    rho_kWithin_nDMP = round(ct_count$estimate, 4),
    p_kWithin_nDMP = ct_count$p.value,
    rho_kWithin_hasDMP = round(ct_binary$estimate, 4),
    p_kWithin_hasDMP = ct_binary$p.value
  ))

  cat(sprintf("  %-12s n=%4d  rho(kW,nDMP)=%+.4f p=%.4f  rho(kW,hasDMP)=%+.4f p=%.4f\n",
              mod, nrow(mod_hubs),
              ct_count$estimate, ct_count$p.value,
              ct_binary$estimate, ct_binary$p.value))
}

metric1_results$padj_count <- p.adjust(metric1_results$p_kWithin_nDMP, method = "BH")
metric1_results$padj_binary <- p.adjust(metric1_results$p_kWithin_hasDMP, method = "BH")

write.table(metric1_results, file.path(outdir, "tables", "S08_kWithin_vs_DMP.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)


# ══════════════════════════════════════════════════════════════════════════════
# METRIC 2: kME vs Methylation Status
# ══════════════════════════════════════════════════════════════════════════════

cat("\n--- Metric 2: Module Membership vs Methylation ---\n")

metric2_results <- data.frame()

for (mod in unique_mods) {
  mod_hubs <- all_hub_genes %>% filter(Module == mod)
  if (nrow(mod_hubs) < 30) next

  mod_hubs$has_dmp_val <- has_dmp[mod_hubs$Gene]
  mod_hubs$n_dmps_val <- n_dmps[mod_hubs$Gene]

  # Wilcoxon: kME for DMP vs non-DMP genes
  dmp_kme <- mod_hubs$kME[mod_hubs$has_dmp_val]
  nodmp_kme <- mod_hubs$kME[!mod_hubs$has_dmp_val]

  if (length(dmp_kme) >= 5 & length(nodmp_kme) >= 5) {
    wt <- wilcox.test(dmp_kme, nodmp_kme)
    wt_p <- wt$p.value
    diff_median <- median(dmp_kme) - median(nodmp_kme)
  } else {
    wt_p <- NA
    diff_median <- NA
  }

  # Spearman: kME vs n_dmps
  ct <- suppressWarnings(
    cor.test(mod_hubs$kME, mod_hubs$n_dmps_val, method = "spearman")
  )

  metric2_results <- rbind(metric2_results, data.frame(
    module = mod,
    n_genes = nrow(mod_hubs),
    median_kME_DMP = round(median(dmp_kme, na.rm = TRUE), 4),
    median_kME_noDMP = round(median(nodmp_kme, na.rm = TRUE), 4),
    diff_median_kME = round(diff_median, 4),
    wilcox_p = wt_p,
    rho_kME_nDMP = round(ct$estimate, 4),
    p_kME_nDMP = ct$p.value
  ))

  cat(sprintf("  %-12s kME: DMP=%.3f noDMP=%.3f diff=%+.3f (Wilcox p=%.4f)  rho(kME,n)=%+.4f\n",
              mod, median(dmp_kme, na.rm = TRUE), median(nodmp_kme, na.rm = TRUE),
              ifelse(is.na(diff_median), 0, diff_median),
              ifelse(is.na(wt_p), 1, wt_p),
              ct$estimate))
}

metric2_results$wilcox_padj <- p.adjust(metric2_results$wilcox_p, method = "BH")
metric2_results$rho_padj <- p.adjust(metric2_results$p_kME_nDMP, method = "BH")

write.table(metric2_results, file.path(outdir, "tables", "S09_kME_vs_methylation.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)


# ══════════════════════════════════════════════════════════════════════════════
# METRIC 3 & 4: TOM Concordance and WMCS
# ══════════════════════════════════════════════════════════════════════════════

cat("\n--- Metrics 3-4: TOM Concordance & WMCS ---\n")
cat("  Loading TOM blocks (memory-intensive)...\n")

# Load TOM blocks
tom_files <- list.files(file.path(outdir, "data"), pattern = "improved_TOM.*\\.RData$",
                        full.names = TRUE)

if (length(tom_files) == 0) {
  cat("  WARNING: No TOM files found. Skipping TOM-based metrics.\n")
  cat("  Make sure 02_improved_wgcna.R was run first with saveTOMs=TRUE\n")
  tom_available <- FALSE
} else {
  tom_available <- TRUE
  cat("  Found", length(tom_files), "TOM block files\n")
}

metric34_results <- data.frame()

if (tom_available) {
  # Process one block at a time to manage memory
  for (tom_file in tom_files) {
    cat("  Loading:", basename(tom_file), "... ")
    load(tom_file)
    tom_obj_name <- ls()[grepl("TOM", ls())]
    if (length(tom_obj_name) == 0) {
      # Try getting the first large matrix
      for (obj in ls()) {
        if (is.matrix(get(obj)) && nrow(get(obj)) > 100) {
          tom_obj_name <- obj
          break
        }
      }
    }

    if (length(tom_obj_name) > 0) {
      TOM_block <- get(tom_obj_name[1])
      cat("dim =", nrow(TOM_block), "x", ncol(TOM_block), "\n")

      # Get gene names for this block (from net$blockGenes)
      block_num <- as.numeric(gsub(".*block-(\\d+)\\.RData", "\\1", basename(tom_file)))
      if (!is.na(block_num) && block_num <= length(net$blockGenes)) {
        block_gene_indices <- net$blockGenes[[block_num]]
        block_genes <- colnames(expr_matrix)[block_gene_indices]

        if (length(block_genes) == nrow(TOM_block)) {
          rownames(TOM_block) <- block_genes
          colnames(TOM_block) <- block_genes

          # For each module represented in this block, compute WMCS
          for (mod in unique_mods) {
            mod_genes_in_block <- intersect(block_genes,
                                            names(module_colors)[module_colors == mod])
            dmp_genes_in_block <- mod_genes_in_block[has_dmp[mod_genes_in_block]]

            if (length(dmp_genes_in_block) >= 10) {
              # Direction vector: +1 hyper, -1 hypo
              dir_vec <- ifelse(meth_direction[dmp_genes_in_block] == "Hyper", 1, -1)
              dir_vec[is.na(dir_vec)] <- 0

              # Pairwise concordance weighted by TOM
              tom_sub <- TOM_block[dmp_genes_in_block, dmp_genes_in_block]
              n_genes_sub <- length(dmp_genes_in_block)

              wmcs_num <- 0
              wmcs_den <- 0
              n_concordant <- 0
              n_total_pairs <- 0

              for (ii in 1:(n_genes_sub - 1)) {
                for (jj in (ii + 1):n_genes_sub) {
                  conc <- dir_vec[ii] * dir_vec[jj]
                  tom_val <- tom_sub[ii, jj]
                  wmcs_num <- wmcs_num + tom_val * conc
                  wmcs_den <- wmcs_den + tom_val
                  if (conc == 1) n_concordant <- n_concordant + 1
                  n_total_pairs <- n_total_pairs + 1
                }
              }

              wmcs <- ifelse(wmcs_den > 0, wmcs_num / wmcs_den, NA)
              pct_concordant <- ifelse(n_total_pairs > 0,
                                        n_concordant / n_total_pairs * 100, NA)

              # Check if this module already has a result (from a previous block)
              existing <- which(metric34_results$module == mod)
              if (length(existing) > 0) {
                # Accumulate
                metric34_results$wmcs_num[existing] <- metric34_results$wmcs_num[existing] + wmcs_num
                metric34_results$wmcs_den[existing] <- metric34_results$wmcs_den[existing] + wmcs_den
                metric34_results$n_concordant[existing] <- metric34_results$n_concordant[existing] + n_concordant
                metric34_results$n_pairs[existing] <- metric34_results$n_pairs[existing] + n_total_pairs
                metric34_results$n_dmp_genes[existing] <- max(metric34_results$n_dmp_genes[existing],
                                                              length(dmp_genes_in_block))
              } else {
                metric34_results <- rbind(metric34_results, data.frame(
                  module = mod,
                  n_dmp_genes = length(dmp_genes_in_block),
                  n_pairs = n_total_pairs,
                  n_concordant = n_concordant,
                  wmcs_num = wmcs_num,
                  wmcs_den = wmcs_den
                ))
              }
            }
          }
        }
      }

      rm(TOM_block)
      gc()
    }
  }

  # Finalize WMCS
  if (nrow(metric34_results) > 0) {
    metric34_results$WMCS <- round(metric34_results$wmcs_num / metric34_results$wmcs_den, 4)
    metric34_results$pct_concordant <- round(metric34_results$n_concordant /
                                               metric34_results$n_pairs * 100, 1)

    cat("\n  WMCS results:\n")
    for (i in 1:nrow(metric34_results)) {
      r <- metric34_results[i, ]
      cat(sprintf("    %-12s WMCS=%+.4f  concordance=%.1f%%  pairs=%d\n",
                  r$module, r$WMCS, r$pct_concordant, r$n_pairs))
    }

    write.table(metric34_results %>% select(module, n_dmp_genes, n_pairs, n_concordant,
                                              pct_concordant, WMCS),
                file.path(outdir, "tables", "S10_WMCS_results.tsv"),
                sep = "\t", row.names = FALSE, quote = FALSE)
  }
}


# ══════════════════════════════════════════════════════════════════════════════
# METRIC 5: Hub-Periphery Methylation Gradient
# ══════════════════════════════════════════════════════════════════════════════

cat("\n--- Metric 5: Hub-Periphery Methylation Gradient ---\n")

metric5_results <- data.frame()

for (mod in unique_mods) {
  mod_hubs <- all_hub_genes %>% filter(Module == mod)
  if (nrow(mod_hubs) < 50) next

  # Sort by kWithin and divide into quintiles
  mod_hubs <- mod_hubs %>% arrange(desc(kWithin))
  mod_hubs$quintile <- cut(seq_len(nrow(mod_hubs)),
                           breaks = quantile(seq_len(nrow(mod_hubs)),
                                             probs = seq(0, 1, 0.2)),
                           labels = c("Q1_hub", "Q2", "Q3", "Q4", "Q5_periph"),
                           include.lowest = TRUE)

  mod_hubs$has_dmp_val <- has_dmp[mod_hubs$Gene]
  mod_hubs$n_dmps_val <- n_dmps[mod_hubs$Gene]

  quintile_summary <- mod_hubs %>%
    group_by(quintile) %>%
    summarise(
      n_genes = n(),
      dmp_rate = mean(has_dmp_val) * 100,
      mean_n_dmps = mean(n_dmps_val),
      .groups = "drop"
    )

  # Trend test: Spearman of quintile rank vs DMP rate
  quintile_summary$rank <- 1:5
  ct <- suppressWarnings(
    cor.test(quintile_summary$rank, quintile_summary$dmp_rate, method = "spearman")
  )

  gradient <- quintile_summary$dmp_rate[1] - quintile_summary$dmp_rate[5]

  metric5_results <- rbind(metric5_results, data.frame(
    module = mod,
    n_genes = nrow(mod_hubs),
    Q1_hub_dmp_pct = round(quintile_summary$dmp_rate[1], 1),
    Q5_periph_dmp_pct = round(quintile_summary$dmp_rate[5], 1),
    gradient = round(gradient, 1),
    trend_rho = round(ct$estimate, 4),
    trend_p = ct$p.value
  ))

  cat(sprintf("  %-12s Q1(hub)=%.1f%% Q5(periph)=%.1f%% gradient=%+.1f rho=%.3f p=%.3f\n",
              mod, quintile_summary$dmp_rate[1], quintile_summary$dmp_rate[5],
              gradient, ct$estimate, ct$p.value))
}

metric5_results$trend_padj <- p.adjust(metric5_results$trend_p, method = "BH")

write.table(metric5_results, file.path(outdir, "tables", "S11_hub_periphery_gradient.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)


# ══════════════════════════════════════════════════════════════════════════════
# METRIC 6: Eigengene-Methylation Burden Correlation
# ══════════════════════════════════════════════════════════════════════════════

cat("\n--- Metric 6: Eigengene-Methylation Burden ---\n")

module_burden <- data.frame()
for (mod in unique_mods) {
  mod_genes <- names(module_colors)[module_colors == mod]
  n_total <- length(mod_genes)
  n_dmp <- sum(has_dmp[mod_genes])
  n_dmr <- sum(has_dmr[mod_genes])
  mean_n_dmps <- mean(n_dmps[mod_genes])

  module_burden <- rbind(module_burden, data.frame(
    module = mod,
    n_genes = n_total,
    n_dmp_genes = n_dmp,
    pct_dmp = round(n_dmp / n_total * 100, 1),
    n_dmr_genes = n_dmr,
    pct_dmr = round(n_dmr / n_total * 100, 1),
    mean_n_dmps = round(mean_n_dmps, 2)
  ))
}

# If we have Wilcoxon ME differences, correlate with methylation burden
me_wilcox_file <- file.path(outdir, "tables", "S07_ME_wilcoxon_improved.tsv")
if (file.exists(me_wilcox_file)) {
  me_wilcox <- read.delim(me_wilcox_file)

  module_burden <- module_burden %>%
    left_join(me_wilcox %>% select(module, cohens_d, direction, p_value),
              by = "module")

  # Correlate |Cohen's d| with methylation burden
  if (sum(!is.na(module_burden$cohens_d)) >= 5) {
    ct <- cor.test(abs(module_burden$cohens_d), module_burden$pct_dmp, method = "spearman")
    cat(sprintf("  |Cohen's d| vs %%DMP: rho=%.4f, p=%.4f (n=%d modules)\n",
                ct$estimate, ct$p.value, sum(!is.na(module_burden$cohens_d))))

    ct2 <- cor.test(abs(module_burden$cohens_d), module_burden$pct_dmr, method = "spearman")
    cat(sprintf("  |Cohen's d| vs %%DMR: rho=%.4f, p=%.4f\n",
                ct2$estimate, ct2$p.value))
  }
}

write.table(module_burden, file.path(outdir, "tables", "S12_module_methylation_burden.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)


# ══════════════════════════════════════════════════════════════════════════════
# REGION-STRATIFIED ANALYSIS
# ══════════════════════════════════════════════════════════════════════════════

cat("\n--- Region-Stratified kWithin vs DMP Analysis ---\n")

region_results <- data.frame()
for (reg in c("Intergenic", "Intron", "Promoter", "Exon")) {
  for (mod in unique_mods) {
    mod_hubs <- all_hub_genes %>% filter(Module == mod)
    if (nrow(mod_hubs) < 30) next

    # Genes with DMPs in this region
    reg_genes <- meth_lookup$gene_id[meth_lookup$primary_region == reg]
    mod_hubs$has_reg_dmp <- mod_hubs$Gene %in% reg_genes

    if (sum(mod_hubs$has_reg_dmp) >= 5) {
      ct <- suppressWarnings(
        cor.test(mod_hubs$kWithin, as.numeric(mod_hubs$has_reg_dmp), method = "spearman")
      )

      region_results <- rbind(region_results, data.frame(
        module = mod,
        region = reg,
        n_genes = nrow(mod_hubs),
        n_with_dmp_region = sum(mod_hubs$has_reg_dmp),
        rho = round(ct$estimate, 4),
        p_value = ct$p.value
      ))
    }
  }
}

if (nrow(region_results) > 0) {
  region_results$padj <- p.adjust(region_results$p_value, method = "BH")

  cat("\n  Region-stratified results (significant only):\n")
  sig_reg <- region_results %>% filter(padj < 0.1)
  if (nrow(sig_reg) > 0) {
    for (i in 1:nrow(sig_reg)) {
      r <- sig_reg[i, ]
      cat(sprintf("    %-12s × %-10s rho=%+.4f padj=%.4f\n",
                  r$module, r$region, r$rho, r$padj))
    }
  } else {
    cat("    None significant at padj < 0.10\n")
  }

  write.table(region_results, file.path(outdir, "tables", "S13_region_stratified_kWithin.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
}


# ══════════════════════════════════════════════════════════════════════════════
# GENE-LEVEL MASTER TABLE
# ══════════════════════════════════════════════════════════════════════════════

cat("\n--- Gene-Level Master Table ---\n")

master <- data.frame(
  gene_id = names(module_colors),
  module = module_colors,
  has_dmp = has_dmp[names(module_colors)],
  n_dmps = n_dmps[names(module_colors)],
  has_dmr = has_dmr[names(module_colors)],
  meth_direction = meth_direction[names(module_colors)],
  primary_region = region_lookup[names(module_colors)],
  stringsAsFactors = FALSE
)

# Add kWithin and kME from hub genes table
master <- master %>%
  left_join(all_hub_genes %>% select(Gene, kWithin, kME, IsHub),
            by = c("gene_id" = "Gene"))

# Add expression (waterfall: cluster → local)
deg_file <- NULL
for (dp in c("Part2_DEGs/data/DEresults_amputated_vs_control.csv",
             "results/02_rnaseq/Part2_DEGs/data/DEresults_amputated_vs_control.csv")) {
  if (file.exists(dp)) { deg_file <- dp; break }
}
if (!is.null(deg_file) && file.exists(deg_file)) {
  deg <- read.csv(deg_file)
  colnames(deg)[1] <- "gene_id"
  master <- master %>%
    left_join(deg %>% select(gene_id, log2FoldChange, padj) %>%
                rename(deg_log2FC = log2FoldChange, deg_padj = padj),
              by = "gene_id")
}

write.table(master, file.path(outdir, "data", "gene_level_master_table.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("  Saved: gene_level_master_table.tsv (", nrow(master), "genes)\n")


# ══════════════════════════════════════════════════════════════════════════════
# VISUALIZATION
# ══════════════════════════════════════════════════════════════════════════════

cat("\n--- Visualization ---\n")

# Plot 1: kWithin vs n_dmps scatter (top 6 modules)
top_mods <- metric1_results %>% arrange(p_kWithin_nDMP) %>% head(6)
plot_data <- all_hub_genes %>%
  filter(Module %in% top_mods$module) %>%
  mutate(n_dmps_val = n_dmps[Gene],
         has_dmp_val = has_dmp[Gene])

p1 <- ggplot(plot_data, aes(x = kWithin, y = n_dmps_val, color = has_dmp_val)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 0.5) +
  facet_wrap(~ Module, scales = "free", ncol = 3) +
  scale_color_manual(values = c("TRUE" = "#E74C3C", "FALSE" = "#BDC3C7")) +
  labs(title = "Intramodular Connectivity vs DMP Count",
       x = "kWithin", y = "Number of DMPs", color = "Has DMP") +
  theme_minimal(base_size = 11)

ggsave(file.path(outdir, "plots", "S05_kWithin_vs_nDMP.png"), p1, width = 12, height = 8, dpi = 300)
ggsave(file.path(outdir, "plots", "S05_kWithin_vs_nDMP.pdf"), p1, width = 12, height = 8)
cat("  Saved: S05_kWithin_vs_nDMP\n")

# Plot 2: Hub-periphery gradient
if (nrow(metric5_results) > 0) {
  gradient_data <- data.frame()
  for (mod in metric5_results$module) {
    mod_hubs <- all_hub_genes %>% filter(Module == mod) %>% arrange(desc(kWithin))
    mod_hubs$quintile <- cut(seq_len(nrow(mod_hubs)),
                             breaks = quantile(seq_len(nrow(mod_hubs)), probs = seq(0, 1, 0.2)),
                             labels = c("Q1\n(Hub)", "Q2", "Q3", "Q4", "Q5\n(Periph)"),
                             include.lowest = TRUE)
    mod_hubs$has_dmp_val <- has_dmp[mod_hubs$Gene]

    qsum <- mod_hubs %>%
      group_by(quintile) %>%
      summarise(dmp_rate = mean(has_dmp_val) * 100, .groups = "drop") %>%
      mutate(module = mod)
    gradient_data <- rbind(gradient_data, qsum)
  }

  p2 <- ggplot(gradient_data, aes(x = quintile, y = dmp_rate, group = module, color = module)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    labs(title = "Hub-Periphery Methylation Gradient",
         subtitle = "DMP rate by connectivity quintile within each module",
         x = "Connectivity Quintile", y = "% Genes with DMPs", color = "Module") +
    theme_minimal(base_size = 12)

  ggsave(file.path(outdir, "plots", "S06_hub_periphery_gradient.png"), p2, width = 10, height = 7, dpi = 300)
  ggsave(file.path(outdir, "plots", "S06_hub_periphery_gradient.pdf"), p2, width = 10, height = 7)
  cat("  Saved: S06_hub_periphery_gradient\n")
}

# Plot 3: Module methylation burden bubble plot
if ("cohens_d" %in% colnames(module_burden)) {
  p3 <- ggplot(module_burden %>% filter(!is.na(cohens_d)),
               aes(x = abs(cohens_d), y = pct_dmp,
                   size = n_genes, color = module)) +
    geom_point(alpha = 0.7) +
    geom_text(aes(label = module), size = 3, vjust = -1, show.legend = FALSE) +
    scale_size_continuous(range = c(3, 15), name = "Module Size") +
    labs(title = "Module Expression Change vs Methylation Burden",
         subtitle = "|Cohen's d| (tail ctrl vs amp) vs % genes with DMPs",
         x = "|Cohen's d| (Expression Difference)",
         y = "% Genes with DMPs") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "right")

  ggsave(file.path(outdir, "plots", "S07_burden_bubble.png"), p3, width = 10, height = 7, dpi = 300)
  ggsave(file.path(outdir, "plots", "S07_burden_bubble.pdf"), p3, width = 10, height = 7)
  cat("  Saved: S07_burden_bubble\n")
}


# ══════════════════════════════════════════════════════════════════════════════
# SUMMARY
# ══════════════════════════════════════════════════════════════════════════════

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("METHYLATION-NETWORK INTEGRATION SUMMARY\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat("METRIC 1 — kWithin vs DMP rate:\n")
sig_m1 <- metric1_results %>% filter(padj_count < 0.05 | padj_binary < 0.05)
if (nrow(sig_m1) > 0) {
  for (i in 1:nrow(sig_m1)) {
    cat(sprintf("  %-12s rho=%+.4f padj=%.4f\n",
                sig_m1$module[i], sig_m1$rho_kWithin_nDMP[i], sig_m1$padj_count[i]))
  }
} else {
  cat("  No modules show significant kWithin-DMP correlation\n")
}

cat("\nMETRIC 2 — kME vs methylation:\n")
sig_m2 <- metric2_results %>% filter(wilcox_padj < 0.05 | rho_padj < 0.05)
if (nrow(sig_m2) > 0) {
  for (i in 1:nrow(sig_m2)) {
    cat(sprintf("  %-12s kME_diff=%+.4f padj=%.4f\n",
                sig_m2$module[i], sig_m2$diff_median_kME[i], sig_m2$wilcox_padj[i]))
  }
} else {
  cat("  No modules show significant kME-methylation relationship\n")
}

if (nrow(metric34_results) > 0) {
  cat("\nMETRIC 3-4 — TOM concordance & WMCS:\n")
  for (i in 1:nrow(metric34_results)) {
    r <- metric34_results[i, ]
    cat(sprintf("  %-12s WMCS=%+.4f concordance=%.1f%%\n",
                r$module, r$WMCS, r$pct_concordant))
  }
}

cat("\nMETRIC 5 — Hub-periphery gradient:\n")
sig_m5 <- metric5_results %>% filter(trend_padj < 0.05)
if (nrow(sig_m5) > 0) {
  for (i in 1:nrow(sig_m5)) {
    cat(sprintf("  %-12s Q1=%.1f%% Q5=%.1f%% gradient=%+.1f\n",
                sig_m5$module[i], sig_m5$Q1_hub_dmp_pct[i],
                sig_m5$Q5_periph_dmp_pct[i], sig_m5$gradient[i]))
  }
} else {
  cat("  No significant hub-periphery gradients detected\n")
}

cat("\nINTERPRETATION GUIDE:\n")
cat("  Metric 1: rho > 0 → highly connected genes get MORE DMPs (hub targeting)\n")
cat("  Metric 2: diff > 0 → DMP genes have HIGHER kME (core targeting)\n")
cat("  Metric 4: WMCS > 0 → co-expressed gene pairs tend to have same meth direction\n")
cat("  Metric 5: gradient > 0 → hubs are MORE methylated than periphery\n")
cat("  If all NS → methylation targeting is INDEPENDENT of network topology\n")
cat("             This is itself an important finding (stochastic regulation)\n")

cat("\n=== Script 03 complete ===\n")
