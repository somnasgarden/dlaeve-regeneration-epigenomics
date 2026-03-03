#!/usr/bin/env Rscript
# =============================================================================
# Blue Module Enrichment Plots (Relaxed Threshold)
# =============================================================================
# The WGCNA pipeline skips plots when padj > 0.05 after correction.
# Blue module has enrichment data but no plots for BP, MF, CC, Reactome.
# This script generates plots using top 15 terms by raw p-value.
# =============================================================================
# Run locally: Rscript local/01_blue_module_enrichment_plots.R
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

theme_set(theme_minimal() +
            theme(text = element_text(size = 12),
                  plot.title = element_text(size = 14, face = "bold"),
                  axis.title = element_text(size = 12),
                  legend.text = element_text(size = 10)))

# ── Paths ─────────────────────────────────────────────────────────────────────

DATA_DIR   <- "results/02_rnaseq/Part5_WGCNA/data"
OUTPUT_DIR <- "results/02_rnaseq/Part5_WGCNA/plots"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Categories to plot (KEGG already has a plot)
categories <- c(
  BP       = "BP_enrichment_blue.tsv",
  MF       = "MF_enrichment_blue.tsv",
  CC       = "CC_enrichment_blue.tsv",
  Reactome = "Reactome_enrichment_blue.tsv"
)

# ── Plot function ─────────────────────────────────────────────────────────────

plot_enrichment <- function(df, category, n_terms = 15) {
  if (nrow(df) == 0) {
    cat("  No terms for", category, "— skipping\n")
    return(invisible(NULL))
  }

  top <- df %>%
    arrange(pvalue) %>%
    head(n_terms) %>%
    mutate(
      neg_log10_p = -log10(pvalue),
      Description = factor(Description, levels = rev(Description)),
      Significant = ifelse(padj < 0.05, "padj < 0.05", "padj >= 0.05")
    )

  p <- ggplot(top, aes(x = neg_log10_p, y = Description, fill = Significant)) +
    geom_col(width = 0.7) +
    scale_fill_manual(
      values = c("padj < 0.05" = "#2166AC", "padj >= 0.05" = "#92C5DE"),
      name = "Significance"
    ) +
    labs(
      title = paste0(category, " Enrichment — Blue Module (D. laeve)"),
      subtitle = paste0("Top ", nrow(top), " terms by raw p-value"),
      x = expression(-log[10](p)),
      y = NULL
    ) +
    theme(
      axis.text.y = element_text(size = 9),
      plot.subtitle = element_text(size = 10, color = "grey40")
    )

  # Save PDF + PNG
  base_name <- paste0("blue_", tolower(category), "_enrichment")
  ggsave(file.path(OUTPUT_DIR, paste0(base_name, ".pdf")),
         p, width = 10, height = 6)
  ggsave(file.path(OUTPUT_DIR, paste0(base_name, ".png")),
         p, width = 10, height = 6, dpi = 300)

  cat("  Saved:", base_name, "\n")
}

# ── Generate plots ────────────────────────────────────────────────────────────

cat("Generating blue module enrichment plots (relaxed threshold)\n")
cat("==========================================================\n\n")

for (cat_name in names(categories)) {
  file_path <- file.path(DATA_DIR, categories[cat_name])

  if (!file.exists(file_path)) {
    cat("  Missing:", file_path, "— skipping\n")
    next
  }

  df <- read.delim(file_path, stringsAsFactors = FALSE)
  cat(cat_name, ":", nrow(df), "terms loaded\n")
  plot_enrichment(df, cat_name)
}

cat("\nDone. Plots saved to:", OUTPUT_DIR, "\n")
