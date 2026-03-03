#!/usr/bin/env Rscript
# =============================================================================
# Part 3: Genomic Annotation, Functional Enrichment & Interpretation
# =============================================================================
# Prerequisites: Run Parts 1 & 2 first
# Note: This script uses your custom GFF for annotation instead of
#       TxDb.Hsapiens. Adjust org.Db if a species-specific one exists.
# =============================================================================

suppressPackageStartupMessages({
  library(bsseq)
  library(GenomicRanges)
  library(rtracklayer)
  library(data.table)
  library(ggplot2)
  library(pheatmap)
  library(dplyr)
  library(tidyr)
  library(scales)
})

# Try loading annotation packages — these are for functional enrichment.
# If your organism has no org.Db, we'll do custom annotation from the GFF only.
has_chipseeker    <- requireNamespace("ChIPseeker", quietly = TRUE)
has_clusterprof   <- requireNamespace("clusterProfiler", quietly = TRUE)

# For a non-model organism you may not have TxDb or org.Db.
# We'll build a TxDb from your GFF directly.
if (requireNamespace("GenomicFeatures", quietly = TRUE)) {
  library(GenomicFeatures)
} else {
  BiocManager::install("GenomicFeatures", update = FALSE, ask = FALSE)
  library(GenomicFeatures)
}

if (has_chipseeker)  library(ChIPseeker)
if (has_clusterprof) library(clusterProfiler)

options(stringsAsFactors = FALSE)
options(scipen = 999)

theme_set(theme_minimal() +
            theme(text = element_text(size = 12),
                  plot.title = element_text(size = 14, face = "bold")))

# >>> EDIT THIS <<<
BASE_DIR <- "/mnt/data/alfredvar/"
RESULTS_DIR <- "/mnt/data/alfredvar/rlopezt/feb_w2_scripts/Preliminary/"

GFF_FILE <- file.path(BASE_DIR,
                      "30-Genoma/31-Alternative_Annotation_EviAnn/derLaeGenome_namesDlasi_v2.fasta.functional_note.pseudo_label.gff")

keep_chr <- paste0("chr", 1:31)

# ── 1. Load Previous Results ────────────────────────────────────────────────

cat("Loading results from Parts 1 & 2...\n")
bs_obj      <- readRDS(file.path(RESULTS_DIR, "bsseq_object.rds"))
dmp_results <- readRDS(file.path(RESULTS_DIR, "dmp_results.rds"))
dmr_results <- readRDS(file.path(RESULTS_DIR, "dmr_results.rds"))

cat(sprintf("  BSseq: %s sites | DMPs: %s | DMRs: %s\n",
            format(nrow(bs_obj), big.mark = ","),
            format(nrow(dmp_results), big.mark = ","),
            format(nrow(dmr_results), big.mark = ",")))

# ── 2. Build TxDb from Your GFF ─────────────────────────────────────────────

cat("\nBuilding TxDb from GFF annotation...\n")
# This creates a transcript database from your custom GFF
# Adjust format if needed ("gff3" or "gtf")
txdb <- makeTxDbFromGFF(GFF_FILE, format = "gff3")
cat("  TxDb created successfully.\n")

# Also load raw GFF for gene-level info
genes_gff <- import(GFF_FILE)
genes_gff <- genes_gff[seqnames(genes_gff) %in% keep_chr]
genes     <- genes_gff[genes_gff$type == "gene"]

# Build gene ID mapping table
gene_info <- data.frame(
  gene_id   = mcols(genes)$ID,
  gene_name = ifelse(!is.null(mcols(genes)$Name), mcols(genes)$Name, mcols(genes)$ID),
  chr       = as.character(seqnames(genes)),
  start     = start(genes),
  end       = end(genes),
  strand    = as.character(strand(genes)),
  stringsAsFactors = FALSE
)

# Classify gene types
gene_info$gene_type <- ifelse(grepl("XLOC.*lncRNA", gene_info$gene_id), "lncRNA",
                              ifelse(grepl("^LOC_", gene_info$gene_id), "protein_coding", "other"))

cat(sprintf("  Gene annotation loaded: %s genes\n", format(nrow(gene_info), big.mark = ",")))
cat("  Gene types:\n")
print(table(gene_info$gene_type))

# ── 3. Custom Annotation Function ───────────────────────────────────────────
# For non-model organisms without ChIPseeker-compatible annoDb,
# we annotate DMPs/DMRs by overlap with genomic features from the GFF.

annotate_with_gff <- function(query_gr, genes_gr, promoter_upstream = 3000,
                              promoter_downstream = 500) {
  # Define genomic feature regions
  promoters_gr <- promoters(genes_gr, upstream = promoter_upstream,
                            downstream = promoter_downstream)
  gene_bodies  <- genes_gr
  
  # Find overlaps
  query_df <- as.data.frame(query_gr)
  query_df$annotation <- "Intergenic"
  query_df$nearest_gene <- NA
  query_df$distance_to_gene <- NA
  
  # Promoter overlaps
  prom_hits <- findOverlaps(query_gr, promoters_gr)
  if (length(prom_hits) > 0) {
    query_df$annotation[queryHits(prom_hits)] <- "Promoter"
    query_df$nearest_gene[queryHits(prom_hits)] <-
      mcols(genes_gr)$ID[subjectHits(prom_hits)]
  }
  
  # Gene body overlaps (only for sites not already in promoters)
  intergenic_idx <- which(query_df$annotation == "Intergenic")
  if (length(intergenic_idx) > 0) {
    body_hits <- findOverlaps(query_gr[intergenic_idx], gene_bodies)
    if (length(body_hits) > 0) {
      orig_idx <- intergenic_idx[queryHits(body_hits)]
      query_df$annotation[orig_idx] <- "Gene Body"
      query_df$nearest_gene[orig_idx] <-
        mcols(genes_gr)$ID[subjectHits(body_hits)]
    }
  }
  
  # For remaining intergenic: find nearest gene
  still_intergenic <- which(query_df$annotation == "Intergenic")
  if (length(still_intergenic) > 0) {
    nearest_idx <- nearest(query_gr[still_intergenic], genes_gr)
    valid       <- !is.na(nearest_idx)
    query_df$nearest_gene[still_intergenic[valid]] <-
      mcols(genes_gr)$ID[nearest_idx[valid]]
    query_df$distance_to_gene[still_intergenic[valid]] <-
      distance(query_gr[still_intergenic[valid]], genes_gr[nearest_idx[valid]])
  }
  
  return(query_df)
}

# ── 4. Annotate DMPs ────────────────────────────────────────────────────────

cat("\n--- Annotating DMPs ---\n")

if (nrow(dmp_results) > 0) {
  dmp_gr <- GRanges(seqnames = dmp_results$chr,
                    ranges = IRanges(start = dmp_results$pos, width = 1),
                    strand = "*",
                    methylation_diff = dmp_results$diff,
                    pvalue = dmp_results$pval,
                    fdr = dmp_results$fdr)
  
  dmp_anno <- annotate_with_gff(dmp_gr, genes)
  
  write.table(dmp_anno, file.path(RESULTS_DIR, "dmps_annotated.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Annotation summary
  anno_summary_dmp <- as.data.frame(table(dmp_anno$annotation))
  colnames(anno_summary_dmp) <- c("Feature", "Count")
  anno_summary_dmp$Percentage <- round(anno_summary_dmp$Count / sum(anno_summary_dmp$Count) * 100, 2)
  
  cat("  DMP annotation distribution:\n")
  print(anno_summary_dmp)
  
  # Annotation pie chart
  p_anno_pie <- ggplot(anno_summary_dmp, aes(x = "", y = Count, fill = Feature)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y") +
    scale_fill_brewer(palette = "Set2") +
    labs(title = "DMP Genomic Feature Distribution") +
    theme_void() +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  
  ggsave(file.path(RESULTS_DIR, "dmp_annotation_pie.pdf"), p_anno_pie,
         width = 8, height = 6)
  
  # Annotation bar chart by methylation direction
  dmp_anno$direction <- ifelse(dmp_anno$methylation_diff > 0, "Hyper", "Hypo")
  
  anno_direction <- dmp_anno %>%
    group_by(annotation, direction) %>%
    summarise(n = n(), .groups = "drop")
  
  p_anno_bar <- ggplot(anno_direction, aes(x = annotation, y = n, fill = direction)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("Hyper" = "#E74C3C", "Hypo" = "#3498DB")) +
    labs(title = "DMP Distribution by Feature and Direction",
         x = "Genomic Feature", y = "Count", fill = "Direction") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(RESULTS_DIR, "dmp_annotation_direction_bar.pdf"), p_anno_bar,
         width = 10, height = 6)
  
  cat("  DMP annotation plots saved.\n")
}

# ── 5. Annotate DMRs ────────────────────────────────────────────────────────

cat("\n--- Annotating DMRs ---\n")

if (nrow(dmr_results) > 0) {
  dmr_gr <- GRanges(seqnames = dmr_results$chr,
                    ranges = IRanges(start = dmr_results$start, end = dmr_results$end))
  
  dmr_anno <- annotate_with_gff(dmr_gr, genes)
  
  # Add back DMR-specific columns
  dmr_anno$nCG        <- dmr_results$nCG
  dmr_anno$diff_methy <- dmr_results$diff.Methy
  dmr_anno$areaStat   <- dmr_results$areaStat
  
  write.table(dmr_anno, file.path(RESULTS_DIR, "dmrs_annotated.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  anno_summary_dmr <- as.data.frame(table(dmr_anno$annotation))
  colnames(anno_summary_dmr) <- c("Feature", "Count")
  anno_summary_dmr$Percentage <- round(anno_summary_dmr$Count / sum(anno_summary_dmr$Count) * 100, 2)
  
  cat("  DMR annotation distribution:\n")
  print(anno_summary_dmr)
  
  p_dmr_pie <- ggplot(anno_summary_dmr, aes(x = "", y = Count, fill = Feature)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y") +
    scale_fill_brewer(palette = "Set2") +
    labs(title = "DMR Genomic Feature Distribution") +
    theme_void() +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  
  ggsave(file.path(RESULTS_DIR, "dmr_annotation_pie.pdf"), p_dmr_pie,
         width = 8, height = 6)
  cat("  DMR annotation plots saved.\n")
}

# ── 6. Gene-Level Summary ───────────────────────────────────────────────────

cat("\n--- Gene-Level Methylation Summary ---\n")

# Combine DMP annotations to get per-gene methylation changes
if (nrow(dmp_results) > 0 && exists("dmp_anno")) {
  gene_meth_summary <- dmp_anno %>%
    filter(!is.na(nearest_gene)) %>%
    group_by(nearest_gene, annotation) %>%
    summarise(
      n_dmps          = n(),
      n_hyper         = sum(methylation_diff > 0, na.rm = TRUE),
      n_hypo          = sum(methylation_diff < 0, na.rm = TRUE),
      mean_diff       = mean(methylation_diff, na.rm = TRUE),
      max_abs_diff    = max(abs(methylation_diff), na.rm = TRUE),
      min_fdr         = min(fdr, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(min_fdr)
  
  # Add gene type info
  gene_meth_summary <- gene_meth_summary %>%
    left_join(gene_info[, c("gene_id", "gene_type")],
              by = c("nearest_gene" = "gene_id"))
  
  write.table(gene_meth_summary,
              file.path(RESULTS_DIR, "gene_level_methylation_summary.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  cat(sprintf("  Genes with DMPs: %s\n",
              format(length(unique(gene_meth_summary$nearest_gene)), big.mark = ",")))
  
  # Top 20 most affected genes
  top_genes <- gene_meth_summary %>%
    group_by(nearest_gene) %>%
    summarise(total_dmps = sum(n_dmps),
              mean_diff  = mean(mean_diff),
              gene_type  = first(gene_type),
              .groups = "drop") %>%
    arrange(desc(total_dmps)) %>%
    head(20)
  
  cat("\n  Top 20 genes by DMP count:\n")
  print(as.data.frame(top_genes))
  
  write.table(top_genes, file.path(RESULTS_DIR, "top_20_genes_by_dmp_count.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Bar plot of top genes
  p_top_genes <- ggplot(top_genes, aes(x = reorder(nearest_gene, total_dmps),
                                       y = total_dmps, fill = mean_diff)) +
    geom_bar(stat = "identity") +
    scale_fill_gradient2(low = "#3498DB", mid = "white", high = "#E74C3C",
                         midpoint = 0, name = "Mean\nMeth Diff") +
    coord_flip() +
    labs(title = "Top 20 Genes by Number of DMPs",
         x = "Gene", y = "Number of DMPs") +
    theme_minimal()
  
  ggsave(file.path(RESULTS_DIR, "top_genes_dmp_barplot.pdf"), p_top_genes,
         width = 12, height = 8)
}

# ── 7. Functional Enrichment (if clusterProfiler available) ──────────────────

cat("\n--- Functional Enrichment Analysis ---\n")

# For non-model organisms, we attempt enrichment if gene IDs can be mapped.
# If your organism has no org.Db, you can:
#   (a) Use a custom GMT file for enricher()
#   (b) Use GO terms from your GFF functional_note
#   (c) Skip this section

if (has_clusterprof && nrow(dmp_results) > 0 && exists("dmp_anno")) {
  
  # Extract unique gene IDs associated with DMPs
  dmp_genes <- unique(dmp_anno$nearest_gene[!is.na(dmp_anno$nearest_gene)])
  cat(sprintf("  Unique genes associated with DMPs: %s\n", length(dmp_genes)))
  
  # --- Option A: If you have functional annotations in the GFF ---
  # Extract GO terms from the GFF 'Dbxref' or 'Ontology_term' attributes if present
  has_go_terms <- FALSE
  
  if ("Ontology_term" %in% colnames(mcols(genes_gff))) {
    cat("  Found Ontology_term annotations in GFF. Building custom GO mapping...\n")
    
    # Build gene-to-GO mapping
    go_entries <- genes_gff[genes_gff$type == "gene" & !is.na(mcols(genes_gff)$Ontology_term)]
    
    if (length(go_entries) > 0) {
      gene2go <- data.frame(
        gene_id = mcols(go_entries)$ID,
        GO      = mcols(go_entries)$Ontology_term,
        stringsAsFactors = FALSE
      )
      
      # If GO terms are semicolon-separated, split them
      gene2go <- gene2go %>%
        mutate(GO = strsplit(GO, "[;,]")) %>%
        unnest(GO) %>%
        mutate(GO = trimws(GO)) %>%
        filter(GO != "" & grepl("^GO:", GO))
      
      if (nrow(gene2go) > 0) {
        has_go_terms <- TRUE
        cat(sprintf("  Built gene-to-GO mapping: %s genes, %s GO terms\n",
                    length(unique(gene2go$gene_id)),
                    length(unique(gene2go$GO))))
        
        # Custom enrichment with enricher
        # Build TERM2GENE and TERM2NAME
        TERM2GENE <- gene2go[, c("GO", "gene_id")]
        
        enrich_result <- enricher(
          gene      = dmp_genes,
          TERM2GENE = TERM2GENE,
          pvalueCutoff  = 0.05,
          pAdjustMethod = "BH",
          minGSSize = 5
        )
        
        if (!is.null(enrich_result) && nrow(enrich_result@result[enrich_result@result$p.adjust < 0.05, ]) > 0) {
          write.table(enrich_result@result,
                      file.path(RESULTS_DIR, "dmp_custom_go_enrichment.txt"),
                      sep = "\t", quote = FALSE, row.names = FALSE)
          
          # Plot top terms
          top_terms <- enrich_result@result %>%
            filter(p.adjust < 0.05) %>%
            arrange(p.adjust) %>%
            head(20)
          
          p_enrich <- ggplot(top_terms, aes(x = Count, y = reorder(Description, Count),
                                            color = p.adjust)) +
            geom_point(aes(size = Count)) +
            scale_color_gradient(low = "red", high = "blue", name = "Adj. p-value") +
            labs(title = "Top Enriched Terms (DMP-associated Genes)",
                 x = "Gene Count", y = "Term") +
            theme_minimal() +
            theme(axis.text.y = element_text(size = 8))
          
          ggsave(file.path(RESULTS_DIR, "dmp_go_enrichment_dotplot.pdf"),
                 p_enrich, width = 14, height = 10)
          cat("  GO enrichment plots saved.\n")
        } else {
          cat("  No significant GO enrichment found.\n")
        }
      }
    }
  }
  
  if (!has_go_terms) {
    cat("  No GO annotations found in GFF.\n")
    cat("  For non-model organisms, consider:\n")
    cat("    - Adding GO annotations via InterProScan/eggNOG-mapper\n")
    cat("    - Using a custom GMT file with enricher()\n")
    cat("    - Manual curation of top genes\n")
  }
  
  # --- Repeat for DMRs if available ---
  if (nrow(dmr_results) > 0 && exists("dmr_anno") && has_go_terms) {
    dmr_genes <- unique(dmr_anno$nearest_gene[!is.na(dmr_anno$nearest_gene)])
    cat(sprintf("\n  Unique genes associated with DMRs: %s\n", length(dmr_genes)))
    
    enrich_dmr <- enricher(
      gene      = dmr_genes,
      TERM2GENE = TERM2GENE,
      pvalueCutoff  = 0.05,
      pAdjustMethod = "BH",
      minGSSize = 5
    )
    
    if (!is.null(enrich_dmr) && nrow(enrich_dmr@result[enrich_dmr@result$p.adjust < 0.05, ]) > 0) {
      write.table(enrich_dmr@result,
                  file.path(RESULTS_DIR, "dmr_custom_go_enrichment.txt"),
                  sep = "\t", quote = FALSE, row.names = FALSE)
      cat("  DMR GO enrichment results saved.\n")
    } else {
      cat("  No significant GO enrichment for DMR genes.\n")
    }
  }
} else {
  cat("  clusterProfiler not available or no DMPs found. Skipping enrichment.\n")
}

# ── 8. Chromosome-Level Summary ─────────────────────────────────────────────

cat("\n--- Chromosome-Level DMP/DMR Distribution ---\n")

if (nrow(dmp_results) > 0) {
  chr_summary <- dmp_results %>%
    mutate(
      direction = ifelse(diff > 0, "Hyper", "Hypo"),
      chr_num   = as.numeric(gsub("chr", "", chr))
    ) %>%
    filter(!is.na(chr_num)) %>%
    group_by(chr, chr_num, direction) %>%
    summarise(count = n(), .groups = "drop") %>%
    arrange(chr_num)
  
  p_chr <- ggplot(chr_summary, aes(x = reorder(chr, chr_num), y = count, fill = direction)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("Hyper" = "#E74C3C", "Hypo" = "#3498DB")) +
    labs(title = "DMP Distribution by Chromosome",
         x = "Chromosome", y = "Number of DMPs", fill = "Direction") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(RESULTS_DIR, "dmp_chromosome_distribution.pdf"), p_chr,
         width = 14, height = 6)
  cat("  Chromosome distribution plot saved.\n")
}

# ── 9. Gene Type Methylation Analysis ────────────────────────────────────────

cat("\n--- Methylation by Gene Type ---\n")

if (exists("gene_meth_summary")) {
  gene_type_meth <- gene_meth_summary %>%
    filter(!is.na(gene_type)) %>%
    group_by(gene_type) %>%
    summarise(
      n_genes       = n_distinct(nearest_gene),
      total_dmps    = sum(n_dmps),
      mean_meth_diff = mean(mean_diff, na.rm = TRUE),
      pct_hyper     = sum(n_hyper) / sum(n_dmps) * 100,
      pct_hypo      = sum(n_hypo)  / sum(n_dmps) * 100,
      .groups = "drop"
    )
  
  cat("  Methylation changes by gene type:\n")
  print(as.data.frame(gene_type_meth))
  
  write.table(gene_type_meth,
              file.path(RESULTS_DIR, "methylation_by_gene_type.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

# ── 10. Final Summary Report ────────────────────────────────────────────────

cat("\n")
cat("================================================================\n")
cat("           METHYLATION ANALYSIS COMPLETE                        \n")
cat("================================================================\n")
cat("Results directory:", RESULTS_DIR, "\n\n")

cat("Part 1 outputs (QC & EDA):\n")
cat("  - bsseq_object.rds\n")
cat("  - global_methylation_stats.txt\n")
cat("  - sample_correlation_heatmap.pdf\n")
cat("  - pca_analysis.pdf\n")
cat("  - methylation_density.pdf\n\n")

cat("Part 2 outputs (Differential Methylation):\n")
cat("  - differentially_methylated_positions.txt\n")
cat("  - differentially_methylated_regions.txt\n")
cat("  - analysis_summary.txt\n")
cat("  - volcano_plot_dmps.pdf\n")
cat("  - manhattan_plot_methylation_differences.pdf\n")
cat("  - genome_wide_methylation.pdf\n")
cat("  - dmp_level_heatmap.pdf / dmr_level_heatmap.pdf\n\n")

cat("Part 3 outputs (Annotation & Enrichment):\n")
cat("  - dmps_annotated.txt / dmrs_annotated.txt\n")
cat("  - dmp_annotation_pie.pdf / dmr_annotation_pie.pdf\n")
cat("  - dmp_annotation_direction_bar.pdf\n")
cat("  - gene_level_methylation_summary.txt\n")
cat("  - top_20_genes_by_dmp_count.txt\n")
cat("  - top_genes_dmp_barplot.pdf\n")
cat("  - dmp_chromosome_distribution.pdf\n")
cat("  - methylation_by_gene_type.txt\n")
cat("  - *_go_enrichment.txt (if GO annotations available)\n")
cat("================================================================\n")