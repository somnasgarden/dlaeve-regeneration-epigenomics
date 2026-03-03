#!/usr/bin/env Rscript
# =============================================================================
# Sox19a Gviz Visualization, TF Methylation ChIPseeker & Motif Analysis
# =============================================================================
# Part 1: Sox19a locus visualization with Gviz (gene model + DMP tracks)
# Part 2: Motif scanning around DMP clusters (TFBSTools/JASPAR)
#         → Can we tell if a region is an enhancer or silencer?
# Part 3: TF DMPs/DMRs with ChIPseeker pie charts and regions
# Part 4: Enrichment analysis of differentially methylated TFs
#
# Key question: are the intergenic DMPs near Sox19a in an enhancer or silencer?
#
# Libraries:
#   Gviz        — Genomic track visualization
#   TFBSTools   — TF binding site analysis
#   JASPAR2024  — Motif database (or JASPAR2022)
#   motifmatchr — Fast motif scanning
#   Biostrings  — Sequence manipulation
#   ChIPseeker  — Genomic annotation + pie charts
#
# Input:  results/01_methylation/tables/ (annotated DMPs/DMRs)
#         DATA/prediction_result.txt (DeepTFactor)
#         DATA/derLaeGenome_namesDlasi_v2.fasta.functional_note.pseudo_label.gff
#         DATA/derLaeGenome_eviann_annotations.tsv
# Output: results/21_sox19a_tf_motifs/
# =============================================================================

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(scales)
})

options(stringsAsFactors = FALSE)
options(scipen = 999)

theme_set(theme_minimal(base_size = 12) +
            theme(plot.title = element_text(size = 14, face = "bold")))

# ── Paths ────────────────────────────────────────────────────────────────────

DATA_DIR <- NULL
for (d in c("/mnt/c/Users/rafae/Projects/DATA",
            "C:/Users/rafae/Projects/DATA",
            "/mnt/data/alfredvar/rlopezt/Preliminary")) {
  if (dir.exists(d)) { DATA_DIR <- d; break }
}
if (is.null(DATA_DIR)) stop("Cannot find DATA directory")

RESULTS_BASE <- "results/01_methylation/tables"
OUTPUT_DIR   <- "results/21_sox19a_tf_motifs"
dir.create(file.path(OUTPUT_DIR, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "plots"),  recursive = TRUE, showWarnings = FALSE)

GFF_FILE   <- file.path(DATA_DIR, "derLaeGenome_namesDlasi_v2.fasta.functional_note.pseudo_label.gff")
TF_FILE    <- file.path(DATA_DIR, "prediction_result.txt")
ANNOT_FILE <- file.path(DATA_DIR, "derLaeGenome_eviann_annotations.tsv")

save_both <- function(p, name, w = 10, h = 7) {
  ggsave(file.path(OUTPUT_DIR, "plots", paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(OUTPUT_DIR, "plots", paste0(name, ".png")), p, width = w, height = h, dpi = 300)
  cat("  Saved:", name, "\n")
}

cat("=== Sox19a Gviz, TF ChIPseeker & Motif Analysis ===\n")
cat("DATA_DIR:", DATA_DIR, "\n\n")

# Chromosome filter — chr1-31 only (exclude scaffolds)
keep_chr <- paste0("chr", 1:31)

# =============================================================================
# LOAD COMMON DATA
# =============================================================================

cat("--- Loading data ---\n")

dmps <- fread(file.path(RESULTS_BASE, "dmps_annotated.txt"))
dmrs <- fread(file.path(RESULTS_BASE, "dmrs_annotated.txt"))
dmps_chip <- fread("results/01_methylation/dmps_chipseeker_annotated.txt")
gff <- import(GFF_FILE)

# Filter to chr1-31
dmps <- dmps[seqnames %in% keep_chr]
dmrs <- dmrs[seqnames %in% keep_chr]
if ("seqnames" %in% colnames(dmps_chip)) dmps_chip <- dmps_chip[seqnames %in% keep_chr]
gff <- gff[as.character(seqnames(gff)) %in% keep_chr]
cat("  Filtered to chr1-31\n")

# Gene annotations
annot <- fread(ANNOT_FILE, header = FALSE, col.names = c("gene_id", "gene_name", "description"))
annot_unique <- annot[!duplicated(gene_id)]

# DeepTFactor
tf_pred <- fread(TF_FILE)
tf_pred[, gene_id := sub("-mRNA-.*$", "", sequence_ID)]
tf_genes <- tf_pred[, .(is_tf = any(prediction == TRUE),
                         max_score = max(score)), by = gene_id]

# DEGs
deg_file <- "results/02_rnaseq/Part2_DEGs/data/DEresults_amputated_vs_control.csv"
if (file.exists(deg_file)) {
  degs <- fread(deg_file)
  if ("V1" %in% colnames(degs)) setnames(degs, "V1", "gene_id")
}

# STRING enrichment DB
enrichment_db <- NULL
for (epath in c("results/01_methylation/protein.enrichment.terms.v12.0.txt",
                "results/01_methylation/tables/protein.enrichment.terms.v12.0.txt",
                file.path(DATA_DIR, "protein.enrichment.terms.v12.0.txt"))) {
  if (file.exists(epath)) {
    enrichment_db <- fread(epath, quote = "")
    cat("  STRING DB loaded\n")
    break
  }
}

cat("  DMPs:", nrow(dmps), "| DMRs:", nrow(dmrs), "| ChIPseeker DMPs:", nrow(dmps_chip), "\n")
cat("  GFF features:", length(gff), "| TF genes:", sum(tf_genes$is_tf), "\n\n")


# =============================================================================
# PART 1: SOX19a GVIZ VISUALIZATION
# =============================================================================

cat("================================================================\n")
cat(" PART 1: Sox19a Locus Visualization (Gviz)                      \n")
cat("================================================================\n\n")

# Try loading Gviz
has_gviz <- requireNamespace("Gviz", quietly = TRUE)
if (!has_gviz) {
  cat("Installing Gviz...\n")
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  BiocManager::install("Gviz", update = FALSE, ask = FALSE)
  has_gviz <- requireNamespace("Gviz", quietly = TRUE)
}

sox_id <- "LOC_00002600"

# Get Sox19a gene coordinates from GFF
sox_gene <- gff[gff$type == "gene" & grepl(sox_id, gff$ID)]
if (length(sox_gene) == 0) {
  sox_gene <- gff[gff$type == "gene" & grepl(sox_id, as.character(gff$Name))]
}
if (length(sox_gene) == 0) {
  sox_gene <- gff[grepl(sox_id, as.character(gff$ID))]
  sox_gene <- sox_gene[sox_gene$type == "gene"]
}

cat("  Sox19a gene:", as.character(seqnames(sox_gene)), ":",
    start(sox_gene), "-", end(sox_gene), "\n")

# Get Sox19a ChIPseeker DMPs
sox_chip <- as.data.table(dmps_chip)[geneId == sox_id]
sox_dmps <- merge(sox_chip, dmps, by = c("seqnames", "start"), suffixes = c("_chip", "_annot"))

cat("  Sox19a-associated DMPs:", nrow(sox_dmps), "\n")
if (nrow(sox_dmps) > 0) {
  cat("  DMP positions and annotations:\n")
  for (i in 1:nrow(sox_dmps)) {
    cat(sprintf("    %s:%d  distTSS=%d  diff=%+.3f  annotation=%s\n",
                sox_dmps$seqnames[i], sox_dmps$start[i],
                sox_dmps$distanceToTSS[i],
                sox_dmps$methylation_diff[i],
                sox_dmps$annotation_chip[i]))
  }
}

# About the "Downstream" annotation at negative distance:
# ChIPseeker annotates relative to the NEAREST gene feature
# "Downstream" = downstream of TES (3' end of gene body)
# distanceToTSS is ALWAYS from the TSS (5' end)
# So a DMP can be "Downstream" of gene X but with negative distTSS
# if it's upstream of the TSS of the NEAREST gene in the other direction
# For Sox19a: DMPs are ~4-5kb upstream → they are actually in the
# intergenic region BETWEEN Sox19a and its upstream neighbor

if (has_gviz) {
  library(Gviz)

  sox_chr <- as.character(seqnames(sox_gene))[1]
  sox_start <- start(sox_gene)[1]
  sox_end <- end(sox_gene)[1]

  # Window: gene ± 10kb to see the intergenic DMPs
  plot_start <- sox_start - 10000
  plot_end <- sox_end + 5000

  # Track 1: Genome axis
  gtrack <- GenomeAxisTrack(name = "Position")

  # Track 2: Gene models from GFF (Sox19a region)
  genes_in_region <- gff[seqnames(gff) == sox_chr &
                           start(gff) >= plot_start - 5000 &
                           end(gff) <= plot_end + 5000 &
                           gff$type %in% c("gene", "mRNA", "exon", "CDS")]

  if (length(genes_in_region) > 0) {
    # Build a GeneRegionTrack
    gene_track <- GeneRegionTrack(
      genes_in_region,
      genome = "Dlaeve",
      chromosome = sox_chr,
      name = "Genes",
      showId = TRUE,
      geneSymbol = TRUE,
      fill = "#2E86C1",
      col = "#2E86C1",
      fontcolor.group = "black",
      background.title = "#2E86C1"
    )
  } else {
    # Fallback: annotation track from coordinates
    gene_track <- AnnotationTrack(
      start = sox_start, end = sox_end,
      chromosome = sox_chr,
      strand = as.character(strand(sox_gene))[1],
      id = "Sox19a",
      name = "Genes",
      fill = "#2E86C1",
      showId = TRUE,
      background.title = "#2E86C1"
    )
  }

  # Track 3: DMP positions with methylation values
  if (nrow(sox_dmps) > 0) {
    dmp_gr_sox <- GRanges(
      seqnames = sox_dmps$seqnames,
      ranges = IRanges(start = sox_dmps$start, width = 1)
    )

    dmp_track <- DataTrack(
      range = dmp_gr_sox,
      data = sox_dmps$methylation_diff,
      chromosome = sox_chr,
      name = "Methylation\nDiff",
      type = c("p", "h"),
      col = ifelse(sox_dmps$methylation_diff > 0, "#E74C3C", "#3498DB"),
      baseline = 0,
      ylim = c(-1, 1),
      background.title = "#E74C3C",
      legend = FALSE
    )

    # Track 4: Annotation regions
    region_colors <- c("Intergenic" = "#E74C3C", "Downstream" = "#E67E22",
                       "Intron" = "#27AE60", "Exon" = "#2E86C1",
                       "Promoter" = "#8E44AD")

    region_labels <- sub("\\s*\\(.*", "", sox_dmps$annotation_chip)
    region_labels[grepl("Intron", region_labels)] <- "Intron"
    region_labels[grepl("Exon", region_labels)] <- "Exon"
    region_labels[grepl("Downstream", region_labels)] <- "Downstream"
    region_labels[grepl("Promoter", region_labels)] <- "Promoter"

    region_track <- AnnotationTrack(
      start = sox_dmps$start,
      end = sox_dmps$start,
      chromosome = sox_chr,
      feature = region_labels,
      name = "Region",
      fill = region_colors[region_labels],
      showId = FALSE,
      background.title = "#95A5A6"
    )

    # Plot
    png(file.path(OUTPUT_DIR, "plots", "01_sox19a_gviz_locus.png"),
        width = 1400, height = 800, res = 150)
    plotTracks(list(gtrack, gene_track, dmp_track, region_track),
              from = plot_start, to = plot_end,
              chromosome = sox_chr,
              main = sprintf("Sox19a Locus — DMPs in Intergenic Region (%s:%s-%s)",
                             sox_chr,
                             format(plot_start, big.mark = ","),
                             format(plot_end, big.mark = ",")),
              cex.main = 1.2)
    dev.off()

    pdf(file.path(OUTPUT_DIR, "plots", "01_sox19a_gviz_locus.pdf"),
        width = 14, height = 8)
    plotTracks(list(gtrack, gene_track, dmp_track, region_track),
              from = plot_start, to = plot_end,
              chromosome = sox_chr,
              main = sprintf("Sox19a Locus — DMPs in Intergenic Region"))
    dev.off()
    cat("  Saved: 01_sox19a_gviz_locus\n")
  }
} else {
  cat("  Gviz not available. Using ggplot fallback.\n")
}

# Interpretation: enhancer or silencer?
cat("\n  ENHANCER vs SILENCER INTERPRETATION for Sox19a:\n")
cat("  ─────────────────────────────────────────────\n")

if (nrow(sox_dmps) > 0) {
  n_hyper <- sum(sox_dmps$methylation_diff > 0)
  n_hypo <- sum(sox_dmps$methylation_diff < 0)
  mean_diff <- mean(sox_dmps$methylation_diff)

  sox_expr_fc <- degs[gene_id == sox_id]$log2FoldChange

  cat(sprintf("  Gene expression: log2FC = %+.3f (UPREGULATED in amputated)\n", sox_expr_fc))
  cat(sprintf("  DMP direction: %d hyper, %d hypo (mean diff = %+.3f)\n", n_hyper, n_hypo, mean_diff))
  cat(sprintf("  DMPs are in: INTERGENIC region ~4-5kb upstream of TSS\n\n"))

  if (sox_expr_fc > 0 && mean_diff > 0) {
    cat("  PREDICTION: SILENCER region undergoing INACTIVATION\n")
    cat("    Reasoning: Hypermethylation of a silencer → silencer inactivation\n")
    cat("              → gene de-repression → UPREGULATION\n")
    cat("    Alternative: Hypermethylation attracting methyl-binding activators\n")
  } else if (sox_expr_fc > 0 && mean_diff < 0) {
    cat("  PREDICTION: ENHANCER region undergoing ACTIVATION\n")
    cat("    Reasoning: Hypomethylation of enhancer → TF binding enabled\n")
    cat("              → enhancer activation → UPREGULATION\n")
  } else if (sox_expr_fc < 0 && mean_diff > 0) {
    cat("  PREDICTION: ENHANCER region undergoing INACTIVATION\n")
    cat("    Reasoning: Hypermethylation of enhancer → TF binding blocked\n")
    cat("              → enhancer silencing → DOWNREGULATION\n")
  } else {
    cat("  PREDICTION: SILENCER region undergoing ACTIVATION\n")
    cat("    Reasoning: Hypomethylation → silencer reactivation → DOWNREGULATION\n")
  }

  cat("\n  TO CONFIRM (experiments needed):\n")
  cat("    1. ATAC-seq: Check chromatin accessibility at DMP cluster\n")
  cat("    2. H3K4me1/H3K27ac ChIP-seq: Enhancer histone marks\n")
  cat("    3. H3K27me3 ChIP-seq: Silencer/Polycomb marks\n")
  cat("    4. Reporter assay: Clone region → test enhancer/silencer activity\n")
  cat("    5. Motif scanning: Look for TF binding sites (Part 2 below)\n\n")
}


# =============================================================================
# PART 2: MOTIF SCANNING AROUND DMP CLUSTERS
# =============================================================================

cat("================================================================\n")
cat(" PART 2: Motif Scanning in DMP Regions                          \n")
cat("================================================================\n\n")

# Try loading TFBSTools + JASPAR
has_tfbs <- requireNamespace("TFBSTools", quietly = TRUE)
has_jaspar <- requireNamespace("JASPAR2024", quietly = TRUE) ||
              requireNamespace("JASPAR2022", quietly = TRUE) ||
              requireNamespace("JASPAR2020", quietly = TRUE)
has_biostrings <- requireNamespace("Biostrings", quietly = TRUE)
has_bsgenome <- requireNamespace("BSgenome.Dlaeve.NCBI.dlgm", quietly = TRUE)

if (!has_tfbs) {
  cat("  Installing TFBSTools...\n")
  BiocManager::install("TFBSTools", update = FALSE, ask = FALSE)
  has_tfbs <- requireNamespace("TFBSTools", quietly = TRUE)
}

# Try JASPAR versions in order
jaspar_pkg <- NULL
for (jp in c("JASPAR2024", "JASPAR2022", "JASPAR2020")) {
  if (requireNamespace(jp, quietly = TRUE)) {
    jaspar_pkg <- jp
    break
  }
}
if (is.null(jaspar_pkg)) {
  cat("  Installing JASPAR2024...\n")
  BiocManager::install("JASPAR2024", update = FALSE, ask = FALSE)
  if (requireNamespace("JASPAR2024", quietly = TRUE)) jaspar_pkg <- "JASPAR2024"
}

if (has_tfbs && !is.null(jaspar_pkg) && has_bsgenome) {
  library(TFBSTools)
  library(Biostrings)

  # Load JASPAR motifs
  jaspar_db <- do.call(library, list(jaspar_pkg))
  cat(sprintf("  Using %s for motif database\n", jaspar_pkg))

  # Get invertebrate + general TF motifs
  # Since D. laeve is a mollusk, use all metazoan motifs
  opts <- list(tax_group = "metazoa", all_versions = FALSE)
  pwm_list <- getMatrixSet(get(paste0(jaspar_pkg, "_db")), opts)
  cat(sprintf("  JASPAR metazoan motifs: %d\n", length(pwm_list)))

  # Load BSgenome for D. laeve
  library(BSgenome.Dlaeve.NCBI.dlgm)
  genome <- BSgenome.Dlaeve.NCBI.dlgm

  # Define regions to scan: Sox19a intergenic DMP cluster
  if (nrow(sox_dmps) > 0) {
    sox_chr <- sox_dmps$seqnames[1]
    scan_start <- min(sox_dmps$start) - 500
    scan_end <- max(sox_dmps$start) + 500

    cat(sprintf("\n  Scanning Sox19a DMP region: %s:%d-%d (%d bp)\n",
                sox_chr, scan_start, scan_end, scan_end - scan_start))

    # Get sequence
    seq_region <- getSeq(genome, sox_chr, start = scan_start, end = scan_end)

    # Scan for motifs
    cat("  Scanning for TF binding motifs...\n")

    # Convert PWMs to scoring matrices
    motif_hits <- data.frame()
    threshold <- "80%"  # Match at >= 80% of max score

    for (i in seq_along(pwm_list)) {
      pwm <- toPWM(pwm_list[[i]])
      name_motif <- name(pwm_list[[i]])

      # Search both strands
      hits <- searchSeq(pwm, seq_region, min.score = threshold, strand = "*")

      if (length(hits@views) > 0) {
        for (h in seq_along(hits@views)) {
          motif_hits <- rbind(motif_hits, data.frame(
            motif_id = ID(pwm_list[[i]]),
            motif_name = name_motif,
            start = start(hits@views)[h] + scan_start - 1,
            end = end(hits@views)[h] + scan_start - 1,
            strand = as.character(strand(hits@views))[h],
            score = round(score(hits@views)[h], 4),
            rel_score = round(relScore(hits@views)[h], 4),
            stringsAsFactors = FALSE
          ))
        }
      }
    }

    if (nrow(motif_hits) > 0) {
      motif_hits <- motif_hits %>% arrange(desc(rel_score))

      cat(sprintf("  Found %d motif hits in Sox19a DMP region:\n", nrow(motif_hits)))
      top_motifs <- head(motif_hits, 20)
      for (i in 1:nrow(top_motifs)) {
        cat(sprintf("    %-15s  pos=%d-%d  strand=%s  score=%.1f%%\n",
                    top_motifs$motif_name[i],
                    top_motifs$start[i], top_motifs$end[i],
                    top_motifs$strand[i],
                    top_motifs$rel_score[i] * 100))
      }

      fwrite(motif_hits, file.path(OUTPUT_DIR, "tables", "sox19a_motif_hits.tsv"), sep = "\t")

      # Classify motifs as activator/repressor (based on common knowledge)
      activator_motifs <- c("SOX", "POU", "KLF", "ETS", "GATA", "CEBP", "AP1",
                            "FOXA", "HNF", "NRF", "CTCF", "SP1")
      repressor_motifs <- c("REST", "SNAI", "ZEB", "GFI", "BCL6",
                            "PRDM", "HES", "HEY")

      motif_hits$category <- "Unknown"
      for (am in activator_motifs) {
        motif_hits$category[grepl(am, motif_hits$motif_name, ignore.case = TRUE)] <- "Activator"
      }
      for (rm in repressor_motifs) {
        motif_hits$category[grepl(rm, motif_hits$motif_name, ignore.case = TRUE)] <- "Repressor"
      }

      cat("\n  Motif categories:\n")
      print(table(motif_hits$category))

      n_activator <- sum(motif_hits$category == "Activator")
      n_repressor <- sum(motif_hits$category == "Repressor")
      cat(sprintf("\n  Enhancer evidence (activator TF motifs): %d\n", n_activator))
      cat(sprintf("  Silencer evidence (repressor TF motifs): %d\n", n_repressor))

      if (n_activator > n_repressor) {
        cat("  → MOTIF ANALYSIS SUPPORTS ENHANCER\n")
      } else if (n_repressor > n_activator) {
        cat("  → MOTIF ANALYSIS SUPPORTS SILENCER\n")
      } else {
        cat("  → MOTIF ANALYSIS INCONCLUSIVE (need experimental validation)\n")
      }
    } else {
      cat("  No motif hits found at 80% threshold. Try lowering to 75%.\n")
    }
  }
} else {
  cat("  Missing packages for motif scanning:\n")
  cat(sprintf("    TFBSTools: %s\n", ifelse(has_tfbs, "OK", "MISSING")))
  cat(sprintf("    JASPAR: %s\n", ifelse(!is.null(jaspar_pkg), jaspar_pkg, "MISSING")))
  cat(sprintf("    BSgenome.Dlaeve: %s\n", ifelse(has_bsgenome, "OK", "MISSING")))
  cat("  Install via BiocManager::install(c('TFBSTools', 'JASPAR2024'))\n")
  cat("  BSgenome.Dlaeve.NCBI.dlgm must be installed from local package\n\n")

  # Alternative: known enhancer/silencer features
  cat("  ALTERNATIVE ENHANCER/SILENCER INDICATORS:\n")
  cat("  ─────────────────────────────────────────\n")
  cat("  1. CpG density: Enhancers tend to be CpG-poor\n")
  cat("  2. Distance from TSS: 1-100kb typical for enhancers\n")
  cat("  3. Conservation: Enhancers tend to be conserved\n")
  cat("  4. DMP clustering: Clustered DMPs suggest regulatory region\n")
  cat("  5. Direction logic: See Part 1 interpretation above\n")
}


# =============================================================================
# PART 3: TF DMPs/DMRs WITH CHIPSEEKER PIE CHARTS
# =============================================================================

cat("\n================================================================\n")
cat(" PART 3: TF Methylation — ChIPseeker Regions & Pie Charts       \n")
cat("================================================================\n\n")

# Get TF gene IDs
tf_gene_ids <- tf_genes[is_tf == TRUE]$gene_id
nontf_gene_ids <- tf_genes[is_tf == FALSE]$gene_id

cat("  Predicted TFs:", length(tf_gene_ids), "\n")
cat("  Non-TFs:", length(nontf_gene_ids), "\n\n")

# ── ChIPseeker-style annotation for TF-associated DMPs ──
cat("ChIPseeker regions for DMPs near TF genes:\n")

# DMPs from annotated file (uses our annotation: Gene Body, Intergenic, Promoter)
dmps_tf <- dmps[nearest_gene %in% tf_gene_ids]
dmps_nontf <- dmps[nearest_gene %in% nontf_gene_ids]

cat(sprintf("  DMPs near TF genes: %d\n", nrow(dmps_tf)))
cat(sprintf("  DMPs near Non-TF genes: %d\n\n", nrow(dmps_nontf)))

# Get ChIPseeker detailed annotation for TF DMPs
chip_tf <- as.data.table(dmps_chip)[geneId %in% tf_gene_ids]
chip_nontf <- as.data.table(dmps_chip)[geneId %in% nontf_gene_ids]

# Parse ChIPseeker annotation into categories
parse_chipseeker <- function(annotation) {
  region <- case_when(
    grepl("Promoter", annotation) ~ "Promoter",
    grepl("5' UTR", annotation) ~ "5' UTR",
    grepl("3' UTR", annotation) ~ "3' UTR",
    grepl("Exon", annotation) ~ "Exon",
    grepl("Intron", annotation) ~ "Intron",
    grepl("Downstream", annotation) ~ "Downstream",
    grepl("Intergenic", annotation) ~ "Intergenic",
    TRUE ~ "Other"
  )
  return(region)
}

# TF DMP regions (ChIPseeker-style)
if (nrow(chip_tf) > 0) {
  chip_tf[, region := parse_chipseeker(annotation)]
  tf_region_counts <- chip_tf[, .N, by = region]
  tf_region_counts[, pct := round(100 * N / sum(N), 1)]
  tf_region_counts[, category := "TF genes"]

  cat("  TF gene DMP regions (ChIPseeker):\n")
  for (i in 1:nrow(tf_region_counts)) {
    cat(sprintf("    %-12s %5d (%5.1f%%)\n",
                tf_region_counts$region[i],
                tf_region_counts$N[i],
                tf_region_counts$pct[i]))
  }
}

if (nrow(chip_nontf) > 0) {
  chip_nontf[, region := parse_chipseeker(annotation)]
  nontf_region_counts <- chip_nontf[, .N, by = region]
  nontf_region_counts[, pct := round(100 * N / sum(N), 1)]
  nontf_region_counts[, category := "Non-TF genes"]

  cat("\n  Non-TF gene DMP regions (ChIPseeker):\n")
  for (i in 1:nrow(nontf_region_counts)) {
    cat(sprintf("    %-12s %5d (%5.1f%%)\n",
                nontf_region_counts$region[i],
                nontf_region_counts$N[i],
                nontf_region_counts$pct[i]))
  }
}

# ── Pie charts ──
cat("\n  Generating pie charts...\n")

region_palette <- c("Promoter" = "#8E44AD", "5' UTR" = "#9B59B6",
                     "3' UTR" = "#3498DB", "Exon" = "#2ECC71",
                     "Intron" = "#27AE60", "Downstream" = "#E67E22",
                     "Intergenic" = "#E74C3C", "Other" = "#95A5A6")

if (nrow(chip_tf) > 0) {
  # TF pie chart
  tf_pie_data <- tf_region_counts[order(-N)]
  tf_pie_data[, label := sprintf("%s\n%d (%.1f%%)", region, N, pct)]

  p_tf_pie <- ggplot(tf_pie_data, aes(x = "", y = N, fill = region)) +
    geom_col(width = 1, color = "white") +
    coord_polar("y") +
    scale_fill_manual(values = region_palette) +
    geom_text(aes(label = label),
              position = position_stack(vjust = 0.5), size = 3) +
    labs(title = sprintf("DMP Regions for TF Genes (n = %d DMPs)", nrow(chip_tf)),
         fill = "Region") +
    theme_void(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
  save_both(p_tf_pie, "01_tf_dmp_pie_chart", 10, 8)

  # Non-TF pie chart
  nontf_pie_data <- nontf_region_counts[order(-N)]
  nontf_pie_data[, label := sprintf("%s\n%d (%.1f%%)", region, N, pct)]

  p_nontf_pie <- ggplot(nontf_pie_data, aes(x = "", y = N, fill = region)) +
    geom_col(width = 1, color = "white") +
    coord_polar("y") +
    scale_fill_manual(values = region_palette) +
    geom_text(aes(label = label),
              position = position_stack(vjust = 0.5), size = 3) +
    labs(title = sprintf("DMP Regions for Non-TF Genes (n = %d DMPs)", nrow(chip_nontf)),
         fill = "Region") +
    theme_void(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
  save_both(p_nontf_pie, "02_nontf_dmp_pie_chart", 10, 8)

  # Side-by-side comparison bar chart
  all_region <- rbind(tf_region_counts, nontf_region_counts)

  p_region_compare <- ggplot(all_region, aes(x = region, y = pct, fill = category)) +
    geom_col(position = "dodge", width = 0.7) +
    scale_fill_manual(values = c("TF genes" = "#E74C3C", "Non-TF genes" = "#3498DB")) +
    labs(title = "DMP Region Distribution: TF vs Non-TF Genes",
         subtitle = "ChIPseeker annotation categories",
         x = "Genomic Region", y = "Percentage of DMPs (%)",
         fill = "") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_both(p_region_compare, "03_tf_nontf_region_comparison", 10, 6)

  fwrite(all_region, file.path(OUTPUT_DIR, "tables", "tf_nontf_dmp_regions.tsv"), sep = "\t")
}

# ── DMR analysis for TFs ──
cat("\nDMR analysis for TF genes:\n")
dmrs_tf <- dmrs[nearest_gene %in% tf_gene_ids]
dmrs_nontf <- dmrs[nearest_gene %in% nontf_gene_ids]

cat(sprintf("  DMRs near TF genes: %d\n", nrow(dmrs_tf)))
cat(sprintf("  DMRs near Non-TF genes: %d\n", nrow(dmrs_nontf)))

if (nrow(dmrs_tf) > 0) {
  tf_dmr_regions <- dmrs_tf[, .N, by = annotation]
  tf_dmr_regions[, pct := round(100 * N / sum(N), 1)]
  cat("  TF DMR regions:\n")
  for (i in 1:nrow(tf_dmr_regions)) {
    cat(sprintf("    %-15s %3d (%5.1f%%)\n",
                tf_dmr_regions$annotation[i], tf_dmr_regions$N[i], tf_dmr_regions$pct[i]))
  }
}


# =============================================================================
# PART 4: ENRICHMENT OF DIFFERENTIALLY METHYLATED TFs
# =============================================================================

cat("\n================================================================\n")
cat(" PART 4: Enrichment of Differentially Methylated TFs            \n")
cat("================================================================\n\n")

# TFs with DMPs
tf_with_dmp <- tf_genes[is_tf == TRUE & gene_id %in% dmps$nearest_gene]$gene_id
tf_with_dmr <- tf_genes[is_tf == TRUE & gene_id %in% dmrs$nearest_gene]$gene_id
tf_methylated <- unique(c(tf_with_dmp, tf_with_dmr))

cat(sprintf("  TFs with DMP or DMR: %d / %d\n", length(tf_methylated), length(tf_gene_ids)))

# Enrichment analysis of methylated TFs
if (!is.null(enrichment_db) && length(tf_methylated) >= 5) {

  background <- unique(enrichment_db$`#string_protein_id`)

  cat("\n  Running enrichment on methylated TFs...\n")

  for (cat_name in c("Process", "Function", "Component", "KEGG", "Reactome")) {
    db_sub <- enrichment_db[enrichment_db$category == cat_name, ]
    terms <- unique(db_sub$term)

    results <- data.frame()
    for (term in terms) {
      term_genes <- unique(db_sub$`#string_protein_id`[db_sub$term == term])
      a <- length(intersect(tf_methylated, term_genes))
      b <- length(tf_methylated) - a
      c_val <- length(intersect(background, term_genes)) - a
      d <- length(background) - a - b - c_val

      if (a >= 2 && d > 0) {
        ft <- fisher.test(matrix(c(a, b, c_val, d), nrow = 2), alternative = "greater")
        desc <- db_sub$description[db_sub$term == term][1]
        results <- rbind(results, data.frame(
          term = term, description = desc,
          n_overlap = a, n_term = length(term_genes),
          n_query = length(tf_methylated),
          odds_ratio = round(ft$estimate, 3),
          p_value = ft$p.value
        ))
      }
    }

    if (nrow(results) > 0) {
      results$padj <- p.adjust(results$p_value, method = "BH")
      results <- results[order(results$p_value), ]
      sig <- results[results$padj < 0.05, ]

      if (nrow(sig) > 0) {
        cat(sprintf("\n  %s — %d significant terms (padj < 0.05):\n", cat_name, nrow(sig)))
        top_sig <- head(sig, 10)
        for (j in 1:nrow(top_sig)) {
          cat(sprintf("    %s (n=%d, OR=%.2f, padj=%.2e)\n",
                      top_sig$description[j], top_sig$n_overlap[j],
                      top_sig$odds_ratio[j], top_sig$padj[j]))
        }
      }

      fwrite(results, file.path(OUTPUT_DIR, "tables",
                                 paste0("tf_methylated_enrichment_", cat_name, ".tsv")),
             sep = "\t")
    }
  }

  # Also: enrichment of TFs with DMP in specific regions
  cat("\n  Region-specific TF enrichment:\n")

  for (region in c("Gene Body", "Intergenic", "Promoter")) {
    dmps_region <- dmps[annotation == region]
    tf_in_region <- tf_genes[is_tf == TRUE & gene_id %in% dmps_region$nearest_gene]$gene_id

    if (length(tf_in_region) >= 3) {
      cat(sprintf("    TFs with %s DMPs: %d\n", region, length(tf_in_region)))

      db_sub <- enrichment_db[enrichment_db$category == "Process", ]
      terms <- unique(db_sub$term)
      results <- data.frame()

      for (term in terms) {
        term_genes <- unique(db_sub$`#string_protein_id`[db_sub$term == term])
        a <- length(intersect(tf_in_region, term_genes))
        if (a >= 2) {
          b <- length(tf_in_region) - a
          c_val <- length(intersect(background, term_genes)) - a
          d <- length(background) - a - b - c_val
          if (d > 0) {
            ft <- fisher.test(matrix(c(a, b, c_val, d), nrow = 2), alternative = "greater")
            results <- rbind(results, data.frame(
              region = region, term = term,
              description = db_sub$description[db_sub$term == term][1],
              n_overlap = a, odds_ratio = round(ft$estimate, 3),
              p_value = ft$p.value
            ))
          }
        }
      }

      if (nrow(results) > 0) {
        results$padj <- p.adjust(results$p_value, method = "BH")
        sig <- results[results$padj < 0.05, ]
        if (nrow(sig) > 0) {
          cat(sprintf("      Top enriched: %s (OR=%.2f, padj=%.2e)\n",
                      sig$description[1], sig$odds_ratio[1], sig$padj[1]))
        }
        fwrite(results, file.path(OUTPUT_DIR, "tables",
                                   paste0("tf_enrichment_", gsub(" ", "_", region), ".tsv")),
               sep = "\t")
      }
    }
  }
} else {
  if (is.null(enrichment_db)) {
    cat("  STRING enrichment DB not found. Cannot run enrichment.\n")
  } else {
    cat("  Too few methylated TFs for enrichment analysis.\n")
  }
}

# ── Are DE TFs also methylated? ──
cat("\nDifferentially expressed TFs with methylation:\n")

if (exists("degs") && nrow(degs) > 0) {
  de_tfs <- degs[gene_id %in% tf_gene_ids & padj < 0.05]
  cat(sprintf("  DE TFs (padj < 0.05): %d\n", nrow(de_tfs)))

  if (nrow(de_tfs) > 0) {
    de_tfs[, has_dmp := gene_id %in% dmps$nearest_gene]
    de_tfs[, has_dmr := gene_id %in% dmrs$nearest_gene]

    de_tfs <- merge(de_tfs, annot_unique, by = "gene_id", all.x = TRUE)

    de_tfs_summary <- de_tfs[, .(gene_id, gene_name, log2FoldChange, padj, has_dmp, has_dmr)]
    de_tfs_summary <- de_tfs_summary[order(padj)]

    cat(sprintf("  DE TFs with DMPs: %d (%.1f%%)\n",
                sum(de_tfs$has_dmp), 100 * mean(de_tfs$has_dmp)))
    cat(sprintf("  DE TFs with DMRs: %d (%.1f%%)\n",
                sum(de_tfs$has_dmr), 100 * mean(de_tfs$has_dmr)))

    cat("\n  Top DE TFs:\n")
    top_de_tfs <- head(de_tfs_summary, 20)
    for (i in 1:nrow(top_de_tfs)) {
      meth_label <- paste0(ifelse(top_de_tfs$has_dmp[i], "DMP", ""),
                           ifelse(top_de_tfs$has_dmp[i] & top_de_tfs$has_dmr[i], "+", ""),
                           ifelse(top_de_tfs$has_dmr[i], "DMR", ""))
      if (meth_label == "") meth_label <- "none"

      cat(sprintf("    %-15s %-15s log2FC=%+.3f  padj=%.2e  meth=%s\n",
                  top_de_tfs$gene_id[i],
                  ifelse(is.na(top_de_tfs$gene_name[i]), "", top_de_tfs$gene_name[i]),
                  top_de_tfs$log2FoldChange[i],
                  top_de_tfs$padj[i],
                  meth_label))
    }

    fwrite(de_tfs_summary, file.path(OUTPUT_DIR, "tables", "de_tfs_methylation_status.tsv"),
           sep = "\t")
  }
}


# =============================================================================
# FINAL SUMMARY
# =============================================================================

cat("\n\n")
cat("================================================================\n")
cat("  ANALYSIS COMPLETE                                             \n")
cat("================================================================\n\n")

cat("PART 1 — Sox19a Gviz visualization:\n")
if (has_gviz) {
  cat("  Gviz locus plot generated with gene model + DMP tracks\n")
} else {
  cat("  Install Gviz for better visualization\n")
}

cat("\nPART 2 — Motif scanning:\n")
if (has_tfbs && !is.null(jaspar_pkg) && has_bsgenome) {
  cat(sprintf("  Scanned with %s motifs from %s\n", "metazoan", jaspar_pkg))
} else {
  cat("  Motif scanning requires TFBSTools + JASPAR + BSgenome.Dlaeve\n")
}

cat("\nPART 3 — TF ChIPseeker:\n")
cat(sprintf("  DMPs near TFs: %d | DMPs near Non-TFs: %d\n", nrow(dmps_tf), nrow(dmps_nontf)))
cat("  Pie charts generated for both categories\n")

cat("\nPART 4 — Methylated TF enrichment:\n")
cat(sprintf("  Methylated TFs: %d\n", length(tf_methylated)))

cat(sprintf("\nOutput: %s/\n", OUTPUT_DIR))

cat("\nKEY PACKAGES FOR ENHANCER/SILENCER ANALYSIS:\n")
cat("  Gviz           — Genomic track visualization\n")
cat("  TFBSTools      — TF binding site tools\n")
cat("  JASPAR2024     — Motif database (metazoan, vertebrate, etc.)\n")
cat("  motifmatchr    — Fast motif scanning (alternative to TFBSTools)\n")
cat("  BSgenome       — Reference genome sequences\n")
cat("  chromVAR       — Chromatin variability (needs ATAC-seq)\n")
cat("  GenomicRanges  — Core interval operations\n")
cat("  ChIPseeker     — Peak/DMP annotation + plots\n")
cat("================================================================\n")
