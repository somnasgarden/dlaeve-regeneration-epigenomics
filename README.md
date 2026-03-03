# Epigenomic Regulation of Tail Regeneration in *Deroceras laeve*

How does DNA methylation control gene expression during tissue regeneration in an invertebrate?

This project integrates whole-genome bisulfite sequencing (WGBS), RNA-seq, and gene regulatory network analysis to answer that question in the gray field slug *Deroceras laeve* — an emerging model for regenerative biology.

---

## Background

Tail amputation in *D. laeve* triggers a regenerative response involving widespread changes in both DNA methylation and gene expression. Unlike in mammals, where promoter methylation is the primary regulatory mechanism, our analysis reveals that **methylation changes in intergenic and intronic regions** — not promoters — drive expression changes during regeneration. This suggests regulation through enhancers, silencers, and distal regulatory elements, challenging the classical promoter-centric model of epigenetic regulation.

## Repository Structure

```
dlaeve-regeneration-epigenomics/
│
├── cluster/                          Scripts that run on the computing cluster (SLURM)
│   ├── 01-methylation/                   WGBS: QC → DMPs/DMRs → annotation → enrichment → TEs
│   ├── 02-rnaseq/                        RNA-seq: normalization → batch correction → DEGs → clustering
│   ├── 03-coexpression/                  WGCNA: module detection, hub genes, trait correlations
│   ├── 04-integration/                   Methylation–expression correlation + module enrichment
│   ├── 05-promoter/                      Promoter CpG island classification (HCP/ICP/LCP)
│   ├── 06-tf-prediction/                 DeepTFactor TF prediction + validation
│   ├── 07-expanded-analysis/             Permutation tests, bootstrap CIs, network null models
│   └── 08-wgcna-sensitivity/             WGCNA parameter sweep + methylation-network integration
│
├── local/                            Scripts for local PC analysis (figures, exploration)
│   ├── 01_blue_module_enrichment_plots.R
│   ├── 02_threshold_sensitivity_analysis.R
│   ├── 03_methylation_expression_predictors.R
│   ├── 04_candidate_genes_for_cloning.R
│   ├── 05_deep_methylation_analysis.R
│   ├── 06_novel_findings_analysis.R
│   ├── 07_sox19a_intergenic_regulation.R
│   ├── 08_morphogenesis_dmp_enrichment.R
│   ├── 09_dmr_spatial_analysis.R
│   ├── 10_methylation_threshold_and_synthesis.R
│   ├── 11_publication_figures.R
│   ├── 12_dmr_deep_analysis.R
│   ├── 13_entropy_and_novel_patterns.R
│   ├── 14_machine_learning.R
│   ├── 15_te_age_methylation.R
│   ├── 16_cross_species_comparison.R
│   ├── 17_te_tf_methylation_landscape.R
│   └── 18_sox19a_gviz_tf_chipseeker.R
│
├── results/                          Output data from cluster + local runs (not tracked)
│
├── figures/                          Publication-ready figures
│   ├── gene-level/                       fig01–fig10: methylation vs expression by region
│   └── module-level/                     fig11–fig19: module DMP enrichment, hub genes
│
├── papers/                           Paper drafts and references
│   ├── 00_draft_outline.md               Findings outline with 15 key results
│   └── 02_comprehensive.md              Full comprehensive paper (~8000 words)
│
├── report.html                      Interactive HTML5 wiki with all figures and analysis
├── CLAUDE.md                         Project context for AI-assisted development
├── .gitignore
└── LICENSE                           MIT
```

## Pipeline Overview

The analysis proceeds in eight stages, each with its own subfolder under `cluster/`:

| Stage | Folder | Input | Output |
|-------|--------|-------|--------|
| 1. Methylation | `01-methylation/` | Bismark CpG reports | 18,754 DMPs, 1,424 DMRs |
| 2. Transcriptomics | `02-rnaseq/` | HTSeq count matrices | Normalized counts, DEGs |
| 3. Co-expression | `03-coexpression/` | Batch-corrected VST | 12 WGCNA modules, hub genes |
| 4. Integration | `04-integration/` | DMPs + DEGs + modules | Meth–expr correlations, module enrichment |
| 5. Promoters | `05-promoter/` | Genome + methylation | HCP/ICP/LCP classification |
| 6. TF prediction | `06-tf-prediction/` | Protein sequences | 4,824 predicted TFs |
| 7. Expanded analysis | `07-expanded-analysis/` | All prior results | Permutation/bootstrap validation |
| 8. WGCNA sensitivity | `08-wgcna-sensitivity/` | RNA-seq + methylation | Parameter sweep, bicor network, topology metrics |

Each folder contains a `run.slurm` template: `sbatch run.slurm <script.R>`

Local scripts (`local/01–18`) perform downstream analysis, figure generation, information-theoretic exploration, TE/TF methylation landscape analysis, and Sox19a locus visualization.

## Experimental Design

| Data type | Samples | Comparison |
|-----------|---------|------------|
| **WGBS** (methylation) | C1, C2 vs A1, A2 | Control vs Amputated tails |
| **RNA-seq** (expression) | C1S1–C4S4 vs T2S6–T4S8 | Control vs Amputated tails |
| **RNA-seq** (multi-tissue) | 6 tissues × 4 conditions | Meta-analysis for co-expression |

Tissues: tail, eye, body wall, head, juvenile, ovotestis
Conditions: control, amputated, irradiated, fungicide

## Key Findings

**1. Methylation direction carries zero information about expression direction.**
Mutual information between methylation direction and expression direction is MI = 0.00009 bits (p = 0.55). The classical "hypermethylation = silencing" model is entirely absent.

**2. Intergenic methylation is the primary regulatory mechanism — not promoters.**
63% of DMPs fall in intergenic/intronic regions. DMPs cluster at 73x the expected rate within 1 kb, forming coherent regulatory domains.

**3. Transcription factors are methylated but expression-protected.**
DeepTFactor CNN-based TFs are enriched for DMPs (OR = 1.607, p = 4.32e-13) and DMRs (OR = 1.645, p = 3.33e-06). However, 0/6 DE TFs have any methylation changes. Top methylated TFs: ZFP2 (30 DMPs), Mafb (25), ATF2 (18). Hox genes: 0/27 with DMPs.

**4. Module-specific methylation erases global correlation.**
Global methylation-expression correlation is null (rho = -0.005, p = 0.80), but within-module correlations reach |r| = 0.81. Opposing module effects cancel at the genome-wide level.

**5. Hub genes are methylated but not differentially expressed.**
Blue module hub genes are enriched for DMRs (OR = 1.97, padj = 0.033) but none is DE. Network centrality positively correlates with methylation responsiveness (rho = 0.176, p = 0.007), contradicting buffering.

**6. Promoter landscape dominated by intermediate CpG content.**
91.6% of genes have Intermediate CpG Promoters (ICP), diverging sharply from mammalian genomes and rendering the CpG island model inapplicable.

**7. Methylation is a binary switch, not a rheostat.**
No dose-response relationship between methylation magnitude and DE probability (Cochran-Armitage z = 1.009, p = 0.31).

## Dependencies

**Bioconductor:** bsseq, DSS, GenomicRanges, rtracklayer, ChIPseeker, DESeq2, WGCNA, clusterProfiler

**CRAN:** data.table, dplyr, tidyr, ggplot2, pheatmap, RColorBrewer

**External:** Bismark (upstream alignment), GENIE3 (regulatory networks), DeepTFactor (TF prediction)

## Genome

All analyses use a custom *D. laeve* genome assembly (`BSgenome.Dlaeve.NCBI.dlgm`) with EviAnn-based annotation. Standard genome databases do not cover this organism.

## Data

Raw sequencing data and large processed objects (.rds, .RData) are not included due to size. See `CLAUDE.md` for expected data paths and file descriptions.

## License

MIT
