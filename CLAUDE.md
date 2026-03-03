# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Researcher

Rafa (rlopezt) — molecular genetics / bioinformatics. Studying epigenetic regulation of tail regeneration in the gastropod mollusk *Deroceras laeve*. Background in molecular genetics, stem cells, bioinformatics, and neuroscience.

## Project summary

This project asks: **how does DNA methylation regulate gene expression during tail regeneration in an invertebrate?**

We integrate whole-genome bisulfite sequencing (WGBS), RNA-seq across multiple tissues and conditions, and gene regulatory network inference to build a mechanistic picture of regenerative epigenomics in *D. laeve*.

## Organism & genome

- **Species**: *Deroceras laeve* (gray field slug, gastropod mollusk)
- **Genome assembly**: `BSgenome.Dlaeve.NCBI.dlgm` (custom)
- **Annotation**: `derLaeGenome_namesDlasi_v2.fasta.functional_note.pseudo_label.gff` (EviAnn)
- This is a non-model invertebrate. **Do NOT assume mammalian regulatory paradigms apply.**

## Experimental design

### Methylation (WGBS)
- Control: C1, C2 (intact tails)
- Amputated: A1, A2 (regenerating tails)
- Input: Bismark CpG_report.txt.gz files

### Transcriptomics (RNA-seq)
- Multi-tissue meta-analysis: tail, eye, bodywall, head, juvenile, ovotestis
- Conditions: control, amputated, irradiated, fungicide
- Tail amputation contrast: C1S1–C4S4 (Control) vs T2S6–T4S8 (Amputated)
- **T1S5 was removed as an outlier during WGCNA**
- Batch-corrected VST counts used for co-expression analysis

## Repository structure

```
cluster/          R scripts + SLURM jobs (8 pipeline stages, 41 files)
local/            R scripts for local analysis (18 scripts)
results/          Output data (untracked, 25 subdirectories)
figures/          Publication-ready figures (fig01–fig19)
papers/           Paper drafts + reference PDFs
report.html       HTML5 wiki-style project report
```

### `cluster/` — Pipeline scripts (run on SLURM)

Each subfolder has a `run.slurm` template: `sbatch run.slurm <script.R>`

#### `cluster/01-methylation/` — WGBS differential methylation (DSS)

| Script | What it does | Key outputs |
|--------|-------------|-------------|
| `01_quality_control.R` | Load Bismark files, build BSseq object, QC plots | `bsseq_object.rds` |
| `02_differential_methylation.R` | DSS statistical testing for DMPs and DMRs | `dmp_results.rds`, `dmr_results.rds` |
| `03_genomic_annotation.R` | Map DMPs/DMRs to promoters, gene body, intergenic | Annotated tables |
| `04_functional_enrichment.R` | ChIPseeker annotation + GO/KEGG/Reactome (STRING-based) | Enrichment tables |
| `05_te_methylation.R` | Vectorized TE–methylation overlap (~100x faster) | TE overlap tables |

#### `cluster/02-rnaseq/` — RNA-seq normalization, DE, clustering

| Script | What it does | Key outputs |
|--------|-------------|-------------|
| `01_normalization.R` | DESeq2 normalization, sample QC, PCA | Normalized count matrices |
| `02_batch_correction.R` | limma `removeBatchEffect` for experiment effects | Batch-corrected VST counts |
| `03_differential_expression.R` | DESeq2 DEGs per condition vs control | DEG tables, union list |
| `04_hierarchical_clustering.R` | Silhouette-optimized hierarchical clustering of DEGs | Cluster assignments |
| `05_kmeans_clustering.R` | K-means clustering + comparison with hierarchical | Cluster assignments, ARI |

#### `cluster/03-coexpression/` — WGCNA (64GB, 8 CPUs)

| Script | What it does | Key outputs |
|--------|-------------|-------------|
| `01_wgcna_module_detection.R` | Signed network, module detection, hub genes, trait correlations, enrichment, preservation | Module assignments, eigengenes, hub lists, enrichment tables |

#### `cluster/04-integration/` — Methylation-expression integration

| Script | What it does | Key outputs |
|--------|-------------|-------------|
| `01_methylation_vs_expression.R` | Gene-level DMP/DMR methylation vs log2FC correlation by region | Scatter plots, correlation tables, significant gene lists |
| `02_module_methylation_enrichment.R` | Fisher's exact test for module DMP enrichment, hub gene methylation status | Module enrichment tables, master summary, heatmaps |

#### `cluster/05-promoter/` — Promoter CpG classification

| Script | What it does | Key outputs |
|--------|-------------|-------------|
| `01_promoter_classification.R` | Weber criteria HCP/ICP/LCP classification, methylation-expression by class | Classification tables, HTML report |

#### `cluster/06-tf-prediction/` — Transcription factor prediction

| Script | What it does | Key outputs |
|--------|-------------|-------------|
| `01_tf_validation.R` | Validate DeepTFactor predictions vs functional annotations | Contingency table, Fisher's test |
| `deeptfactor/` | DeepTFactor Python tool (CNN model, environment, runner) | `prediction_result.txt` (73,425 proteins) |

#### `cluster/07-expanded-analysis/` — Statistical validation (64GB)

| Script | What it does | Key outputs |
|--------|-------------|-------------|
| `01_permutation_enrichment.R` | Permutation tests (10K) for module DMP/DMR enrichment | Empirical p-values |
| `02_network_null_models.R` | Degree-preserving null models for network topology | Null distribution tables |
| `03_bootstrap_correlations.R` | Bootstrap CIs for module-region correlations | Confidence interval tables |

#### `cluster/08-wgcna-sensitivity/` — WGCNA parameter sweep (64GB, 8 CPUs)

| Script | What it does | Key outputs |
|--------|-------------|-------------|
| `01_parameter_sensitivity.R` | 5 variance filters x full power sweep, outlier detection | Best parameters, R² plots |
| `02_improved_wgcna.R` | Bicor network with optimized parameters, module comparison | Improved modules, TOM, hub genes |
| `03_methylation_network_integration.R` | Topology-methylation metrics (WMCS, hub gradient, TOM concordance) | Network-level tables, master gene table |
| `04_enrichment_across_parameters.R` | 7 filters × 9 powers with enrichment tracking, sulfation emergence, connectivity analysis | Enrichment shift tables, sulfation tracking, ARI heatmap |

### `local/` — Local PC scripts (18 analysis scripts)

| Script | What it does | Output dir |
|--------|-------------|-----------|
| `01_blue_module_enrichment_plots.R` | Blue module GO plots | results/02_rnaseq/ |
| `02_threshold_sensitivity_analysis.R` | DMP threshold sensitivity | results/06_threshold_sensitivity/ |
| `03_methylation_expression_predictors.R` | Module x region correlations | results/10_meth_expr_predictors/ |
| `04_candidate_genes_for_cloning.R` | Candidate gene scoring | results/11_candidate_genes/ |
| `05_deep_methylation_analysis.R` | 10-part deep methylation dive | results/07_deep_analysis/ |
| `06_novel_findings_analysis.R` | 7 novel findings | results/12_novel_findings/ |
| `07_sox19a_intergenic_regulation.R` | Sox19a case study | results/13_intergenic_regulation/ |
| `08_morphogenesis_dmp_enrichment.R` | Category enrichment | results/14_morphogenesis_enrichment/ |
| `09_dmr_spatial_analysis.R` | DMR spatial analysis | results/08_dmr_spatial/ |
| `10_methylation_threshold_and_synthesis.R` | Threshold + synthesis | results/15_synthesis/ |
| `11_publication_figures.R` | Publication figures | results/22_publication_figures/ |
| `12_dmr_deep_analysis.R` | DMR deep analysis, TF tests | results/09_dmr_deep_analysis/ |
| `13_entropy_and_novel_patterns.R` | Information theory, entropy | results/16_entropy_novel/ |
| `14_machine_learning.R` | ML prediction of DE from methylation features (RF, GBM, nnet) | results/17_machine_learning/ |
| `15_te_age_methylation.R` | TE methylation by evolutionary age (Kimura → My) | results/18_te_age_methylation/ |
| `16_cross_species_comparison.R` | Cross-species methylation comparison for regeneration | results/19_cross_species/ |
| `17_te_tf_methylation_landscape.R` | TE genomic context, genomation meta-gene profile, TF methylation | results/20_te_tf_landscape/ |
| `18_sox19a_gviz_tf_chipseeker.R` | Sox19a Gviz locus, motif scanning (JASPAR), TF ChIPseeker pies, methylated TF enrichment | results/21_sox19a_tf_motifs/ |

### `papers/` — Paper drafts and references

| File | What it is |
|------|-----------|
| `00_draft_outline.md` | Findings outline (15 key results, figure plan, candidate genes) |
| `02_comprehensive.md` | Full paper (~8000 words, 10 figures, complete methods) |
| `Evidence that direct inhibition of .pdf` | Reference PDF (~23 MB) |
| `The chromosome-level genome assembly of the slug...pdf` | *D. laeve* genome assembly paper (~14 MB) |

### `results/` — Analysis output data (untracked, 25 subdirectories)

| Subfolder | Source | Contents |
|-----------|--------|----------|
| `01_methylation/` | cluster/01 | BSseq object, DMP/DMR tables, ChIPseeker, enrichment, TE analysis |
| `02_rnaseq/` | cluster/02-03 | Part1–Part5_WGCNA: counts, DEGs, clusters, modules |
| `03_integration/` | cluster/04 | MXT figures, tables, motif analysis BED files |
| `04_promoter/` | cluster/05 | HCP/ICP/LCP classification, methylation stats |
| `05_tf_prediction/` | cluster/06 | DeepTFactor prediction_result.txt, protein FASTA |
| `06_threshold_sensitivity/` | local/02 | DMP threshold sensitivity tables |
| `07_deep_analysis/` | local/05 | 15 deep analysis tables |
| `08_dmr_spatial/` | local/09 | DMR spatial clustering |
| `09_dmr_deep_analysis/` | local/12 | DMR enrichment, TF tests, large DMRs |
| `10_meth_expr_predictors/` | local/03 | Module x region Spearman correlations |
| `11_candidate_genes/` | local/04 | Candidate gene scoring tables |
| `12_novel_findings/` | local/06 | Novel findings tables |
| `13_intergenic_regulation/` | local/07 | Sox19a intergenic regulation |
| `14_morphogenesis_enrichment/` | local/08 | Category DMP enrichment |
| `15_synthesis/` | local/10 | Threshold + key statistics |
| `16_entropy_novel/` | local/13 | Shannon entropy, binary switch, four-layer model |
| `17_machine_learning/` | local/14 | ML prediction of DE from methylation |
| `18_te_age_methylation/` | local/15 | TE methylation by evolutionary age |
| `19_cross_species/` | local/16 | Cross-species methylation comparison |
| `20_te_tf_landscape/` | local/17 | TE genomic context, TF DMP/DMR enrichment (combined) |
| `20_te_methylation/` | local/17 (split) | TE-DMP genomic context plots + tables |
| `21_sox19a_tf_motifs/` | local/18 | Sox19a locus, TF ChIPseeker pies (combined) |
| `21_tf_methylation/` | local/17 (split) | TF methylation rates, module enrichment, top TFs |
| `22_publication_figures/` | local/11 | Publication-ready figure set |
| `22_sox19a_locus/` | local/18 (split) | Sox19a locus plots + tables |

### `figures/` — Publication-ready figures (19 files)

- `gene-level/` — fig01–fig10: methylation vs expression by region, quadrant density, significant genes
- `module-level/` — fig11–fig19: module DMP enrichment, Fisher tests, hub gene status, heatmaps

### `report.html` — HTML5 multi-page wiki report

Interactive 9-page report with JavaScript page switching, sidebar navigation with sub-section links, keyboard arrow navigation, and lightbox figure viewer. Pages: Home, Key Findings, Methylation, RNA-seq, Integration, Advanced, Context, Figures, Technical.

## Data paths

### Local reference data (read-only — NEVER write or modify)

`C:\Users\rafae\Projects\DATA\`

| File / folder | Description |
|---------------|-------------|
| `A1.CpG_report.txt`, `A2.CpG_report.txt` | Bismark CpG reports — amputated samples |
| `C1.CpG_report.txt`, `C2.CpG_report.txt` | Bismark CpG reports — control samples |
| `derLaeGenome_namesDlasi_v2.fasta.functional_note.pseudo_label.gff` | *D. laeve* GFF annotation (EviAnn) |
| `derLaeGenome_eviann_annotations.tsv` | EviAnn annotations table |
| `collapsed_te_age_data.tsv` | Collapsed transposable element age data |
| `genie3_top500k.tsv` | GENIE3 top 500k regulatory predictions |
| `STRG0A31YWK.protein.orthology.v12.0.txt` | Protein orthology file |
| `prediction_result.txt` | DeepTFactor TF predictions (73,425 proteins, 4,824 TFs) |
| `counts_HTseq_EviAnn/` | HTSeq count matrices (EviAnn annotation) |

### Cluster paths (remote)

```
BASE_DIR  = /mnt/data/alfredvar/
METH_DIR  = /mnt/data/alfredvar/rlopezt/Preliminary/
TRANS_DIR = /mnt/data/alfredvar/rlopezt/CorrelationMatrix/
MXT_DIR   = /mnt/data/alfredvar/rlopezt/MXT/
GFF_FILE  = /mnt/data/alfredvar/30-Genoma/31-Alternative_Annotation_EviAnn/
            derLaeGenome_namesDlasi_v2.fasta.functional_note.pseudo_label.gff
```

Scripts use a **waterfall loading pattern**: try local `DATA` path first → cluster path → rebuild from upstream. See "Architecture: dual-path system" for details.

## Key biological findings

### 1. MI = 0 between methylation direction and expression direction
Mutual information is 0.00009 bits (p = 0.55). The classical "hypermethylation = silencing" rule is entirely absent genome-wide.

### 2. Intergenic methylation is the main regulatory mechanism
63% of DMPs in intergenic/intronic regions. DMPs cluster at 73x expected rate within 1 kb. Regulation through enhancers and silencers, NOT promoters.

### 3. TF methylation paradox (revised four-layer model)
DeepTFactor CNN TFs are **ENRICHED** for DMPs (OR = 1.607, p = 4.32e-13) and DMRs (OR = 1.645, p = 3.33e-06). 14.8% of TFs have DMPs vs 9.8% of non-TFs. However, 0/6 DE TFs have any methylation changes — methylation targets TFs but does NOT change their expression. Old keyword-based finding (OR = 0.77 depleted) was an artifact of incomplete TF identification.

### 4. Module-specific methylation erases global correlation
Global rho = -0.005 (NS), but within-module |r| up to 0.81. Yellow module sole DMP-enriched (68% thresholds). Yellow sole hyper-dominant; 9/10 others hypo-dominant.

### 5. Hub gene methylation paradox
Blue hub genes enriched for DMRs (OR = 1.97, padj = 0.033) but NONE are DE. kME positively correlates with methylation response (rho = 0.176, p = 0.007) — contradicts buffering.

### 6. Promoter landscape dominated by ICPs
91.6% ICP, 6.5% HCP, 1.9% LCP. Unprecedented in any characterized organism. Renders CpG island model inapplicable.

### 7. Binary switch, not rheostat
Cochran-Armitage z = 1.009, p = 0.31 (no dose-response). Methylation is qualitative (which genes/modules), not quantitative (how much).

## Current analysis frontier

The core analysis is complete. All 18 local scripts have been run. Three paper drafts exist in `papers/`. A multi-page HTML5 wiki report (`report.html`) summarizes the full project.

### New local analyses (scripts 17–18, completed Mar 2026)
- **TE/TF methylation landscape** (local/17): 51.1% DMPs overlap TEs (enriched intergenic OR=1.305, p=3.10e-17). DeepTFactor TFs ENRICHED for DMPs (OR=1.607, p=4.32e-13). Top TFs: ZFP2 (30 DMPs), Mafb (25), ATF2 (18). 0/6 DE TFs are methylated.
- **Sox19a Gviz + TF ChIPseeker** (local/18): 13 hypermethylated DMPs ~4-5kb upstream of Sox19a. Prediction: silencer inactivation. TF promoter DMPs 25% vs non-TF 20.4%. 6 DE TFs identified (Scrt1, REPTOR-BP, Y-box, Catenin, sox8, CEBPA), none methylated.

### Completed local analyses (scripts 14–16, added Feb 2025)
- **Machine learning** (local/14): RF, GBM, nnet prediction of DE from methylation features
- **TE age methylation** (local/15): TE methylation as function of evolutionary age (Kimura → My)
- **Cross-species comparison** (local/16): D. laeve vs other species methylation patterns

### Cluster scripts to run
- `cluster/07-expanded-analysis/` — Permutation, bootstrap, network null models (validation)
- `cluster/08-wgcna-sensitivity/` — WGCNA parameter sweep (4 scripts: parameter sensitivity, improved WGCNA with bicor, methylation-network integration, enrichment across parameters with sulfation tracking)

### Future experiments
- Motif analysis (BED files prepared in `results/03_integration/Motif_Analysis/`)
- ATAC-seq to confirm enhancer accessibility at DMP clusters
- Time-series WGBS for temporal dynamics
- 5-azacytidine treatment to test methylation requirement
- Sox19a upstream regulatory region reporter assay

## Running the pipeline

### On the cluster (SLURM)
Each pipeline stage has a `run.slurm` template. Submit jobs with:
```bash
cd cluster/01-methylation/
sbatch run.slurm 01_quality_control.R
```
Memory requirements: 32GB for most stages, 64GB + 8 CPUs for WGCNA (`03-coexpression/`).

### Locally (R scripts)
```bash
Rscript local/01_blue_module_enrichment_plots.R
```
Local scripts expect data in `C:\Users\rafae\Projects\DATA\` (WSL: `/mnt/c/Users/rafae/Projects/DATA/`).

### DeepTFactor (Python)
```bash
cd cluster/06-tf-prediction/deeptfactor/
conda env create -f environment.yml   # Python 3.6, PyTorch 1.2, CUDA 10.0
conda activate deeptfactor
python tf_running.py --input <fasta> --output <outdir>
```

## Pipeline dependency graph

Scripts within each pipeline run sequentially — each loads outputs from the previous step. Cross-pipeline dependencies exist for integration.

```
01-methylation:  01_qc → 02_dm → 03_annot → 04_enrich → 05_te
                   ↓ bsseq_object.rds   ↓ dmp/dmr_results.rds
                                         ↓
02-rnaseq:  01_norm → [02_batch] → 03_de → 04_hier → 05_kmeans
              ↓ Part1_objects.RData  ↓ Part1b  ↓ Part2  ↓ Part3  ↓ Part4
                                                                   ↓
03-coexpression:  01_wgcna ← Part4_objects.RData
                    ↓ Part5_WGCNA_objects.RData
                    ↓
04-integration:  01_meth_vs_expr ← methylation + rnaseq + wgcna
                 02_module_enrich ← methylation + wgcna

05-promoter:     01_promoter_class  (semi-independent, reads raw cluster data)
06-tf-prediction: 01_tf_validation  (standalone, uses pre-computed DeepTFactor)
```

**Batch correction is optional**: `03_differential_expression.R` has a commented toggle to load from either `Part1_QC/` (no batch correction) or `Part1b_BatchCorrection/`.

## Architecture: dual-path system

### Methylation pipeline — absolute cluster paths
Scripts hardcode `RESULTS_DIR <- "/mnt/data/alfredvar/rlopezt/feb_w2_scripts/Preliminary/"` and save `.rds` files there. Each script redefines `BASE_DIR` and `RESULTS_DIR` at the top.

### RNA-seq pipeline — relative Part-numbered directories
Scripts create sequentially numbered output folders in the working directory:
- `Part1_QC/{plots,data}/` → `Part1_objects.RData`
- `Part1b_BatchCorrection/{plots,data}/` → `Part1b_objects.RData`
- `Part2_DEGs/{plots,data}/` → `Part2_objects.RData`
- `Part3_Hierarchical/{plots,data}/` → `Part3_objects.RData`
- `Part4_Kmeans/{plots,data}/` → `Part4_objects.RData`
- `Part5_WGCNA/{inputs,plots,data}/` → `Part5_WGCNA_objects.RData`

### Integration pipeline — MXT directory on cluster
Uses `MXT_DIR` with subdirectories: `Tables/`, `PDFs/`, `PNGs/`, `RDS_Objects/`.

### Local scripts — relative to repo root
Load from `results/02_rnaseq/Part5_WGCNA/data/`, output to `results/02_rnaseq/Part5_WGCNA/plots/`.

## Technical conventions

### Self-installing packages
Scripts auto-install missing dependencies at startup via BiocManager. The standard header pattern:
```r
essential_packages <- c("bsseq", "DSS", ...)
for (pkg in essential_packages) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
}
suppressPackageStartupMessages({ library(bsseq); ... })
```

### Global options (set in every script)
```r
options(stringsAsFactors = FALSE)
options(scipen = 999)
```

### Plotting
- Dual save via `save_both()` and `save_pheatmap_both()` helpers: every figure saved as PDF + PNG (300 DPI)
- Consistent ggplot2 theme: `theme_minimal()`, 12pt text, 14pt bold titles
- Figure naming: `fig##_descriptive_name.png` or `MXT_##_descriptive_name.png`

### Analysis parameters object
Scripts define a standard `analysis_params` list at the top:
```r
analysis_params <- list(
  min_coverage = 5, min_samples = 2, min_diff = 0.1,
  fdr_threshold = 0.05, min_dmr_cpgs = 3, min_dmr_length = 50
)
```

### WGCNA threading
WGCNA scripts call `enableWGCNAThreads(nThreads = 8)` — must match `#SBATCH --cpus-per-task=8`.

### Key R packages
- `bsseq`, `DSS` — BS-seq analysis
- `GenomicRanges`, `rtracklayer` — Genomic intervals
- `ChIPseeker` — Region annotation
- `DESeq2` — Differential expression
- `WGCNA` — Co-expression networks
- `GENIE3` — Regulatory network inference
- `Gviz` — Genomic track visualization
- `TFBSTools`, `JASPAR2024` — Motif scanning
- `genomation` — Meta-gene methylation profiles
- `ggplot2`, `pheatmap` — Visualization

## CRITICAL: What NOT to assume

1. **Do NOT use human/mouse databases** for annotation — this is a non-model invertebrate
2. **Do NOT assume promoter methylation matters most** — intergenic regions are where the action is
3. **Do NOT discard intergenic DMPs as noise** — they are likely the most functionally relevant
4. **Do NOT assume methylation → silencing** — the relationship is complex and context-dependent
5. **The genome is custom** — standard databases (UCSC, Ensembl) don't cover this organism
