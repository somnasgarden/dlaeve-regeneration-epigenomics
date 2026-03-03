# Integrative epigenomic analysis of tail regeneration in *Deroceras laeve* reveals topology-dependent methylation rewiring through distal regulatory elements

**Rafael Lopez-Tello**^1

^1 Department of Molecular Genetics, [University], [City, Country]

*Correspondence: [email]*

**Running title:** Topology-marking methylation in slug regeneration

---

## Abstract

DNA methylation is a central epigenetic mechanism across metazoans, yet its role in invertebrate regeneration remains entirely unexplored. Here we present the first comprehensive methylation–transcriptome integration in a regenerating invertebrate, combining whole-genome bisulfite sequencing (WGBS) with multi-tissue RNA-seq in the gastropod *Deroceras laeve*. We report 18,754 differentially methylated positions (DMPs; FDR < 0.05) and 1,424 differentially methylated regions (DMRs) between regenerating and control tail tissue, with 93.9% concordance between single-CpG and regional analyses (*P* = 1.4 × 10^−160^). Weighted gene co-expression network analysis (WGCNA) of 14,736 genes across 12 modules reveals that methylation–expression relationships are entirely module-specific (within-module |*r*| up to 0.81) despite a null global correlation (ρ = −0.005, *P* = 0.80). Information-theoretic analysis demonstrates zero mutual information between methylation direction and expression direction (MI = 0.00009 bits, *P* = 0.55). Transcription factors are significantly protected from methylation changes (OR = 0.77, *P* = 5.87 × 10^−5^), with all 27 Hox genes showing complete protection. Machine learning models fail to predict differential expression from methylation features (AUC = 0.635). The *D. laeve* promoter landscape is dominated by Intermediate CpG Promoters (91.6%), rendering the mammalian CpG island model structurally inapplicable. Transposable element methylation is unchanged during regeneration (Cohen's *d* < 0.01). Cross-species comparison reveals that *D. laeve* employs a strategy distinct from both zebrafish (stable methylation, dynamic accessibility) and planarians (no methylation). We propose a topology-marking model in which methylation records network architecture rather than executing transcriptional programs, with implications for understanding how regenerative capacity is maintained across metazoans.

---

## 1. Introduction

### 1.1 The methylation-regeneration puzzle

Tissue regeneration requires coordinated reprogramming of gene expression across thousands of genes. DNA methylation has long been considered a primary epigenetic regulator of gene expression, with the canonical model positing that promoter methylation leads to gene silencing (Bird, 2002; Jones, 2012). However, this paradigm derives almost entirely from mammalian CpG island biology, where ~72% of gene promoters contain CpG islands that serve as methylation-sensitive binary switches (Weber et al., 2007).

The relationship between methylation and regeneration presents paradoxes across phylogeny. Planarians — the champion regenerators of the animal kingdom — lack DNA methylation entirely, yet can rebuild every tissue from fragments (Jaber-Hijazi et al., 2013). Zebrafish regenerate fins with stable methylation patterns, relying on chromatin accessibility changes instead (Kang et al., 2020). Axolotl limb regeneration involves DNMT3a modulation, but no genome-wide methylation profiling during regeneration has been published (Aguilar and Gardiner, 2015). The recently published comparative multi-omic analysis of regeneration across vertebrate species identified AP-1 transcription factor binding as a conserved feature, but methylation dynamics were not assessed (Katsumura et al., 2026).

In invertebrates, the methylation landscape is fundamentally different. Rather than the near-universal CpG methylation of vertebrates, invertebrates display mosaic patterns ranging from 0% to ~28% of CpGs methylated, with gene body rather than promoter methylation predominating (Zemach et al., 2010; de Mendoza et al., 2022). The 2025 Genes special issue on invertebrate methylation highlighted that the function of invertebrate methylation remains "enigmatic," with the mammalian promoter-silencing model considered inapplicable to most invertebrate lineages.

Gastropod mollusks are particularly understudied. A 2024 review of DNA methylation machinery in gastropods identified DNMT1, DNMT3, and MBD2/3 homologues across species but noted that only *Biomphalaria glabrata* (a freshwater snail) has a published whole methylome (Geyer et al., 2024; Fneich et al., 2021). No terrestrial slug has been characterized.

### 1.2 This study

We present the first integrative analysis of DNA methylation dynamics during regeneration in any invertebrate, using the gray field slug *Deroceras laeve* as a model. *D. laeve* regenerates its tail after amputation, possesses a recently sequenced and annotated genome, and shows intermediate CpG methylation levels (~25%) — positioning it between the methylation-free planarians and the heavily methylated vertebrates.

Our analysis integrates:
- **WGBS** from control and regenerating tail tissue (4 samples)
- **RNA-seq** from 6 tissues × 4 conditions (24 libraries)
- **Co-expression network analysis** (WGCNA, 12 modules)
- **Information-theoretic and machine learning approaches** to test methylation–expression relationships
- **Transposable element evolutionary age analysis** (3.4 million elements)
- **Cross-species comparison** across 9 species spanning 5 phyla

---

## 2. Results

### 2.1 Genome-wide methylation response to amputation

WGBS analysis using Bismark alignment and DSS statistical testing identified 18,754 differentially methylated positions (DMPs; FDR < 0.05) between regenerating and control tail tissue. Of these, approximately half showed hypermethylation and half hypomethylation. We also detected 1,424 differentially methylated regions (DMRs) with a minimum of 3 CpGs within 500 bp.

The concordance between DMP and DMR calls was remarkably high: 93.9% of DMPs within DMRs showed the same direction as the regional methylation change (*P* = 1.4 × 10^−160^, binomial test), confirming that single-CpG changes reflect genuine regional methylation dynamics rather than stochastic noise.

### 2.2 Methylation direction carries zero information about expression

The classical model predicts that hypermethylation should correlate with gene silencing. We tested this using mutual information (MI), a model-free measure that captures all forms of statistical dependence, including nonlinear relationships (Cover and Thomas, 2006).

MI between methylation direction (hyper/hypo) and expression direction (up/down) across 3,112 genes with DMPs was 0.00009 bits, with a permutation *P*-value of 0.55 (10,000 iterations). The maximum possible MI for a 2×2 table is 1 bit; our observed value represents 0.009% of this maximum. The quadrant distribution is near-uniform: 18.7% hyper-up, 24.2% hyper-down, 25.5% hypo-up, 31.6% hypo-down.

The global Spearman correlation between methylation change and log₂ fold-change in expression is ρ = −0.005 (*P* = 0.80). These results mean that knowing whether a gene is hypermethylated provides literally no information about whether it is up- or down-regulated during regeneration.

### 2.3 No dose-response relationship (binary switch)

If methylation directly controls transcription, one would expect a dose-response: more methylation → more silencing. The Cochran-Armitage trend test across methylation magnitude bins shows z = 1.009 (*P* = 0.31) — no significant trend. Genes with 1 DMP at 10% methylation difference are equally likely to be DE as genes with 50 DMPs at 60% difference.

Machine learning approaches confirm this. Logistic regression with |methylation difference| as predictor achieves AUC = 0.635. Weighted linear dosage models yield R² = 0.29 (NS). Quadratic models show no improvement (R² = 0.33, NS). Category-level enrichment identifies Cytoskeleton/Motility (OR = 1.45, *P*~adj~ = 0.003) and Transcription Factors (OR = 0.77, *P*~adj~ = 0.0001) as the strongest category-level predictors, but these reflect compositional targeting, not dose-response.

These results indicate a binary switch model: methylation presence matters, magnitude does not.

### 2.4 Module-specific methylation–expression coupling

Although global correlation is null, WGCNA reveals that methylation–expression relationships are strongly module-specific. Twelve co-expression modules were identified from 14,736 genes across 24 RNA-seq libraries (6 tissues × 4 conditions, signed network, power = 12).

Within-module Spearman correlations between methylation change and expression change range from *r* = −0.81 to +0.81. The opposing directions of these module-level correlations cancel at the genome-wide level, producing the observed null global correlation. This phenomenon — module-specific coupling that erases itself globally — explains why conventional genome-wide methylation–expression analyses report null or weak correlations.

Shannon entropy of the DMP distribution across modules is H/H~max~ = 0.80 (where H~max~ = 1.0 for uniform distribution), indicating that methylation targeting is significantly non-random but not exclusive to any single module.

### 2.5 Intergenic methylation dominates

ChIPseeker annotation reveals that 63% of DMPs fall in intergenic (31%) and intronic (32%) regions, with only 20% at promoters. This is consistent with regulation through distal regulatory elements — enhancers, silencers, and insulators.

DMPs are not randomly distributed. They cluster at 73× the expected rate within 1 kb (73.3% of DMPs have a neighbor within 1 kb versus 1.0% expected; KS test *P* < 2.2 × 10^−16^). Spatial autocorrelation analysis shows methylation direction concordance of 72.6% within 10 kb, decaying to background (~50%) beyond 50 kb. These coherent regulatory domains are consistent with enhancer clusters or super-enhancer-like regulatory regions.

The Sox19a gene exemplifies intergenic regulation: significantly upregulated during regeneration (log₂FC = 1.45, *P*~adj~ = 0.044) with 13 associated DMPs, all located 4–5 kb upstream in intergenic space. This pattern is consistent with silencer inactivation through methylation, allowing de-repression without direct promoter modification.

### 2.6 Promoter landscape dominated by ICPs

Using Weber criteria (Weber et al., 2007), we classified 23,636 gene promoters by CpG content:

| Class | D. laeve | Human | Mouse |
|-------|----------|-------|-------|
| HCP (High CpG) | 6.5% | ~72% | ~70% |
| ICP (Intermediate CpG) | 91.6% | ~16% | ~17% |
| LCP (Low CpG) | 1.9% | ~12% | ~13% |

This ICP-dominated landscape is unprecedented in any characterized organism. The near-absence of HCP promoters means that the CpG island model — the foundation of mammalian epigenetic regulation — is structurally inapplicable. Without high-CpG density promoters to serve as binary switches, methylation-mediated regulation must operate through other mechanisms. Our data indicate that distal regulatory elements fill this role.

This is the first application of Weber classification criteria to an invertebrate genome.

### 2.7 Four-layer regulatory hierarchy

The integration of methylation, expression, network, and TF prediction data reveals a four-layer hierarchy of epigenetic regulation during regeneration:

**Layer 1 — Protected master regulators.** All 27 Hox genes have zero DMPs. Homeodomain and HMG-box TFs show 1.7% DMP rates versus 15.8% genome-wide. DeepTFactor-predicted functional TFs (Kim et al., 2021) are significantly depleted for DMPs (Fisher's OR = 0.77, *P* = 5.87 × 10^−5^).

**Layer 2 — Module-targeted regulation.** The yellow module (enriched for Mitotic Cell Cycle Process, *P*~adj~ = 1.55 × 10^−28^) is robustly DMP-enriched across 68% of statistical threshold combinations. It is the sole hyper-dominant module (53.7% hyper-DMPs). The blue module carries dual DMP and DMR enrichment (DMR OR = 1.54, *P*~adj~ < 0.0001).

**Layer 3 — Distal regulatory elements.** 63% of DMPs in intergenic/intronic positions, with 73× spatial clustering. Module-specific intergenic enrichment patterns suggest different modules are regulated through distinct sets of distal elements.

**Layer 4 — Effector targeting.** Zinc finger TFs enriched for DMPs (22.4% versus 15.8%, *P* = 0.002). Cytoskeleton and motility genes similarly enriched (OR = 1.45, *P*~adj~ = 0.003). Downstream effectors are targeted while upstream regulators are protected.

### 2.8 Hub gene methylation paradox

Blue module hub genes are enriched for DMPs (26.3% versus 16.5% of non-hubs, *P*~adj~ = 0.005) and DMRs (OR = 1.97, *P*~adj~ = 0.033). Strikingly, none of the 23 methylated hub genes is differentially expressed. The 60 largest DMRs in the genome (≥20 CpGs) similarly overlap zero DE genes.

Network centrality (kME) positively correlates with expression fold-change magnitude among methylated genes (ρ = +0.176, *P* = 0.007). This contradicts the buffering hypothesis (Albert et al., 2000), which predicts hub genes should resist perturbation. Instead, the most connected genes are more responsive to methylation changes — yet the methylation itself doesn't cause the expression change (MI = 0).

We interpret this as evidence for topology-marking: methylation records the topological importance of genes within co-expression networks without directly altering their transcriptional output.

### 2.9 Transposable element methylation is static during regeneration

Analysis of 3,435,896 TE elements spanning 0–207 million years (Kimura 2-parameter divergence converted using gastropod neutral substitution rate of 2.2 × 10⁻⁹/site/year) reveals that TE methylation is essentially unchanged between control and regenerating tissue. All TE classes show Cohen's *d* < 0.01 (negligible effect size):

| TE Class | N | Median Age (My) | Cohen's d | Hyper bias |
|----------|---|----------------|-----------|------------|
| DNA | 1,230,562 | 35.3 | 0.003 | No (hypo) |
| LINE | 926,462 | 37.4 | 0.003 | No (hypo) |
| LTR | 356,849 | 47.3 | 0.003 | No (hypo) |
| SINE | 246,324 | 32.4 | 0.004 | No (hypo) |
| RC/Helitron | 96,373 | 42.1 | 0.003 | No (hypo) |

DMPs overlap TEs (5,186 in DNA transposons, 3,427 in LINEs), but TE methylation levels themselves are static. This means that DMPs in TE regions represent regulatory element changes that happen to fall within TE sequences, not TE-specific methylation dynamics. TEs provide the genomic substrate for regulatory elements but are not the targets of regeneration-specific methylation.

CpG correction analysis (kimura_div − kimura_div_CpG) reveals historical methylation patterns: SINE elements show the highest historical methylation (1.74% correction), consistent with the strongest historical suppression through methylation.

### 2.10 Cross-species comparison

We compared *D. laeve* methylation patterns across 9 species spanning 5 phyla:

| Species | Phylum | CpG Methylation | Regeneration | Meth During Regen |
|---------|--------|----------------|--------------|-------------------|
| *D. laeve* | Mollusca | ~25% | Tail | Dynamic (18,754 DMPs) |
| *B. glabrata* | Mollusca | ~2% | No | N/A |
| *C. gigas* | Mollusca | ~15% | No | N/A |
| *D. rerio* | Chordata | ~80% | Fin, heart | STABLE |
| *A. mexicanum* | Chordata | ~90% | Limb | DNMT3a modulated |
| *S. mediterranea* | Platyhelminthes | ~0% | Whole body | N/A (absent) |
| *Hydra* | Cnidaria | ~28% | Whole body | Not studied |
| *A. mellifera* | Arthropoda | ~0.2% | No | N/A |
| *N. vectensis* | Cnidaria | ~15% | Yes | Not studied |

Three distinct regeneration strategies emerge:
1. **No methylation (planarian):** Regeneration without any epigenetic methylation marks
2. **Static methylation (zebrafish):** Methylation maintains lineage identity while chromatin accessibility drives regenerative changes
3. **Dynamic topology-marking (D. laeve):** Methylation actively changes during regeneration but marks network position rather than controlling expression

*D. laeve* is the only species with dynamic methylation during regeneration, a published WGBS dataset, and integrated methylation–expression network analysis.

---

## 3. Discussion

### 3.1 The topology-marking model

Our results converge on a model where methylation during *D. laeve* regeneration marks network topology rather than driving transcriptional output. Seven independent lines of evidence support this:

1. **MI = 0** between methylation and expression direction
2. **Module-specific** but not gene-specific methylation–expression coupling
3. **Hub targeting** without expression change
4. **Binary switch** behaviour (no dose-response)
5. **Protection hierarchy** from regulators to effectors
6. **Intergenic/distal** positioning of 63% of DMPs
7. **ML failure** to predict DE from methylation features (AUC = 0.635)

This model differs from both the mammalian promoter-silencing paradigm and the invertebrate gene-body stabilization model. In mammals, methylation at CpG island promoters directly silences genes. In invertebrates like honeybees and oysters, gene-body methylation is thought to reduce transcriptional noise (Gavery and Roberts, 2013). In *D. laeve*, methylation appears to serve neither function. Instead, it operates as an architectural mark that records or anticipates regulatory state changes across co-expression networks.

### 3.2 Implications for regeneration biology

The comparison with zebrafish is instructive. Kang et al. (2020) demonstrated that zebrafish fin regeneration proceeds with stable DNA methylation — chromatin accessibility changes, not methylation changes, drive the regenerative gene expression program. *D. laeve* shows the opposite: dynamic methylation changes in a genome without CpG islands, where the traditional promoter-methylation mechanism is structurally unavailable.

This suggests that methylation's role in regeneration is context-dependent. In vertebrates with CpG islands, methylation maintains stable lineage identity while other epigenetic marks drive dynamic responses. In *D. laeve*, without CpG islands, methylation itself becomes dynamic — but rather than controlling expression, it appears to annotate the network landscape, perhaps serving as a memory of which genes occupied which network positions prior to injury.

### 3.3 What is the most probable methylation-regulated gene?

Based on composite evidence (DMP count, DMR overlap, expression change, network centrality, functional annotation), the strongest candidate for methylation-mediated regulation is Sox19a. This SRY-related HMG box transcription factor shows:
- Significant upregulation (log₂FC = 1.45, *P*~adj~ = 0.044)
- 13 DMPs all in intergenic space (4–5 kb upstream)
- Consistent with silencer inactivation through methylation
- Known role in developmental regulation across metazoans

However, the MI = 0 result means that even Sox19a's regulation is not through the classical methylation → silencing pathway. The intergenic DMPs likely affect a silencer element, with methylation reducing its repressive activity — a mechanism that produces the opposite of the textbook prediction (methylation → activation).

### 3.4 Honest limitations

1. **Sample size.** Only 2 biological replicates per WGBS condition limits statistical power and prevents methylation heterogeneity analysis (PDR, epipolymorphism).

2. **Single time point.** We capture one snapshot of regeneration. Time-course WGBS would reveal temporal dynamics.

3. **No chromatin data.** Without ATAC-seq or ChIP-seq, we infer enhancer function from methylation alone. ATAC-seq would directly identify open regulatory elements.

4. **No causal validation.** DNMT inhibitor experiments, CRISPR perturbation, or reporter assays would test whether observed methylation changes are causal.

5. **No physical interaction data.** Hi-C or Capture-C would confirm that intergenic DMPs physically interact with target gene promoters.

6. **Weber criteria in invertebrates.** The HCP/ICP/LCP classification was designed for mammalian genomes. The 91.6% ICP result may partly reflect a need for invertebrate-calibrated thresholds.

7. **Genome quality.** Custom genome assembly without repeat masking or N-filtering in promoter analysis. However, the overall patterns are robust to threshold variation (sensitivity analysis across 25 threshold combinations).

### 3.5 Future directions

- **ATAC-seq** on regenerating tail tissue to identify enhancers at DMP clusters
- **Time-course WGBS** to reveal temporal methylation dynamics
- **5-azacytidine treatment** to test methylation requirement for regeneration
- **Co-methylation WGCNA** to compare methylation and expression modules
- **MethNet-style analysis** to link intergenic DMPs to target genes within 1 Mb windows
- **Methylation heterogeneity analysis** (Metheor tool) to distinguish uniform shifts from cell-type mixing

---

## 4. Methods

### 4.1 Organism and experimental design

*Deroceras laeve* (gray field slug) adults were subjected to tail amputation. Regenerating tail tissue was collected alongside control (intact) tail tissue. WGBS: 2 control (C1, C2) and 2 amputated (A1, A2) samples. RNA-seq: 6 tissues (tail, eye, body wall, head, juvenile, ovotestis) × 4 conditions (control, amputated, irradiated, fungicide), 24 libraries after outlier removal (T1S5 excluded).

### 4.2 WGBS analysis

Bismark v0.23+ (Krueger and Andrews, 2011) was used for alignment to the custom *D. laeve* genome (BSgenome.Dlaeve.NCBI.dlgm). Pipeline choice validated by the 2025 Pipeline Olympics (Nucleic Acids Research), which confirmed Bismark as having the best memory footprint among WGBS pipelines. BSseq objects constructed in R using the bsseq package (Hansen et al., 2012). DMPs called using DSS (Feng et al., 2014): beta-binomial model with arcsine link function and Wald test, FDR < 0.05. DMRs: areaStat threshold, minimum 3 CpGs within 500 bp, minimum length 50 bp. DSS was chosen because it provides proper dispersion estimation with as few as 2 replicates, outperforming alternatives (methylKit, Fisher's exact) in the 2021 DMR benchmarking study (Wu et al., 2021).

### 4.3 RNA-seq analysis

HTSeq count matrices were normalized using DESeq2 (Love et al., 2014), which uses negative binomial generalized linear models with dispersion shrinkage. Batch effects between experiments were removed using limma::removeBatchEffect (Ritchie et al., 2015) applied to VST-transformed counts for visualization and WGCNA, while keeping batch as a covariate in the DESeq2 model for differential expression testing. Differential expression: Wald test, FDR < 0.05.

### 4.4 Co-expression network analysis

WGCNA (Langfelder and Horvath, 2008) with signed adjacency matrix, soft-thresholding power selected for scale-free topology (R² > 0.85), blockwise module detection (maximum block size 20,000), minimum module size 30, merge cut height 0.25. 12 modules detected from 14,736 genes. Hub genes defined as top 10% by module membership (kME). WGCNA is the standard for co-expression network analysis (>15,000 citations) and has been applied to invertebrate transcriptomics in oysters, honeybees, and corals.

### 4.5 Genomic annotation

ChIPseeker (Yu et al., 2015) for annotation of DMPs and DMRs relative to genomic features, using the *D. laeve* EviAnn GFF annotation. This approach mirrors the DMRichR FAIRification framework for non-model organisms (Bioinformatics Advances, 2025).

### 4.6 Information theory and statistical testing

Mutual information: MI(X;Y) = Σ p(x,y) log₂[p(x,y) / p(x)p(y)], computed for X ∈ {hyper, hypo}, Y ∈ {up, down}. Permutation null distribution from 10,000 label shuffles.

Shannon entropy: H = −Σ pᵢ log₂(pᵢ), computed for DMP distribution across modules.

Cochran-Armitage trend test for monotonic trend in DE proportion across methylation magnitude bins. Fisher's exact tests with Benjamini-Hochberg correction for enrichment analyses. Spearman correlations throughout.

### 4.7 Machine learning

Logistic regression with |methylation difference| as predictor (baseline). Weighted linear and quadratic models for dosage response. Category-level enrichment using Fisher's exact tests on functional categories.

### 4.8 Transposable element analysis

TE coordinates and Kimura divergence from RepeatMasker output. Kimura-to-million-years conversion: Time(My) = kimura_div / (2 × 2.2 × 10⁻⁹ × 10⁶), using gastropod-specific neutral substitution rate (Kimura, 1980). CpG correction computed as kimura_div − kimura_div_CpG, reflecting historical methylation via C→T deamination at methylated CpGs. Cohen's *d* for effect size between conditions.

### 4.9 Promoter classification

Weber criteria (Weber et al., 2007): 500 bp sliding windows scanned from −2 kb to +500 bp relative to TSS. HCP: any window with CpG O/E ≥ 0.75 AND GC% ≥ 55%. LCP: no window with CpG O/E ≥ 0.48. ICP: remainder. 23,636 genes classified.

---

## 5. Data Availability

Raw sequencing data deposited in NCBI SRA [accession TBD]. Complete analysis code available at https://github.com/somnasgarden/dlaeve-regeneration-epigenomics. The repository includes 16 local analysis scripts, 8 cluster pipeline stages, and a comprehensive HTML report (report.html).

---

## 6. Acknowledgments

[TBD]

---

## 7. References

Aguilar, C. and Gardiner, D.M. (2015). DNA methylation dynamics regulate the formation of a regenerative wound epidermis during axolotl limb regeneration. *PLoS ONE* 10, e0118649.

Albert, R., Jeong, H., and Barabási, A.-L. (2000). Error and attack tolerance of complex networks. *Nature* 406, 378–382.

Bird, A. (2002). DNA methylation patterns and epigenetic memory. *Genes Dev.* 16, 6–21.

Cover, T.M. and Thomas, J.A. (2006). *Elements of Information Theory*, 2nd edn. Wiley.

de Mendoza, A., Poppe, D., Buckberry, S., et al. (2021). The emergence of the brain non-CpG methylation system in vertebrates. *Nat. Ecol. Evol.* 5, 369–378.

de Mendoza, A. et al. (2022). Genome-scale, base-resolution DNA methylation profiles across 580 animal species. *Nat. Commun.* 13, 5464.

Feng, H., Conneely, K.N., and Wu, H. (2014). A Bayesian hierarchical model to detect differentially methylated loci from single nucleotide resolution sequencing data. *Nucleic Acids Res.* 42, e69.

Fneich, S. et al. (2021). 5-Methyl-cytosine and 5-hydroxy-methyl-cytosine in the genome of *Biomphalaria glabrata*. *Epigenetics Chromatin* 14, 5.

Gavery, M.R. and Roberts, S.B. (2013). Predominant intragenic methylation is associated with gene expression characteristics in a bivalve mollusc. *PeerJ* 1, e215.

Geyer, K.K. et al. (2024). DNA methylation machinery in gastropod mollusks. *Life* 14, 537.

Hansen, K.D., Langmead, B., and Irizarry, R.A. (2012). BSmooth: from whole genome bisulfite sequencing reads to differentially methylated regions. *Genome Biol.* 13, R83.

Jaber-Hijazi, F. et al. (2013). Planarian MBD2/3 is required for adult stem cell pluripotency independently of DNA methylation. *Dev. Biol.* 384, 141–153.

Jones, P.A. (2012). Functions of DNA methylation: islands, start sites, gene bodies and beyond. *Nat. Rev. Genet.* 13, 484–492.

Kang, J. et al. (2020). Modulation of tissue repair by regeneration enhancer elements. *Genome Biol.* 21, 178.

Katsumura, T. et al. (2026). Comparative multi-omic analysis of fin and limb regeneration. *Nat. Commun.* 17, 68801.

Kim, G.B. et al. (2021). DeepTFactor: a deep learning-based tool for the prediction of transcription factors. *Proc. Natl Acad. Sci. USA* 118, e2021171118.

Kimura, M. (1980). A simple method for estimating evolutionary rates of base substitutions through comparative studies of nucleotide sequences. *J. Mol. Evol.* 16, 111–120.

Krueger, F. and Andrews, S.R. (2011). Bismark: a flexible aligner and methylation caller for bisulfite-seq applications. *Bioinformatics* 27, 1571–1572.

Langfelder, P. and Horvath, S. (2008). WGCNA: an R package for weighted correlation network analysis. *BMC Bioinformatics* 9, 559.

Levy-Jurgenson, A., Tekpli, X., Kristensen, V.N., and Yakhini, Z. (2020). Predicting gene expression from DNA methylation and transcript structure. *Nucleic Acids Res.* 48, 11292–11303.

Love, M.I., Huber, W., and Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biol.* 15, 550.

Pipeline Olympics (2025). Comprehensive benchmarking of WGBS computational pipelines. *Nucleic Acids Res.* 53, gkaf970.

Ritchie, M.E. et al. (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Res.* 43, e47.

Weber, M. et al. (2007). Distribution, silencing potential and evolutionary impact of promoter DNA methylation in the human genome. *Nat. Genet.* 39, 457–466.

Wu, H. et al. (2021). Benchmarking tools for detection of differentially methylated regions. *Brief. Bioinform.* 22, bbab168.

Yu, G., Wang, L.G., and He, Q.Y. (2015). ChIPseeker: an R/Bioconductor package for ChIP peak annotation, comparison and visualization. *Bioinformatics* 31, 2382–2383.

Zemach, A., McDaniel, I.E., Silva, P., and Zilberman, D. (2010). Genome-wide evolutionary analysis of eukaryotic DNA methylation. *Science* 328, 916–919.

Zhang, Y. et al. (2024). MethNet: integrating methylation and expression data to identify cis regulatory elements. *Nat. Commun.* 15, 50380.

---

## Figure Legends

**Figure 1. Methylation direction carries zero information about expression change.** (a) Mutual information between methylation direction (hyper/hypo) and expression direction (up/down): MI = 0.00009 bits (*P* = 0.55, 10,000 permutations). (b) Quadrant distribution showing near-uniform allocation across hyper-up (18.7%), hyper-down (24.2%), hypo-up (25.5%), and hypo-down (31.6%). (c) Cochran-Armitage trend test: no dose-response (z = 1.009, *P* = 0.31). (d) ML model comparison: all approaches fail to predict DE from methylation features (best AUC = 0.635).

**Figure 2. Intergenic methylation in a CpG-island-free genome.** (a) DMP genomic distribution: 31% intergenic, 32% intronic, 20% promoter, 17% other. (b) Promoter CpG classification: 91.6% ICP, 6.5% HCP, 1.9% LCP. (c) DMP spatial clustering: 73× enrichment within 1 kb. (d) Distance-decay of methylation concordance from 72.6% (10 kb) to 50% (>50 kb).

**Figure 3. Four-layer regulatory hierarchy.** (a) TF family DMP rates: Hox 0/27, homeodomain/HMG 1.7%, genome-wide 15.8%. (b) Module DMP enrichment: yellow consistently enriched (68% of thresholds). (c) Sox19a case study: intergenic DMPs with upregulation. (d) Zinc finger effectors enriched versus core TFs protected. (e) Schematic of the four-layer model.

**Figure 4. Hub gene paradox.** (a) Blue module hub DMP enrichment: 26.3% versus 16.5% (*P*~adj~ = 0.005). (b) kME versus |log₂FC|: positive correlation (ρ = +0.176, *P* = 0.007). (c) Largest DMRs: zero DE overlap. (d) TE methylation: static during regeneration (all *d* < 0.01).

**Figure 5. Cross-species comparison.** (a) CpG methylation levels across 9 species. (b) Three regeneration strategies: absent, static, topology-marking. (c) Feature comparison heatmap across regenerating species. (d) Paradigm shift: all 7 mammalian assumptions rejected in *D. laeve*.

**Figure 6. Module-level analysis.** (a) Eigengene heatmap across all samples. (b) Module-trait correlations. (c) Within-module methylation–expression correlations. (d) DMP enrichment across modules.

**Figure 7. Entropy and information theory.** (a) Shannon entropy of DMP distribution (H/H~max~ = 0.80). (b) Jensen-Shannon divergence between module methylation profiles. (c) Region-specific DE rates. (d) Network buffering contradicted: kME versus methylation response.

**Figure 8. TE evolutionary age analysis.** (a) TE age landscape spanning 0–207 My by class. (b) CpG correction revealing historical methylation. (c) DMP distribution by TE class. (d) Methylation levels unchanged across TE ages and conditions.

**Figure 9. Machine learning and predictive analysis.** (a) Logistic regression AUC = 0.635. (b) Dosage model comparison (linear R² = 0.29, quadratic R² = 0.33). (c) Category-level feature importance. (d) Module feature analysis.

**Figure 10. Topology-marking model synthesis.** (a) Complete model: methylation marks network architecture. (b) Comparison with mammalian and invertebrate paradigms. (c) Global versus within-module correlation showing cancellation effect. (d) Regulatory domain model with spatial clustering.

---

## Supplementary Tables

**Table S1.** DMP statistics (18,754 DMPs): coordinates, methylation difference, FDR, genomic annotation, module assignment.
**Table S2.** DMR statistics (1,424 DMRs): coordinates, areaStat, number of CpGs, mean methylation difference.
**Table S3.** Module DMP enrichment (Fisher's exact tests) across 25 threshold combinations.
**Table S4.** TF family methylation analysis with DeepTFactor validation.
**Table S5.** Hub gene methylation status for all 12 modules.
**Table S6.** Shannon entropy, KL divergence, and Jensen-Shannon divergence metrics.
**Table S7.** Machine learning model performance comparison.
**Table S8.** TE age analysis by class (3.4M elements).
**Table S9.** Cross-species methylation comparison (9 species, 5 phyla).
**Table S10.** Paradigm contrast table: mammalian versus *D. laeve*.
**Table S11.** Candidate genes ranked by composite evidence score.
**Table S12.** Promoter CpG classification (23,636 genes).
**Table S13.** DMP spatial clustering statistics.
**Table S14.** Module methylation direction asymmetry.
**Table S15.** Four-layer model evidence matrix.
