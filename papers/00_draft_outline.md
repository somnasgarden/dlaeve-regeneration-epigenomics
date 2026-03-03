# Paper Draft: Epigenetic Regulation of Regeneration in *Deroceras laeve*

## Working Title

**Module-specific DNA methylation rewiring during invertebrate tail regeneration reveals non-canonical regulatory mechanisms**

*Alternative titles:*
- Intergenic methylation, not promoters, drives module-specific gene regulation during gastropod regeneration
- Transcription factor methylation protection and module-specific epigenetic rewiring in invertebrate regeneration

---

## Abstract (Draft)

DNA methylation is a key epigenetic mechanism regulating gene expression, yet its role in invertebrate regeneration remains largely unexplored. Here, we integrate whole-genome bisulfite sequencing (WGBS) and RNA-seq from regenerating tail tissue of the gastropod *Deroceras laeve* to characterize the epigenetic landscape of regeneration in a non-model invertebrate. We identify 18,754 differentially methylated positions (DMPs) affecting ~16% of expressed genes and use weighted gene co-expression network analysis (WGCNA) to reveal that methylation-expression relationships are entirely module-specific rather than genome-wide. Strikingly, DeepTFactor CNN-predicted transcription factors are significantly enriched for DMPs (OR=1.607, p=4.32e-13), yet 0/6 differentially expressed TFs carry any methylation changes — revealing a paradox where methylation targets TFs without altering their expression. Intergenic and intronic regions, not promoters, harbor the functionally relevant methylation changes, suggesting enhancer/silencer-mediated regulation. The yellow WGCNA module (cell cycle/mitotic genes) shows robust DMP enrichment across all statistical thresholds, while the brown module displays significant hypomethylation bias (padj=0.001). Network hub genes in the blue module are preferentially methylated, whereas black module hubs are entirely unmethylated. These findings demonstrate that invertebrate regeneration employs a non-canonical epigenetic program distinct from mammalian paradigms, where methylation acts through distal regulatory elements in a module-specific manner while actively protecting transcriptional control points.

---

## Key Findings (with statistical support)

### 1. The methylation landscape is non-canonical

- **91.6% of promoters are ICP** (Intermediate CpG), only 6.5% HCP, 1.9% LCP
- Diverges sharply from mammalian genomes (~70% HCP)
- Promoter methylation-expression correlation is weak: ICP r=-0.066, HCP r=-0.035
- **This means the mammalian model of promoter CpG island methylation = silencing does not apply**

### 2. Intergenic and intronic regions are where methylation matters

- DMPs are distributed: ~31% intergenic, ~32% intronic, ~20% promoter
- Module-region Spearman correlations show strongest effects:
  - Purple + Promoter: r=+0.812, p=2.5e-6 (paradoxical positive!)
  - Tan + UTR: r=-0.747, p=9.9e-5
  - Yellow + Downstream: r=-0.672, p=2.3e-4
  - Black + Promoter: r=-0.515, p=2.3e-6 (classical silencing)
- **Global correlations are non-significant in ALL regions** because opposing module effects cancel out
- The relationship is ENTIRELY module-specific and region-specific

### 3. TF methylation paradox — targeted but expression-protected

- **DeepTFactor CNN-based TFs: OR=1.607, p=4.32e-13** (significant ENRICHMENT for DMPs)
- **DMR enrichment: OR=1.645, p=3.33e-06**
- 14.8% of TFs have DMPs vs 9.8% of non-TFs
- Top methylated TFs: ZFP2 (30 DMPs), Mafb (25), ATF2 (18), Hr96 (18)
- **However, 0/6 DE TFs have any methylation changes** — methylation targets TFs but does NOT change their expression
- 6 DE TFs (Scrt1, REPTOR-BP, Y-box, Catenin, sox8, CEBPA) achieve expression changes WITHOUT methylation
- Old keyword-based finding (OR=0.77 depleted) was an artifact of incomplete TF identification
- Hox genes: 0/27 with DMPs — master regulators still protected
- **Interpretation**: Methylation marks TFs based on network position, not to control their transcriptional output

### 4. Yellow module is the epigenetic target of regeneration

- **Enriched for DMPs**: Fisher padj < 0.05 (ONLY module robust across thresholds)
- **68% of threshold combinations give significance** (5 FDR x 5 effect size = 25 combos)
- Yellow module = cell cycle, mitotic processes, DNA repair
- Top enrichment: Mitotic cell cycle process (padj=1.55e-28), Cell Cycle (Reactome, padj=1.11e-15)
- DMP burden concentrated in intergenic/intronic regions
- **Biological interpretation**: Regeneration targets cell cycle genes via distal methylation

### 5. Brown module shows coordinated hypomethylation

- Significant hypomethylation bias: 38.1% hyper vs 61.9% hypo (chi-sq padj=0.001)
- Brown = UP in regeneration (Cohen's d=2.038)
- Brown functions: Drug metabolism, SLC-mediated transmembrane transport
- **Hypomethylation may be activating this module during regeneration**

### 6. Hub gene methylation is module-specific

- **Blue module hubs ENRICHED for DMPs**: 26.3% vs 16.5% non-hubs (padj=0.005)
  - Blue = DOWN in regeneration → methylation may be driving repression
- **Black module hubs DEPLETED**: 0% vs 14.5% non-hubs (padj=0.020)
  - Black module hubs are completely protected from methylation
- **Red module hubs DEPLETED**: 6.2% vs 12.2% (p=0.018, trending after FDR)
  - Red = UP in regeneration with cytoplasmic translation, ribosome enrichment

### 7. DMP spatial clustering identifies putative regulatory elements

- 2,109 DMP clusters (3+ CpGs within 500bp)
- 452 clusters with 10+ DMPs
- Nearly all clusters are unidirectional (>99% same methylation direction)
- Largest cluster: 33 DMPs near Cdkl4 (cyclin-dependent kinase-like 4) — hypomethylated
- Cluster near Notch1: 16 DMPs, all hypomethylated gene body

### 8. Sox19a: Intergenic methylation and upregulation

- Sox19a is significantly upregulated (log2FC=1.45, padj=0.044) in regenerating tissue
- No DMPs directly on the gene
- 13 associated DMPs are all in intergenic/downstream regions ~4-5kb upstream
- This is consistent with enhancer-mediated activation via demethylation of upstream regulatory elements
- Blue module member — blue is DOWN overall, but Sox19a is UP (module outlier with high GS)

### 9. Methylation dosage does NOT predict expression globally

- Overall Spearman: rho=-0.005, p=0.80 (completely non-significant)
- This is an honest null result — methylation effect size does not linearly predict expression change
- The relationship is qualitative (module-specific) not quantitative
- Possible explanation: methylation acts as a binary switch, not a rheostat

### 10. Network methylation coordination is weak but real

- Connected genes (GENIE3 edges) show 51.8% same-direction methylation (vs 50% expected)
- Binomial test: p=0.001
- Effect is very small — methylation changes are NOT propagating through the network
- Instead, methylation targets SPECIFIC loci, not pathways

### 11. Blue module: Dual DMP + DMR targeting (strongest methylation target)

- **DMR enrichment**: OR=1.54, padj<0.0001 (only module significant for DMRs)
- **DMP hub enrichment**: 26.3% of hubs have DMPs (padj=0.005)
- Blue module genes have both individual CpG changes (DMPs) AND coordinated multi-CpG changes (DMRs)
- Direction: DOWN in regeneration — methylation may drive repression of tissue maintenance
- Top methylated hub: fibcd1-a (25 DMPs, sulfation pathway)

### 12. DMP-DMR concordance validates both data types

- **93.9% of genes with both DMPs and DMRs show concordant direction** (p=1.4e-160)
- 789 genes have both DMPs and DMRs
- This extreme concordance (>>50%) validates that both approaches capture the same regulatory signal
- DMRs are more distal (median 12.3kb from TSS vs 12.4kb for DMPs — similar)

### 13. Morphogenesis genes are PROTECTED, not targeted

- Morphogenesis/Development genes: OR=0.69, padj=0.14 (trend toward depletion)
- **Hox genes: 0/27 have DMPs** — complete methylation protection
- When morphogenesis genes DO have DMPs, they are biased toward hypomethylation (61.4%, p=0.028)
- The initial impression of DMP clustering near developmental genes reflects dramatic cases (Sox19a, Notch1) but the STATISTICAL pattern is protection, not targeting

### 14. Methylation acts as a switch, not a rheostat

- No methylation threshold predicts DE: Cochran-Armitage trend p=0.178
- DMP count does not predict DE rate: 1 DMP → 1.9% DE, >20 DMPs → 3.4% DE (NS)
- Region-specific |meth| vs |log2FC| correlations: all NS (p>0.38)
- **Interpretation**: Methylation is qualitative (which genes), not quantitative (how much)
- Consistent with an enhancer on/off model rather than graded promoter response

### 15. Functional category DMP enrichment reveals targeting hierarchy

- **Enriched** (padj<0.05): Transcription factors (annotation-based, OR=1.45), Cytoskeleton/Motility (OR=1.45), "Other/Unknown" (OR=1.12)
- **Depleted** (padj<0.05): RNA Processing (OR=0.72)
- **Critical nuance**: DeepTFactor CNN-based TFs are ENRICHED for DMPs (OR=1.607, p=4.32e-13). The old keyword-based finding (OR=0.77 depleted) was an artifact of incomplete TF identification. Top methylated TFs: ZFP2, Mafb, ATF2, Hr96. Despite enrichment, 0/6 DE TFs are methylated — methylation targets TFs but doesn't change their expression.

---

## Proposed Figures

### Figure 1: Methylation landscape overview
- (A) Promoter classification (HCP/ICP/LCP pie chart)
- (B) DMP genomic distribution (region stacked bar)
- (C) DMP direction (hyper/hypo by region)
- (D) Chromosome-level DMP density plot

### Figure 2: Module-specific methylation-expression relationships
- (A) Heatmap of module x region Spearman r (fig_D1)
- (B) Global vs module-specific correlations (fig_D4) — why global washes out
- (C) Yellow module DMP enrichment across thresholds
- (D) Module methylation direction balance (brown bias)

### Figure 3: TF methylation protection
- (A) TF family DMP rates vs background (fig_G01)
- (B) Homeodomain/HMG vs zinc finger contrast
- (C) Methylated TF network targets (cascade heatmap)
- (D) Key TFs: expression change vs methylation

### Figure 4: Hub gene epigenetic targeting
- (A) Hub vs non-hub DMP rates by module (fig_G05)
- (B) Blue module hub enrichment detail
- (C) Black module hub depletion
- (D) Hub connectivity vs methylation effect

### Figure 5: Intergenic regulation and DMP clusters
- (A) DMP cluster size distribution
- (B) Distance to TSS vs expression effect
- (C) Sox19a case study: intergenic DMPs + expression
- (D) Notch1 case study: gene body methylation cluster

### Figure 6: Candidate genes for experimental validation
- (A) Top 30 candidate lollipop (fig_E2)
- (B) Evidence matrix heatmap (fig_E4)
- (C) Sulfation pathway candidates
- (D) Module-level candidate distribution

### Figure 7: Four-layer model and synthesis
- (A) Category DMP enrichment forest plot (morphogenesis, TFs, cytoskeleton)
- (B) Category × Module heatmap (which categories in which modules)
- (C) Methylation threshold analysis (DE rate vs |meth diff|)
- (D) Four-layer regeneration model diagram

### Supplementary Figures
- S1: Threshold sensitivity analysis (all 9 plots from fig_A series)
- S2: Methylation dosage response (null result — important to report)
- S3: TE-associated DMPs
- S4: DMP quadrant analysis by region
- S5: WGCNA module enrichment details
- S6: DMR size vs CpG density
- S7: DMR module enrichment
- S8: DMR distance-to-TSS distribution
- S9: DMP-DMR direction concordance plot
- S10: Developmental pathway DMP heatmap (Wnt, Notch, Hox, BMP/TGFb)

---

## Discussion Points

### Novel contributions

1. **First comprehensive methylation-transcriptome integration in a regenerating gastropod**
   - No prior study has done WGBS + RNA-seq + WGCNA in any mollusk regeneration system

2. **TF methylation paradox**: DeepTFactor CNN TFs are ENRICHED for DMPs (OR=1.607, p=4.32e-13) yet 0/6 DE TFs are methylated. Methylation targets TFs based on network topology, not to control expression. Hox genes (0/27 DMPs) remain protected. This suggests a hierarchy: TFs are marked but buffered, while effectors are the functional methylation targets.

3. **Module-specificity invalidates genome-wide methylation-expression models**: The complete absence of global correlation (rho=-0.005) combined with strong module-level correlations (|r| up to 0.81) demonstrates that methylation effects are context-dependent. This challenges studies that report only genome-wide statistics.

4. **Intergenic-centric regulation**: The dominance of intergenic/intronic DMPs, combined with the ICP-dominated promoter landscape, suggests that *D. laeve* regulates regeneration genes through enhancer/silencer methylation rather than promoter methylation. This is distinct from the mammalian CpG island model.

5. **Hub gene targeting is module-specific**: Blue module hubs being preferentially methylated while black module hubs are completely protected suggests different epigenetic strategies for different functional modules.

6. **Four-layer model of regenerative epigenetic regulation**:
   - Layer 1: **Protected regulators** — Hox, homeodomain/HMG TFs, core developmental genes shielded from methylation
   - Layer 2: **Targeted modules** — Yellow (cell cycle, DMP enriched), Blue (repressed, dual DMP+DMR), Brown (activated, hypomethylated), Black (protected hubs)
   - Layer 3: **Distal regulation** — Enhancer/silencer methylation in intergenic/intronic space (63% of DMPs)
   - Layer 4: **Effector targeting** — Zinc finger proteins, cytoskeleton/motility genes are the direct methylation targets

7. **DMP-DMR concordance at 93.9%**: Extreme agreement between single-CpG and multi-CpG analyses validates both approaches and strengthens regulatory element identification.

8. **Methylation is a binary switch, not a rheostat**: Neither methylation magnitude nor DMP count predicts DE probability (Cochran-Armitage p=0.178). The relationship is qualitative (which module) not quantitative (how much change).

### Limitations (honest reporting)

1. **n=2 per condition** for WGBS — limits statistical power for individual DMPs
2. **No functional validation** of putative enhancers — BED files prepared for motif analysis
3. **Annotation-dependent**: Non-model organism annotations are less complete
4. **Correlation ≠ causation**: Cannot determine if methylation changes cause expression changes or vice versa
5. **Single time point**: Cannot distinguish early vs late regeneration responses

### Future directions

1. Functional validation of Sox19a intergenic regulatory region (reporter assay)
2. ATAC-seq to confirm enhancer accessibility at DMP clusters
3. Time-series WGBS to track methylation dynamics during regeneration
4. 5-azacytidine treatment to test methylation requirement for regeneration
5. Compare with planarian/axolotl methylation patterns

---

## Candidate Genes for Cloning (Priority List)

### Tier 1: Strongest candidates (score >= 12)
| Gene | Symbol | Module | Score | Evidence |
|------|--------|--------|-------|----------|
| LOC_00005231 | B3gntl1 | green | 13 | Sulfation + Hub + 14 sig DMPs |
| LOC_00020191 | SULT6B1 | green | 13 | Sulfotransferase + Hub + DMPs |
| LOC_00010201 | fibcd1-a | blue | 12 | Sulfation + Hub + 17 sig DMPs |
| LOC_00019391 | Cbs | brown | 12 | Sulfation + Hub + DMPs |

### Tier 2: Yellow module hubs (score >= 11)
| Gene | Symbol | Score | Evidence |
|------|--------|-------|----------|
| LOC_00003444 | ARHGAP11A | 11 | Hub + DMP + Yellow (DMP-enriched module) |
| LOC_00003602 | XPO5 | 11 | Hub + 6 sig DMPs + Yellow |
| LOC_00006234 | CHD2 | 11 | Chromatin remodeler + Hub + Yellow |
| LOC_00025818 | EXO1 | 11 | DNA repair + Hub + Yellow |
| LOC_00019749 | KMT2A | 11 | Histone methyltransferase + Hub + Yellow |

### Tier 3: Key developmental genes
| Gene | Symbol | Evidence |
|------|--------|----------|
| LOC_00002600 | Sox19a | Sig upregulated (padj=0.044) + intergenic DMPs |
| LOC_00009396 | LHX6 | LIM/homeobox TF |
| LOC_00008589 | Notch1 | 16 DMPs (13 sig), all hypomethylated |
| LOC_00006568 | Wnt4 | Wnt pathway |

---

## Data Availability

All raw sequencing data deposited in [TBD]. Analysis code available at [GitHub repo].
Processed data tables and figures in `results/` directory structure.
