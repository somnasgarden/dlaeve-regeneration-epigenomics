# Bisulfite Alignment Pipeline

Scripts to preprocess and align whole-genome bisulfite sequencing (WGBS) reads using Bismark. The pipeline processes 4 paired-end samples (A1r, A2r, C1, C2r) through quality control, alignment, deduplication, and methylation extraction.

## Pipeline Overview

Raw FASTQ → Cutadapt → Trimmomatic → Bismark Alignment → Coverage/QC → Deduplication → Methylation Extraction → Cytosine Report
---

## Scripts

### `00_cutadapt.sh` — Hard-clip 5′ bases
Removes the first 12 bases from each read using `cutadapt -u 12`. This accounts for random priming bias common in WGBS library preparation. Processes R1 and R2 files for all 4 samples independently.

**Output:** `*_cut.fq.gz` files in `03_cutadapt/`

---

### `01_trimmomatic.sh` — Adapter trimming and quality filtering
Runs `trimmomatic PE` on the cutadapt-clipped reads to remove Illumina TruSeq3 adapters and low-quality bases. Parameters: `LEADING:30 TRAILING:30 SLIDINGWINDOW:4:30 MINLEN:120`. Runs FastQC on all paired outputs for quality assessment.

**Output:** `trimmed_*_R{1,2}.paired.fq.gz` in `05_trimmomatic/`

---

### `01_prepare_bismark_genome.sh` — Genome indexing
Prepares the reference genome for bisulfite alignment using `bismark_genome_preparation`. Runs with `--parallel 4` (2 internal processes × 4 threads = 8 CPUs total) and generates genomic composition statistics.

**Output:** Bismark/Bowtie2 index files in the genome folder (`50-Metilacion/Genoma/`)

---

### `02_alignment_bismark.sh` — Bismark alignment
Aligns all trimmed paired-end reads to the prepared genome using `bismark --parallel 4` with Bowtie2. Automatically collects all R1/R2 paired files from the trimmomatic output directory using `find`.

**Output:** BAM files in `06_bismark_alignment/`

---

### `03_aln_coverage_report.sh` — Nucleotide coverage report (array job)
SLURM array job (1–4) running `bam2nuc` per sample to generate nucleotide coverage and composition statistics from each alignment BAM. Useful for QC of bisulfite conversion efficiency.

**Output:** Per-sample reports in `07_methylation_extraction_and_reports/<SAMPLE_ID>_output/`

---

### `04_deduplicate_bams.sh` — PCR duplicate removal (array job)
SLURM array job (1–4) running `deduplicate_bismark` per sample to remove PCR duplicates introduced during library amplification, which is critical for accurate methylation quantification.

**Output:** Deduplicated BAMs in `08_deduplication/`

---

### `05_extract_methylation.sh` — Methylation extraction
Runs `bismark_methylation_extractor` on all deduplicated BAMs. Uses `--parallel 8`, `--ignore_r2 2` (skips the first 2 bases of R2 to avoid end-repair bias), and generates gzipped CpG/non-CpG coverage files and a bedGraph.

**Output:** Methylation call files and cytosine reports in `09_methylation_calls/`

---

### `06_standalone_coverage2cytosine.sh` — Cytosine report (standalone)
A flexible, CLI-driven wrapper around `coverage2cytosine` that can be called individually per sample. Accepts `--input`, `--dir`, `--out`, and `--genome` arguments. Useful for re-running or debugging a single sample without resubmitting the full pipeline.

**Usage:**
```bash
sbatch 06_standalone_coverage2cytosine.sh \
  -i 09_methylation_calls/C2_paired_bismark_bt2_pe.deduplicated.bismark.cov.gz \
  -d 09_methylation_calls \
  -o C2
```

**Output:** Full cytosine-context report (`.CX_report.txt.gz`) in the specified output directory.

---

## Directory Structure

```
.
├── 01.RawData/            # Raw FASTQ inputs
├── 03_cutadapt/           # Hard-clipped reads
├── 05_trimmomatic/        # Adapter-trimmed paired reads
├── 06_bismark_alignment/  # Aligned BAMs
├── 07_methylation_extraction_and_reports/  # bam2nuc QC reports
├── 08_deduplication/      # Deduplicated BAMs
├── 09_methylation_calls/  # Final methylation output
├── 00_Genoma/             # Reference genome and Bismark index
└── slurm_logs/            # SLURM stdout/stderr logs
```

## Samples

| Sample ID | Condition |
|-----------|-----------|
| A1r       | Amputated tail, replicate 1 |
| A2r       | Amputated tail, replicate 2 |
| C1        | Control, replicate 1 |
| C2r       | Control, replicate 2 |

## Software Versions

| Tool | Version |
|------|---------|
| Bismark | 0.25.1 |
| Bowtie2 | 2.5.4 |
| Samtools | 1.22.1 |
| Cutadapt | 5.1 |
| Trimmomatic | 0.33 |
| TrimGalore | 0.6.10 |
| FastQC | 0.11.3 |
