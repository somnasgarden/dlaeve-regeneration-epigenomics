#!/usr/bin/env bash
#SBATCH --job-name=methylation
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err
#SBATCH --mem=32G

set -euo pipefail

# Run each script in order:
#   sbatch run.slurm 01_quality_control.R
#   sbatch run.slurm 02_differential_methylation.R
#   sbatch run.slurm 03_genomic_annotation.R
#   sbatch run.slurm 04_functional_enrichment.R
#   sbatch run.slurm 05_te_methylation.R

RSCRIPT="${1:?Usage: sbatch run.slurm <script.R>}"
Rscript "$RSCRIPT"
