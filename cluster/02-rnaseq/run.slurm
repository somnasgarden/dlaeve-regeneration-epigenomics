#!/usr/bin/env bash
#SBATCH --job-name=rnaseq
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err
#SBATCH --mem=32G

set -euo pipefail

# Run each script in order:
#   sbatch run.slurm 01_normalization.R
#   sbatch run.slurm 02_batch_correction.R
#   sbatch run.slurm 03_differential_expression.R
#   sbatch run.slurm 04_hierarchical_clustering.R
#   sbatch run.slurm 05_kmeans_clustering.R

RSCRIPT="${1:?Usage: sbatch run.slurm <script.R>}"
Rscript "$RSCRIPT"
