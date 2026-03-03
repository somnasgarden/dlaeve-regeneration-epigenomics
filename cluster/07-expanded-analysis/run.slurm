#!/usr/bin/env bash
#SBATCH --job-name=expanded
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --time=08:00:00

set -euo pipefail

# Expanded analyses — require cluster for:
#   - GENIE3 500K edge loading (memory)
#   - Permutation testing (CPU time)
#   - Bootstrap confidence intervals (CPU time)
#
# Run each script:
#   sbatch run.slurm 01_permutation_enrichment.R
#   sbatch run.slurm 02_network_null_models.R
#   sbatch run.slurm 03_bootstrap_correlations.R

RSCRIPT="${1:?Usage: sbatch run.slurm <script.R>}"
Rscript "$RSCRIPT"
