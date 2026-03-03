#!/usr/bin/env bash
#SBATCH --job-name=mxt-integration
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err
#SBATCH --mem=32G

set -euo pipefail

# Run each script in order:
#   sbatch run.slurm 01_methylation_vs_expression.R
#   sbatch run.slurm 02_module_methylation_enrichment.R

RSCRIPT="${1:?Usage: sbatch run.slurm <script.R>}"
Rscript "$RSCRIPT"
