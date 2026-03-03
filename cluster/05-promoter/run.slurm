#!/usr/bin/env bash
#SBATCH --job-name=promoter
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err
#SBATCH --mem=32G

set -euo pipefail

#   sbatch run.slurm 01_promoter_classification.R

RSCRIPT="${1:?Usage: sbatch run.slurm <script.R>}"
Rscript "$RSCRIPT"
