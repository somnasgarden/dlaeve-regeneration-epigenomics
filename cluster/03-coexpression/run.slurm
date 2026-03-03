#!/usr/bin/env bash
#SBATCH --job-name=wgcna
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

set -euo pipefail

# WGCNA needs more memory and threads for TOM computation
#   sbatch run.slurm 01_wgcna_module_detection.R

RSCRIPT="${1:?Usage: sbatch run.slurm <script.R>}"
Rscript "$RSCRIPT"
