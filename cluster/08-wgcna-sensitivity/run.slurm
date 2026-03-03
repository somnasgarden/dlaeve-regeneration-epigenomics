#!/usr/bin/env bash
#SBATCH --job-name=wgcna-sens
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00

set -euo pipefail

# WGCNA Sensitivity & Methylation-Network Integration
#
# Run scripts sequentially (each depends on the previous):
#   sbatch run.slurm 01_parameter_sensitivity.R    # ~3-5 hours
#   sbatch run.slurm 02_improved_wgcna.R           # ~1-2 hours
#   sbatch run.slurm 03_methylation_network_integration.R  # ~2-4 hours
#
# Independent exploration (does NOT depend on 01-03):
#   sbatch run.slurm 04_enrichment_across_parameters.R     # ~6-10 hours
#   Tests 7 variance filters × 9 powers, runs enrichment per config,
#   tracks sulfation/GAG term emergence, connectivity vs enrichment

RSCRIPT="${1:?Usage: sbatch run.slurm <script.R>}"
Rscript "$RSCRIPT"
