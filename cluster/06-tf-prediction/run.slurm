#!/usr/bin/env bash
#SBATCH --job-name=tf-prediction
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err
#SBATCH --mem=16G

set -euo pipefail

# DeepTFactor prediction (Python):
#   sbatch run.slurm deeptfactor/tf_running.py
#
# TF validation analysis (R):
#   sbatch run.slurm 01_tf_validation.R

SCRIPT="${1:?Usage: sbatch run.slurm <script>}"

if [[ "$SCRIPT" == *.py ]]; then
    python "$SCRIPT"
elif [[ "$SCRIPT" == *.R ]]; then
    Rscript "$SCRIPT"
fi
