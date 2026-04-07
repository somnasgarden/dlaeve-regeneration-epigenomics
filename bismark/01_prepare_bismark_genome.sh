#!/bin/bash

#SBATCH --job-name=bismark_prep      # Job name
#SBATCH --output=obismark_prep_%j.out  # Standard output log
#SBATCH --error=obismark_prep_%j.err   # Standard error log
#SBATCH --nodes=1                   # Run on a single node
#SBATCH --ntasks=1                  # Run a single task
#SBATCH --cpus-per-task=8           # Request 8 CPUs for this task
#SBATCH --mem=100G                  # Request 100GB of memory

#  Bismark runs 2 indexing processes in parallel (top/bottom strand).
#  --parallel 4, which gives 4 threads *to each* of those processes.
# Total CPUs = 2 processes * 4 threads/process = 8 CPUs.


echo "Starting Bismark genome preparation..."
echo "Job ID: $SLURM_JOB_ID"
echo "Running on host: $(hostname)"
echo "Allocated CPUs: $SLURM_CPUS_PER_TASK"
echo "Allocated Memory: $SLURM_MEM_PER_NODE"

# --- Load Required Modules ---
module load bamtools/2.5.1 bowtie2/2.5.4

export PATH=/mnt/.../80-scripts/81-bin/Bismark-0.25.1:$PATH
export PERL5LIB=$HOME/.perl/lib:$PERL5LIB

bismark_genome_preparation \
  --verbose \
  /mnt/data/../50-Genoma/51-Metilacion/Genoma/ \
  --parallel 4 --genomic_composition

echo "Bismark preparation finished with exit code $?."
