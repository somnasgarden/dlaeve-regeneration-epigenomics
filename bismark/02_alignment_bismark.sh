#!/bin/bash

#SBATCH --job-name=bismark_aln      # Job name
#SBATCH --output=oalnbismark_%j.out  # Standard output log
#SBATCH --error=oalnbismark_%j.err   # Standard error log
#SBATCH --nodes=1                   # Run on a single node
#SBATCH --ntasks=1                  # Run a single task
#SBATCH --cpus-per-task=20           # Request 8 CPUs for this task
#SBATCH --mem=100G                  # Request 100GB of memory


echo "Starting Bismark alignment..."
echo "Job ID: $SLURM_JOB_ID"
echo "Running on host: $(hostname)"
echo "Allocated CPUs: $SLURM_CPUS_PER_TASK"
echo "Allocated Memory: $SLURM_MEM_PER_NODE"

# --- Load Required Modules ---
# Uncomment and update these lines if Bismark or Bowtie2 are
# loaded as modules on your HPC cluster.
#
module load bowtie2/2.5.4 samtools/1.22.1
export PATH=/mnt/data/.../80-scripts/81-bin/Bismark-0.25.1:$PATH
export PERL5LIB=$HOME/.perl/lib:$PERL5LIB

fastq_dir=$(pwd)/"05_trimmomatic/"
reads_1_pat="R1.paired.fq.gz"
reads_2_pat="R2.paired.fq.gz"
reads1="$(find -L ${fastq_dir} -name "*$reads_1_pat" -type f | sort | awk -v OFS=, '{$1=$1;print}' | paste -sd,)"
reads2="$(find -L ${fastq_dir} -name "*$reads_2_pat" -type f | sort | awk -v OFS=, '{$1=$1;print}' | paste -sd,)"

# --- Run the Command ---
bismark --parallel 4 -o 06_bismark_alignment 00_Genoma -1 $reads1 -2 $reads2

echo "Bismark alignment finished with exit code $?."
