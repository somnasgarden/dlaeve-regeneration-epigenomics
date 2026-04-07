#!/bin/bash
# Run job through bash shell
#
# Your job name
#SBATCH -J wil_slurm
#
# Use current working directory
#SBATCH --chdir=./
#
#stdout and stderr
#SBATCH --output=%x.o%j
#
#
# If modules are needed, source modules environment (Do not delete the next line):
. /etc/profile.d/modules.sh
#
# Add any modules you might require:
module load trimgalore/0.6.10
module load fastqc/0.11.3
module load anaconda3/2025.06
source activate cutadapt-5.1
#
# here.
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#
# Write your commands in the next line

cutadapt -u 12 -o A1r_1_cut.fq.gz ../01.RawData/A1r/A1r_1.fq.gz
cutadapt -u 12 -o A1r_2_cut.fq.gz ../01.RawData/A1r/A1r_2.fq.gz
cutadapt -u 12 -o A2r_1_cut.fq.gz ../01.RawData/A2r/A2r_1.fq.gz
cutadapt -u 12 -o A2r_2_cut.fq.gz ../01.RawData/A2r/A2r_2.fq.gz
cutadapt -u 12 -o C1_1_cut.fq.gz ../01.RawData/C1/C1_1.fq.gz
cutadapt -u 12 -o C1_2_cut.fq.gz ../01.RawData/C1/C1_2.fq.gz
cutadapt -u 12 -o C2r_1_cut.fq.gz ../01.RawData/C2r/C2r_1.fq.gz
cutadapt -u 12 -o C2r_2_cut.fq.gz ../01.RawData/C2r/C2r_2.fq.gz
