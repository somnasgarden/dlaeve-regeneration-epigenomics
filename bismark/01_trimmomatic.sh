#!/bin/bash
# 
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
module load trimmomatic/0.33
#
# 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#
# Write your commands in the next line

trimmomatic PE \
   ../03_cutadapt/A1r_1_cut.fq.gz \
   ../03_cutadapt/A1r_2_cut.fq.gz \
   ./trimmed_A1_R1.paired.fq.gz \
   ./trimmed_A1_R1.unpaired.fq.gz \
   ./trimmed_A1_R2.paired.fq.gz \
   ./trimmed_A1_R2.unpaired.fq.gz \
   ILLUMINACLIP:TruSeq3_Illumina.fa:2:30:10 LEADING:30 TRAILING:30 SLIDINGWINDOW:4:30 MINLEN:120

trimmomatic PE \
   ../03_cutadapt/A2r_1_cut.fq.gz \
   ../03_cutadapt/A2r_2_cut.fq.gz \
   ./trimmed_A2_R1.paired.fq.gz \
   ./trimmed_A2_R1.unpaired.fq.gz \
   ./trimmed_A2_R2.paired.fq.gz \
   ./trimmed_A2_R2.unpaired.fq.gz \
   ILLUMINACLIP:TruSeq3_Illumina.fa:2:30:10 LEADING:30 TRAILING:30 SLIDINGWINDOW:4:30 MINLEN:120

trimmomatic PE \
   ../03_cutadapt/C1_1_cut.fq.gz \
   ../03_cutadapt/C1_2_cut.fq.gz \
   ./trimmed_C1_R1.paired.fq.gz \
   ./trimmed_C1_R1.unpaired.fq.gz \
   ./trimmed_C1_R2.paired.fq.gz \
   ./trimmed_C1_R2.unpaired.fq.gz \
   ILLUMINACLIP:TruSeq3_Illumina.fa:2:30:10 LEADING:30 TRAILING:30 SLIDINGWINDOW:4:30 MINLEN:120

trimmomatic PE \
   ../03_cutadapt/C2r_1_cut.fq.gz \
   ../03_cutadapt/C2r_2_cut.fq.gz \
   ./trimmed_C2_R1.paired.fq.gz \
   ./trimmed_C2_R1.unpaired.fq.gz \
   ./trimmed_C2_R2.paired.fq.gz \
   ./trimmed_C2_R2.unpaired.fq.gz \
   ILLUMINACLIP:TruSeq3_Illumina.fa:2:30:10 LEADING:30 TRAILING:30 SLIDINGWINDOW:4:30 MINLEN:120

fastqc *paired.fq.gz
