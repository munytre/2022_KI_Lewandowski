#!/bin/bash
 
#SBATCH -A snic2022-22-143
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 72:00:00
#SBATCH -J fastq_dump

# Error report
set -euo pipefail

# Load required tools (Detection)
module load bioinfo-tools
module load sratools/2.10.9 

# Starting time
echo -e "\n`date` Start"

# Working directory
wd="$1"
cd "${wd}"

# Download
echo -e "\n`date` Start download"
fastq-dump SRR13782529_GSM5106143_02_E4_15_H09_Mus_musculus_RNA-Seq.sra --gzip --split-files -v

echo -e "\n`date` Done with download"