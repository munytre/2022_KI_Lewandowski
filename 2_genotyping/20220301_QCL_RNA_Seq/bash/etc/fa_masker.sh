#!/bin/bash
 
#SBATCH -A snic2022-22-143
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 72:00:00
#SBATCH -J fasta_masker

# Error report
set -euo pipefail

# Load required tools)
module load bioinfo-tools
module load BEDTools/2.29.2

# Starting time
echo -e "\n`date` Start"

# Working directory
wd="$1"
cd "${wd}"

# Input
input_file="$2"
bed_file="$3"

# Output
output_file="$4"

# Masking
echo -e "\n`date` Start Masking"
bedtools maskfasta -fi ${input_file} -bed ${bed_file} -fo ${output_file}

echo -e "\n`date` Done with Masking"