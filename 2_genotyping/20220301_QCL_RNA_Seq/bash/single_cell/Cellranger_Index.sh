#!/bin/bash
 
#SBATCH -A snic2022-22-143
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 24:00:00
#SBATCH -J cellranger_index

# Error report
set -euo pipefail

# Starting time
echo -e "\n`date` Start"

# Working directory
wd="$1"
cd "${wd}"

# Set resources
ref="$2"
gtf="$3"

# Set output
genomename="$4"

# Load required tools
module load bioinfo-tools
module load cellranger/6.1.2

# Indexing
echo -e "\n`date` Start indexing"
cellranger mkref --genome=${genomename} \
    --fasta=${ref} \
    --genes=${gtf}

echo -e "\n`date` Done with indexing (Cellranger)"