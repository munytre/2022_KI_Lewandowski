#!/bin/bash
 
#SBATCH -A snic2022-22-143
#SBATCH -p core
#SBATCH -n 12
#SBATCH --mem 64G
#SBATCH -t 12:00:00
#SBATCH -J star_index

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
readlength="$4"

# Set output
genomedir="$5"

# Load required tools
module load bioinfo-tools
module load star/2.7.9a

# Indexing
echo -e "\n`date` Start indexing"
STAR \
    --runMode genomeGenerate \
    --runThreadN 12 \
    --genomeDir $genomedir \
    --genomeFastaFiles $ref \
    --sjdbOverhang $readlength \
    --sjdbGTFfile $gtf

echo -e "\n`date` Done with indexing (STAR)"