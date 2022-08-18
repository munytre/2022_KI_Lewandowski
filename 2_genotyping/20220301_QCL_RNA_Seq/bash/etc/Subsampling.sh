#!/bin/bash
 
#SBATCH -A snic2022-22-143
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 12:00:00
#SBATCH -J subsampling_parameters

# Error report
set -euo pipefail

# Starting time
echo -e "\n`date` Start"

# Working directory
wd="$1"
cd "${wd}"

# Load required tools (Detection)
module load bioinfo-tools
module load seqtk/1.2-r101

# Unzipping
echo -e "\n`date` Unzipping"
gzip -dc SRR8375274/SRR8375274_1.fastq.gz > $SNIC_TMP/SRR8375274_1.fastq
gzip -dc SRR8375274/SRR8375274_2.fastq.gz > $SNIC_TMP/SRR8375274_2.fastq

# Subsampling with seqtk
echo -e "\n`date` Subsampling"
seqtk sample -s747 $SNIC_TMP/SRR8375274_1.fastq 1000000 > SRR8375274_1_sub.fastq
seqtk sample -s747 $SNIC_TMP/SRR8375274_1.fastq 1000000 > SRR8375274_2_sub.fastq

# Deleting raw.fastq
echo -e "\n`date` Deleting unzipped files"
rm $SNIC_TMP/SRR8375274_1.fastq
rm $SNIC_TMP/SRR8375274_2.fastq

# Gzipping subsamples
echo -e "\n`date` Gzipping subsamples"
gzip SRR8375274_1_sub.fastq
gzip SRR8375274_2_sub.fastq

echo -e "\n`date` Done with subsampling"