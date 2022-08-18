#!/bin/bash
 
#SBATCH -A snic2022-22-143
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 72:00:00
#SBATCH -J 10X_gen

# Change these settings for different runs! Change SBATCH --array to the number
# of samples you have in shell

# Error report
set -euo pipefail

# Working directory
wd="$1"

# Set resources
ref_Cellranger="$2"

# Load required tools
echo -e "\n`date` Load modules"
module load bioinfo-tools
module load samtools/1.14
module load sratools/2.10.9
module load cellranger/6.1.2

# Assign names for arrays
cd "${wd}"
names_1=($(cat jobs1))
selected_sample_1=${names_1[$((SLURM_ARRAY_TASK_ID-1))]}
names_2=($(cat jobs2))
selected_sample_2=${names_2[$((SLURM_ARRAY_TASK_ID-1))]}
names_3=($(cat jobs3))
selected_sample_3=${names_3[$((SLURM_ARRAY_TASK_ID-1))]}
echo -e "Selected samples: ${selected_sample_1}, ${selected_sample_2}, ${selected_sample_3}"

links_1=($(cat ftp_links1.txt))
ftp_link_1=${links_1[$((SLURM_ARRAY_TASK_ID-1))]}
links_2=($(cat ftp_links2.txt))
ftp_link_2=${links_2[$((SLURM_ARRAY_TASK_ID-1))]}
links_3=($(cat ftp_links3.txt))
ftp_link_3=${links_3[$((SLURM_ARRAY_TASK_ID-1))]}

exp_names=($(cat exp_name.txt))
experiment=${exp_names[$((SLURM_ARRAY_TASK_ID-1))]}
echo -e "Selected experiment: ${experiment}"

sample_nrs_1=($(cat sample_nr1.txt))
sample_number_1=${sample_nrs_1[$((SLURM_ARRAY_TASK_ID-1))]}
sample_nrs_2=($(cat sample_nr2.txt))
sample_number_2=${sample_nrs_2[$((SLURM_ARRAY_TASK_ID-1))]}
sample_nrs_3=($(cat sample_nr3.txt))
sample_number_3=${sample_nrs_3[$((SLURM_ARRAY_TASK_ID-1))]}
echo -e "Selected sample numbers: ${sample_number_1}, ${sample_number_2}, ${sample_number_3}"

# Generate temporary folder structure
mkdir -p $SNIC_TMP/{processed,data}

# Start acquiring files
echo -e "\n`date` Downloading ${experiment} with prefetch"
cd "$SNIC_TMP/data/"
mkdir -p ${experiment}
cd ${experiment}
prefetch -p $selected_sample_1
prefetch -p $selected_sample_2
prefetch -p $selected_sample_3

# Convert sra to fastq in parallel
echo -e "\n`date` Start conversion of ${experiment} with xargs and fastq-dump"
ls */*.sra | xargs -P 3 -n 1 fastq-dump --gzip --split-files

# Remove folders with .sra files
rm -r $selected_sample_1
rm -r $selected_sample_2
rm -r $selected_sample_3

# Rename fastqs to Cellranger-compatible files
# For lane 2
mv "${selected_sample_1}_1.fastq.gz" "${experiment}_${sample_number_1}_R1_001.fastq.gz"
mv "${selected_sample_1}_2.fastq.gz" "${experiment}_${sample_number_1}_R2_001.fastq.gz"
mv "${selected_sample_1}_3.fastq.gz" "${experiment}_${sample_number_1}_I1_001.fastq.gz"

# For lane 2
mv "${selected_sample_2}_1.fastq.gz" "${experiment}_${sample_number_2}_R1_001.fastq.gz"
mv "${selected_sample_2}_2.fastq.gz" "${experiment}_${sample_number_2}_R2_001.fastq.gz"
mv "${selected_sample_2}_3.fastq.gz" "${experiment}_${sample_number_2}_I1_001.fastq.gz"

# For lane 2
mv "${selected_sample_3}_1.fastq.gz" "${experiment}_${sample_number_3}_R1_001.fastq.gz"
mv "${selected_sample_3}_2.fastq.gz" "${experiment}_${sample_number_3}_R2_001.fastq.gz"
mv "${selected_sample_3}_3.fastq.gz" "${experiment}_${sample_number_3}_I1_001.fastq.gz"

# cd to processed directory
cd "$SNIC_TMP/processed/"

# Start counting
cellranger count --id=${experiment} \
    --fastqs="$SNIC_TMP/data/${experiment}" \
    --transcriptome=${ref_Cellranger} \
    --no-bam \
    --localcores=16 \
    --localmem=96

# Copy output to wd
mkdir -p ${wd}/${experiment}
cp -r $SNIC_TMP/processed/${experiment}/outs ${wd}/${experiment}

# List all used modules in run
echo -e "\n`date` SessionInfo"
module list

echo -e "\n`date` Done with download"