#!/bin/bash
 
#SBATCH -A snic2022-22-143
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 96:00:00
#SBATCH -J GATK_Filter

# Change these settings for different runs! Change SBATCH --array to the number
# of samples you have in shell

# Error report
set -euo pipefail

# Working directory
wd="$1"

ref_gen="$2"

# Load required tools
echo -e "\n`date` Load modules"
module load bioinfo-tools
module load TrimGalore/0.6.1
module load star/2.7.9a
module load samtools/1.14
module load htseq/0.12.4
module load GATK/4.2.0.0
module load picard/2.23.4
module load R/4.1.1
module load R_packages/4.1.1

# Assign names for arrays
cd "${wd}"
names=($(cat jobs))
selected_sample=${names[$((SLURM_ARRAY_TASK_ID-1))]}
echo -e "Selected sample: ${selected_sample}"

cd ${wd}/processed/${selected_sample}/GATK
gatk VariantFiltration \
    -R ${ref_gen} \
    -V ${selected_sample}_19_44500000_45000000.vcf.gz \
    -O ${selected_sample}_19_44500000_45000000_filtered.vcf.gz \
    --filter-name "FS" \
	--filter "FS > 30.0" \
	--filter-name "QD" \
	--filter "QD < 2.0"

# List all used modules in run
echo -e "\n`date` SessionInfo"
module list

echo -e "\n`date` Finished!"