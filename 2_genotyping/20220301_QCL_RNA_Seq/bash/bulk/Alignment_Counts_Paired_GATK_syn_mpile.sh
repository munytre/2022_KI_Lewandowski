#!/bin/bash
 
#SBATCH -A snic2022-22-143
#SBATCH -p core
#SBATCH -n 12
#SBATCH -t 96:00:00
#SBATCH -J GATK_Counts

# Change these settings for different runs! Change SBATCH --array to the number
# of samples you have in shell

# Error report
set -euo pipefail

# Working directory
wd="$1"

# Set resources
ref_STAR="$2"
ref_gtf="$3"
ref_gen="$4"

# Synapse login
synapse_user="$5"
synapse_api="$6"

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
module load synapseclient/2.3.1

# Assign names for arrays
cd "${wd}"
names=($(cat jobs))
selected_sample=${names[$((SLURM_ARRAY_TASK_ID-1))]}
echo -e "Selected sample: ${selected_sample}"
links1=($(cat syn_links1.txt))
syn_link1=${links1[$((SLURM_ARRAY_TASK_ID-1))]}
links2=($(cat syn_links2.txt))
syn_link2=${links2[$((SLURM_ARRAY_TASK_ID-1))]}

# Generate temporary folder structure
mkdir -p $SNIC_TMP/{processed,data}

# Download data with synapse
echo -e "\n`date` Downloading ${selected_sample} with Synapse"
cd "$SNIC_TMP/data/"
mkdir -p ${selected_sample}
cd ${selected_sample}
synapse login --rememberMe -u ${synapse_user} -p ${synapse_api}
synapse get ${syn_link1}
synapse get ${syn_link2}
mv ${selected_sample}_*_R1_001.fastq.gz ${selected_sample}_1.fastq.gz
mv ${selected_sample}_*_R2_001.fastq.gz ${selected_sample}_2.fastq.gz

# Create processed dir with underlying dirs for sample
cd "$SNIC_TMP/processed/"
mkdir -p "${selected_sample}"

# Run TrimGalore on the fastq
echo -e "\n`date` Filtering and trimming ${selected_sample} ..."
trim_galore --cores 8 \
	--trim-n \
	--fastqc \
	--gzip \
	--paired \
	"$SNIC_TMP/data/${selected_sample}/${selected_sample}_1.fastq.gz" \
	"$SNIC_TMP/data/${selected_sample}/${selected_sample}_2.fastq.gz" \
	-o "$SNIC_TMP/processed/${selected_sample}/trimgalore"

# Copy QC files to wd
echo -e "\n`date` Copying QC files for ${selected_sample}"
mkdir -p ${wd}/processed/${selected_sample}/trimgalore/
cp $SNIC_TMP/processed/${selected_sample}/trimgalore/*.html ${wd}/processed/${selected_sample}/trimgalore/

### Alignment and counts ###
#Run STAR
echo -e "\n`date` Mapping ${selected_sample} with STAR"
STAR --genomeDir ${ref_STAR} \
    --runThreadN 12 \
    --twopassMode Basic \
    --readFilesIn \
    "$SNIC_TMP/processed/${selected_sample}/trimgalore/${selected_sample}_1_val_1.fq.gz" \
    "$SNIC_TMP/processed/${selected_sample}/trimgalore/${selected_sample}_2_val_2.fq.gz" \
    --readFilesCommand zcat \
    --outFileNamePrefix "$SNIC_TMP/processed/${selected_sample}/STAR/${selected_sample}" \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes All \
    --limitBAMsortRAM 70000000000

# Index by samtools
echo -e "\n`date` Indexing ${selected_sample} with samtools"
samtools index -@ 12 "$SNIC_TMP/processed/${selected_sample}/STAR/${selected_sample}Aligned.sortedByCoord.out.bam"

# Copy log files to wd
echo -e "\n`date` Copying log files for ${selected_sample}"
mkdir -p ${wd}/processed/${selected_sample}/STAR
cp -R $SNIC_TMP/processed/${selected_sample}/STAR/*.out ${wd}/processed/${selected_sample}/STAR

### GATK - RNAseq version ###
cd $SNIC_TMP/processed/${selected_sample}/STAR

# Add read groups for HaplotypeCaller
echo -e "\n`date` Adding read groups for ${selected_sample} with AddOrReplaceReadGroups"
java -Xmx70G -XX:ParallelGCThreads=12 -jar $PICARD_ROOT/picard.jar AddOrReplaceReadGroups I=${selected_sample}Aligned.sortedByCoord.out.bam \
    O=${selected_sample}RGs.bam \
    RGID=1 \
    RGLB=lib1 \
    RGPL=illumina \
    RGPU=unit1 \
    RGSM=20

# Index after adding read groups
echo -e "\n`date` Re-index ${selected_sample} with samtools"
samtools index -@ 12 ${selected_sample}RGs.bam

# Remove duplicated reads
echo -e "\n`date` Removing duplicated reads from ${selected_sample} with MarkDuplicates"
java -Xmx70G -XX:ParallelGCThreads=12 -jar $PICARD_ROOT/picard.jar MarkDuplicates I=${selected_sample}RGs.bam \
    O=${selected_sample}RemovedDups.bam \
    M=${selected_sample}RemovedDups.txt \
    REMOVE_DUPLICATES=TRUE

# Sort BAM after MarkDuplicates
echo -e "\n`date` Sorting ${selected_sample} after removal of duplicated reads with SortSam"
java -Xmx70G -XX:ParallelGCThreads=12 -jar $PICARD_ROOT/picard.jar SortSam I=${selected_sample}RemovedDups.bam \
    O=${selected_sample}RemovedDups_sorted.bam \
    SORT_ORDER=coordinate \
    CREATE_INDEX=TRUE

# Acquiring coverage with samtools mpileup
echo -e "\n`date` Acquiring coverage for APOE in ${selected_sample} with samtools mpileup"
samtools mpileup ${selected_sample}RemovedDups_sorted.bam \
    -f ${ref_gen} \
    -r 19:44905791-44909393 \
    -o ${selected_sample}_mpileup.txt

# Copy relevant files from GATK pipeline
echo -e "\n`date` Copying pileup output for ${selected_sample} to wd"
mkdir -p ${wd}/processed/${selected_sample}/samtools
cp $SNIC_TMP/processed/${selected_sample}/STAR/${selected_sample}_mpileup.txt ${wd}/processed/${selected_sample}/samtools

# List all used modules in run
echo -e "\n`date` SessionInfo"
module list

echo -e "\n`date` Finished!"