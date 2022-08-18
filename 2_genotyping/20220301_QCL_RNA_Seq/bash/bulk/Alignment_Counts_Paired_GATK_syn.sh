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
names=($(cat jobs_nounderscore))
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
mv ${selected_sample}.r1.fastq.gz ${selected_sample}_1.fastq.gz
mv ${selected_sample}.r2.fastq.gz ${selected_sample}_2.fastq.gz

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

# Counting by htseq-count
echo -e "\n`date` Counting raw expression values of ${selected_sample} with htseq-count"
mkdir -p ${wd}/processed/${selected_sample}/htseq_count/
htseq-count -n 12 \
    --order pos \
    --stranded reverse \
    --counts_output "${wd}/processed/${selected_sample}/htseq_count/${selected_sample}.count" \
    "$SNIC_TMP/processed/${selected_sample}/STAR/${selected_sample}Aligned.sortedByCoord.out.bam" \
    "${ref_gtf}"

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

# Split reads on junctions (as described in best practices GATK - RNAseq)
echo -e "\n`date` Split reads on junctions for ${selected_sample} with SplitNCigarReads"
gatk --java-options "-Xmx70G -XX:ParallelGCThreads=12" SplitNCigarReads -R ${ref_gen} \
    -I ${selected_sample}RemovedDups_sorted.bam \
    -O ${selected_sample}RemovedDups_sorted_split.bam

# Generate table for QC Score recalculation
echo -e "\n`date` Generate table for QC Score recalculation for ${selected_sample} with BaseRecalibrator"
gatk --java-options "-Xmx70G -XX:ParallelGCThreads=12" BaseRecalibrator -I ${selected_sample}RemovedDups_sorted_split.bam \
    -R ${ref_gen} \
    --known-sites /home/munytre/RESOURCES/Homo_sapiens.GRCh38/VCF/dbsnp_146_non_chr.hg38.vcf.gz \
    -O ${selected_sample}_BQSR.table

# Recalculate QC Scores for all reads with BQSR
echo -e "\n`date` Recalculating QC Scored for ${selected_sample} with ApplyBQSR"
gatk --java-options "-Xmx70G -XX:ParallelGCThreads=12" ApplyBQSR -R ${ref_gen} \
    -I ${selected_sample}RemovedDups_sorted_split.bam \
    --bqsr-recal-file ${selected_sample}_BQSR.table \
    -O ${selected_sample}BQSR.bam

echo -e "\n`date` Generate 2nd table for QC Score recalculation for ${selected_sample} with BaseRecalibrator"
gatk --java-options "-Xmx70G -XX:ParallelGCThreads=12" BaseRecalibrator -I ${selected_sample}BQSR.bam \
    -R ${ref_gen} \
    --known-sites /home/munytre/RESOURCES/Homo_sapiens.GRCh38/VCF/dbsnp_146_non_chr.hg38.vcf.gz \
    -O ${selected_sample}_BQSR2.table

# Generate BQSR QC plots
echo -e "\n`date` Generate BQSR plots for ${selected_sample} with AnalyzeCovariates"
gatk --java-options "-Xmx70G -XX:ParallelGCThreads=12" AnalyzeCovariates -before ${selected_sample}_BQSR.table \
    -after ${selected_sample}_BQSR2.table \
    -plots ${selected_sample}.pdf

# Generate VCFs
echo -e "\n`date` Generate VCF for ${selected_sample} with HaplotypeCaller"
gatk --java-options "-Xmx70G -XX:ParallelGCThreads=12" HaplotypeCaller -R ${ref_gen} \
    -I ${selected_sample}BQSR.bam \
    -O ${selected_sample}_19_44500000_45000000.vcf.gz \
    --native-pair-hmm-threads 12 \
    --dbsnp /home/munytre/RESOURCES/Homo_sapiens.GRCh38/VCF/dbsnp_146_non_chr.hg38.vcf.gz \
    -L 19:44500000-45000000

# Filtering VCF (According to best practices GATK for RNAseq)
echo -e "\n`date` Filtering VCF for ${selected_sample} with VariantFiltration"
gatk VariantFiltration \
    -R ${ref_gen} \
    -V ${selected_sample}_19_44500000_45000000.vcf.gz \
    -O ${selected_sample}_19_44500000_45000000_filtered.vcf.gz \
    --filter-name "FS" \
	--filter "FS > 30.0" \
	--filter-name "QD" \
	--filter "QD < 2.0"

# Copy relevant files from GATK pipeline
echo -e "\n`date` Copying VCF and BQRS plots for ${selected_sample} to wd"
mkdir -p ${wd}/processed/${selected_sample}/GATK
cp $SNIC_TMP/processed/${selected_sample}/STAR/*.txt ${wd}/processed/${selected_sample}/GATK
cp $SNIC_TMP/processed/${selected_sample}/STAR/*.pdf ${wd}/processed/${selected_sample}/GATK
cp $SNIC_TMP/processed/${selected_sample}/STAR/*.vcf* ${wd}/processed/${selected_sample}/GATK

# List all used modules in run
echo -e "\n`date` SessionInfo"
module list

echo -e "\n`date` Finished!"