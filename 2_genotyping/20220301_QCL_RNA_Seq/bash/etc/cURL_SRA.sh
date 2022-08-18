#!/bin/bash
 
#SBATCH -A snic2022-22-143
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 72:00:00
#SBATCH -J download_sra

# Error report
set -euo pipefail

# Starting time
echo -e "\n`date` Start"

# Working directory
wd="$1"
cd "${wd}"

# Download
echo -e "\n`date` Start download"
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/srr/SRR137/029/SRR13782529 -o SRR13782529_GSM5106143_02_E4_15_H09_Mus_musculus_RNA-Seq.sra
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/srr/SRR137/030/SRR13782530 -o SRR13782530_GSM5106143_02_E4_15_H09_Mus_musculus_RNA-Seq.sra
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/srr/SRR137/031/SRR13782531 -o SRR13782531_GSM5106143_02_E4_15_H09_Mus_musculus_RNA-Seq.sra
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/srr/SRR137/032/SRR13782532 -o SRR13782532_GSM5106144_03_E3_10_H10_Mus_musculus_RNA-Seq.sra
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/srr/SRR137/033/SRR13782533 -o SRR13782533_GSM5106144_03_E3_10_H10_Mus_musculus_RNA-Seq.sra
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/srr/SRR137/034/SRR13782534 -o SRR13782534_GSM5106144_03_E3_10_H10_Mus_musculus_RNA-Seq.sra

echo -e "\n`date` Done with download"