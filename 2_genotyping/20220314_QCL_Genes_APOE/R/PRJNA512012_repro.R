###--------------------------------------------------------------------------###
### Project: ApoE in ALS                                                       -
### Purpose: Reproduction analysis (PRJNA512012)                               -
### Authors: QC Lin                                                            -
### Last Edit: 2022-08-15, QC Lin                                              -
### Contact: q.c.lin@students.uu.nl or sebastian.lewandowski@ki.se             -
###--------------------------------------------------------------------------###

### Dependencies and general settings ------------------------------------------
library(tidyverse)

### Count files ----------------------------------------------------------------
# Load in data
SRR8375392 <- read.delim("20220301_QCL_RNA_Seq/PRJNA512012/processed/SRR8375392/htseq_count/SRR8375392.count",
                               header=FALSE)
SRR8375325 <- read.delim("20220301_QCL_RNA_Seq/PRJNA512012/processed/SRR8375325/htseq_count/SRR8375325.count",
                               header=FALSE)

SRR8375392_repro <- read.delim("20220301_QCL_RNA_Seq/PRJNA512012_repro/processed/SRR8375392/htseq_count/SRR8375392.count",
                               header=FALSE)
SRR8375325_repro <- read.delim("20220301_QCL_RNA_Seq/PRJNA512012_repro/processed/SRR8375325/htseq_count/SRR8375325.count",
                               header=FALSE)

# See if files are equal
setequal(SRR8375392,SRR8375392_repro)
setequal(SRR8375325,SRR8375325_repro)
setequal(SRR8375325,SRR8375392)

### VCF files ------------------------------------------------------------------
# Load in data
SRR8375392_vcf_repro <- read.table("20220301_QCL_RNA_Seq/PRJNA512012_repro/processed/SRR8375392/GATK/SRR8375392_19_44500000_45000000_filtered.vcf.gz",
                                   quote="\"")
SRR8375325_vcf_repro <- read.table("20220301_QCL_RNA_Seq/PRJNA512012_repro/processed/SRR8375325/GATK/SRR8375325_19_44500000_45000000_filtered.vcf.gz",
                                   quote="\"")

SRR8375392_vcf <- read.table("20220301_QCL_RNA_Seq/PRJNA512012/processed/SRR8375392/GATK/SRR8375392_19_44500000_45000000_filtered.vcf.gz",
                                   quote="\"")
SRR8375325_vcf <- read.table("20220301_QCL_RNA_Seq/PRJNA512012/processed/SRR8375325/GATK/SRR8375325_19_44500000_45000000_filtered.vcf.gz",
                                   quote="\"")

# See if files are equal
setequal(SRR8375392_vcf,SRR8375392_vcf_repro)
setequal(SRR8375325_vcf,SRR8375325_vcf_repro)
setequal(SRR8375325_vcf,SRR8375392_vcf)