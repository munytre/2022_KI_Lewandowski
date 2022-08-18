###--------------------------------------------------------------------------###
### Project: ApoE in ALS                                                       -
### Purpose: Reproduction analysis (ROSMAP1)                                   -
### Authors: QC Lin                                                            -
### Last Edit: 2022-08-15, QC Lin                                              -
### Contact: q.c.lin@students.uu.nl or sebastian.lewandowski@ki.se             -
###--------------------------------------------------------------------------###

### Dependencies and general settings ------------------------------------------
library(tidyverse)

### Count files ----------------------------------------------------------------
# Load in data
`583_120522` <- read.delim("20220301_QCL_RNA_Seq/ROSMAP_batch1/processed/583_120522/htseq_count/583_120522.count",
                               header=FALSE)
`331_120501` <- read.delim("20220301_QCL_RNA_Seq/ROSMAP_batch1/processed/331_120501/htseq_count/331_120501.count",
                               header=FALSE)

`583_120522_repro` <- read.delim("20220301_QCL_RNA_Seq/ROSMAP_batch1_repro/processed/583_120522/htseq_count/583_120522.count",
                               header=FALSE)
`331_120501_repro` <- read.delim("20220301_QCL_RNA_Seq/ROSMAP_batch1_repro/processed/331_120501/htseq_count/331_120501.count",
                               header=FALSE)

# See if files are equal
setequal(`583_120522`,`583_120522_repro`)
setequal(`331_120501`,`331_120501_repro`)
setequal(`583_120522`,`331_120501`)

### VCF files ------------------------------------------------------------------
# Load in data
`583_120522_vcf_repro` <- read.table("20220301_QCL_RNA_Seq/ROSMAP_batch1_repro/processed/583_120522/GATK/583_120522_19_44500000_45000000_filtered.vcf.gz",
                                   quote="\"")
`331_120501_vcf_repro` <- read.table("20220301_QCL_RNA_Seq/ROSMAP_batch1_repro/processed/331_120501/GATK/331_120501_19_44500000_45000000_filtered.vcf.gz",
                                   quote="\"")

`583_120522_vcf` <- read.table("20220301_QCL_RNA_Seq/ROSMAP_batch1/processed/583_120522/GATK/583_120522_19_44500000_45000000_filtered.vcf.gz",
                                   quote="\"")
`331_120501_vcf` <- read.table("20220301_QCL_RNA_Seq/ROSMAP_batch1/processed/331_120501/GATK/331_120501_19_44500000_45000000_filtered.vcf.gz",
                                   quote="\"")

# See if files are equal
setequal(`583_120522_vcf`,`583_120522_vcf_repro`)
setequal(`331_120501_vcf`,`331_120501_vcf_repro`)
setequal(`583_120522_vcf`,`331_120501_vcf`)
