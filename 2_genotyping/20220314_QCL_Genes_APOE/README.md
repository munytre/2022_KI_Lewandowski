# 20220314_QCL_Genes_APOE
## General information
This repository contains scripts used to analyse APOE allele variants in  datasets involving ALS and AD.

Scripts written in R are located in the R folder. <br>
Resources folder contains metadata and other resources needed to perform the R scripts.
## R
### FilteringROSMAP.R
 - Script used to stochastically  filter ROSMAP samples to act as validation samples
### APOE_SNPs.R
 - Script to assign genotype to SNPs and assess genotype quality
### mpileup.R
- Script to assess read depth for SNPs
### DESeq2_Human_ALS.R
- Script containing DESeq2 analysis and PCAs for PRJNA512012
### PRJNA512012_repro.R
- Script to assess reproducibility for analysed samples in PRJNA512012
### ROSMAP1_repro.R
- Script to assess reproducibility for analysed samples in ROSMAP batch 1
### Survival_PRJNA512012.R
- Script to perform Cox Proportional Hazard Models on genotyped samples of PRJNA512012
## Contact
For questions don't hesitate to send me a message at:
q.c.lin@students.uu.nl