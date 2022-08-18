###--------------------------------------------------------------------------###
### Project: ApoE in ALS                                                       -
### Purpose: ROSMAP filtering script                                           -
### Authors: QC Lin                                                            -
### Last Edit: 2022-08-15, QC Lin                                              -
### Contact: q.c.lin@students.uu.nl or sebastian.lewandowski@ki.se             -
###--------------------------------------------------------------------------###

### Dependencies and general settings ------------------------------------------
library(tidyverse)
set.seed(747)

### Load data ------------------------------------------------------------------
ROSMAP_clinical <- read.csv("20220314_QCL_Genes_APOE/ROSMAP/ROSMAP_clinical.csv")
ROSMAP_biospecimen_metadata <- read.csv("20220314_QCL_Genes_APOE/ROSMAP/ROSMAP_biospecimen_metadata.csv")
ROSMAP_assay_rnaSeq_metadata <- read.csv("20220314_QCL_Genes_APOE/ROSMAP/ROSMAP_assay_rnaSeq_metadata.csv")

ROSMAP <- left_join(ROSMAP_assay_rnaSeq_metadata, ROSMAP_biospecimen_metadata, by = "specimenID")
ROSMAPClin <- left_join(ROSMAP, ROSMAP_clinical, by = "individualID")
ROSMAPClin_1 <- filter(ROSMAPClin,
                       notes == "data contribution batch 1")
ROSMAPClin_1 <- filter(ROSMAPClin_1,
                       assay == "rnaSeq")
ROSMAPClin_1_APOE <- filter(ROSMAPClin_1,
                            !is.na(apoe_genotype))

# Subset samples for APOE genotype
ROSMAP_APOE22 <- filter(ROSMAPClin_1_APOE,
                        apoe_genotype == "22")
ROSMAP_APOE23 <- filter(ROSMAPClin_1_APOE,
                        apoe_genotype == "23")
ROSMAP_APOE24 <- filter(ROSMAPClin_1_APOE,
                        apoe_genotype == "24")
ROSMAP_APOE33 <- filter(ROSMAPClin_1_APOE,
                        apoe_genotype == "33")
ROSMAP_APOE34 <- filter(ROSMAPClin_1_APOE,
                        apoe_genotype == "34")
ROSMAP_APOE44 <- filter(ROSMAPClin_1_APOE,
                        apoe_genotype == "44")

# Sampling ---------------------------------------------------------------------
ROSMAP_APOE22_R <- sample(c(ROSMAP_APOE22$specimenID), 5, replace=FALSE)
ROSMAP_APOE23_R <- sample(c(ROSMAP_APOE23$specimenID), 30, replace=FALSE)
ROSMAP_APOE24_R <- sample(c(ROSMAP_APOE24$specimenID), 5, replace=FALSE)
ROSMAP_APOE33_R <- sample(c(ROSMAP_APOE33$specimenID), 40, replace=FALSE)
ROSMAP_APOE34_R <- sample(c(ROSMAP_APOE34$specimenID), 35, replace=FALSE)
ROSMAP_APOE44_R <- sample(c(ROSMAP_APOE44$specimenID), 5, replace=FALSE)

ROSMAP_selected_samples <- c(ROSMAP_APOE22_R,
                             ROSMAP_APOE23_R,
                             ROSMAP_APOE24_R,
                             ROSMAP_APOE33_R,
                             ROSMAP_APOE34_R,
                             ROSMAP_APOE44_R)

write(sort(ROSMAP_selected_samples),
      "~/Desktop/20220314_QCL_Genes_APOE/Final/120_jobs.txt",
      sep = "\n")

filtered_samples <- filter(ROSMAPClin_1_APOE, specimenID %in% ROSMAP_selected_samples) %>%
  dplyr::select(specimenID, apoe_genotype)
write_delim(filtered_samples,
      "~/Desktop/20220314_QCL_Genes_APOE/Final/120_samples_genotype.txt")
