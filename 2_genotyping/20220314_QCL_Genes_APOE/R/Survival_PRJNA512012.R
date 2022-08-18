###--------------------------------------------------------------------------###
### Project: ApoE in ALS                                                       -
### Purpose: Survival analysis of genotyped APOE in PRJNA512012                -
### Authors: QC Lin                                                            -
### Last Edit: 2022-08-15, QC Lin                                              -
### Contact: q.c.lin@students.uu.nl or sebastian.lewandowski@ki.se             -
###--------------------------------------------------------------------------###

### Dependencies and general settings ------------------------------------------
library(survival)
library(survminer)
library(GEOquery)
library(tidyverse)
library(GGally)
library(broom)
library(ggbeeswarm)
library(ggpubr)

### Load data ------------------------------------------------------------------
# Load in metadata
gse <- getGEO(filename = "20220301_QCL_GSE124439/metadata/GSE124439_series_matrix.txt")
metadata <- gse@phenoData@data
colnames(metadata) <- sub(".", "_", colnames(metadata), fixed=TRUE)
colnames(metadata) <- sub(":", "_", colnames(metadata), fixed=TRUE)
colnames(metadata) <- sub(" ", "_", colnames(metadata), fixed=TRUE)

metadata_excel1 <- read.delim("20220301_QCL_GSE124439/metadata/GSE124439_sheet1.txt",
                              skip = 1,
                              nrows = 176)
colnames(metadata_excel1)[1] <- "title"

metadata_excel2 <- read.delim("20220301_QCL_GSE124439/metadata/GSE124439_sheet2.txt",
                              skip = 1,
                              nrows = 77)
colnames(metadata_excel2)[1] <- "subject_id_ch1"
metadata_excel2 <- metadata_excel2[,c(1,5:8)]

metadata <- left_join(metadata, metadata_excel1, by = "title")
metadata <- left_join(metadata, metadata_excel2, by = "subject_id_ch1")
metadata <- metadata[,-c(10:43)]

SRR_sample <- read.delim("20220301_QCL_GSE124439/metadata/PRJNA512012_tsv.txt")
SRR_sample <- SRR_sample[,c(2,4)]
colnames(SRR_sample) <- c("SRR", "title")

metadata_SRR_included <- left_join(SRR_sample,
                                   metadata,
                                   by = "title")

# Load in genotype data
genotype_PRJNA512012 <- read.csv("20220314_QCL_Genes_APOE/Final/VCF_output_PRJNA512012.txt",
                                 sep="")

# Join metadata and genotype data
pre_metadata_SRR_included <- left_join(metadata_SRR_included,
                                       genotype_PRJNA512012,
                                       by = c("SRR" = "sample"))
# Change names for readability
metadata_SRR_included <- pre_metadata_SRR_included %>%
  dplyr::select(SRR,
                source_name_ch1,
                Tissue,
                Gender,
                sample_group_ch1,
                subject_id_ch1,
                C9orf72.repeat.expansion.,
                ALS.NMF.subtype..,
                Age.of.Onset,
                Site.of.Onset,
                Disease.Duration..months.,
                Age.of.Death,
                APOE) %>%
  `colnames<-`(c("SRR",
                 "Source",
                 "Tissue",
                 "Gender",
                 "Subject_Group",
                 "Individual",
                 "C9orf72_Expansion",
                 "Subtype",
                 "Age_Onset",
                 "Site_Onset",
                 "Disease_Duration",
                 "Age_Death",
                 "APOE"))

# Change column to correct datatype
metadata_SRR_included$APOE <- as.factor(metadata_SRR_included$APOE)
metadata_SRR_included <- metadata_SRR_included %>%
  mutate(APOE4_State = if_else(APOE %in% c(34,44),
                               "APOE4pos",
                               "APOE4neg"))
metadata_SRR_included$Source <- as.factor(metadata_SRR_included$Source)
metadata_SRR_included$Tissue <- as.factor(metadata_SRR_included$Tissue)
metadata_SRR_included$Gender <- as.factor(metadata_SRR_included$Gender)
metadata_SRR_included$Subject_Group <- as.factor(metadata_SRR_included$Subject_Group)
metadata_SRR_included$Individual <- as.factor(metadata_SRR_included$Individual)
metadata_SRR_included$C9orf72_Expansion <- as.factor(metadata_SRR_included$C9orf72_Expansion)
metadata_SRR_included$Subtype <- as.factor(metadata_SRR_included$Subtype)
metadata_SRR_included$Site_Onset <- as.factor(metadata_SRR_included$Site_Onset)
metadata_SRR_included$APOE4_State <- as.factor(metadata_SRR_included$APOE4_State)

# Remove redundant data
rm(gse,
   genotype_PRJNA512012,
   metadata,metadata_excel1,
   metadata_excel2,
   SRR_sample,
   pre_metadata_SRR_included)

# Update rownames and filter for ALS-only samples
rownames(metadata_SRR_included) <- metadata_SRR_included$SRR
metadata_SRR_included_ALS <- metadata_SRR_included[(rownames(metadata_SRR_included) %in% metadata_SRR_included$SRR[metadata_SRR_included$Subject_Group == "ALS Spectrum MND"]),]

# Load in survival data
# NYGC
survival_data <- metadata_SRR_included_ALS %>%
  select(Gender,
         Subject_Group,
         Individual,
         Age_Onset,
         Site_Onset,
         Disease_Duration,
         Age_Death,
         APOE,
         APOE4_State,
         Subtype) %>%
  unique()

survival_data <- survival_data[!is.na(survival_data$Disease_Duration),]
survival_data$APOE <- factor(survival_data$APOE, levels = c(33,23,24,34,44))
survival_data$Status <- !is.na(survival_data$Age_Death)
survival_data$Status <- as.integer(as.logical(survival_data$Status))
survival_data <- survival_data %>%
  mutate(onset_bulbar = Site_Onset %in% "Bulbar")

# surv_obj <- coxph(formula = Surv(time = Disease_Duration, event = Status) ~ onset_bulbar + Age_Onset + APOE, data = survival_data)
# ggforest(surv_obj,
#          data = survival_data)
# cox.zph(surv_obj)

# UTRECHT
als_beads_utrecht <- readRDS(file = "20220131_QCL_Survival_Analysis/generated_data_R/apoeALS.rds")
survival_data_utrecht <- als_beads_utrecht %>%
  select(gender,
         class,
         sample_name,
         age_at_onset,
         onset_bulbar,
         survival_months,
         apoe,
         died)
survival_data_utrecht$Gender <- if_else(survival_data_utrecht$gender == "M",
                                        "Male",
                                        "Female")
survival_data_utrecht$APOE <- if_else(survival_data_utrecht$apoe == "ALS Apo-e3e3",
                                      33,
                                      34,
                                      missing = NULL)
survival_data_utrecht <- survival_data_utrecht %>%
  mutate(APOE4_State = if_else(APOE == c(34,44),
                               "APOE4pos",
                               "APOE4neg",
                               missing = NULL))
survival_data_utrecht$Status <- as.integer(as.logical(survival_data_utrecht$died))

# Join NYGC and UTRECHT data
NYCG <- survival_data %>%
  select(Individual,
         Gender,
         Disease_Duration,
         Status,
         Onset_Bulbar = onset_bulbar,
         Age_Onset,
         APOE,
         APOE4_State)
UTRECHT <- survival_data_utrecht %>%
  select(Individual = sample_name,
         Gender,
         Disease_Duration = survival_months,
         Status,
         Onset_Bulbar = onset_bulbar,
         Age_Onset = age_at_onset,
         APOE,
         APOE4_State)

UTRECHT$Gender <- as.factor(UTRECHT$Gender)
UTRECHT$Disease_Duration <- as.integer(UTRECHT$Disease_Duration)
UTRECHT$Age_Onset <- as.integer(UTRECHT$Age_Onset)

UTRECHT$APOE <- as.factor(UTRECHT$APOE) %>%
  factor(levels = c(33,23,24,34,44))
UTRECHT$APOE4_State <- as.factor(UTRECHT$APOE4_State)


joined_cohorts <- list(New_York = NYCG,
                       Utrecht = UTRECHT) %>%
  bind_rows(.id = "cohort")

### Cox Proportional Hazard Models ---------------------------------------------
# Joined (APOE4_State)
surv_obj <- coxph(formula = Surv(time = Disease_Duration, event = Status) ~ Onset_Bulbar + Age_Onset + APOE4_State + Gender + strata(cohort), data = joined_cohorts)
ggcoef_model(surv_obj,
             add_reference_rows = T,
             exponentiate = T,
             variable_labels = c(
               Onset_Bulbar = "Bulbar Onset",
               Age_Onset = "Age of Onset",
               APOE4_State = "APOE4 Status"
             )) +
  labs(x = "HR") +
#  scale_x_continuous(limits = c(0,3.5)) +
  ggtitle("Bulbar Onset + Age of Onset + APOE4 Status + strata(Cohort)")
cox.zph(surv_obj)

# Joined (APOE)
surv_obj <- coxph(formula = Surv(time = Disease_Duration, event = Status) ~ Onset_Bulbar + Age_Onset + APOE + Gender + strata(cohort), data = joined_cohorts)
ggcoef_model(surv_obj,
             add_reference_rows = T,
             exponentiate = T,
             variable_labels = c(
               Onset_Bulbar = "Bulbar Onset",
               Age_Onset = "Age of Onset",
               APOE = "APOE Genotype"
             )) +
  labs(x = "HR") +
#  scale_x_continuous(limits = c(0,3.5)) +
  ggtitle("Bulbar Onset + Age of Onset + APOE Genotype + strata(Cohort)")
cox.zph(surv_obj)

# Joined (APOE:AoO)
surv_obj <- coxph(formula = Surv(time = Disease_Duration, event = Status) ~ Age_Onset*APOE + Onset_Bulbar, data = UTRECHT)
ggcoef_model(surv_obj,
             add_reference_rows = T,
             exponentiate = T,
             variable_labels = c(
               Onset_Bulbar = "Bulbar Onset",
               Age_Onset = "Age of Onset",
               APOE = "APOE Genotype"
             )) +
  labs(x = "HR") +
  #  scale_x_continuous(limits = c(0,3.5)) +
  ggtitle("Age of Onset:APOE Genotype + Bulbar Onset")
cox.zph(surv_obj)

# Utrecht (APOE4_State)
surv_obj <- coxph(formula = Surv(time = Disease_Duration, event = Status) ~ Onset_Bulbar + Age_Onset + APOE4_State + Gender, data = UTRECHT)
ggcoef_model(surv_obj,
             add_reference_rows = T,
             exponentiate = T,
             variable_labels = c(
               Onset_Bulbar = "Bulbar Onset",
               Age_Onset = "Age of Onset",
               APOE4_State = "APOE4 Status"
             )) +
  labs(x = "HR") +
#  scale_x_continuous(limits = c(0,3.5)) +
  ggtitle("Bulbar Onset + Age of Onset + APOE4 Status")
cox.zph(surv_obj)

# Utrecht (APOE)
surv_obj <- coxph(formula = Surv(time = Disease_Duration, event = Status) ~ Onset_Bulbar + Age_Onset + APOE + Gender, data = UTRECHT)
ggcoef_model(surv_obj,
             add_reference_rows = T,
             exponentiate = T,
             variable_labels = c(
               Onset_Bulbar = "Bulbar Onset",
               Age_Onset = "Age of Onset",
               APOE = "APOE Genotype"
             )) +
  labs(x = "HR") +
#  scale_x_continuous(limits = c(0,3.5)) +
  ggtitle("Bulbar Onset + Age of Onset + APOE Genotype")
cox.zph(surv_obj)

# NYCG (APOE4_State)
surv_obj <- coxph(formula = Surv(time = Disease_Duration, event = Status) ~ Onset_Bulbar + Age_Onset + APOE4_State, data = NYCG)
ggcoef_model(surv_obj,
             add_reference_rows = T,
             exponentiate = T,
             variable_labels = c(
               Onset_Bulbar = "Bulbar Onset",
               Age_Onset = "Age of Onset",
               APOE4_State = "APOE4 Status"
             )) +
  labs(x = "HR") +
#  scale_x_continuous(limits = c(0,3.5)) +
  ggtitle("Bulbar Onset + Age of Onset + APOE4 Status")
cox.zph(surv_obj)

# NYCG (APOE)
surv_obj <- coxph(formula = Surv(time = Disease_Duration, event = Status) ~ Onset_Bulbar + Age_Onset + APOE, data = NYCG)
ggcoef_model(surv_obj,
             add_reference_rows = T,
             exponentiate = T,
             variable_labels = c(
               Onset_Bulbar = "Bulbar Onset",
               Age_Onset = "Age of Onset",
               APOE = "APOE Genotype"
             )) +
  labs(x = "HR") +
#  scale_x_continuous(limits = c(0,3.5)) +
  ggtitle("Bulbar Onset + Age of Onset + APOE Genotype")
cox.zph(surv_obj)

# NYCG (SUBTYPE INCLUDED - APOE4_State)
surv_obj <- coxph(formula = Surv(time = Disease_Duration, event = Status) ~ onset_bulbar + Age_Onset + APOE4_State + Subtype, data = survival_data[survival_data$Disease_Duration < 80,])
ggforest(surv_obj, data = survival_data)
ggcoef_model(surv_obj,
             add_reference_rows = T,
             exponentiate = T,
             variable_labels = c(
               onset_bulbar = "Bulbar Onset",
               Age_Onset = "Age of Onset",
               APOE4_State = "APOE4 Status"
             )) +
  labs(x = "HR") +
  #  scale_x_continuous(limits = c(0,3.5)) +
  ggtitle("Bulbar Onset + Age of Onset + APOE4 Status + ALS Subtype")
cox.zph(surv_obj)

# NYCG (SUBTYPE INCLUDED - APOE)
surv_obj <- coxph(formula = Surv(time = Disease_Duration, event = Status) ~ onset_bulbar + Age_Onset + APOE + Subtype, data = survival_data[survival_data$Disease_Duration < 80,])
ggforest(surv_obj, data = survival_data)
ggcoef_model(surv_obj,
             add_reference_rows = T,
             exponentiate = T,
             variable_labels = c(
               onset_bulbar = "Bulbar Onset",
               Age_Onset = "Age of Onset",
               APOE = "APOE Genotype"
             )) +
  labs(x = "HR") +
  #  scale_x_continuous(limits = c(0,3.5)) +
  ggtitle("Bulbar Onset + Age of Onset + APOE Genotype + ALS Subtype")
cox.zph(surv_obj)

# NYCG (SUBTYPE STRATA - APOE4_State)
surv_obj <- coxph(formula = Surv(time = Disease_Duration, event = Status) ~ onset_bulbar + Age_Onset + APOE4_State + strata(Subtype), data = survival_data)
ggcoef_model(surv_obj,
             add_reference_rows = T,
             exponentiate = T,
             variable_labels = c(
               onset_bulbar = "Bulbar Onset",
               Age_Onset = "Age of Onset",
               APOE4_State = "APOE4 Status"
             )) +
  labs(x = "HR") +
  #  scale_x_continuous(limits = c(0,3.5)) +
  ggtitle("Bulbar Onset + Age of Onset + APOE4 Status + ALS Subtype + Strata(Subtype)")
cox.zph(surv_obj)

# NYCG (SUBTYPE STRATA - APOE)
surv_obj <- coxph(formula = Surv(time = Disease_Duration, event = Status) ~ onset_bulbar + Age_Onset + APOE + strata(Subtype), data = survival_data)
ggcoef_model(surv_obj,
             add_reference_rows = T,
             exponentiate = T,
             variable_labels = c(
               onset_bulbar = "Bulbar Onset",
               Age_Onset = "Age of Onset",
               APOE = "APOE Genotype"
             )) +
  labs(x = "HR") +
  #  scale_x_continuous(limits = c(0,3.5)) +
  ggtitle("Bulbar Onset + Age of Onset + APOE Genotype + ALS Subtype + Strata(Subtype)")
cox.zph(surv_obj)

### Tests ----------------------------------------------------------------------
# # NYCG (BULBAR AND LIMB SUBSET) Data
# test_NYCG_site <- NYCG[rownames(survival_data[survival_data$Site_Onset %in% c("Bulbar", "Limb"),]),]
# rownames(survival_data[survival_data$Site_Onset %in% c("Bulbar", "Limb"),])
# 
# surv_obj <- coxph(formula = Surv(time = Disease_Duration, event = Status) ~ Onset_Bulbar + Age_Onset + APOE + Gender, data = test_NYCG_site)
# ggcoef_model(surv_obj,
#              add_reference_rows = T,
#              exponentiate = T,
#              variable_labels = c(
#                Onset_Bulbar = "Bulbar Onset",
#                Age_Onset = "Age of Onset",
#                APOE = "APOE Genotype"
#              )) +
#   labs(x = "HR") +
# #  scale_x_continuous(limits = c(0,5)) +
#   ggtitle("Bulbar Onset + Age of Onset + APOE Genotype")
# cox.zph(surv_obj)
# 
# # NYCG (BULBAR AND LIMB SUBSET - APOE4_State)
# surv_obj <- coxph(formula = Surv(time = Disease_Duration, event = Status) ~ Onset_Bulbar + Age_Onset + APOE4_State + Gender, data = test_NYCG_site)
# ggcoef_model(surv_obj,
#              add_reference_rows = T,
#              exponentiate = T,
#              variable_labels = c(
#                Onset_Bulbar = "Bulbar Onset",
#                Age_Onset = "Age of Onset",
#                APOE4_State = "APOE4 Status"
#              )) +
#   labs(x = "HR") +
#   scale_x_continuous(limits = c(0,3.5)) +
#   ggtitle("Bulbar Onset + Age of Onset + APOE4 Status")
# cox.zph(surv_obj)
# 
# # NYCG (DISEASE DURATION < 80 SUBSET) Data 
# NYCG_dd_80 <- NYCG[NYCG$Disease_Duration < 80,]
# 
# # NYCG (DISEASE DURATION < 80 SUBSET - APOE)
# surv_obj <- coxph(formula = Surv(time = Disease_Duration, event = Status) ~ Onset_Bulbar + Age_Onset + APOE + Gender, data = NYCG_dd_80)
# ggforest(surv_obj)
# ggcoef_model(surv_obj,
#              add_reference_rows = T,
#              exponentiate = T,
#              variable_labels = c(
#                Onset_Bulbar = "Bulbar Onset",
#                Age_Onset = "Age of Onset",
#                APOE = "APOE Genotype"
#              )) +
#   labs(x = "HR") +
#   #  scale_x_continuous(limits = c(0,5)) +
#   ggtitle("Bulbar Onset + Age of Onset + APOE Genotype")
# cox.zph(surv_obj)
# 
# # NYCG (DISEASE DURATION < 80 SUBSET - APOE4_State)
# surv_obj <- coxph(formula = Surv(time = Disease_Duration, event = Status) ~ Onset_Bulbar + Age_Onset + APOE4_State, data = NYCG_dd_80)
# ggforest(surv_obj)
# ggcoef_model(surv_obj,
#              add_reference_rows = T,
#              exponentiate = T,
#              variable_labels = c(
#                Onset_Bulbar = "Bulbar Onset",
#                Age_Onset = "Age of Onset",
#                APOE4_State = "APOE4 Status"
#              )) +
#   labs(x = "HR") +
# #  scale_x_continuous(limits = c(0,3.5)) +
#   ggtitle("Bulbar Onset + Age of Onset + APOE4 Status")
# cox.zph(surv_obj)
# 
# # UGLY PLOT?
# ggplot(NYCG,
#        aes(Status, Disease_Duration)) +
#   geom_beeswarm() +
#   geom_boxplot()

### Correlation ----------------------------------------------------------------
# Correlation between Age of Onset and Disease Duration between genotypes
ggscatter(UTRECHT[UTRECHT$Disease_Duration < 80,],
          "Disease_Duration",
          "Age_Onset",
          add = "reg.line",
          conf.int = TRUE,
          cor.coef = TRUE,
          cor.method = "spearman")

ggscatter(UTRECHT[UTRECHT$APOE4_State == "APOE4neg" & UTRECHT$Disease_Duration < 80,],
          "Disease_Duration",
          "Age_Onset",
          add = "reg.line",
          conf.int = TRUE,
          cor.coef = TRUE,
          cor.method = "spearman",)

ggscatter(UTRECHT[UTRECHT$APOE4_State == "APOE4pos" & UTRECHT$Disease_Duration < 80,],
          "Disease_Duration",
          "Age_Onset",
          add = "reg.line",
          conf.int = TRUE,
          cor.coef = TRUE,
          cor.method = "spearman",)

ggqqplot(UTRECHT,
         "Disease_Duration")
ggqqplot(UTRECHT,
         "Age_Onset")


# gse_extra <- getGEO(filename = "GSE116622_series_matrix.txt")
# metadata_extra <- gse_extra@phenoData@data
# colnames(metadata_extra) <- sub(".", "_", colnames(metadata_extra), fixed=TRUE)
# colnames(metadata_extra) <- sub(":", "_", colnames(metadata_extra), fixed=TRUE)
# colnames(metadata_extra) <- sub(" ", "_", colnames(metadata_extra), fixed=TRUE)
