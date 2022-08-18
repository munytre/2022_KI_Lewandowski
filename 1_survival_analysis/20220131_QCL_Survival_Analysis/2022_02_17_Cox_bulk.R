###--------------------------------------------------------------------------###
### Project: ApoE in ALS                                                       -
### Purpose: Multivariate analyses with cutoff groups                          -
### Authors: FS Sanders, QC Lin                                                -
### Last Edit: 2022-02-21, QC Lin                                              -
### Contact: q.c.lin@students.uu.nl or sebastian.lewandowski@ki.se             -
###--------------------------------------------------------------------------###
### Comment:                                                                   -
### The following script is used for creating the final three cox models for   -
### the SPP1 and COL6a1 analyses. Should be working if ran all at once, for    -
### any questions contact me through Sebastian. I removed other frailty        -
### testing and cohort specific models from this script. After every cox       -
### figure is finished Schoenfeld test results are printed in the console.     -
### ~ F.S.Sanders                                                              -
###--------------------------------------------------------------------------###

### Dependencies and general settings ------------------------------------------

# Load libraries
library(tidyverse)
library(gridExtra)
library(survminer)
library(survival)
library(broom)
library(rlang)
library(GGally)

# Set working directory and seed for reproducibility
setwd(".")
set.seed(747)

# Create output directory for pdfs
ifelse(!dir.exists(file.path(getwd(), "output")),
       dir.create(file.path(getwd(), "output_R")),
       FALSE)

# Helper function for method of approach (generation of Cox models)
.get_data <- function(fit, data = NULL, complain = TRUE) {
  if(is.null(data)){
    if (complain)
      warning ("The `data` argument is not provided. Data will be extracted from model fit.")
    data <- eval(fit$call$data)
    if (is.null(data))
      stop("The `data` argument should be provided either to ggsurvfit or survfit.")
  }
  data
}

### Load data ------------------------------------------------------------------

# Load from rds objects
bead_info_utrecht <- readRDS(file = "generated_data_R/bead_info_utrecht.rds")
full_beads_utrecht <- readRDS(file = "generated_data_R/als.ma.apoe.rds")
als_beads_utrecht <- readRDS(file = "generated_data_R/apoeALS.rds")
# For the following two, need to run "2020_02_24_Cutpoint_Analyses_Publ.R" first
bead_info_utrecht_pvals_loose <- read.csv("output_R/bead_info_utrecht_pvals_loose.csv")
bead_info_utrecht_pvals_e34 <- read.csv("output_R/bead_info_utrecht_pvals_strict_e34.csv")

### Prep data (and generate cutoff dataframes) ---------------------------------

# Change names for readability
columns_of_interest <- c("survobj",
                         c(as.character(bead_info_utrecht$Bead_name)),
                         "gender",
                         "onset_bulbar",
                         "sampling_age",
                         "gap",
                         "apoe",
                         "cohort")

# Gather data and implement APOE genotype
allcoxdf <- als_beads_utrecht %>%
  select(all_of(columns_of_interest)) %>%
  mutate(survobj,
            Gender = gender,
            Onset = factor(onset_bulbar, labels = c("Thorasic/Spinal", "Bulbar")),
            Sampling_Age = sampling_age,
            Sampling_delay = gap,
            ApoE = factor(apoe),
            cohort,
            .keep = "unused") %>%
  filter(!is.na(ApoE)) %>%
  filter(ApoE == "ALS Apo-e3e3" | ApoE == "ALS Apo-e3e4") %>%
  droplevels()

# Making cutoffs for beads of interest
# for (i in c(as.character(bead_info_utrecht$Bead_name))){
#   single.cut <- surv_cutpoint(als_beads_utrecht,
#                               time = "survival_months",
#                               event = "died",
#                               variables = i,
#                               minprop = 0.1)
#   
#   max <- single.cut[[i]]
#   cut <- max$estimate
#   column_level <- paste0(i,"_Level")
#   
#   allcoxdf <- allcoxdf %>% mutate(level = case_when(
#     allcoxdf[i] <= cut ~ "Lower", 
#     allcoxdf[i] > cut ~ "Upper"))
#   allcoxdf$level <- as.factor(allcoxdf$level)
#   
#   names(allcoxdf)[names(allcoxdf) == "level"] <- column_level
#   }

### Multivariate Cox (Bulk - Continuous) ---------------------------------------
for (protein in c(as.character(bead_info_utrecht$Bead_name))){
  # Building the model
  res.clu <- expr(coxph(survobj ~ ApoE*!!sym(protein) + Onset + Sampling_Age + Gender +
                          strata(Sampling_delay), data = allcoxdf))
  res_cox <- as.formula(res.clu)
  summary_res.clu <- summary(res_cox)
  # Schoenfeld residuals test
  zph_res.clu <- cox.zph(res_cox)
  cox_results <- list(summary_res.clu,zph_res.clu)
  # Write files
  capture.output(cox_results, file = paste0("output_R/cox_multi_interaction/cox_results_", protein, ".txt"))
  # Output forest plots
  forestplot <- ggcoef_model(res_cox,
                             add_reference_rows = F,
                             exponentiate = T) +
    labs(x = "HR")
  cairo_pdf(paste0("output_R/cox_multi_interaction/forest_", protein, ".pdf"),
      width = 10,
      height = 7)
  grid.draw(forestplot)
  dev.off()
}

### Multivariate Cox (Genotyped groups) ----------------------------------------
# For APOEe3e3
e33coxdf <- allcoxdf %>%
  filter(ApoE == "ALS Apo-e3e3") %>%
  droplevels()

for (protein in c(as.character(bead_info_utrecht$Bead_name))){
  # Building the model
  protein_level <- paste0(protein,"_Level")
  res.clu <- expr(coxph(survobj ~ !!sym(protein) + Onset + Sampling_Age + Gender +
                          strata(Sampling_delay), data = e33coxdf))
  res_cox <- as.formula(res.clu)
  summary_res.clu <- summary(res_cox)
  # Schoenfeld residuals test
  zph_res.clu <- cox.zph(res_cox)
  cox_results <- list(summary_res.clu,zph_res.clu)
  # Write files
  capture.output(cox_results, file = paste0("output_R/cox_multi_e33/cox_results_e33_", protein, ".txt"))
  # Output forest plots
  forestplot <- ggcoef_model(res_cox,
                             add_reference_rows = F,
                             exponentiate = T) +
    labs(x = "HR")
  cairo_pdf(paste0("output_R/cox_multi_e33/forest_e33_", protein, ".pdf"),
            width = 10,
            height = 7)
  grid.draw(forestplot)
  dev.off()
}

# For APOEe3e4
e34coxdf <- allcoxdf %>%
  filter(ApoE == "ALS Apo-e3e4") %>%
  droplevels()

for (protein in c(as.character(bead_info_utrecht$Bead_name))){
  # Building the model
  protein_level <- paste0(protein,"_Level")
  res.clu <- expr(coxph(survobj ~ !!sym(protein) + Onset + Sampling_Age + Gender +
                          strata(Sampling_delay), data = e34coxdf))
  res_cox <- as.formula(res.clu)
  summary_res.clu <- summary(res_cox)
  # Schoenfeld residuals test
  zph_res.clu <- cox.zph(res_cox)
  cox_results <- list(summary_res.clu,zph_res.clu)
  # Write files
  capture.output(cox_results, file = paste0("output_R/cox_multi_e34/cox_results_e34_", protein, ".txt"))
  # Output forest plots
  forestplot <- ggcoef_model(res_cox,
                             add_reference_rows = F,
                             exponentiate = T) +
    labs(x = "HR")
  cairo_pdf(paste0("output_R/cox_multi_e34/forest_e34_", protein, ".pdf"),
            width = 10,
            height = 7)
  grid.draw(forestplot)
  dev.off()
}

### Multivariate Cox (Bulk - Categorical) --------------------------------------
# for (protein in c(as.character(bead_info_utrecht$Bead_name))){
#   # Building the model
#   protein_level <- paste0(protein,"_Level")
#   res.clu <- expr(coxph(survobj ~ ApoE + !!sym(protein_level) + Onset + Sampling_Age + !!sym(protein_level):ApoE +
#                   strata(Sampling_delay), data = allcoxdf))
#   res_cox <- as.formula(res.clu)
#   summary_res.clu <- summary(res_cox)
#   # Schoenfeld residuals test
#   zph_res.clu <- cox.zph(res_cox)
#   cox_results <- list(summary_res.clu,zph_res.clu)
#   # Write files
#   capture.output(cox_results, file = paste0("output_R/cox_multi_result_test/cox_results_", protein, ".txt"))
# }

### Multivariate Cox (Bulk - TEST) ---------------------------------------------
# for (protein in c(as.character(bead_info_utrecht$Bead_name))){
#   # Building the model
#   res.clu <- expr(coxph(survobj ~ ApoE + ApoE:!!sym(protein) + Onset + Sampling_Age + Gender +
#                           strata(Sampling_delay), data = allcoxdf))
#   res_cox <- as.formula(res.clu)
#   summary_res.clu <- summary(res_cox)
#   # Schoenfeld residuals test
#   zph_res.clu <- cox.zph(res_cox)
#   cox_results <- list(summary_res.clu,zph_res.clu)
#   # Write files
#   capture.output(cox_results, file = paste0("output_R/cox_multi_interaction/test_cox_results_", protein, ".txt"))
#   # Output forest plots
#   forestplot <- ggcoef_model(res_cox,
#                              add_reference_rows = F,
#                              exponentiate = T) +
#     labs(x = "HR")
#   cairo_pdf(paste0("output_R/cox_multi_interaction/test_forest_", protein, ".pdf"),
#             width = 10,
#             height = 7)
#   grid.draw(forestplot)
#   dev.off()
# }
# 
### Multivariate Cox (Bulk - TEST_NOINTERACTION) -------------------------------
# for (protein in c(as.character(bead_info_utrecht$Bead_name))){
#   # Building the model
#   res.clu <- expr(coxph(survobj ~ ApoE + !!sym(protein) + Onset + Sampling_Age + Gender +
#                           strata(Sampling_delay), data = allcoxdf))
#   res_cox <- as.formula(res.clu)
#   summary_res.clu <- summary(res_cox)
#   # Schoenfeld residuals test
#   zph_res.clu <- cox.zph(res_cox)
#   cox_results <- list(summary_res.clu,zph_res.clu)
#   # Write files
#   capture.output(cox_results, file = paste0("output_R/cox_multi_interaction/test_noint_cox_results_", protein, ".txt"))
#   # Output forest plots
#   forestplot <- ggcoef_model(res_cox,
#                              add_reference_rows = F,
#                              exponentiate = T) +
#     labs(x = "HR")
#   cairo_pdf(paste0("output_R/cox_multi_interaction/test_noint_forest_", protein, ".pdf"),
#             width = 10,
#             height = 7)
#   grid.draw(forestplot)
#   dev.off()
# }
# 
# res.clu <- coxph(survobj ~ ApoE*B.33_LAMP1_ab24170 + Onset + Sampling_Age + Gender +
#                    strata(Sampling_delay), data = allcoxdf)
# ggadjustedcurves(res.clu, data = allcoxdf)
