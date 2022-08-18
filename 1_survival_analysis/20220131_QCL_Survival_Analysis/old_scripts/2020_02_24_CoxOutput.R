###--------------------------------------------------------------------------###
### Project: ApoE in ALS                                                       -
### Purpose: Multivariate analyses with cutoff groups                          -
### Authors: FS Sanders, QC Lin                                                -
### Last Edit: 2022-02-09, QC Lin                                              -
### Contact: q.c.lin@students.uu.nl or sebastian.lewandowski@ki.se             -
###--------------------------------------------------------------------------###
### Comment:                                                                   -
### The following script is used for creating the final three coxmodels for    -
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
library(grid)
library(forestmodel)
library(rlang)

# Set working directory and seed for reproducibility
setwd("~/Desktop/20220131_QCL_Survival_Analysis/")
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

# General parameters for forest plot
cpositions <- c(0.02, 0.22, 0.4)
fontsize <-  0.7
noDigits <- 2

# Set fontsize of annotation in forest plots
annot_size_mm <- fontsize *
  as.numeric(convertX(unit(theme_get()$text$size, "pt"), "mm"))

### Load data ------------------------------------------------------------------

# Load from rds objects
bead_info_utrecht_pvals_loose <- read.csv("~/Desktop/20220131_QCL_Survival_Analysis/output_R/bead_info_utrecht_pvals_loose.csv")
bead_info_utrecht_pvals_e34 <- read.csv("~/Desktop/20220131_QCL_Survival_Analysis/output_R/bead_info_utrecht_pvals_strict_e34.csv")
bead_info_utrecht <- readRDS(file = "generated_data_R/bead_info_utrecht.rds")
full_beads_utrecht <- readRDS(file = "generated_data_R/als.ma.apoe.rds")
als_beads_utrecht <- readRDS(file = "generated_data_R/apoeALS.rds")

### Generate cutoff dataframes -------------------------------------------------

# Change names for readability
columns_of_interest <- c("survobj",
                         c(bead_info_utrecht_pvals_e34$Bead_name),
                         "gender",
                         "onset_bulbar",
                         "sampling_age",
                         "gap",
                         "apoe",
                         "cohort")

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
for (i in c(bead_info_utrecht_pvals_e34$Bead_name)){
  single.cut <- surv_cutpoint(als_beads_utrecht,
                              time = "survival_months",
                              event = "died",
                              variables = i,
                              minprop = 0.1)
  
  max <- single.cut[[i]]
  cut <- max$estimate
  column_level <- paste0(i,"_Level")
  
  allcoxdf <- allcoxdf %>% mutate(level = case_when(
    allcoxdf[i] <= cut ~ "Lower", 
    allcoxdf[i] > cut ~ "Upper"))
  
  names(allcoxdf)[names(allcoxdf) == "level"] <- column_level
  }


##
single.cut_33 <- surv_cutpoint(filter(als_beads_utrecht, apoe == "ALS Apo-e3e3"),
                            time = "survival_months",
                            event = "died",
                            variables = c("B.183_AIF1_HPA049234"),
                            minprop = 0.1)
max183_33 <- single.cut_33[["B.183_AIF1_HPA049234"]]
cut183_33 <- max183_33$estimate

single.cut_34 <- surv_cutpoint(filter(als_beads_utrecht, apoe == "ALS Apo-e3e4"),
                            time = "survival_months",
                            event = "died",
                            variables = c("B.183_AIF1_HPA049234"),
                            minprop = 0.1)
max183_34 <- single.cut_34[["B.183_AIF1_HPA049234"]]
cut183_34 <- max183_34$estimate

single.cut <- surv_cutpoint(als_beads_utrecht,
                               time = "survival_months",
                               event = "died",
                               variables = c("B.183_AIF1_HPA049234"),
                               minprop = 0.1)
max183 <- single.cut[["B.183_AIF1_HPA049234"]]
cut183 <- max183$estimate
##
# single.cut <- surv_cutpoint(combinedALS,
#                             time = "survival_months",
#                             event = "died", 
#                             variables = c("RawB.123"),
#                             minprop = 0.1)
# max123 <- single.cut[["RawB.123"]]
# cut123 <- max123$estimate
# single.cut <- surv_cutpoint(combinedALS,
#                             time = "survival_months",
#                             event = "died",
#                             variables = c("RawB.94"),
#                             minprop = 0.0)
# max94 <- single.cut[["RawB.94"]]
# cut94 <- max94$estimate

# Stratifying on earlier establish cutoffs.
##
#e33
allcoxdf <- allcoxdf %>% mutate(AIF1_Level_e33 = case_when(
  allcoxdf$AIF1 <= cut183_33 ~ "Lower", 
  allcoxdf$AIF1 > cut183_33 ~ "Upper"))
#e34
allcoxdf <- allcoxdf %>% mutate(AIF1_Level_e34 = case_when(
  allcoxdf$AIF1 <= cut183_34 ~ "Lower", 
  allcoxdf$AIF1 > cut183_34 ~ "Upper"))
#non
allcoxdf <- allcoxdf %>% mutate(AIF1_Level = case_when(
  allcoxdf$AIF1 <= cut183 ~ "Lower", 
  allcoxdf$AIF1 > cut183 ~ "Upper"))
##
# allcoxdf <- allcoxdf %>% mutate(SPP1_Level = case_when(
#   allcoxdf$SPP1 <= cut53 ~ "Lower", 
#   allcoxdf$SPP1 > cut53 ~ "Upper"))
# allcoxdf <- allcoxdf %>% mutate(COL6A1_Level = case_when(
#   allcoxdf$COL6A1 <= cut123 ~ "Lower", 
#   allcoxdf$COL6A1 > cut123 ~ "Upper"))
# allcoxdf <- allcoxdf %>% mutate(NEFL_Level = case_when(
#   allcoxdf$NEFL <= cut94 ~ "Lower", 
#   allcoxdf$NEFL > cut94 ~ "Upper"))

### Multivariate Cox (Bulk) ----------------------------------------------------
for (protein in c(bead_info_utrecht_pvals_e34$Bead_name)){
  protein_level <- paste0(protein,"_Level")
  res.clu <- expr(coxph(survobj ~ !!sym(protein_level):ApoE + Onset + Sampling_Age +
                     strata(Sampling_delay), data = allcoxdf))
  res_cox <- as.formula(res.clu)
  summary_res.clu <- summary(res_cox)
  zph_res.clu <- cox.zph(res_cox)
  cox_results <- list(summary_res.clu,zph_res.clu)
  capture.output(cox_results, file = paste0("output_R/cox_multi_result/cox_results_", protein, ".txt"))
}

(!!as.symbol(var))
# Model creation for the covariates
res.clu <- coxph(survobj ~ B.183_AIF1_HPA049234_Level*ApoE + Onset + Sampling_Age +
                   strata(Sampling_delay), data = allcoxdf)
summary_res.clu <- summary(res.clu)
zph_res.clu <- cox.zph(res.clu)
cox_results <- list(summary_res.clu,zph_res.clu)
capture.output(cox_results, file = paste0("output_R/cox_results_", protein, ".txt"))
forestmodel::forest_model(res.clu)

# Schoenfeld residuals test
cat("Multivariate Cox model finished, proportionality for covariates are tested. See below.\n")
cox.zph(res.clu)

publish(res_cox)

### Continuous Univariate Cox --------------------------------------------------

# Title for forest plot
main <-  "Univariate Cox Models: Hazard ratio"
# Covariates
covariates <- c("SPP1", "COL6A1", "NEFL", "Gender", "Onset", "Sampling_Age")

# Model creation for each covariate
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('survobj ~', x, '+ strata(Sampling_delay) + cluster(cohort)')))
univ_models <- lapply(univ_formulas, function(x){coxph(x, data = allcoxdf)})

# Extract values from models
univ_results <- lapply(univ_models,
                       function(y){ 
                         coef <- as.data.frame(tidy(y))
                         x <- summary(y)
                         N <- x$n
                         p.value<-signif(x$wald["pvalue"], digits=3)
                         wald.test<-signif(x$wald["test"], digits=2)
                         estimate <-signif(x$coef[1], digits=3);#coeficient beta
                         estimate.1 <-signif(x$coef[2], digits=3);#exp(beta)
                         conf.low <- signif(x$conf.int[,"lower .95"], 3)
                         conf.high <- signif(x$conf.int[,"upper .95"],3)
                         conf.low.1 <- signif(x$conf.int[,"lower .95"], 3)
                         conf.high.1 <- signif(x$conf.int[,"upper .95"],3)
                         ci <- paste0("(", 
                                      conf.low.1, "-", conf.high.1, ")")
                         res<-c(N, p.value, estimate, estimate.1, conf.low,
                                conf.high, conf.low.1, conf.high.1, ci)
                         names(res)<-c( "N","p.value", "estimate", "estimate.1", "conf.low",
                                        "conf.high","conf.low.1", "conf.high.1", "ci")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res <- t(as.data.frame(univ_results, check.names = FALSE))
res <- tibble::rownames_to_column(as.data.frame(res), "var")

# Correct for group size
datalist <- list()
covariates <- c("SPP1_Level", "COL6A1_Level", "NEFL_Level", "Gender", "Onset")
res$N <- as.numeric(as.character(res$N))
for (i in 1:length(covariates)){
  var <- covariates[i]
  adf <- as.data.frame(table(allcoxdf[, var]))
  datalist[[i]] <- cbind(var = var, adf, pos = 1:nrow(adf))
}
total <- do.call(rbind,datalist)
n <- subset(total, pos == 2) %>% 
  select(Freq)
res$N[4:5] <- n[4:5,]
res$N <- paste0("(N=",res$N,")")

# Change datatypes
res$N <- as.character(res$N)
res$p.value <- as.numeric(as.character(res$p.value))
res$estimate <- as.numeric(as.character(res$estimate))
res$estimate.1 <- as.character(res$estimate.1)
res$conf.low <- as.numeric(as.character(res$conf.low))
res$conf.high <- as.numeric(as.character(res$conf.high))
res$conf.low.1 <- as.character(res$conf.low.1)
res$conf.high.1 <- as.character(res$conf.high.1)
res$ci <- as.character(res$ci)

# Add stars (p.value) to dataframe
res$stars <- paste0(round(res$p.value, 3), " ",
                    ifelse(res$p.value < 0.05, "*",""),
                    ifelse(res$p.value < 0.01, "*",""),
                    ifelse(res$p.value < 0.001, "*",""))
res$stars[which(res$p.value < 0.001)] = "<0.001 ***"

# Flip order
res <- res[nrow(res):1, ]

# Setting plot sizes and scales
rangeb <- range(log(res$conf.low), log(res$conf.high), na.rm = TRUE)
breaks <- axisTicks(rangeb/2, log = TRUE, nint = 7)
rangeplot <- rangeb
# Make plot twice as wide as needed to create space for annotations
rangeplot[1] <- rangeplot[1] - diff(rangeb)
# Increase white space on right for p-vals
rangeplot[2] <- rangeplot[2] + .15 * diff(rangeb)
width <- diff(rangeplot)
# Y-coordinates for labels
y_variable <- rangeplot[1] +  cpositions[1] * width
y_nlevel <- rangeplot[1]  +  cpositions[2] * width
y_cistring <- rangeplot[1]  +  cpositions[3] * width
y_stars <- rangeb[2]
x_annotate <- seq_len(nrow(res))

# Change column names for readability
res$var[res$var == "Onset"] <- "Onset Bulbar"
res$var[res$var == "Gender"] <- "Gender Male"
res$var[res$var == "Sampling_Age"] <- "Sampling Age"

# Build forest plot
p <- ggplot(res, aes(seq_along(var), exp(estimate))) +
  geom_rect(aes(xmin = seq_along(var) - .5, xmax = seq_along(var) + .5,
                ymin = exp(rangeplot[1]), ymax = exp(rangeplot[2]),
                fill = ordered(seq_along(var) %% 2 + 1))) +
  scale_fill_manual(values = c("#FFFFFF33", "#00000033"), guide = "none") +
  geom_point(pch = 15, size = 4) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.15) +
  geom_hline(yintercept = 1, linetype = 3) +
  coord_flip(ylim = exp(rangeplot)) +
  ggtitle(main) +
  scale_y_log10(
    name = "",
    labels = sprintf("%g", breaks),
    expand = c(0.02, 0.02),
    breaks = breaks) +
  theme_minimal() +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position = "none",
        panel.border=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  xlab("") +
  annotate(geom = "text", x = x_annotate, y = exp(y_variable),
           label = res$var, fontface = "bold", hjust = 0,
           size = annot_size_mm) +
  annotate(geom = "text", x = x_annotate, y = exp(y_nlevel),
           label = res$N, fontface = "italic", hjust = 0,
           size = annot_size_mm) +
  annotate(geom = "text", x = x_annotate, y = exp(y_cistring),
           label = res$estimate.1, size = annot_size_mm) +
  annotate(geom = "text", x = x_annotate, y = exp(y_cistring),
           label = res$ci, size = annot_size_mm,
           vjust = 2,  fontface = "italic") +
  annotate(geom = "text", x = x_annotate, y = exp(y_stars),
           label = res$p.value, size = annot_size_mm,
           hjust = -0.2,  fontface = "italic") +
  annotate(geom = "text", x = 0.5, y = exp(y_variable),
           label = paste0("Univariate Cox Models. Continuous data."),
           size = annot_size_mm, hjust = 0, vjust = 1.2,  fontface = "italic")

# Switch off clipping for p-vals, bottom annotation
gt <- ggplot_gtable(ggplot_build(p))
gt$layout$clip[gt$layout$name == "panel"] <- "off"
grid.draw(gt)

# Print forest plot
pdf(file = "output_R/2020_Univ_Cox_Continuous.pdf", width = 10, height = 7)
ggpubr::as_ggplot(gt)
dev.off()

# Schoenfeld residuals tests
cat("Univariate Cox model finished, proportionality for all univariate models are tested seperately. See below.\n")
unitest <- coxph(survobj ~ SPP1 + strata(Sampling_delay) + cluster(cohort), data = allcoxdf)
cox.zph(unitest)
unitest2 <- coxph(survobj ~ COL6A1+ strata(Sampling_delay) + cluster(cohort), data = allcoxdf)
cox.zph(unitest2)
unitest3 <- coxph(survobj ~ NEFL + strata(Sampling_delay) + cluster(cohort), data = allcoxdf)
cox.zph(unitest3)

### Continuous Univariate Cox with Cutoff --------------------------------------

# Title for forest plot
main <-  "Univariate Cox Models: Hazard ratio"
# Covariates
covariates <- c("SPP1_Level", "COL6A1_Level", "NEFL_Level", "Gender", "Onset", "Sampling_Age")

# Model creation for each covariate
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('survobj ~', x, '+ strata(Sampling_delay) + cluster(cohort)')))
univ_models <- lapply(univ_formulas, function(x){coxph(x, data = allcoxdf)})

# Extract values from models
univ_results <- lapply(univ_models,
                       function(y){ 
                         coef <- as.data.frame(tidy(y))
                         # conf.low <- signif(coef$conf.low, digits = 3)
                         # conf.high <- signif(coef$conf.high, digits = 3)
                         x <- summary(y)
                         N <- x$n
                         p.value <- signif(x$wald["pvalue"], digits=3)
                         wald.test <- signif(x$wald["test"], digits=2)
                         estimate <- signif(x$coef[1], digits=3);#coeficient beta
                         estimate.1 <-signif(x$coef[2], digits=3);#exp(beta)
                         conf.low <- signif(x$conf.int[,"lower .95"], 3)
                         conf.high <- signif(x$conf.int[,"upper .95"],3)
                         conf.low.1 <- signif(x$conf.int[,"lower .95"], 3)
                         conf.high.1 <- signif(x$conf.int[,"upper .95"],3)
                         ci <- paste0("(", 
                                      conf.low.1, "-", conf.high.1, ")")
                         res <- c(N, p.value, estimate, estimate.1, conf.low,
                                  conf.high, conf.low.1, conf.high.1, ci)
                         names(res) <- c( "N","p.value", "estimate", "estimate.1", "conf.low",
                                          "conf.high","conf.low.1", "conf.high.1", "ci")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res <- t(as.data.frame(univ_results, check.names = FALSE))
res <- tibble::rownames_to_column(as.data.frame(res), "var")

# Correct for group size
datalist <- list()
covariates <- c("SPP1_Level", "COL6A1_Level", "NEFL_Level", "Gender", "Onset")
res$N <- as.numeric(as.character(res$N))
for (i in 1:length(covariates)){
  var <- covariates[i]
  adf <- as.data.frame(table(allcoxdf[, var]))
  datalist[[i]] <- cbind(var = var, adf, pos = 1:nrow(adf))
}
total <- do.call(rbind,datalist)
n <- subset(total, pos == 2) %>% 
  select(Freq)
res$N[1:5] <- n[1:5,]
res$N <- paste0("(N=",res$N,")")

# Change datatypes
res$p.value <- as.numeric(as.character(res$p.value))
res$estimate <- as.numeric(as.character(res$estimate))
res$estimate.1 <- as.character(res$estimate.1)
res$conf.low <- as.numeric(as.character(res$conf.low))
res$conf.high <- as.numeric(as.character(res$conf.high))
res$conf.low.1 <- as.character(res$conf.low.1)
res$conf.high.1 <- as.character(res$conf.high.1)
res$ci <- as.character(res$ci)

# Add stars (p.value) to dataframe
res$stars <- paste0(round(res$p.value, 3), " ",
                    ifelse(res$p.value < 0.05, "*",""),
                    ifelse(res$p.value < 0.01, "*",""),
                    ifelse(res$p.value < 0.001, "*",""))
res$stars[which(res$p.value < 0.001)] = "<0.001 ***"

# Flip order
res <- res[nrow(res):1, ]

# Setting plot sizes and scales
rangeb <- range(log(res$conf.low), log(res$conf.high), na.rm = TRUE)
breaks <- axisTicks(rangeb/2, log = TRUE, nint = 7)
rangeplot <- rangeb
# Make plot twice as wide as needed to create space for annotations
rangeplot[1] <- rangeplot[1] - diff(rangeb)
# Increase white space on right for p-vals
rangeplot[2] <- rangeplot[2] + .15 * diff(rangeb)
width <- diff(rangeplot)
# Y-coordinates for labels
y_variable <- rangeplot[1] +  cpositions[1] * width
y_nlevel <- rangeplot[1]  +  cpositions[2] * width
y_cistring <- rangeplot[1]  +  cpositions[3] * width
y_stars <- rangeb[2]
x_annotate <- seq_len(nrow(res))

# CHange column names for readability
res$var[res$var == "Onset"] <- "Onset Bulbar"
res$var[res$var == "Gender"] <- "Gender Male"
res$var[res$var == "SPP1_Level"] <- "High level SPP1"
res$var[res$var == "COL6A1_Level"] <- "High level COL6A1"
res$var[res$var == "NEFL_Level"] <- "High level NEFL"
res$var[res$var == "Sampling_Age"] <- "Sampling Age"

# Build forest plot
p <- ggplot(res, aes(seq_along(var), exp(estimate))) +
  geom_rect(aes(xmin = seq_along(var) - .5, xmax = seq_along(var) + .5,
                ymin = exp(rangeplot[1]), ymax = exp(rangeplot[2]),
                fill = ordered(seq_along(var) %% 2 + 1))) +
  scale_fill_manual(values = c("#FFFFFF33", "#00000033"), guide = "none") +
  geom_point(pch = 15, size = 4) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.15) +
  geom_hline(yintercept = 1, linetype = 3) +
  coord_flip(ylim = exp(rangeplot)) +
  ggtitle(main) +
  scale_y_log10(
    name = "",
    labels = sprintf("%g", breaks),
    expand = c(0.02, 0.02),
    breaks = breaks) +
  theme_minimal() +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position = "none",
        panel.border=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  xlab("") +
  annotate(geom = "text", x = x_annotate, y = exp(y_variable),
           label = res$var, fontface = "bold", hjust = 0,
           size = annot_size_mm) +
  annotate(geom = "text", x = x_annotate, y = exp(y_nlevel),
           label = res$N, fontface = "italic", hjust = 0,
           size = annot_size_mm) +
  annotate(geom = "text", x = x_annotate, y = exp(y_cistring),
           label = res$estimate.1, size = annot_size_mm) +
  annotate(geom = "text", x = x_annotate, y = exp(y_cistring),
           label = res$ci, size = annot_size_mm,
           vjust = 2,  fontface = "italic") +
  annotate(geom = "text", x = x_annotate, y = exp(y_stars),
           label = res$p.value, size = annot_size_mm,
           hjust = -0.2,  fontface = "italic") +
  annotate(geom = "text", x = 0.5, y = exp(y_variable),
           label = paste0("Univariate Cox Models.\nHigh protein level categories based on cutpoints."),
           size = annot_size_mm, hjust = 0, vjust = 1.2,  fontface = "italic")

# Switch off clipping for p-vals, bottom annotation
gt <- ggplot_gtable(ggplot_build(p))
gt$layout$clip[gt$layout$name == "panel"] <- "off"
grid.draw(gt)

# Print forest plot
pdf(file = "output_R/2020_Univ_Thresh_Cox.pdf", width = 10, height = 7)
ggpubr::as_ggplot(gt)
dev.off()

# Schoenfeld residuals tests
cat("Univariate Cox model (Threshold) finished, proportionality for all univariate models are tested seperately. See below.\n")
unitest <- coxph(survobj ~ SPP1_Level + strata(Sampling_delay)+ cluster(cohort), data = allcoxdf)
cox.zph(unitest)
unitest2 <- coxph(survobj ~ COL6A1_Level + strata(Sampling_delay)+ cluster(cohort), data = allcoxdf)
cox.zph(unitest2)
unitest3 <- coxph(survobj ~ NEFL_Level + strata(Sampling_delay)+ cluster(cohort), data = allcoxdf)
cox.zph(unitest3)

### Continuous Multivariate Cox ------------------------------------------------

# Title for forest plot
main <-  "Hazard Ratio, Continuous Protein Data"

# Model creation for the covariates
res.clu <- coxph(survobj ~ SPP1 +
                   Gender + Onset + Sampling_Age + strata(Sampling_delay) +
                   cluster(cohort), data = allcoxdf)

# Extract values from models
x <- summary(res.clu)
N <- x$n
coef <- as.data.frame(x$coefficients)
conf <- as.data.frame(x$conf.int)
p.value <- signif(coef["Pr(>|z|)"], digits = 3)
estimate <- signif(coef[1], digits=3);#coeficient beta
estimate.1 <- signif(coef[2], digits=3);#exp(beta)
conf.low <- signif(conf[3], 3)
conf.high <- signif(conf[4],3)
# Binding together and renaming
res <- cbind(N, p.value, estimate, estimate.1, conf.low, conf.high)
res <- tibble::rownames_to_column(as.data.frame(res), "var")
res <- res %>% transmute(var,
                         N,
                         p.value = res[,"Pr(>|z|)"],
                         estimate = coef,
                         estimate.1 = exp(coef),
                         conf.low = res[,"lower .95"],
                         conf.high = res[,"upper .95"],
                         ci = paste0("(", 
                                     conf.low, "-", conf.high, ")"))
names(res)<-c( "var","N","p.value", "estimate", "estimate.1", 
               "conf.low", "conf.high", "ci")

# Correct for group size
datalist <- list()
covariates <- c("SPP1","Gender","Onset","Sampling_Age")
data  <- .get_data(res.clu, data = allcoxdf)
for (i in 1:length(covariates)){
  var <- covariates[i]
  adf <- as.data.frame(table(data[, var]))
  datalist[[i]] <- cbind(var = var, adf, pos = 1:nrow(adf))
}
total <- do.call(rbind,datalist)
n <- subset(total, pos == 2) %>% 
  select(Freq)
res$N[2:3] <- n[2:3,]
res$N <- paste0("(N=",res$N,")")

# Add stars (p.value) to dataframe
res$stars <- paste0(round(res$p.value, 3), " ",
                    ifelse(res$p.value < 0.05, "*",""),
                    ifelse(res$p.value < 0.01, "*",""),
                    ifelse(res$p.value < 0.001, "*",""))
res$stars[which(res$p.value < 0.001)] = "<0.001 ***"

# Change datatypes
res$N <- as.character(res$N)
res$p.value <- as.character(res$p.value)
res$estimate <- as.numeric(as.character(res$estimate))
res$estimate.1 <- as.character(round(res$estimate.1,2))
res$conf.low <- as.numeric(as.character(res$conf.low))
res$conf.high <- as.numeric(as.character(res$conf.high))
res$ci <- as.character(res$ci)

# Flip order
res <- res[nrow(res):1, ]

# Setting plot sizes and scales
rangeb <- range(log(res$conf.low), log(res$conf.high), na.rm = TRUE)
breaks <- axisTicks(rangeb/2, log = TRUE, nint = 7)
breaks <- breaks[1:6]
rangeplot <- rangeb
# Make plot twice as wide as needed to create space for annotations
rangeplot[1] <- rangeplot[1] - diff(rangeb)
# Increase white space on right for p-vals
rangeplot[2] <- rangeplot[2] + .15 * diff(rangeb)
width <- diff(rangeplot)
# Y-coordinates for labels
y_variable <- rangeplot[1] +  cpositions[1] * width
y_nlevel <- rangeplot[1]  +  cpositions[2] * width
y_cistring <- rangeplot[1]  +  cpositions[3] * width
y_stars <- rangeb[2]
x_annotate <- seq_len(nrow(res))

# Change column names for readability
res$var[res$var == "OnsetBulbar"] <- "Onset Bulbar"
res$var[res$var == "GenderM"] <- "Gender Male"
res$var[res$var == "SPP1"] <- "SPP1"
res$var[res$var == "Sampling_Age"] <- "Sampling Age"

# Build forest plot
p <- ggplot(res, aes(seq_along(var), exp(estimate))) +
  geom_rect(aes(xmin = seq_along(var) - .5, xmax = seq_along(var) + .5,
                ymin = exp(rangeplot[1]), ymax = exp(rangeplot[2]),
                fill = ordered(seq_along(var) %% 2 + 1))) +
  scale_fill_manual(values = c("#FFFFFF33", "#00000033"), guide = "none") +
  geom_point(pch = 15, size = 4) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.15) +
  geom_hline(yintercept = 1, linetype = 3) +
  coord_flip(ylim = exp(rangeplot)) +
  ggtitle(main) +
  scale_y_log10(
    name = "",
    labels = sprintf("%g", breaks),
    expand = c(0.02, 0.02),
    breaks = breaks) +
  theme_minimal() +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position = "none",
        panel.border=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  xlab("") +
  annotate(geom = "text", x = x_annotate, y = exp(y_variable),
           label = res$var, fontface = "bold", hjust = 0,
           size = annot_size_mm) +
  annotate(geom = "text", x = x_annotate, y = exp(y_nlevel),
           label = res$N, fontface = "italic", hjust = 0,
           size = annot_size_mm) +
  annotate(geom = "text", x = x_annotate, y = exp(y_cistring),
           label = res$estimate.1, size = annot_size_mm) +
  annotate(geom = "text", x = x_annotate, y = exp(y_cistring),
           label = res$ci, size = annot_size_mm,
           vjust = 2,  fontface = "italic") +
  annotate(geom = "text", x = x_annotate, y = exp(y_stars),
           label = res$p.value, size = annot_size_mm,
           hjust = -0.2,  fontface = "italic") +
  annotate(geom = "text", x = 0.5, y = exp(y_variable),
           label = paste0("Multivariate coxmodels, includes corrections for sampling delay and cohort identity."),
           size = annot_size_mm, hjust = 0, vjust = 1.2,  fontface = "italic")

# Switch off clipping for p-vals, bottom annotation
gt <- ggplot_gtable(ggplot_build(p))
gt$layout$clip[gt$layout$name == "panel"] <- "off"
grid.draw(gt)

# Print model
pdf(file = "output_R/2020_Multiv_Cox_Cont.pdf", width = 10, height = 7)
ggpubr::as_ggplot(gt)
dev.off()

# Schoenfeld residuals test
cat("Multivariate Cox model finished, proportionality for covariates are tested. See below.\n")
cox.zph(res.clu)

### Continuous Multivariate Cox with Cutoff --------------------------------

# Title for forest plot
main <-  "Hazard ratio"

# Model creation for the covariates
res.clu <- coxph(survobj ~ SPP1_Level + COL6A1_Level + NEFL_Level +
                   Gender + Onset + Sampling_Age + strata(Sampling_delay) +
                   cluster(cohort), data = allcoxdf)

# Extract values from models
x <- summary(res.clu)
N <- x$n
coef <- as.data.frame(x$coefficients)
conf <- as.data.frame(x$conf.int)
p.value <- signif(coef["Pr(>|z|)"], digits = 3)
estimate <-signif(coef[1], digits=3);#coeficient beta
estimate.1 <-signif(coef[2], digits=3);#exp(beta)
conf.low <- signif(conf[3], 3)
conf.high <- signif(conf[4],3)
# Binding together and renaming
res<-cbind(N, p.value, estimate, estimate.1, conf.low, conf.high)
res <- tibble::rownames_to_column(as.data.frame(res), "var")
res <- res %>% transmute(var,
                         N,
                         p.value = res[,"Pr(>|z|)"],
                         estimate = coef,
                         estimate.1 = exp(coef),
                         conf.low = res[,"lower .95"],
                         conf.high = res[,"upper .95"],
                         ci = paste0("(", 
                                     conf.low, "-", conf.high, ")"))
names(res)<-c( "var","N","p.value", "estimate", "estimate.1", 
               "conf.low", "conf.high", "ci")

# Correct for group size
datalist <- list()
covariates <- c("SPP1_Level", "COL6A1_Level", "NEFL_Level", "Gender", "Onset")
data  <- .get_data(res.clu, data = allcoxdf)
for (i in 1:length(covariates)){
  var <- covariates[i]
  adf <- as.data.frame(table(data[, var]))
  datalist[[i]] <- cbind(var = var, adf, pos = 1:nrow(adf))
}
total <- do.call(rbind,datalist)
n <- subset(total, pos == 2) %>% 
  select(Freq)
res$N[1:5] <- n[1:5,]
res$N <- paste0("(N=",res$N,")")

# Add stars (p.value) to dataframe
res$stars <- paste0(round(res$p.value, 3), " ",
                    ifelse(res$p.value < 0.05, "*",""),
                    ifelse(res$p.value < 0.01, "*",""),
                    ifelse(res$p.value < 0.001, "*",""))
res$stars[which(res$p.value < 0.001)] = "<0.001 ***"

# Change datatypes
res$N <- as.character(res$N)
res$p.value <- as.character(res$p.value)
res$estimate <- as.numeric(as.character(res$estimate))
res$estimate.1 <- as.character(round(res$estimate.1,2))
res$conf.low <- as.numeric(as.character(res$conf.low))
res$conf.high <- as.numeric(as.character(res$conf.high))
res$ci <- as.character(res$ci)

# Flip order
res <- res[nrow(res):1, ]

# Setting plot sizes and scales
rangeb <- range(log(res$conf.low), log(res$conf.high), na.rm = TRUE)
breaks <- axisTicks(rangeb/2, log = TRUE, nint = 7)
rangeplot <- rangeb
# Make plot twice as wide as needed to create space for annotations
rangeplot[1] <- rangeplot[1] - diff(rangeb)
# Increase white space on right for p-vals
rangeplot[2] <- rangeplot[2] + .15 * diff(rangeb)
width <- diff(rangeplot)
# Y-coordinates for labels
y_variable <- rangeplot[1] +  cpositions[1] * width
y_nlevel <- rangeplot[1]  +  cpositions[2] * width
y_cistring <- rangeplot[1]  +  cpositions[3] * width
y_stars <- rangeb[2]
x_annotate <- seq_len(nrow(res))

# Change column names for readability
res$var[res$var == "OnsetBulbar"] <- "Onset Bulbar"
res$var[res$var == "GenderM"] <- "Gender Male"
res$var[res$var == "SPP1_LevelUpper"] <- "High level SPP1"
res$var[res$var == "COL6A1_LevelUpper"] <- "High level COL6A1"
res$var[res$var == "NEFL_LevelUpper"] <- "High level NEFL"
res$var[res$var == "Sampling_Age"] <- "Sampling Age"

# Build forest plot
p <- ggplot(res, aes(seq_along(var), exp(estimate))) +
  geom_rect(aes(xmin = seq_along(var) - .5, xmax = seq_along(var) + .5,
                ymin = exp(rangeplot[1]), ymax = exp(rangeplot[2]),
                fill = ordered(seq_along(var) %% 2 + 1))) +
  scale_fill_manual(values = c("#FFFFFF33", "#00000033"), guide = "none") +
  geom_point(pch = 15, size = 4) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.15) +
  geom_hline(yintercept = 1, linetype = 3) +
  coord_flip(ylim = exp(rangeplot)) +
  ggtitle(main) +
  scale_y_log10(
    name = "",
    labels = sprintf("%g", breaks),
    expand = c(0.02, 0.02),
    breaks = breaks) +
  theme_minimal() +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position = "none",
        panel.border=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  xlab("") +
  annotate(geom = "text", x = x_annotate, y = exp(y_variable),
           label = res$var, fontface = "bold", hjust = 0,
           size = annot_size_mm) +
  annotate(geom = "text", x = x_annotate, y = exp(y_nlevel),
           label = res$N, fontface = "italic", hjust = 0,
           size = annot_size_mm) +
  annotate(geom = "text", x = x_annotate, y = exp(y_cistring),
           label = res$estimate.1, size = annot_size_mm) +
  annotate(geom = "text", x = x_annotate, y = exp(y_cistring),
           label = res$ci, size = annot_size_mm,
           vjust = 2,  fontface = "italic") +
  annotate(geom = "text", x = x_annotate, y = exp(y_stars),
           label = res$p.value, size = annot_size_mm,
           hjust = -0.2,  fontface = "italic") +
  annotate(geom = "text", x = 0.5, y = exp(y_variable),
           label = paste0("Multivariate coxmodels, includes corrections for sampling delay and cohort identity."),
           size = annot_size_mm, hjust = 0, vjust = 1.2,  fontface = "italic")

# Switch off clipping for p-vals, bottom annotation
gt <- ggplot_gtable(ggplot_build(p))
gt$layout$clip[gt$layout$name == "panel"] <- "off"
grid.draw(gt)

# Print model
pdf(file = "output_R/2020_Multiv_Cox_Thresholds.pdf", width = 10, height = 7)
ggpubr::as_ggplot(gt)
dev.off()

# Schoenfeld residuals test
cat("Multivariate Cox model (Threshold) finished, proportionality for covariates are tested. See below.\n")
cox.zph(res.clu)
