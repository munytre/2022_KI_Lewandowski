###--------------------------------------------------------------------------###
### Project: ApoE in ALS                                                       -
### Purpose: Cutpoint analyses script for combined cohorts                     -
### Authors: FS Sanders, QC Lin                                                -
### Last Edit: 2022-02-09, QC Lin                                              -
### Contact: q.c.lin@students.uu.nl or sebastian.lewandowski@ki.se             -
###--------------------------------------------------------------------------###

### Dependencies and general settings ------------------------------------------

# Load libraries
library(tidyverse)
library(gridExtra)
library(survminer)
library(ggbeeswarm)
library(survival)

# Set working directory and seed for reproducibility
setwd(".")
set.seed(747)

# Create output directory for pdfs
ifelse(!dir.exists(file.path(getwd(), "output_R")),
       dir.create(file.path(getwd(), "output_R")),
       FALSE)
ifelse(!dir.exists(file.path(getwd(), "output_R_nonstrat")),
       dir.create(file.path(getwd(), "output_R_nonstrat")),
       FALSE)
ifelse(!dir.exists(file.path(getwd(), "output_R_e33")),
       dir.create(file.path(getwd(), "output_R_e33")),
       FALSE)
ifelse(!dir.exists(file.path(getwd(), "output_R_e34")),
       dir.create(file.path(getwd(), "output_R_e34")),
       FALSE)

# Run helper function created within the package (group separation)
.dichotomize <- function(x, cutpoint, labels = c("low", "high")){
  grps <- x
  grps[x <= cutpoint] = labels[1]
  grps[x > cutpoint] = labels[2]
  
  grps
}

### Cutpoint analyses ----------------------------------------------------------

# Function for cutpoint analysis
Cutpoint_analyses <- function(df1,df2,protein,output_folder,stratification){
  # Cutpoint and density
  print(protein)
  single.cut <- surv_cutpoint(df1,
                              time = "survival_months",
                              event = "died",
                              variables = protein,
                              minprop = 0.1)
  max_stat <- single.cut[[protein]]
  cutpoint <- as.numeric(max_stat$estimate)
  print(cutpoint)
  # Create the data frame to plot
  p_data <- data.frame(
    stats = max_stat$stats,
    cuts = max_stat$cuts,
    grps = .dichotomize(max_stat$cuts, cutpoint)
  )
  # Variables for plotting
  vline_df <- data.frame(x1 = cutpoint, x2 = cutpoint,
                         y1 = 0, y2 = max(p_data$stats))
  down <- floor(min(df2[[protein]]))
  top <- ceiling(max(df2[[protein]]))
  limits <- c(down,top)
  posx <- top-0.3
  cutpoint_label <- paste0("Cutpoint: ", round(cutpoint,2))
  x1 <- y1 <- x2 <- y2 <- NULL
  posmidy <- 1
  
  # Cutpoint plot
  m <- ggplot(data = p_data, mapping=aes_string("cuts", "stats")) +
    geom_point(aes_string(color = "grps"), shape = 19, size = 0.5) +
    geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2),
                 data = vline_df, linetype = "dashed", size = 0.5) +
    ggplot2::annotate("text", x = posx, y=posmidy,
                      label = cutpoint_label, size = 3) +
    labs(y = "Standardized Log-Rank Statistic", x = paste0(" MFI")) +
    scale_x_continuous(limits = limits)+ coord_flip()
  # Apply cutpoints
  controls <- subset(df2, df2$class == "Control" | df2$class == "Neuro control")
  controls$class <- "Controls"
  controls$level <- ifelse(controls[[protein]] <= cutpoint, "Low", "High")
  df1$level <- ifelse(df1[[protein]] <= cutpoint, "Low", "High")
  df1$class <- "ALS"
  
  # Disease boxplot
  a <- ggboxplot(df1,
                 x= "class",
                 y= protein,
                 ylab = paste0(protein," MFI"),
                 xlab = "ALS Patients",
                 outlier.shape = " ") +
    geom_hline(yintercept = cutpoint, linetype = 2) +
    geom_beeswarm(aes(color = df1[[protein]] < cutpoint), show.legend = FALSE) +
    scale_y_continuous(limits = limits)
  # Controls boxplot
  b <- ggboxplot(controls,
                 x= "cohort",
                 y= protein,
                 ylab = paste0(protein," MFI"),
                 xlab = "Controls",
                 outlier.shape = " ") +
    geom_hline(yintercept = cutpoint, linetype = 2) +
    geom_beeswarm(aes(color = controls[[protein]] < cutpoint), show.legend = FALSE) +
    scale_y_continuous(limits = limits)
  # Density plot
  e <- ggplot(data = df1, aes(x=df1[[protein]], fill=level))+
    geom_density( alpha=0.3) +
    geom_vline(xintercept = cutpoint, linetype = 2) +
    scale_x_continuous(limits = limits) +
    labs(x= paste0(protein," MFI"), y="Density") +
    coord_flip()
  # Custom legend showing number of censored can be checked with the following:
  namedf <- "Combined"
  t <- table(df1$level)
  high <- paste0("High: ",t[1])
  low <- paste0("Low: ",t[2])
  highlow <- c(high,low)
  fit <- survfit(Surv(survival_months, died) ~ level, df1)
  km <- ggsurvplot(fit,
                   data = df1,
                   palette = c("#ff3333","#0099ff" ),
                   conf.int = FALSE,
                   legend.title = paste0(protein," KM ",namedf),
                   legend=c(0.8,0.8),
                   pval = TRUE,
                   legend.labs = highlow,
                   pval.size = 4)
  kmp <- km$plot
  risk <- km$table
  lay <- rbind(c(1,2,3,4,5,5),
               c(1,2,3,4,5,5))
  # c(NA,NA,NA,6,6,6))
  plotlist <<- list(b,a,m,e,kmp)
  
  # Show survival ratio for level=high (dead/alive)
  grp_high <- subset(df1, df1$level == "High")
  n_high <-table(grp_high$died)
  alive_high <- n_high[1]
  dead_high <- n_high[2]
  cat(paste0("High\n",
             "Alive: ", alive_high,
             ",Dead: ", dead_high, "\n"))
  cat("-----\n")
  # Show survival ratio for level=low (dead/alive)
  grp_low <- subset(df1, df1$level == "Low")
  n_low <- table(grp_low$died)
  alive_low <- n_low[1]
  dead_low <- n_low[2]
  cat(paste0("Low\n",
             "Alive: ", alive_low,
             ",Dead: ", dead_low, "\n"))
  cat("-----\n")
  # Show survival ratio for all (dead/alive)
  n <- table(df1$died)
  alive <- n[1]
  dead <- n[2]
  cat(paste0("Total ALS: ", alive + dead, "\n",
             "Alive: ", alive,
             ",Dead: ", dead, "\n"))
  cat("-----\n")
  # Show amount of controls
  cat(paste0("Total controls: ", nrow(controls), "\n"))
  
  pdf(file = paste0(output_folder,protein,stratification,"_cutpoint_",namedf,".pdf"),
      width = 16,height = 6)
  gridExtra::grid.arrange(grobs = plotlist, layout_matrix = lay)
  dev.off()
  pval <- surv_pvalue(fit, data = df1)
  return(pval[,2])
}


### Load data ------------------------------------------------------------------

# Load from rds objects
bead_info_utrecht <- readRDS(file = "generated_data_R/bead_info_utrecht.rds")
full_beads_utrecht <- readRDS(file = "generated_data_R/als.ma.apoe.rds")
als_beads_utrecht <- readRDS(file = "generated_data_R/apoeALS.rds")

### Genotype separation --------------------------------------------------------

genotypeALS <- subset(als_beads_utrecht,
                      als_beads_utrecht$apoe == "ALS Apo-e3e3" | als_beads_utrecht$apoe == "ALS Apo-e3e4" )
genotypeALSe33 <- subset(als_beads_utrecht,
                        als_beads_utrecht$apoe == "ALS Apo-e3e3")
genotypeALSe34 <- subset(als_beads_utrecht,
                         als_beads_utrecht$apoe == "ALS Apo-e3e4")


### Results --------------------------------------------------------------------

# All beads
bead_info_utrecht$pval_nonstrat <- map_dbl(as.character(bead_info_utrecht$Bead_name),
                                           ~Cutpoint_analyses(als_beads_utrecht,
                                                              full_beads_utrecht,
                                                              .x,
                                                              "output_R_nonstrat/",
                                                              "_nonstrat"))
bead_info_utrecht$pval_e33 <- map_dbl(as.character(bead_info_utrecht$Bead_name),
                                           ~Cutpoint_analyses(genotypeALSe33,
                                                              full_beads_utrecht,
                                                              .x,
                                                              "output_R_e33/",
                                                              "_e33"))
bead_info_utrecht$pval_e34 <- map_dbl(as.character(bead_info_utrecht$Bead_name),
                                           ~Cutpoint_analyses(genotypeALSe34,
                                                              full_beads_utrecht,
                                                              .x,
                                                              "output_R_e34/",
                                                              "_e34"))

# Write pvals to csv
write_csv(bead_info_utrecht,
          file = "output_R/bead_info_utrecht_pvals.csv")

### P-Val filter ----------------------------------------------------------------

# Make new variable to work with
pvals <- bead_info_utrecht

# Loose filtering (nonstrat no sign, e33 or e34 sign)
pval_interest <- pvals %>%
  filter(pval_nonstrat > 0.05) %>%
  filter(pval_e33 <= 0.05 | pval_e34 <= 0.05)

# Strict filtering (nonstrat no sign, e33 and e34 sign)
pval_strict <- pvals %>%
  filter(pval_nonstrat > 0.05) %>%
  filter(pval_e33 <= 0.05 & pval_e34 <= 0.05)

# Strict filtering (nonstrat and e34 no sign, e33 sign)
pval_strict_e33 <- pvals %>%
  filter(pval_nonstrat > 0.05) %>%
  filter(pval_e33 <= 0.05 & pval_e34 > 0.05)

# Strict filtering (nonstrat and e33 no sign, e34 sign)
pval_strict_e34 <- pvals %>%
  filter(pval_nonstrat > 0.05) %>%
  filter(pval_e33 > 0.05 & pval_e34 <= 0.05)

# Write pvals to csv
write_csv(pval_interest,
          file = "output_R/bead_info_utrecht_pvals_loose.csv")
# Write pvals to csv
write_csv(pval_strict,
          file = "output_R/bead_info_utrecht_pvals_strict.csv")
# Write pvals to csv
write_csv(pval_strict_e33,
          file = "output_R/bead_info_utrecht_pvals_strict_e33.csv")
# Write pvals to csv
write_csv(pval_strict_e34,
          file = "output_R/bead_info_utrecht_pvals_strict_e34.csv")


### Folkert's analyses (OLD) ---------------------------------------------------

# # Non-cohort specific
# Cutpoint_analyses(combinedALS,
#                   cohorts,
#                   "RawB.53")
# Cutpoint_analyses(combinedALS,
#                   cohorts,
#                   "RawB.123")
# 
# # Cohort specific
# cohortB <- cohorts %>%
#   filter(cohort %in% "Leuven")
# combinedALSB <- combinedALS %>%
#   filter(cohort %in% "Leuven")
# Cutpoint_analyses(combinedALSB,
#                   cohortB,
#                   "RawB.53")
