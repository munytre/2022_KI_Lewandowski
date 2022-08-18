###--------------------------------------------------------------------------###
### Project: ApoE in ALS                                                       -
### Purpose: Readout Mouse                                                     -
### Authors: QC Lin                                                            -
### Last Edit: 2022-08-18, QC Lin                                              -
### Contact: q.c.lin@students.uu.nl or sebastian.lewandowski@ki.se             -
###--------------------------------------------------------------------------###

### Dependencies and general settings ------------------------------------------

# Load libraries
library(tidyverse)
library(viridis)

### Functions ------------------------------------------------------------------
load_data <- function(path_data,path_metadata,n_celltypes){
  # Load in proportion data
  READOUT_all <- read.csv(path_data,
                          row.names=1)
  # Prepare data in long format
  cell_ratio <- READOUT_all
  cell_ratio$names <- rownames(cell_ratio)
  cell_ratio_tidy <- pivot_longer(cell_ratio,
                                  cols = 1:all_of(n_celltypes),
                                  names_to ="celltype",
                                  values_to = "ratio")
  # Prepare metadata
  load(path_metadata)
  sample_age     = sample_data[,"Age"]
  sample_genotype  = sample_data[,"Genotype"]
  metadata     = data.frame(age=sample_age,genotype=sample_genotype,names=sample_data$ID_REF)
  rownames(metadata) = metadata$names
  # Merge metadata and proportion data
  cell_ratio_tidy_joined <- left_join(cell_ratio_tidy,
                                      metadata,
                                      by = "names")
  return(cell_ratio_tidy_joined)
}

facet_plot_age <- function(data){
  ggplot(data,
         aes(x = celltype,
             y = ratio,
             fill = genotype)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    labs(x = "Celltypes",
         y = "Estimated proportion",
         fill = "Genotype") +
    facet_grid(rows = vars(factor(age, levels = c("28 days",
                                                  "42 days",
                                                  "56 days",
                                                  "70 days",
                                                  "98 days",
                                                  "112 days",
                                                  "126 days"))))
}

facet_plot_age_y_sqrt <- function(data){
  ggplot(data,
         aes(x = celltype,
             y = ratio,
             fill = genotype)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    labs(x = "Celltypes",
         y = "Estimated proportion",
         fill = "Genotype") +
    facet_grid(rows = vars(factor(age, levels = c("28 days",
                                                  "42 days",
                                                  "56 days",
                                                  "70 days",
                                                  "98 days",
                                                  "112 days",
                                                  "126 days")))) +
    scale_y_sqrt()
}

facet_plot_genotype <- function(data){
  ggplot(data,
         aes(x = factor(age, levels = c("28 days",
                                        "42 days",
                                        "56 days",
                                        "70 days",
                                        "98 days",
                                        "112 days",
                                        "126 days")),
             y = ratio,
             fill = genotype)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          strip.text.x = element_text(size = 5)) +
    labs(x = "Celltypes",
         y = "Estimated proportion",
         fill = "Genotype") +
    facet_grid(genotype ~ celltype)
}

facet_plot_genotype_y_sqrt <- function(data){
  ggplot(data,
         aes(x = factor(age, levels = c("28 days",
                                        "42 days",
                                        "56 days",
                                        "70 days",
                                        "98 days",
                                        "112 days",
                                        "126 days")),
             y = ratio,
             fill = genotype)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          strip.text.x = element_text(size = 5)) +
    labs(x = "Celltypes",
         y = "Estimated proportion",
         fill = "Genotype") +
    facet_grid(genotype ~ celltype) +
    scale_y_sqrt()
}

### Load data ------------------------------------------------------------------
# Mean single cell reference
Mouse_mean_data <- load_data("20220531_Oscars_tool/READOUT/READOUT_all.csv",
                             "Data/sample_data.Rda",
                             11)
# Individual single cell reference
Mouse_7238_data <- load_data("20220531_Oscars_tool/READOUT/READOUT_Mouse_7238.csv",
                             "Data/sample_data.Rda",
                             11)

### Plots ----------------------------------------------------------------------
facet_plot_age(Mouse_mean_data)
facet_plot_genotype(Mouse_mean_data)
facet_plot_age(Mouse_7238_data)
facet_plot_genotype(Mouse_7238_data)
facet_plot_age_y_sqrt(Mouse_7238_data)
facet_plot_genotype_y_sqrt(Mouse_7238_data)

### Old code -------------------------------------------------------------------
# READOUT_all <- read.csv("20220531_Oscars_tool/READOUT/READOUT_all.csv",
#                         row.names=1)
# # Prepare data in long format
# cell_ratio <- READOUT_all
# cell_ratio$names <- rownames(cell_ratio)
# cell_ratio_tidy <- pivot_longer(cell_ratio,
#                                 cols = 1:11,
#                                 names_to ="celltype",
#                                 values_to = "ratio")
# 
# load("Data/sample_data.Rda")
# sample_age     = sample_data[,"Age"]
# sample_genotype  = sample_data[,"Genotype"]
# metadata     = data.frame(age=sample_age,genotype=sample_genotype,names=sample_data$ID_REF)
# rownames(metadata) = metadata$names
# 
# cell_ratio_tidy_joined <- left_join(cell_ratio_tidy,
#                                     metadata,
#                                     by = "names")