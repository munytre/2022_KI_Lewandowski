###--------------------------------------------------------------------------###
### Project: ApoE in ALS                                                       -
### Purpose: Readout Human Anterior Horn                                       -
### Authors: QC Lin                                                            -
### Last Edit: 2022-08-18, QC Lin                                              -
### Contact: q.c.lin@students.uu.nl or sebastian.lewandowski@ki.se             -
###--------------------------------------------------------------------------###

### Dependencies and general settings ------------------------------------------

# Load libraries
library(tidyverse)
library(viridis)

### Functions ------------------------------------------------------------------
load_data <- function(path_data,n_celltypes){
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
  # Generate metadata
  metadata <- as_tibble(cell_ratio$names) %>%
    `colnames<-`("names")
  metadata$condition <- metadata %>%
    separate(names,
             into = c("location", "condition", "other")) %>%
    pull("condition")
  # Merge metadata and proportion data
  cell_ratio_tidy_joined <- left_join(cell_ratio_tidy,
                                      metadata,
                                      by = "names")
  return(cell_ratio_tidy_joined)
}

simple_plot <- function(data){
  ggplot(data,
         aes(x = celltype,
             y = ratio,
             fill = condition)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    labs(x = "Celltypes",
         y = "Estimated proportion",
         fill = "Condition")
}

simple_plot_y_sqrt <- function(data){
  ggplot(data,
         aes(x = celltype,
             y = ratio,
             fill = condition)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    labs(x = "Celltypes",
         y = "Estimated proportion",
         fill = "Condition") +
    scale_y_sqrt()
}

### Load data ------------------------------------------------------------------
# Mean single cell reference
AH_mean_data <- load_data("20220531_Oscars_tool/READOUT/READOUT_HUMAN_AH.csv",
                          10)
# Individual single cell reference
AH_7238_data <- load_data("20220531_Oscars_tool/READOUT/READOUT_HUMAN_AH_7238.csv",
                          11)
### Plots ----------------------------------------------------------------------
simple_plot(AH_mean_data)
simple_plot(AH_7238_data)
simple_plot_y_sqrt(AH_7238_data)
