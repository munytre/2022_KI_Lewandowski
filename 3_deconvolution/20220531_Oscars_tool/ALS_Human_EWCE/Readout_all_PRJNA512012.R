###--------------------------------------------------------------------------###
### Project: ApoE in ALS                                                       -
### Purpose: Readout PRJNA512012                                               -
### Authors: QC Lin                                                            -
### Last Edit: 2022-08-18, QC Lin                                              -
### Contact: q.c.lin@students.uu.nl or sebastian.lewandowski@ki.se             -
###--------------------------------------------------------------------------###

### Dependencies and general settings ------------------------------------------

# Load libraries
library(tidyverse)
library(viridis)
library(ggpubr)
library(rstatix)

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
  # Load metadata
  metadata_PRJNA512012_2 <- read.delim(path_metadata)
  metadata_PRJNA512012 <- read.delim(path_metadata)[,c(1,4,13)] %>%
    `colnames<-`(c("names",
                   "Gender",
                   "APOE"))
  # Merge metadata and proportion data
  cell_ratio_tidy_joined <- left_join(cell_ratio_tidy,
                                      metadata_PRJNA512012,
                                      by = "names")
  cell_ratio_tidy_joined$APOE <- as.factor(cell_ratio_tidy_joined$APOE)
  cell_ratio_tidy_joined_33_34 <- cell_ratio_tidy_joined %>%
    filter(APOE %in% c(33,34))
  return(cell_ratio_tidy_joined_33_34)
}

simple_plot <- function(data){
  ggplot(data,
         aes(x = celltype,
             y = ratio,
             fill = APOE)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    labs(x = "Celltypes",
         y = "Estimated proportion") +
    scale_fill_manual(name = "APOE",
                      values = c("orange","brown"))
}

simple_plot_y_sqrt <- function(data){
  ggplot(data,
         aes(x = celltype,
             y = ratio,
             fill = APOE)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    labs(x = "Celltypes",
         y = "Estimated proportion") +
    scale_fill_manual(name = "APOE",
                      values = c("orange","brown")) +
    scale_y_sqrt()
}

plot_gendersplit <- function(data,output_stats,output_normality_test){
  # Prepare order of x-axis in plot
  data$celltype <- factor(data$celltype,
                          levels = c("Astrocytes",
                                     "Interneurons",
                                     "Microglia",
                                     "Oligodendrocyte.Precursor",
                                     "Oligodendrocytes",
                                     "Pericytes",
                                     "Pyramidal.Neurons",
                                     "Vascular.and.Leptomeningeal.Cells",
                                     "Vascular.Endothelial",
                                     "Vascular.Smooth.Muscle.Cell",
                                     "X.none."))
  # Designing plot
  bxp <- ggboxplot(data,
                   x = "celltype",
                   y = "ratio",
                   fill = "APOE") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    labs(x = "Celltypes",
         y = "Estimated proportion") +
    scale_fill_manual(name = "APOE",
                      values = c("orange","brown")) +
    facet_grid(rows = vars(Gender))
  # Statistics
  normality_test <- data %>%
    group_by(Gender, celltype) %>%
    shapiro_test(ratio)
  stat.test <- data %>%
    group_by(Gender, celltype) %>%
    wilcox_test(ratio ~ APOE) %>%
    add_significance()
  # Write statistics to disk
  write_csv(normality_test,
            file = output_normality_test)
  write_csv(stat.test,
            file = output_stats)
  # Set coordinates for significance bars
  stat.test <- stat.test %>%
    add_xy_position(x = "celltype", dodge = 0.8)
  # Plot with significance bars
  plot <- bxp + 
    stat_pvalue_manual(
      stat.test,
      label = "p.signif",
      tip.length = 0.01,
      hide.ns = TRUE) +
    scale_y_continuous(expand = expansion(mult = c(0.01, 0.1)))
  return(plot)
}

### Load data ------------------------------------------------------------------
# Mean single cell reference - MC
PRJNA512012_mean_data <- load_data("20220531_Oscars_tool/READOUT/READOUT_Human_all_ALS_PRJNA512012.csv",
                                   "20220531_Oscars_tool/ALS_Human_EWCE/Data/Metadata_PRJNA512012_ALS_MC.csv",
                                   10)
# Individual single cell reference - MC
PRJNA512012_7238_data <- load_data("20220531_Oscars_tool/READOUT/READOUT_HUMAN_PRJNA512012_7238.csv",
                                   "20220531_Oscars_tool/ALS_Human_EWCE/Data/Metadata_PRJNA512012_ALS_MC.csv",
                                   11)
# Mean single cell reference - FC
PRJNA512012_mean_data_FC <- load_data("20220531_Oscars_tool/READOUT/READOUT_PRJNA512012_FC.csv",
                                      "20220531_Oscars_tool/ALS_Human_EWCE/Data/Metadata_PRJNA512012_ALS_FC.csv",
                                      10)
# Individual single cell reference - FC
PRJNA512012_7238_data_FC <- load_data("20220531_Oscars_tool/READOUT/READOUT_PRJNA512012_FC_7238.csv",
                                      "20220531_Oscars_tool/ALS_Human_EWCE/Data/Metadata_PRJNA512012_ALS_FC.csv",
                                      11)

### Plots ----------------------------------------------------------------------
# MC
simple_plot(PRJNA512012_mean_data)
simple_plot(PRJNA512012_7238_data)
simple_plot_y_sqrt(PRJNA512012_7238_data)
plot_gendersplit(PRJNA512012_mean_data,
                 "20220531_Oscars_tool/ALS_Human_EWCE/Output_READOUT_Human/Stats_PRJNA512012_Celltype_AvgSS_GenderSplit.txt",
                 "20220531_Oscars_tool/ALS_Human_EWCE/Output_READOUT_Human/Normality_PRJNA512012_Celltype_AvgSS_GenderSplit.txt")
plot_gendersplit(PRJNA512012_7238_data,
                 "20220531_Oscars_tool/ALS_Human_EWCE/Output_READOUT_Human/Stats_PRJNA512012_Celltype_7238SS_GenderSplit.txt",
                 "20220531_Oscars_tool/ALS_Human_EWCE/Output_READOUT_Human/Normality_PRJNA512012_Celltype_7238SS_GenderSplit.txt")

# FC
simple_plot(PRJNA512012_mean_data_FC)
simple_plot(PRJNA512012_7238_data_FC)
simple_plot_y_sqrt(PRJNA512012_7238_data_FC)
plot_gendersplit(PRJNA512012_mean_data_FC,
                 "20220531_Oscars_tool/ALS_Human_EWCE/Output_READOUT_Human/Stats_PRJNA512012_Celltype_AvgSS_GenderSplit_FC.txt",
                 "20220531_Oscars_tool/ALS_Human_EWCE/Output_READOUT_Human/Normality_PRJNA512012_Celltype_AvgSS_GenderSplit_FC.txt")
plot_gendersplit(PRJNA512012_7238_data_FC,
                 "20220531_Oscars_tool/ALS_Human_EWCE/Output_READOUT_Human/Stats_PRJNA512012_Celltype_7238SS_GenderSplit_FC.txt",
                 "20220531_Oscars_tool/ALS_Human_EWCE/Output_READOUT_Human/Normality_PRJNA512012_Celltype_7238SS_GenderSplit_FC.txt")

### Age of Onset comparison ----------------------------------------------------
# MC
meta_MC <- PRJNA512012_mean_data[c(1,4,5,6)] %>%
  unique()
ggplot(meta_MC,
       aes(x = Age_Onset, color = Gender, linetype = APOE)) +
  geom_density()
meta_MC_comb <- meta_MC %>%
  unite(group, c(Gender,APOE))
ggplot(meta_MC_comb,
       aes(x = Age_Onset, after_stat(count), fill = group, color = group)) +
  geom_density(position = "fill") +
  scale_fill_manual(breaks = c("Female_33","Female_34","Male_33","Male_34"),
                    values = c("#F8766D","red","#00BFC4","blue")) +
  scale_colour_manual(breaks = c("Female_33","Female_34","Male_33","Male_34"),
                      values = c("#F8766D","red","#00BFC4","blue"))

# FC
meta_FC <- PRJNA512012_mean_data_FC[c(1,4,5,6)] %>%
  unique()
ggplot(meta_FC,
       aes(x = Age_Onset, color = Gender, linetype = APOE)) +
  geom_density()
ggplot(meta_FC,
       aes(x = Age_Onset, color = Gender, linetype = APOE)) +
  geom_freqpoly()
meta_FC_comb <- meta_FC %>%
  unite(group, c(Gender,APOE))
ggplot(meta_FC_comb,
       aes(x = Age_Onset, after_stat(count), fill = group, color = group)) +
  geom_density(position = "fill") +
  scale_fill_manual(breaks = c("Female_33","Female_34","Male_33","Male_34"),
                    values = c("#F8766D","red","#00BFC4","blue")) +
  scale_colour_manual(breaks = c("Female_33","Female_34","Male_33","Male_34"),
                      values = c("#F8766D","red","#00BFC4","blue"))

### Old code -------------------------------------------------------------------
# # # READOUT_all <- read.csv("20220531_Oscars_tool/READOUT/READOUT_Human_all_ALS_PRJNA512012.csv", row.names=1)
# # READOUT_all <- read.csv("20220531_Oscars_tool/READOUT/READOUT_HUMAN_PRJNA512012_7238.csv", row.names=1)
# # # 
# # cell_ratio <- READOUT_all
# # cell_ratio$names <- rownames(cell_ratio)
# # # 
# # # # cell_ratio_tidy <- pivot_longer(cell_ratio, cols = 1:10, names_to ="celltype", values_to = "ratio")
# # cell_ratio_tidy <- pivot_longer(cell_ratio, cols = 1:11, names_to ="celltype", values_to = "ratio")
# # # 
### LOAD THE SAMPLE ANNOTATIONS ####
# # metadata_PRJNA512012_2 <- read.delim("20220531_Oscars_tool/ALS_Human_EWCE/Data/Metadata_PRJNA512012_ALS_MC.csv")
# # metadata_PRJNA512012 <- read.delim("20220531_Oscars_tool/ALS_Human_EWCE/Data/Metadata_PRJNA512012_ALS_MC.csv")[,c(1,4,13)] %>%
# #   `colnames<-`(c("names","Gender","APOE"))
# # cell_ratio_tidy_joined <- left_join(cell_ratio_tidy,
# #                                     metadata_PRJNA512012,
# #                                     by = "names")
# # cell_ratio_tidy_joined$APOE <- as.factor(cell_ratio_tidy_joined$APOE)
# # cell_ratio_tidy_joined_33_34 <- cell_ratio_tidy_joined %>%
# #   filter(APOE %in% c(33,34))
# 
# ggplot(cell_ratio_tidy_joined_33_34,
#        aes(x = celltype,
#            y = ratio,
#            fill = APOE)) +
#   geom_boxplot() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#   labs(x = "Celltypes",
#        y = "Estimated proportion") +
#   scale_fill_manual(name = "APOE",
#                     values = c("orange","brown"))
#   # +
#   # coord_flip() +
#   # scale_y_sqrt()
# 
# cell_ratio_tidy_joined_33_34$celltype <- factor(cell_ratio_tidy_joined_33_34$celltype,
#                                                 levels = c("Astrocytes",
#                                                            "Interneurons",
#                                                            "Microglia",
#                                                            "Oligodendrocyte.Precursor",
#                                                            "Oligodendrocytes",
#                                                            "Pericytes",
#                                                            "Pyramidal.Neurons",
#                                                            "Vascular.and.Leptomeningeal.Cells",
#                                                            "Vascular.Endothelial",
#                                                            "Vascular.Smooth.Muscle.Cell",
#                                                            "X.none."))
# bxp <- ggboxplot(cell_ratio_tidy_joined_33_34,
#           x = "celltype",
#            y = "ratio",
#            fill = "APOE") +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#   labs(x = "Celltypes",
#        y = "Estimated proportion") +
#   scale_fill_manual(name = "APOE",
#                     values = c("orange","brown")) +
#   facet_grid(rows = vars(Gender))
# 
# library(rstatix)
# cell_ratio_tidy_joined_33_34 %>%
#   group_by(Gender, celltype) %>%
#   shapiro_test(ratio)
# 
# stat.test <- cell_ratio_tidy_joined_33_34 %>%
#   group_by(Gender, celltype) %>%
#   wilcox_test(ratio ~ APOE) %>%
#   # adjust_pvalue(method = "BH") %>%
#   add_significance()
# stat.test
# 
# # write_csv(stat.test,
#           # file = "20220531_Oscars_tool/ALS_Human_EWCE/Output_READOUT_Human/Stats_PRJNA512012_Celltype_AvgSS_GenderSplit.txt")
# # write_csv(stat.test,
#           # file = "20220531_Oscars_tool/ALS_Human_EWCE/Output_READOUT_Human/Stats_PRJNA512012_Celltype_7238SS_GenderSplit.txt")
# 
# # stat.test_female <- cell_ratio_tidy_joined_33_34[cell_ratio_tidy_joined_33_34$Gender == "Female",] %>%
# #   group_by(celltype) %>%
# #   t_test(ratio ~ APOE) %>%
# #   adjust_pvalue(method = "BH") %>%
# #   add_significance()
# # 
# # stat.test_male <- cell_ratio_tidy_joined_33_34[cell_ratio_tidy_joined_33_34$Gender == "Male",] %>%
# #   group_by(celltype) %>%
# #   t_test(ratio ~ APOE) %>%
# #   adjust_pvalue(method = "BH") %>%
# #   add_significance()
# 
# stat.test <- stat.test %>%
#   add_xy_position(x = "celltype", dodge = 0.8)
# 
# bxp + 
#   stat_pvalue_manual(
#     stat.test,
#     label = "p.signif",
#     tip.length = 0.01,
#     hide.ns = TRUE) +
#   scale_y_continuous(expand = expansion(mult = c(0.01, 0.1)))
