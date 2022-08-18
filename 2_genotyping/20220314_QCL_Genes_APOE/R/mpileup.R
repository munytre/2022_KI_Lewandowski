###--------------------------------------------------------------------------###
### Project: ApoE in ALS                                                       -
### Purpose: Read depth analysis (mpileup)                                     -
### Authors: QC Lin                                                            -
### Last Edit: 2022-08-15, QC Lin                                              -
### Contact: q.c.lin@students.uu.nl or sebastian.lewandowski@ki.se             -
###--------------------------------------------------------------------------###

### Dependencies and general settings ------------------------------------------
library(tidyverse)
library(ggpubr)

### Custom analysis function ---------------------------------------------------
# rs7412
rs7412 <- function(location){
  count_files_PRJ <- list.files(path = location,
                                pattern = "mpileup.txt", 
                                full.names = T, 
                                recursive = T)
  count_names_PRJ <- list.files(path = location,
                                pattern = "mpileup.txt",
                                full.names = F, 
                                recursive = T)
  count_names_PRJ <- gsub("\\/.*","",count_names_PRJ)
  counts_PRJ <- purrr::map(count_files_PRJ,
                           ~ read.delim(.,header = FALSE, quote = ""))
  names(counts_PRJ) <- count_names_PRJ
  PRJ <- lapply(counts_PRJ, function(x) x[(names(x) %in% c("V2", "V4"))])
  PRJ_SNP <- PRJ %>%
    purrr::reduce(dplyr::full_join, by = "V2")
  PRJ_SNP_rs7412 <- filter(PRJ_SNP,
                           V2 == "44908822")
  rownames(PRJ_SNP_rs7412) <- PRJ_SNP_rs7412[,1]
  PRJ_SNP_rs7412 <- PRJ_SNP_rs7412[,-1]
  colnames(PRJ_SNP_rs7412) <- count_names_PRJ
  PRJ_SNP_rs7412[is.na(PRJ_SNP_rs7412)] <- 0
  return(PRJ_SNP_rs7412)
}

# rs429358
rs429358 <- function(location){
  count_files_PRJ <- list.files(path = location,
                                pattern = "mpileup.txt", 
                                full.names = T, 
                                recursive = T)
  count_names_PRJ <- list.files(path = location,
                                pattern = "mpileup.txt",
                                full.names = F, 
                                recursive = T)
  count_names_PRJ <- gsub("\\/.*","",count_names_PRJ)
  counts_PRJ <- purrr::map(count_files_PRJ,
                           ~ read.delim(.,header = FALSE, quote = ""))
  names(counts_PRJ) <- count_names_PRJ
  PRJ <- lapply(counts_PRJ, function(x) x[(names(x) %in% c("V2", "V4"))])
  PRJ_SNP <- PRJ %>%
    purrr::reduce(dplyr::full_join, by = "V2")
  PRJ_SNP_rs429358 <- filter(PRJ_SNP,
                           V2 == "44908684")
  rownames(PRJ_SNP_rs429358) <- PRJ_SNP_rs429358[,1]
  PRJ_SNP_rs429358 <- PRJ_SNP_rs429358[,-1]
  colnames(PRJ_SNP_rs429358) <- count_names_PRJ
  PRJ_SNP_rs429358[is.na(PRJ_SNP_rs429358)] <- 0
  return(PRJ_SNP_rs429358)
}

### Load data and analysis  ----------------------------------------------------
# Load rs7412
PRJNA512012 <- rs7412("20220301_QCL_RNA_Seq/PRJNA512012_mpileup/processed/")
ROSMAP1 <- rs7412("20220301_QCL_RNA_Seq/ROSMAP_batch1_mpileup/processed/")
ROSMAP3 <- rs7412("20220301_QCL_RNA_Seq/ROSMAP_batch3_mpileup/processed/")

PRJNA512012_t <- as.data.frame(t(PRJNA512012))
ROSMAP1_t <- as.data.frame(t(ROSMAP1))
ROSMAP3_t <- as.data.frame(t(ROSMAP3))
rs7412_list <- list(PRJNA512012 = PRJNA512012_t, ROSMAP1 = ROSMAP1_t, ROSMAP3 = ROSMAP3_t)

rs7412_coverage <- bind_rows(rs7412_list,
          .id = "source") %>%
  `colnames<-`(c("source", "coverage"))

compare_means(coverage ~ source,  data = rs7412_coverage)
my_comparisons <- list( c("PRJNA512012", "ROSMAP1"), c("PRJNA512012", "ROSMAP3"), c("ROSMAP1", "ROSMAP3") )

# Plot results
ggboxplot(rs7412_coverage,
          x = "source",
          y = "coverage",
          fill = "source") +
  yscale("log10") +
  stat_compare_means(comparisons = list(c("PRJNA512012", "ROSMAP1"),
                                        c("PRJNA512012", "ROSMAP3"),
                                        c("ROSMAP1", "ROSMAP3"))) +
  stat_compare_means(label.y = 5) +
  scale_x_discrete(labels=c("PRJNA512012" = "PRJNA512012",
                            "ROSMAP1" = "ROSMAP batch 1",
                            "ROSMAP3" = "ROSMAP batch 3")) +
  labs(x = "Dataset",
       y = "Coverage") +
  theme_gray() +
  theme(legend.position = "none")

# Load rs429358
PRJNA512012_2 <- rs429358("20220301_QCL_RNA_Seq/PRJNA512012_mpileup/processed/")
ROSMAP1_2 <- rs429358("20220301_QCL_RNA_Seq/ROSMAP_batch1_mpileup/processed/")
ROSMAP3_2 <- rs429358("20220301_QCL_RNA_Seq/ROSMAP_batch3_mpileup/processed/")

PRJNA512012_t_2 <- as.data.frame(t(PRJNA512012_2))
ROSMAP1_t_2 <- as.data.frame(t(ROSMAP1_2))
ROSMAP3_t_2 <- as.data.frame(t(ROSMAP3_2))
rs429358_list <- list(PRJNA512012 = PRJNA512012_t_2, ROSMAP1 = ROSMAP1_t_2, ROSMAP3 = ROSMAP3_t_2)

rs429358_coverage <- bind_rows(rs429358_list,
                             .id = "source") %>%
  `colnames<-`(c("source", "coverage"))

# Plot results
ggboxplot(rs429358_coverage,
          x = "source",
          y = "coverage",
          fill = "source") +
  yscale("log10") +
  stat_compare_means(comparisons = list(c("PRJNA512012", "ROSMAP1"),
                                        c("PRJNA512012", "ROSMAP3"),
                                        c("ROSMAP1", "ROSMAP3"))) +
  stat_compare_means(label.y = 5) +
  scale_x_discrete(labels=c("PRJNA512012" = "PRJNA512012",
                            "ROSMAP1" = "ROSMAP batch 1",
                            "ROSMAP3" = "ROSMAP batch 3")) +
  labs(x = "Dataset",
       y = "Coverage") +
  theme_gray() +
  theme(legend.position = "none")
