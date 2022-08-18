###--------------------------------------------------------------------------###
### Project: ApoE in ALS                                                       -
### Purpose: Genotype assignment                                               -
### Authors: QC Lin                                                            -
### Last Edit: 2022-08-15, QC Lin                                              -
### Contact: q.c.lin@students.uu.nl or sebastian.lewandowski@ki.se             -
###--------------------------------------------------------------------------###

### Dependencies and general settings ------------------------------------------
# Load libraries
library(tidyverse)
library(ggbeeswarm)
library(viridis)
library(VariantAnnotation)
library(ggpubr)

# Set working directory and seed for reproducibility
setwd(".")
set.seed(747)

### ---- GENERAL ----
### Load in general resources --------------------------------------------------
protein_coding_genes <- scan("20220314_QCL_Genes_APOE/Protein_coding_genes.txt",
                             what = "")

### ---- PERCENTAGE OF READS ON APOE ----
### Load data - ROSMAP3_APOE ---------------------------------------------------
count_files <- list.files(path = "20220301_QCL_RNA_Seq/ROSMAP/processed",
                          pattern = ".count", 
                          full.names = T, 
                          recursive = T)

count_names <- list.files(path = "20220301_QCL_RNA_Seq/ROSMAP/processed",
                          pattern = ".count", 
                          full.names = F, 
                          recursive = T)
count_names <- gsub("\\/.*","",count_names)
counts <- purrr::map(count_files,
                     ~ read.delim(.,header = FALSE))
names(counts) <- count_names

colnames <- c("gene", "count")
counts <- lapply(counts, `names<-`, colnames)
counts_ROSMAP3 <- bind_rows(counts, .id = 'id')


extra <- c(tail(counts$RISK_105_redo$gene, n = 5))
counts_ROSMAP3 <- filter(counts_ROSMAP3,
                         !(gene %in% extra))
# counts_ROSMAP3 <- filter(counts_ROSMAP3,
#                          gene %in% protein_coding_genes)

total_gene_reads <- counts_ROSMAP3 %>% group_by(id) %>%
  summarise(sum = sum(count))
APOE_ROSMAP3 <- filter(counts_ROSMAP3,
                       gene %in% "ENSG00000130203")

APOE_ROSMAP3 <- left_join(APOE_ROSMAP3,
                          total_gene_reads,
                          by = "id")
APOE_ROSMAP3$percentage_APOE <- APOE_ROSMAP3$count/APOE_ROSMAP3$sum*100

### Load data - ROSMAP1_APOE ---------------------------------------------------
count_files <- list.files(path = "20220301_QCL_RNA_Seq/ROSMAP_batch1/processed",
                          pattern = ".count", 
                          full.names = T, 
                          recursive = T)
count_names <- list.files(path = "20220301_QCL_RNA_Seq/ROSMAP_batch1/processed",
                          pattern = ".count", 
                          full.names = F, 
                          recursive = T)
count_names <- gsub("\\/.*","",count_names)
counts <- purrr::map(count_files,
                     ~ read.delim(.,header = FALSE))
names(counts) <- count_names

colnames <- c("gene", "count")
counts <- lapply(counts, `names<-`, colnames)
counts_ROSMAP1 <- bind_rows(counts, .id = 'id')

counts_ROSMAP1 <- filter(counts_ROSMAP1,
                         !(gene %in% extra))
# counts_ROSMAP1 <- filter(counts_ROSMAP1,
#                          gene %in% protein_coding_genes)

total_gene_reads <- counts_ROSMAP1 %>% group_by(id) %>%
  summarise(sum = sum(count))
APOE_ROSMAP1 <- filter(counts_ROSMAP1,
                       gene %in% "ENSG00000130203")
APOE_ROSMAP1 <- left_join(APOE_ROSMAP1,
                          total_gene_reads,
                          by = "id")
APOE_ROSMAP1$percentage_APOE <- APOE_ROSMAP1$count/APOE_ROSMAP1$sum*100

### Load data - HUMAN_ALS_APOE -------------------------------------------------
count_files <- list.files(path = "20220301_QCL_RNA_Seq/PRJNA512012/processed",
                          pattern = ".count", 
                          full.names = T, 
                          recursive = T)
count_names <- list.files(path = "20220301_QCL_RNA_Seq/PRJNA512012/processed",
                          pattern = ".count", 
                          full.names = F, 
                          recursive = T)
count_names <- gsub("\\/.*","",count_names)
counts <- purrr::map(count_files,
                     ~ read.delim(.,header = FALSE))
names(counts) <- count_names

colnames <- c("gene", "count")
counts <- lapply(counts, `names<-`, colnames)
counts_PRJNA512012 <- bind_rows(counts, .id = 'id')

counts_PRJNA512012 <- filter(counts_PRJNA512012,
                             !(gene %in% extra))
# counts_PRJNA512012 <- filter(counts_PRJNA512012,
#                              gene %in% protein_coding_genes)

total_gene_reads <- counts_PRJNA512012 %>% group_by(id) %>%
  summarise(sum = sum(count))
APOE_PRJNA512012 <- filter(counts_PRJNA512012,
                           gene %in% "ENSG00000130203")
APOE_PRJNA512012 <- left_join(APOE_PRJNA512012,
                              total_gene_reads,
                              by = "id")
APOE_PRJNA512012$percentage_APOE <- APOE_PRJNA512012$count/APOE_PRJNA512012$sum*100

### Combined percentages -------------------------------------------------------
APOE <- list(APOE_ROSMAP3 = APOE_ROSMAP3,
             APOE_PRJNA512012 = APOE_PRJNA512012,
             APOE_ROSMAP1 = APOE_ROSMAP1)
APOE <- bind_rows(APOE,
                  .id = "dataset")

ggplot(APOE,
       aes(x = dataset, y = percentage_APOE)) +
  geom_violin(alpha = 0.5) +
  geom_beeswarm(aes(colour = dataset)) +
  #scale_color_viridis(option = "turbo",
  #                    discrete = TRUE) +
  labs(x = NULL,
       y = "Reads counted for APOE by HTSeq-Counts (%)") +
  theme(legend.position="none") +
  scale_x_discrete(labels=c("APOE_PRJNA512012" = "PRJNA512012", "APOE_ROSMAP1" = "ROSMAP batch 1",
                            "APOE_ROSMAP3" = "ROSMAP batch 3"))

ggplot(APOE,
       aes(x = dataset, y = count)) +
  geom_boxplot(aes(fill = dataset)) +
  #scale_fill_viridis(option = "turbo",
  #                    discrete = TRUE) +
  scale_y_log10() +
  labs(x = NULL,
       y = "Absolute read count") +
  theme(legend.position="none") +
  scale_x_discrete(labels=c("APOE_PRJNA512012" = "PRJNA512012", "APOE_ROSMAP1" = "ROSMAP batch 1",
                          "APOE_ROSMAP3" = "ROSMAP batch 3"))

### ---- UNIQUELY MAPPED READS STAR ----
### Load data - ROSMAP3_APOE ---------------------------------------------------
ROSMAP3_STAR <- read.delim("20220301_QCL_RNA_Seq/ROSMAP/stats/Uniquely_mapped_reads_STAR.txt",
                           header = FALSE)
colnames(ROSMAP3_STAR) <- c("sample","reads")
ROSMAP3_STAR$est_full_length <- ROSMAP3_STAR$reads*300

### Load data - ROSMAP1_APOE ---------------------------------------------------
ROSMAP1_STAR <- read.delim("20220301_QCL_RNA_Seq/ROSMAP_batch1/stats/Uniquely_mapped_reads_STAR.txt",
                           header = FALSE)
colnames(ROSMAP1_STAR) <- c("sample","reads")
ROSMAP1_STAR$est_full_length <- ROSMAP1_STAR$reads*200

### Load data - HUMAN_ALS_APOE -------------------------------------------------
PRJNA512012_STAR <- read.delim("20220301_QCL_RNA_Seq/PRJNA512012/stats/Uniquely_mapped_reads_STAR.txt",
                               header = FALSE)
colnames(PRJNA512012_STAR) <- c("sample","reads")
PRJNA512012_STAR$est_full_length <- PRJNA512012_STAR$reads*250

### Combine datasets - STAR ----------------------------------------------------
combined_datasets <- list(ROSMAP3_STAR = ROSMAP3_STAR,
                          PRJNA512012_STAR = PRJNA512012_STAR,
                          ROSMAP1_STAR = ROSMAP1_STAR) %>%
  bind_rows(.id = "dataset")

ggplot(combined_datasets,
       aes(x = dataset, y = reads)) +
  geom_boxplot(aes(fill = dataset)) +
  #scale_fill_viridis(option = "viridis",
  #                   discrete = TRUE) +
  labs(x = NULL,
       y = "Uniquely mapped reads by STAR") +
  theme(legend.position="none") +
  scale_x_discrete(labels=c("PRJNA512012_STAR" = "PRJNA512012", "ROSMAP1_STAR" = "ROSMAP batch 1",
                            "ROSMAP3_STAR" = "ROSMAP batch 3"))

ggplot(combined_datasets,
       aes(x = dataset, y = est_full_length)) +
  geom_boxplot(aes(fill = dataset)) +
  #scale_fill_viridis(option = "viridis",
  #                   discrete = TRUE) +
  labs(x = NULL,
       y = "Estimated Nucleotides Covered") +
  theme(legend.position="none") +
  scale_y_log10() +
  scale_x_discrete(labels=c("PRJNA512012_STAR" = "PRJNA512012", "ROSMAP1_STAR" = "ROSMAP batch 1",
                            "ROSMAP3_STAR" = "ROSMAP batch 3"))

### ---- VCFs - HUMAN_ALS (FS & QD) ----
### Load data - HUMAN_ALS ------------------------------------------------------
vcf_files <- list.files(path = "20220301_QCL_RNA_Seq/PRJNA512012/processed",
                        pattern = "_filtered.vcf.gz.tbi", 
                        full.names = T,
                        recursive = T)
vcf_files <- gsub(".tbi", "",vcf_files)
vcf_names <- list.files(path = "20220301_QCL_RNA_Seq/PRJNA512012/processed",
                        pattern = "_filtered.vcf.gz.tbi", 
                        full.names = F,
                        recursive = T)
vcf_names <- gsub("\\/.*","",vcf_names)
variant_calls_PRJNA512012 <- purrr::map(vcf_files,
                                        readVcf)
names(variant_calls_PRJNA512012) <- vcf_names

### FS & QD for region 19:44500000-45000000 ------------------------------------
# FS
variant_calls_FS_PRJNA512012 <- purrr::map(variant_calls_PRJNA512012,
                                           ~ as_tibble(info(.x)$FS)) %>%
  bind_rows(.id = "sample")
# QD
variant_calls_QD_PRJNA512012 <- purrr::map(variant_calls_PRJNA512012,
                                           ~ as_tibble(info(.x)$QD)) %>%
  bind_rows(.id = "sample")

### FS & QD for APOE - rs7412|rs429358 -----------------------------------------
variant_calls_info_PRJNA512012 <- purrr::map(variant_calls_PRJNA512012,
                                             ~ info(.x) %>%
                                               subset(rownames(.x) %in% c("rs7412","rs429358"))) %>%
  keep(~ nrow(.) > 0)
# FS
variant_calls_apoe_FS_PRJNA512012 <- purrr::map(variant_calls_info_PRJNA512012,
                                                ~ as_tibble(.$FS)) %>%
  bind_rows(.id = "sample")
# QD
variant_calls_apoe_QD_PRJNA512012 <- purrr::map(variant_calls_info_PRJNA512012,
                                                ~ as_tibble(.$QD)) %>%
  bind_rows(.id = "sample")

# DP
variant_calls_apoe_DP_PRJNA512012 <- purrr::map(variant_calls_info_PRJNA512012,
                                                ~ as_tibble(.$DP)) %>%
  bind_rows(.id = "sample")


### ---- VCFs - ROSMAP3 (FS & QD) ----
### Load data - ROSMAP3 ------------------------------------------------------
vcf_files <- list.files(path = "20220301_QCL_RNA_Seq/ROSMAP/processed",
                        pattern = "_filtered.vcf.gz.tbi", 
                        full.names = T,
                        recursive = T)
vcf_files <- gsub(".tbi", "",vcf_files)
vcf_names <- list.files(path = "20220301_QCL_RNA_Seq/ROSMAP/processed",
                        pattern = "_filtered.vcf.gz.tbi", 
                        full.names = F,
                        recursive = T)
vcf_names <- gsub("\\/.*","",vcf_names)
variant_calls_ROSMAP3 <- purrr::map(vcf_files,
                                    readVcf)
names(variant_calls_ROSMAP3) <- vcf_names

### FS & QD for region 19:44500000-45000000 ------------------------------------
# FS
variant_calls_FS_ROSMAP3 <- purrr::map(variant_calls_ROSMAP3,
                                       ~ as_tibble(info(.x)$FS)) %>%
  bind_rows(.id = "sample")
# QD
variant_calls_QD_ROSMAP3 <- purrr::map(variant_calls_ROSMAP3,
                                       ~ as_tibble(info(.x)$QD)) %>%
  bind_rows(.id = "sample")



### FS & QD for APOE - rs7412|rs429358 -----------------------------------------
variant_calls_info_ROSMAP3 <- purrr::map(variant_calls_ROSMAP3,
                                         ~ info(.x) %>%
                                           subset(rownames(.x) %in% c("rs7412","rs429358"))) %>%
  keep(~ nrow(.) > 0)
# FS
variant_calls_apoe_FS_ROSMAP3 <- purrr::map(variant_calls_info_ROSMAP3,
                                            ~ as_tibble(.$FS)) %>%
  bind_rows(.id = "sample")
# QD
variant_calls_apoe_QD_ROSMAP3 <- purrr::map(variant_calls_info_ROSMAP3,
                                            ~ as_tibble(.$QD)) %>%
  bind_rows(.id = "sample")
# DP
variant_calls_apoe_DP_ROSMAP3 <- purrr::map(variant_calls_info_ROSMAP3,
                                            ~ as_tibble(.$DP)) %>%
  bind_rows(.id = "sample")

### ---- VCFs - ROSMAP1 (FS & QD) ----
### Load data - ROSMAP1 ------------------------------------------------------
vcf_files <- list.files(path = "20220301_QCL_RNA_Seq/ROSMAP_batch1/processed",
                        pattern = "_filtered.vcf.gz.tbi", 
                        full.names = T,
                        recursive = T)
vcf_files <- gsub(".tbi", "",vcf_files)
vcf_names <- list.files(path = "20220301_QCL_RNA_Seq/ROSMAP_batch1/processed",
                        pattern = "_filtered.vcf.gz.tbi", 
                        full.names = F,
                        recursive = T)
vcf_names <- gsub("\\/.*","",vcf_names)
variant_calls_ROSMAP1 <- purrr::map(vcf_files,
                                    readVcf)
names(variant_calls_ROSMAP1) <- vcf_names

### FS & QD for region 19:44500000-45000000 ------------------------------------
# FS
variant_calls_FS_ROSMAP1 <- purrr::map(variant_calls_ROSMAP1,
                                       ~ as_tibble(info(.x)$FS)) %>%
  bind_rows(.id = "sample")
# QD
variant_calls_QD_ROSMAP1 <- purrr::map(variant_calls_ROSMAP1,
                                       ~ as_tibble(info(.x)$QD)) %>%
  bind_rows(.id = "sample")


### FS & QD for APOE - rs7412|rs429358 -----------------------------------------
variant_calls_info_ROSMAP1 <- purrr::map(variant_calls_ROSMAP1,
                                         ~ info(.x) %>%
                                           subset(rownames(.x) %in% c("rs7412","rs429358"))) %>%
  keep(~ nrow(.) > 0)
# FS
variant_calls_apoe_FS_ROSMAP1 <- purrr::map(variant_calls_info_ROSMAP1,
                                            ~ as_tibble(.$FS)) %>%
  bind_rows(.id = "sample")
# QD
variant_calls_apoe_QD_ROSMAP1 <- purrr::map(variant_calls_info_ROSMAP1,
                                            ~ as_tibble(.$QD)) %>%
  bind_rows(.id = "sample")
# DP
variant_calls_apoe_DP_ROSMAP1 <- purrr::map(variant_calls_info_ROSMAP1,
                                            ~ as_tibble(.$DP)) %>%
  bind_rows(.id = "sample")
# BQRS
variant_calls_apoe_BQSR_ROSMAP1 <- purrr::map(variant_calls_info_ROSMAP1,
                                            ~ as_tibble(.$BaseQRankSum)) %>%
  bind_rows(.id = "sample")

### ---- PLOTS ----
density_standard <- function(dataset){
  ggplot(dataset,
         aes(x = value)) +
    geom_density()
}
density_standard(variant_calls_FS_PRJNA512012) +
  geom_density(data = variant_calls_FS_ROSMAP3) +
  geom_density(data = variant_calls_FS_ROSMAP1)

ggplot(variant_calls_FS_PRJNA512012,
       aes(x = value)) +
  geom_area(stat = "bin", binwidth = 1, alpha = 0.25, fill = "green") +
  geom_area(data = variant_calls_FS_ROSMAP3,
            stat = "bin", binwidth = 1, alpha = 0.5, fill = "purple") +
  geom_area(data = variant_calls_FS_ROSMAP1,
            stat = "bin", binwidth = 1, alpha = 0.25, fill = "red") +
  labs(x = "FisherStrand",
       y = "Counts") +
  scale_y_continuous(trans = "pseudo_log")

ggplot(variant_calls_QD_PRJNA512012,
       aes(x = value)) +
  geom_area(stat = "bin", binwidth = 1, alpha = 0.25, fill = "green") +
  geom_area(data = variant_calls_QD_ROSMAP3,
            stat = "bin", binwidth = 1, alpha = 0.5, fill = "purple") +
  geom_area(data = variant_calls_QD_ROSMAP1,
            stat = "bin", binwidth = 1, alpha = 0.25, fill = "red") +
  labs(x = "QualityByDepth",
       y = "Counts")

ggplot(variant_calls_apoe_DP_PRJNA512012,
       aes(x = value)) +
  geom_histogram(stat = "bin", binwidth = 100,alpha = 0.25, fill = "green") +
  geom_histogram(data = variant_calls_apoe_DP_ROSMAP3,
                stat = "bin", binwidth = 100, alpha = 0.5, fill = "purple") +
  geom_histogram(data = variant_calls_apoe_DP_ROSMAP1,
                 stat = "bin", binwidth = 100, alpha = 0.25, fill = "red") +
  labs(x = "Read Depth",
       y = "Counts")

DP <- list(PRJNA512012 = variant_calls_apoe_DP_PRJNA512012,
           ROSMAP3 = variant_calls_apoe_DP_ROSMAP3,
           ROSMAP1 = variant_calls_apoe_DP_ROSMAP1) %>%
  bind_rows(.id = "set")

ggplot(DP,
       aes(x = set,
           y = value)) +
  geom_violin(alpha = 0.5) +
  geom_beeswarm(aes(colour = set)) +
  scale_y_log10() +
  #scale_color_viridis(option = "turbo",
  #                    discrete = TRUE) +
  labs(x = "Dataset",
       y = "Read Depth") +
  theme(legend.position = "none") +
  scale_x_discrete(labels=c("PRJNA512012" = "PRJNA512012", "ROSMAP1" = "ROSMAP batch 1",
                            "ROSMAP3" = "ROSMAP batch 3"))

### ---- SNP CALLING ----
### Genotyping ROSMAP1 ---------------------------------------------------------
variant_calls_ROSMAP1_filtered <- map(variant_calls_ROSMAP1,
                                      ~ subset(.x, FILTER == "PASS"))

variant_calls_info_ROSMAP1_filtered <- purrr::map(variant_calls_ROSMAP1_filtered,
                                         ~ info(.x) %>%
                                           subset(rownames(.x) %in% c("rs7412","rs429358"))) %>%
  keep(~ nrow(.) > 0)

variant_calls_info_ROSMAP1_filtered <- purrr::map(variant_calls_info_ROSMAP1_filtered,
                                                  ~ as.data.frame(.x) %>%
                                                    rownames_to_column)

rsID_ROSMAP1_filtered <- purrr::map(variant_calls_info_ROSMAP1_filtered,
                                    ~ dplyr::select(.x,rowname,AF)) %>%
  bind_rows(.id = "sample")

rsID_ROSMAP1_filtered$AF <- unlist(rsID_ROSMAP1_filtered$AF)

APOE_genotype <- rsID_ROSMAP1_filtered %>%
  group_by(sample) %>%
  pivot_wider(names_from = "rowname", values_from = "AF") %>%
  mutate(APOE = case_when(rs429358 == 1.0 ~ "44",
                          rs7412 == 1.0 ~ "22",
                          rs429358 == 0.5 && rs7412 == 0.5 ~ "24",
                          rs429358 == 0.5 ~ "34",
                          rs7412 == 0.5 ~ "23"
  ))

table(APOE_genotype$APOE)
samples_genotype <- read.csv("20220314_QCL_Genes_APOE/Final/120_samples_genotype.txt", sep="") %>%
  `colnames<-`(c("sample","apoe_genotype"))
test <- left_join(samples_genotype,
                  APOE_genotype,
                  by = "sample")
test$APOE[is.na(test$APOE)] <- 33
test <- mutate(test,
               correct = APOE == apoe_genotype)

test <- test %>%
  mutate(type = case_when(correct == TRUE & APOE == "33" ~ "TN",
                          correct == TRUE & APOE != "33" ~ "TP",
                          correct == FALSE & APOE == "33" ~ "FN",
                          correct == FALSE & APOE != "33" ~ "FP"))
# write.table(test,
#             file = "20220314_QCL_Genes_APOE/Final/VCF_output.txt",
#             row.names = FALSE,
#             quote = FALSE)

table(test$type)

#   mutate(prediction = if_else(apoe_genotype == "33",
#                               0,
#                               1))
# prediction(test$APOE, test$apoe_genotype)

test_DP_gen <- left_join(test_DP, test, by = "sample")

ggplot(test_DP_gen,
       aes(x = set,
           y = value)) +
  geom_violin(alpha = 0.5) +
  geom_beeswarm(aes(col = set)) +
  geom_beeswarm(aes(col = correct)) +
  scale_y_log10() +
  scale_color_viridis(option = "turbo",
                      discrete = TRUE) +
  labs(x = "Dataset",
       y = "Read Depth") +
  theme(legend.position = "none")

### Genotyping HUMAN ALS -------------------------------------------------------
variant_calls_HUMAN_ALS_filtered <- map(variant_calls_PRJNA512012,
                                      ~ subset(.x, FILTER == "PASS"))

variant_calls_info_HUMAN_ALS_filtered <- purrr::map(variant_calls_HUMAN_ALS_filtered,
                                                  ~ info(.x) %>%
                                                    subset(rownames(.x) %in% c("rs7412","rs429358"))) %>%
  keep(~ nrow(.) > 0)

variant_calls_info_HUMAN_ALS_filtered <- purrr::map(variant_calls_info_HUMAN_ALS_filtered,
                                                  ~ as.data.frame(.x) %>%
                                                    rownames_to_column)

rsID_HUMAN_ALS_filtered <- purrr::map(variant_calls_info_HUMAN_ALS_filtered,
                                    ~ dplyr::select(.x,rowname,AF)) %>%
  bind_rows(.id = "sample")

rsID_HUMAN_ALS_filtered$AF <- unlist(rsID_HUMAN_ALS_filtered$AF)

APOE_genotype_HUMAN_ALS <- rsID_HUMAN_ALS_filtered %>%
  group_by(sample) %>%
  pivot_wider(names_from = "rowname", values_from = "AF") %>%
  mutate(APOE = case_when(rs429358 == 1.0 ~ "44",
                          rs7412 == 1.0 ~ "22",
                          rs429358 == 0.5 && rs7412 == 0.5 ~ "24",
                          rs429358 == 0.5 ~ "34",
                          rs7412 == 0.5 ~ "23"
  ))
table(APOE_genotype_HUMAN_ALS$APOE)

names_HUMAN_ALS <- as.data.frame(count_names) %>%
  `colnames<-`("sample")

APOE_genotype_HUMAN_ALS_table <- left_join(names_HUMAN_ALS,
                                           APOE_genotype_HUMAN_ALS,
                                           by = "sample")
APOE_genotype_HUMAN_ALS_table$APOE[is.na(APOE_genotype_HUMAN_ALS_table$APOE)] <- 33
table(APOE_genotype_HUMAN_ALS_table$APOE)

# write.table(APOE_genotype_HUMAN_ALS_table,
#             file = "20220314_QCL_Genes_APOE/Final/VCF_output_PRJNA512012.txt",
#             row.names = FALSE,
#             quote = FALSE)

