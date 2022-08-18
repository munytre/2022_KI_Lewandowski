###--------------------------------------------------------------------------###
### Project: ApoE in ALS                                                       -
### Purpose: DESeq2 to explore gender effects in GSE124439                     -
### Authors: QC Lin                                                            -
### Last Edit: 2022-08-12, QC Lin                                              -
### Contact: q.c.lin@students.uu.nl or sebastian.lewandowski@ki.se             -
###--------------------------------------------------------------------------###

### Dependencies and general settings ------------------------------------------

# Load libraries
library(GEOquery)
library(tidyverse)
library(DESeq2)
library(ggpubr)
library(rstatix)

# Set working directory and seed for reproducibility
setwd(".")
set.seed(747)

### Custom PCA function --------------------------------------------------------
plotPCA.custom <- function (object, intgroup = "condition", ntop = 500, returnData = FALSE){
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, 
                                               drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }
  else {
    colData(object)[[intgroup]]
  }
  d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3], PC4 = pca$x[, 4], group = group, 
                  intgroup.df, name = colnames(object))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:2]
    return(d)
  }
  ggplot(data = d, aes_string(x = "PC2", y = "PC3", color = "group")) + 
    geom_point(size = 3) +
    xlab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
    ylab(paste0("PC3: ", round(percentVar[3] * 100), "% variance")) +
    coord_fixed()
}

### Load data ------------------------------------------------------------------
# Load in metadata (1)
gse <- getGEO(filename = "metadata/GSE124439_series_matrix.txt")
metadata <- gse@phenoData@data
colnames(metadata) <- sub(".", "_", colnames(metadata), fixed=TRUE)
colnames(metadata) <- sub(":", "_", colnames(metadata), fixed=TRUE)
colnames(metadata) <- sub(" ", "_", colnames(metadata), fixed=TRUE)

# Split disease vs control (optional)
# metadata <- metadata %>%
# mutate(condition = characteristics_ch1_3 %in% c("sample group: ALS Spectrum MND"))
# metadata$condition[metadata$condition == TRUE] <- "ALS"
# metadata$condition[metadata$condition == FALSE] <- "Control"

# Load in metadata (2)
metadata_excel1 <- read.delim("metadata/GSE124439_sheet1.txt",
                              skip = 1,
                              nrows = 176)
colnames(metadata_excel1)[1] <- "title"
metadata_excel2 <- read.delim("metadata/GSE124439_sheet2.txt",
                              skip = 1,
                              nrows = 77)
colnames(metadata_excel2)[1] <- "subject_id_ch1"
metadata_excel2 <- metadata_excel2[,c(1,5:8)]
metadata <- left_join(metadata, metadata_excel1, by = "title")
metadata <- left_join(metadata, metadata_excel2, by = "subject_id_ch1")
metadata <- metadata[,-c(10:43)]

# Load in metadata (3)
SRR_sample <- read.delim("metadata/PRJNA512012_tsv.txt")
SRR_sample <- SRR_sample[,c(2,4)]
colnames(SRR_sample) <- c("SRR", "title")

# Load in sex chromosomal gene names
# Y
y_genes <- scan(file = "sex_linked_genes/Y_specific.txt",
                what = "",
                quiet = TRUE)
# X & Y
x_y_genes <- scan(file = "sex_linked_genes/XY.txt",
                  what = "",
                  quiet = TRUE)

# Load in count files
count_files <- list.files(path = "data",
                          pattern = ".txt", 
                          full.names = T, 
                          recursive = T)
counts <- purrr::map(count_files,
                     read.delim)
names(counts) <- row.names(metadata)

count_matrix <- counts %>%
  purrr::reduce(dplyr::left_join, by = "gene.TE")

count_matrix_nonfiltered <- count_matrix

rownames(count_matrix_nonfiltered) <- count_matrix_nonfiltered[,1]
count_matrix_nonfiltered <- count_matrix_nonfiltered[,-1]
colnames(count_matrix_nonfiltered) <- row.names(metadata)
count_matrix_nonfiltered <- as.matrix(count_matrix_nonfiltered)

# # Sex-linked genes filtering
# count_matrix <- filter(count_matrix,
#                        !(gene.TE %in% y_genes))

rownames(count_matrix) <- count_matrix[,1]
count_matrix <- count_matrix[,-1]
colnames(count_matrix) <- row.names(metadata)
count_matrix <- as.matrix(count_matrix)

### Sex Chromosome count analysis ----------------------------------------------
genes_oscar_1 <- c("AKAP17A", "ASMT", "ASMTL", "CD99", "CRLF2", "CSF2RA",
                   "DHRSX", "GTPBP6", "IL3RA", "IL9R", "P2RY8", "PLCXD1",
                   "PPP2R3B", "SHOX", "SLC25A6", "VAMP7", "WASH6P", "ZBED1")
genes_oscar_2 <- c("AKAP17A", "ASMT", "ASMTL", "CD99", "CRLF2", "CSF2RA",
                   "DHRSX", "GTPBP6", "IL3RA", "IL9R", "P2RY8", "PLCXD1",
                   "PPP2R3B", "SHOX", "SLC25A6", "VAMP7", "WASH6P", "ZBED1")
genes_oscar_3 <- c("AMELX", "DDX3X", "EIF1AX", "KDM5C", "KDM6A", "NLGN4X",
                  "PCDH11X", "RPS4X", "TBL1X", "TGIF2LX", "TMSB4X", "USP9X",
                  "VCX3A", "VCX3B", "VCX", "VCX2", "ZFX")
genes_oscar_4 <- c("AMELY", "DDX3Y", "EIF1AY", "KDM5D", "UTY", "NLGN4Y",
                  "PCDH11Y", "RPS4Y1", "RPS4Y2", "TBL1Y", "TGIF2LY", "TMSB4Y",
                  "USP9Y", "VCY", "VCY1B", "ZFY")
genes_oscar_5 <- c("BPY2", "BPY2B", "BPY2C", "CDY1", "CDY1B", "CDY2A", "CDY2B",
                   "DAZ1", "DAZ2", "DAZ3", "DAZ4", "HSFY1", "HSFY2", "PRY",
                   "PRY2", "PRYP3", "RBMY1A1", "RBMY1B", "RBMY1D", "RBMY1E",
                   "RBMY1F", "RBMY1J", "SRY", "TSPY1", "TSPY10", "TSPY2",
                   "TSPY3", "TSPY4", "TSPY8", "TSPY9P")
genes_oscar_total <- unique(c(genes_oscar_1,
                              genes_oscar_2,
                              genes_oscar_3,
                              genes_oscar_4,
                              genes_oscar_5))

# Load in transcript lengths of genes of interest
genes_oscar_length <- read.delim("APOE_homology/genes_oscar_length.txt") %>%
  group_by(Gene.name) %>%
  summarise(avg = mean(Transcript.length..including.UTRs.and.CDS.))

counts_oscar <- subset(count_matrix_nonfiltered,
                       rownames(count_matrix_nonfiltered) %in% genes_oscar_total) %>%
  as.data.frame()

# Calculate TPMs
counts_oscar$Gene.name <- rownames(counts_oscar)
counts_oscar <- left_join(counts_oscar, genes_oscar_length, by = "Gene.name")
x <- counts_oscar %>%
  select(c(1:176)) / counts_oscar$avg
tpm_oscar <- t( t(x) * 1e6 / colSums(x) )
rownames(tpm_oscar) <- counts_oscar$Gene.name

tpm_oscar_df <- as.data.frame(t(tpm_oscar))
joined_tpm <- bind_cols(tpm_oscar_df, metadata)

joined_tpm_long <- pivot_longer(joined_tpm,
                                cols = 1:78,
                                names_to = "Gene",
                                values_to = "TPM")

# Assign genes to their respective gene group
joined_tpm_long$group <- NA
joined_tpm_long$group[joined_tpm_long$Gene %in% genes_oscar_1 & joined_tpm_long$Gender == "Female"] <- "1"
joined_tpm_long$group[joined_tpm_long$Gene %in% genes_oscar_2 & joined_tpm_long$Gender == "Male"] <- "1"
joined_tpm_long$group[joined_tpm_long$Gene %in% genes_oscar_3] <- "3"
joined_tpm_long$group[joined_tpm_long$Gene %in% genes_oscar_4] <- "4"
joined_tpm_long$group[joined_tpm_long$Gene %in% genes_oscar_5] <- "5"

# Statistical testing
pwc <- joined_tpm_long %>%
  group_by(group) %>%
  pairwise_wilcox_test(TPM ~ Gender,
              p.adjust.method = "BH")

# Add position to Wilcox test for ggplot
pwc <- pwc %>% add_xy_position(x = "group", group = "Gender", y.trans = log10)

# Plotting
ggboxplot(joined_tpm_long,
          x = "group",
          y = "TPM",
          color = "Gender") +
  yscale("log10", .format = TRUE) +
  stat_pvalue_manual(pwc, step.increase = 0.1) +
  labs(
    caption = get_pwc_label(pwc)
  )

# Statistical testing
pwt <- joined_tpm_long %>%
  group_by(group) %>%
  pairwise_t_test(TPM ~ Gender,
         p.adjust.method = "BH")

# Add position to T-test test for ggplot
pwt <- pwt %>% add_xy_position(x = "group", group = "Gender", y.trans = log10)

# Plotting
ggboxplot(joined_tpm_long,
          x = "group",
          y = "TPM",
          color = "Gender") +
  yscale("log10", .format = TRUE) +
  stat_pvalue_manual(pwt, step.increase = 0.1) +
  labs(
    caption = get_pwc_label(pwt)
  )
