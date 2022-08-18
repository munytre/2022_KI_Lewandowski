###--------------------------------------------------------------------------###
### Project: ApoE in ALS                                                       -
### Purpose: Genotype assignment                                               -
### Authors: QC Lin                                                            -
### Last Edit: 2022-08-15, QC Lin                                              -
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

###GENERAL ----
###Functions -------------------------------------------------------------------
plotPCA.custom <- function (object, intgroup = "condition", ntop = 500, returnData = FALSE) {
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
test_match_order <- function(x,y) {
  
  if (isTRUE(all.equal(x,y))) print('Perfect match in same order')
  
  if (!isTRUE(all.equal(x,y)) && isTRUE(all.equal(sort(x),sort(y)))) print('Perfect match in wrong order')
  
  if (!isTRUE(all.equal(x,y)) && !isTRUE(all.equal(sort(x),sort(y)))) print('No match')
}

###DATA LOAD AND ASSEMBLY ----
### Metadata -------------------------------------------------------------------
# Load in metadata (1)
gse <- getGEO(filename = "20220301_QCL_GSE124439/metadata/GSE124439_series_matrix.txt")
metadata <- gse@phenoData@data
colnames(metadata) <- sub(".", "_", colnames(metadata), fixed=TRUE)
colnames(metadata) <- sub(":", "_", colnames(metadata), fixed=TRUE)
colnames(metadata) <- sub(" ", "_", colnames(metadata), fixed=TRUE)

# Load in metadata (2)
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

# Load in metadata (3)
SRR_sample <- read.delim("20220301_QCL_GSE124439/metadata/PRJNA512012_tsv.txt")
SRR_sample <- SRR_sample[,c(2,4)]
colnames(SRR_sample) <- c("SRR", "title")

metadata_SRR_included <- left_join(SRR_sample,
                                  metadata,
                                  by = "title")

# Load in genotype data
genotype_PRJNA512012 <- read.csv("20220314_QCL_Genes_APOE/Final/VCF_output_PRJNA512012.txt",
                                 sep="")

pre_metadata_SRR_included <- left_join(metadata_SRR_included,
                                   genotype_PRJNA512012,
                                   by = c("SRR" = "sample"))

# Rename columns for readability
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

# Reformat data to correct format
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
# write.table(metadata_SRR_included,
#             file = "20220531_Oscars_tool/ALS_Human_EWCE/Data/Metadata_PRJNA512012.csv",
#             quote = F,
#             sep = "\t",
#             row.names = F)

### List of sex-linked genes ---------------------------------------------------
x_y_genes <- scan(file = "20220301_QCL_GSE124439/sex_linked_genes/XY_ENS.txt",
                  what = "",
                  quiet = TRUE)

### Count files ----------------------------------------------------------------
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

### Matrix ---------------------------------------------------------------------
count_matrix <- counts %>%
  purrr::reduce(dplyr::left_join, by = "gene")

count_matrix_nonfiltered <- count_matrix

rownames(count_matrix_nonfiltered) <- count_matrix_nonfiltered[,1]
count_matrix_nonfiltered <- count_matrix_nonfiltered[,-1]
colnames(count_matrix_nonfiltered) <- count_names

count_matrix_nonfiltered <- as.matrix(count_matrix_nonfiltered)

#REMEMBER, DO NOT TO USE "gene" with filter
count_matrix <- filter(count_matrix,
                       !(gene %in% x_y_genes))

rownames(count_matrix) <- count_matrix[,1]
count_matrix <- count_matrix[,-1]
colnames(count_matrix) <- count_names
count_matrix <- as.matrix(count_matrix)
rownames(metadata_SRR_included) <- metadata_SRR_included$SRR

other_neuro_disorders <- rownames(dplyr::filter(metadata_SRR_included, Subject_Group %in% "Other Neurological Disorders"))

### Original matrix and metadata -----------------------------------------------
# As back-up
count_matrix_combined_ctrl <- count_matrix
metadata_SRR_included_combined_ctrl <- metadata_SRR_included

### Pre-Filtering (Outliers) ---------------------------------------------------
#Removal of category Other Neurological Disorders
count_matrix_nonfiltered_ALSNNC <- count_matrix_nonfiltered[,!(colnames(count_matrix_nonfiltered) %in% other_neuro_disorders)]
metadata_SRR_included_nonfiltered_ALSNNC <- metadata_SRR_included[!(rownames(metadata_SRR_included) %in% other_neuro_disorders),]

#Check if order and naming is still the same
test_match_order(colnames(count_matrix_nonfiltered_ALSNNC),
                 metadata_SRR_included_nonfiltered_ALSNNC$SRR)
all(rownames(metadata_SRR_included_nonfiltered_ALSNNC) == colnames(count_matrix_nonfiltered_ALSNNC))

#DESeq2 (ALS vs NNC - To find technical problematic sample)
dds <- DESeqDataSetFromMatrix(countData = count_matrix_nonfiltered_ALSNNC,
                              colData = metadata_SRR_included_nonfiltered_ALSNNC,
                              design = ~ Subject_Group) %>%
  DESeq()

#Variance Stabilizing Transformation pre-PCA
vsd <- vst(dds, blind=TRUE)
#PCAs with variables of interests (PC12)
plotPCA(vsd,
        intgroup="Tissue")
plotPCA(vsd,
        intgroup="Subject_Group")
plotPCA(vsd,
        intgroup="Gender")
plotPCA(vsd,
        intgroup="Subtype")
plotPCA(vsd,
        intgroup="APOE")
plotPCA(vsd,
        intgroup="APOE4_State")
#PCAs with variables of interests (PC2)
plotPCA.custom(vsd,
               intgroup="Gender")
plotPCA.custom(vsd,
               intgroup="APOE")
plotPCA.custom(vsd,
               intgroup="Tissue")
plotPCA.custom(vsd,
               intgroup="APOE4_State")
#Extract the outlier
outlier.PCA <- plotPCA.custom(vsd,
                              intgroup="Gender")
outlier_sample <- rownames(outlier.PCA$data[outlier.PCA$data$Gender == "Female" & outlier.PCA$data$PC3 > 0,])

# Remove possible technical outlier
count_matrix_technical <- count_matrix[,!(colnames(count_matrix) %in% outlier_sample)]
metadata_SRR_included_technical <- metadata_SRR_included[!(rownames(metadata_SRR_included) %in% outlier_sample),]
test_match_order(colnames(count_matrix_technical),
                 metadata_SRR_included_technical$SRR)
all(rownames(metadata_SRR_included_technical) == colnames(count_matrix_technical))

#Results
res_technical <- results(dds,
                      contrast = c("Subject_Group", "ALS Spectrum MND", "Non-Neurological Control"))
resOrdered_technical <- res_technical[order(res_technical$padj),]
summary(resOrdered_technical)

### Pre-filtering (Sample Specific) --------------------------------------------
# Only ALS (no technical outlier)
count_matrix_ALS <- count_matrix[,(colnames(count_matrix) %in% metadata_SRR_included$SRR[metadata_SRR_included$Subject_Group == "ALS Spectrum MND"])]
metadata_SRR_included_ALS <- metadata_SRR_included[(rownames(metadata_SRR_included) %in% metadata_SRR_included$SRR[metadata_SRR_included$Subject_Group == "ALS Spectrum MND"]),]
#Check if order and naming is still the same
test_match_order(colnames(count_matrix_ALS),
                 metadata_SRR_included_ALS$SRR)
all(rownames(metadata_SRR_included_ALS) == colnames(count_matrix_ALS))

## Only ALS - Frontal Cortex (no NAs)
count_matrix_ALS_FC <- count_matrix_ALS[,colnames(count_matrix_ALS) %in% metadata_SRR_included_ALS$SRR[metadata_SRR_included_ALS$Tissue == "Frontal Cortex"]] # & !is.na(metadata_SRR_included_ALS$Age_Onset) & metadata_SRR_included_ALS$Age_Onset <= 70 & metadata_SRR_included_ALS$Age_Onset >= 60]
metadata_SRR_included_ALS_FC <- metadata_SRR_included_ALS[rownames(metadata_SRR_included_ALS) %in% metadata_SRR_included_ALS$SRR[metadata_SRR_included_ALS$Tissue == "Frontal Cortex"],] # & !is.na(metadata_SRR_included_ALS$Age_Onset) & metadata_SRR_included_ALS$Age_Onset <= 70 & metadata_SRR_included_ALS$Age_Onset >= 60,]
#Check if order and naming is still the same
test_match_order(colnames(count_matrix_ALS_FC),
                 metadata_SRR_included_ALS_FC$SRR)
all(rownames(metadata_SRR_included_ALS_FC) == colnames(count_matrix_ALS_FC))

## Only ALS - Motor Cortex
count_matrix_ALS_MC <- count_matrix_ALS[,!(colnames(count_matrix_ALS) %in% metadata_SRR_included_ALS$SRR[metadata_SRR_included_ALS$Tissue == "Frontal Cortex"])]
metadata_SRR_included_ALS_MC <- metadata_SRR_included_ALS[!(rownames(metadata_SRR_included_ALS) %in% metadata_SRR_included_ALS$SRR[metadata_SRR_included_ALS$Tissue == "Frontal Cortex"]),]
#Check if order and naming is still the same
test_match_order(colnames(count_matrix_ALS_MC),
                 metadata_SRR_included_ALS_MC$SRR)
all(rownames(metadata_SRR_included_ALS_MC) == colnames(count_matrix_ALS_MC))

## Only ALS - Motor Cortex (Lateral & no NAs)
count_matrix_ALS_MC_lateral <- count_matrix_ALS_MC[,colnames(count_matrix_ALS_MC) %in% metadata_SRR_included_ALS_MC$SRR[metadata_SRR_included_ALS_MC$Tissue == "Motor Cortex (Lateral)" & !is.na(metadata_SRR_included_ALS_MC$Age_Onset)]]
metadata_SRR_included_ALS_MC_lateral <- metadata_SRR_included_ALS_MC[rownames(metadata_SRR_included_ALS_MC) %in% metadata_SRR_included_ALS_MC$SRR[metadata_SRR_included_ALS_MC$Tissue == "Motor Cortex (Lateral)" & !is.na(metadata_SRR_included_ALS_MC$Age_Onset)],]
#Check if order and naming is still the same
test_match_order(colnames(count_matrix_ALS_MC_lateral),
                 metadata_SRR_included_ALS_MC_lateral$SRR)
all(rownames(metadata_SRR_included_ALS_MC_lateral) == colnames(count_matrix_ALS_MC_lateral))

## Only Motor Cortex
count_matrix_MC <- count_matrix_technical[,!(colnames(count_matrix_technical) %in% metadata_SRR_included_technical$SRR[metadata_SRR_included_technical$Tissue == "Frontal Cortex"])]
metadata_SRR_included_MC <- metadata_SRR_included_technical[!(rownames(metadata_SRR_included_technical) %in% metadata_SRR_included_technical$SRR[metadata_SRR_included_technical$Tissue == "Frontal Cortex"]),]
#Check if order and naming is still the same
test_match_order(colnames(count_matrix_MC),
                 metadata_SRR_included_MC$SRR)
all(rownames(metadata_SRR_included_MC) == colnames(count_matrix_MC))

###DESEQ2 ----
### Running DESeq2 (ALS vs NNC - ALS Only) -------------------------------------
dds_ALS <- DESeqDataSetFromMatrix(countData = count_matrix_ALS,
                                  colData = metadata_SRR_included_ALS,
                                  design = ~ APOE) %>%
  DESeq()

table(metadata_SRR_included_ALS$APOE)

vsd_ALS <- vst(dds_ALS, blind=TRUE)
plotPCA(vsd_ALS,
        intgroup="Subject_Group")
plotPCA(vsd_ALS,
        intgroup="Gender")
plotPCA(vsd_ALS,
        intgroup="Tissue")
plotPCA(vsd_ALS,
        intgroup="Subtype")
plotPCA(vsd_ALS,
        intgroup="APOE")
plotPCA(vsd_ALS,
        intgroup="APOE4_State")
plotPCA.custom(vsd_ALS,
               intgroup="Gender")
plotPCA.custom(vsd_ALS,
               intgroup="APOE")
plotPCA.custom(vsd_ALS,
               intgroup="APOE4_State")

res_ALS <- results(dds_ALS,
                   contrast = c("APOE", "44", "33"))
resOrdered_ALS <- res_ALS[order(res_ALS$padj),]
summary(res_ALS)
#Subset the gene names that are significant and have a log2FC > 1 or -1
res_de_genes_ALS_DOWN <- rownames(subset(resOrdered_ALS,
                                         resOrdered_ALS$padj < 0.1 & resOrdered_ALS$log2FoldChange < 0))
res_de_genes_ALS_UP <- rownames(subset(resOrdered_ALS,
                                         resOrdered_ALS$padj < 0.1 & resOrdered_ALS$log2FoldChange > 0))
resSign_ALS_DOWN <- as.data.frame(subset(resOrdered_ALS)[res_de_genes_ALS_DOWN,])
resSign_ALS_UP <- as.data.frame(subset(resOrdered_ALS)[res_de_genes_ALS_UP,])
# write.table(res_de_genes_ALS, file = "ALS.txt", row.names = F, col.names = F, quote = F)

### Running DESeq2 (ALS vs NNC - ALS Only_FC) ----------------------------------
ggplot(metadata_SRR_included_ALS_FC,
       aes(x=Age_Onset)) +
  geom_density(binwidth=1)

dds_ALS_FC <- DESeqDataSetFromMatrix(countData = count_matrix_ALS_FC,
                                     colData = metadata_SRR_included_ALS_FC,
                                     design = ~ APOE) %>%
  DESeq()

table(metadata_SRR_included_ALS_FC$APOE)

#PCA
vsd_ALS_FC <- vst(dds_ALS_FC, blind=TRUE)
plotPCA(vsd_ALS_FC,
        intgroup="Subject_Group")
plotPCA(vsd_ALS_FC,
        intgroup="Gender")
plotPCA(vsd_ALS_FC,
        intgroup="Tissue")
plotPCA(vsd_ALS_FC,
        intgroup="Subtype")
plotPCA(vsd_ALS_FC,
        intgroup="APOE")
plotPCA(vsd_ALS_FC,
        intgroup="APOE4_State")
plotPCA.custom(vsd_ALS_FC,
               intgroup="Gender")
plotPCA.custom(vsd_ALS_FC,
               intgroup="APOE")
plotPCA.custom(vsd_ALS_FC,
               intgroup="APOE4_State")
#Results
res_ALS_FC <- results(dds_ALS_FC,
                      contrast = c("APOE", "34", "33"))
resOrdered_ALS_FC <- res_ALS_FC[order(res_ALS_FC$padj),]
summary(resOrdered_ALS_FC)
#Subset the gene names that are significant and have a log2FC > 1 or -1
res_de_genes_ALS_FC_DOWN <- rownames(subset(resOrdered_ALS_FC,
                                       resOrdered_ALS_FC$padj < 0.1 & resOrdered_ALS_FC$log2FoldChange < 0))
res_de_genes_ALS_FC_UP <- rownames(subset(resOrdered_ALS_FC,
                                            resOrdered_ALS_FC$padj < 0.1 & resOrdered_ALS_FC$log2FoldChange > 0))
resSign_ALS_FC_DOWN <- as.data.frame(subset(resOrdered_ALS_FC)[res_de_genes_ALS_FC_DOWN,])
resSign_ALS_FC_UP <- as.data.frame(subset(resOrdered_ALS_FC)[res_de_genes_ALS_FC_UP,])
# write.table(resSign, file = "DESeq2_tables/FC_APOE_AoObinned_44vs33_6070.txt", quote = FALSE)

### Running DESeq2 (ALS vs NNC - ALS Only_MC) ----------------------------------
## Only ALS - Motor Cortex
count_matrix_ALS_MC <- count_matrix_ALS[,!colnames(count_matrix_ALS) %in% metadata_SRR_included_ALS$SRR[metadata_SRR_included_ALS$Tissue == "Frontal Cortex"]] # & !is.na(metadata_SRR_included_ALS$Age_Onset) & metadata_SRR_included_ALS$Age_Onset <= 70 & metadata_SRR_included_ALS$Age_Onset >= 40]
metadata_SRR_included_ALS_MC <- metadata_SRR_included_ALS[!rownames(metadata_SRR_included_ALS) %in% metadata_SRR_included_ALS$SRR[metadata_SRR_included_ALS$Tissue == "Frontal Cortex"],] # & !is.na(metadata_SRR_included_ALS$Age_Onset) & metadata_SRR_included_ALS$Age_Onset <= 70 & metadata_SRR_included_ALS$Age_Onset >= 40,]
#Check if order and naming is still the same
test_match_order(colnames(count_matrix_ALS_MC),
                 metadata_SRR_included_ALS_MC$SRR)
all(rownames(metadata_SRR_included_ALS_MC) == colnames(count_matrix_ALS_MC))

table(metadata_SRR_included_ALS_MC$APOE)

dds_ALS_MC <- DESeqDataSetFromMatrix(countData = count_matrix_ALS_MC,
                                     colData = metadata_SRR_included_ALS_MC,
                                     design = ~ APOE) %>%
  DESeq()

vsd_ALS_MC <- vst(dds_ALS_MC, blind=TRUE)
plotPCA(vsd_ALS_MC,
        intgroup="Subject_Group")
plotPCA(vsd_ALS_MC,
        intgroup="Tissue")
plotPCA(vsd_ALS_MC,
        intgroup="Gender")
plotPCA(vsd_ALS_MC,
        intgroup="Subtype")
plotPCA(vsd_ALS_MC,
        intgroup="APOE")
plotPCA(vsd_ALS_MC,
        intgroup="APOE4_State")
plotPCA(vsd_ALS_MC,
        intgroup="Disease_Duration")
plotPCA(vsd_ALS_MC,
        intgroup="Age_Death")
plotPCA.custom(vsd_ALS_MC,
               intgroup="Gender")
plotPCA.custom(vsd_ALS_MC,
               intgroup="APOE")
plotPCA.custom(vsd_ALS_MC,
               intgroup="APOE4_State")

#Results
res_ALS_MC <- results(dds_ALS_MC,
                      contrast = c("APOE", "34", "33"))
resOrdered_ALS_MC <- res_ALS_MC[order(res_ALS_MC$padj),]
summary(resOrdered_ALS_MC)
#Subset the gene names that are significant and have a log2FC > 1 or -1
res_de_genes_ALS_MC_DOWN <- rownames(subset(resOrdered_ALS_MC,
                                            resOrdered_ALS_MC$padj < 0.1 & resOrdered_ALS_MC$log2FoldChange < 0))
res_de_genes_ALS_MC_UP <- rownames(subset(resOrdered_ALS_MC,
                                          resOrdered_ALS_MC$padj < 0.1 & resOrdered_ALS_MC$log2FoldChange > 0))
resSign_ALS_MC_DOWN <- as.data.frame(subset(resOrdered_ALS_MC)[res_de_genes_ALS_MC_DOWN,])
resSign_ALS_MC_UP <- as.data.frame(subset(resOrdered_ALS_MC)[res_de_genes_ALS_MC_UP,])
# write.table(resSign, file = "DESeq2_tables/MC_APOE_34vs33_4070_up01_log1.txt", quote = FALSE)

###FOR DECONVOLUTION### ----
### Make matrix
deconvo_matrix <- as.data.frame(count_matrix_ALS_MC)
deconvo_matrix$initial_alias <- rownames(deconvo_matrix)
# write.table(rownames(deconvo_matrix),
#             file = "20220531_Oscars_tool/ENS_DECONVO.txt",
#             quote = F,
#             row.names = F,
#             col.names = F)

# Filter orthologs
ENS_GENEID_GPROFILER <- read.csv("20220531_Oscars_tool/ENS_GENEID_GPROFILER.csv")[,c(1,3)]
deconvo_matrix_withGN <- inner_join(deconvo_matrix, 
                                    ENS_GENEID_GPROFILER,
                                    by = "initial_alias")
deconvo_matrix_withGN <- filter(deconvo_matrix_withGN,
                                name != "None") 
One2One_orthos <- read.delim("20220531_Oscars_tool/One2One_orthos.csv",
                             header=FALSE) %>%
  `colnames<-`(c("name","MGI.symbol"))

# Final matrix
deconvo_final <- inner_join(One2One_orthos,
                            deconvo_matrix_withGN,
                            by = "name")
deconvo_final <- deconvo_final[!duplicated(deconvo_final[ , "name"]),]
rownames(deconvo_final) <- deconvo_final$MGI.symbol
deconvo_final <- deconvo_final[,-c(1,2,81)]
# write.table(deconvo_final,
#             file = "20220531_Oscars_tool/ALS_Human_EWCE/Data/tt_PRJNA512012_MC_ALS.csv",
#             quote = F,
#             sep = "\t")

# write.table(metadata_SRR_included_ALS_MC,
#             file = "20220531_Oscars_tool/ALS_Human_EWCE/Data/Metadata_PRJNA512012_ALS_MC.csv",
#             quote = F,
#             sep = "\t",
#             row.names = F)

###DECONVOLUTION EXTRACTION ENDS###

###DESEQ2 (Extra) ----
### Running DESeq2 (ALS vs NNC - ALS Only_MC) ----------------------------------
dds_MC <- DESeqDataSetFromMatrix(countData = count_matrix_MC,
                                     colData = metadata_SRR_included_MC,
                                     design = ~ APOE4_State) %>%
  DESeq()

vsd_MC <- vst(dds_MC, blind=TRUE)
plotPCA(vsd_MC,
        intgroup="Subject_Group")
plotPCA(vsd_MC,
        intgroup="Gender")
plotPCA(vsd_MC,
        intgroup="Subtype")
plotPCA(vsd_MC,
        intgroup="APOE")
plotPCA(vsd_MC,
        intgroup="APOE4_State")
plotPCA(vsd_MC,
        intgroup="Disease_Duration")
plotPCA(vsd_MC,
        intgroup="Age_Death")
plotPCA.custom(vsd_MC,
               intgroup="Gender")
plotPCA.custom(vsd_MC,
               intgroup="APOE")
plotPCA.custom(vsd_MC,
               intgroup="APOE4_State")
t1 <-plotPCA.custom_rainbow(vsd_MC,
                       intgroup="Disease_Duration")
t2 <- t1$data[!(rownames(t1$data) %in% c("SRR8375314","SRR8375401","SRR8375289")),]

t3 <- t1
t3$data <- t2 
t3

ggplot(metadata_SRR_included,
       aes(x=Disease_Duration, y=APOE, color = Tissue)) +
  geom_point()

### Running DESeq2 (ALS vs NNC - ALS MC Lateral) -------------------------------
metadata_SRR_included_ALS_MC_lateral <- metadata_SRR_included_ALS_MC_lateral %>%
  mutate(agegroup = case_when(Age_Onset >= 55 & Age_Onset <= 75 ~ 'Mid',
                              Age_Onset < 55  ~ 'Low',
                              Age_Onset > 75 ~ 'High'))
metadata_SRR_included_ALS_MC_lateral <- metadata_SRR_included_ALS_MC_lateral %>%
  mutate(OnsetBulbar = case_when(Site_Onset == "Bulbar" ~ 'Yes',
                                 Site_Onset != "Bulbar" ~ 'No'))
dds_ALS_MC_lateral <- DESeqDataSetFromMatrix(countData = count_matrix_ALS_MC_lateral,
                                             colData = metadata_SRR_included_ALS_MC_lateral,
                                             design = ~ APOE + agegroup + Gender + OnsetBulbar) %>%
  DESeq()

vsd_ALS_MC_lateral <- vst(dds_ALS_MC_lateral, blind=TRUE)
plotPCA(vsd_ALS_MC_lateral,
        intgroup="Subject_Group")
plotPCA(vsd_ALS_MC_lateral,
        intgroup="Gender")
plotPCA(vsd_ALS_MC_lateral,
        intgroup="Subtype")
plotPCA(vsd_ALS_MC_lateral,
        intgroup="APOE")
plotPCA(vsd_ALS_MC_lateral,
        intgroup="APOE4_State")
plotPCA.custom(vsd_ALS_MC_lateral,
               intgroup="Gender")
plotPCA.custom(vsd_ALS_MC_lateral,
               intgroup="APOE")
plotPCA.custom(vsd_ALS_MC_lateral,
               intgroup="APOE4_State")

res_ALS_MC_lateral <- results(dds_ALS_MC_lateral,
                              contrast = c("APOE", "44", "33"))
resOrdered_ALS_MC_lateral <- res_ALS_MC_lateral[order(res_ALS_MC_lateral$padj),]
summary(resOrdered_ALS_MC_lateral)
#Subset the gene names that are significant and have a log2FC > 1 or -1
res_de_genes_ALS_MC_lateral <- rownames(subset(resOrdered_ALS_MC_lateral,
                                               resOrdered_ALS_MC_lateral$padj < 0.1))
resSign <- as.data.frame(subset(resOrdered_ALS_MC_lateral)[res_de_genes_ALS_MC_lateral,])
# write.table(resSign, file = "DESeq2_tables/MC_lateral_APOE_AoO_44vs33.txt", quote = FALSE)

### Running DESeq2 (ALS vs NNC - technical removed) ----------------------------
dds_t <- DESeqDataSetFromMatrix(countData = count_matrix_technical,
                                colData = metadata_SRR_included_technical,
                                design = ~ Subject_Group) %>%
  DESeq()

vsd_t <- vst(dds_t, blind=TRUE)
plotPCA(vsd_t,
        intgroup="Subject_Group")
plotPCA(vsd_t,
        intgroup="Gender")
plotPCA(vsd_t,
        intgroup="Subtype")
plotPCA(vsd_t,
        intgroup="APOE")
plotPCA(vsd_t,
        intgroup="APOE4_State")
plotPCA.custom(vsd_t,
               intgroup="Gender")
plotPCA.custom(vsd_t,
               intgroup="APOE")
plotPCA.custom(vsd_t,
               intgroup="APOE4_State")

res_t <- results(dds_t,
                 contrast = c("Subject_Group", "ALS Spectrum MND", "Non-Neurological Control"))
resOrdered_t <- res_t[order(res_t$padj),]
#Subset the gene names that are significant and have a log2FC > 2 or -2
res_de_genes_t_up <- rownames(subset(resOrdered_t,
                                     resOrdered_t$padj < 0.1 & resOrdered_t$log2FoldChange > 0))
res_de_genes_t_down <- rownames(subset(resOrdered_t,
                                       resOrdered_t$padj < 0.1 & resOrdered_t$log2FoldChange < 0))
vsd_t_nonblind <- vst(dds_t, blind=FALSE)
summary(resOrdered_t)
resSigUp <- as.data.frame(subset(resOrdered_t)[res_de_genes_t_up,])
resSigDown <- as.data.frame(subset(resOrdered_t)[res_de_genes_t_down,])
rownames(resSigUp)

write.table(resSigUp, file = "testD", quote = FALSE)

### Extra code -----------------------------------------------------------------
# library(viridis)
# plotPCA.custom_rainbow <- function (object, intgroup = "condition", ntop = 500, returnData = FALSE) {
#   rv <- rowVars(assay(object))
#   select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
#                                                      length(rv)))]
#   pca <- prcomp(t(assay(object)[select, ]))
#   percentVar <- pca$sdev^2/sum(pca$sdev^2)
#   if (!all(intgroup %in% names(colData(object)))) {
#     stop("the argument 'intgroup' should specify columns of colData(dds)")
#   }
#   intgroup.df <- as.data.frame(colData(object)[, intgroup, 
#                                                drop = FALSE])
#   group <- if (length(intgroup) > 1) {
#     factor(apply(intgroup.df, 1, paste, collapse = " : "))
#   }
#   else {
#     colData(object)[[intgroup]]
#   }
#   d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3], PC4 = pca$x[, 4], group = group, 
#                   intgroup.df, name = colnames(object))
#   if (returnData) {
#     attr(d, "percentVar") <- percentVar[1:2]
#     return(d)
#   }
#   ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group")) + 
#     geom_point(size = 3) +
#     xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
#     ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
#     coord_fixed() +
#     scale_color_viridis(option = "viridis",
#                         discrete = FALSE)
# }
# plotPCA.custom_rainbow(vsd_ALS_MC,
#                        intgroup="Disease_Duration")