###--------------------------------------------------------------------------###
### Project: ApoE in ALS                                                       -
### Purpose: PRJAN512012 deconvolution data extraction                         -
### Authors: QC Lin                                                            -
### Last Edit: 2022-08-12, QC Lin                                              -
### Contact: q.c.lin@students.uu.nl or sebastian.lewandowski@ki.se             -
###--------------------------------------------------------------------------###

### Dependencies and general settings ------------------------------------------

# Load libraries
library(GEOquery)
library(tidyverse)
# Set seed
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
###Metadata --------------------------------------------------------------------
# Load metadata (1)
gse <- getGEO(filename = "20220301_QCL_GSE124439/metadata/GSE124439_series_matrix.txt")
metadata <- gse@phenoData@data
colnames(metadata) <- sub(".", "_", colnames(metadata), fixed=TRUE)
colnames(metadata) <- sub(":", "_", colnames(metadata), fixed=TRUE)
colnames(metadata) <- sub(" ", "_", colnames(metadata), fixed=TRUE)

# Load metadata (2)
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

# Load metadata (3)
SRR_sample <- read.delim("20220301_QCL_GSE124439/metadata/PRJNA512012_tsv.txt")
SRR_sample <- SRR_sample[,c(2,4)]
colnames(SRR_sample) <- c("SRR", "title")

metadata_SRR_included <- left_join(SRR_sample,
                                   metadata,
                                   by = "title")

# Load genotype data
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

###List of sex-linked genes ----------------------------------------------------
x_y_genes <- scan(file = "20220301_QCL_GSE124439/sex_linked_genes/XY_ENS.txt",
                  what = "",
                  quiet = TRUE)

###Count files -----------------------------------------------------------------
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

###Matrix ----------------------------------------------------------------------
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

###Pre-filtering (Sample Specific) ---------------------------------------------
# Only ALS (no technical outlier)
count_matrix_ALS <- count_matrix[,(colnames(count_matrix) %in% metadata_SRR_included$SRR[metadata_SRR_included$Subject_Group == "ALS Spectrum MND"])]
metadata_SRR_included_ALS <- metadata_SRR_included[(rownames(metadata_SRR_included) %in% metadata_SRR_included$SRR[metadata_SRR_included$Subject_Group == "ALS Spectrum MND"]),]
#Check if order and naming is still the same
test_match_order(colnames(count_matrix_ALS),
                 metadata_SRR_included_ALS$SRR)
all(rownames(metadata_SRR_included_ALS) == colnames(count_matrix_ALS))

###Running DESeq2 (ALS vs NNC) -------------------------------------------------
## Only ALS - Motor Cortex
count_matrix_ALS_MC <- count_matrix_ALS[,!colnames(count_matrix_ALS) %in% metadata_SRR_included_ALS$SRR[metadata_SRR_included_ALS$Tissue == "Frontal Cortex"]] # & !is.na(metadata_SRR_included_ALS$Age_Onset) & metadata_SRR_included_ALS$Age_Onset <= 70 & metadata_SRR_included_ALS$Age_Onset >= 40]
metadata_SRR_included_ALS_MC <- metadata_SRR_included_ALS[!rownames(metadata_SRR_included_ALS) %in% metadata_SRR_included_ALS$SRR[metadata_SRR_included_ALS$Tissue == "Frontal Cortex"],] # & !is.na(metadata_SRR_included_ALS$Age_Onset) & metadata_SRR_included_ALS$Age_Onset <= 70 & metadata_SRR_included_ALS$Age_Onset >= 40,]
#Check if order and naming is still the same
test_match_order(colnames(count_matrix_ALS_MC),
                 metadata_SRR_included_ALS_MC$SRR)
all(rownames(metadata_SRR_included_ALS_MC) == colnames(count_matrix_ALS_MC))

table(metadata_SRR_included_ALS_MC$APOE)

## Only ALS - Frontal Cortex
count_matrix_ALS_FC <- count_matrix_ALS[,colnames(count_matrix_ALS) %in% metadata_SRR_included_ALS$SRR[metadata_SRR_included_ALS$Tissue == "Frontal Cortex"]] # & !is.na(metadata_SRR_included_ALS$Age_Onset) & metadata_SRR_included_ALS$Age_Onset <= 70 & metadata_SRR_included_ALS$Age_Onset >= 40]
metadata_SRR_included_ALS_FC <- metadata_SRR_included_ALS[rownames(metadata_SRR_included_ALS) %in% metadata_SRR_included_ALS$SRR[metadata_SRR_included_ALS$Tissue == "Frontal Cortex"],] # & !is.na(metadata_SRR_included_ALS$Age_Onset) & metadata_SRR_included_ALS$Age_Onset <= 70 & metadata_SRR_included_ALS$Age_Onset >= 40,]
#Check if order and naming is still the same
test_match_order(colnames(count_matrix_ALS_FC),
                 metadata_SRR_included_ALS_FC$SRR)
all(rownames(metadata_SRR_included_ALS_FC) == colnames(count_matrix_ALS_FC))

table(metadata_SRR_included_ALS_FC$APOE)

### For deconvolution ####
deconvo_matrix <- as.data.frame(count_matrix_ALS_MC)
deconvo_matrix$initial_alias <- rownames(deconvo_matrix)
# write.table(rownames(deconvo_matrix),
#             file = "20220531_Oscars_tool/ENS_DECONVO.txt",
#             quote = F,
#             row.names = F,
#             col.names = F)
ENS_GENEID_GPROFILER <- read.csv("20220531_Oscars_tool/resources/ENS_GENEID_GPROFILER.csv")[,c(1,3)]
deconvo_matrix_withGN <- inner_join(deconvo_matrix, 
                                    ENS_GENEID_GPROFILER,
                                    by = "initial_alias")
deconvo_matrix_withGN <- filter(deconvo_matrix_withGN,
                                name != "None") 
One2One_orthos <- read.delim("20220531_Oscars_tool/resources/One2One_orthos.csv",
                             header=FALSE) %>%
  `colnames<-`(c("name","MGI.symbol"))
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

###DECONVO ENDS###