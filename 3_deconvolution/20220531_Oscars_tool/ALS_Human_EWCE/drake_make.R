#############################################################
## GSE18920 --- Splicing array data from human spinal cord ##
#############################################################
library(limma)
library(One2One)
library(tidyr)
library(EWCE)
library(R.utils)
library(drake)
library(ggplot2)
library(reshape)
sourceDirectory("R")
#source("/Users/natske/Datasets that are too large to store elsewhere/SOD1 Spinal Cord (GSE18597)/generate.bootstrap.plots.for.transcriptome_SOD1.r")

dir.create("Results",showWarnings = FALSE)
dir.create("Results/Tables",showWarnings = FALSE)
dir.create("Results/Figures",showWarnings = FALSE)

load(file="Data/Tidy/ctd_OligosNCortex_(woLev)_thresh(0)_trim(0).rda")

plan <- drake_plan(
  path = "20220531_Oscars_tool/ALS_Human_EWCE/",
  data_spinal = load_als_data(path),
  annot_spinal = load_als_annot(path,data_spinal),
  tt_spinal = run_als_diffExp_analysis(annot_spinal,data_spinal),
  ewce_spinal = run_ewce(ctd,tt_spinal),
  catchOut = generate.bootstrap.plots.for.transcriptome(sct_data=ctd,tt=tt_spinal,thresh=250,annotLevel=1,reps=1000,full_results=ewce_spinal,listFileName="Human ALS Spinal")
)

config <- drake_config(plan)
vis_drake_graph(config)

make(plan)

### Extract expression data ----------------------------------------------------
# For mean expression data (single cell)
mean_ss <- ctd[[1]][["mean_exp"]]

# For bulk expression data (Anterior Horn)
library(tidyverse)
loadd(ctd,ewce_spinal,tt_spinal)
data_spinal = load_als_data("20220531_Oscars_tool/ALS_Human_EWCE/")
name_filter <- tt_spinal[,c(1,8)]
data_spinal$HGNC.symbol <- rownames(data_spinal)
bulk_AH <- left_join(name_filter,
                     data_spinal,
                     by = "HGNC.symbol")
rownames(bulk_AH) <- bulk_AH$MGI.symbol
tt_AH <- bulk_AH[3:22]
