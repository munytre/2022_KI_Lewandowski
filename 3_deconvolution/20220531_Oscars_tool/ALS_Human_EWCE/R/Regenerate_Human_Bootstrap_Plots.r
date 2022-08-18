regenerate_human_bootstrap_plots <- function(x){
  library(EWCE)
  load("/Users/natske/Datasets that are too large to store elsewhere/Disease Transcriptomes/ALS/Data/Tidy/ctd_OligosNCortex_(woLev)_thresh(0)_trim(0).rda")
  #setwd("/Users/natske/Datasets that are too large to store elsewhere/Disease Transcriptomes/ALS")
  tt_spinal = read.csv(file="Results/Tables/tt_HumanSpinal.csv")
  #load(ewce_spinal,file="EWCE_SPINAL.Rda")
  #load(file="Results/Tables/EWCE_SPINAL.Rda")
  ewce_spinal = read.csv(file="Results/Tables/EWCE_SPINAL.csv")
  #write.csv(ewce_spinal$joint_results,file="EWCE_SPINAL.csv")
  pdf(file="Fig_EWCE_HumanSpinal.pdf",width=6,height=5)
  ewce.plot(ewce_spinal$joint_results)
  dev.off()
  
  #source("/Users/natske/Google Drive/EWCE/R/generate.bootstrap.plots.for.transcriptome.r")
  source("/Users/natske/Datasets that are too large to store elsewhere/SOD1 Spinal Cord (GSE18597)/generate.bootstrap.plots.for.transcriptome_SOD1.r")
  #load("/Users/natske/Datasets that are too large to store elsewhere/Disease Transcriptomes/ALS/EWCE_SPINAL.Rda")
  generate.bootstrap.plots.for.transcriptome(sct_data=celltype_data,tt=tt_spinal,thresh=250,annotLevel=1,reps=1000,full_results=ewce_spinal,listFileName="Human ALS Spinal")
}