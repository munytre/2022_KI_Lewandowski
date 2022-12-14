regenerate_bootstrap_plots <- function(){
    load("/Users/natske/Datasets that are too large to store elsewhere/SOD1 Spinal Cord (GSE18597)/tt.rda")
    setwd("/Users/natske/Datasets that are too large to store elsewhere/SOD1 Spinal Cord (GSE18597)/")
    ages = c(28,42,56,70,98,112,126)
    tt_list = list(tt_28,tt_42,tt_56,tt_70,tt_98,tt_112,tt_126)
    source("/Users/natske/Datasets that are too large to store elsewhere/SOD1 Spinal Cord (GSE18597)/generate.bootstrap.plots.for.transcriptome_SOD1.r")
    load("/Users/natske/Datasets that are too large to store elsewhere/SOD1 Spinal Cord (GSE18597)/celltype_data_OligosNCortex_(woLev)_thresh(0)_trim(0).rda")
    
    
    #for(thresh in c(100,150,200,250,300,350)){
    for(thresh in c(250)){
        load(file=sprintf("SOD1_EWCE_RES_Thresh%s.rda",thresh))
        for(ageI in 1:length(ages)){
            print(sprintf("Age: %s",ages[ageI]))
            tt=tt_list[[ageI]]
            colnames(tt)[colnames(tt)=="ID"] = "MGI.symbol"
            generate.bootstrap.plots.for.transcriptome(sct_data=celltype_data,tt=tt,thresh=thresh,annotLevel=1,reps=1000,full_results=combined_ewce[[ageI]],listFileName=sprintf("%sdays_MOUSE",ages[ageI]))
        }
    }

}