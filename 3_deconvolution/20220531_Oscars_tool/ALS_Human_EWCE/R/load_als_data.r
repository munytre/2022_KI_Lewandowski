load_als_data <- function(path){
  setwd(path)
  data_spinal = read.csv(sprintf("%sData/Raw/GSE18920_extended_for_R_wGeneNames.txt",path),stringsAsFactors = FALSE,sep="\t")
  rownames(data_spinal) = gsub(" ","",data_spinal$NAME)
  data_spinal = data_spinal[,-1]
  data_spinal = data_spinal[,grep("AH",colnames(data_spinal))] #<--- THIS IS THE LINE I SHOULD HAVE HAD PREVIOUSLY
  return(data_spinal)
}