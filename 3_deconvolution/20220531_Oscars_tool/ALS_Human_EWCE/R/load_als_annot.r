load_als_annot <- function(path,data_spinal){
  setwd(path)
  annot_spinal_status = rep("ALS",dim(data_spinal)[2])
  annot_spinal_status[grep("CTRL",colnames(data_spinal))] = "CTRL"
  annot_spinal_gender = rep("Male",dim(data_spinal)[2])
  annot_spinal_gender[grep("_F",colnames(data_spinal))] = "Female"
  annot_spinal_age = as.numeric(as.character(gsub("M|F","",gsub("No.*","",gsub(".*_","",colnames(data_spinal))))))
  annot_spinal = data.frame(status=annot_spinal_status,gender=annot_spinal_gender,age=annot_spinal_age)
  #annot_spinal_tissue = gsub("_.*","",colnames(data_spinal))
  return(annot_spinal)
}
