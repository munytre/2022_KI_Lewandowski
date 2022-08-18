run_als_diffExp_analysis <- function(annot,data){
  m2h = One2One::ortholog_data_Mouse_Human$orthologs_one2one %>% dplyr::select(human.symbol,mouse.symbol) %>% dplyr::rename(HGNC.symbol=human.symbol,MGI.symbol=mouse.symbol)
  mod  = model.matrix(~annot$gender+factor(annot$status,levels=c("CTRL","ALS")))
  colnames(mod)[2:3] = c("Gender","ALS")
  tt_spinal = prep.tt(data,mod,m2h,"tt_HumanSpinal.csv",coef="ALS")
  return(tt_spinal)
}

