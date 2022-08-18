run_ewce <- function(ctd,tt,thresh = 250){
  ewce_spinal = ewce_expression_data(ctd,annotLevel=1,tt=tt,sortBy="t",thresh=thresh,reps=10000,ttSpecies = "human",sctSpecies="mouse")
  ewce_spinal$joint_results$Q = p.adjust(ewce_spinal$joint_results$p)
  ewce_spinal$joint_results[order(ewce_spinal$joint_results$sd_from_mean),]
  
  save(ewce_spinal,file="Results/ewce_spinal.rda")
  write.csv(ewce_spinal$joint_results,file="Results/Tables/EWCE_SPINAL.csv")
  pdf(file="Results/Figures/Fig_EWCE_HumanSpinal.pdf",width=6,height=5)
  print(ewce_plot(ewce_spinal$joint_results))
  dev.off()
  
  return(ewce_spinal)
}