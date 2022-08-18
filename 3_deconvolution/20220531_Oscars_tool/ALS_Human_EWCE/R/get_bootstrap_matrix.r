get_bootstrap_matrix <- function(nReps,mouse.hits,mouse.bg,cell_types,sct_data,annotLevel){
  combinedGenes = unique(c(mouse.hits, mouse.bg))
  exp_mats = list()
  for(cc in cell_types){ # For each celltype...
    exp_mats[[cc]] = matrix(0,nrow=nReps,ncol=length(mouse.hits)) # Create an empty matrix, with a row for each bootstrap replicate
    rownames(exp_mats[[cc]]) = sprintf("Rep%s",1:nReps)
  }
  
  for(s in 1:nReps){
    bootstrap_set = sample(combinedGenes,length(mouse.hits))
    ValidGenes = rownames(sct_data[[annotLevel]]$specificity)[rownames(sct_data[[annotLevel]]$specificity) %in% bootstrap_set]
    
    expD = sct_data[[annotLevel]]$specificity[ValidGenes,]
    
    for(cc in cell_types){
      exp_mats[[cc]][s,] = sort(expD[,cc])
    }
  }
  return(exp_mats)
}