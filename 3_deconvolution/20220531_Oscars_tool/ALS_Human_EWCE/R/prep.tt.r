prep.tt <- function(exp_data,mod,m2h,title,coef){
  fit = lmFit(exp_data,mod)
  eb = eBayes(fit)
  tt = topTable(eb, coef=coef, adjust="BH",number=1000000)	
  tt2 = cbind(tt,HGNC.symbol=rownames(tt))
  tt2$HGNC.symbol=as.character(tt2$HGNC.symbol)
  tt3 = merge(tt2,m2h,by="HGNC.symbol")
  tt4 = tt3[order(tt3$P.Value),]
  write.csv(tt4,file=sprintf("Results/Tables/%s",title))
  return(tt4)
}