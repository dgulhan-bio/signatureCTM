get_error_frob<-function(genomes, signatures, exposures){
  if(sum(exposures[,1]) == 1){
    for(i in 1:ncol(exposures)){
      exposures[,i] = sum(genomes[,i])*exposures[,i]
    }
  }
  reco <- signatures %*% exposures
  diff <- genomes - reco
  error <- sqrt(sum(diag(t(diff) %*% (diff))))
  error <- error/sqrt(sum(diag(t(genomes) %*% (genomes))))
  return(error)
}