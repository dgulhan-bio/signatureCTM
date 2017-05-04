genomes_rep<-rep(0, 12192)
for(i in 1:127){
    first = as.integer((i-1)*96+1)
    last = as.integer(i*96)
    print(first)
    print(last)
    genomes_rep[first:last] <- t(original[,i])
}

signatures_rep<-matrix(0, 12192, 635)
for(i in 1:127){
    first = as.integer((i-1)*96+1)
    last = as.integer(i*96)
    print(first)
    print(last)
    signatures_rep[first:last,as.integer((i-1)*5+1):as.integer(i*5)] <- signature_mtrx
}

exposures_flat<-coef(nnls(signatures_rep, genomes_rep))

exposures_nnls <-matrix(0,5,127)
for(i in 1:127){
  exposures_nnls[,i] <- exposures_flat[as.integer((i-1)*5+1):as.integer(i*5)]
  exposures_nnls[,i]<-exposures_nnls[,i]/sum(exposures_nnls[,i])
}

barplot(exposures_nnls, col = c(3,4,5,6,7))

 genomes_nnls <-signature_mtrx %*% exposures_nnls
 diff_nnls <- original - genomes_nnls
 error_nnls<-sum(diag(diff_nnls %*% t(diff_nnls)))
