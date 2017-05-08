get_exposures_simultaneous<-function(){

  source('match_signatures.R')
  source('get_error_frob.R')
  library('nnls')
  load('100run.Rda')

  signatures_nmf<-as.matrix(read.csv('nmfprocesses.dat',header = FALSE))
  exposures_nmf<-as.matrix(read.csv('nmfexposures.dat',header = FALSE))
  genomes_nmf<-as.matrix(read.csv('originalgenomes.dat',header = FALSE))
  match_index<-matchSignatures(signatures_nmf, signatures_min)

  exposures_nmf_resorted<-exposures_nmf
  for(i in 1:5){
    exposures_nmf_resorted[i,] = exposures_nmf[match_index[[i]],]
  }

  for(i in 1:127){
    sum = sum(original[,i])
    exposures_min[,i] = sum*exposures_min[,i]  
  }

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
    signatures_rep[first:last,as.integer((i-1)*5+1):as.integer(i*5)] <- signatures_min
  }

  exposures_flat<-coef(nnls(signatures_rep, genomes_rep))

  exposures_nnls <-matrix(0,5,127)
  for(i in 1:127){
    exposures_nnls[,i] <- exposures_flat[as.integer((i-1)*5+1):as.integer(i*5)]
    #exposures_nnls[,i]<-exposures_nnls[,i]/sum(exposures_nnls[,i])
  }

  png('exposuresComparison.png', width=6000, height=1500, res=300)
  par(mfrow=c(1,3)) 
  barplot(exposures_nnls, col = c(3,4,5,6,7))
  mtext(text='decompose in ctm',side=3)
  barplot(exposures_min, col = c(3,4,5,6,7))
  mtext(text='ctm algo', side=3)
  barplot(exposures_nmf_resorted, col = c(3,4,5,6,7))
  mtext(text='nmf algo',side=3)
  dev.off()

  error_nnls<-getErrorFrob(original, signatures_min, exposures_nnls)
  error_ctm<-getErrorFrob(original, signatures_min, exposures_min)
  error_nmf<-getErrorFrob(genomes_nmf, signatures_nmf, exposures_nmf)

  print(sprintf('nnls %.3f ctm %.3f nmf %.3f', error_nnls, error_ctm, error_nmf))
}