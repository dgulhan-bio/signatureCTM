library('tm')
library('topicmodels')
library('ggplot2')

#docs <- Corpus(DirSource("/Users/dgulhan/Park/CTM/workingDir/documentTest"))
#table<-read.table('testrep.txt')
load('frame3BaseAging.Rda')
ds <- DataframeSource(frame3BaseAging)
corpus <- Corpus(ds)

vocab_list <- read.table('vocabSorted.txt')
vocab_mtrx<-array(unlist(vocab_list), dim = c(nrow(vocab_list), ncol(vocab_list), length(vocab_list)))

dtm<-DocumentTermMatrix(corpus)
inspect(dtm)

#dtm<-DocumentTermMatrix(docs, list(dictionary = vocab_mtrx))

signatures_min <- matrix(1, 96, 5)
exposures_min <- matrix(1, 5, 127)

min = 1

signatures_nmf<-read.csv('nmfresult.dat')

for(i in 1:1000){
  
  ctm_tmp <- CTM(dtm, 5, method = "VEM")

  beta<-slot(ctm_tmp, 'beta')
  signature_mtrx <- exp(t(beta))
  signatures_df<-data.frame(signature_mtrx)
  rownames(signatures_df)<-vocab_mtrx
  names(signature_mtrx)<-vocab_mtrx

  gamma<-slot(ctm_tmp, 'gamma')
  exposure_mtrx<-t(gamma)
  exposures_df<-data.frame(exposure_mtrx)
  terms<-slot(ctm_tmp, 'terms')
 
  reco <- signature_mtrx %*% exposure_mtrx
  original<-matrix(0,96,127)
  matrix_dtm<-as.matrix(dtm)

  for (j in 1:127){
    original[,j] = matrix_dtm[j,]
  }
  
  for(j in 1:127){
    original[,j] = original[,j]/sum(original[,j])
  }

#  for (k in 1:5){
#    print(k)
#    png(sprintf('bars_num%d_totalnum%d.png',k,5), width=6000, height=1500, res=300)

#    barplot(signature_mtrx[,k], col=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6), las=2, ylim=c(0, 0.25),names.arg = vocab_mtrx)
#    dev.off()
# }
  diff <- original - reco
  error <-sum(diag(diff %*% t(diff)))
  inspect(dtm)
 
  
  if(error<min){
    min = error
    signatures_min = signature_mtrx
    exposures_min = exposure_mtrx
  }
  #print(terms)
}


for (k in 1:5){
  print(k)
  png(sprintf('bars_num%d_totalnum%d.png',k,5), width=6000, height=1500, res=300)

  barplot(signature_mtrx[,k], col=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6), las=2, ylim=c(0, 0.25),names.arg = vocab_mtrx)
  dev.off()
}

save(original, reco, error, signatures_min, exposures_min ,file="1000run.Rda")

print(error)
 