library('tm')
library('topicmodels')
library('ggplot2')
source('get_error_frob.R')

nsig_max = 10

#docs <- Corpus(DirSource("/Users/dgulhan/Park/CTM/workingDir/documentTest"))
#table <- read.table('testrep.txt')
load('craig/split_out_documents_2/frame3BaseAging.Rda')
ds <- DataframeSource(frame3BaseAging)
corpus <- Corpus(ds)

vocab_list <- read.table('vocabSorted2.txt')
vocab_mtrx <- array(unlist(vocab_list), dim = c(nrow(vocab_list), ncol(vocab_list), length(vocab_list)))

dtm <- DocumentTermMatrix(corpus)
inspect(dtm)

ntype <- dim(dtm)[[2]]
ngenome <- dim(dtm)[[1]]

total_columns <- 0
for(isig in 2:nsig_max){
  total_columns <- total_columns + isig
}

min_signatures_array = matrix(0, ntype, total_columns)
min_exposures_array = matrix(0, total_columns, ngenome)

total_columns <- 0

min_error = rep(nsig_max)

for(isig in 2:nsig_max){
 
  #dtm<-DocumentTermMatrix(docs, list(dictionary = vocab_mtrx))

  signatures_min <- matrix(1, ntype, isig)
  exposures_min <- matrix(1, isig, ngenome)

  min = 1
  for(i in 1:100){
  
    ctm_tmp <- CTM(dtm, isig, method = "VEM")

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
    matrix_dtm<-as.matrix(dtm)
    original<-matrix(0,dim(matrix_dtm)[[2]],dim(matrix_dtm)[[1]])
  
    for (j in 1:ngenome){
      original[,j] = matrix_dtm[j,]
    }

    originalNorm<-original
  
    for(j in 1:ngenome){
      originalNorm[,j] = originalNorm[,j]/sum(originalNorm[,j])
    }

    error <- getErrorFrob(original, signature_mtrx, exposure_mtrx)
  
    if(error<min){
      min = error
      print(min)
      signatures_min = signature_mtrx
      exposures_min = exposure_mtrx
    }
    #print(terms)
  }
 

  for (k in 1:isig){
    min_signatures_array[,total_columns + k] = signatures_min[, k]
    min_exposures_array[total_columns + k,] = exposures_min[k, ]
    min_error[k] = min
    print(k)
    png(sprintf('bars_num%d_totalnum%d.png',k,isig), width=6000, height=1500, res=300)

    barplot(signatures_min[,k], col=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6), las=2, ylim=c(0, 0.25),names.arg = vocab_mtrx)
    dev.off() 
  }
  total_columns <- total_columns + isig
}

save(min_signatures_array, originalNorm, original,  min_error, file="single_double_100run.Rda")

 