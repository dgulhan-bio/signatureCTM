rm(list=ls())
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

min_signatures_array <- matrix(0, ntype, total_columns)
clusters_array <- matrix(0, ntype, total_columns)
min_exposures_array <- matrix(0, total_columns, ngenome)
variance_all <- rep(0, total_columns)
total_columns <- 0

min_error = rep(nsig_max)

for(isig in 2:nsig_max){
 
  #dtm<-DocumentTermMatrix(docs, list(dictionary = vocab_mtrx))

  signatures_min <- matrix(1, ntype, isig)
  exposures_min <- matrix(1, isig, ngenome)

  min = 1
  all_signatures <- matrix(0, isig*100, ntype)
  
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

    error <- get_error_frob(original, signature_mtrx, exposure_mtrx)
    dim(all_signatures)
    dim(signature_mtrx)
    all_signatures[integer((i-1)*isig+1):(i*isig), ] <- t(signature_mtrx)
    
    if(error<min){
      min = error
      print(min)
      signatures_min = signature_mtrx
      exposures_min = exposure_mtrx
    }
    #print(terms)
  } 
   
  fit <- kmeans(all_signatures, isig) # 5 cluster solution
  aggregate(all_signatures, by=list(fit$cluster), FUN=mean)
  all_signatures <- data.frame(all_signatures, fit$cluster)
  
  num_entry_in_cluster = rep(0,isig)
  mean_clusters = matrix(0, ntype, isig) 
  variance_clusters = rep(0,isig)

  for (k in 1:isig){
    mean_cluster = rep(0,ntype)
    mean_cluster <- colMeans(all_signatures[all_signatures$fit.cluster==k,])[1:ntype,]
    print(mean_cluster)
    variance <- sum(apply(all_signatures[all_signatures$fit.cluster==k,1:ngenome],2,var))
    #for(j in 1:500){
    #  if( all_signature$fit.cluster[j] == k ){
    #    mean_cluster <- mean_cluster + all_signatures[j,1:ntype]
    #    num_entry_in_cluster[k] = num_entry_in_cluster[k] + 1
    #  }
    #}
    mean_cluster = mean_cluster/num_entry_in_cluster[k]
    variance_clusters[k] <- variance
    mean_clusters[,k] <- t(mean_cluster)
  }

  for (k in 1:isig){
    min_signatures_array[,total_columns + k] = signatures_min[, k]
    clusters_array[total_columns + k] = mean_clusters[, k] 
    variance_all[total_columns + k] = variance_clusters[k]
    min_exposures_array[total_columns + k,] = exposures_min[k, ]
    min_error[k] = min
    print(k)
    png(sprintf('bars_num%d_totalnum%d.png',k,isig), width=6000, height=1500, res=300)

    barplot(signatures_min[,k], col=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6), las=2, ylim=c(0, 0.25),names.arg = vocab_mtrx)
    dev.off() 
  }
  total_columns <- total_columns + isig
}

save(min_signatures_array, clusters_array, originalNorm, original,  min_error, variance_all, file="single_double_100run.Rda")

 