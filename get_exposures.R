sort_signatures <- function(signatures, measure){
  sig_nrow = nrow(signatures)
  sig_ncolumn = ncol(signatures)
  
  enlarged = matrix(0, sig_ncolumn, sig_nrow + 1) 
  
  enlarged[,1:sig_nrow] = t(signatures)
  enlarged[,sig_nrow+1] = measure
   
  signature_frame <- data.frame(enlarged)
  colnames(signature_frame)[[sig_nrow+1]] = 'measure'
  sorted_matrix <- enlarged[order(signature_frame$measure, decreasing = TRUE),]
  signatures_sorted <- sorted_matrix[,1:sig_nrow]
  return(t(signatures_sorted))  
}

get_exposures <- function(){

  source('match_signatures.R')
  source('get_error_frob.R')
  source('calculate_bic.R')
  library('nnls')
  load('100run.Rda')
  
  ngenome = dim(original)[[2]]
  nsig = dim(signatures_min)[[2]]
  ntype = dim(signatures_min)[[1]]
  
  #convert fractional exposures to absolute exposures
  for(i in 1:ngenome){
    sum = sum(original[,i])
    exposures_min[,i] = sum*exposures_min[,i]  
  }
  
  exposures_nnls <- matrix(0, nsig, ngenome)   
  
  #for each genome find optimal number of signatures and the corresponding exposures
  for(igenome in 1:ngenome){
    this_genome <- original[,igenome]
    similarity <- rep(0, nsig)
  
    for(i in 1:nsig){
      similarity[[i]] <- this_genome %*% signatures_min[,i] / sqrt(this_genome %*% this_genome * signatures_min[,i] %*% signatures_min[,i])
    }

    #sort signatures such that the most similar to this genome comes first
    signatures_sorted <- sort_signatures(signatures_min, similarity)
 
    min_exposure <- rep(0,1)
    min_bic = 10000
    min_error = 10000
    min_signature <- matrix(0,2,1)
    
    #check all possible number of total signatures
    for(i in 1:nsig){
      #get subset of signatures acc. to the total number 
      signature_subset <- matrix(0, ntype, i) 
      for(j in 1:i){
       signature_subset[,j] <- signatures_sorted[,j]
      }
      
      #calculate exposure
      fit_nnls<-nnls(signature_subset, this_genome)
      this_exposure <- coef(fit_nnls)
      remnant <- signature_subset %*% this_exposure - this_genome
      #print(dim(remnant))
      frob_error <- sqrt(sum(diag(t(remnant) %*% remnant)))/sqrt(sum(diag(this_genome %*% this_genome)))
      
      #calculate bic, the min bic is chosen to be the optimal number of signatures
      bic <- calculate_bic(signature_subset, this_exposure, this_genome)
      if(bic < min_bic && min_error - frob_error > 0.01){
        min_bic = bic
        min_error = frob_error
        min_exposure <- this_exposure
        min_signature <- signature_subset      
      }
    }
   print(sprintf('genome %d, num sig = %d, frob error = %.3f, bic = %.3f', igenome, length(min_exposure), frob_error, bic))
    
    #fill in info of best exposures for the plot  
    for(i in 1:nsig){
      for(j in 1:length(min_exposure)){
        comparison <- (min_signature[,j] %*% signatures_min[,i]/sqrt(min_signature[,j] %*% min_signature[,j] * signatures_min[,i] %*% signatures_min[,i])) 
        if( comparison == 1 ){
          #print(dim(exposures_nnls))
          #print(j)
          #print(igenome)
          exposures_nnls[i, igenome] = min_exposure[[j]]
        }
      }
    }
    #return(this_exposure)
  }
  print(sprintf('frob error all = %.3f', get_error_frob(original, signatures_min, exposures_nnls)))
  png('exposuresComparison.png', width=6000, height=1500, res=300)
  par(mfrow=c(1,3))
  barplot(exposures_nnls, col = c(3,4,5,6,7))
  mtext(text='decompose in ctm',side=3)
  barplot(exposures_min, col = c(3,4,5,6,7))
  dev.off()
}
