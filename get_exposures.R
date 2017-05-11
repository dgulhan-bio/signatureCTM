savePlot <- function(myPlot,i) {
        pdf(sprintf("myPlot%d.pdf",i))
        print(myPlot)
        dev.off()
}

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

get_exposures <- function(add_nmf_result = TRUE){

  source('match_signatures.R')
  source('get_error_frob.R')
  source('calculate_bic.R')
  library('nnls')
  library('ggplot2')
  load('100run.Rda')
  
  barcols = rep(0,192)
  barcols[1:16] = 1
  barcols[17:32] = 2
  barcols[33:48] = 3
  barcols[49:64] = 4
  barcols[65:80] = 5
  barcols[81:96] = 6
  barcols[97:112] = 1
  barcols[113:128] = 2
  barcols[129:144] = 3
  barcols[145:160] = 4
  barcols[161:176] = 5
  barcols[177:192] = 6
  print(barcols)
  print(length(barcols))
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
    min_remnant <- rep(0,96)
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
      reco <- signature_subset %*% this_exposure
      reco[ reco < 1 ] = 0
      remnant <- reco - this_genome
      
      #print(dim(remnant))
      frob_error <- sqrt(sum(diag(t(remnant) %*% remnant)))/sqrt(sum(diag(this_genome %*% this_genome)))
      
      #calculate bic, the min bic is chosen to be the optimal number of signatures
      bic <- calculate_bic(signature_subset, this_exposure, this_genome)
      if(bic < min_bic && min_error - frob_error > 0.01){
        min_bic = bic
        min_error = frob_error
        min_exposure <- this_exposure
        min_signature <- signature_subset      
        min_remnant <- remnant 
      }
    }

    combine_remnant_original = matrix(0, 192, 3)
    combine_remnant_original[1:96,1] = t(remnant)
    combine_remnant_original[97:192,1] = t(this_genome)
    combine_remnant_original[1:96,2] = c(1:96)
    combine_remnant_original[97:192,2] = c(1:96)
    combine_remnant_original[1:96,3] = rep(3,96)
    combine_remnant_original[97:192,3] = rep(4,96)

    df_rem_ori <- data.frame(combine_remnant_original)
    colnames(df_rem_ori) <- c('remnant', 'number', 'before')
    #print(df_rem_ori)
    myplot<-ggplot(df_rem_ori, aes(x = number, y = remnant, fill = before)) + geom_bar(stat = "identity")
    savePlot(myplot,igenome)

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
  plot_width = 4000
  if(add_nmf_result) plot_width = 6000
  png('exposuresComparison.png', width = plot_width, height = 1500, res = 300)
  par(mfrow = c(1,3))
  barplot(exposures_nnls, col = c(3,4,5,6,7))
  mtext(text = 'nnls to decompose ctm sig',side = 3)
  barplot(exposures_min, col = c(3,4,5,6,7))
  mtext(text = 'result ctm',side = 3)
  if(add_nmf_result){
    source('match_signatures.R')
    source('get_error_frob.R')

    signatures_nmf<-as.matrix(read.csv('nmfprocesses.dat',header = FALSE))
    exposures_nmf<-as.matrix(read.csv('nmfexposures.dat',header = FALSE))
    genomes_nmf<-as.matrix(read.csv('originalgenomes.dat',header = FALSE))
    match_index<-matchSignatures(signatures_nmf, signatures_min)

    exposures_nmf_resorted<-exposures_nmf
    for(i in 1:5){
      exposures_nmf_resorted[i,] = exposures_nmf[match_index[[i]],]
    }
    barplot(exposures_nmf_resorted, col = c(3,4,5,6,7))
    mtext(text='result nmf sanger',side=3)
  }

  legend("topright", legend=c("Signature1", "Signature2", "Signature3", "Signature4", "Signature5"), fill = c(3,4,7,6,5), cex = 1.5)

  dev.off()

  save(exposures_nnls, original, signatures_min, file = "exposures_genome_by_genome.Rda")
}
