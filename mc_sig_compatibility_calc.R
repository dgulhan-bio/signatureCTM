library('ggplot2')
source('match_signatures.R')

mc_sig_compatibility_calc <- function(){
  load('out_genome_sig_asso_genome_by_genome.Rda')
  ntype = dim(genome_sig_assoc)[[1]]
  ngenome = dim(genome_sig_assoc)[[2]]
  nsig = dim(genome_sig_assoc)[[3]]
    
  total_snvs = rep(0, ngenome)
  
  for(i in 1:ngenome){
    total_snvs[[i]] <- sum(genome_sig_assoc[,i,])
  }
  
  load('exposures_genome_by_genome.Rda')

  compatibility <- matrix(0, ngenome, nsig)
  fraction <- matrix(0, ngenome, nsig)
  for(i in 1:ngenome){
    normalized_mc <- genome_sig_assoc[,i,]/total_snvs[[i]]
    print(dim(normalized_mc))
    index <- match_signatures(signatures_min, normalized_mc)
    for(j in 1:nsig){
      if( index[[j]] == -1 ){
        compatibility[i,j] = 0
        fraction[i,j] = 0
        next
      }
      compatibility[i,j] <- normalized_mc[,index[[j]]] %*% signatures_min[,j] / sqrt(normalized_mc[,index[[j]]] %*% normalized_mc[,index[[j]]] * signatures_min[,j] %*% signatures_min[,j] )
      fraction[i,j] <- sum(normalized_mc[,j])
    }
  }

  save(genome_sig_assoc, fraction, compatibility, file = 'out_mc_frac_comp.Rda')
  
}