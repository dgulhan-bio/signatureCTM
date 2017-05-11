MC_decomposition <- function(genomes, exposures, signatures){
  ntype <- dim(genomes)[[1]]
  ngenome <- dim(genomes)[[2]]
  nsig <- dim(signatures)[[2]]
  print(sprintf('ntype %d ngenome %d nsig %d',ntype, ngenome, nsig))
  reco <- signatures %*% exposures
 
  MAX <- max(max(genomes), max(reco))
  genome_sig_assoc = array(0, c(96, 127, 5))
 
  for(i in 1:ntype){
    for(j in 1:ngenome){
      this_snv_num <- genomes[i,j]
      if(this_snv_num == 0) next
      
      # print(sprintf('type = %d genome = %d snv_num = %d', i, j, this_snv_num))
      total = 0
      for(isnv in 1:this_snv_num){
        #print(total)
        trial = 0
        while(TRUE){
          dice_sig <- floor(nsig*runif(1, 0, 1)) + 1
          #print(sprintf('dice_sig %d', dice_sig))
          number_prob <- signatures[i, dice_sig] * exposures[dice_sig, j]
          dice_num <- runif(1, 0, MAX)
          if(number_prob >= dice_num){
            # print(sprintf('dice %.3f  reco %.3f', dice_num, number_prob))
            genome_sig_assoc[i, j, dice_sig] = genome_sig_assoc[i, j, dice_sig] + 1
            total = total + 1
            break
          }
          trial = trial + 1
        }
      }
      print(sprintf('type = %d genome = %d snv_num = %d total = %d', i, j, this_snv_num, total))
    }
  }
  save(genome_sig_assoc, file = 'out_genome_sig_asso_genome_by_genome.Rda')
  test <- rep(0,96)
  for(i in 1:nsig){
    test = test + genome_sig_assoc[, 1, i]
  }
  print(test - genomes[,1])
}

#run_decomposition<-function(){
  load('exposures_genome_by_genome.Rda')
  MC_decomposition(original, exposures_nnls, signatures_min)
#}