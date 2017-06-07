match_signatures<-function(signatures, signatures_new){
  index <- rep(0,5)
  largest_similarity <- rep(0,5)
  nsig <- dim(signatures)[[2]]
  
  for(i in 1:nsig){
    max <- 0
    if(sum(signatures_new[,i]) == 0){
      index[[i]] <- (-1)
      next
    }
    for(j in 1:nsig){
      similarity <- signatures_new[,i]%*%signatures[,j]/sqrt(signatures_new[,i]%*%signatures_new[,i]*signatures[,j]%*%signatures[,j])
      print(sprintf('%d %d %.2f %.6f %.6f',i, j, similarity, (signatures_new[,i]%*%signatures_new[,i]), (signatures[,j]%*%signatures[,j])))
      if(similarity > max){
        max <- similarity
        index[[i]] <- j
      }
    }
    largest_similarity[[i]] <- max
  }
  return(index)
}