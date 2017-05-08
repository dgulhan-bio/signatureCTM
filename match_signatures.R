match_signatures<-function(signatures, signatures_min){
  index = rep(0,5)
  largest_similarity = rep(0,5)
  for(i in 1:nsig){
    max = 0
    for(j in 1:nsig){
      similarity = signatures_min[,i]%*%signatures[,j]/sqrt(signatures_min[,i]%*%signatures_min[,i]*signatures[,j]%*%signatures[,j])
      print(sprintf('%d %d %.2f',i, j, similarity))
      if(similarity > max){
        max = similarity
        index[[i]] = j
      }
    }
    largest_similarity[[i]] = max
  }
  return(index)
}