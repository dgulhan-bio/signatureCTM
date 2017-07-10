library('tm')
library('topicmodels')
library('ggplot2')
source('get_error_frob.R')
source('class_defs.R')
source('make_dtm_from_vcf.R')
  
setMethod("set.mine", signature(object = "sign.mine"), function(object, input.dtm, nsig.min, nsig.max, nsig.step = as.integer(1), method = 'ctm'){
  object@nsig.min <- as.integer(nsig.min)
  object@nsig.max <- as.integer(nsig.max)
  object@nsig.step <- as.integer(nsig.step)

  object@method <- method
  object@input.dtm <- input.dtm
  object@types <- input.dtm$dimnames$Terms
  number.of.miners <- (nsig.max - nsig.min)/(nsig.step) + 1
  object@num.miners <- as.integer(number.of.miners)
  object@miners <- lapply( rep("sign.miner", number.of.miners), FUN = new )
  index <- 1

  for(nsig_miner in seq(from = nsig.min, to = nsig.max, by = nsig.step)){
    object@miners[[index]]@nsig <- as.integer(nsig_miner)
    index <- index + 1
  }
  return(object)
})

setMethod("calc.frob.error", signature(object = "sign.miner"), function(object, input.matrix){
 exposures <- object@exps
 signatures <- object@signs
 if(sum(exposures[,1]) == 1){
    for(i in 1:ncol(exposures)){
      exposures[,i] = sum(input.matrix[,i])*exposures[,i]
    }
  }
 
  reco <- signatures %*% exposures

  diff <- input.matrix - reco

  error <- sqrt(sum(diag(t(diff) %*% (diff))))
  error <- error/sqrt(sum(diag(t(input.matrix) %*% (input.matrix))))
  return(error)
})

setMethod("run.mine", signature(object = "sign.mine"), function(object){
  index <- 1
  input.matrix <- t(as.matrix(object@input.dtm))
  for(isig in seq(from = object@nsig.min, to = object@nsig.max, by = object@nsig.step)){
    min.err <- 1000
    min.miner <- new('sign.mine')
    print(sprintf('signature %d', isig))
    for(iter in 1:100){
      if(object@method == "ctm"){
        miner.method <- CTM(object@input.dtm, isig, method = "VEM", control = list( cg = (list (iter.max = 1000L, tol = 10^-4))))
        object@miners[[index]]@method <- new('topic.mod.obj', ctm = miner.method)
      }else if(object@method == "lda"){
        miner.method <- LDA(object@input.dtm, isig, method = "Gibbs")
        object@miners[[index]]@method <- new('topic.mod.obj', lda = miner.method)
      }else stop(sprintf('invalid method %s', object@method))

      beta<-slot(miner.method, 'beta')
      gamma<-slot(miner.method, 'gamma')
      object@miners[[index]]@signs <- exp(t(beta))
      object@miners[[index]]@exps <- t(gamma)
      error <- calc.frob.error(object@miners[[index]], input.matrix)
      object@miners[[index]]@frob.err <- error

      if(error < min.err){
        min.err <- error
        min.miner <- object@miners[[index]]
      }
    }
    object@miners[[index]] <- min.miner
    index <- index + 1
  }
  return(object)
})

setMethod("calc.perplexity.miner", signature(object = "sign.miner"), function(object, new.data, method.desc){
  if(method.desc == "ctm") object@perplexity <- perplexity(object@method@ctm, new.data)
  else if(method.desc == "lda") object@perplexity <- perplexity(object@method@lda, new.data) 
  else stop(sprintf('invalid method %s', method.desc))
  return(object)
})

setMethod("calc.perplexity", signature(object = "sign.mine"), function(object, new.data){
  for(index in 1:object@num.miners){
    object@miners[[index]] <- calc.perplexity.miner(object@miners[[index]], new.data, object@method)
  }
  return(object)
})


setMethod("calc.bic.miner", signature(object = "sign.miner"), function(object, input.dtm, method.desc){
  if(method.desc == "ctm") method <- object@method@ctm
  else if(method.desc == "lda") method <- object@method@lda
  else stop(sprintf('invalid method %s', method.desc))
   
  ngenome <- input.dtm$nrow
  ntype <- input.dtm$ncol
  nsig <- object@nsig

  bic <- (2*as.numeric(logLik(method))) - (ngenome + ntype)*log(ngenome)*nsig
  object@bic <- bic
  return(object)
})

setMethod("calc.bic", signature(object = "sign.mine"), function(object){
  for(index in 1:object@num.miners){
    object@miners[[index]] <- calc.bic.miner(object@miners[[index]], object@input.dtm, object@method)
  }
  return(object)
})

find_signatures <- function(){

  nsig_max = 6
  nsig_min = 2
  
  input_file_list <- 'test.txt'
  dtm <- make_dtm_from_vcf_list(input_file_list)
  
  signature_mine <- new("sign.mine")
  signature_mine <- set.mine(signature_mine, input.dtm = dtm, nsig.min = nsig_min, nsig.max = nsig_max)
  signature_mine <- run.mine(signature_mine)
  print('calculating perplexity')
  signature_mine <- calc.perplexity(signature_mine, dtm)
  print('calculating bic')
  signature_mine <- calc.bic(signature_mine)
 
  print(signature_mine@num.miners)
  for(index in 1:signature_mine@num.miners){
    miner <- signature_mine@miners[[index]]
    print(sprintf('nsig: %d, error: %.3f, perplexity: %.3f, bic %.0f', miner@nsig, miner@frob.err, miner@perplexity, miner@bic))
  } 

}
