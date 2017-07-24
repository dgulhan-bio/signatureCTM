library('tm')
library('topicmodels')
library('ggplot2')
#source('get_error_frob.R')
source('class_defs.R')
source('make_dtm_from_vcf.R')
  
setMethod("set.stack", signature(object = "sign.stack"), function(object, input.dtm, nsig.min, nsig.max, nsig.step = as.integer(1), method = 'ctm'){
  object@nsig.min <- as.integer(nsig.min)
  object@nsig.max <- as.integer(nsig.max)
  object@nsig.step <- as.integer(nsig.step)

  object@method <- method
  object@input.dtm <- input.dtm
  object@types <- input.dtm$dimnames$Terms
  number.of.insts <- (nsig.max - nsig.min)/(nsig.step) + 1
  object@num.insts <- as.integer(number.of.insts)
  object@insts <- lapply( rep("sign.inst", number.of.insts), FUN = new )
  index <- 1

  for(nsig_inst in seq(from = nsig.min, to = nsig.max, by = nsig.step)){
    object@insts[[index]]@nsig <- as.integer(nsig_inst)
    index <- index + 1
  }
  return(object)
})

setMethod("calc.frob.error", signature(object = "sign.inst"), function(object, input.matrix){
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

setMethod("run.calc", signature(object = "sign.stack"), function(object){
  index <- 1
  input.matrix <- t(as.matrix(object@input.dtm))
  for(isig in seq(from = object@nsig.min, to = object@nsig.max, by = object@nsig.step)){
    min.err <- 1000
    min.inst <- new('sign.stack')
    print(sprintf('signature %d', isig))
    for(iter in 1:100){
      if(object@method == "ctm"){
        inst.method <- CTM(object@input.dtm, isig, method = "VEM", control = list( cg = (list (iter.max = 1000L, tol = 10^-4))))
        object@insts[[index]]@method <- new('topic.mod.obj', ctm = inst.method)
      }else if(object@method == "lda"){
        inst.method <- LDA(object@input.dtm, isig, method = "Gibbs")
        object@insts[[index]]@method <- new('topic.mod.obj', lda = inst.method)
      }else stop(sprintf('invalid method %s', object@method))

      beta<-slot(inst.method, 'beta')
      gamma<-slot(inst.method, 'gamma')
      object@insts[[index]]@signs <- exp(t(beta))
      object@insts[[index]]@exps <- t(gamma)
      error <- calc.frob.error(object@insts[[index]], input.matrix)
      object@insts[[index]]@frob.err <- error

      if(error < min.err){
        min.err <- error
        min.inst <- object@insts[[index]]
      }
    }
    object@insts[[index]] <- min.inst
    index <- index + 1
  }
  return(object)
})

setMethod("calc.perplexity.inst", signature(object = "sign.inst"), function(object, new.data, method.desc){
  if(method.desc == "ctm") object@perplexity <- perplexity(object@method@ctm, new.data)
  else if(method.desc == "lda") object@perplexity <- perplexity(object@method@lda, new.data) 
  else stop(sprintf('invalid method %s', method.desc))
  return(object)
})

setMethod("calc.perplexity", signature(object = "sign.stack"), function(object, new.data){
  for(index in 1:object@num.insts){
    object@insts[[index]] <- calc.perplexity.inst(object@insts[[index]], new.data, object@method)
  }
  return(object)
})


setMethod("calc.bic.inst", signature(object = "sign.inst"), function(object, input.dtm, method.desc){
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

setMethod("calc.bic", signature(object = "sign.stack"), function(object){
  for(index in 1:object@num.insts){
    object@insts[[index]] <- calc.bic.inst(object@insts[[index]], object@input.dtm, object@method)
  }
  return(object)
})

find_signatures <- function(){

  nsig_max = 6
  nsig_min = 2
  
  input_file_list <- 'test.txt'
  dtm <- make_dtm_from_vcf_list(input_file_list)
  
  signature_stack <- new("sign.stack")
  signature_stack <- set.stack(signature_stack, input.dtm = dtm, nsig.min = nsig_min, nsig.max = nsig_max)
  signature_stack <- run.calc(signature_stack)
  print('calculating perplexity')
  signature_stack <- calc.perplexity(signature_stack, dtm)
  print('calculating bic')
  signature_stack <- calc.bic(signature_stack)
 
  print(signature_stack@num.insts)
  for(index in 1:signature_stack@num.insts){
    inst <- signature_stack@insts[[index]]
    print(sprintf('nsig: %d, error: %.3f, perplexity: %.3f, bic %.0f', inst@nsig, inst@frob.err, inst@perplexity, inst@bic))
  } 

}
