library('tm')
library('topicmodels')
library('ggplot2')
library('fields')
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

setMethod("cluster.signs", signature(object = "sign.stack"), function(object, by = "median"){
  nsig <- object@insts[[1]]@nsig
  ntype <- dim(object@insts[[1]]@signs)[[1]]
  ngenome <- dim(object@insts[[1]]@exps)[[2]]
 
  matrix.all.inst.signs <- rep(0, 0, ntype) 
  matrix.all.inst.exps <- rep(0, 0, ngenome) 

  perplexity <- rep(0, object@nsig.max - object@nsig.min + 1)
  bic <- rep(0, object@nsig.max - object@nsig.min + 1)  

  for(inst in 1:(object@nsig.max - object@nsig.min + 1)){
 
    if(object@insts[[inst]]@nsig != nsig) stop('clustering works on iterations on same number of signs')
 
    signatures.this <- object@insts[[inst]]@signs
    exposures.this <- object@insts[[inst]]@exps
    perplexity[[inst]] <- calc.perplexity.inst(object@insts[[inst]],  object@input.dtm, object@method)@perplexity
    bic[[inst]] <- calc.bic.inst(object@insts[[inst]], object@input.dtm, object@method)@bic

    matrix.all.inst.signs <- rbind(matrix.all.inst.signs, t(signatures.this))
    matrix.all.inst.exps <- rbind(matrix.all.inst.exps, exposures.this)
    rm(signatures.this, exposures.this)
  }    

  mean.perplexity <-  mean(perplexity)
  median.perplexity <- median(perplexity)
  sd.perplexity <- sd(perplexity)
  q1.perplexity <- median(perplexity[perplexity < median.perplexity])
  q3.perplexity <- median(perplexity[perplexity > median.perplexity])

  mean.bic <- mean(bic)
  median.bic <- median(bic)
  sd.bic <- sd(bic)

  q1.bic <- median(bic[bic < median.bic])
  q3.bic <- median(bic[bic > median.bic])

  sign.clusters <- kmeans(matrix.all.inst.signs, nsig)
  median.signs <- matrix(0, ntype, nsig) 
  median.exps <- matrix(0, nsig, ngenome) 
  mean.signs <- matrix(0, ntype, nsig)
  mean.exps <- matrix(0, nsig, ngenome)

  for(isig in 1:nsig){
    for(itype in 1:ntype){
      median.signs[itype, isig] <- median(matrix.all.inst.signs[sign.clusters$cluster == isig, itype])
      mean.signs[itype, isig] <- mean(matrix.all.inst.signs[sign.clusters$cluster == isig, itype]) 
    }
    for(igenome in 1:ngenome){
      median.exps[isig, igenome] <- median(matrix.all.inst.exps[sign.clusters$cluster == isig, igenome])
      mean.exps[isig, igenome] <- mean(matrix.all.inst.exps[sign.clusters$cluster == isig, igenome])
    }
  }
  
  #calculate silhouette width
  dist.for.all <- rdist(matrix.all.inst.signs)
   
  silhou.ave <- rep(0, nsig)
  silhou.sd <- rep(0, nsig)
  for(isig in 1:nsig){
    print(sprintf('isig = %d', isig))

    si.vec <- rep(0, dim(dist.for.all)[[1]])
    ai.vec <- rep(0, dim(dist.for.all)[[1]])
    bi.vec <- rep(0, dim(dist.for.all)[[1]])
 
    indices.clus <- which(sign.clusters$cluster == isig)
    print(indices.clus)
    
    for(iele in indices.clus){
      clus.dist <- dist.for.all[sign.clusters$cluster == isig, iele]
  
      ai <- (length(indices.clus))*mean(clus.dist)/(length(indices.clus) + 1)
      ai.vec[[iele]] <- ai

      other.dist <- dist.for.all[-indices.clus, iele]
     
      bi <- min(other.dist)
      bi.vec[[iele]] <- bi
 
      si <- 0
      if(ai < bi) si <- (1 - ai/bi)
      else if(ai == bi) si <- 0
      else si <- (bi/ai - 1)
  
      si.vec[[iele]] <- si
    } 
    print(si.vec)
    print(si.vec)
    silhou.ave[[isig]] <- mean(si.vec)   
    silhou.sd[[isig]] <- sd(si.vec) 
  }

  cluster.centers <- new("sign.inst")
  cluster.centers@signs <- (median.signs)
  cluster.centers@exps <- (median.exps)
  cluster.centers@perp.mean <- mean.perplexity
  cluster.centers@perp.median <- median.perplexity
  cluster.centers@perp.sd <- sd.perplexity
  cluster.centers@perp.q1 <- q1.perplexity
  cluster.centers@perp.q3 <- q3.perplexity
  cluster.centers@bic.mean <- median.bic
  cluster.centers@bic.median <- median.bic
  cluster.centers@bic.sd <- sd.bic
  cluster.centers@bic.q1 <- q1.bic
  cluster.centers@bic.q3 <- q3.bic
  cluster.centers@silhou.width <- silhou.ave
  cluster.centers@silhou.width.sd <- silhou.sd
  cluster.centers@nsig <- nsig
  
  return(cluster.centers)
})

setMethod("run.calc", signature(object = "sign.stack"), function(object, cluster = FALSE, n.iter = 100, by = "median"){
  index <- 1
  input.matrix <- t(as.matrix(object@input.dtm))
  for(isig in seq(from = object@nsig.min, to = object@nsig.max, by = object@nsig.step)){
    min.err <- 1000
    min.inst <- new("sign.inst")
    print(sprintf("signature %d", isig))
    if(cluster){
      stack.tmp <- new("sign.stack")
      stack.tmp <- set.stack(stack.tmp, input.dtm = dtm, nsig.min = 1, nsig.max = n.iter)
    }
    for(iter in 1:n.iter){
      inst.tmp <- new("sign.inst")
      inst.tmp@nsig <- isig
      if(object@method == "ctm"){
        inst.method <- CTM(object@input.dtm, isig, method = "VEM", control = list( cg = (list (iter.max = 1000L, tol = 10^-4)), seed = iter + 2000))
        inst.tmp@method <- new('topic.mod.obj', ctm = inst.method)
      }else if(object@method == "lda"){
        inst.method <- LDA(object@input.dtm, isig, method = "Gibbs", control = list(seed = iter + 3000))
        inst.tmp@method <- new('topic.mod.obj', lda = inst.method)
      }else stop(sprintf('invalid method %s', object@method))

      beta<-slot(inst.method, 'beta')
      gamma<-slot(inst.method, 'gamma')
      inst.tmp@signs <- exp(t(beta))
      inst.tmp@exps <- t(gamma)
      error <- calc.frob.error(inst.tmp, input.matrix)
      inst.tmp@frob.err <- error

      object@insts[[index]]@frob.err <- error
      if(cluster){
        stack.tmp@insts[[iter]] <- inst.tmp
      }
      if(error < min.err){
        min.err <- error
        min.inst <- inst.tmp
      }
    }

    if(cluster){
      cluster.center <- cluster.signs(stack.tmp, by)
      object@insts[[index]] <- cluster.center
      print(dim(input.matrix))
      object@insts[[index]]@frob.err <- calc.frob.error(cluster.center, input.matrix)
      rm(stack.tmp)
    }else{
      object@insts[[index]] <- min.inst
    }
    index <- index + 1
  }  
  rm(input.matrix)
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

  nsig_max = 3
  nsig_min = 2
  
  input_file_list <- 'test_small.txt'
  dtm <- make_dtm_from_vcf_list(input_file_list)
  
  signature_stack <- new("sign.stack")
  signature_stack <- set.stack(signature_stack, input.dtm = dtm, nsig.min = nsig_min, nsig.max = nsig_max)
  signature_stack <- run.calc(signature_stack, cluster = TRUE)
#  print('calculating perplexity')
#  signature_stack <- calc.perplexity(signature_stack, dtm)
#  print('calculating bic')
#  signature_stack <- calc.bic(signature_stack)
 
  save(signature_stack, file = "test_output_dim_2_6_sil.Rda") 
  print(signature_stack@num.insts)
  for(index in 1:signature_stack@num.insts){
    inst <- signature_stack@insts[[index]]
    print(sprintf('nsig: %d, error: %.3f, perplexity: %.3f, bic %.0f', inst@nsig, inst@frob.err, inst@perplexity, inst@bic))
  } 

}
