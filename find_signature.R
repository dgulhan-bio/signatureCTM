library('tm')
library('topicmodels')
library('ggplot2')
library('fields')
library('nnls')
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

setMethod("calc.frob.error", signature(object = "sign.inst"), function(object, input.matrix, by = "median"){
 exposures <- object@exps
 signatures <- object@signs

 if(by == "mean"){ 
   exposures <- object@mean.exps
   signatures <- object@mean.signs
 }else if(by == "q1"){
   exposures <- object@q1.exps
   signatures <- object@q1.signs
 }else if(by == "q3"){
   exposures <- object@q3.exps
   signatures <- object@q3.signs
 }

 if(round(sum(exposures[,1]), digit = 0) == 1){
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

setMethod("cluster.signs", signature(object = "sign.stack"), function(object){
  nsig <- object@insts[[1]]@nsig
  ntype <- dim(object@insts[[1]]@signs)[[1]]
  ngenome <- dim(object@insts[[1]]@exps)[[2]]
 
  matrix.all.inst.signs <- matrix(0, 0, ntype) 
  matrix.all.inst.exps <- matrix(0, 0, ngenome) 

  perplexity <- rep(0, object@nsig.max - object@nsig.min + 1)
  bic <- rep(0, object@nsig.max - object@nsig.min + 1)  
  print(sprintf('total insts: %d', object@nsig.max - object@nsig.min + 1))

  for(inst in 1:(object@nsig.max - object@nsig.min + 1)){ 
    if(object@insts[[inst]]@nsig != nsig) stop('clustering works on iterations on same number of signs')
    if(object@insts[[inst]]@frob.err > 0.4) next
    
    signatures.this <- object@insts[[inst]]@signs
    exposures.this <- object@insts[[inst]]@exps
    perplexity[[inst]] <- calc.perplexity.inst(object@insts[[inst]],  object@input.dtm, object@method)
    bic[[inst]] <- calc.bic.inst(object@insts[[inst]], object@input.dtm, object@method)@bic
    matrix.all.inst.signs <- rbind(matrix.all.inst.signs, t(signatures.this))
    matrix.all.inst.exps <- rbind(matrix.all.inst.exps, exposures.this)
    rm(signatures.this, exposures.this)
  }    

  sign.clusters <- kmeans(matrix.all.inst.signs, nsig, nstart = 20)

  #calculate silhouette width
  dist.for.all <- rdist(matrix.all.inst.signs)
   
  silhou.ave <- rep(0, nsig)
  silhou.sd <- rep(0, nsig)
  silhou.mean <- rep(0, nsig)
  silhou.median <- rep(0, nsig)
  silhou.q1 <- rep(0, nsig)
  silhou.q3 <- rep(0, nsig)

#  if(exists('indices.filtered')) rm(indices.filtered)
   
 
  for(isig in 1:nsig){
 
    indices.clus <- which(sign.clusters$cluster == isig)

    si.vec <- rep(0, length(indices.clus))
    ai.vec <- rep(0, length(indices.clus))
    bi.vec <- rep(0, length(indices.clus))

    count <- 1
    
    for(iele in indices.clus){
      clus.dist <- dist.for.all[sign.clusters$cluster == isig, iele]
  
      ai <- (length(indices.clus))*mean(clus.dist)/(length(indices.clus) + 1)
      ai.vec[[count]] <- ai

      other.dist <- dist.for.all[-indices.clus, iele]
     
      bi <- min(other.dist)
      bi.vec[[count]] <- bi
 
      si <- 0
      if(ai < bi) si <- (1 - ai/bi)
      else if(ai == bi) si <- 0
      else si <- (bi/ai - 1)
  
      si.vec[[count]] <- si

      count <- count + 1
    } 
    
#    indices.gt.q1 <- (si.vec > median(si.vec))
#    indices.gt.q1 <- (si.vec > median(si.vec[si.vec < median(si.vec)]))
#    indices.gt.q1 <- (si.vec > median(si.vec[si.vec > median(si.vec)]))
    
#    cropped.si <- si.vec[indices.gt.q1]
    
#    print(indices.clus[indices.gt.q1])
#    if(exists('indices.filtered')) indices.filtered <- c(indices.filtered, indices.clus[indices.gt.q1])
#    else indices.filtered <- indices.clus[indices.gt.q1]

#    print(indices.gt.q1)
    
    #print('si.vec')
    #print(si.vec)
    #print('ai.vec')
    #print(ai.vec)
    silhou.ave[[isig]] <- mean(si.vec)   
    silhou.sd[[isig]] <- sd(si.vec) 
    silhou.median[[isig]] <- median(si.vec)
    silhou.q1[[isig]] <- median(si.vec[si.vec < median(si.vec)])
    silhou.q3[[isig]] <- median(si.vec[si.vec > median(si.vec)])
  }

  save(matrix.all.inst.signs, matrix.all.inst.exps, file = sprintf("sign_stack_simul_lda_nsig%d.Rda",nsig))

#  print('indices.filtered')
#  print(indices.filtered)
#  indices.filtered <- sort(indices.filtered)
#  perplexity <- perplexity[unique(round(indices.filtered/nsig))]
#  bic <- bic[unique(round(indices.filtered/nsig))]
  
#  print('indices')
#  print(indices.filtered)
  print('matrix all')
  print(dim(matrix.all.inst.signs))
#  matrix.all.inst.signs <- matrix.all.inst.signs[indices.filtered, ]
  print(dim(matrix.all.inst.signs))
  print('silhou.q1')
  print(silhou.q1)
  print('perplexity')
  print(perplexity)

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

  print('second clustering')
  sign.clusters <- kmeans(matrix.all.inst.signs, nsig)
  median.signs <- matrix(0, ntype, nsig) 
  mean.signs <- matrix(0, ntype, nsig)
  q1.signs <- matrix(0, ntype, nsig)
  q3.signs <- matrix(0, ntype, nsig) 
  median.exps <- matrix(0, nsig, ngenome) 
  mean.exps <- matrix(0, nsig, ngenome)
  q1.exps <- matrix(0, nsig, ngenome)
  q3.exps <- matrix(0, nsig, ngenome) 
  
#  indices.closest.to.median <- rep(0, nsig)
#  indices.closest.to.mean <- rep(0, nsig)
#  indices.closest.to.q1 <- rep(0, nsig)
#  indices.closest.to.q3 <- rep(0, nsig)
  
  for(isig in 1:nsig){
    matrix.inst <- matrix.all.inst.signs[sign.clusters$cluster == isig,]
    
    for(itype in 1:ntype){
      matrix.inst.type <- matrix.inst[, itype]
      median.signs[itype, isig] <- median(matrix.inst.type)
      mean.signs[itype, isig] <- mean(matrix.inst.type)       
      q1.signs[itype, isig] <- median(matrix.inst.type[matrix.inst.type < median.signs[itype, isig]])
      q3.signs[itype, isig] <- median(matrix.inst.type[matrix.inst.type > median.signs[itype, isig]])
    }
    
    dist.to.median <- apply(matrix.all.inst.signs, 1, function(x)sqrt(sum((x - median.signs[,isig])^2)))
    dist.to.mean <- apply(matrix.all.inst.signs, 1, function(x)sqrt(sum((x - mean.signs[,isig])^2)))
    dist.to.q1 <- apply(matrix.all.inst.signs, 1, function(x)sqrt(sum((x - q1.signs[,isig])^2)))
    dist.to.q3 <- apply(matrix.all.inst.signs, 1, function(x)sqrt(sum((x - q3.signs[,isig])^2)))
     
    print(sprintf('min.dist.to.mean %f', min(dist.to.median)))
    
#    indices.closest.to.median[[isig]] <- which(dist.to.median == min(dist.to.median))
#    indices.closest.to.mean[[isig]] <- which(dist.to.mean == min(dist.to.mean))
#    indices.closest.to.q1[[isig]] <- which(dist.to.q1 == min(dist.to.q1))
#    indices.closest.to.q3[[isig]] <- which(dist.to.q3 == min(dist.to.q3))
   
#    median.exps[isig, ] <- matrix.all.inst.exps[indices.closest.to.median[[isig]],]
#    mean.exps[isig, ] <- matrix.all.inst.exps[indices.closest.to.mean[[isig]],]
#    q1.exps[isig, ] <- matrix.all.inst.exps[indices.closest.to.q1[[isig]],]
#    q3.exps[isig, ] <- matrix.all.inst.exps[indices.closest.to.q3[[isig]],]

#    for(igenome in 1:ngenome){
#      matrix.inst.genome <- matrix.all.inst.exps[sign.clusters$cluster == isig, igenome]
#      median.exps[isig, igenome] <- median(matrix.inst.genome)
#      mean.exps[isig, igenome] <- mean(matrix.inst.genome)
#      q1.exps[isig, igenome] <- median(matrix.inst.genome[matrix.inst.genome < median.exps[isig, igenome]])
#      q3.exps[isig, igenome] <- median(matrix.inst.genome[matrix.inst.genome > median.exps[isig, igenome]])
#    }
  }

  for(igenome in 1:ngenome){
    median.exps[, igenome] <- coef(nnls(median.signs, object@input.matrix[,igenome]))
    mean.exps[, igenome] <- coef(nnls(mean.signs, object@input.matrix[,igenome]))
    q1.exps[, igenome] <- coef(nnls(q1.signs, object@input.matrix[,igenome]))
    q3.exps[, igenome] <- coef(nnls(q3.signs, object@input.matrix[,igenome]))
  }



  cluster.centers <- new("sign.inst")
  cluster.centers@signs <- median.signs
  cluster.centers@exps <- median.exps
  cluster.centers@mean.signs <- mean.signs
  cluster.centers@mean.exps <- mean.exps
  cluster.centers@q1.signs <- q1.signs
  cluster.centers@q1.exps <- q1.exps
  cluster.centers@q3.signs <- q3.signs
  cluster.centers@q3.exps <- q3.exps
  cluster.centers@perp.mean <- mean.perplexity
  cluster.centers@perp.median <- median.perplexity
  cluster.centers@perp.sd <- sd.perplexity
  cluster.centers@perp.q1 <- q1.perplexity
  cluster.centers@perp.q3 <- q3.perplexity
  cluster.centers@bic.mean <- mean.bic
  cluster.centers@bic.median <- median.bic
  cluster.centers@bic.sd <- sd.bic
  cluster.centers@bic.q1 <- q1.bic
  cluster.centers@bic.q3 <- q3.bic
  cluster.centers@silhou.width <- silhou.ave
  cluster.centers@silhou.width.sd <- silhou.sd
  cluster.centers@silhou.width.mean <- silhou.median
  cluster.centers@silhou.width.q1 <- silhou.q1
  cluster.centers@silhou.width.q3 <- silhou.q3
  cluster.centers@nsig <- nsig
  
  return(cluster.centers)
})

setMethod("run.calc", signature(object = "sign.stack"), function(object, cluster = FALSE, n.iter = 200){
  index <- 1
  input.matrix <- t(as.matrix(object@input.dtm))
  object@input.matrix <- input.matrix
  print(object@method)

  for(isig in seq(from = object@nsig.min, to = object@nsig.max, by = object@nsig.step)){
    min.err <- 1000
    min.inst <- new("sign.inst")
    print(sprintf("signature %d", isig))
    if(cluster){
      stack.tmp <- new("sign.stack")
      stack.tmp <- set.stack(stack.tmp, input.dtm = object@input.dtm, method = object@method, nsig.min = 1, nsig.max = n.iter)
      stack.tmp@input.matrix <- input.matrix
    }
    for(iter in 1:n.iter){
      inst.tmp <- new("sign.inst")
      inst.tmp@nsig <- isig
      if(object@method == "ctm"){
        inst.method <- CTM(object@input.dtm, isig, method = "VEM", control = list( cg = (list (iter.max = 500L, tol = 10^-4)), seed = iter + 2000))
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
      save(stack.tmp, file = "test.Rda")
      cluster.center <- cluster.signs(stack.tmp)
      object@insts[[index]] <- cluster.center
      object@insts[[index]]@frob.err <- calc.frob.error(cluster.center, input.matrix)
      object@insts[[index]]@frob.err.mean <- calc.frob.error(cluster.center, input.matrix, by = "mean")
      object@insts[[index]]@frob.err.q1 <- calc.frob.error(cluster.center, input.matrix, by = "q1")
      object@insts[[index]]@frob.err.q3 <- calc.frob.error(cluster.center, input.matrix, by = "q3")
    }

    object@min.insts[[index]] <- min.inst
    object@min.insts[[index]]@frob.err <- calc.frob.error(min.inst, input.matrix)
    object@min.insts[[index]]@frob.err.mean <- object@min.insts[[index]]@frob.err
    object@min.insts[[index]]@frob.err.q1 <- object@min.insts[[index]]@frob.err
    object@min.insts[[index]]@frob.err.q3 <- object@min.insts[[index]]@frob.err
    object@min.insts[[index]]@perp.median <- calc.perplexity.inst(min.inst,  object@input.dtm, object@method)
    print(calc.perplexity.inst(min.inst,  object@input.dtm, object@method))
    object@min.insts[[index]]@perp.mean <- object@min.insts[[index]]@perp.median
    object@min.insts[[index]]@perp.q1 <- object@min.insts[[index]]@perp.q1
    object@min.insts[[index]]@perp.q3 <- object@min.insts[[index]]@perp.q3
    object@min.insts[[index]]@bic <- calc.bic.inst(min.inst, object@input.dtm, object@method)
    object@min.insts[[index]]@bic.mean <- object@min.insts[[index]]@bic
    object@min.insts[[index]]@bic.q1 <- object@min.insts[[index]]@bic
    object@min.insts[[index]]@bic.q3 <- object@min.insts[[index]]@bic
    object@min.insts[[index]]@silhou.width <- 1
    object@min.insts[[index]]@silhou.width.sd <- 0
    object@min.insts[[index]]@silhou.width.mean <- 1
    object@min.insts[[index]]@silhou.width.q1 <- 1
    object@min.insts[[index]]@silhou.width.q3 <- 1


    #}else{
    #  object@insts[[index]] <- min.inst
    #}
    index <- index + 1
  }  


  rm(input.matrix)
  rm(min.inst) 
  rm(inst.tmp)
  if(cluster) rm(stack.tmp)
  return(object)
})

setMethod("calc.perplexity.inst", signature(object = "sign.inst"), function(object, new.data, method.desc){
  perplexity <- 0
  if(method.desc == "ctm") perplexity <- perplexity(object@method@ctm, new.data)
  else if(method.desc == "lda") perplexity <- perplexity(object@method@lda, new.data) 
  else stop(sprintf('invalid method %s', method.desc))
  return(perplexity)
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

  return(bic)
})

setMethod("calc.bic", signature(object = "sign.stack"), function(object){
  for(index in 1:object@num.insts){
    object@insts[[index]] <- calc.bic.inst(object@insts[[index]], object@input.dtm, object@method)
  }
  return(object)
})

find_signatures <- function(){

  nsig_max = 9
  nsig_min = 2
  
#  input_file_list <- 'test.txt'
#  dtm <- make_dtm_from_vcf_list(input_file_list)

  load('kidney_dtm_ncontext_3.Rda')

#  dtm <- make_dtm_from_vcf_list(input_file_list, ncontext = 5)
  
  signature_stack <- new("sign.stack") 
  str(dtm)

  signature_stack <- set.stack(signature_stack, input.dtm = dtm, nsig.min = nsig_min, nsig.max = nsig_max, method = "lda") 
  signature_stack <- run.calc(signature_stack, cluster = TRUE) 
  # print('calculatingperplexity') 
  # signature_stack <- calc.perplexity(signature_stack,dtm) 
  # print('calculating bic') 
  # signature_stack <- calc.bic(signature_stack)
 
  save(signature_stack, file = "test_output_kidney_dim_2_9_lda_sil_min_context3.Rda") 
  print(signature_stack@num.insts)
  for(index in 1:signature_stack@num.insts){
    inst <- signature_stack@insts[[index]]
    print(sprintf('nsig: %d, error: %.3f, perplexity: %.3f, bic %.0f', inst@nsig, inst@frob.err, inst@mean.perp, inst@bic))
  } 
  rm(signature_stack)
}
