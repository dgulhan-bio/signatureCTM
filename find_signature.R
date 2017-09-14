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
  print('end of set stack')
  return(object)
})

#getters used in apply
setMethod(".same.dim", signature(object = "sign.inst"), function(object, nsig.comp){
  return(object@nsig != nsig.comp)
})

setMethod("get.signs", signature(object = "sign.inst"), function(object, by = "inst"){
  return(object@signs[by])
})

setMethod("get.exps", signature(object = "sign.inst"), function(object, by = "inst"){
  return(object@exps[by])
})

setMethod("get.frob.err", signature(object = "sign.inst"), function(object, by = "inst"){
  return(object@frob.err[by])
})

#looks at the list of signatures, median, mean, q1, q3 for clusters
#or just inst for single iterations and returns the error as a list
#with the same names
setMethod("calc.frob.error", signature(object = "sign.inst"), function(object, input.matrix){
  by <- names(object@signs)
  error.list <- list()
  
  for(icalc in length(by)){
    exposures <- object@exps[[icalc]]
    signatures <- object@signs[[icalc]]

    if(round(sum(exposures[,1]), digit = 0) == 1){
      for(i in 1:ncol(exposures)){
        exposures[,i] = sum(input.matrix[,i])*exposures[,i]
      }
    }
 
    reco <- signatures %*% exposures
    diff <- input.matrix - reco
    error <- sqrt(sum(diag(t(diff) %*% (diff))))
    error <- error/sqrt(sum(diag(t(input.matrix) %*% (input.matrix))))

    error.list[[by[[icalc]]]] <- error
  }
  print(error.list)
  return(error.list)
})


setMethod("cluster.signs", signature(object = "sign.stack"), function(object){
  nsig <- object@insts[[1]]@nsig
  ntype <- dim(object@insts[[1]]@signs$inst)[[1]]
  ngenome <- dim(object@insts[[1]]@exps$inst)[[2]]
  ninst <- object@nsig.max - object@nsig.min + 1
  
  perplexity <- rep(0, ninst)
  bic <- rep(0, ninst)
  print(sprintf('total insts: %d', ninst))

  if(sum(sapply(object@insts, .same.dim, nsig.comp = nsig)) != 0) stop('clustering works on iterations on same number of sign\
s')
  perplexity <- sapply(object@insts, calc.perplexity.inst, new.data = object@input.dtm, method.desc = object@method)
  bic <- sapply(object@insts, calc.bic.inst, input.dtm = object@input.dtm, method.desc = object@method)

  list.all.inst.signs <- sapply(object@insts, get.signs)
  list.all.inst.exps <- sapply(object@insts, get.exps)
 
  matrix.all.inst.signs <- matrix(unlist(list.all.inst.signs), byrow = TRUE, nrow = nsig*ninst)
  matrix.all.inst.exps <- matrix(unlist(list.all.inst.exps), byrow = TRUE, nrow = nsig*ninst)
  save(matrix.all.inst.signs, file = 'test.Rda')
  rm(list.all.inst.signs)
  rm(list.all.inst.exps)

  error.vector <- sapply(object@insts, get.frob.err, by = "inst")
  error.vector <- as.vector(sapply(error.vector, rep, times = nsig))
#  matrix.all.inst.signs <- matrix.all.inst.signs[!which(error.vector > 0.5),] 
#  matrix.all.inst.exps <- matrix.all.inst.exps[!which(error.vector > 0.5),]
  
  sign.clusters <- kmeans(matrix.all.inst.signs, centers = nsig, nstart = 20)
  #calculate silhouette width
  dist.for.all <- rdist(matrix.all.inst.signs)
  sd.silhou <- rep(0, nsig)
  mean.silhou <- rep(0, nsig)
  median.silhou <- rep(0, nsig)
  q1.silhou <- rep(0, nsig)
  q3.silhou <- rep(0, nsig)
  for(isig in 1:nsig){
    indices.clus <- which(sign.clusters$cluster == isig)
    dist.same.clus <- dist.for.all[indices.clus,]
    dist.other.clus <- dist.for.all[-indices.clus,]

    ai.vec <- (length(indices.clus)*apply(dist.same.clus, 2, mean)/(length(indices.clus) + 1))[indices.clus]
    bi.vec <- apply(dist.other.clus, 2, min)[indices.clus]
    f.sil <- function(x,y){
      if(x < y) return( 1- x/y)
      else if(x == y) return(0)
      else return(y/x - 1)
    }
    si.vec <- mapply(f.sil, ai.vec, bi.vec)

    mean.silhou[[isig]] <- mean(si.vec)
    sd.silhou[[isig]] <- sd(si.vec)
    median.silhou[[isig]] <- median(si.vec)
    q1.silhou[[isig]] <- median(si.vec[si.vec < median(si.vec)])
    q3.silhou[[isig]] <- median(si.vec[si.vec > median(si.vec)])
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
  median.signs <- matrix(0, ntype, nsig)
  mean.signs <- matrix(0, ntype, nsig)
  q1.signs <- matrix(0, ntype, nsig)
  q3.signs <- matrix(0, ntype, nsig)
  median.exps <- matrix(0, nsig, ngenome)
  mean.exps <- matrix(0, nsig, ngenome)
  q1.exps <- matrix(0, nsig, ngenome)
  q3.exps <- matrix(0, nsig, ngenome)
  for(isig in 1:nsig){
    matrix.inst <- matrix.all.inst.signs[sign.clusters$cluster == isig,]
    median.signs[,isig] <- apply(matrix.inst, 2, median)
    mean.signs[,isig] <- apply(matrix.inst, 2, mean)
    q1.signs[,isig] <- apply(matrix.inst , 2, function(x){ median(x[x < median(x)])})
    q3.signs[,isig] <- apply(matrix.inst , 2, function(x){ median(x[x > median(x)])})
  }
  
  median.exps <- apply(object@input.matrix, 2, function(x, signs) coef(nnls(signs, x)), signs = median.signs)
  mean.exps <- apply(object@input.matrix, 2, function(x, signs) coef(nnls(signs, x)), signs = mean.signs)
  q1.exps <- apply(object@input.matrix, 2, function(x, signs) coef(nnls(signs, x)), signs = q1.signs)
  q3.exps <- apply(object@input.matrix, 2, function(x, signs) coef(nnls(signs, x)), signs = q3.signs)
  cluster.centers <- new("sign.inst")
  signs.l <- list(median = median.signs,
                mean = mean.signs,
                q1 = q1.signs,
                q3 = q3.signs )
  exps.l <- list(median = median.exps,
                 mean = mean.exps,
                 q1 = q1.exps,
                 q3 = q3.exps )
  perp.l <- list(median = median.perplexity,
                 mean = mean.perplexity,
                 q1 = q1.perplexity,
                 q3 = q3.perplexity)
  bic.l <- list(median = median.bic,
                mean = mean.bic,
                q1 = q1.bic,
                q3 = q3.bic)
  silhou.l <- list(median = median.silhou,
                mean = mean.silhou,
                q1 = q1.silhou,
                q3 = q3.silhou,
                sd = sd.silhou)

  cluster.centers@signs <- signs.l
  cluster.centers@exps <- exps.l
  cluster.centers@perp <- perp.l
  cluster.centers@bic <- bic.l
  cluster.centers@silhou.width <- silhou.l
  cluster.centers@nsig <- nsig

  return(cluster.centers)
})

#runs one topic model calculation for a single inst
#saves the signatures exposures and error 
setMethod(".run.iter", signature(object = "sign.inst"), function(object, n.iter, input.dtm, input.matrix, isig, method){
  seed.inst <- object@iter + 2000

  object <- new("sign.inst")
  object@nsig <- isig
  if(method == "ctm"){
    inst.method <- CTM(input.dtm, isig, method = "VEM", control = list( cg = (list (iter.max = 500L, tol = 10^-4)), seed =  seed.inst)) #add object@iter before 200
    object@method <- new('topic.mod.obj', ctm = inst.method)
  }else if(method == "lda"){
    inst.method <- LDA(input.dtm, isig, method = "Gibbs", control = list(seed = seed.inst)) #add object@iter before 300
    object@method <- new('topic.mod.obj', lda = inst.method)
  }else stop(sprintf('invalid method %s', object@method))
  beta<-slot(inst.method, 'beta')
  gamma<-slot(inst.method, 'gamma')
  object@signs <- list( inst = exp(t(beta)))
  object@exps <- list( inst = t(gamma))

  error <- calc.frob.error(object, input.matrix)
  object@frob.err <- error

  return(object)
})

setMethod(".get.min.inst", signature(object = "sign.stack"), function(object){
  error.vector <- unlist(sapply(object@insts, get.frob.err))
  min.ind <- which(min(error.vector) == error.vector)
  return(object@insts[[min.ind[[1]]]])
})

setMethod("run.calc", signature(object = "sign.stack"), function(object, cluster = FALSE, n.iter = 400){
  index <- 1
  input.matrix <- t(as.matrix(object@input.dtm))
  object@input.matrix <- input.matrix
 
  for(isig in seq(from = object@nsig.min, to = object@nsig.max, by = object@nsig.step)){
    print(sprintf("signature %d", isig))
    
    stack.tmp <- new("sign.stack")
    stack.tmp <- set.stack(stack.tmp, input.dtm = object@input.dtm, method = object@method, nsig.min = 1, nsig.max = n.iter)
   
    stack.tmp@input.matrix <- input.matrix
    stack.tmp@insts <- lapply(seq_along(stack.tmp@insts), function(y, i){y[[i]]@iter <- i; return(y[[i]])}, y = stack.tmp@insts) 
    
    #when you do this step multi core change n.iter
    #split it into groups rather than per iteration
    stack.tmp@insts <- lapply(stack.tmp@insts, .run.iter, n.iter = n.iter, input.dtm = object@input.dtm, input.matrix = input.matrix, isig = isig, method = object@method)
    min.inst <- .get.min.inst(stack.tmp)
    
    if(cluster){
      save(stack.tmp, file = "test.Rda")
      cluster.center <- cluster.signs(stack.tmp)
      object@insts[[index]] <- cluster.center
      object@insts[[index]]@frob.err <- calc.frob.error(cluster.center, input.matrix)
    }

    object@min.insts[[index]] <- min.inst
    object@min.insts[[index]]@frob.err <- calc.frob.error(min.inst, input.matrix)
    object@min.insts[[index]]@perp <- list(inst = calc.perplexity.inst(min.inst,  object@input.dtm, object@method))
    object@min.insts[[index]]@bic <- list(inst = calc.bic.inst(min.inst, object@input.dtm, object@method))
    object@min.insts[[index]]@silhou.width <- list(inst = 1)
    index <- index + 1
  }  

  rm(input.matrix)
  rm(min.inst) 
  rm(stack.tmp)
  return(object)
})

setMethod("calc.perplexity.inst", signature(object = "sign.inst"), function(object, new.data, method.desc){
  perplexity <- 0
  if(method.desc == "ctm") perplexity <- perplexity(object@method@ctm, new.data)
  else if(method.desc == "lda") perplexity <- perplexity(object@method@lda, new.data) 
  else stop(sprintf('invalid method %s', method.desc))
  return(perplexity)
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

find_signatures <- function(){

  nsig_max = 9
  nsig_min = 2
  
#  input_file_list <- 'test.txt'
#  dtm <- make_dtm_from_vcf_list(input_file_list)

  load('kidney_dtm_ncontext_3.Rda')

#  dtm <- make_dtm_from_vcf_list(input_file_list, ncontext = 5)
  
  signature_stack <- new("sign.stack") 
  
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
