library(tm)
library(topicmodels)

setClass("topic.mod.obj", representation(ctm = "CTM_VEM", lda = "LDA_Gibbs"))

setClass("sign.inst", representation(nsig = "integer",
                                      signs= "matrix", 
                                      exps = "matrix", 
                                      mean.signs= "matrix", 
                                      mean.exps = "matrix", 
                                      q1.signs= "matrix", 
                                      q1.exps = "matrix", 
                                      q3.signs= "matrix", 
                                      q3.exps = "matrix", 
                                      method = "topic.mod.obj",
                                      perplexity = "numeric",  
                                      perp.mean = "numeric",
                                      perp.median = "numeric",
                                      perp.sd = "numeric",
                                      perp.q1 = "numeric",
                                      perp.q3 = "numeric",
                                      bic = "numeric", 
                                      bic.mean = "numeric",
                                      bic.median = "numeric",
                                      bic.sd = "numeric",
                                      bic.q1 = "numeric",
                                      bic.q3 = "numeric",
                                      frob.err = "numeric", 
                                      frob.err.mean = "numeric", 
                                      frob.err.q1 = "numeric", 
                                      frob.err.q3 = "numeric", 
                                      silhou.width = "vector",
                                      silhou.width.sd = "vector",
                                      silhou.width.mean = "vector",
                                      silhou.width.q1 = "vector",
                                      silhou.width.q3 = "vector",
                                      choice.met = "character"))


setGeneric("calc.frob.error", def = function(object, input.matrix, by = "median"){
  standardGeneric("calc.frob.error")
})

setGeneric("calc.perplexity.inst", function(object, new.data, method.desc){
  standardGeneric("calc.perplexity.inst")
})

setGeneric("calc.bic.inst", function(object, input.dtm, method.desc){
  standardGeneric("calc.bic.inst")
})

setGeneric(".same.dim", function(object, nsig.comp){
  standardGeneric(".same.dim")
}

setGeneric("get.signs", function(object)){
  standardGeneric("get.signs")
}

setGeneric("get.exps", function(object)){
  standardGeneric("get.exps")
}

setGeneric("get.frob.err", function(object)){
  standardGeneric("get.frob.error")
}


setClass("sign.stack", representation(nsig.min = "integer", 
                                     nsig.step = "integer", 
                                     nsig.max = "integer", 
                                     method = "character", 
                                     input.dtm = "DocumentTermMatrix", 
                                     input.matrix = "matrix",
                                     types = "vector", 
                                     insts = "list",
                                     min.insts = "list",
                                     num.insts = "integer"))

setGeneric("run.calc", def = function(object, cluster = FALSE, n.iter = 100){
  standardGeneric("run.calc")
})

setGeneric("set.stack", def = function(object, input.dtm, nsig.min, nsig.max, nsig.step = as.integer(1), method = "ctm"){
  standardGeneric("set.stack")
})

setGeneric("types", signature(object = "sign.stack"), def = function(object){
  standardGeneric("types")
})

setGeneric("genomes", signature(object = "sign.stack"), def = function(object){
  standardGeneric("genomes")
})

setGeneric("calc.perplexity", function(object, new.data){
  standardGeneric("calc.perplexity")
})

setGeneric("calc.bic", function(object){
  standardGeneric("calc.bic")
})

setGeneric("cluster.signs", function(object){
  standardGeneric("cluster.signs")
})