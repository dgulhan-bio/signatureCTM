library(tm)
library(topicmodels)

setClass("topic.mod.obj", representation(ctm = "CTM_VEM", lda = "LDA_VEM"))

setClass("sign.miner", representation(nsig = "integer",
                                      signs= "matrix", 
                                      exps = "matrix", 
                                      method = "topic.mod.obj",
                                      perplexity = "numeric",  
                                      bic = "numeric", 
                                      frob.err = "numeric", 
                                      silhou.width = "vector",
                                      choice.met = "character"))

setGeneric("set.miners", def = function(object, nsig, signs, exps){
  standardGeneric("set.miners")
})

setGeneric("calc.frob.error", def = function(object, input.matrix){
  standardGeneric("calc.frob.error")
})

setGeneric("calc.perplexity.miner", function(object, new.data, method.desc){
  standardGeneric("calc.perplexity.miner")
})

setGeneric("calc.bic.miner", function(object, input.dtm, method.desc){
  standardGeneric("calc.bic.miner")
})

setClass("sign.mine", representation(nsig.min = "integer", 
                                     nsig.step = "integer", 
                                     nsig.max = "integer", 
                                     method = "character", 
                                     input.dtm = "DocumentTermMatrix", 
                                     types = "vector", 
                                     miners = "list",
                                     num.miners = "integer"))

setGeneric("run.mine", def = function(object){
  standardGeneric("run.mine")
})

setGeneric("set.mine", def = function(object, input.dtm, nsig.min, nsig.max, nsig.step = as.integer(1), method = "ctm"){
  standardGeneric("set.mine")
})

setGeneric("types", signature(object = "sign.mine"), def = function(object){
  standardGeneric("types")
})

setGeneric("genomes", signature(object = "sign.mine"), def = function(object){
  standardGeneric("genomes")
})

setGeneric("calc.perplexity", function(object, new.data){
  standardGeneric("calc.perplexity")
})

setGeneric("calc.bic", function(object){
  standardGeneric("calc.bic")
})