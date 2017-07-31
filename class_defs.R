library(tm)
library(topicmodels)

setClass("topic.mod.obj", representation(ctm = "CTM_VEM", lda = "LDA_VEM"))

setClass("sign.inst", representation(nsig = "integer",
                                      signs= "matrix", 
                                      exps = "matrix", 
                                      method = "topic.mod.obj",
                                      perplexity = "numeric",  
                                      bic = "numeric", 
                                      frob.err = "numeric", 
                                      silhou.width = "vector",
                                      choice.met = "character"))


setGeneric("calc.frob.error", def = function(object, input.matrix){
  standardGeneric("calc.frob.error")
})

setGeneric("calc.perplexity.inst", function(object, new.data, method.desc){
  standardGeneric("calc.perplexity.inst")
})

setGeneric("calc.bic.inst", function(object, input.dtm, method.desc){
  standardGeneric("calc.bic.inst")
})

setClass("sign.stack", representation(nsig.min = "integer", 
                                     nsig.step = "integer", 
                                     nsig.max = "integer", 
                                     method = "character", 
                                     input.dtm = "DocumentTermMatrix", 
                                     types = "vector", 
                                     insts = "list",
                                     num.insts = "integer"))

setGeneric("run.calc", def = function(object, cluster = FALSE, n.iter = 100, by = "median"){
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

setGeneric("cluster.signs", function(object, by = "median"){
  standardGeneric("cluster.signs")
})