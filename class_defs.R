library(tm)
library(topicmodels)

setClass("topic.mod.obj", representation(ctm = "CTM_VEM", lda = "LDA_Gibbs"))

setClass("sign.inst", representation(nsig = "integer",
                                     signs= "list", 
                                     exps = "list",  
                                     method = "topic.mod.obj",
                                     perp = "list", 
                                     bic = "list",
                                     frob.err = "list", 
                                     silhou.width = "list"))

setGeneric("calc.frob.error", def = function(object, input.matrix){
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
})

setGeneric("get.signs", function(object, by = "inst"){
  standardGeneric("get.signs")
})

setGeneric("get.exps", function(object, by = "inst"){
  standardGeneric("get.exps")
})

setGeneric("get.frob.err", function(object, by = "inst"){
  standardGeneric("get.frob.err")
})


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