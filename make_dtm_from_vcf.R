library(slam)
library(VariantAnnotation)
library(GenomicRanges)
library(BSgenome)

library(BSgenome.Hsapiens.UCSC.hg19)

##flips a single base
.flip_base <- function(base){
  switch(base, 
  a = {
    return('t')
  },
  c = {
    return('g')    
  },
  g = {
    return('c')
  },
  t = {
    return('a')
  },{
    stop(sprintf('invalid base input %s', base))
  })
}

##flips the bases in a combination to the opposite strand
.flip_strand <- function(comb_in){
  comb_out <- rep('', length(comb_in))
  for(isnv in 1:length(comb_in)){
    snv <- comb_in[[isnv]]
    ref <- substr(snv, 1, 1)
    if( ref == 'c' || ref == 't'){
      comb_out[[isnv]] <- snv
      next
    }else{
      ref <- .flip_base(ref)
      alt <- .flip_base(substr(snv, 2, 2))
      ext <- substr(snv, 3, nchar(snv))

      #every base after the first two belongs to context
      for(iext in 1:nchar(ext)){
        replacement <- .flip_base(substr(ext, iext, iext))
        if(!exists('flip_ext')) flip_ext <- replacement 
        else flip_ext <- paste0(flip_ext, replacement)
      }

      comb_out[[isnv]] <- paste0(ref, alt, flip_ext)
      rm(flip_ext)
    }
  }

  return(comb_out)
}

##given nmer context, ref, alt information calculates 
##counts for each type of snv
.convert_seq_to_vector <- function(context, ref_vector, alt_vector, types, nstrand = 1){
  
  context <- tolower(context)
  ref_vector <- tolower(ref_vector)
  alt_vector <- tolower(alt_vector)
  nsnv <- length(context)

  if(nsnv != length(ref_vector) || nsnv != length(alt_vector))
    stop('dimensions of context ref alt are different')
  if(nsnv == 0) 
    stop('empty snv array')
  
 
  range <- nchar(context[[1]])
  if(range %% 2 == 0)
    stop('context should be an odd number')
  combined <- paste0(ref_vector, alt_vector) 
  combined <- paste0(combined, substr(context, 1, floor(range/2)))
  combined <- paste0(combined, substr(context, ceiling(range/2)+1, range))

  if(nstrand == 1) combined <- .flip_strand(combined)
  

  count_vector <- rep(0, length(types))

  #binary search for counting the occurances of each type
  for(snv in combined){
    L <- 1
    R <- length(types)
    m <- 1
    while(L <= R){
      m <- floor((L + R)/2)
      if(types[[m]] < snv) L <- m+1
      else if (types[[m]] > snv) R <- m-1
      else{
        count_vector[[m]] <- count_vector[[m]] + 1
        break
      }
    }
  }
  return(count_vector)
}

##given nmer size and strand choice returns 
## an array of types of snvs
.make_type <- function(ncontext = 3, nstrand = 1){
  components <- c('a', 'c', 'g', 't')
  base_in <- ''
  if(nstrand == 1) base_in <- c('c', 't')
  else base_in <- components  
#  types <- rep('', 6*4^(ncontext -1)*nstrand )
 
  nextend <- (ncontext - 1)
  index <- 1
  
  #first construct types without context
  #first factor is count of possible alt
  #second factor is count of possible ref
  types <- rep ('', (length(components) -1) * length(base_in)) 

  for(base in base_in){
    alt_out <- components[!grepl(base, components)]
    for(alt in alt_out){
      types[[index]] <- paste0(base, alt)
      index <- index + 1
    }
  }

  types <- rep(types, length(components)^nextend) 
  types <- sort(types)
 
  for(iextend in 1:nextend){
    
    length_rep <- length(components)^(nextend - iextend)

    nrep <- length(components)^(iextend - 1) * (length(components) -1) * length(base_in)

    for(irep in 1:nrep){
      for(component in components){
        if(!exists("types_tmp")) types_tmp <- rep(component, length_rep)
        else types_tmp <- c(types_tmp, rep(component, length_rep))
      }
    }
    for(index in 1:length(types)){
      types[[index]] <- paste0(types[[index]], types_tmp[[index]])
    }
    rm(types_tmp)
  }

  types <- sort(types)
  return(types)
}

##converts a matrix to a document term matrix
make_dtm_from_matrix <- function(input_matrix, sample_names, types){
  dtm <- as.simple_triplet_matrix(t(input_matrix))
  dtm$dimnames$Docs <- as.character(1:length(sample_names))
  #dtm$dimnames$Docs <- sample_names
  #print(types)
  dtm$dimnames$Terms <- types

  class(dtm) <- c("DocumentTermMatrix", "simple_triplet_matrix")
  attr(dtm, "weighting") <- c("term frequency", "tf")
  return(dtm)
}

##converts a single vcf to a vector of counts of types
make_vector_from_vcf <- function(vcf_file, ref_genome = BSgenome.Hsapiens.UCSC.hg19, types, ncontext = 3, nstrand = 1){
  #get vcf obj using genomic ranges
  vcf <- readVcf(vcf_file)
  gr <- granges(vcf)  
  
  #get a range of ncontext around snv
  gr_context <- resize(gr, ncontext, fix = 'center')
  
  #read the context around snv from the reference
  seq_start <- start(gr_context)
  seq_end <- end(gr_context)

  chrom_nums <- paste0('chr',as.vector(seqnames(gr_context)))
  context_seq <- getSeq(ref_genome, names = chrom_nums, start = seq_start, end = seq_end, as.character = TRUE)

  #get ref and alt
  ref_vector <- as.character(ref(vcf))
  alt_vector <- as.character(unlist(alt(vcf)))

  #convert snv arrays into counts specified by types
  count_vector <- .convert_seq_to_vector(context_seq, ref_vector, alt_vector, types, nstrand)
  return(count_vector)
}

##given a list of vcf files returns a matrix
make_matrix_from_vcf_list <- function(vcf_file_list, ref_genome = BSgenome.Hsapiens.UCSC.hg19, ncontext = 3, nstrand = 1){
  vcf_files_table <- read.table(vcf_file_list)
  vcf_files <- as.character(vcf_files_table$V1)
  matrix_snvs <- matrix(0, 6*4^(ncontext - 1)*nstrand, length(vcf_files)) 
  types <- .make_type(ncontext, nstrand)
  for(ifile in 1:length(vcf_files)){
    count_vector <- make_vector_from_vcf(vcf_files[[ifile]], ref_genome, types, ncontext, nstrand)
    matrix_snvs[, ifile] <- count_vector
  }  
  rownames(matrix_snvs) <- types
  colnames(matrix_snvs) <- vcf_files
  return(matrix_snvs)
}

##given a list of vcf files returns a document term matrix
make_dtm_from_vcf_list <- function(vcf_file_list, ref_genome = BSgenome.Hsapiens.UCSC.hg19, ncontext = 3, nstrand = 1){
  vcf_files_table <- read.table(vcf_file_list)
  vcf_files <- as.character(vcf_files_table$V1)
  matrix_snvs <- matrix(0, 6*4^(ncontext - 1)*nstrand, length(vcf_files))
  types <- .make_type(ncontext, nstrand)
  for(ifile in 1:length(vcf_files)){
    count_vector <- make_vector_from_vcf(vcf_files[[ifile]], ref_genome, types, ncontext, nstrand)
    matrix_snvs[, ifile] <- count_vector
  }  
  print('end of make dtm from vcf list')
  #print(types)
  return(make_dtm_from_matrix(matrix_snvs, vcf_files, types))
}

