library(slam)
library(VariantAnnotation)
library(GenomicRanges)
library(BSgenome)

library(BSgenome.Hsapiens.UCSC.hg19)

make_dtm_from_matrix <- function(input_matrix, sample_names, types){
  dtm <- as.simple_triplet_matrix(input_matrix)
  dtm$dimnames$Docs <- sample_names
  dtm$dimnames$Terms <- types
  class(dtm) <- c("DocumentTermMatrix", "simple_triplet_matrix")
  attr(dtm, "weighting") <- c("term frequency", "tf")
}

.convert_to_vector <- function(context, ref, alt, types){
  nsnv <- length(context)
  context <- tolower(context)
  ref <- tolower(ref)
  alt <- tolower(alt)
  
  if(nsnv != length(ref) || nsnv != length(alt) || nsnv != length(types))
    stop('dimensions of context ref alt are different')
  if(nsnv == 0) 
    stop('empty snv array')
  
  range <- nchar(context[[1]])
  if(nchar %% 2 == 0)
    stop('context should be an odd number')
  
  combined <- lapply(ref, paste0, alt) 
  combined <- lapply(combined, paste0, substr(context, 1, floor(nchar/2))
  combined <- lapply(combined, paste0, substr(context, range - floor(nchar/2), range) 
  
  count_vector <- rep(0, length(types))

  #binary search for counting the occurances of each type
  for(snv in combined){
    L <- 0
    R <- nsnv - 1
    m <- 0
    while(L < R){
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

.make_type <- function(ncontext = 3, nstrand = 1){
  components <- c('a', 'c', 'g', 't')
  if(nstrand == 1) base_in <- c('c', 't')
  else base_in <- components  
  types <- rep('', 3*4^(context -1)*nstrand ))
  index <- 1
  for(base in basein){
    for(alt in components[!grepl('a',components)]){
      for(5prime in components){
        for(3prime in components){     
          types[[index]] paste0(base, alt, 5prime, 3prime)
          index <- index + 1
        }
      }
    }
  }
  types <- sort(types)
}

make_vector_from_vcf <- function(vcf_file, ref_genome = BSgenome.Hsapiens.UCSC.hg19, ncontext = 3){
  #get vcf obj using genomic ranges
  vcf <- readVcf(vcf_file)
  gr <- granges(vcf)  
  
  #get a range of ncontext around snv
  gr_context <- resize(gr, ncontext, fix = 'center')
  
  #read the context around snv from the reference
  seq_start <- gr_context@ranges@start ranges
  seq_end <- gr_context@ranges@start + gr_context@ranges@width - 1
  chrom_nums<-paste0('chr',as.vector(gr_context@seqnames), seq = '')
  context_seq <- getSeq(ref_genome, chrom_nums, start = seq_start, end = seq_end, as.character = TRUE)
  
  #convert snv arrays into counts specified by types
  count_array <- .convert_to_vector(context_seq, ref, alt)
  return(count_array)
}

make_dtm_from_vcf_list <- function(vcf_file_list, ref_genome = BSgenome.Hsapiens.UCSC.hg19, ncontext = 3, nstrand = 1){
  vcf_files_table <- read.table(vcf_file_list)
  vcf_files <- vcf_files_table$V1
  matrix_snvs <- matrix(0, 3*4^(ncontext - 1)*nstrand, 0)
  for(ifile in 1:length(vcf_files)){
    matrix_snvs <- cbind(matrix_snvs, make_dtm_from_vcf(vcf_files[[ifile]], ref_genome))
  }  
  types <- .make_type(ncontext, nstrand)
  make_dtm_from_matrix(matrix_snvs, vcf_files, types)
}