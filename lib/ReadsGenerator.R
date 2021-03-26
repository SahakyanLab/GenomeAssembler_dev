###########################################################################################
###### This file takes the randomly generated genome and generates reads of varying  ######
###### sizes to mimic the DNA sonication process; the flow is the following 		     ######
###### 1. assign uniform or non-uniform breakage probabilities to a k-mer            ######
######    length across the whole genome                                             ######
###### 2. generate reads biased based on these probabilities                         ######
###########################################################################################

ReadsGenerator <- function(len = len, G.cont = G.cont, C.cont = C.cont, 
                           A.cont = A.cont, k.mer = 2, NCPU = NCPU,
                           length.of.read = length.of.read, # length of each read
                           multiplier = multiplier, # read generation number multiplier
                           prob.type = "non-uniform", # "uniform" or "non-uniform"
                           # ind = i, 
                           k = k,
                           seed = seed
                           ){
  set.seed(seed)
  source("../lib/RandomSeqGenerator.R")
  packages <- function(x){
    for(i in x){
      if(!require(i, character.only = TRUE)){
        install.packages(i, dependencies = TRUE)
        library(i, character.only = TRUE)
      } else {
        library(i, character.only = TRUE)
      }
    }
  }
  suppressMessages(packages(c("gtools", "pbapply", "foreach", "doRNG", "doParallel")))
  
  # generate random genome sequence 
  randSeq <- RandomSeqGenerator(len = len, G.cont = G.cont, C.cont = C.cont, 
                                A.cont = A.cont, seed = seed)
  
  reverseComp <- function(x){
    paste(unname(unlist(map[unlist(
      lapply(strsplit(x, NULL), rev))])), collapse = "")
  }
  
  # add dummy base to each end of the random sequence
  randSeq <- c("N", randSeq, "N")
  bases   <- c("A", "G", "C", "T")
  # obtain all possible (4**k) k-mer permutations
  two.mers <- permutations(n = length(bases), v = bases, 
                           r = 2, repeats.allowed = TRUE)
  # concatenate k-mers into one string
  two.mers <- sapply(1:dim(two.mers)[1], function(i){
    paste(two.mers[i,], collapse = "")
  })
  
  # initialise character mapping list to obtain reverse complements
  map  <- list("A" = "T", "T" = "A", "C" = "G", "G" = "C")
  # initialise probability vector of length (4**k)/2 
  if(prob.type=="uniform"){
    prob <- rep(1/length(two.mers), length(two.mers)/2)
    prob <- c(prob, mean(prob)) # avg prob to end-bases
    prob <- prob/sum(prob) # normalise all probs to sum of one
    prob <- prob/2
  }
  if (prob.type=="non-uniform"){
    prob <- runif(n = 10, min = 0, max = 1)
    prob <- c(prob, mean(prob)) # avg prob to end-bases
    prob <- prob/sum(prob) # normalise all probs to sum of one
    prob <- prob/2
  }
  
  # add dummy k-mers from each end of randSeq onto k-mer vector
  two.mers <- c(two.mers,
                paste(randSeq[1:2], collapse = ""),
                paste(tail(randSeq, n = 2), collapse = ""))
  
  # initialise reference data frame with k-mer and associated probabilities
  ref.prob <- data.frame("fwd.kmer" = two.mers,
                         "prob"     = NA)
  # pre-fill probabilities for dummy k-mers
  dummy.ind <- grep("N", ref.prob$fwd.kmer)
  ref.prob$prob[dummy.ind] <- tail(prob, n = 1)
  prob <- prob[1:(length(prob)-1)]
  
  # associate each k-mer string with its corresponding probability value
  i = 1
  while(length(prob)>0){
    # obtain index of forward k-mer string
    fwd.ind      <- match(ref.prob$fwd.kmer[i], ref.prob$fwd.kmer)
    
    if(is.na(ref.prob$prob[i])){
      # obtain reverse complement string
      rev.comp     <- reverseComp(two.mers[i])
      rev.ind      <- match(rev.comp, ref.prob$fwd.kmer)
      
      # replace NAs with probability values 
      ref.prob[fwd.ind, "prob"] <- prob[1]
      ref.prob[rev.ind, "prob"] <- prob[1]
      prob <- prob[-1] 
    } 
    i = i+1
  }
  
  # obtain all k-mers of the randomly generated string
  end  <- length(randSeq)-2+1
  kmer <- substring(paste(randSeq, collapse = ""),
                    first = 1:end, last = (1:end)+2-1)
  # initialise data frame with all k-mers from ref.seq and associated probs
  kmer.from.seq <- data.frame("kmer" = kmer,
                              "prob" = NA)
  # replace NAs with probability values from reference data frame
  kmer.ind           <- match(kmer.from.seq$kmer, ref.prob$fwd.kmer)
  kmer.from.seq$prob <- ref.prob$prob[kmer.ind]
  
  # sample 10k breakpoint positions at a time
  len.breakpoint.positions <- length(randSeq)*1000
  len.sampling             <- 10000
  
  print(paste("Sampling",  len.breakpoint.positions, 
              "breakpoint positions..."), quote = F)
  breakpoint.positions <- pbreplicate(n = len.breakpoint.positions/len.sampling, {
    sampling <- kmer.from.seq[sample(x = nrow(kmer.from.seq), size = len.sampling,
                                     replace = TRUE, prob = kmer.from.seq$prob),]
    sampling <- round(as.integer(rownames(sampling)), 0)
  }, simplify = TRUE)
  print("Breakpoint positions sampled!", quote = F)
  
  # generate all reads
  read.len    <- length.of.read
  three.sd    <- (read.len*0.1)*3 # 1.SD = 10% of read length
  upper.limit <- read.len+three.sd
  lower.limit <- read.len-three.sd
  max.runs    <- read.len*multiplier
  end         <- length(randSeq)-1
  # set-up cluster for parallel computation
  NCPU <- NCPU
  cl <- makeCluster(NCPU)
  registerDoParallel(cl)
  registerDoRNG(seed = seed)
  
  # Helper function to combine results and return list of vectors from foreach
  comb <- function(x, ...) {
    lapply(seq_along(x),
           function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
  }
  
  print("Obtaining all reads...", quote = F)
  reads <- foreach(i = 1:max.runs, 
                   .combine = 'comb', .multicombine = TRUE,
                   .init = list(list(), list()),
                   .inorder = FALSE)%dopar%{
   sampling.points <- sample(x = breakpoint.positions, size = 2,
                             replace = FALSE)
   read <- paste(randSeq[(min(sampling.points)+1):max(sampling.points)],
                 collapse = "")
   if(nchar(read) >= (lower.limit) & nchar(read) <= upper.limit){
     out.put <- read
     read.length <- max(sampling.points)-min(sampling.points)
     return(list(out.put, read.length))
   }
  }
  stopCluster(cl)
  
  reads.strings <- reads[[1]]
  reads.strings <- reads.strings[which(nchar(reads.strings)>4)]
  reads.freq    <- unlist(reads[[2]])
  coverage <- (length(reads.strings)*read.len)/length(randSeq)
  print(paste(length(reads.strings), " reads obtained!",
              "Coverage is: ", ceiling(coverage)), quote = F)
  
  # save generated reads into text file for use in genome assembly
  file.reads <- file("../data/reads.txt")
  writeLines(unlist(reads.strings), file.reads, sep = "\n")
  close(file.reads)
  
  # write.csv(unlist(reads.strings),
  #           file = paste0("../data/reads-",ind,".csv"),
  #           row.names = FALSE)
   
  # remove dummy bases on each end of the randomly generated sequence
  randSeq  <- randSeq[2:(length(randSeq)-1)]
  ref.prob <- ref.prob[1:(dim(ref.prob)[1]-2),]

  # write.csv(paste(randSeq, collapse = ""),
  #           file = paste0("../data/ref-seq-",ind,".csv"),
  #           row.names = FALSE)
  
  # generating k-mers from the reads
  reads <- unlist(reads.strings)
  # obtain all k-mers per read
  kmers <- lapply(1:length(reads), function(i){
    sequence <- unlist(strsplit(reads[i], split = ""))
    end      <- length(sequence)-k+1
    kmer     <- substring(paste(sequence, collapse = ""),
                          first = 1:end, last = (1:end)+k-1)
  })
  kmers   <- unlist(kmers)
  # obtain weighted k-mers
  kmer.df <- as.data.frame(table(kmers))
  kmers   <- unlist(lapply(kmer.df$kmers, as.character))
  newList <- list(kmers, k)

  return(list(randSeq, ref.prob, length(reads.strings), ceiling(coverage),
              reads.freq, kmer.from.seq, newList))
}
#
# setwd("/Volumes/Paddy_Backup/ProjectBoard_Patrick/02-Proof_of_principle/lib/")
#
# 
# out <- ReadsGenerator(len = 1000, G.cont = 0.25, C.cont = 0.25,
#                       A.cont = 0.25, NCPU = 4,
#                       length.of.read = 100,
#                       multiplier = 53, k = 9,
#                       prob.type = "non-uniform", seed = 1234)