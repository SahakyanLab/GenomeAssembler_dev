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
                           prob.type = "non-uniform" # "uniform" or "non-uniform"
                           ){
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
  suppressMessages(packages(c("gtools", "pbapply", "foreach", "doSNOW")))
  
  # generate random genome sequence 
  randSeq <- RandomSeqGenerator(len = len, G.cont = G.cont, C.cont = C.cont, 
                                A.cont = A.cont)
  
  # add dummy base to each end of the random sequence
  randSeq <- c("N", randSeq, "N")
  bases   <- c("A", "G", "C", "T")
  # obtain all possible (4**k) k-mer permutations
  two.mers <- permutations(n = length(bases), v = bases, 
                           r = k.mer, repeats.allowed = TRUE)
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
    prob <- runif(n = length(two.mers)/2, min = 0, max = 1)
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
  # associate each k-mer string with its corresponding probability value
  for(i in 1:(length(prob)-1)){
    # obtain index of forward k-mer string
    fwd.ind      <- match(ref.prob$fwd.kmer[i], ref.prob$fwd.kmer)
    
    # obtain reverse complement string
    split.string <- strsplit(two.mers[i], split = "")[[1]]
    rev.comp     <- unname(unlist(map[split.string]))
    rev.comp     <- paste(rev.comp, collapse = "")
    rev.ind      <- match(rev.comp, ref.prob$fwd.kmer)
    
    # replace NAs with probability values 
    ref.prob[fwd.ind, "prob"] <- prob[i]
    ref.prob[rev.ind, "prob"] <- prob[i]
  }
  
  # append final prob value to the dummy k-mers
  end <- dim(ref.prob)[1]
  ref.prob[(end-1):end, "prob"] <- tail(prob, n = 1)
  
  # obtain all k-mers of the randomly generated string
  end  <- length(randSeq)-k.mer+1
  kmer <- substring(paste(randSeq, collapse = ""),
                    first = 1:end, last = (1:end)+k.mer-1)
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
  cl   <- makeSOCKcluster(NCPU)
  registerDoSNOW(cl)
  
  # setup text progress bar for use in foreach loop
  # supported by the doSNOW package
  pb       <- txtProgressBar(min = 1, max = max.runs, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts     <- list(progress = progress)
  
  # Helper function to combine results and return list of vectors from foreach
  comb <- function(x, ...) {
    lapply(seq_along(x),
           function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
  }
  
  print("Obtaining all reads...", quote = F)
  reads <- foreach(i = 1:max.runs, 
                   .combine = 'comb', .multicombine = TRUE,
                   .init = list(list(), list()),
                   .inorder = FALSE, .options.snow = opts)%dopar%{
   sampling.points <- sample(x = breakpoint.positions, size = 2,
                             replace = FALSE)
   read <- paste(randSeq[(min(sampling.points)+1):max(sampling.points)],
                 collapse = "")
   if(nchar(read) >= (lower.limit) & nchar(read) <= upper.limit){
     out.put <- read
     return(list(out.put, sampling.points))
   }
  }
  close(pb)
  stopCluster(cl)
  
  reads.strings <- reads[[1]]
  reads.strings <- reads.strings[which(nchar(reads.strings)>4)]
  reads.freq    <- reads[[2]]
  reads.freq    <- reads.freq[which(nchar(reads.freq)>4)]
  coverage <- (length(reads.strings)*read.len)/length(randSeq)
  print(paste(length(reads.strings), " reads obtained!",
              "Coverage is: ", ceiling(coverage)), quote = F)
  
  # save generated reads into text file for use in genome assembly
  file.reads <- file("../data/reads.txt")
  writeLines(unlist(reads.strings), file.reads, sep = "\n")
  close(file.reads)
  
  # remove dummy bases on each end of the randomly generated sequence
  randSeq  <- randSeq[2:(length(randSeq)-1)]
  ref.prob <- ref.prob[1:(dim(ref.prob)[1]-2),]

  return(list(randSeq, ref.prob, length(reads.strings), ceiling(coverage), 
              reads.freq, kmer.from.seq))
}

# setwd("/Volumes/Paddy_Backup/ProjectBoard_Patrick/02-Proof_of_principle/lib/")
# 
# out <- ReadsGenerator(len = 1000, G.cont = 0.25, C.cont = 0.25, 
#                A.cont = 0.25, NCPU = 4, 
#                length.of.read = 10,
#                multiplier = ((50*1000)/10)*50, 
#                prob.type = "non-uniform")
# 
# out <- ReadsGenerator(len = 1000, G.cont = 0.25, C.cont = 0.25, 
#                       A.cont = 0.25, NCPU = 4, 
#                       length.of.read = 50,
#                       multiplier = ((50*1000)/50)*10, 
#                       prob.type = "non-uniform")
# 
#  multiplier = ((coverage * genome.length)/average.read.length)*10
# 