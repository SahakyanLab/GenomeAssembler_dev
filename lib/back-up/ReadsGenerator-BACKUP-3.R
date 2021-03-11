# This file randomly splits the sequences based 
# on the distribution of the read sizes

ReadsGenerator <- function(len = len, G.cont = G.cont, C.cont = C.cont, 
                           A.cont = A.cont, k.mer = 2, NCPU = NCPU, 
                           length.of.read = length.of.read, # length of each read
                           multiplier = multiplier, # read generation number multiplier
                           prob.type = "non-uniform" # "uniform" or "non-uniform"
){
  source("../lib/RandomSeqGenerator.R")
  # load required packages or install if not installed
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
  suppressMessages(packages(c("gtools", "foreach", "pbapply", "doSNOW")))
  
  # generate random genome sequence 
  randSeq <- RandomSeqGenerator(len = len, G.cont = G.cont, C.cont = C.cont, 
                                A.cont = A.cont)
  bases <- c("A", "G", "C", "T")
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
  } else if (prob.type=="non-uniform"){
    prob <- runif(n = length(two.mers)/2, min = 0, max = 1)
    prob <- prob/sum(prob)
    prob <- prob/2
  }
  
  # initialise data frame with k-mer and associated probabilities
  # reference table 
  ref.prob <- data.frame("fwd.kmer" = two.mers,
                         "prob"     = NA)
  # associate each k-mer string with its corresponding probability value
  for(i in 1:length(prob)){
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
  
  # obtain all k-mers of the randomly generated string
  end  <- length(randSeq)-k.mer+1
  kmer <- substring(paste(randSeq, collapse = ""),
                    first = 1:end, last = (1:end)+k.mer-1)
  # initialise data frame with k-mer and associated probabilities
  # table for all k-mers from randomly generated sequence
  kmer.from.seq <- data.frame("kmer" = kmer,
                              "prob" = NA)
  # replace NAs with probability values from reference data frame
  kmer.ind <- match(kmer.from.seq$kmer, ref.prob$fwd.kmer)
  kmer.from.seq$prob <- ref.prob$prob[kmer.ind]
  
  # sample 10k breakpoint positions at a time
  len.breakpoint.positions <- length(randSeq)*1000
  len.sampling <- 10000
  
  print(paste("Sampling",  len.breakpoint.positions, "breakpoint positions..."), quote = F)
  breakpoint.positions <- pbreplicate(n = len.breakpoint.positions/len.sampling, {
    sampling <- kmer.from.seq[sample(x = nrow(kmer.from.seq), size = len.sampling,
                                     replace = TRUE, prob = kmer.from.seq$prob),]
    sampling <- round(as.numeric(rownames(sampling)), 0)
  }, simplify = TRUE)
  print("Breakpoint positions sampled!", quote = F)
  
  # generate all reads
  read.len    <- length.of.read
  three.sd    <- (read.len*0.1)*3 # 1 SD = 10% of read length
  upper.limit <- read.len+three.sd
  lower.limit <- read.len-three.sd
  max.runs    <- read.len*multiplier
  end         <- length(randSeq)-1
  # set-up cluster for parallel computation
  NCPU <- NCPU
  cl   <- makeSOCKcluster(NCPU)
  registerDoSNOW(cl)
  
  # setup text progress bar for use in foreach loop
  # pb       <- txtProgressBar(min = 1, max = max.runs, style = 3)
  # progress <- function(n) setTxtProgressBar(pb, n)
  # opts     <- list(progress = progress)
  # 
  # print("Obtaining all reads...", quote = F)
  # reads <- foreach(i=1:max.runs, .combine = "rbind",
  #                  .inorder = FALSE, .options.snow=opts)%dopar%{
  #  sampling.point <- sample(x = breakpoint.positions, size = 1,
  #                           replace = FALSE)
  #  if((sampling.point >= lower.limit) & (sampling.point <= upper.limit)){
  #    # check if breakpoint is within acceptable distance from start
  #    read <- paste(randSeq[1:sampling.point],
  #                  collapse = "")
  #    if(nchar(read) >= (lower.limit) & nchar(read) <= upper.limit){
  #      out.put <- read
  #      if(rbinom(n = 1, size = 1, prob = 0.5)==1){
  #        return(out.put)
  #      }
  #    }
  #  }
  #  if(((length(randSeq)-sampling.point) >= lower.limit) &
  #     ((length(randSeq)-sampling.point) <= upper.limit)){
  #    # check if breakpoint is within acceptable distance to the end
  #    read <- paste(randSeq[sampling.point:length(randSeq)],
  #                  collapse = "")
  #    if(nchar(read) >= (lower.limit) & nchar(read) <= upper.limit){
  #      out.put <- read
  #      if(rbinom(n = 1, size = 1, prob = 0.5)==1){
  #        return(out.put)
  #      }
  #    }
  #  }
  #  if(sampling.point > upper.limit &
  #     ((length(randSeq)-sampling.point) > upper.limit)){
  #    sampling.point.two <- sample(x = breakpoint.positions, size = 1,
  #                                 replace = FALSE)
  #    read <- paste(randSeq[(sampling.point+1):sampling.point.two],
  #                  collapse = "")
  #    if(nchar(read) >= (lower.limit) & nchar(read) <= upper.limit){
  #      out.put <- read
  #      if(rbinom(n = 1, size = 1, prob = 0.5)==1){
  #        return(out.put)
  #      }
  #    }
  #  }
  # }
  # close(pb)
  # stopCluster(cl)
  # coverage <- (length(reads)*read.len)/length(randSeq)
  # print(paste(length(reads), " reads obtained!",
  #             "Coverage is: ", ceiling(coverage)), quote = F)
  
  ################################################################################
  ################################## VERSION 2 ###################################
  ################################################################################
  # setup text progress bar for use in foreach loop
  # supported by the doSNOW package
  pb       <- txtProgressBar(min = 1, max = (max.runs/2), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts     <- list(progress = progress)
  
  print("Obtaining all reads...", quote = F)
  
  # Test helper function
  comb <- function(x, ...) {
    lapply(seq_along(x),
           function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
  }
  
  reads.one <- foreach(i=1:(max.runs/2), 
                       .combine='comb', .multicombine=TRUE,
                       .init=list(list(), list()),
                       # .combine = "rbind",
                       .inorder = FALSE, .options.snow = opts)%dopar%{
   sampling.points <- sample(x = breakpoint.positions, size = 2,
                             replace = FALSE)
   read <- paste(randSeq[(min(sampling.points)+1):max(sampling.points)],
                 collapse = "")
   if(nchar(read) >= (lower.limit) & nchar(read) <= upper.limit){
     if(rbinom(n = 1, size = 1, prob = 0.5)==1){
       out.put <- read
       # return(out.put)
       return(list(out.put, sampling.points))
     }
   }
  }
  close(pb)
  stopCluster(cl)
  
  NCPU <- NCPU
  cl   <- makeSOCKcluster(NCPU)
  registerDoSNOW(cl)
  reads.two <- foreach(i=1:(max.runs/2), 
                       .combine='comb', .multicombine=TRUE,
                       .init=list(list(), list()),
                       # .combine = "rbind",
                       .inorder = FALSE, .options.snow = opts)%dopar%{
   sampling.point <- sample(x = breakpoint.positions, size = 1,
                            replace = FALSE)
   if((sampling.point >= lower.limit) & (sampling.point <= upper.limit)){
     # check if breakpoint is within acceptable distance from start
     read <- paste(randSeq[1:sampling.point],
                   collapse = "")
     if(nchar(read) >= (lower.limit) & nchar(read) <= upper.limit){
       if(rbinom(n = 1, size = 1, prob = 0.5)==1){
         out.put <- read
         # return(out.put)
         return(list(out.put, sampling.point))
       }
     }
   } else if(((length(randSeq)-sampling.point) >= lower.limit) &
             ((length(randSeq)-sampling.point) <= upper.limit)){
     # check if breakpoint is within acceptable distance to the end
     read <- paste(randSeq[sampling.point:length(randSeq)],
                   collapse = "")
     if(nchar(read) >= (lower.limit) & nchar(read) <= upper.limit){
       if(rbinom(n = 1, size = 1, prob = 0.5)==1){
         out.put <- read
         # return(out.put)
         return(list(out.put, sampling.point))
       }
     }
   }
  }
  close(pb)
  stopCluster(cl)
  # reads <- c(reads.one, reads.two)
  reads      <- c(reads.one[[1]], reads.two[[1]])
  reads      <- reads[which(nchar(reads)>4)]
  reads.freq <- c(reads.one[[2]], reads.two[[2]])
  reads.freq <- reads.freq[which(nchar(reads.freq)>4)]
  cat("\n")
  coverage <- (length(reads)*read.len)/length(randSeq)
  print(paste(length(reads), " reads obtained!",
              "Coverage is: ", ceiling(coverage)), quote = F)
  ################################################################################
  ################################################################################
  # save generated reads into text file for use in genome assembly
  file.reads <- file("../data/reads.txt")
  writeLines(unlist(reads), file.reads, sep = "\n")
  close(file.reads)
  
  return(list(randSeq, ref.prob, length(reads), ceiling(coverage), 
              reads.freq, kmer.from.seq))
}

# setwd("/Volumes/Paddy_Backup/ProjectBoard_Patrick/02-Proof_of_principle/lib/")
# randSeq <- ReadsGenerator(len = 600, G.cont = 0.3, C.cont = 0.3,
#                           A.cont = 0.2, k.mer = 2, NCPU = 2,
#                           length.of.read = 100,
#                           multiplier = 50,
#                           prob.type = "non-uniform")