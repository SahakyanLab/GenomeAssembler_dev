# This file randomly splits the sequences based 
# on the distribution of the read sizes

ReadsGenerator <- function(len = len, G.cont = G.cont, C.cont = C.cont, 
                           A.cont = A.cont, NCPU = NCPU, k.mer = k.mer,
                           prob.type = "non-uniform" # "uniform" or "non-uniform"
                           ){
  source("RandomSeqGenerator.R")
  # load required packages or install if not installed
  packages <- function(x){
    for(i in x){
      if(!require(i, character.only = TRUE)){
        install.packages(i, dependencies = TRUE)
        suppressMessages(library(i, character.only = TRUE))
      } else {
        suppressMessages(library(i, character.only = TRUE))
      }
    }
  }
  suppressMessages(packages(c("gtools", "pbapply", "doSNOW")))
  
  # generate random genome sequence 
  randSeq <- RandomSeqGenerator(len = 1000, G.cont = 0.3, C.cont = 0.3, A.cont = 0.2)
  bases <- c("A", "G", "C", "T")
  # obtain all possible (4**k) k-mer permutations
  two.mers <- permutations(n = length(bases), v = bases, 
                           r = k.mer, repeats.allowed = TRUE)
  # concatenate k-mers into one string
  two.mers <-sapply(1:dim(two.mers)[1], function(i){
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
  kmer.from.seq <- data.frame("kmer" = kmer,
                              "prob" = NA)
  # replace NAs with probability values from reference data frame
  kmer.ind <- match(kmer.from.seq$kmer, ref.prob$fwd.kmer)
  kmer.from.seq$prob <- ref.prob$prob[kmer.ind]
  
  # sample 100 breakpoint positions at a time
  len.breakpoint.positions <- length(randSeq)*100000
  len.sampling <- 100
  
  print("sampling breakpoint positions...", quote = F)
  breakpoint.positions <- pbsapply(1:(len.breakpoint.positions/len.sampling), function(i){
    sampling <- kmer.from.seq[sample(x = nrow(kmer.from.seq), size = len.sampling, 
                                     replace = TRUE, prob = kmer.from.seq$prob),]
    sampling <- round(as.numeric(rownames(sampling)), 0) 
  })
  print("breakpoint positions sampled!", quote = F)
  
  # generate all reads
  read.len    <- 100
  avg         <- read.len
  three.sd    <- (read.len*0.1)*3
  upper.limit <- avg+three.sd
  lower.limit <- avg-three.sd
  max.runs    <- read.len*5000
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
  
  print("Obtaining all reads...", quote = F)
  reads <- foreach(i=1:max.runs, .combine = "rbind", 
                   .inorder = FALSE, .options.snow=opts)%dopar%{
    sampling.points <- sample(x = breakpoint.positions, size = 2, 
                              replace = FALSE)
    if(sampling.points[1]==1 | sampling.points[2]==1){
     # for start of genome string
     read <- paste(randSeq[min(sampling.points):max(sampling.points)], 
                   collapse = "")
    } else if((sampling.points[1]==end | sampling.points[2]==end)){
     # for end of genome string
     read <- paste(randSeq[(min(sampling.points)+1):(max(sampling.points)+1)], 
                   collapse = "")
    } else {
     read <- paste(randSeq[(min(sampling.points)+1):max(sampling.points)], 
                   collapse = "")
    }
    if(nchar(read) >= (lower.limit) & nchar(read) <= upper.limit){
     out.put <- read
     return(out.put)
    }              
  }
  close(pb)
  stopCluster(cl)
  print("Reads obtained!", quote = F)
  
  # save generated reads into text file for use in genome assembly
  file.reads <- file("../data/reads.txt")
  writeLines(reads, file.reads, sep = "\n")
  close(file.reads)
  
  return(list(randSeq, ref.prob))
}

# randSeq <- ReadsGenerator(len = 1000, G.cont = 0.3, C.cont = 0.3,
#                           A.cont = 0.2, NCPU = 2, k.mer = 2, prob.type = "non-uniform")
