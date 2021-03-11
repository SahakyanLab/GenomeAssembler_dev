# This file randomly splits the sequences based 
# on the distribution of the read sizes

ReadSizeGenerator <- function(len = len, G.cont = G.cont, C.cont = C.cont, A.cont = A.cont){
  # generate random sequence 
  randSeq <- RandomSeqGenerator(len = len, 
                                G.cont = G.cont, C.cont = C.cont, A.cont = A.cont)
  
  # define total number of reads and average read length 
  total.reads  <- length(randSeq)*100
  # read.length  <- 10
  # read.length  <- 150
  read.length  <- 250
  # total.reads  <- ceiling((coverage*length(randSeq))/read.length)
  if(length(randSeq)<=read.length){
    stop("Choose a bigger length value!")
  }
  
  # generate all reads
  reads <- lapply(1:total.reads, function(i){
    rand.read.size <- round(rnorm(n = 1, mean = read.length, sd = 1), 0)
    max.start.pos  <- length(randSeq)-rand.read.size
    rand.start.pos <- round(runif(n = 1, min = 1, max = max.start.pos), 0)
    read           <- paste(randSeq[rand.start.pos:(rand.start.pos+rand.read.size)], 
                            collapse = "")
  })
  reads <- unlist(reads)
  
  # save all reads in a text file
  file.reads <- file("../data/reads.txt")
  writeLines(reads, file.reads, sep = "\n")
  close(file.reads)
  return(randSeq)
}