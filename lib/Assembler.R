###########################################################################################
###### This function takes all the required scripts together to assemble a genome    ######
###### the flow is the following													                           ######
###### 1. generate random sequence                                                   ######
###### 2. assign uniform or non-uniform breakage probabilities to a k-mer            ######
######    length across the whole genome                                             ######
###### 3. generate reads biased based on these probabilities                         ######
###### 4. de novo assemble genome based on generated reads                           ######
###### 5. obtain alignment scores for each de novo generated assembly against        ######
######    reference sequence                                                         ######
###### 6. obtain alignment scores for each de novo generated assembly ONLY           ######
######    taking into account the breakage probabilities and alignment against reads ######
###########################################################################################

string.reconstruction <- function(len = len, G.cont = G.cont, 
                                  C.cont = C.cont, A.cont = A.cont, 
                                  k = k, NCPU = NCPU, 
                                  timer = timer, # max.time per DBG assembly
                                  length.of.read = length.of.read, # length of each read
                                  multiplier = multiplier, # read generation number multiplier
                                  prob.type = prob.type # "uniform" or "non-uniform"
                                  ){
  # measure execution time 
  start.time <- Sys.time()
  
  # generate random genome sequence
  randSeq <- ReadsGenerator(len = len, G.cont = G.cont, C.cont = C.cont, 
                            A.cont = A.cont, k.mer = 2, NCPU = NCPU,
                            length.of.read = length.of.read,
                            multiplier = multiplier, 
                            prob.type = prob.type)
  
  # traverse all paths in the rooted tree and construct the Eulerian paths
  output <- eulerian.path(get.balance.count(
    debrujin.graph.from.kmers(kmer.composition(k = k))), timer = timer)
  
  # obtain most plausible sequence based on breakage probability scores
  breakage.output <- breakage.prob.scoring(randSeq[[2]], output)
  # calculate the alignment scores per reconstructed genome vs. reference sequence
  align <- alignment.scoring(path = output, randSeq[[1]], NCPU = NCPU)
  
  # measure execution time
  end.time   <- Sys.time()
  time.taken <- end.time-start.time

  RESULTS                   <- NULL
  RESULTS$ref.seq           <- paste(randSeq[[1]], collapse = "") # reference sequence
  RESULTS$alignment.seq     <- align[[1]] # highest-score assembled sequence
  RESULTS$alignment.scores  <- align[[2]] # all calculated alignment scores
  RESULTS$alignment.res     <- align[[3]] # all calculated alignment scores in a df
  RESULTS$breakage.best.seq <- breakage.output[[1]] # highest breakage prob score sequence
  RESULTS$prob.ref.table    <- breakage.output[[2]] # all breakage probability scores
  RESULTS$prob.table        <- randSeq[[2]] # table of probs to bias reads generation
  RESULTS$kmer.prob.seq     <- randSeq[[6]] # all k-mers from genome and associated probs
  RESULTS$total.reads       <- randSeq[[3]] # total number of reads generated
  RESULTS$reads.freq        <- unlist(randSeq[[5]]) # read positions in the genome
  RESULTS$coverage          <- randSeq[[4]] # total coverage 
  RESULTS$execution.time    <- time.taken   # execution time to run all functions
  
  return(RESULTS)
}