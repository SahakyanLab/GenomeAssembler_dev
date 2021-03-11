###########################################################################################
###### This function aligns each generated de novo genome with all the reads         ######
###### and calculates an alignment score based on the generated breakage             ######
###### probabilities and the frequency of occurrence/overlap with that specific      ######
###### de novo assembled genome sequence                                             ######
###########################################################################################

breakage.prob.scoring <- function(randSeq, output){
  if(length(output)==0){
    # if DBG assembly took too long, it will stop at X-minutes, where X is
    # the upper time limit for exploration; in such a case, it will generate
    # no assembled result, hence no pathway obtained and we return NULL
    return()
  } else {
    # import all generatd reads 
    reads <- scan("../data/reads.txt", what = character(), sep = "", quiet = TRUE)
    # create copy of probability reference table to be reset in each iteration
    prob.ref.copy <- randSeq
    # initialise vector to store the sums of breakpoint probability occurrences
    prob.ref.sums <- vector(mode = "numeric", length = length(output))
    pb <- txtProgressBar(min = 1, max = length(output), style = 3)
    sums.index = 1
    
    print("Counting k-mer breakpoints...", quote = F)
    for(path in 1:length(output)){
      prob.ref <- prob.ref.copy
      # counter column for number of occurrence of k-mer breakages
      prob.ref$freq <- 0
      out.put <- strsplit(output[[path]], split = "")[[1]]
      for(i in 1:length(reads)){
        read.result <- gregexpr(pattern = reads[i], output[[path]], perl=TRUE)
        read.start  <- read.result[[1]][1]
        if(read.start!=(-1)){
          # extract attributes from read and reference sequence overlap
          read.length <- attr(read.result[[1]], 'match.length')
          read.end    <- read.start+read.length-1
          
          # obtain the k-mers that are broken on either end of the read
          if(read.start==1){
            # if k-mer is at the start of the genome
            kmer.start <- paste(out.put[read.start:(read.start+1)], collapse = "")
          } else {
            kmer.start <- paste(out.put[(read.start-1):read.start], collapse = "")
          }
          
          if(length(out.put)==read.end){
            # if k-mer is at the end of the genome
            kmer.end <- paste(out.put[(read.end-1):read.end], collapse = "")
          } else {
            kmer.end <- paste(out.put[read.end:(read.end+1)], collapse = "")
          }
          
          # obtain index of the above k-mers
          start.ind <- match(kmer.start, prob.ref$fwd.kmer)
          end.ind   <- match(kmer.end, prob.ref$fwd.kmer)
          
          # increase counter of given k-mer occurrence in data frame
          prob.ref$freq[start.ind] <- prob.ref$freq[start.ind]+1
          prob.ref$freq[end.ind]   <- prob.ref$freq[end.ind]+1  
        }
        if(i==length(reads)){
          prob.ref.sums[sums.index] <- sum(prob.ref$prob*prob.ref$freq)
          sums.index = sums.index+1 
        }
      }
      setTxtProgressBar(pb, path)
    }
    close(pb)
    print("All k-mer breakpoints counted!", quote = F)
    
    sorted.sums <- match(sort(prob.ref.sums, decreasing = TRUE), prob.ref.sums)
    return(list(output[sorted.sums[1]][[1]], prob.ref.sums))
  }
}
