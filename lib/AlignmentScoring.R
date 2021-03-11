###########################################################################################
###### This function aligns each generated de novo genome with the reference         ######
###### sequence and calculates a variety of alignment scores 						             ######
###### including highe-order kmer, levenshtein, local-global, global-local,  		     ######
###### local, global and a final score that sums all results into account			       ######
###########################################################################################

alignment.scoring <- function(path, randSeq, NCPU = NCPU){
  if(length(path)==0){
    # if DBG assembly took too long, it will stop at X-minutes, where X is
    # the upper time limit for exploration; in such a case, it will generate
    # no assembled result, hence no pathway obtained and we return NULL
    return()
  } else {
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
    suppressPackageStartupMessages(
      packages(c("Biostrings", "foreach", "itertools", "doSNOW")))
    
    # set-up cluster for parallel computation
    cl   <- makeSOCKcluster(NCPU)
    registerDoSNOW(cl)
    
    # remove strings which contain NAs
    if(length(grep("NA", path))>0){
      path <- path[-grep("NA", path)]
    }
    
    rand.seq.combined <- paste(randSeq, collapse = "")
    ref.seq <- DNAString(rand.seq.combined)
    # obtain all higher-order k-mers in the reference genome
    k = 25
    kmer.count <- function(seq, k){
      sequence <- unlist(strsplit(as.character(seq), split = ""))
      end      <- length(sequence)-k+1
      kmer     <- substring(paste(sequence, collapse = ""),
                            first = 1:end, last = (1:end)+k-1)
      kmer     <- unlist(kmer)
      return(kmer)
    }
    kmers.ref <- kmer.count(as.character(ref.seq), k)
    
    # setup text progress bar for use in foreach loop
    # supported by the doSNOW package
    pb       <- txtProgressBar(min = 1, max = length(path), style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts     <- list(progress = progress)
    
    # calculate alignment scores of all paths and return only the sequence with max score
    print("Calculating alignment scores...", quote = F)
    alignment.runs <- foreach(i = 1:length(path), 
                              .combine = "rbind", 
                              .packages = c("foreach", "Biostrings"), 
                              .inorder = TRUE, .options.snow = opts)%dopar%{
      seq.new <- DNAString(path[[i]])
      mat     <- nucleotideSubstitutionMatrix(match    = 1, 
                                              mismatch = -3, baseOnly = TRUE)
      gap.opening   = 5 
      gap.extension = 2
      # global pairwise alignment
      globalAlign  <- 
        pairwiseAlignment(ref.seq, seq.new, 
                          type = "global", substitutionMatrix = mat, 
                          gapOpening = gap.opening, gapExtension = gap.extension, 
                          scoreOnly = TRUE)
      # local pairwise alignment
      localAlign   <- 
        pairwiseAlignment(ref.seq, seq.new, 
                          type = "local", substitutionMatrix = mat,
                          gapOpening = gap.opening, gapExtension = gap.extension, 
                          scoreOnly = TRUE)
      # global-local pairwise alignment
      global.localAlign <-
        pairwiseAlignment(ref.seq, seq.new,
                          type = "global-local", substitutionMatrix = mat,
                          gapOpening = gap.opening, gapExtension = gap.extension, 
                          scoreOnly = TRUE)
      # local-global pairwise alignment
      local.globalAlign <-
        pairwiseAlignment(ref.seq, seq.new,
                          type = "local-global", substitutionMatrix = mat,
                          gapOpening = gap.opening, gapExtension = gap.extension, 
                          scoreOnly = TRUE)
      # higher order k-meric signature counting
      kmers.new  <- kmer.count(as.character(seq.new), k)
      kmer.score <- ceiling(length(
        match(which(!is.na(match(kmers.ref, kmers.new))==TRUE),
              kmers.new))/length(kmers.ref)*100)
      # levenshtein distance; the lower the levenshtein score, the better;
      # hence need to subtract from final score
      ls.score <- as.numeric(stringDist(
        c(as.character(ref.seq), as.character(seq.new)),
        method = "levenshtein"))
      
      finalAlign <- globalAlign+localAlign+global.localAlign+
        local.globalAlign+kmer.score-ls.score
      # create data frame for all the outputs
      out.put <- data.frame(score        = finalAlign,
                            local        = localAlign,
                            global       = globalAlign,
                            global.local = global.localAlign,
                            local.global = local.globalAlign,
                            kmer         = kmer.score,
                            levenshtein  = ls.score,
                            seq          = as.character(seq.new)
                            )
      return(out.put)
    }
    close(pb)
    stopCluster(cl)  
    print("Alignment scores calculated!", quote = F)
    return(list(alignment.runs$seq[match(max(alignment.runs$score), 
                                         alignment.runs$score)], 
                alignment.runs$score, alignment.runs))
  }
}