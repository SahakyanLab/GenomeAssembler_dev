###########################################################################################
###### This function takes all the required scripts together to assemble a genome    ######
###### and runs a correlation coefficient analysis on various alignment metrics      ######
###### the flow is the following                                                     ######
###### 1. generate random sequence                                                   ######
###### 2. assign uniform or non-uniform breakage probabilities to a k-mer            ######
######    length across the whole genome                                             ######
###### 3. generate reads biased based on these probabilities                         ######
###### 4. de novo assemble genome based on generated reads                           ######
###### 5. obtain alignment scores for each de novo generated assembly against        ######
######    reference sequence                                                         ######
###### 6. obtain alignment scores for each de novo generated assembly ONLY           ######
######    taking into account the breakage probabilities and alignment against reads ######
###### 7. conduct analysis to see performance of the alignment scoring metrics       ######
###########################################################################################

setwd("/Volumes/Paddy_Backup/ProjectBoard_Patrick/02-Proof_of_principle/scripts/")
source("../lib/ReadsGenerator.R")
source("../lib/DeBruijnGraph.R")
source("../lib/BreakageProbScoring.R")
source("../lib/AlignmentScoring.R")
source("../lib/Assembler.R")

metrics.analysis <- function(
  len = len, G.cont = G.cont,
  C.cont = C.cont, A.cont = A.cont,
  k = k, NCPU = NCPU,
  max.runs = max.runs, # number of iterations
  timer = timer, # max.time per DBG assembly
  length.of.read = length.of.read, # length of each read
  multiplier = multiplier, # read generation number multiplier
  prob.type = prob.type, # "uniform" or "non-uniform"
  analysis = analysis # "r-squared" or "breakpoints" or "breakpoints.performance"
  ){
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
  suppressMessages(packages(c("foreach", "doSNOW")))
  
  # set-up cluster for parallel computation
  NCPU <- NCPU
  cl   <- makeSOCKcluster(NCPU)
  registerDoSNOW(cl)
  
  # setup text progress bar for use in foreach loop
  # supported by the doSNOW package
  pb       <- txtProgressBar(min = 1, max = max.runs, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts     <- list(progress = progress)
  
  if(analysis=="r-squared"){
    print("Running genome assembler and analysing key metrics...", quote = F)
    out.put <- foreach(i=1:max.runs, .combine = "rbind",
                       .export = ls(globalenv()),
                       .inorder = FALSE, .options.snow = opts)%dopar%{
    results <- string.reconstruction(len = len, G.cont = G.cont,
                                    C.cont = C.cont, A.cont = A.cont,
                                    k = k, NCPU = NCPU,
                                    timer = timer,
                                    length.of.read = length.of.read,
                                    multiplier = multiplier,
                                    prob.type = prob.type)
    
    possibleError <- tryCatch(
     suppressWarnings(summary(lm(results$alignment.res$global ~
                                   results$alignment.res$local+
                                   results$alignment.res$kmer))$r.squared),
     error = function(e){e})
    
    if(identical(length(results$alignment.scores), 
                 length(results$prob.ref.table)) & 
      !inherits(possibleError, "error")){
     # obtain correlation coefficients
     # suppress warnings in cases where number of results are very low and column values
     # are identical; would yield a perfect fit in the linear regression calculation;
     # remove all such cases at the very end
     rsq.kmer <-
       suppressWarnings(summary(lm(results$alignment.res$global ~
                                     results$alignment.res$local+
                                     results$alignment.res$kmer))$r.squared)
     rsq.levenshtein <-
       suppressWarnings(summary(lm(results$alignment.res$global ~
                                     results$alignment.res$local+
                                     results$alignment.res$levenshtein))$r.squared)
     rsq.local.global <-
       suppressWarnings(summary(lm(results$alignment.res$global ~
                                     results$alignment.res$local+
                                     results$alignment.res$local.global))$r.squared)
     rsq.global.vs.local <-
       suppressWarnings(summary(lm(results$alignment.res$global ~
                                     results$alignment.res$local))$r.squared)
     
     # Check if highest scores match between alignment and breakage probabilities
     scores.df  <- data.frame("alignment.scores" = results$alignment.scores,
                              "prob.ref.sums"    = results$prob.ref.table)
     counter = 0
     if(identical(
       scores.df[match(max(scores.df$alignment.scores),
                       scores.df$alignment.scores), "prob.ref.sums"],
       max(scores.df$prob.ref.sums)
     )==TRUE){
       counter = 1
     }
    } else {
     # if DBG is too large and assembly takes too long, we stop the run at
     # timer minutes; and save as NA to be omitted later
     rsq.kmer            = NA
     rsq.levenshtein     = NA
     rsq.local.global    = NA
     rsq.global.vs.local = NA
     counter             = NA
    }
    df <- data.frame(r.squared.kmer            = rsq.kmer,
                     r.squared.levenshtein     = rsq.levenshtein,
                     r.squared.local.global    = rsq.local.global,
                     r.squared.global.vs.local = rsq.global.vs.local,
                     goodness.of.predict       = counter,
                     sequence.length           = nchar(results$ref.seq),
                     kmer                      = k,
                     coverage                  = results$coverage)
    return(df)
    }
    print("Metrics analysis complete!", quote = F)
    close(pb)
    stopCluster(cl)
    out.put <- na.omit(out.put)
    return(out.put)
  }
  
  # ref.seq <- vector(mode= "character", length = max.runs)
  if(analysis=="breakpoints"){
    print("Running genome assembler and analysing key metrics...", quote = F)
    out.put <- foreach(i=1:max.runs, .combine = "rbind",
                       .export = ls(globalenv()),
                       .inorder = FALSE, .options.snow = opts)%dopar%{
     results <- string.reconstruction(len = len, G.cont = G.cont,
                                      C.cont = C.cont, A.cont = A.cont,
                                      k = k, NCPU = NCPU,
                                      timer = timer,
                                      length.of.read = length.of.read,
                                      multiplier = multiplier,
                                      prob.type = prob.type)
     
     possibleError <- tryCatch(
       suppressWarnings(summary(lm(results$alignment.res$global ~
                                     results$alignment.res$local+
                                     results$alignment.res$kmer))$r.squared),
       error = function(e){e})
     
     if(identical(length(results$alignment.scores), 
                  length(results$prob.ref.table)) & 
        !inherits(possibleError, "error")){
       
       # df <- data.frame(prob.table        = scale(results$prob.ref.table),
       #                  total.align.score = scale(results$alignment.res$score),
       #                  local             = scale(results$alignment.res$local),
       #                  global            = scale(results$alignment.res$global),
       #                  global.local      = scale(results$alignment.res$global.local),
       #                  local.global      = scale(results$alignment.res$local.global),
       #                  kmer              = scale(results$alignment.res$kmer),
       #                  levenshtein       = scale(results$alignment.res$levenshtein),
       #                  all.sequences     = results$alignment.res$seq,
       #                  ref.seq           = rep(results$ref.seq, length(results$prob.ref.table))
       #                  )
       
       df <- data.frame(prob.table        = results$prob.ref.table,
                        total.align.score = results$alignment.res$score,
                        local             = results$alignment.res$local,
                        global            = results$alignment.res$global,
                        global.local      = results$alignment.res$global.local,
                        local.global      = results$alignment.res$local.global,
                        kmer              = results$alignment.res$kmer,
                        levenshtein       = results$alignment.res$levenshtein,
                        all.sequences     = results$alignment.res$seq,
                        ref.seq           = rep(results$ref.seq, length(results$prob.ref.table))
       )
       return(df)
     }
    }
    print("Metrics analysis complete!", quote = F)
    close(pb)
    stopCluster(cl)
    out.put <- na.omit(out.put)
    # data frame of correlation coefficients between the breakage probability scores
    # and all other scoring metrics 
    rsq.df <- data.frame(rsq.local        = summary(lm(out.put$prob.table ~
                                                         out.put$local))$r.squared,
                         rsq.global       = summary(lm(out.put$prob.table ~
                                                         out.put$global))$r.squared,
                         rsq.global.local = summary(lm(out.put$prob.table ~
                                                         out.put$global.local))$r.squared,
                         rsq.local.global = summary(lm(out.put$prob.table ~
                                                         out.put$local.global))$r.squared,
                         rsq.kmer         = summary(lm(out.put$prob.table ~
                                                         out.put$kmer))$r.squared,
                         rsq.levenshtein  = summary(lm(out.put$prob.table ~
                                                        out.put$levenshtein))$r.squared)
    return(list(out.put, rsq.df))
    return(out.put)
  }
  
  if(analysis=="breakpoints.performance"){
    print("Running genome assembler and calculating performance of breakpoint scores...", quote = F)
    # Helper function to combine results and return list of vectors from foreach
    comb <- function(x, ...) {
      lapply(seq_along(x),
             function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
    }
    
    out.put <- foreach(i = 1:max.runs, 
                       # .combine = "rbind",
                       .export = ls(globalenv()),
                       .combine = 'comb', .multicombine = TRUE,
                       .init = list(list(), list()),
                       .inorder = FALSE, .options.snow = opts)%dopar%{
      results <- string.reconstruction(len = len, G.cont = G.cont,
                                       C.cont = C.cont, A.cont = A.cont,
                                       k = k, NCPU = NCPU,
                                       timer = timer,
                                       length.of.read = length.of.read,
                                       multiplier = multiplier,
                                       prob.type = prob.type)
      
      df.norm <- data.frame(
       local             = results$alignment.res$local,
       global            = results$alignment.res$global,
       global.local      = results$alignment.res$global.local,
       local.global      = results$alignment.res$local.global,
       kmer              = results$alignment.res$kmer,
       levenshtein       = results$alignment.res$levenshtein)
    
      # check if lengths of alignment scores are the same
      if(identical(length(results$alignment.scores), 
                  length(results$prob.ref.table))){
       # check index of the highest score in each obtained alignment scores
       max.ind.probs <- match(max(results$prob.ref.table), results$prob.ref.table)
       max.ind.align <- match(max(results$alignment.res$score), results$alignment.res$score)
       ind.alignment <- c(max.ind.probs, max.ind.align)
       
       # difference scores
       diff.df <- data.frame(
         local        = df.norm$local[ind.alignment[1]]-df.norm$local[ind.alignment[2]],
         global       = df.norm$global[ind.alignment[1]]-df.norm$global[ind.alignment[2]],
         global.local = df.norm$global.local[ind.alignment[1]]-df.norm$global.local[ind.alignment[2]],
         local.global = df.norm$local.global[ind.alignment[1]]-df.norm$local.global[ind.alignment[2]],
         kmer         = df.norm$kmer[ind.alignment[1]]-df.norm$kmer[ind.alignment[2]],
         levenshtein  = df.norm$levenshtein[ind.alignment[1]]-df.norm$levenshtein[ind.alignment[2]]
       )
       
       df.levenshtein <- data.frame(
         levenshtein = df.norm$levenshtein[ind.alignment[1]])
       
       return(list(diff.df, df.levenshtein))
      }
    }
    print("Breakpoint scores performances calculated!", quote = F)
    close(pb)
    stopCluster(cl)
    
    diff.df        <- na.omit(do.call("rbind", out.put[[1]]))
    df.levenshtein <- na.omit(do.call("rbind", out.put[[2]]))
    df.levenshtein <- data.frame(
      levenshtein = mean(df.levenshtein$levenshtein),
      read.length = read.length)
    
    # scoring matrix
    mat <- matrix(NA, nrow = 2, ncol = 6)
    for(col in 1:dim(mat)[2]){
      mat[1, col] <- ceiling(length(which(diff.df[, col]>=0))/dim(diff.df)[1]*100)
      mat[2, col] <- 100-(ceiling(length(which(diff.df[, col]>=0))/dim(diff.df)[1]*100))
      if(col==6){
        mat[1, col] <- ceiling(length(which(diff.df[, col]<=0))/dim(diff.df)[1]*100)
        mat[2, col] <- 100-(ceiling(length(which(diff.df[, col]<=0))/dim(diff.df)[1]*100))
      }
    }
    colnames(mat) <- c("local", "global", "global-local", "local-global", "kmer", "levenshtein")
    rownames(mat) <- c("breakpoints.scores", "alignment.scores")
    
    return(list(mat, df.levenshtein))
  }
}

NCPU        = 4 
max.runs    = 20
length      = 1000
kmer        = 9
read.length = 15
multiplier  = 10000

# Run the assembly algorithm for max.runs number of iterations and conduct
# analysis on alignment scores metrics
results <- metrics.analysis(len = length, G.cont = 0.25, C.cont = 0.25,
                            A.cont = 0.25, timer = 1,
                            k = kmer, NCPU = NCPU,
                            max.runs = max.runs,
                            length.of.read = read.length,
                            multiplier = multiplier,
                            prob.type = "non-uniform",
                            analysis = "breakpoints.performance") 

# Run the assembly algorithm 
results <- string.reconstruction(len = length, G.cont = 0.25, C.cont = 0.25,
                                 A.cont = 0.25, timer = 1,
                                 k = kmer, NCPU = NCPU,
                                 length.of.read = read.length,
                                 multiplier = multiplier,
                                 prob.type = "non-uniform")

##################################################################
# write.csv(results, file = "../data/correlation_coefficient_analysis.csv",
#           row.names = FALSE)
##################################################################
