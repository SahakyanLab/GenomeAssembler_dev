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
suppressMessages(packages(c("foreach", "ggExtra", "ggplot2", "cowplot", "Rcpp",
                            "stringr", "doRNG", "dplyr", "tidyr", "doParallel", "itertools")))

setwd("/Volumes/Paddy_Backup/ProjectBoard_Patrick/02-Proof_of_principle/scripts/")
source("../lib/ReadsGenerator.R")
source("../lib/DeBruijnGraph.R")
source("../lib/BreakageProbScoring.R")
source("../lib/AlignmentScoring.R")
source("../lib/Assembler.R")
# sourceCpp(file = '../lib/edlibFunction.cpp')

##########################################################################################
pearson.correlation <- function(
  len = len, G.cont = G.cont,
  C.cont = C.cont, A.cont = A.cont,
  k = k, NCPU = NCPU,
  max.runs = max.runs, # number of iterations
  timer = timer, # max.time per DBG assembly
  length.of.read = length.of.read, # length of each read
  multiplier = multiplier, # read generation number multiplier
  prob.type = prob.type, # "uniform" or "non-uniform"
  seed = seed
){
  # set-up cluster for parallel computation
  cl <- makeCluster(NCPU)
  registerDoParallel(cl)
  
  # Declare that parallel RNG should be used for in a parallel foreach() call.
  # %dorng% will still result in parallel processing; uses %dopar% internally.
  registerDoRNG(seed = seed)
  
  comb <- function(x, ...) {
    lapply(seq_along(x),
           function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
  }
  
  print("Running genome assembler and pearson correlation coefficients...", quote = F)
  out.put <- foreach(i = 1:max.runs, 
                     .combine = 'comb', .multicombine = TRUE,
                     .init = list(list(), list(), list(), list()),
                     .export = ls(globalenv()),
                     .inorder = TRUE)%dopar%{
   results <- string.reconstruction(len = len, G.cont = G.cont,
                                    C.cont = C.cont, A.cont = A.cont,
                                    k = k, NCPU = NCPU,
                                    timer = timer,
                                    length.of.read = length.of.read,
                                    multiplier = multiplier,
                                    prob.type = prob.type,
                                    seed = seed)
   
   # check if lengths of alignment scores are the same
   if(identical(length(results$alignment.scores),
                length(results$prob.ref.table))){
     df        <- cbind(results$prob.ref.table, results$alignment.res[1:6])
     
     # data frame of the breakage probability scores and the length of each 
     # assembled genome sequence
     seq.dist.df <- data.frame(prob.table        = results$prob.ref.table,
                               assembled.seq     = results$prob.ref.seq,
                               total.align.score = results$alignment.res$score)
     
     # check index of the highest score in each obtained alignment scores
     max.ind.probs <- match(max(results$prob.ref.table), results$prob.ref.table)
     max.ind.align <- match(max(results$alignment.res$score), results$alignment.res$score)
     ind.alignment <- c(max.ind.probs, max.ind.align)
     
     # difference scores
     diff.df <- data.frame(
       local        = results$alignment.res$local[ind.alignment[1]]-results$alignment.res$local[ind.alignment[2]],
       global       = results$alignment.res$global[ind.alignment[1]]-results$alignment.res$global[ind.alignment[2]],
       global.local = results$alignment.res$global.local[ind.alignment[1]]-results$alignment.res$global.local[ind.alignment[2]],
       local.global = results$alignment.res$local.global[ind.alignment[1]]-results$alignment.res$local.global[ind.alignment[2]],
       kmer         = results$alignment.res$kmer[ind.alignment[1]]-results$alignment.res$kmer[ind.alignment[2]])
     
     # read lengths distribution
     df.reads <- data.frame(reads = results$reads.freq)
   }
   return(list(df, seq.dist.df, diff.df, df.reads))
  }
  stopCluster(cl)

  rsq.df      <- do.call("rbind", out.put[[1]])
  seq.dist.df <- do.call("rbind", out.put[[2]])
  diff.df     <- do.call("rbind", out.put[[3]])
  df.reads    <- do.call("rbind", out.put[[4]])
  
  rsq.df      <- na.omit(rsq.df)
  seq.dist.df <- na.omit(seq.dist.df)
  diff.df     <- na.omit(diff.df)
  df.reads    <- na.omit(df.reads)
  
  # data frame of correlation coefficients between the breakage probability scores
  # and all other scoring metrics 
  colnames(rsq.df) <- c("probs", "alignment.score", "local", "global", 
                        "global.local", "local.global", "kmer")
  rsq.df.sum <- suppressWarnings(as.data.frame(t(cor(rsq.df)[1,]))[3:7])
  
  # scoring matrix for percent of correctness 
  mat <- matrix(NA, nrow = 2, ncol = 5)
  for(col in 1:dim(mat)[2]){
    mat[1, col] <- ceiling(length(which(diff.df[, col]>=0))/dim(diff.df)[1]*100)
    mat[2, col] <- 100-(ceiling(length(which(diff.df[, col]>=0))/dim(diff.df)[1]*100))
    if(col==6){
      mat[1, col] <- ceiling(length(which(diff.df[, col]<=0))/dim(diff.df)[1]*100)
      mat[2, col] <- 100-(ceiling(length(which(diff.df[, col]<=0))/dim(diff.df)[1]*100))
    }
  }
  rownames(mat) <- c("breakpoints.scores", "alignment.scores")
  # convert to stacked ggplot friendly format
  dat <- as.data.frame(t(mat))
  
  datm <- dat %>% 
    mutate(ind = factor(row_number())) %>%  
    gather(variable, value, -ind)
  
  datm$ind <- rep(c("local", "global", "global-local", "local-global", "kmer"), 2)
  
  return(list(rsq.df, rsq.df.sum, seq.dist.df, diff.df, datm, df.reads))
}

##########################################################################################

multi.regression <- function(
  len = len, G.cont = G.cont,
  C.cont = C.cont, A.cont = A.cont,
  k = k, NCPU = NCPU,
  max.runs = max.runs, # number of iterations
  timer = timer, # max.time per DBG assembly
  length.of.read = length.of.read, # length of each read
  multiplier = multiplier, # read generation number multiplier
  prob.type = prob.type, # "uniform" or "non-uniform"
  seed = seed
){
  # set-up cluster for parallel computation
  cl <- makeCluster(NCPU)
  registerDoParallel(cl)
  
  # Declare that parallel RNG should be used for in a parallel foreach() call.
  # %dorng% will still result in parallel processing; uses %dopar% internally.
  registerDoRNG(seed = seed)
  
  print("Running genome assembler and analysing key metrics...", quote = F)
  out.put <- foreach(i = 1:max.runs, 
                     .combine = 'rbind', .export = ls(globalenv()),
                     .inorder = FALSE)%dopar%{
   results <- string.reconstruction(len = len, G.cont = G.cont,
                                    C.cont = C.cont, A.cont = A.cont,
                                    k = k, NCPU = NCPU,
                                    timer = timer,
                                    length.of.read = length.of.read,
                                    multiplier = multiplier,
                                    prob.type = prob.type,
                                    seed = seed)
    
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
     rsq.local.global    = NA
     rsq.global.vs.local = NA
     counter             = NA
    }
    df <- data.frame(r.squared.kmer           = rsq.kmer,
                    r.squared.local.global    = rsq.local.global,
                    r.squared.global.vs.local = rsq.global.vs.local,
                    goodness.of.predict       = counter,
                    sequence.length           = nchar(results$ref.seq),
                    kmer                      = k,
                    coverage                  = results$coverage)
    return(df)
  }
  print("Metrics analysis complete!", quote = F)
  stopCluster(cl)
  out.put <- na.omit(out.put)
  return(out.put)
}
##########################################################################################
# import human genome atgc composition as data frame
human.atgc.df <- read.csv(file = "../data/HumanGenomeATGC.csv")

# define parameters for assembly functions
NCPU        = 4 
max.runs    = 5
length      = 1000
kmer        = 9
readLengthSeq <- c(15, 30, 60, 90, 100, 120)
multiplierSeq <- c(15000, 1800, 60, 90, 53, 32)

# readLengthSeq <- c(60)
# multiplierSeq <- c(80)

for(i in 1:length(readLengthSeq)){
  results <- pearson.correlation(
    len = length, 
    G.cont = human.atgc.df$G.cont, 
    C.cont = human.atgc.df$C.cont,
    A.cont = human.atgc.df$A.cont, 
    timer = 1,
    k = kmer, NCPU = NCPU,
    max.runs = max.runs,
    length.of.read = readLengthSeq[i],
    multiplier = multiplierSeq[i],
    prob.type = "non-uniform", 
    seed = 1935)

  if(i==1){
    cor.df <- data.frame(results[[2]])
  } else {
    cor.df <- rbind(cor.df, results[[2]])
  }
  
  if(i==length(readLengthSeq)){
    write.csv(cor.df, file = "../data/correlation-coefficient-df.csv",
              row.names = FALSE)
  }
  # Plots of assembled sequence length distribution vs. breakage prob scores
  sorted.len  <- sort(results[[3]]$assembled.seq, decreasing = TRUE)
  sorted.len  <- unique(sorted.len)[1]
  top.results <- as.numeric(rownames(results[[3]][which(results[[3]]$assembled.seq==sorted.len),]))
  
  p.seq <- ggplot(results[[3]], 
              aes(x = prob.table, y = assembled.seq)) + geom_point() +
    geom_point(data=results[[3]][top.results,], 
               aes(x = prob.table, y = assembled.seq), colour="red") +
    labs(x = "Breakage Probability Scores", y = "Assembled Sequence Length")
  
  p.seq <- ggMarginal(p.seq, margins = "y", type="density",
                      colour = "black", fill = "black", alpha = .3)
  
  p.align.probs <- ggplot(results[[3]], 
                        aes(x = prob.table, y = total.align.score)) + geom_point() + 
    geom_point(data=results[[3]][top.results,], 
               aes(x = prob.table, y = total.align.score), colour="red") +
    labs(x = "Breakage Probability Scores", y = "Reference Alignment Scores")
  
  combined.plots <- plot_grid(p.seq, NULL, p.align.probs, 
                              nrow = 1,
                              rel_widths = c(1, 0.05, 1),
                              labels = c("A", "", "B"))
  
  title <- ggdraw() + draw_label(
    paste0("Seq Length: ",length,"; ",
           "Kmer: ",kmer,"; ",
           "Read Length: ",readLengthSeq[i]), fontface='bold')
  
  # rel_heights values control title margins
  combined.plots <- plot_grid(title, combined.plots, ncol=1, rel_heights=c(0.1, 1))
  
  ggsave(combined.plots, width=10, height=6,
         file=paste0("../figures/BP-Score_vs_Seq-Length/",
                     "Len_",length, "-",
                     "Kmer",kmer,"-",
                     "ReadLength",readLengthSeq[i],".pdf"))
  
  # Plots of correlation coefficients of breakage probability scores vs.
  # alignment scoring metrics
  p1 <- ggplot(results[[1]], aes(x = probs, y = local)) + geom_point() +
    geom_text(data = data.frame(xpos = -Inf, ypos =  Inf,
                                annotateText = paste0("r = ", round(results[[2]]$local,4)),
                                hjustvar = 0, vjustvar = 1), 
              aes(x = xpos, y = ypos, hjust = hjustvar, vjust = vjustvar, label = annotateText),
              fontface = "bold") + geom_smooth(method='lm', formula= y~x) + 
    labs(x = "Breakage Probability", y = paste(str_to_title(names(results[[2]])[1]), "Alignment"))
  
  p2 <- ggplot(results[[1]], aes(x = probs, y = global)) + geom_point() +
    geom_text(data = data.frame(xpos = -Inf, ypos =  Inf,
                                annotateText = paste0("r = ", round(results[[2]]$global,4)),
                                hjustvar = 0, vjustvar = 1), 
              aes(x = xpos, y = ypos, hjust = hjustvar, vjust = vjustvar, label = annotateText),
              fontface = "bold") + geom_smooth(method='lm', formula= y~x) + 
    labs(x = "Breakage Probability", y = paste(str_to_title(names(results[[2]])[2]), "Alignment"))
  
  p3 <- ggplot(results[[1]], aes(x = probs, y = global.local)) + geom_point() +
    geom_text(data = data.frame(xpos = -Inf, ypos =  Inf,
                                annotateText = paste0("r = ", round(results[[2]]$global.local,4)),
                                hjustvar = 0, vjustvar = 1), 
              aes(x = xpos, y = ypos, hjust = hjustvar, vjust = vjustvar, label = annotateText),
              fontface = "bold") + geom_smooth(method='lm', formula= y~x) + 
    labs(x = "Breakage Probability",
         y = paste(str_to_title(gsub("[.]", "-", names(results[[2]])[3])), "Alignment"))
  
  p4 <- ggplot(results[[1]], aes(x = probs, y = local.global)) + geom_point() +
    geom_text(data = data.frame(xpos = -Inf, ypos =  Inf,
                                annotateText = paste0("r = ", round(results[[2]]$local.global,4)),
                                hjustvar = 0, vjustvar = 1), 
              aes(x = xpos, y = ypos, hjust = hjustvar, vjust = vjustvar, label = annotateText),
              fontface = "bold") + geom_smooth(method='lm', formula= y~x) + 
    labs(x = "Breakage Probability", 
         y = paste(str_to_title(gsub("[.]", "-", names(results[[2]])[4])), "Alignment"))
  
  p5 <- ggplot(results[[1]], aes(x = probs, y = kmer)) + geom_point() +
    geom_text(data = data.frame(xpos = -Inf, ypos =  Inf,
                                annotateText = paste0("r = ", round(results[[2]]$kmer,4)),
                                hjustvar = 0, vjustvar = 1), 
              aes(x = xpos, y = ypos, hjust = hjustvar, vjust = vjustvar, label = annotateText),
              fontface = "bold") + geom_smooth(method='lm', formula= y~x) + 
    labs(x = "Breakage Probability", y = paste(str_to_title(names(results[[2]])[5]), "Alignment"))
  
  # p6 <- ggplot(results[[1]], aes(x = probs, y = levenshtein)) + geom_point() +
  #   geom_text(data = data.frame(xpos = -Inf, ypos =  Inf,
  #                               annotateText = paste0("r = ", round(results[[2]]$levenshtein,4)),
  #                               hjustvar = 0, vjustvar = 1),
  #             aes(x = xpos, y = ypos, hjust = hjustvar, vjust = vjustvar, label = annotateText),
  #             fontface = "bold") + geom_smooth(method='lm', formula= y~x) +
  #   labs(x = "Breakage Probability", y = paste(str_to_title(names(results[[2]])[6]), "Alignment"))
  
  title <- ggdraw() + draw_label(
    paste0("Correlation Coefficients of Breakage Probability Scores vs. Alignment Scores",
           "\n",
           "Seq Length: ",length,"; ",
           "Kmer: ",kmer,"; ",
           "Read Length: ",readLengthSeq[i]), fontface='bold')
  
  combined.cor.plots <- plot_grid(p1, p2, p5,
                                  p4, p3,
                                  nrow = 2)
  
  # rel_heights values control title margins
  combined.cor.plots <- plot_grid(title, combined.cor.plots, ncol = 1, rel_heights = c(0.1, 1))
  
  ggsave(combined.cor.plots, width=10, height=6,
         file=paste0("../figures/Correlations/",
                     "Len_",length, "-",
                     "Kmer",kmer,"-",
                     "ReadLength",readLengthSeq[i],".pdf"))
  
  # Stacked barplot showing the percent of correctness of breakpoints vs. 
  # each alignment metric; we define "correctness" as a breakpoint score that is
  # equal to or greater than the alignment score
  percent.correctness <- ggplot(results[[5]], aes(x = ind, y = value, fill = variable)) + 
    geom_bar(position="stack", stat="identity") + 
    labs(x = "Alignment Scores Metrics", y = "Percent of correctness", 
         title = "Percent of correctness between breakpoints and alignment scores")
  
  ggsave(percent.correctness, width=10, height=6,
         file=paste0("../figures/Alignment_scores/",
                     "Len_",length, "-",
                     "Kmer",kmer,"-",
                     "ReadLength",readLengthSeq[i],".pdf"))
  
  # density plot of read lengths
  p <- ggplot(results[[6]], aes(x=reads)) + 
    geom_density() + lims(x = c(min(results[[6]])*0.75,max(results[[6]])*1.25)) + 
    labs(x = "Read Length", y = "Density", 
         title = "Distribution of read lengths") 
  
  ggsave(p, width=10, height=6,
         file=paste0("../figures/ReadLength/ReadLength-dist-",
                     "Len_",length, "-",
                     "Kmer",kmer,"-",
                     "ReadLength",readLengthSeq[i],".pdf"))
}
# 
# readLengthSeq <- c(15, 30, 60, 90, 100, 120)
# multiplierSeq <- c(15000, 1800, 60, 90, 53, 32)
# i = 5
# results <- string.reconstruction(
#   len = length, 
#   G.cont = human.atgc.df$G.cont, 
#   C.cont = human.atgc.df$C.cont,
#   A.cont = human.atgc.df$A.cont, 
#   timer = 1,
#   k = kmer, NCPU = NCPU,
#   length.of.read = readLengthSeq[i],
#   multiplier = multiplierSeq[i],
#   prob.type = "non-uniform", 
#   seed = 9264)
# 
# # Plots of correlation coefficients of breakage probability scores vs.
# # alignment scoring metrics
# df.probs <- data.frame(prob.table = results$prob.ref.table, 
#                        total.align.score = results$alignment.scores)
# 
# p.align.probs <- ggplot(df.probs, 
#                         aes(x = prob.table, y = total.align.score)) + geom_point() 
# p.align.probs
# 
# df <- cbind(results$prob.ref.table, results$alignment.res)
# colnames(df) <- c("probs", "alignment.score", "local", "global", 
#                   "global.local", "local.global", "kmer")
# p1 <- ggplot(df, aes(x = probs, y = local)) + geom_point()
# p2 <- ggplot(df, aes(x = probs, y = global)) + geom_point()
# p3 <- ggplot(df, aes(x = probs, y = global.local)) + geom_point()
# p4 <- ggplot(df, aes(x = probs, y = local.global)) + geom_point()
# p5 <- ggplot(df, aes(x = probs, y = kmer)) + geom_point()
#  
# title <- ggdraw() + draw_label(
#   paste0("Correlation Coefficients of Breakage Probability Scores vs. Alignment Scores",
#          "\n",
#          "Seq Length: ",length,"; ",
#          "Kmer: ",kmer,"; ",
#          "Read Length: ",readLengthSeq[i]), fontface='bold')
# 
# combined.cor.plots <- plot_grid(p1, p2, p5,
#                                 p4, p3,
#                                 nrow = 2)
# 
# combined.cor.plots
# 9264
# 5392
# 9999
# 111000
# results <- multi.regression(
#   len = length, G.cont = 0.25, C.cont = 0.25,
#   A.cont = 0.25, timer = 1,
#   k = kmer, NCPU = NCPU,
#   max.runs = max.runs,
#   length.of.read = read.length,
#   multiplier = multiplier,
#   prob.type = "non-uniform")
# 