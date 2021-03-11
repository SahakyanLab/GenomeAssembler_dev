# This file assembles the reads together using the 
# De Bruijn Graph Assembly algorithm approach

kmer.composition <- function(k){
  # import all generatd reads 
  reads <- scan("../data/reads.txt", what = character(), sep = "", quiet = TRUE)
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
  print("Weighted k-mers obtained!", quote = FALSE)
  return(newList)
}

debrujin.graph.from.kmers <- function(newList){
  # re-define objects from previous functions' output
  patterns <- newList[[1]]
  k        <- newList[[2]]
  # define prefix for each k-mer
  prefix <- sapply(1:length(patterns), function(i){
    substring(patterns[i], first = 1, last = k-1)
  })
  # define suffix for each k-mer
  suffix <- sapply(1:length(patterns), function(i){
    substring(patterns[i], first = 2, last = k)
  })
  
  # Initialise dictionary
  dict <- vector(mode = "list", length = length(unique(prefix)))
  
  # set keys as prefix
  names(dict) <- unique(prefix)
  
  # fill dict with pre- and suffixes
  for(i in 1:length(patterns)){
    if(match(prefix[i], names(dict))==i){
      dict[i] <- suffix[i]
    } else {
      ind <- match(prefix[i], names(dict))
      dict[[ind]] <- c(dict[[ind]], suffix[i])
    }
  }
  print("De Brjuin graph obtained from k-mers!", quote = FALSE)
  newList <- list(dict, prefix, suffix)
  return(newList)
}

get.balance.count <- function(newList){
  # re-define objects from previous functions' output
  dict   <- newList[[1]]
  prefix <- newList[[2]]
  suffix <- newList[[3]]
  
  # Initialise dictionary to count the 
  # number of in-degree and out-degree per vertex 
  # in an Eulerian trail
  balanced.count <- vector(mode = "list", length = length(dict))
  
  # set keys as prefix
  names(balanced.count) <- names(dict)
  
  # set values to zero
  balanced.count <- lapply(balanced.count, function(i) {
    balanced.count[i] <- 0
  })
  
  prefix.table <- as.data.frame(table(prefix))
  suffix.table <- as.data.frame(table(suffix))

  for(i in 1:length(balanced.count)){
    balanced.count[which(as.character(
      prefix.table$prefix)[i] == names(balanced.count))] <- -prefix.table$Freq[i]
  }
  # fill dictionary with the total number of in-degree and out-degree
  # per vertex in the Eulerian trail
  for(i in 1:dim(suffix.table)[1]){
    possibleError <- tryCatch(
      balanced.count[which(as.character(suffix.table$suffix)[i] == names(balanced.count))][[1]],
      error = function(e){e}
      )
    
    if(!inherits(possibleError, "error")){
      val <- balanced.count[which(as.character(
        suffix.table$suffix)[i] == names(balanced.count))][[1]] + suffix.table$Freq[i]
      balanced.count[which(as.character(
        suffix.table$suffix)[i] == names(balanced.count))] <- val 
    } else {
      list1          <- as.list(1)
      names(list1)   <- as.character(suffix.table$suffix)[i]
      balanced.count <- c(balanced.count, list1)
      rm(list1)
    }
  }
  print("Eulerian path balanced!", quote = FALSE)
  newList <- list(newList[[1]], suffix, balanced.count)
  return(newList)
}

genome.construction <- function(path){
  # initialise vector to save genome path
  genome <- character()
  # reconstruct genome path for a given string
  len <- nchar(path[1])
  genome <- c(substring(path[1], first = 1:len, last = 1:len),
              substring(path[2:length(path)], first = len, last = len))
  genome <- paste(genome, collapse = "")
  return(genome)
} 

eulerian.path <- function(newList){
  # re-define objects from previous functions' output
  dict           <- newList[[1]]
  suffix         <- newList[[2]]
  balanced.count <- newList[[3]]
  # create copies for use in the pathway explorations
  dict.copy <- dict
  
  # obtain all possible starting points of the genome
  start.point <- list()
  for(i in 1:length(balanced.count)){
    if(balanced.count[[i]]==-1){
      start.point <- c(start.point, names(balanced.count)[i])
    }
  }
  
  # initialise list of reconstructed genome paths with the highest score
  pathways <- list(); pathways
  # initialise text progress bar
  print("Traversing all paths in the rooted tree...", quote = F)
  pb <- txtProgressBar(min = 1, max = length(start.point)+1, style = 3)
  start.point.index = 1

    while(start.point.index<=length(start.point)){
    # initialise list to store number of branching points per step
    branch.point <- list()
    # initialise list to save the reconstructed genome path
    all.paths <- list()
    # initialise current path with starting point
    path <- start.point[[start.point.index]]
    main.index = 1
    continue = TRUE
    while(continue){
      if(main.index==1){
        ### iteration one ###
        i = 1
        while(TRUE){
          last.item <- tail(path, n = 1)
          all.items <- dict[last.item][[1]]
          # store the number of branching points per step
          if(length(all.items)>0){
            branch.point[[i]] <- seq(1,length(all.items),1) 
          }
          # always traverse the tree along the left-most path
          new.item  <- all.items[1]
          if(length(new.item)>0){
            path <- c(path, new.item)       
            if(length(dict[last.item][[1]])>1){
              dict[last.item][[1]] <- dict[last.item][[1]][-which(new.item==dict[last.item][[1]])]
            } else {
              dict[last.item] <- NULL
              }
            } else {
              # save reconstructed genome path into list
              all.paths <- c(all.paths, genome.construction(path))
              # pop reverse path from the list 
              # until latest branching points first index
              for(l in length(branch.point):1){
                if(length(branch.point[[l]])>1){
                  branch.point[[l]] <- branch.point[[l]][-branch.point[[l]][1]]
                  break
                } else {
                  branch.point[[l]] <- NULL
                }
              }
              main.index = main.index+1
              break
            }
          i = i+1
        }
      } else {
        ### iteration two and onward ###
        while(length(branch.point)>=0){
          # reset dict
          dict <- dict.copy
          # re-initialise current path with starting point
          path <- start.point[[start.point.index]]
          i = 1
          while(TRUE){
            last.item <- tail(path, n = 1)
            all.items <- dict[last.item][[1]]
            if(i<=length(branch.point)){
              # always traverse the tree along the left-most path
              new.item  <- all.items[branch.point[[i]][1]]
            } else {
              new.item  <- all.items[1]
              # store the number of branching points per step
              if(length(all.items)>0){
                branch.point[[i]] <- seq(1,length(all.items),1)
              }
            }
            if(length(new.item)>0){
              path <- c(path, new.item)
              if(length(dict[last.item][[1]])>1){
                dict[last.item][[1]] <- dict[last.item][[1]][-which(new.item==dict[last.item][[1]])]
              } else {
                dict[last.item] <- NULL
                }
              } else {
                # save reconstructed genome path into list
                all.paths <- c(all.paths, genome.construction(path))
                # pop reverse path from the list 
                # until latest branching points first index
                for(l in length(branch.point):1){
                  if(length(branch.point[[l]])>1){
                    branch.point[[l]] <- branch.point[[l]][-branch.point[[l]][1]]
                    break
                  } else {
                    branch.point[[l]] <- NULL
                  }
                }
                break
              }
            i = i+1
          }
          # if no more branches to traverse with current starting point
          # save all generated strings from current starting point
          # to another list; restart everything with
          # different starting point
          if(length(branch.point)==0){
            pathways <- c(pathways, all.paths)
            start.point.index = start.point.index+1
            continue = FALSE
            break
          }
        }
      }
    }
    setTxtProgressBar(pb, start.point.index)
  }
  close(pb)
  print("Eulerian paths constructed!", quote = FALSE)
  return(pathways)
}

alignment.scoring <- function(path, randSeq, NCPU = NCPU){
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
  suppressMessages(packages(c("Biostrings", 
                              "doParallel", "foreach", "itertools")))
  
  # set-up cluster for parallel computation
  cl <- makeCluster(NCPU)
  registerDoParallel(cl)
  
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
  
  # calculate alignment scores of all paths and return only the sequence with max score
  print("Calculating alignment scores...", quote = F)
  alignment.runs <- foreach(chunk = isplitVector(1:length(path), 
                                               chunks = ceiling(length(path)/NCPU)),
                            .packages = c("foreach", "Biostrings"), .combine = "rbind")%dopar%{
    foreach(i = chunk, .combine = "rbind", .inorder = FALSE)%do%{
      seq.new <- DNAString(path[[i]])
      mat     <- nucleotideSubstitutionMatrix(match = 1, 
                                              mismatch = -3, baseOnly = TRUE)
      # global pairwise alignment
      globalAlign  <- 
        pairwiseAlignment(ref.seq, seq.new, 
                          type = "global", substitutionMatrix = mat, 
                          gapOpening = 5, gapExtension = 2, scoreOnly = TRUE)
      # local pairwise alignment
      localAlign   <- 
        pairwiseAlignment(ref.seq, seq.new, 
                          type = "local", substitutionMatrix = mat,
                          gapOpening = 5, gapExtension = 2, scoreOnly = TRUE)
      # global-local pairwise alignment
      global.localAlign <- 
        pairwiseAlignment(ref.seq, seq.new, 
                          type = "global-local", substitutionMatrix = mat,
                          gapOpening = 5, gapExtension = 2, scoreOnly = TRUE)
      # local-global pairwise alignment
      local.globalAlign <- 
        pairwiseAlignment(ref.seq, seq.new, 
                          type = "local-global", substitutionMatrix = mat,
                          gapOpening = 5, gapExtension = 2, scoreOnly = TRUE)
      # higher order k-meric signature counting
      kmers.new  <- kmer.count(as.character(seq.new), k)
      kmer.score <- ceiling(length(match(which(!is.na(match(kmers.ref, kmers.new))==TRUE), 
                                         kmers.new))/length(kmers.ref)*100)
      # levenshtein distance; the lower the levenshtein score, the better;
      # hence need to subtract from final score
      ls.score <- as.numeric(stringDist(c(as.character(ref.seq), as.character(seq.new)), 
                                        method = "levenshtein"))
      finalAlign <- globalAlign+localAlign+global.localAlign+local.globalAlign+kmer.score-ls.score
      out.put <- data.frame(score = finalAlign,
                            seq = as.character(seq.new))
      return(out.put)
    }
  }
  stopCluster(cl)  
  print("Alignment scores calculated!", quote = F)
  return(list(alignment.runs$seq[match(max(alignment.runs$score), alignment.runs$score)], 
              alignment.runs$score))
  
  # return(list(alignment.runs$seq[match(max(alignment.runs$score), alignment.runs$score)],
  #             max(alignment.runs$score)))
}

