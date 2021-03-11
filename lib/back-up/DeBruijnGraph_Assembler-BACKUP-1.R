# This file assembles the reads together using the 
# Paired De Bruijn Graph Assembly algorithm approach

# string.reconstruction <- function(k){
#   return(genome.path.construction(
#     eulerian.path(
#       get.balance.count(
#         debrujin.graph.from.kmers(
#           kmer.composition(k)
#           )
#         )
#       )
#     ))
# }

kmer.composition <- function(k){
  reads <- scan("../data/reads.txt", what = character(), sep = "", quiet = TRUE)
  kmers <- lapply(1:length(reads), function(i){
    sequence <- unlist(strsplit(reads[i], split = ""))
    end      <- length(sequence)-k+1
    kmer     <- substring(paste(sequence, collapse=""),
                          first=1:end, last=(1:end)+k-1)
  })
  kmers   <- unlist(kmers)
  kmer.df <- as.data.frame(table(kmers))
  kmers   <- unlist(lapply(kmer.df$kmers, as.character))
  newList <- list(kmers, k)
  return(newList)
}

debrujin.graph.from.kmers <- function(newList){
  # re-define objects from previous functions' output
  patterns <- newList[[1]]
  k        <- newList[[2]]
  
  prefix <- sapply(1:length(patterns), function(i){
    substring(patterns[i], first = 1, last = k-1)
  })
  
  suffix <- sapply(1:length(patterns), function(i){
    substring(patterns[i], first = 2, last = k)
  })
  
  # Initialise dictionary
  dict <- vector(mode="list", length=length(unique(prefix)))
  
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
  print("De Brjuin graph obtained from kmers!", quote = FALSE)
  newList <- list(dict, prefix, suffix)
  return(newList)
}

get.balance.count <- function(newList){
  # re-define objects from previous functions' output
  dict   <- newList[[1]]
  prefix <- newList[[2]]
  suffix <- newList[[3]]
  
  # Initialise dictionary
  balanced.count <- vector(mode="list", length=length(dict))
  
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
  print("Eulerian path is balanced!", quote = FALSE)
  newList <- list(newList[[1]], suffix, balanced.count)
  return(newList)
}

eulerian.path <- function(newList){
  # re-define objects from previous functions' output
  dict           <- newList[[1]]
  suffix         <- newList[[2]]
  balanced.count <- newList[[3]]
  # create copies for use in the pathway explorations
  dict.copy           <- dict
  suffix.copy         <- suffix
  balanced.count.copy <- balanced.count
  
  # obtain all possible starting points of the genome
  start.point <- list()
  for(i in 1:length(balanced.count)){
    if(balanced.count[[i]]==-1){
      start.point <- c(start.point, names(balanced.count)[i])
    }
  }
  
  # obtain all possible branching points and number of options per branching point
  tot.ind <- unlist(sapply(1:length(dict), function(i){
    if(length(dict[[i]])>1){
      length(dict[[i]])
    }
  }))
  
  # initialise list of pathway combinations
  pathways <- vector(mode = "list", length(start.point)); pathways
  
  # loop over combinations
  combs <- expand.grid(lapply(1:length(tot.ind), function(i){1:tot.ind[i]})); combs
  
  # initialise empty lists
  combs.paths <- vector(mode="list", length=dim(combs)[1]); combs.paths
  
  # loop over all possible pathway combinations
  len.start.points <- 1
  while(len.start.points<=length(pathways)){
    rows <- 1
    while(rows<=dim(combs)[1]){
      dict           <- dict.copy
      suffix         <- suffix.copy
      balanced.count <- balanced.count.copy
      
      path <- start.point[[len.start.points]]
      
      ### Loop over individual columns per row ### 
      i=1
      while(i<=dim(combs)[2]+1){
        last.item <- tail(path, n=1)
        
        if(length(dict[last.item][[1]])>1){
          new.item  <- dict[last.item][[1]][combs[rows,i]]
          i = i+1
        } else {
          new.item  <- dict[last.item][[1]][1]
        }
        
        if(length(new.item)>0){
          path    <- c(path, new.item)
        
          if(length(dict[last.item][[1]])>1){
            dict[last.item][[1]] <- dict[last.item][[1]][-which(new.item==dict[last.item][[1]])]
          } else {
            dict[last.item] <- NULL
          }
        } else {
          combs.paths[[rows]] <- path
          rows = rows+1
          break
        }
      }
    }
    pathways[[len.start.points]] <- c(pathways[[len.start.points]], combs.paths)
    len.start.points = len.start.points+1
  }
  print("Eulerian path is constructed!", quote = FALSE)
  return(pathways)
}

genome.path.construction <- function(path, randSeq){
  source("Pairwise_Alignment.R")
  
  combs.paths <- path
  # initialise vector
  all.genomes <- vector(mode = "list", length = length(combs.paths))
  # obtain all reconstructed genomes and save as characters
  for(combs.paths.len in 1:length(combs.paths)){
    for(i in 1:length(combs.paths[[combs.paths.len]])){
      genome <- character()
      len    <- nchar(combs.paths[[combs.paths.len]][[i]][1])
      genome <- c(substring(combs.paths[[combs.paths.len]][[i]][1], first=1:len, last=1:len), 
                  substring(combs.paths[[combs.paths.len]][[i]][2:length(combs.paths[[combs.paths.len]][[i]])], 
                            first=len, last=len))
      all.genomes[[combs.paths.len]] <- c(all.genomes[[combs.paths.len]], list(genome))
    }
  }
  
  # compare length of reference genome and generated genomes
  compare.len <- sapply(1:length(combs.paths), function(i){
    length(randSeq)-sapply(all.genomes[[i]], length)
  })
  
  # return row index of the genomes with the smallest differences 
  # to the reference genome 
  smallest.diff <- apply(compare.len, 2, which.min); smallest.diff
  
  # obtain optimal alignment scores
  scores <- sapply(1:length(smallest.diff), function(i){
    aligner(seq1 = paste(all.genomes[[i]][[smallest.diff[i]]], collapse = ""),
            seq2 = paste(randSeq, collapse = ""), score = TRUE)
  })
  
  # return re-contructed genome path of best match
  total.matches <- max(scores)
  best.match    <- all.genomes[[which.max(scores)]][[smallest.diff[which.max(scores)]]]
  best.match    <- paste(best.match, collapse = "")
  
  print("Full genome path is reconstructed!", quote = FALSE)
  return(list(best.match, paste0(total.matches,"% similarity")))
}

randSeq <- ReadSizeGenerator(len = 100, G.cont = 0.3, C.cont = 0.3, A.cont = 0.2)
# paste(randSeq, collapse = "")
genome.path.construction(eulerian.path(get.balance.count(debrujin.graph.from.kmers(kmer.composition(6)))), randSeq)

