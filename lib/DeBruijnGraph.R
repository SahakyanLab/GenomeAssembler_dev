###########################################################################################
###### This file assembles the breakage probability-biased reads together into a     ######
###### fully assembled genome using the De Bruin Graph Assembly algorithm approach;  ######
###### the k-mers are weighted in the graph; sources for this approach are below 	   ######
###### 1. https://www.cs.jhu.edu/~langmea/resources/lecture_notes/assembly_dbg.pdf 	 ######
###### 2. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4719071/						           ######
###########################################################################################

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
  
  # Initialise dictionary to count the number of in-degree and out-degree 
  # per vertex in an Eulerian path
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

eulerian.path <- function(newList, 
                          timer = timer # limit time to run DBG assembly attempt
                          ){
  # re-define objects from previous functions' output
  dict           <- newList[[1]]
  suffix         <- newList[[2]]
  balanced.count <- newList[[3]]
  # create a copy of the dictionary for use in the pathway explorations
  dict.copy <- dict
  
  # obtain all possible starting points of the genome
  start.point <- list()
  for(i in 1:length(balanced.count)){
    if(balanced.count[[i]]==-1){
      start.point <- c(start.point, names(balanced.count)[i])
    }
  }
  
  # initialise list of reconstructed genome paths with the highest score
  pathways <- list()
  # initialise text progress bar
  pb <- txtProgressBar(min = 1, max = length(start.point)+1, style = 3)
  start.point.index = 1
  
  print("Traversing all paths in the rooted tree...", quote = F)
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
        #####################
        ### iteration one ###
        #####################
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
        ################################
        ### iteration two and onward ###
        ################################
        time1 <- Sys.time()
        while(length(branch.point)>=0){
          # reset dict
          dict <- dict.copy
          # re-initialise current path with starting point
          path <- start.point[[start.point.index]]
          i = 1
          while(TRUE){
            time2 <- Sys.time()
            if(as.numeric(difftime(time2, time1 , units = 'mins'))>=timer){
              # If the DBG is too large, the execution time is too long
              # Hence, skip current iteration and restart everything
              # with a new sequence
              print("De Bruijn Graph too large. Skip to next iteration.", 
                    quote = FALSE)
              return()
            }
            
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
          # if no more branches to traverse with current starting point save all
          # generated strings from current starting point to another list; 
          # restart everything with different starting point
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