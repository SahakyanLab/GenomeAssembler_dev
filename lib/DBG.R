DBG <- R6::R6Class(
    classname = "DBG",
    inherit = GenerateReads,
    public = list(
        #' @field sequencing_reads Character vector of all sequencing reads.
        sequencing_reads = NULL,

        #' @field read_kmers Character vector of all k-mers extracted
        #'  from sequencing reads.
        read_kmers = NULL,

        #' @field dbg_kmer Numeric vector of the k-mer sizes for use in the 
        #'  de bruijn graph assembly process.
        dbg_kmer = 9,

        #' @field pathways Character vector of all de novo assembled genomes.
        pathways = NULL,

        #' @field results Data.table of assembled solutions and their scores.
        results = NULL,

        initialize = function(seq_len, read_len, G_cont, C_cont, 
                              A_cont, kmer, dbg_kmer, seed, action){
            if(!missing(dbg_kmer)) self$dbg_kmer <- dbg_kmer
            super$initialize(
                seq_len = seq_len,
                read_len = read_len,
                G_cont = G_cont,
                C_cont = C_cont,
                A_cont = A_cont,
                kmer = kmer,
                seed = seed,
                action = action
            )
        },

        #' @description
        #' Run de novo genome assembly process with de bruijn graph data structure.
        #' @return None.
        run_assembler = function(){
            start.time <- Sys.time()
            cur.msg <- "Running a proof of principle simulation"
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(cur.msg, l, "\n", sep = "")

            self$get_reads()
            private$get_kmers_from_reads()
            private$dbg_from_kmers()
            private$get_balance_count()
            private$get_eulerian_path()
            private$compare_break_scores()

            # time taken for full processing for this experiment
            final.t <- Sys.time() - start.time
            cat(paste(c(rep("-", 70), "\n"), collapse = ""))
            cat("Final time taken:", signif(final.t[[1]], digits = 3), 
                attr(final.t, "units"), "\n")
            cat(paste(c(rep("-", 70), "\n"), collapse = ""))
        }
    ), 
    private = list(
        #' @field dict List. key-value pairs of prefix and suffix for each k-mer.
        dict = NULL,

        #' @field prefix Character vector. Prefix for each k-mer.
        prefix = NULL,

        #' @field suffix Character vector. Suffix for each k-mer.
        suffix = NULL,

        #' @field balanced_count List. Counts the numbers of in-degree and out-degree 
        #'  per vertex of the Eulerian path. Keys are k-mer prefixes and values are
        #'  the in/out-degree counts.
        balanced_count = NULL,

        #' @description
        #' Extract kmers from sonicated sequencing reads.
        #' @return None.
        get_kmers_from_reads = function(){
            # progress message
            t1 <- Sys.time()
            cur.msg <- "Extracting k-mers from sequencing reads"
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))

            reads <- fread(
                "../data/reads/reads.txt", 
                header = FALSE,
                showProgress = FALSE
            )
            self$sequencing_reads <- unlist(reads$V1)

            # obtain all k-mers per read
            read.kmers <- lapply(1:length(self$sequencing_reads), function(i){
                sequence <- unlist(strsplit(self$sequencing_reads[i], split = ""))
                end <- length(sequence)-self$dbg_kmer+1
                dbg.kmer <- substring(
                    text = paste(sequence, collapse = ""),
                    first = 1:end, 
                    last = (1:end)+self$dbg_kmer-1)
                return(dbg.kmer)
            })
            read.kmers <- unlist(read.kmers)

            # obtain k-mers
            kmer.df <- as.data.table(table(read.kmers))
            self$read_kmers <- unlist(
                lapply(kmer.df$read.kmers, as.character),
                use.names = FALSE
            )

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")
        },

        #' @description
        #' Create de bruijn graph from k-mers
        #' @return None.
        dbg_from_kmers = function(){
            # progress message
            t1 <- Sys.time()
            cur.msg <- "Create de bruijn graph from k-mers"
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))

            # define prefix for each k-mer
            private$prefix <- sapply(1:length(self$read_kmers), function(i){
                substring(
                    text = self$read_kmers[i], 
                    first = 1, 
                    last = self$dbg_kmer-1
                )
            })
            # define suffix for each k-mer
            private$suffix <- sapply(1:length(self$read_kmers), function(i){
                substring(
                    text = self$read_kmers[i], 
                    first = 2, 
                    last = self$dbg_kmer
                )
            })
            
            # Initialise dictionary
            private$dict <- vector(
                mode = "list", 
                length = length(unique(private$prefix))
            )
            
            # set keys as prefix
            names(private$dict) <- unique(private$prefix)
            
            # fill dict with pre- and suffixes
            for(i in 1:length(self$read_kmers)){
                if(match(private$prefix[i], names(private$dict)) == i){
                    private$dict[i] <- private$suffix[i]
                } else {
                    ind <- match(private$prefix[i], names(private$dict))
                    private$dict[[ind]] <- c(private$dict[[ind]], private$suffix[i])
                }
            }

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")
        },
        
        #' @description
        #' Count the number of in-degree and out-degree per vertex
        #' in the Eulerian path.
        #' @return None.
        get_balance_count = function(){
            # progress message
            t1 <- Sys.time()
            cur.msg <- "Count in/out-degrees per vertex in Eulerian path"
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))

            private$balanced_count <- vector(
                mode = "list", 
                length = length(private$dict)
            )
            
            # set keys as prefix
            names(private$balanced_count) <- names(private$dict)
            
            # set values to zero
            private$balanced_count <- lapply(private$balanced_count, function(i) {
                private$balanced_count[i] <- 0
            })
            prefix.table <- as.data.table(table(private$prefix))
            suffix.table <- as.data.table(table(private$suffix))

            for(i in 1:length(private$balanced_count)){
                private$balanced_count[which(as.character(
                    prefix.table$V1)[i] == names(private$balanced_count)
                )] <- -prefix.table$N[i]
            }
            
            # fill dictionary with the total number of in-degree and out-degree
            # per vertex in the Eulerian trail
            for(i in 1:dim(suffix.table)[1]){
                possibleError <- tryCatch(
                    private$balanced_count[which(
                        as.character(suffix.table$V1)[i] == names(private$balanced_count)
                    )][[1]],
                    error = function(e){e}
                )
                
                if(!inherits(possibleError, "error")){
                    val <- private$balanced_count[which(as.character(
                        suffix.table$V1)[i] == names(private$balanced_count)
                    )][[1]] + suffix.table$N[i]
                    private$balanced_count[which(as.character(
                        suffix.table$V1)[i] == names(private$balanced_count)
                    )] <- val 
                } else {
                    list1 <- as.list(1)
                    names(list1) <- as.character(suffix.table$V1)[i]
                    private$balanced_count <- c(private$balanced_count, list1)
                    rm(list1)
                }
            }

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")
        },

        #' @description
        #' Traverse along all branches of the rooted tree.
        #' @return None.
        get_eulerian_path = function(){
            t1 <- Sys.time()
            cur.msg <- "Traversing along the de bruijn graph"
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))
            # cat(cur.msg, l, "\n", sep = "")

            # create a copy of the dictionary for use in the pathway explorations
            dict <- copy(private$dict)
            timer <- 100

            # obtain all possible starting points of the genome
            start.point <- names(private$balanced_count)[which(
                private$balanced_count == -1
            )]

            # initialise list of reconstructed genome paths with the highest score
            pathways <- list()
            # initialise text progress bar
            # pb <- txtProgressBar(min = 1, max = length(start.point)+1, style = 3)
            start.point.index = 1
            
            while(start.point.index<=length(start.point)){
                # initialise list to store number of branching points per step
                branch.point <- list()
                # initialise list to save the reconstructed genome path
                all.paths <- list()
                # initialise current path with starting point
                path <- start.point[start.point.index]
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
                                all.paths <- c(all.paths, private$reconstruct_genome(path))
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
                            dict <- copy(private$dict)
                            # re-initialise current path with starting point
                            path <- start.point[start.point.index]
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
                                        dict[last.item][[1]] <- 
                                            dict[last.item][[1]][-which(
                                                new.item==dict[last.item][[1]]
                                        )]
                                    } else {
                                        dict[last.item] <- NULL
                                    }
                                } else {
                                    # save reconstructed genome path into list
                                    all.paths <- c(all.paths, private$reconstruct_genome(path))
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
                # setTxtProgressBar(pb, start.point.index)
            }
            # close(pb)
            self$pathways <- unlist(pathways, use.names = FALSE)

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")
        },

        #' @description
        #' Reconstruct genome from assembly
        #' @return Character vector of reconstructed genome.
        reconstruct_genome = function(path){
            # initialise vector to save genome path
            genome <- character()
            # reconstruct genome path for a given string
            len <- nchar(path[1])
            genome <- c(
                substring(text = path[1], first = 1:len, last = 1:len),
                substring(text = path[2:length(path)], first = len, last = len)
            )
            genome <- paste(genome, collapse = "")
            return(genome)
        },

        #' @description
        #' Aligns de novo solutions with sequencing reads. Calculates an 
        #' alignment score based on the breakage probabilities.
        #' @return None.
        compare_break_scores = function(){     
            t1 <- Sys.time()
            cur.msg <- "Calculating alignment scores based on breakage probabilities"
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(cur.msg, l, "\n", sep = "")
                   
            res <- pbapply::pblapply(1:length(self$pathways), function(path){
                # counter column for number of occurrence of k-mer breakages
                self$df_prob$freq <- 0
                self$df_prob$bp.count <- 0
                for(i in 1:length(self$sequencing_reads)){
                    read.result <- gregexpr(
                        pattern = self$sequencing_reads[i], 
                        text = self$pathways[path]
                    )
                    read.start <- read.result[[1]][1]
                    if(read.start != (-1)){
                        # extract attributes from read and reference sequence overlap
                        read.length <- attr(read.result[[1]], 'match.length')
                        read.end <- read.start+read.length-1
                        
                        # obtain the k-mers that are broken on either end of the read
                        if(read.start == 1){
                            # if k-mer is at the start of the genome
                            kmer.start <- substring(
                                text = self$pathways[path],
                                first = 1,
                                last = self$kmer
                            )
                        } else {
                            kmer.start <- substring(
                                text = self$pathways[path],
                                first = read.start-self$kmer/2,
                                last = read.start+self$kmer/2-1
                            )
                        }
                        
                        if(nchar(self$pathways[path]) == read.end){
                            # if k-mer is at the end of the genome
                            kmer.end <- substring(
                                text = self$pathways[path],
                                first = read.end-self$kmer/2+1
                            )
                            if(!grepl(pattern = "N", x = kmer.end)){
                                kmer.end <- paste0(c(
                                    kmer.end,
                                    rep("N", self$kmer/2)
                                ), collapse = "")
                            }
                        } else {
                            kmer.end <- substring(
                                text = self$pathways[path],
                                first = read.end-self$kmer/2+1,
                                last = read.end+self$kmer/2
                            )
                        }
                        
                        # obtain index of the above k-mers
                        if(kmer.start %in% self$df_prob$kmer){
                            start.ind <- match(kmer.start, self$df_prob$kmer) 
                        } else {
                            kmer.start <- paste(Biostrings::reverseComplement(
                                Biostrings::DNAStringSet(kmer.start)
                            ))
                            start.ind <- match(kmer.start, self$df_prob$kmer)
                        }
                        
                        if(kmer.end %in% self$df_prob$kmer){
                            end.ind <- match(kmer.end, self$df_prob$kmer)
                        } else {
                            kmer.end <- paste(Biostrings::reverseComplement(
                                Biostrings::DNAStringSet(kmer.end)
                            ))
                            end.ind <- match(kmer.end, self$df_prob$kmer)
                        }
                        
                        if(!is.na(start.ind)){
                            # increase counter of given k-mer occurrence in data frame
                            self$df_prob$freq[start.ind] <- self$df_prob$freq[start.ind]+1
                            # count number of inferred breakpoints
                            self$df_prob$bp.count[start.ind] <- self$df_prob$bp.count[start.ind]+2
                        }
                        if(!is.na(end.ind)){
                            # increase counter of given k-mer occurrence in data frame
                            self$df_prob$freq[end.ind] <- self$df_prob$freq[end.ind]+1
                            # count number of inferred breakpoints
                            self$df_prob$bp.count[end.ind] <- self$df_prob$bp.count[end.ind]+2
                        }
                    }

                    if(i == length(self$sequencing_reads)){
                        bp.score <<- sum(self$df_prob$prob*self$df_prob$freq, na.rm = TRUE)
                        kmer.count <<- sum(self$df_prob$bp.count, na.rm = TRUE)
                        kmer.break <<- sum(self$df_prob$freq, na.rm = TRUE)
                    }
                }
                return(data.table(
                    denovo.result = self$pathways[path],
                    denovo.len = nchar(self$pathways[path]),
                    bp.score = bp.score,
                    # bp.score.norm = bp.score/kmer.count,
                    bp.score.norm = bp.score/nchar(self$pathways[path]),
                    kmer.count = kmer.count,
                    kmer.break = kmer.break
                ))
            })
            res <- rbindlist(res)

            # only accept the top 50% of the de novo assembled solutions
            setorder(res, -bp.score.norm)
            self$results <- res[complete.cases(res)]
        }
    )
)