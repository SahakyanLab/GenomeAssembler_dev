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

        #' @field dbg_kmer_table Data.table of k-mer counts generated from 
        #'  sequencing read and their weightings based on the frequency.
        dbg_kmer_table = NULL,

        #' @field pathways Character vector of all de novo assembled genomes.
        pathways = NULL,

        #' @field results Data.table of assembled solutions and their scores.
        results = NULL,

        initialize = function(seq_len, read_len, G_cont, C_cont, A_cont, 
                              kmer, dbg_kmer, seed, action, uniform_prob){
            super$initialize(
                seq_len = seq_len,
                read_len = read_len,
                G_cont = G_cont,
                C_cont = C_cont,
                A_cont = A_cont,
                kmer = kmer,
                dbg_kmer = dbg_kmer,
                seed = seed,
                action = action,
                uniform_prob = uniform_prob
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
                end <- nchar(self$sequencing_reads[i])-self$dbg_kmer+1
                return(substring(
                    self$sequencing_reads[i],
                    first = 1:end, 
                    last = (1:end)+self$dbg_kmer-1
                ))
            })
            read.kmers <- unlist(read.kmers, use.names = FALSE)

            # obtain k-mers
            kmer.df <- as.data.table(table(read.kmers))
            kmer.df[, N.weights := N/sum(N, na.rm = TRUE)]
            self$dbg_kmer_table <- kmer.df
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
            # cat(paste0(cur.msg, l))
            cat(cur.msg, l, "\n", sep = "")

            # test_dict = list(
            #     "NGTG" = "GTGT",
            #     "GTGT" = "TGTA",
            #     "TGTA" = c("GTAC", "GTAT"),
            #     "GTAC" = "TACC",
            #     "GTAT" = "TATC",
            #     "TACC" = "ACCT",
            #     "TATC" = c("ATCG", "ATCC"),
            #     "ATCC" = "TCCG"
            # )

            # ############################
            # initialise list to store number of branching points per step
            # dict=copy(test_dict)
            # queue_dict <- list()
            # queue_path <- list()
            # queue_i <- list()
            # already_traversed <- test_dict[lengths(test_dict) > 1]
            # ############################

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
            }

            # obtain all possible starting points of the genome
            start.point <- names(private$balanced_count)[which(
                private$balanced_count == -1
            )]

            # traverse all paths
            all.paths <- pbapply::pblapply(1:length(start.point), function(x){
                # create a copy of the dictionary for use in the pathway explorations
                dict <- copy(private$dict)
                queue_dict <- list()
                queue_path <- list()
                queue_i <- list()
                already_traversed <- private$dict[lengths(private$dict) > 1]

                if(length(already_traversed) > 0){
                    already_traversed_index <- lapply(1:length(already_traversed), function(ind){
                        temp <- as.data.table(already_traversed[[ind]])
                        temp[, key := names(already_traversed[ind])]
                        setnames(temp, c("value", "key"))
                        return(temp)
                    })
                    already_traversed_index <- rbindlist(already_traversed_index)
                    already_traversed_index[, index := 0]
                    setcolorder(already_traversed_index, c("key", "value", "index"))
                }
                
                # initialise list to save the reconstructed genome path
                all.paths <- list()
                # initialise current path with starting point
                path <- start.point[x]
                i <- 1

                while(TRUE){
                    last.item <- tail(path, n = 1)
                    all.items <- dict[last.item][[1]]
                    cur.item <- all.items[1]

                    if(length(all.items) > 0){
                        path <- c(path, cur.item)
                        if(length(all.items) > 1){
                            dict[last.item][[1]] <- 
                                dict[last.item][[1]][-which(
                                    cur.item == dict[last.item][[1]]
                                )]
                            
                            if(already_traversed_index[(last.item == key) & (cur.item == value), "index"] == 0){
                                # save list
                                already_traversed[last.item][[1]] <- 
                                    already_traversed[last.item][[1]][-which(
                                        cur.item == already_traversed[last.item][[1]]
                                    )]
                                already_traversed <- already_traversed[lengths(already_traversed) > 0]
                            } else {
                                already_traversed[last.item][[1]] <- dict[last.item][[1]]
                            }
                            already_traversed_index[(last.item == key) & (cur.item == value), "index"] <- i

                            if(!is.null(already_traversed[last.item][[1]])){
                                queue_dict <- c(queue_dict, list(dict))
                                queue_dict[[1]][last.item][[1]] <- c(dict[last.item][[1]], cur.item)
                                queue_path <- c(queue_path, list(head(path, n = -1)))
                                queue_i <- c(queue_i, list(i))
                            }
                        } else {
                            dict[last.item] <- NULL
                        }
                    } else {
                        # save reconstructed genome path into list
                        all.paths <- c(all.paths, reconstruct_genome(path))

                        # reset loop based on current queue
                        list.len <- length(queue_i)
                        if(list.len > 0){
                            i <- queue_i[[list.len]]
                            dict <- queue_dict[[list.len]]
                            path <- queue_path[[list.len]]

                            queue_dict <- queue_dict[-list.len]
                            queue_path <- queue_path[-list.len]
                            queue_i <- queue_i[-list.len]
                        } else {
                            break
                        }
                    }
                    i <- i+1
                }
                return(unlist(all.paths, use.names = FALSE))
            })
            self$pathways <- unlist(all.paths, use.names = FALSE)
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
                    bp.score.norm = bp.score/kmer.count,
                    # bp.score.norm = bp.score/nchar(self$pathways[path]),
                    kmer.count = kmer.count,
                    kmer.break = kmer.break
                ))
            })
            res <- rbindlist(res)

            # only accept the top 50% of the de novo assembled solutions
            setorder(res, -bp.score)
            self$results <- res[complete.cases(res)]
        }
    )
)
