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

        #' @field contigs Character vector of all contigs generated through 
        #'   traversal of the de bruijn graph.
        contigs = NULL,

        #' @field scaffolds Character vector of all de novo assembled genomes.
        scaffolds = NULL,

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
            private$wdbg_from_kmers()
            private$get_branching_nodes()
            private$get_contigs()
            private$get_many_scaffolds()
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

        #' @field branching_nodes Character vector of all nodes that have a branch.
        branching_nodes = NULL,

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
            self$read_kmers <- unlist(read.kmers, use.names = FALSE)

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")
        },

        #' @description
        #' Create weighted de bruijn graph from k-mers
        #' @return None.
        wdbg_from_kmers = function(){
            # progress message
            t1 <- Sys.time()
            cur.msg <- "Create weighted de bruijn graph from k-mers"
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
            for(i in 1:length(private$prefix)){
                if(is.null(private$dict[[private$prefix[i]]])){
                    private$dict[[private$prefix[i]]] <- private$suffix[i]
                } else {
                    if(!private$suffix[i] %in% private$dict[[private$prefix[i]]]){
                        private$dict[[private$prefix[i]]] <- c(
                            private$dict[[private$prefix[i]]], 
                            private$suffix[i]
                        )
                    }
                }
            }

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")
        },
        
        #' @description
        #' Find the nodes that have more than one outgoing node.
        #' in the Eulerian path.
        #' @return None.
        get_branching_nodes = function(){
            # progress message
            t1 <- Sys.time()
            cur.msg <- "Count in/out-degrees per vertex in Eulerian path"
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))

            to.fill <- unique(c(unique(private$prefix), unique(private$suffix)))
            private$balanced_count <- vector(
                mode = "list", 
                length = length(to.fill)
            )
            
            # set keys as prefix
            names(private$balanced_count) <- to.fill

            # set values to zero
            private$balanced_count <- lapply(private$balanced_count, function(i){
                private$balanced_count[i] <- c(0, 0)
            })

            # get in-degree and out-degree per node
            for(i in 1:length(private$dict)){
                node <- names(private$dict[i])
                edges <- private$dict[[i]]

                # out-degree 
                private$balanced_count[[node]][2] <- 
                    private$balanced_count[[node]][2]+length(edges)

                # in-degree
                for(edge in edges){
                    private$balanced_count[[edge]][1] <- 
                        private$balanced_count[[edge]][1]+1
                }
            }
            private$branching_nodes <- sapply(private$balanced_count, function(node){
                return((node[1] == 1) & (node[2] == 1))
            })
            private$branching_nodes <- names(which(private$branching_nodes == FALSE))
            private$branching_nodes <- private$branching_nodes[which(
                private$branching_nodes %in% names(private$dict)
            )]

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")
        },

        #' @description
        #' Generate contigs by the paths from branching nodes
        #' @return None.
        get_contigs = function(){
            # progress message
            t1 <- Sys.time()
            cur.msg <- "Generating contigs from branching nodes."
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))

            contigs <- character()
            for(node in private$branching_nodes){
                for(edge in private$dict[[node]]){
                    path <- character()
                    current.node <- edge
                    while(!current.node %in% private$branching_nodes){
                        if(is.null(private$dict[[current.node]])) break
                        path <- c(path, current.node)
                        current.node <- private$dict[[current.node]]
                    }
                    path <- c(node, path, current.node)
                    contigs <- c(contigs, private$reconstruct_path(path))
                }
            }
            self$contigs <- unique(contigs)

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")
        },

        #' @description
        #' Assemble contigs into scaffolds.
        #' @param contigs Character vector of all contigs generated from branching nodes.
        #' @return Character vector of contigs.
        assemble_contigs = function(contigs){
            for(kmer in (self$dbg_kmer-1):1){
                len_changed <- TRUE
                while(len_changed){
                    temp <- length(contigs)
                    for(i in 1:length(contigs)){
                        if(contigs[i] == "") next
                        for(j in 1:length(contigs)){
                            if(contigs[i] != contigs[j]){
                                suffix <- private$substr_suffix(x = contigs[i], n = kmer)
                                prefix <- private$substr_prefix(x = contigs[j], n = kmer)
                                if(suffix == prefix){
                                    contigs[i] <- private$reconstruct_contig(
                                        first_contig = contigs[i], 
                                        second_contig = contigs[j], 
                                        n = kmer
                                    )
                                    contigs[j] <- ""
                                }
                            }
                        }
                    }
                    contigs <- contigs[contigs != ""]
                    if(temp == length(contigs)){
                        len_changed <- FALSE
                    }
                }
            }
            return(contigs)
        },

        #' @description
        #' Traverse as many paths as possible as combining the various contigs into scaffolds.
        #' @return None.
        get_many_scaffolds = function(){
            # progress message
            t1 <- Sys.time()
            cur.msg <- "Generating scaffolds from contigs"
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(cur.msg, l, "\n", sep = "")

            set.seed(private$seed)
            contig.matrix <- replicate(n = 3000, sample(self$contigs))
            contig.matrix <- cbind(
                contig.matrix, 
                matrix(c(self$contigs, sort(self$contigs)), ncol = 2)
            )
            contig.results <- pbapply::pblapply(1:ncol(contig.matrix), function(x){
                return(private$assemble_contigs(contigs = contig.matrix[,x]))
            })
            contig.results <- unlist(contig.results, use.names = FALSE)
            self$scaffolds <- unique(contig.results)
            # self$scaffolds <- self$scaffolds[which(
            #     nchar(self$scaffolds) <= length(self$genome_seq)
            # )]
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

            res <- pbapply::pblapply(self$scaffolds, function(path){
                read.matches <- stringr::str_locate_all(
                    string = path,
                    pattern = self$sequencing_reads
                )
                read.matches.ind <- which(lengths(read.matches) > 0)

                if(length(read.matches.ind) > 0){
                    df.read.matches <- as.data.table(do.call(rbind, read.matches[read.matches.ind]))
                    reads <- rep(
                        self$sequencing_reads[read.matches.ind], 
                        lengths(read.matches[read.matches.ind])/2
                    )
                    path.len <- nchar(path)
                    df.read.matches[, read := reads]
                    df.read.matches[, any.ns := stringr::str_detect(string = read, pattern = "N")]
                    df.read.matches[, `:=`(
                        left.kmer.start = ifelse(any.ns & (start <= self$kmer/2), 
                            1, start-self$kmer/2
                        ),
                        left.kmer.end = ifelse(any.ns & (start <= self$kmer/2), ifelse(
                            start < self$kmer/2, self$kmer,
                            start+self$kmer/2-1), start+self$kmer/2-1
                        ),
                        right.kmer.start = ifelse(any.ns & (start > self$kmer/2), ifelse(
                            path.len-end < self$kmer/2, path.len-self$kmer+1, 
                            end-self$kmer/2+1), end-self$kmer/2+1
                        ),
                        right.kmer.end = ifelse(any.ns & (start > self$kmer/2), ifelse(
                            path.len-end < self$kmer/2, path.len, 
                            end+self$kmer/2), end+self$kmer/2
                        )
                    )]
                    df.read.matches[, `:=`(
                        kmer.start = substring(
                            text = path, 
                            first = left.kmer.start, 
                            last = left.kmer.end
                        ),
                        kmer.end = substring(
                            text = path, 
                            first = right.kmer.start, 
                            last = right.kmer.end
                        ),
                        read = NULL
                    )]

                    to.keep <- c("kmer.start", "kmer.end")
                    kmers <- c(df.read.matches$kmer.start, df.read.matches$kmer.end)

                    # obtain index of forward and reverse complement k-mers
                    fwd.ind <- match(kmers, self$df_prob$kmer)
                    rev.ind <- is.na(fwd.ind)
                    fwd.ind <- fwd.ind[!is.na(fwd.ind)]
                    if(any(rev.ind)){
                        rev.comp <- paste(Biostrings::reverseComplement(
                            Biostrings::DNAStringSet(kmers[which(rev.ind)])
                        ))
                        rev.ind <- match(rev.comp, self$df_prob$kmer)
                    }
                    probs.to.extract <- c(fwd.ind, rev.ind)
                    probs.to.extract <- probs.to.extract[probs.to.extract != 0]
                    probs.to.extract <- probs.to.extract[complete.cases(probs.to.extract)]
                    
                    if(length(probs.to.extract) > 0){
                        # breakpoint score = prob of breaking k-mer * freq of k-mer broken
                        probs.to.extract <- as.data.table(table(probs.to.extract))
                        setnames(probs.to.extract, c("ind", "N"))
                        probs.to.extract[, names(probs.to.extract) := lapply(.SD, as.numeric)]
                        probs.to.extract[, prob := self$df_prob$prob[ind]]
                        probs.to.extract[, bp.score := prob*N]
                        bp.score <- sum(probs.to.extract$bp.score, na.rm = TRUE)

                        # nr of k-mers broken
                        kmer.breaks <- nrow(probs.to.extract)
                    } else {
                        bp.score <- 0
                        kmer.breaks <- 0
                    }
                } else {
                    bp.score <- 0
                    kmer.breaks <- 0
                }
                return(list(bp.score, kmer.breaks))
            })
            denovo.len <- nchar(self$scaffolds)
            bp.score <- sapply(res, `[[`, 1)
            kmer.breaks <- sapply(res, `[[`, 2)
            df.bp.score <- data.table(
                denovo.result = self$scaffolds,
                denovo.len = denovo.len,
                bp.score = bp.score,
                # bp.score.norm = bp.score/kmer.breaks,
                bp.score.norm = bp.score/denovo.len,
                kmer.breaks = kmer.breaks
            )
            setorder(df.bp.score, -bp.score)
            self$results <- df.bp.score[complete.cases(df.bp.score)]
        },

        #' @description
        #' Append second contig with first contig.
        #' @param first_contig Character vector of a contig.
        #' @param second_contig Character vector of another contig.
        #' @param n Numeric vector. Concatenating first_contig and second_contig
        #'  by their overlapping prefix and suffix k-mers.
        #' @return Character vector.
        reconstruct_contig = function(first_contig, second_contig, n){
            second.contig <- substring(
                text = second_contig,
                first = n+1,
                last = nchar(second_contig)
            )
            return(paste0(first_contig, second.contig, collapse = ""))
        },

        #' @description
        #' Helper function to get prefix of contig.
        #' @param x Character vector. A contig.
        #' @param n Numeric vector. Length of k-mer to take substring of x.
        #' @return Character vector.
        substr_prefix = function(x, n){
            return(substring(
                text = x, 
                first = 1,
                last = n
            ))
        },

        #' @description
        #' Helper function to get suffix of contig.
        #' @param x Character vector. A contig.
        #' @param n Numeric vector. Length of k-mer to take substring of x.
        #' @return Character vector.
        substr_suffix = function(x, n){
            return(substring(
                text = x, 
                first = nchar(x)-n+1, 
                last = nchar(x)
            ))
        },

        #' @description
        #' Reconstruct path from traversal of graph.
        #' @param path Character vector of a reconstructed genome path.
        #' @return Character vector.
        reconstruct_path = function(path){
            # initialise vector to save genome path
            genome <- character()
            # reconstruct genome path for a given string
            len <- nchar(path[1])
            genome <- c(
                substring(
                    text = path[1], 
                    first = 1:len, 
                    last = 1:len
                ),
                substring(
                    text = path[2:length(path)], 
                    first = len, 
                    last = len
                )
            )
            genome <- paste(genome, collapse = "")
            return(genome)
        }
    )
)
