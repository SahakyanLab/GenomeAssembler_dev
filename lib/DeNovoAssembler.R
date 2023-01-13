DeNovoAssembler <- R6::R6Class(
    classname = "DeNovoAssembler",
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

        initialize = function(seq_len, read_len, kmer, dbg_kmer, ncpu, seed, action){
            if(!missing(ncpu)) private$ncpu <- ncpu
            super$initialize(
                seq_len = seq_len,
                read_len = read_len,
                kmer = kmer,
                dbg_kmer = dbg_kmer,
                seed = seed,
                action = action
            )
        },

        #' @description
        #' Run de novo genome assembly process with de bruijn graph data structure.
        #' @param ind Numeric vector. Subsets the sequence from all sampled ref sequences.
        #' @return None.
        run_assembler = function(ind = 1){
            start.time <- Sys.time()
            cur.msg <- "Running a proof of principle simulation"
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(cur.msg, l, "\n", sep = "")

            self$get_reads(ind = ind)
            private$get_kmers_from_reads()
            private$wdbg_from_kmers()
            private$get_branching_nodes()
            private$get_contigs()
            private$get_many_scaffolds()
            private$compare_break_scores()
            private$save_results()

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

        #' @field ncpu Numeric vector of the number of cores to use 
        #'  for parallel execution.
        ncpu = 2,

        #' @description
        #' Extract kmers from sonicated sequencing reads.
        #' @return None.
        get_kmers_from_reads = function(){
            # progress message
            t1 <- Sys.time()
            cur.msg <- "Extracting k-mers from sequencing reads"
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))

            # read one
            read.one <- fread("../data/reads/read_one.txt", header = FALSE, showProgress = FALSE)
            self$sequencing_reads$read_one <- unlist(read.one$V1)

            # read two
            read.two <- fread("../data/reads/read_two.txt", header = FALSE, showProgress = FALSE)
            read.two <- unlist(read.two$V1)
            read.two <- paste(Biostrings::reverseComplement(Biostrings::DNAStringSet(read.two)))
            self$sequencing_reads$read_two <- read.two

            # obtain all k-mers per read
            read.kmers <- lapply(1:length(self$sequencing_reads$read_one), function(i){
                end <- nchar(self$sequencing_reads$read_one[i])-self$dbg_kmer+1
                return(substring(
                    self$sequencing_reads$read_one[i],
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
        #' Traverse as many paths as possible by combining the various contigs into scaffolds.
        #' @return None.
        get_many_scaffolds = function(){
            # progress message
            t1 <- Sys.time()
            cur.msg <- "Generating scaffolds from contigs"
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")

            N.fac <- factorial(length(self$contigs))
            N.default <- 20000
            N.rep <- which.min(c(N.default, N.fac))
            if(N.rep == 1){
                contig.matrix <- replicate(n = N.default, sample(self$contigs))
                contig.matrix <- split(contig.matrix, col(contig.matrix))
                contig.matrix <- unname(contig.matrix)
                contig.matrix <- unique(contig.matrix)
                contig.matrix <- c(
                    contig.matrix,
                    list(sort(self$contigs), self$contigs)
                )
            } else if(N.rep == 2){
                contig.matrix <- combinat::permn(self$contigs)
            }

            if(private$ncpu > 1){
                cat(paste0(cur.msg, l))
                `%op%` <- `%dopar%`
            } else {
                cat(cur.msg, l, "\n", sep = "")
                `%op%` <- `%do%`
                pb <- txtProgressBar(
                    min = 1, 
                    max = length(contig.matrix), 
                    style = 3
                )
            }

            # set-up cluster for parallel computation
            cl <- makeCluster(private$ncpu)
            registerDoParallel(cl)
            contig.results <- foreach(i = 1:length(contig.matrix), 
                                .combine="rbind",
                                .export=c(ls(globalenv()), "private", "self"),
                                .packages=c("foreach", "data.table", "dplyr", "R6"),
                                .inorder=TRUE)%op%{
                DeNovoAssembler$parent_env <- environment()
                if(private$ncpu < 2) setTxtProgressBar(pb, i)
                return(private$assemble_contigs(contigs = contig.matrix[[i]]))
            } %>% 
            suppressWarnings() # ignore 'already exporting variables' warning from foreach
            stopImplicitCluster()
            stopCluster(cl)
            if(private$ncpu == 1) close(pb)

            contig.results <- as.character(contig.results)
            contig.results <- unlist(contig.results, use.names = FALSE)
            self$scaffolds <- unique(contig.results)

            # retain only top 30% assemblies by length
            temp <- data.table(scaffolds = self$scaffolds, len = nchar(self$scaffolds))
            setorder(temp, -len)
            self$scaffolds <- temp$scaffolds[1:ceiling(nrow(temp)*0.3)]

            if(private$ncpu > 1){
                # time taken for full processing for this experiment
                total.time <- Sys.time() - t1
                cat("DONE! --", signif(total.time[[1]], 2), 
                    attr(total.time, "units"), "\n")
            }
        },

        #' @description
        #' Summary statistics of aligning all de novo solutions with sequencing reads.
        #' @return None.
        compare_break_scores = function(){
            t1 <- Sys.time()
            cur.msg <- "Calculating alignment scores based on breakage probabilities"
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")

            if(private$ncpu > 1){
                cat(paste0(cur.msg, l))
                `%op%` <- `%dopar%`
            } else {
                cat(cur.msg, l, "\n", sep = "")
                `%op%` <- `%do%`
                pb <- txtProgressBar(
                    min = 1, 
                    max = length(self$scaffolds), 
                    style = 3
                )
            }
            # set-up cluster for parallel computation
            cl <- makeCluster(private$ncpu)
            registerDoParallel(cl)
            res <- foreach(i = 1:length(self$scaffolds), 
                                .combine="rbind",
                                .export=c(ls(globalenv()), "private", "self"),
                                .packages=c("foreach", "data.table", "dplyr", "R6"),
                                .inorder=TRUE)%op%{
                DeNovoAssembler$parent_env <- environment()
                if(private$ncpu < 2) setTxtProgressBar(pb, i)
                return(private$analyse_results(path = self$scaffolds[i]))
            } %>% 
            suppressWarnings() # ignore 'already exporting variables' warning from foreach
            stopImplicitCluster()
            stopCluster(cl)
            if(private$ncpu == 1) close(pb)

            if(all(c("matrix", "array") == class(res))){
                bp.score <- as.numeric(res[,1])
                kmer.breaks <- as.numeric(res[,2])
            } else {
                bp.score <- res[[1]]
                kmer.breaks <- res[[2]]
            }
            denovo.len <- nchar(self$scaffolds)
            df.bp.score <- data.table(
                denovo.result = self$scaffolds,
                denovo.len = denovo.len,
                bp.score = bp.score,
                bp.score.norm.by.len = bp.score/denovo.len,
                bp.score.norm.by.break.freq = bp.score/kmer.breaks,
                kmer.breaks = kmer.breaks
            )
            setorder(df.bp.score, -bp.score)
            df.bp.score <- df.bp.score[complete.cases(df.bp.score)]
            df.bp.score <- df.bp.score[kmer.breaks != 0]

             # Levenshtein distance calculations
            lev.mat <- matrix(
                data = c(
                    df.bp.score$denovo.result,
                    rep(paste(self$genome_seq, collapse = ""), 
                        length(df.bp.score$denovo.result))
                ),
                ncol = 2
            )
            lev.dist <- as.numeric(LevenshteinLoop(lev.mat))
            df.bp.score[, lev.dist.vs.true := lev.dist]
            self$results <- df.bp.score
            self$dbg_summary$nr_of_solutions <- nrow(self$results)

            if(private$ncpu > 1){
                # time taken for full processing for this experiment
                total.time <- Sys.time() - t1
                cat("DONE! --", signif(total.time[[1]], 2), 
                    attr(total.time, "units"), "\n")
            }
        },

        #' @description
        #' Aligns de novo solutions with sequencing reads. Calculates an 
        #' alignment score based on the breakage probabilities.
        #' @param path Character vector of a genome sequence (one of the de novo solutions).
        #' @return None.
        analyse_results = function(path){
            get.score <- function(sequencing.reads, is.read.two){
                read.matches <- stringr::str_locate_all(
                    string = path,
                    pattern = sequencing.reads
                )
                read.matches.ind <- which(lengths(read.matches) > 0)
                if(length(read.matches.ind) > 0){
                    df.read.matches <- as.data.table(do.call(rbind, read.matches[read.matches.ind]))
                    reads <- rep(
                        sequencing.reads[read.matches.ind], 
                        lengths(read.matches[read.matches.ind])/2
                    )
                    path.len <- nchar(path)
                    df.read.matches[, read := reads]
                    df.read.matches[, any.ns := stringr::str_detect(string = read, pattern = "N")]
                    
                    if(is.read.two){
                        # read.two
                        df.read.matches[, `:=`(
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
                            kmer.end = substring(
                                text = path, 
                                first = right.kmer.start, 
                                last = right.kmer.end
                            ),
                            read = NULL
                        )]
                    } else {
                        # read one
                        df.read.matches[, `:=`(
                            left.kmer.start = ifelse(any.ns & (start <= self$kmer/2), 
                                1, start-self$kmer/2
                            ),
                            left.kmer.end = ifelse(any.ns & (start <= self$kmer/2), ifelse(
                                start < self$kmer/2, self$kmer,
                                start+self$kmer/2-1), start+self$kmer/2-1
                            )
                        )]
                        df.read.matches[, `:=`(
                            kmer.start = substring(
                                text = path, 
                                first = left.kmer.start, 
                                last = left.kmer.end
                            ),
                            read = NULL
                        )]
                    }

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

                        # true breakage probabilities
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
            }
            read.one.scores <- get.score(
                sequencing.reads = self$sequencing_reads$read_one, 
                is.read.two = FALSE
            )
            read.two.scores <- get.score(
                sequencing.reads = self$sequencing_reads$read_two, 
                is.read.two = TRUE
            )
            return(list(
                read.one.scores[[1]]+read.two.scores[[1]],
                read.one.scores[[2]]+read.two.scores[[2]]
            ))
        },

        #' @description
        #' Save results into csv files.
        #' @return None.
        save_results = function(){
            # progress message
            t1 <- Sys.time()
            cur.msg <- "Saving solutions as csv file and assembly stats as RDS file"
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))

            # as RDS file
            self$dbg_summary$genome_seq <- paste(self$genome_seq, collapse = "")

            dir.create(
                path = "../data/results",
                showWarnings = FALSE,
                recursive = TRUE
            )
            saveRDS(
                self$dbg_summary,
                file = paste0(
                    "../data/results/AssemblyStats_kmer-", self$kmer, 
                    "_SeqLen-", length(self$genome_seq), "_SeqSeed-", private$seed,
                    "_ReadLen-", private$read_len,
                    ".RData"
                )
            )

            fwrite(
                self$results,
                file = paste0(
                    "../data/results/SolutionsTable_kmer-", self$kmer, 
                    "_SeqLen-", length(self$genome_seq), "_SeqSeed-", private$seed,
                    "_ReadLen-", private$read_len,
                    ".csv"
                )
            )

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")
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