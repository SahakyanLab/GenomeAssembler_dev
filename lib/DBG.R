DBG <- R6::R6Class(
    classname = "DBG",
    inherit = GenerateReads,
    public = list(
        #' @field read_kmers Character vector of all k-mers extracted
        #'  from sequencing reads.
        read_kmers = NULL,

        #' @field dbg_kmer Numeric vector of the k-mer sizes for use in the 
        #'  de bruijn graph assembly process.
        dbg_kmer = 9,

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
            self$get_reads()
            private$get_kmers_from_reads()
            private$dbg_from_kmers()
            private$get_balance_count()
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
            reads <- unlist(reads$V1)

            # obtain all k-mers per read
            read.kmers <- lapply(1:length(reads), function(i){
                sequence <- unlist(strsplit(reads[i], split = ""))
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
        }
    )
)