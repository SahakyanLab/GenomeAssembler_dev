DBG <- R6::R6Class(
    classname = "DBG",
    inherit = GenerateReads,
    public = list(
        #' @field read_kmers Character vector of all k-mers extracted
        #'  from sequencing reads.
        read_kmers = NULL,

        #' @field dbg_kmer Numeric vector of the k-mer sizes for use in the 
        #' de bruijn graph assembly process.
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
        }
    ), 
    private = list(
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
        }
    )
)