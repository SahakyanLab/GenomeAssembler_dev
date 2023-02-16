DeNovoAssembler <- R6::R6Class(
    classname = "DeNovoAssembler",
    inherit = GenerateReads,
    public = list(
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

        initialize = function(seq_len, read_len, kmer, dbg_kmer, 
                              seed, ind, action, reads_only, cores){
            if(!missing(reads_only)) private$reads_only <- reads_only
            if(!missing(cores)) private$cores <- cores
            super$initialize(
                seq_len = seq_len,
                read_len = read_len,
                kmer = kmer,
                dbg_kmer = dbg_kmer,
                seed = seed,
                action = action,
                ind = ind
            )
        },

        #' @description
        #' Run de novo genome assembly process with de bruijn graph data structure.
        #' @param bins Numeric vector of how many breaks to create when plotting results.
        #' @return None.
        run_assembler = function(bins = 10){
            start.time <- Sys.time()
            cur.msg <- paste0("Running a proof of principle simulation ", 
                              "for experiment: ", private$ind, "/", 1000)
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(cur.msg, l, "\n", sep = "")

            private$load_all_ref_seq()
            self$get_reads()
            if(!private$reads_only){
                private$get_kmers_from_reads()
                private$get_assemblies()
                private$save_results()
                private$quick_plots(bins = bins)
            }

            # time taken for full processing for this experiment
            final.t <- Sys.time() - start.time
            cat(paste(c(rep("-", 70), "\n"), collapse = ""))
            cat("Final time taken:", signif(final.t[[1]], digits = 3), 
                attr(final.t, "units"), "\n")
            cat(paste(c(rep("-", 70), "\n"), collapse = ""))
        }
    ), 
    private = list(
        #' @field reads_only Boolean. If TRUE, will only generate sequencing reads and 
        #'  not assembly the genome.
        reads_only = FALSE,

        #' @field cores Numeric vector of the number of cores to use in omp.h
        cores = 1,

        #' @description
        #' Extract kmers from sonicated sequencing reads.
        #' @return None.
        get_kmers_from_reads = function(){
            # progress message
            t1 <- Sys.time()
            cur.msg <- "Extracting k-mers from sequencing reads"
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))

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
        get_assemblies = function(){
            # progress message
            t1 <- Sys.time()
            cur.msg <- "Running DBG de novo genome assembler"
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))

            # assemble contigs
            self$scaffolds <- assemble_contigs(
                read_kmers = self$read_kmers, 
                dbg_kmer = self$dbg_kmer,
                seed = private$seed
            )

            # calculate breakage scores
            self$results <- calc_breakscore(
                path = self$scaffolds, 
                sequencing_reads = self$sequencing_reads$read_one, 
                true_solution = self$dbg_summary$genome_seq,
                kmer = self$kmer, 
                bp_kmer = self$kmer_ref$kmer,
                bp_prob = self$kmer_ref$prob,
                num_threads = private$cores
            )
            self$results <- as.data.table(self$results)

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")
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

            dir.create(
                path = paste0("../data/results/exp_", private$ind, "/"),
                showWarnings = FALSE,
                recursive = TRUE
            )
            saveRDS(
                self$dbg_summary,
                file = paste0(
                    "../data/results/exp_", private$ind, 
                    "/AssemblyStats",
                    "_SeqLen-", private$seq_len, 
                    "_SeqSeed-", private$seed,
                    "_ReadLen-", private$read_len,
                    "_kmer-", self$kmer, 
                    ".RData"
                )
            )
            fwrite(
                self$results,
                file = paste0(
                    "../data/results/exp_", private$ind, 
                    "/SolutionsTable",
                    "_SeqLen-", private$seq_len, 
                    "_SeqSeed-", private$seed,
                    "_ReadLen-", private$read_len,
                    "_kmer-", self$kmer, 
                    ".csv"
                )
            )

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")
        },

        #' @description
        #' Generate some quick plots after assembly results are obtained.
        #' @param bins Numeric vector of how many breaks to create when plotting results.
        #' @return None.
        quick_plots = function(bins){
            # true break score vs. levenshtein distance
            binnings <- cut(
                self$results$lev_dist_vs_true, 
                breaks = seq(
                    from = 0, 
                    to = max(self$results$lev_dist_vs_true, 
                             na.rm = TRUE), 
                    length.out = bins
                ), 
                include.lowest = TRUE
            )
            self$results[, bins := binnings]
            mytitle <- paste0(
                "Breakage probability scores vs. binned ",
                "levenshtein distance to the true solution"
            )
            mysubtitle <- paste0(
                "Kmer: ", self$kmer, ". ",
                "Sequence length: ", private$seq_len, ". ",
                "Sequence set.seed: ", private$seed, ". ",
                "Experiment: ", private$ind, ". ",
                "Nr. of solutions: ", nrow(self$results)
            )
            pdf(
                file = paste0(
                    "../figures/exp_", private$ind, 
                    "/Breakscore-vs-Levdist",
                    "_SeqLen-", private$seq_len, 
                    "_SeqSeed-", private$seed,
                    "_ReadLen-", private$read_len,                    
                    "_kmer-", self$kmer, 
                    "_exp-", private$ind,
                    ".pdf"
                ),
                width = 19, 
                height = 5
            )
            par(mfrow = c(1,3))
            par(mar=c(7,6,4,1))
            boxplot(
                bp_score ~ bins, 
                data = self$results, 
                frame = FALSE,
                col = "lightblue",
                border = "black",
                xlab = "",
                ylab = "Actual",
                main = "",
                las = 3
            )
            mtext(line=2.2, at=-0.07, adj=0, cex=1.1, mytitle)
            mtext(line=1, at=-0.07, adj=0, cex=0.9, mysubtitle)
            boxplot(
                bp_score_norm_by_len ~ bins, 
                data = self$results, 
                frame = FALSE,
                col = "lightblue",
                border = "black",
                xlab = "",
                ylab = "Normalised by length",
                main = "",
                las = 3
            )
            boxplot(
                bp_score_norm_by_break_freqs ~ bins, 
                data = self$results, 
                frame = FALSE,
                col = "lightblue",
                border = "black",
                xlab = "",
                ylab = "Normalised by nr of breaks",
                main = "",
                las = 3
            )
            plot.saved <- dev.off()
        }
    )
)