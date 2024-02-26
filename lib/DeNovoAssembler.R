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

        #' @field Table in Data.table format of all k-mers, break probabilities and counts.
        table_read_kmer_prob = NULL,

        initialize = function(seq_len, read_len, kmer, dbg_kmer, seed, 
                              only_kmers_from_reads, ind, action, industry_standard, 
                              coverage_target, plot_results, save_read_files){
            if(!missing(only_kmers_from_reads)) private$only_kmers_from_reads <- only_kmers_from_reads
            if(!missing(industry_standard)) private$industry_standard <- industry_standard
            if(!missing(plot_results)) private$plot_results <- plot_results
            if(!missing(save_read_files)) private$save_read_files <- save_read_files
            
            super$initialize(
                seq_len = seq_len,
                read_len = read_len,
                coverage_target = coverage_target,
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
        #' @param total_iters total number of experiments to run.
        #' @return None.
        run_assembler = function(bins = 6, total_iters = 100){
            start.time <- Sys.time()
            cur.msg <- paste0("Running a p.o.p simulation with read len: ", private$read_len,  
                              ". Experiment: ", private$ind, "/", total_iters)
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(cur.msg, l, "\n", sep = "")

            self$get_reads()
            if(!private$only_kmers_from_reads){
                if(private$industry_standard){
                    private$get_assemblies_industry_standard()
                } else {
                    private$get_kmers_from_reads()
                    private$get_assemblies()

                    if(private$plot_results) private$quick_plots(bins = bins)
                }
                self$results <- private$score_solutions()
                private$save_results()
            } else {
                private$count_read_kmers()
            }

            # remove dir contents if industry standard approach is used
            # if(private$industry_standard & (private$ind == total_iters)){
            if(private$ind == total_iters){
                dirs.to.delete <- list.files(
                    path = "../data/reads",
                    pattern = "exp_",
                    full.names = TRUE
                )
                unlink(paste0(dirs.to.delete), recursive = TRUE)
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
        #' @field save_read_files Boolean. If TRUE, will save all ultrasonicated reads as files.
        save_read_files = TRUE,

        #' @field only_kmers_from_reads Boolean. If TRUE, will only get kmers from ultrasonicated reads.
        only_kmers_from_reads = FALSE,

        #' @field industry_standard Boolean. If TRUE, will use industry standard algorithms for experiments.
        industry_standard = FALSE,

        #' @field plot_results Boolean. If TRUE, will plot results for each experimental run. 
        plot_results = FALSE,

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
        #' Extract kmers of breakage length from sonicated sequencing reads.
        #' @return None.
        count_read_kmers = function(){
            # progress message
            t1 <- Sys.time()
            cur.msg <- "Extracting k-mers from sequencing reads"
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))

            # obtain all k-mers per read
            read.kmers <- lapply(1:length(self$sequencing_reads$read_one), function(i){
                end <- nchar(self$sequencing_reads$read_one[i])-self$kmer+1
                return(substring(
                    self$sequencing_reads$read_one[i],
                    first = 1:end, 
                    last = (1:end)+self$kmer-1
                ))
            })
            all_read_kmers <- unlist(read.kmers, use.names = FALSE)
            
            self$read_kmers <- table(all_read_kmers)
            self$read_kmers <- as.data.table(self$read_kmers)
            setnames(self$read_kmers, c("kmer", "count"))

            # find associated probability of breaking
            count_vals <- self$read_kmers$count[match(
                self$df_prob[[paste0("kmer_", self$kmer)]]$kmer, self$read_kmers$kmer
            )]
            self$table_read_kmer_prob <- copy(self$df_prob[[paste0("kmer_", self$kmer)]])
            self$table_read_kmer_prob[, count := count_vals]
            self$table_read_kmer_prob[, count := ifelse(is.na(count), 0, count_vals)]

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")
        },

        #' @description
        #' Create weighted de bruijn graph from k-mers using industry standard algorithms.
        #' @return None.
        get_assemblies_industry_standard = function(){
            # progress message
            t1 <- Sys.time()
            cur.msg <- "Running velvet de novo genome assembler"
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))

            # run velvet algorithm
            # velvet h
            velveth.command <- paste0(
                "../tests/velvet/velveth ",
                
                "../data/reads/exp_", private$ind, " ", 
                
                self$dbg_kmer, " ",
                
                "-shortPaired -fasta -separate ", 
                "../data/reads/exp_", private$ind, "/read_1",
                "_SeqLen-", format(private$seq_len, scientific = FALSE),
                "_SeqSeed-", private$seed,
                "_ReadLen-", private$read_len,
                "_DBGKmer-", self$dbg_kmer,
                ".fasta ",

                "../data/reads/exp_", private$ind, "/read_2",
                "_SeqLen-", format(private$seq_len, scientific = FALSE),
                "_SeqSeed-", private$seed,
                "_ReadLen-", private$read_len,
                "_DBGKmer-", self$dbg_kmer,
                ".fasta "
            )
            suppressWarnings(system(
                velveth.command, intern = TRUE,
                ignore.stdout = TRUE, ignore.stderr = TRUE
            ))

            # velvet g
            velvetg.command <- paste0(
                "../tests/velvet/velvetg ",
                
                "../data/reads/exp_", private$ind, " ",
                
                "-exp_cov auto ",
                "-cov_cutoff auto ",
                "-scaffolding yes"
            )
            suppressWarnings(system(
                velvetg.command, intern = TRUE,
                ignore.stdout = TRUE, ignore.stderr = TRUE
            ))

            # read contig files
            self$results <- Biostrings::readDNAStringSet(paste0(
                "../data/reads/exp_", private$ind,
                "/contigs.fa"
            )) %>% paste()

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")
        },

        #' @description
        #' Create weighted de bruijn graph from k-mers.
        #' @return None.
        get_assemblies = function(){
            # progress message
            t1 <- Sys.time()
            cur.msg <- "Running DBG de novo genome assembler"
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))

            # get contigs
            self$contigs <- get_contigs(
                read_kmers = self$read_kmers, 
                dbg_kmer = self$dbg_kmer,
                seed = private$seed
            )

            # brute-force assemble contigs
            self$results <- assemble_contigs(
                contig_matrix = self$contigs,
                dbg_kmer = self$dbg_kmer
            )

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
                    "_SeqLen-", format(private$seq_len, scientific = FALSE),
                    "_SeqSeed-", private$seed,
                    "_ReadLen-", private$read_len,
                    "_DBGKmer-", self$dbg_kmer,
                    "_kmer-", self$kmer, 
                    "_IndustryModel-", private$industry_standard,
                    ".RData"
                )
            )
            fwrite(
                self$results,
                file = paste0(
                    "../data/results/exp_", private$ind, 
                    "/SolutionsTable",
                    "_SeqLen-", format(private$seq_len, scientific = FALSE),
                    "_SeqSeed-", private$seed,
                    "_ReadLen-", private$read_len,
                    "_DBGKmer-", self$dbg_kmer,
                    "_kmer-", self$kmer, 
                    "_IndustryModel-", private$industry_standard,
                    ".csv"
                ),
                showProgress = FALSE
            )

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")
        },

        #' @description
        #' Score each contig solution using random and calcualted k-meric breakage probabilities.
        #' @return Evaluated solutions in Data.Table format.
        score_solutions = function(){
            # progress message
            t1 <- Sys.time()
            cur.msg <- "Evaluating each de novo assembled solution"
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))

            results <- lapply(c(FALSE, TRUE), function(random_prob){
                if(random_prob){
                    bp_prob <<- rep(
                        1 / length(self$df_prob$all$prob),
                        length(self$df_prob$all$prob)
                    )
                } else {
                    bp_prob <<- self$df_prob$all$prob
                }

                # calculate breakage scores
                breakscore_results <- calc_breakscore(
                    path = self$results, 
                    sequencing_reads = self$sequencing_reads$read_one, 
                    true_solution = self$dbg_summary$genome_seq,
                    kmer = self$kmer, 
                    bp_kmer = self$df_prob$all$kmer,
                    bp_prob = bp_prob
                )
                breakscore_results <- as.data.table(breakscore_results)

                # sort by break score
                setorder(breakscore_results, -bp_score)
                if(private$industry_standard){
                    results <- breakscore_results[breakscore_results$path_prob_dist_startpos != -1]
                } else {
                    results <- breakscore_results
                    results[, path_prob_dist_startpos := 0]
                }

                if(nrow(results) > 0){
                    # Kulback-Leibler Divergence
                    kl.mat <- sapply(1:nrow(results), function(x){
                        # statistical test of similarity
                        a_solution <- results$path_prob_dist[[x]]

                        # if(private$industry_standard){
                        #     # Tests require both vectors to be of the same length.
                        #     # Thus, truncate the longer one.
                        #     start.pos <- results$path_prob_dist_startpos[x]+1 # 0-index in cpp
                        #     true_solution_prob_truncated <- self$kmer_from_seq[start.pos:(
                        #         start.pos+nchar(results$sequence[x])-1
                        #     )]

                        #     # Discard edge probabilities
                        #     first.ind <- self$kmer/2
                        #     last.ind <- length(a_solution)-first.ind

                        #     if(last.ind < 0){
                        #         PQ <- NA
                        #     } else {
                        #         a_solution_truncated <- a_solution[first.ind:last.ind]
                        #         true_solution_prob_truncated <- true_solution_prob_truncated[first.ind:last.ind]
                        #     }
                        # } else {
                        #     # Tests require both vectors to be of the same length.
                        #     # Thus, truncate the longer one.
                        #     min_length <- min(length(a_solution), length(self$kmer_from_seq))
                        #     true_solution_prob_truncated <- self$kmer_from_seq[1:min_length]
                        #     a_solution_truncated <- a_solution[1:min_length]
                        # }

                        # # remove any nas if present
                        # to_discard <- which(is.na(a_solution_truncated) | is.na(true_solution_prob_truncated))
                        # if(length(to_discard) > 0){
                        #     a_solution_truncated <- a_solution_truncated[-to_discard]
                        #     true_solution_prob_truncated <- true_solution_prob_truncated[-to_discard]
                        # }

                        # # normalise distributions to sum to one.
                        # a_solution_adj <- a_solution_truncated / sum(a_solution_truncated)
                        # true_solution_prob_truncated_adj <- true_solution_prob_truncated / sum(true_solution_prob_truncated)

                        # # Kulback-Leibler Divergence: divergence of P's distribution from Q's
                        # PQ <- rbind(a_solution_adj, true_solution_prob_truncated_adj)  
                        # PQ <- philentropy::KL(PQ) %>% suppressMessages()
                        # PQ <- unname(PQ)

                        # remove any nas if present
                        to_discard <- which(is.na(a_solution))
                        if(length(to_discard) > 0) a_solution <- a_solution[-to_discard]

                        ks_test <- ks.test(
                            x = a_solution,
                            y = self$kmer_from_seq, 
                            alternative = "two.sided"
                        ) %>% suppressWarnings()
                        KS <- unname(ks_test$statistic)

                        return(KS)
                    })
                    results[, stat_test_KS := kl.mat]
                    all.ks.stat <- results$stat_test_KS
                    
                    # as we have many contigs, we need to see if all contigs cover the ref genome
                    full_genome_granges <- GRanges(
                        seqnames = "seq_1",
                        ranges = IRanges(start = 1, end = private$seq_len)
                    )

                    solutions_coverage <- tibble(
                        seqnames = "seq_1",
                        start = results$path_prob_dist_startpos,
                        end = results$path_prob_dist_startpos + results$sequence_len
                        ) %>% 
                        plyranges::as_granges() %>% 
                        reduce(.)
                    uncovered_regions <- setdiff(full_genome_granges, solutions_coverage)
                    covered_region <- (1 - (sum(width(uncovered_regions)) / private$seq_len)) * 100
                } else {
                    results <- breakscore_results
                    all.ks.stat <- NA
                    covered_region <- 0
                }

                # only need to keep the first row because we are summing all contigs anyways
                results[, `:=`(
                    stat_test_KS = all.ks.stat, 
                    path_prob_dist = NULL,
                    path_prob_dist_startpos = NULL,
                    contig_frac_len = covered_region
                )]

                return(results)
            })

            c_results <- dplyr::inner_join(
                x = as_tibble(results[[1]]),
                y = as_tibble(results[[2]]),
                by = c(
                    "sequence", "sequence_len", 
                    "kmer_breaks", "contig_frac_len",
                    "lev_dist_vs_true"
                ),
                suffix = c("_true", "_random")
            )

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")
            
            return(as.data.table(c_results))
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
                "Sequence length: ", format(private$seq_len, scientific = FALSE), ". ",
                "Sequence set.seed: ", private$seed, ". ",
                "Experiment: ", private$ind, ". ",
                "Nr. of solutions: ", nrow(self$results)
            )
            pdf(
                file = paste0(
                    "../figures/exp_", private$ind, 
                    "/Breakscore-vs-Levdist",
                    "_SeqLen-", format(private$seq_len, scientific = FALSE),
                    "_SeqSeed-", private$seed,
                    "_ReadLen-", private$read_len,   
                    "_DBGKmer-", self$dbg_kmer,                 
                    "_kmer-", self$kmer, 
                    "_IndustryModel-", private$industry_standard,
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