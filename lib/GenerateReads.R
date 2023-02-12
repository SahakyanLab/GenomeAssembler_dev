GenerateReads <- R6::R6Class(
    classname = "GenerateReads",
    public = list(
        #' @field genome_seq Character vector. Randomly generated reference genome.
        genome_seq = NULL,

        #' @field df_prob Data.table of kmers, probability ratios and enrichment
        #'  and depletion z-scores obtained through analysis of 
        #'  ultrasonication experiments.
        df_prob = NULL,

        #' @field kmer Numeric vector. Length of kmer for breakage probability. 
        #'  Kmer can only take values of c(4,6,8).
        kmer = 4,

        #' @field dbg_kmer Numeric vector of the k-mer sizes for use in the 
        #'  de bruijn graph assembly process.
        dbg_kmer = 9,

        #' @field dbg_summary List of the summary statistics of the assembly process.
        dbg_summary = NULL,

        #' @field sequencing_reads Character vector of all sequencing reads.
        sequencing_reads = NULL,

        initialize = function(seq_len, read_len, kmer, dbg_kmer, seed, ind, action){
            if(!missing(seq_len)) private$seq_len <- seq_len
            if(!missing(read_len)) private$read_len <- read_len
            if(!missing(kmer)) self$kmer <- kmer
            if(!missing(dbg_kmer)) self$dbg_kmer <- dbg_kmer
            if(!missing(seed)) private$seed <- seed
            if(!missing(action)) private$action <- action
            if(!missing(ind)) private$ind <- ind

            private$get_prob_values()
        },

        #' @description
        #' Load reference sequence from which we will sample
        #' @return None.
        sample_ref_genome = function(){
            # progress message
            t1 <- Sys.time()
            cur.msg <- "Sampling masked hg38 reference genome"
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))

            genome <- BSgenome.Hsapiens.UCSC.hg38.masked::BSgenome.Hsapiens.UCSC.hg38.masked
            genome.attr <- attr(genome, "seqinfo")
            genome.lens <- attr(genome.attr, "seqlengths")[1:22]

            # sample 1k sequences
            N.sample <- 3000
            which.chr <- sample(x = 1:22, size = N.sample, replace = TRUE)
            sample.dt <- data.table(seqnames = 1:22, len = genome.lens)
            sample.dt <- sample.dt[match(which.chr, sample.dt$seqnames)]
            start.pos <- sapply(1:nrow(sample.dt), function(x){
                sample(x = 1:(sample.dt$len[x]-private$seq_len), size = 1)
            })
            sample.dt[, `:=`(
                start = start.pos, 
                width = private$seq_len,
                len = NULL
            )]
            setorder(sample.dt, seqnames)
            sample.dt[, seqnames := paste0("chr", seqnames)]
            seq.names <- paste0(sample.dt$seqnames, "_", sample.dt$start)

            # extract sequences
            sample.dt <- plyranges::as_granges(sample.dt)
            ref.sequences <- Biostrings::getSeq(genome, sample.dt)
            names(ref.sequences) <- seq.names
            ref.sequences <- unique(ref.sequences)

            # filter results
            any.ns <- gregexpr(pattern = "N", text = ref.sequences)
            to.keep <- which(sapply(any.ns, `[`, 1) == -1)
            ref.sequences <- ref.sequences[to.keep]

            # randomly select 1k sequences
            ref.sequences <- sample(ref.sequences, size = 1000, replace = FALSE)
            
            # save results
            dir.create(
                path = "../data/refseq/",
                showWarnings = FALSE
            )
            Biostrings::writeXStringSet(
                x = ref.sequences, 
                filepath = paste0(
                    "../data/refseq/SampledRefGenome",
                    "_SeqLen-", private$seq_len, 
                    "_SeqSeed-", private$seed,
                    ".fasta"
                ), 
                format = "fasta"
            )

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")
        },

        #' @description
        #' Generates random reference sequence and simulates the sonication process
        #' to produce sequencing reads biased by the intrinsic breakage probabilities.
        #' @return None.
        get_reads = function(){
            private$extract_genome_seq()
            private$generate_sequencing_reads()
        }
    ),
    private = list(
        #' @field seq_len Numeric vector. Length of randomly generated reference sequence.
        seq_len = NULL,

        #' @field read_len Numeric vector. Length of sequencing reads.
        read_len = 100,

        #' @field action Character vector. Breakage bias using information from probability ratios
        #'  or enrichment/depletion z-scores. Choose from c("ratio", "zscore").
        action = "ratio",

        #' @field cols_to_keep Character vector. Columns to keep from df_prob.
        cols_to_keep = NULL,

        #' @field seed Numeric vector. Seed used for reproducibility. Set outside of the R6 environment.
        seed = NULL,

        #' @field genome BSgenome object of the human reference genome hg38.
        genome = NULL,

        #' @field ind Numeric vector of the subsetted reference sequence.
        ind = NULL,

        #' @description
        #' Get probability values for ultrasonication experiments.
        #' @return None.
        get_prob_values = function(){
            df.prob <- fread(
                file = paste0("../data/QueryTable/", 
                              "QueryTable_kmer-", 
                              self$kmer, ".csv"),
                showProgress = FALSE
            )
            if(private$action == "zscore"){
                df.prob[, zscore := (
                    zscore-min(zscore, na.rm = TRUE)
                )/(max(zscore, na.rm = TRUE)-min(zscore, na.rm = TRUE)
                )]
            }
            private$cols_to_keep <- to.keep <- c("kmer", private$action)
            self$df_prob <- df.prob[, ..to.keep]
        },

        #' @description
        #' Load in all reference sequences that was sampled from hg38 masked version.
        #' @return None.
        load_all_ref_seq = function(){
            private$genome <- Biostrings::readDNAStringSet(
                filepath = paste0(
                    "../data/refseq/SampledRefGenome",
                    "_SeqLen-", private$seq_len, 
                    "_SeqSeed-", private$seed,
                    ".fasta"
                )
            )
        },

        #' @description
        #' Generate random sequence for use as the reference genome.
        #' @return None.
        extract_genome_seq = function(){
            # progress message
            t1 <- Sys.time()
            cur.msg <- "Extracting reference genome"
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))
            
            self$genome_seq <- paste(private$genome[private$ind])
            self$genome_seq <- strsplit(self$genome_seq, split = "")[[1]]

            # extract base contents
            all.letters <- Biostrings::letterFrequency(
                private$genome[private$ind], 
                letters = "ACGT", 
                OR = 0
            )
            all.letters <- all.letters/width(private$genome[private$ind])
            self$dbg_summary$base_composition <- signif(all.letters[1,], digits = 4)

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")
        },

        #' @description
        #' Generate sequencing reads by chopping up the randomly generated sequence
        #' biased by the intrinsic probabilities of breaking.
        #' @return None.
        generate_sequencing_reads = function(){
            # progress message
            t1 <- Sys.time()
            cur.msg <- "Generating sequencing reads"
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))

            # # Average probability based on non-N dummy bases
            # # add dummy base to start of the random sequence
            # start.of.genome <- genome_seq[1:(kmer-1)]
            # dummy.bases.start <- c(rep("N", (kmer-1)), start.of.genome)
            # dummy.bases.start <- sapply(1:(kmer-1), function(x){
            #     paste0(dummy.bases.start[x:(kmer+x-1)], collapse = "")
            # }, USE.NAMES = FALSE)
            # sapply(1:length(dummy.bases.start), function(x){
            #     kmer.to.find <- stringr::str_remove_all(
            #         string = dummy.bases.start[x],
            #         pattern = "N"
            #     )
            #     similar.kmers <- stringr::str_locate(
            #         string = df.prob$kmer, 
            #         pattern = kmer.to.find
            #     )
            #     kmer.matches <- which(
            #         similar.kmers[, "start"] == (kmer-nchar(kmer.to.find)+1)
            #     )
            #     avg.prob <- as.numeric(colMeans(
            #         df.prob[kmer.matches, ..action], 
            #         na.rm = TRUE
            #     ))
            #     return(avg.prob)
            # })
            self$genome_seq <- c(
                rep("N", (self$kmer/2-1)), 
                self$genome_seq, 
                rep("N", (self$kmer/2-1))
            )

            # add dummy k-mers from each end of randSeq onto k-mer vector
            dummy.base.start <- sapply(1:(self$kmer-1), function(x){
                paste0(self$genome_seq[x:(self$kmer+x-1)], collapse = "")
            }, USE.NAMES = FALSE)
            dummy.base.end <- tail(self$genome_seq, n = (2*self$kmer-2))
            dummy.base.end <- sapply(1:(self$kmer-1), function(x){
                paste0(dummy.base.end[x:(self$kmer+x-1)], collapse = "")
            }, USE.NAMES = FALSE)
            dummy.bases <- c(dummy.base.start, dummy.base.end)

            # assign avg prob to dummy k-mers
            prob <- c(
                self$df_prob[, -c("kmer")],
                rep(colMeans(
                    self$df_prob[, -c("kmer")], 
                    na.rm = TRUE
                ),length(dummy.bases))
            )
            prob <- unlist(prob, use.names = FALSE)
            prob.norm <- prob/sum(prob, na.rm = TRUE)

            dummy.bases.df <- data.table(
                kmer = dummy.bases,
                action = tail(
                    prob, 
                    n = length(dummy.bases)
                )
            )
            setnames(dummy.bases.df, private$cols_to_keep)
            self$df_prob <- rbind(self$df_prob, dummy.bases.df)
            self$df_prob[, prob.norm := prob.norm]
            to.keep <- c("kmer", "prob.norm")
            self$df_prob <- self$df_prob[, ..to.keep]
            setnames(self$df_prob, c("kmer", "prob"))

            # obtain all k-mers of the randomly generated string
            end <- length(self$genome_seq)-self$kmer+1
            genome_seq_kmers <- substring(
                paste(self$genome_seq, collapse = ""),
                first = 1:end, 
                last = (1:end)+self$kmer-1
            )

            # initialise data frame with all k-mers from ref.seq 
            # and associated probs
            kmer.from.seq <- data.table(
                kmer = genome_seq_kmers,
                prob = NA
            )

            # replace NAs with probability values from reference data frame
            fwd.ind <- match(kmer.from.seq$kmer, self$df_prob$kmer)

            # reverse complement kmer index
            rev.comp.kmer.ind <- which(is.na(fwd.ind))
            rev.comp.prob.ind <- match(
                paste(Biostrings::reverseComplement(
                    Biostrings::DNAStringSet(
                        kmer.from.seq$kmer[rev.comp.kmer.ind]
                    )
                )),
                self$df_prob$kmer
            )
            kmer.from.seq$prob[rev.comp.kmer.ind] <- 
                self$df_prob$prob[rev.comp.prob.ind]

            # fwd kmer index
            fwd.prob.ind <- fwd.ind[which(!is.na(fwd.ind))]
            fwd.kmer.ind <- which(!is.na(fwd.ind))
            kmer.from.seq$prob[fwd.kmer.ind] <- self$df_prob$prob[fwd.prob.ind]

            # plot breakage probability over the full genome sequence
            mytitle <- "Kmeric breakage probabilities of genome sequence"
            mysubtitle <- paste0(
                "Kmer: ", self$kmer, ". ",
                "Sequence length: ", length(self$genome_seq), ". ",
                "Sequence set.seed: ", private$seed
            )
            dir.create(
                path = paste0("../figures/exp_", private$ind, "/"),
                showWarnings = FALSE,
                recursive = TRUE
            )
            pdf(
                file = paste0(
                    "../figures/exp_", private$ind,
                    "/BreakProb-vs-GenomeSeq", 
                    "_SeqLen-", private$seq_len, 
                    "_SeqSeed-", private$seed,                    
                    "_kmer-", self$kmer, 
                    ".pdf"
                ),
                width = 11, 
                height = 8
            )
            barplot(
                kmer.from.seq$prob,
                main = "",
                xlab = paste0("Genome sequence (", 
                              self$kmer, "-mer sliding window ", 
                              "by 1 nt)"),
                ylab = "Breakage probabilitiy",
                col = "grey"
            )
            mtext(line=2.2, at=-0.07, adj=0, cex=1.1, mytitle)
            mtext(line=1, at=-0.07, adj=0, cex=0.9, mysubtitle)
            plot.saved <- dev.off()

            # sample multiple breakpoint positions
            len.sampling <- private$seq_len*10
            bp.start.pos <- sample(
                x = nrow(kmer.from.seq),
                size = len.sampling, 
                replace = TRUE,
                prob = kmer.from.seq$prob
            )
            # size selection / discard out-of-bounce reads
            ind.to.keep <- which(
                (bp.start.pos+private$read_len-1) <= length(self$genome_seq)
            )
            bp.start.pos <- bp.start.pos[ind.to.keep]

            # plot breakage probability over the full genome sequence
            mytitle <- paste0(
                "Probability-weighted sampling of genomic sequence ", 
                "positions to replicate ultrasonication process"
            )
            pdf(
                file = paste0(
                    "../figures/exp_", private$ind,
                    "/UltrasonicatedGenome", 
                    "_SeqLen-", private$seq_len, 
                    "_SeqSeed-", private$seed,                    
                    "_kmer-", self$kmer, 
                    ".pdf"
                ),
                width = 11, 
                height = 8
            )
            hist(
                bp.start.pos, 
                breaks = 300,
                main = "",
                xlab = "Genomic sequence position",
                xlim = c(0, private$seq_len)
            )
            mtext(line=2.2, at=-0.07, adj=0, cex=1.1, mytitle)
            mtext(line=1, at=-0.07, adj=0, cex=0.9, mysubtitle)
            plot.saved <- dev.off()

            # synthetic sequencing by synthesis for paired-end sequencing
            # read two information will be used for the breakage scoring and 
            # overall evaluation of the de novo assembly results
            read.one <- IRanges(start = bp.start.pos, width = private$read_len)

            # extract nucleotide bases
            bio.ref.seq <- Biostrings::DNAStringSet(
                paste(self$genome_seq, collapse = "")
            )            
            # get sequencing read for read one
            self$sequencing_reads$read_one <- read.one <- substring(
                text = bio.ref.seq,
                first = start(read.one),
                last = end(read.one)
            )
            read.length <- nchar(read.one)
            self$dbg_summary$coverage <- signif(
                (length(read.one)*private$read_len)/length(self$genome_seq),
                digits = 3
            )
            self$dbg_summary$nr_of_reads <- length(read.one)

            # # remove dummy bases on each end of the randomly generated sequence
            # self$genome_seq <- self$genome_seq[
            #     self$kmer:(length(self$genome_seq)-self$kmer+1)]
            # self$df_prob <- head(self$df_prob, n = -(2*self$kmer-2))

            # save generated reads into text file for use in genome assembly
            dir.create(
                path = paste0("../data/reads/exp_", private$ind, "/"),
                showWarnings = FALSE,
                recursive = TRUE
            )
            file.reads <- file(paste0(
                "../data/reads/exp_", private$ind, 
                "/read_1",
                "_SeqLen-", private$seq_len, 
                "_SeqSeed-", private$seed,
                "_ReadLen-", private$read_len,
                ".txt"
            ))
            writeLines(
                read.one,
                file.reads, 
                sep = "\n"
            )
            close(file.reads)

            # fasta file for read one
            read.one <- Biostrings::DNAStringSet(read.one)
            names(read.one) <- paste0("seq-", 1:length(read.one), ":", private$read_len)
            Biostrings::writeXStringSet(
                x = read.one, 
                filepath = paste0(
                    "../data/reads/exp_", private$ind, 
                    "/read_1",
                    "_SeqLen-", private$seq_len, 
                    "_SeqSeed-", private$seed,
                    "_ReadLen-", private$read_len,
                    ".fasta"
                ),
                format = "fasta"
            )

            # fasta file for read one
            read.two <- Biostrings::DNAStringSet(
                Biostrings::reverseComplement(read.one)
            )
            names(read.two) <- paste0("seq-", 1:length(read.two), ":", private$read_len)
            Biostrings::writeXStringSet(
                x = read.two, 
                filepath = paste0(
                    "../data/reads/exp_", private$ind, 
                    "/read_2",
                    "_SeqLen-", private$seq_len, 
                    "_SeqSeed-", private$seed,
                    "_ReadLen-", private$read_len,
                    ".fasta"
                ),
                format = "fasta"
            )

            genome.seq <- self$dbg_summary$genome_seq <- paste(self$genome_seq, collapse = "")
            genome.seq <- Biostrings::DNAStringSet(genome.seq)
            names(genome.seq) <- "seq-1"
            Biostrings::writeXStringSet(
                x = genome.seq, 
                filepath = paste0(
                    "../data/reads/exp_", private$ind, 
                    "/ref",
                    "_SeqLen-", private$seq_len, 
                    "_SeqSeed-", private$seed,
                    "_ReadLen-", private$read_len,
                    ".txt"
                ),
                format = "fasta"
            )

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")
        }
    )
)