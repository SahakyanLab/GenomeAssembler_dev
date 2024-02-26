GenerateReads <- R6::R6Class(
    classname = "GenerateReads",
    public = list(
        #' @field genome_seq Character vector. Randomly generated reference genome.
        genome_seq = NULL,

        #' @field df_prob Data.table of kmers, probability ratios and enrichment
        #'  and depletion z-scores obtained through analysis of 
        #'  ultrasonication experiments.
        df_prob = NULL,

        #' @field kmer_from_seq Data.Table of kmers and probabilities of true solution.
        kmer_from_seq = NULL,

        #' @field kmer_ref Data.table of all kmers and probabilities.
        kmer_ref = NULL,

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

        initialize = function(seq_len, read_len, coverage_target, 
                              kmer, dbg_kmer, seed, ind, action){
            if(!missing(seq_len)) private$seq_len <- seq_len
            if(!missing(read_len)) private$read_len <- read_len
            if(!missing(coverage_target)) private$coverage_target <- coverage_target
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
            cur.msg <- "Sampling genomic positions from the T2T human genome"
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))

            genome <- BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0

            # get the start and end positions of each chromosome
            refseq.table <- as.data.frame(genome@seqinfo)
            refseq.table <- refseq.table[grepl(
                pattern = "^([1-9]|1[0-9]|2[0-2])$", 
                x = rownames(refseq.table)
            ),]
            refseq.table <- data.table(
                seqnames = rownames(refseq.table),
                len = refseq.table$seqlengths
            )

            # sample sequences
            N.sample <- 1000
            sample.dt <- refseq.table[sample(x = .N, size = N.sample, replace = TRUE)]
            sample.dt[, `:=`(
                start = sapply(len, function(l) sample(x = 1:(l-1), size = 1)),
                width = private$seq_len,
                len = NULL
            )]
            setorder(sample.dt, seqnames, start)
            sample.dt[, seqnames := paste0("chr", seqnames)]
            seq.names <- paste0(sample.dt$seqnames, "_", sample.dt$start)

            # extract sequences
            sample.dt <- plyranges::as_granges(sample.dt)

            if(!any(grepl(pattern = "^chr", x = seqnames(genome)))){
                chr_names <- paste0("chr", seqnames(genome))
                seqnames(genome@seqinfo) <- seqnames(genome) <- chr_names
            }
            ref.sequences <- Biostrings::getSeq(genome, sample.dt)
            names(ref.sequences) <- seq.names
            ref.sequences <- unique(ref.sequences)
            
            # save results
            dir.create(
                path = "../data/refseq/",
                showWarnings = FALSE
            )
            Biostrings::writeXStringSet(
                x = ref.sequences, 
                filepath = paste0(
                    "../data/refseq/SampledRefGenome",
                    "_SeqLen-", format(private$seq_len, scientific = FALSE), 
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
            private$load_all_ref_seq()
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

        #' @field coverage_target Numeric vector of the desired sequencing coverage using 
        #' the Lander/Waterman equation.
        coverage_target = 10,

        #' @description
        #' Get probability values for ultrasonication experiments.
        #' @return None.
        get_prob_values = function(){            
            df_probs <- lapply(seq(2,8,2), function(k){
                dat <- fread(
                    file = paste0("../data/QueryTable/QueryTable_kmer-", k, ".csv"),
                    showProgress = FALSE
                )

                # replace NAs with the lowest probabilities in the table
                nas <- which(is.na(dat$prob))
                if(any(nas)){
                    lowest.probs <- min(dat$prob, na.rm = TRUE)
                    dat$prob[nas] <- lowest.probs
                }
                return(dat)
            })
            df_probs <- rbindlist(df_probs)

            # normalise probability vector to sum to one
            df_probs[, `:=`(
                prob = prob / sum(prob, na.rm = TRUE),
                kmer_len = nchar(kmer)
            )]

            self$df_prob <- split(df_probs, df_probs$kmer_len)
            names(self$df_prob) <- paste0("kmer_", seq(2,8,2))
            self$df_prob <- lapply(self$df_prob, function(x) x[, kmer_len := NULL])
            
            self$df_prob$all <- df_probs
            self$df_prob$all[, kmer_len := NULL]
        },

        #' @description
        #' Load in all reference sequences that was sampled from hg38 masked version.
        #' @return None.
        load_all_ref_seq = function(){ 
            ref.genome.file <- paste0(
                "../data/refseq/SampledRefGenome", 
                "_SeqLen-", format(private$seq_len, scientific = FALSE), 
                "_SeqSeed-", private$seed,
                ".fasta"
            )

            # sample from human genome if file doesn't exist 
            if(!file.exists(ref.genome.file)) assembler$sample_ref_genome()

            # load in the genome
            private$genome <- Biostrings::readDNAStringSet(filepath = ref.genome.file)
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
            self$dbg_summary$base_composition <- all.letters/width(private$genome[private$ind])
            # self$dbg_summary$base_composition <- signif(all.letters[1,], digits = 4)

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
            fwd.ind <- match(kmer.from.seq$kmer, self$df_prob[[paste0("kmer_", self$kmer)]]$kmer)
            self$kmer_from_seq <- kmer.from.seq$prob <- self$df_prob[[paste0("kmer_", self$kmer)]]$prob[fwd.ind]

            if(private$plot_results){
                # plot breakage probability over the full genome sequence
                mytitle <- "Kmeric breakage probabilities of genome sequence"
                mysubtitle <- paste0(
                    "Kmer: ", self$kmer, ". ",
                    "Sequence length: ", format(private$seq_len, scientific = FALSE), ". ",
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
                        "_SeqLen-", format(private$seq_len, scientific = FALSE),
                        "_SeqSeed-", private$seed,                    
                        "_kmer-", self$kmer, 
                        "_IndustryModel-", private$industry_standard,
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
            }

            # sample multiple breakpoint positions
            len.sampling <- ceiling(private$coverage_target*length(self$genome_seq)/private$read_len)
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
            if(private$plot_results){
                mytitle <- paste0(
                    "Probability-weighted sampling of genomic sequence ", 
                    "positions to replicate ultrasonication process"
                )
                pdf(
                    file = paste0(
                        "../figures/exp_", private$ind,
                        "/UltrasonicatedGenome", 
                        "_SeqLen-", format(private$seq_len, scientific = FALSE),
                        "_SeqSeed-", private$seed,
                        "_kmer-", self$kmer, 
                        "_IndustryModel-", private$industry_standard,
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
                    xlim = c(0, private$seq_len),
                    density = FALSE
                )
                mtext(line=2.2, at=-0.07, adj=0, cex=1.1, mytitle)
                mtext(line=1, at=-0.07, adj=0, cex=0.9, mysubtitle)
                plot.saved <- dev.off()
            }

            # #' Simulate ultrasonication of genome to generate fragments of varying sizes
            # broken_fragments <- matrix(c(
            #     bp.start.pos[seq(from = 1, to = length(bp.start.pos)-1, 2)],
            #     bp.start.pos[seq(from = 2, to = length(bp.start.pos), 2)]
            # ), ncol = 2)
            # broken_fragments <- t(apply(broken_fragments, 1, sort))
            # colnames(broken_fragments) <- c("start", "end")

            # #' Simulate sequencing by synthesis where the user specifies the number of 
            # #' cycles to add nucleotides for synthesis. Read two is the reverse complement of 
            # #' the 3' of the broken fragment.
            # read.one <- IRanges(
            #     start = broken_fragments[,"start"], 
            #     width = private$read_len
            # )
            # read.two <- IRanges(
            #     start = broken_fragments[,"end"] - private$read_len, 
            #     width = private$read_len
            # )
            
            #' Old version:
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

            # save generated reads into text file for use in genome assembly
            dir.create(
                path = paste0("../data/reads/exp_", private$ind, "/"),
                showWarnings = FALSE,
                recursive = TRUE
            )

            # fasta file for read one
            read.one <- Biostrings::DNAStringSet(read.one)
            refgenome.name <- names(private$genome[private$ind])
            refgenome.chrname <- stringr::str_extract(
                string = refgenome.name,
                pattern = "[^_]+"
            )
            refgenome.start.pos <- stringr::str_extract(
                string = refgenome.name, pattern = "(?<=_).*"
            ) %>% as.numeric()

            refgenome.start.pos.read1 <- paste0(
                refgenome.chrname, "_",
                # #' After:
                # refgenome.start.pos+broken_fragments[,"start"], "_", 
                # refgenome.start.pos+broken_fragments[,"start"]+private$read_len
                
                # ' Before:
                refgenome.start.pos+bp.start.pos, "_", 
                refgenome.start.pos+bp.start.pos+private$read_len
            )
            names(read.one) <- paste0(
                refgenome.start.pos.read1, ":0_",
                1:length(read.one), "/1"
            )
            if(private$save_read_files){
                Biostrings::writeXStringSet(
                    x = read.one, 
                    filepath = paste0(
                        "../data/reads/exp_", private$ind, 
                        "/read_1",
                        "_SeqLen-", format(private$seq_len, scientific = FALSE),
                        "_SeqSeed-", private$seed,
                        "_ReadLen-", private$read_len,
                        "_DBGKmer-", self$dbg_kmer,
                        ".fasta"
                    ),
                    format = "fasta"
                )
            }

            #' Before: 
            # fasta file for read one
            read.two <- Biostrings::DNAStringSet(
                Biostrings::reverseComplement(read.one)
            )
            names(read.two) <- paste0(
                refgenome.start.pos.read1, ":0_",
                1:length(read.two), "/2"
            )
            if(private$save_read_files){
                Biostrings::writeXStringSet(
                    x = read.two, 
                    filepath = paste0(
                        "../data/reads/exp_", private$ind, 
                        "/read_2",
                        "_SeqLen-", format(private$seq_len, scientific = FALSE),
                        "_SeqSeed-", private$seed,
                        "_ReadLen-", private$read_len,
                        "_DBGKmer-", self$dbg_kmer,
                        ".fasta"
                    ),
                    format = "fasta"
                )
            }

            # reference genome in fasta format
            genome.seq <- self$dbg_summary$genome_seq <- paste(self$genome_seq, collapse = "")
            genome.seq <- Biostrings::DNAStringSet(genome.seq)
            names(genome.seq) <- "seq-1"

            if(private$save_read_files){
                Biostrings::writeXStringSet(
                    x = genome.seq, 
                    filepath = paste0(
                        "../data/reads/exp_", private$ind, 
                        "/ref",
                        "_SeqLen-", format(private$seq_len, scientific = FALSE),
                        "_SeqSeed-", private$seed,
                        "_ReadLen-", private$read_len,
                        "_DBGKmer-", self$dbg_kmer,
                        ".fasta"
                    ),
                    format = "fasta"
                )
            }

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")
        }
    )
)
