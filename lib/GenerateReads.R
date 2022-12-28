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

        initialize = function(seq_len, read_len, G_cont, C_cont, 
                              A_cont, kmer, seed, action){
            if(!missing(seq_len)) private$seq_len <- seq_len
            if(!missing(read_len)) private$read_len <- read_len
            if(!missing(G_cont)) private$G_cont <- G_cont
            if(!missing(C_cont)) private$C_cont <- C_cont
            if(!missing(A_cont)) private$A_cont <- A_cont
            if(!missing(kmer)) self$kmer <- kmer
            if(!missing(seed)) private$seed <- seed
            if(!missing(action)) private$action <- action

            private$get_prob_values()
        },

        #' @description
        #' Generates random reference sequence and simulates the sonication process
        #' to produce sequencing reads biased by the intrinsic breakage probabilities.
        #' @return None.
        get_reads = function(){
            private$generate_genome_seq()
            private$generate_sequencing_reads()
        }
    ),
    private = list(
        #' @field seq_len Numeric vector. Length of randomly generated reference sequence.
        seq_len = 1000,

        #' @field read_len Numeric vector. Length of sequencing reads.
        read_len = 100,

        #' @field G_cont Numeric vector. Amount of guanine content in seq_len.
        G_cont = 0.25,

        #' @field C_cont Numeric vector. Amount of cytosine content in seq_len.
        C_cont = 0.25,

        #' @field A_cont Numeric vector. Amount of adenine content in seq_len.
        A_cont = 0.1,

        #' @field seed Numeric vector. Set the seed for reproducibility.
        seed = 1234,

        #' @field action Character vector. Breakage bias using information from probability ratios
        #'  or enrichment/depletion z-scores. Choose from c("ratio", "zscore").
        action = "ratio",

        #' @field cols_to_keep Character vector. Columns to keep from df_prob.
        cols_to_keep = character(),

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
            private$cols_to_keep <- to.keep <- c("kmer", private$action)
            self$df_prob <- df.prob[, ..to.keep]
        },

        #' @description
        #' Generate random sequence for use as the reference genome.
        #' @return None.
        generate_genome_seq = function(){
            # progress message
            t1 <- Sys.time()
            cur.msg <- "Generating random sequence"
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))

            set.seed(private$seed)
            private$G_cont <- round(private$G_cont*private$seq_len)
            private$C_cont <- round(private$C_cont*private$seq_len)
            private$A_cont <- round(private$A_cont*private$seq_len)
            T_cont <- private$seq_len-(private$G_cont+private$C_cont+private$A_cont)
            
            # Sum of base content not equal to seq_len of genome
            base::stopifnot(identical(
                private$seq_len, 
                (private$G_cont+private$C_cont+private$A_cont+T_cont)
            ))

            genome_seq <- c(
                rep("G", private$G_cont), 
                rep("C", private$C_cont),
                rep("A", private$A_cont), 
                rep("T", T_cont)
            )
            
            # random sampling of base positions in the sequence
            self$genome_seq <- sample(
                x = genome_seq, 
                size = private$seq_len, 
                replace = FALSE
            )

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
                rep("N", (self$kmer-1)), 
                self$genome_seq, 
                rep("N", (self$kmer-1))
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
                )/length(dummy.bases), 
                length(dummy.bases))
            )
            prob <- unlist(prob, use.names = FALSE)

            dummy.bases.df <- data.table(
                kmer = dummy.bases,
                action = tail(
                    prob, 
                    n = length(dummy.bases)
                )
            )
            setnames(dummy.bases.df, private$cols_to_keep)
            self$df_prob <- rbind(self$df_prob, dummy.bases.df) 
            setnames(self$df_prob, c("kmer", "prob"))

            # obtain all k-mers of the randomly generated string
            end  <- length(self$genome_seq)-self$kmer+1
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

            # sample 10k breakpoint positions at a time
            len.breakpoint.positions <- private$seq_len*1000
            len.sampling <- 10000
            breakpoint.positions <- replicate(
                n = len.breakpoint.positions/len.sampling, {
                sampling <- sample(
                    x = nrow(kmer.from.seq), 
                    size = len.sampling,
                    replace = TRUE, 
                    prob = kmer.from.seq$prob
                )
            }, simplify = TRUE)

            # generate all reads
            three.sd    <- (private$read_len*0.1)*3 # 1.SD = 10% of read length
            upper.limit <- private$read_len+three.sd
            lower.limit <- private$read_len-three.sd
            max.runs    <- len.sampling
            end         <- length(self$genome_seq)-1

            sampling.points <- sample(
                x = breakpoint.positions, 
                size = max.runs, 
                replace = FALSE
            )
            start <- seq(from = 1, to = max.runs, by = 2)
            end <- seq(from = 2, to = max.runs, by = 2)
            mat <- matrix(c(sampling.points[start], sampling.points[end]), nrow = 2)
            mat <- apply(mat, 2, sort)
            mat.ir <- IRanges(start = mat[1,]+1, end = mat[2,])
            bio.ref.seq <- Biostrings::DNAStringSet(paste(self$genome_seq, collapse = ""))
            read <- substring(
                text = bio.ref.seq,
                first = start(mat.ir),
                last = end(mat.ir)
            )
            read.length <- nchar(read)
            match.read <- which(read.length >= (lower.limit) & read.length <= (upper.limit))
            reads.strings <- read[match.read]
            reads.freq <- read.length[match.read]
            coverage <- (length(reads.strings)*private$read_len)/length(self$genome_seq)

            # remove dummy bases on each end of the randomly generated sequence
            self$genome_seq <- self$genome_seq[
                self$kmer:(length(self$genome_seq)-self$kmer+1)]
            self$df_prob <- head(self$df_prob, n = -(2*self$kmer-2))

            # save generated reads into text file for use in genome assembly
            dir.create(
                path = "../data/reads/",
                showWarnings = FALSE,
                recursive = TRUE
            )
            file.reads <- file("../data/reads/reads.txt")
            writeLines(
                unlist(reads.strings), 
                file.reads, 
                sep = "\n"
            )
            close(file.reads)

            # fasta file
            all.reads <- unlist(reads.strings)
            all.reads <- Biostrings::DNAStringSet(all.reads)
            names(all.reads) <- paste0("seq-", 1:length(all.reads))
            Biostrings::writeXStringSet(
                x = all.reads, 
                filepath = "../data/reads/toassemble.fasta", 
                format = "fasta"
            )

            rand.seq.combined <- paste(self$genome_seq, collapse = "")
            rand.seq.combined <- Biostrings::DNAStringSet(rand.seq.combined)
            names(rand.seq.combined) <- "seq-1"
            Biostrings::writeXStringSet(
                x = rand.seq.combined, 
                filepath = "../data/reads/ref.fasta", 
                format = "fasta"
            )

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")
        }
    )
)