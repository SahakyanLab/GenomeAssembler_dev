# load dependencies
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(data.table)))
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))
suppressPackageStartupMessages(suppressWarnings(library(pbapply)))
suppressPackageStartupMessages(suppressWarnings(library(Biostrings)))
pbapply::pboptions(char = "=")

args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])
setwd(my.path)

seed = 1234
break_kmer = 8
total_iters = 200
industry_standard = TRUE
coverage_target = 40
seq_len = 50000

loop_grid <- tibble(
    read_lens = c(12,14,16,18,20,25,40),
    all_dbg_kmers = c(11,13,13,15,17,19,37)
)

##################################################################################
#' To demonstrate:
#' 1. Is there a GC content dependency?
##################################################################################
base_dir <- paste0("IndustryModel_", industry_standard)

# load in the genome
ref_genome_file <- paste0(
    "../data/refseq/SampledRefGenome", 
    "_SeqLen-", format(seq_len, scientific = FALSE), 
    "_SeqSeed-", seed,
    ".fasta"
)
refseq <- Biostrings::readDNAStringSet(filepath = ref_genome_file)
refseq <- refseq[1:total_iters]

# count GC content
letter_counts <- Biostrings::letterFrequency(
    refseq, letters = "ACGT", OR = 0
)
letter_counts_norm <- letter_counts / width(refseq)
gc_content <- letter_counts_norm[, "G"]+letter_counts_norm[, "C"]

all_res <- pbapply::pblapply(1:nrow(loop_grid), function(x){
    files <- list.files(
        path = "../data/results",
        pattern = paste0(
            "SolutionsTable_SeqLen-", format(seq_len, scientific = FALSE), 
            "_SeqSeed-", seed, 
            "_ReadLen-", loop_grid$read_lens[x], 
            "_DBGKmer-", loop_grid$all_dbg_kmers[x],  
            "_kmer-", break_kmer, 
            "_IndustryModel-", industry_standard,
            ".csv"
        ),
        full.names = TRUE,
        recursive = TRUE
    )
    files <- stringr::str_sort(files, numeric = TRUE)
    files <- files[1:total_iters]

    df <- lapply(1:length(files), function(i){
        df <- fread(
            files[i], 
            select = c(
                "sequence_len",
                "kmer_breaks",
                "bp_score_norm_by_break_freqs_true",
                "bp_score_norm_by_len_true",
                "bp_score_true",
                "lev_dist_vs_true",
                "stat_test_KS_true"
            ),
            showProgress = FALSE
        )
        df[, gc_content := gc_content[i]]
        return(colMeans(df, na.rm = TRUE))
    })
    df <- do.call(rbind, df)
    df <- as.data.table(df)
    df[, `:=`(
        read_len = loop_grid$read_lens[x],
        dbg_kmer = loop_grid$all_dbg_kmers[x]
    )]
    return(df)
})

df_all_res <- rbindlist(all_res)
unique_read_len <- unique(df_all_res$read_len)
unique_read_len <- paste0("Read len: ", sort(unique_read_len, decreasing = TRUE))

df_all_res <- as_tibble(df_all_res) %>% 
    dplyr::rename(bp_score_norm_by_break_freqs = bp_score_norm_by_break_freqs_true) %>% 
    dplyr::rename(bp_score_norm_by_len = bp_score_norm_by_len_true) %>% 
    dplyr::rename(bp_score = bp_score_true) %>% 
    dplyr::rename(stat_test_KS = stat_test_KS_true) %>% 
    dplyr::mutate(
        read_len = paste0("Read len: ", read_len),
        read_len = factor(read_len, levels = unique_read_len)
    ) %>% 
    tidyr::drop_na()

# GC content vs. breakage scores
colfun = colorRampPalette(c(
    "white","blue","skyblue",
    "chartreuse3","green","yellow",
    "orange","red","darkred"
))
y <- df_all_res$gc_content
nbin <- 150
fac <- 70
size <- 2
ylim <- c(0.3,0.6)

pdf(
    file = paste0(
            "../figures/", base_dir, "/",
            "GC_content-vs-Breakscore",
            "_IndustryModel-", industry_standard, 
            ".pdf"
    ),
    height=10, width=14
)
par(mfrow=c(2,2), mar=c(5, 5, 1, 1), bty="L")

# first plot
x <- df_all_res$bp_score
smoothScatter(
    x=x, y=y, 
    nrpoints=50,
    nbin=nbin,
    bandwidth=c(diff(range(x))/fac, diff(range(y))/fac),
    xlab="Breakage score", ylab="",
    ylim=ylim,
    colramp=colfun,
    cex.lab=size, cex.axis=size, cex.main=size
)

# second plot
x <- df_all_res$stat_test_KS
smoothScatter(
    x=x, y=y, 
    nrpoints=50,
    nbin=nbin,
    bandwidth=c(diff(range(x))/fac, diff(range(y))/fac),
    xlab="KS statistic", ylab="",
    ylim=ylim,
    colramp=colfun,
    cex.lab=size, cex.axis=size, cex.main=size
)

# third plot
x <- df_all_res$bp_score_norm_by_break_freqs * 1e4
smoothScatter(
    x=x, y=y, 
    nrpoints=50,
    nbin=nbin,
    bandwidth=c(diff(range(x))/fac, diff(range(y))/fac),
    xlab = expression("Breakage score norm. by break freq., x10"^-4*""),
    ylab="GC content",
    ylim=ylim,
    colramp=colfun,
    cex.lab=size, cex.axis=size, cex.main=size
)

# third plot
x <- df_all_res$bp_score_norm_by_len * 1e3
smoothScatter(
    x=x, y=y, 
    nrpoints=50,
    nbin=nbin,
    bandwidth=c(diff(range(x))/fac, diff(range(y))/fac),
    xlab = expression("Breakage score norm. by seq. length, x10"^3*""),
    ylab="",
    ylim=ylim,
    colramp=colfun,
    cex.lab=size, cex.axis=size, cex.main=size
)

plot.saved <- dev.off()