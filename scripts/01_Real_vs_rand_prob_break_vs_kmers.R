# load dependencies
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(data.table)))
suppressPackageStartupMessages(suppressWarnings(library(Rcpp)))
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))
suppressPackageStartupMessages(suppressWarnings(library(pbapply)))
suppressPackageStartupMessages(suppressWarnings(library(Biostrings)))
suppressPackageStartupMessages(suppressWarnings(library(plyranges)))
suppressPackageStartupMessages(suppressWarnings(library(BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0)))
pbapply::pboptions(char = "=")

args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])
setwd(my.path)

source("../lib/GenerateReads.R")
source("../lib/DeNovoAssembler.R")

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

x <- 1
kmer_tables <- lapply(seq(2,8,2), function(break_kmer){
    assembler <- DeNovoAssembler$new(
        seq_len = seq_len,
        read_len = loop_grid$read_lens[x],
        coverage_target = coverage_target,
        kmer = break_kmer,
        dbg_kmer = loop_grid$all_dbg_kmers[x],
        seed = seed,
        ind = 1,
        only_kmers_from_reads = TRUE,
        save_read_files = FALSE,
        industry_standard = industry_standard,
        plot_results = FALSE
    )
    assembler$run_assembler(total_iters = 1)
    
    df <- assembler$table_read_kmer_prob

    label <- paste0(
        "Kmer: ", break_kmer, "\n",
        "R^2:  ", signif(cor(df$prob, df$count)^2, 2)
    )
    df[, label := label]

    return(df)
})
kmer_tables <- rbindlist(kmer_tables)

# plot kmer counts vs. kmer breakage probabilities
p1 <- as_tibble(kmer_tables) %>%
    dplyr::mutate(label = as.factor(label)) %>% 
    ggplot(aes(x = prob, y = count)) + 
    geom_point() + 
    facet_wrap(vars(label), nrow = 1, scales = "free") + 
    theme_bw() + 
    theme_classic() + 
    theme(
        axis.text.x = element_blank(),
        text = element_text(size = 15)
    ) + 
    labs(
        x = "Breakage probability",
        y = "Count of kmers"
    )

ggsave(
    filename = "../figures/KmerCounts_vs_Prob.png",
    plot = p1,
    height = 4, width = 12
)