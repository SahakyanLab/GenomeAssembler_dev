# load dependencies
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(stringr)))
suppressPackageStartupMessages(suppressWarnings(library(data.table)))
suppressPackageStartupMessages(suppressWarnings(library(pbapply)))
suppressPackageStartupMessages(suppressWarnings(library(Biostrings)))
suppressPackageStartupMessages(suppressWarnings(library(plyranges)))
suppressPackageStartupMessages(suppressWarnings(library(BSgenome.Hsapiens.UCSC.hg38.masked)))
suppressPackageStartupMessages(suppressWarnings(library(Rcpp)))
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))
pbo <- pbapply::pboptions(type = "txt", char = "=")

setwd("/Users/paddy/Documents/DPhil/github_repos/GenomeAssembler_dev/lib/")
source("../lib/GenerateReads.R")
source("../lib/DeNovoAssembler.R")
Rcpp::sourceCpp("../lib/DeNovoAssembler.cpp")

seed = 1234
seq_len = 1000
i = 1

set.seed(seed = seed)
assembler <- DeNovoAssembler$new(
    seq_len = seq_len,
    read_len = 30,
    kmer = 8,
    dbg_kmer = 9,
    seed = seed,
    action = "ratio",
    ind = i,
    reads_only = FALSE,
    cores = 10
)
# assembler$sample_ref_genome()
assembler$run_assembler(bins = 10)

assembler$results

# dim(assembler$results)
# df <- copy(assembler$results)
# df[, bins := cut(df$lev.dist.vs.true, breaks = 10)]
# boxplot(
#     bp.score ~ bins, 
#     data = df, 
#     frame = FALSE,
#     col = "lightblue",
#     border = "black",
#     xlab = "Levenshtein distance vs. true genome",
#     ylab = "Breakage probability score",
#     main = ""
# )

set.seed(seed = 1234)
private=self=NULL
private$seq_len = 1000
private$read_len = 50
self$kmer = 8
self$dbg_kmer = 9
private$action = "ratio"
private$ncpu = 4
private$seed = 1234
private$ind=ind=1

cp ../../data/reads/exp_1/read_1_SeqLen-1000_SeqSeed-1234_ReadLen-20.fasta .
cp ../../data/reads/exp_1/read_2_SeqLen-1000_SeqSeed-1234_ReadLen-20.fasta .
perl shuffleReads_fasta.pl read_1_SeqLen-1000_SeqSeed-1234_ReadLen-20.fasta read_2_SeqLen-1000_SeqSeed-1234_ReadLen-20.fasta input.fasta
velveth Assem 9 -shortPaired -fasta
velvetg Assem
rm read_1_*; rm read_2_*; rm input.fasta

####################################################################################################################################################################################
#' Assembly programmes
#' @spades
#' spades.py -1 ../reads/read_1.fasta -2 ../reads/read_2.fasta -k 11 -o ./spades --only-assembler

#' @minia
#' ../../tests/minia/bin/minia -in ../reads/read_1.fasta -kmer-size 11 -no-bulge-removal -no-tip-removal -no-ec-removal -nb-cores 1 -nb-glue-partitions 200 -out minia/minia_assembly
#' ../../tests/minia/bin/minia -in ../reads/toassemble.fasta -kmer-size 11 -nb-cores 1 -nb-glue-partitions 200 -out minia/minia_assembly
####################################################################################################################################################################################