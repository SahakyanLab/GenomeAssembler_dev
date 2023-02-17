# read arguments from job submission
args <- commandArgs(trailingOnly = TRUE)
# my.path <- as.character(args[1])
read.len <- as.numeric(args[1])
sample.ref <- as.logical(as.character(args[2]))

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

my.path="/Users/paddy/Documents/DPhil/github_repos/GenomeAssembler_dev/scripts/"
setwd(my.path)
source("../lib/GenerateReads.R")
source("../lib/DeNovoAssembler.R")
Rcpp::sourceCpp("../lib/DeNovoAssembler.cpp")

seed = 1234
seq_len = 1000
cores = 10
read.len = 15
sample.ref = FALSE
kmer = 8
total_iters = 100

for(i in 1:total_iters){
    set.seed(seed = seed)
    assembler <- DeNovoAssembler$new(
        seq_len = seq_len,
        read_len = read.len,
        kmer = kmer,
        dbg_kmer = 9,
        seed = seed,
        action = "ratio",
        ind = i,
        cores = cores
    )
    if(sample.ref & i == 1) assembler$sample_ref_genome()
    assembler$run_assembler(total_iters = total_iters)
}