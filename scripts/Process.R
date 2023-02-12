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
    read_len = 20,
    kmer = 8,
    dbg_kmer = 9,
    seed = seed,
    action = "ratio",
    ind = i,
    reads_only = FALSE
)
assembler$sample_ref_genome()
assembler$run_assembler(bins = 10)