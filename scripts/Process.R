# load dependencies
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(stringr)))
suppressPackageStartupMessages(suppressWarnings(library(data.table)))
suppressPackageStartupMessages(suppressWarnings(library(Biostrings)))

setwd("/Users/paddy/Documents/DPhil/github_repos/GenomeAssembler_dev/lib/")
source("../lib/GenerateReads.R")
source("../lib/DBG.R")

assembler = DBG$new(
    seq_len = 1000,
    read_len = 100,
    G_cont = 0.25,
    C_cont = 0.25,
    A_cont = 0.1,
    kmer = 4,
    dbg_kmer = 9,
    seed = 1234,
    action = "ratio"
)
assembler$run_assembler()

assembler$dbg_summary