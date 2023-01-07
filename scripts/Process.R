# load dependencies
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(stringr)))
suppressPackageStartupMessages(suppressWarnings(library(data.table)))
suppressPackageStartupMessages(suppressWarnings(library(pbapply)))
suppressPackageStartupMessages(suppressWarnings(library(Biostrings)))
pbo <- pbapply::pboptions(type = "txt", char = "=")

setwd("/Users/paddy/Documents/DPhil/github_repos/GenomeAssembler_dev/lib/")
source("../lib/GenerateReads.R")
source("../lib/DBG.R")
assembler <- DBG$new(
    seq_len = 2000,
    read_len = 10,
    G_cont = 0.20,
    C_cont = 0.15,
    A_cont = 0.20,
    kmer = 8,
    dbg_kmer = 9,
    seed = 1,
    action = "ratio",
    uniform_prob = FALSE
)
assembler$run_assembler()

assembler$results
assembler$results[1,]
paste(assembler$genome_seq, collapse = "")
# assembler$results[denovo.len > 1000]

# assembler$results[, bp.score.norm := bp.score/kmer.break]
# setorder(assembler$results, -bp.score.norm)
# assembler$results

private=self=NULL
private$seq_len = 2000
private$read_len = 10
private$G_cont = 0.20
private$C_cont = 0.15
private$A_cont = 0.20
self$kmer = 8
self$dbg_kmer = 9
private$seed = 1
private$action = "ratio"

# seq_reads <- GenerateReads$new(
#     seq_len = 1000,
#     read_len = 10,
#     G_cont = 0.25,
#     C_cont = 0.10,
#     A_cont = 0.25,
#     kmer = 6,
#     dbg_kmer = 11,
#     seed = 1,
#     action = "ratio",
#     uniform_prob = FALSE
# )
# seq_reads$get_reads()
# seq_reads$dbg_summary

####################################################################################################################################################################################
#' Assembly programmes
#' @spades
#' spades.py -s ../reads/toassemble.fasta --isolate -k 7 --threads 1 -o ./spades --only-assembler 

#' @minia
#' ../../tests/minia/bin/minia -in ../reads/toassemble.fasta -kmer-size 11 -no-bulge-removal -no-tip-removal -no-ec-removal -nb-cores 1 -nb-glue-partitions 200 -out minia/minia_assembly
#' ../../tests/minia/bin/minia -in ../reads/toassemble.fasta -kmer-size 11 -nb-cores 1 -nb-glue-partitions 200 -out minia/minia_assembly
####################################################################################################################################################################################