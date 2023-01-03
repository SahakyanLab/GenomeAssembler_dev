# load dependencies
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
# suppressPackageStartupMessages(suppressWarnings(library(stringr)))
suppressPackageStartupMessages(suppressWarnings(library(data.table)))
suppressPackageStartupMessages(suppressWarnings(library(pbapply)))
suppressPackageStartupMessages(suppressWarnings(library(Biostrings)))
pbo <- pbapply::pboptions(type = "txt", char = "=")

setwd("/Users/paddy/Documents/DPhil/github_repos/GenomeAssembler_dev/lib/")
source("../lib/GenerateReads.R")
seq_reads <- GenerateReads$new(
    seq_len = 1000,
    read_len = 50,
    G_cont = 0.25,
    C_cont = 0.25,
    A_cont = 0.1,
    kmer = 6,
    seed = 1,
    action = "ratio",
    uniform_prob = FALSE
)
seq_reads$get_reads()

source("../lib/DBG.R")
assembler <- DBG$new(
    seq_len = 1000,
    read_len = 50,
    G_cont = 0.25,
    C_cont = 0.25,
    A_cont = 0.1,
    kmer = 6,
    dbg_kmer = 11,
    seed = 1,
    action = "ratio",
    uniform_prob = FALSE
)
assembler$run_assembler()

print(dim(assembler$results))
paste(assembler$genome_seq, collapse = "")
assembler$results

# assembler$results[, bp.score.norm := bp.score/kmer.break]
# setorder(assembler$results, -bp.score.norm)
# assembler$results

setorder(assembler$results, -denovo.len)
assembler$results

# private=self=NULL
# private$seq_len = 20
# private$read_len = 10
# private$G_cont = 0.25
# private$C_cont = 0.25
# private$A_cont = 0.1
# self$kmer = 4
# self$dbg_kmer = 5
# private$seed = 1
# private$action = "ratio"