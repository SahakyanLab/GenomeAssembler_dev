# load dependencies
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
# suppressPackageStartupMessages(suppressWarnings(library(stringr)))
suppressPackageStartupMessages(suppressWarnings(library(data.table)))
suppressPackageStartupMessages(suppressWarnings(library(pbapply)))
suppressPackageStartupMessages(suppressWarnings(library(Biostrings)))
pbo <- pbapply::pboptions(type = "txt", char = "=")

setwd("/Users/paddy/Documents/DPhil/github_repos/GenomeAssembler_dev/lib/")
source("../lib/GenerateReads.R")
source("../lib/DBG.R")

assembler <- DBG$new(
    seq_len = 1000,
    read_len = 50,
    G_cont = 0.25,
    C_cont = 0.25,
    A_cont = 0.1,
    kmer = 8,
    dbg_kmer = 11,
    seed = 1,
    action = "ratio"
)
assembler$run_assembler()

print(dim(assembler$results))
paste(assembler$genome_seq, collapse = "")
head(assembler$results, n=5)

setorder(assembler$results, -bp.score)
head(assembler$results, n=5)

# private=self=NULL
# private$seq_len = 100
# private$read_len = 50
# private$G_cont = 0.25
# private$C_cont = 0.25
# private$A_cont = 0.1
# self$kmer = 6
# self$dbg_kmer = 11
# private$seed = 1
# private$action = "ratio"