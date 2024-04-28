#!/bin/bash

pwd="$(pwd)/"

# use in-house simple DBG model to run de novo genome assembly programme
Rscript 02_Real_vs_rand_prob_own.R $pwd

# no relationship between k-mer frequency and breakage probabilities
Rscript 01_Real_vs_rand_prob_break_vs_kmers.R $pwd

# no dependency on GC content
Rscript 03_GC_content_dependency.R $pwd

# Uncomment the below line to optionally run using velvet de novo assembler
# Rscript 00_Real_vs_rand_prob_velvet.R $pwd