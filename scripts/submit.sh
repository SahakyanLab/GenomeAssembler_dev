#!/bin/bash

pwd="$(pwd)/"

# use velvet to run de novo genome assembly programme
Rscript 00_Real_vs_rand_prob_velvet.R $pwd 

# no relationship between k-mer frequency and breakage probabilities
Rscript 01_Real_vs_rand_prob_break_vs_kmers.R $pwd

# no dependency on GC content
Rscript 03_GC_content_dependency.R $pwd

# uncomment the below line to optinally run a simple self-implementation of de novo assemblers
# Rscript 02_Real_vs_rand_prob_own.R $pwd 