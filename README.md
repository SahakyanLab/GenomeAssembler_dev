# De Novo Genome Assembly
This repository is a toy model for assembling a _de novo_ genome sequence by taking
into account the breakage probabilities associated with the randomly generated
reference sequence. Below are the sequential steps of this model.

## Random Sequence Assembler
A genome sequence is assembled based on the four DNA bases of A, T, G and C. The
parameters are the four DNA base contents and the desired length of the sequence.

## Reads Generator
This function takes the randomly generated genome and generates reads of varying sizes
to mimic the DNA sonication process. The workflow is the following:
1. Assign uniform or non-uniform breakage probabilities to a k-mer length across the whole Genome (default is k-mer of 2)
2. Randomly sample 10k breakpoint positions at a time
3. Generate reads biased based on these probabilities. The accepted reads will lie within 3 standard deviations from the input read length.

## De Bruijn Graph Assembler
This function assembles the breakage probability-biased reads together into a fully
assembled genome using the De Bruijn Graph Assembly algorithm approach. The k-mers are weighted in the graph. Sources for the weighted k-mer approach are found
[here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4719071/) and
[here](https://www.cs.jhu.edu/~langmea/resources/lecture_notes/assembly_dbg.pdf).

The workflow of the assembler is the following:

1. Obtain weighted k-mers from all the reads
2. Create a dictionary of prefix and suffix from each k-mer to obtain a de bruijn Graph
3. Count the number of in-degree and out-degree to check if the Eulerian path is balanced
4. Traverse along all branches of the rooted tree and save all de novo assembled genomes

## Alignment Scorings
Two functions exist: (1) pairwise alignment scoring between each de novo assembled genome
and the reference sequence which is the randomly generated sequence we generated at the start.
(2) alignment scoring of each de novo assembled genome based on breakage probabilities
and frequency of overlaps with reads.

### Pairwise Alignment Scorings
This function aligns each generated de novo genome with the reference sequence and
calculates a variety of alignment scores, including higher-order k-meric signature,
levenshtein distance, local-global, global-local, local, global and a final score that sums all results together.

### Breakpoint Alignment Scorings
This function aligns each generated de novo genome with all the reads and calculates
an alignment score based on the generated breakage probabilities and the frequency
of occurrence/overlap with that specific de novo assembled genome sequence.

## Assembler
This function takes all the required scripts together to assemble a genome. The workflow
is the following:

1. Generate random sequence                                        
2. Assign uniform or non-uniform breakage probabilities to a k-mer length across the whole genome
3. Generate reads biased based on these probabilities                         
4. De novo assemble genome based on generated reads                           
5. Obtain pairwise alignment scores for each de novo generated assembly against reference sequence                                                         
6. Obtain alignment scores for each de novo generated assembly ONLY taking into account the breakage probabilities and alignment against reads

## Metrics Analysis
This function runs the assembler function and then runs calculations on the
alignment scoring metrics. It calculates the following:

1. Correlation coefficients of each pairwise alignment scoring metric versus the breakage probability alignment scores
2. The percent of correctness of breakpoints versus each alignment metric, where we defined the "correctness" as a breakpoint score that is equal to or greater than the pairwise alignment score
3. The breakpoint alignment scoring performance by measuring the average levenshtein distance for various read lengths
4. Multi-regression analysis of local versus global versus the other pairwise alignment scoring metrics
