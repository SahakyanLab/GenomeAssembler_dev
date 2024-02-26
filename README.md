# De Novo Genome Assembly
This repository is a proof of concept for incorporating k-meric breakage probabilities in evaluating *de novo* genome assembled solutions.

## Setup

Clone the project:

```
git clone https://github.com/SahakyanLab/GenomeAssembler_dev.git
```

Please follow the instructions below on how to acquire the public datasets, setup the directory stucture, and software necessary to run all the studies from the publication.  At the end of this `README` file, you can find two separate bash script commands that runs the majority of the setup and runs the calculations sequentially. 

## Software requirements
The resource-demanding computations were performed on a single NVIDIA RTX A6000 GPU with 40GB RAM. The developed workflows and analyses employed the [R programming language 4.3.2](https://www.r-project.org/).

Please run the below script to install the latest versions of the R packages necessary to perform the calculations and analyses. 

```bash
bash ./setup/install_packages.sh
```

Please also download and install the below software.

### Edlib

* Please clone the repo from [this link](https://github.com/Martinsos/edlib) (Edlib >= 1.2.7). Place the [edlib.h](https://github.com/Martinsos/edlib/tree/master/edlib/include) and [edlib.cpp](https://github.com/Martinsos/edlib/tree/master/edlib/src) into [lib/edlib/](https://github.com/SahakyanLab/GenomeAssembler_dev/tree/master/lib/edlib/) folder.

### kseq.h

Please download the kseq.h file from [this link](https://github.com/attractivechaos/klib). Place this into [lib/](https://github.com/SahakyanLab/GenomeAssembler_dev/tree/master/lib/) folder.

### phmap.hpp via gtl

Please clone the repo from [this link](https://github.com/greg7mdp/gtl). Place the contents of gtl into [lib/](https://github.com/SahakyanLab/GenomeAssembler_dev/tree/master/lib/) folder.

### `velvet` *de novo* genome assembly algorithm

The `velvet` *de novo* genome assembly algorithm was used in this study to assemble the simulated short sequencing reads and is based on a de Bruijn graph. Originally published [here](http://www.genome.org/cgi/doi/10.1101/gr.074492.107), the GitHub repository can be [found here](https://github.com/dzerbino/velvet.git).

Please clone the repo and place the contents into [lib/](https://github.com/SahakyanLab/GenomeAssembler_dev/tree/master/lib/) folder. Then follow the installation manual.

## Workflow

### Samping of genomic sequences

To sample genomic sequences for reference, we took the latest telomere-to-telomere (T2T) human genome assembly version. We sampled 1k unique sequences across all autosomes, each with a 50kb segment. The random sampling is done with a fixed seed of 1234.

### Sequencing reads generator
This function takes the subsampled T2T genome assembly segment and generates sequencing reads of fixed size to mimic the sequencing read generation of NGS platforms that specify the length by the number of cycles. 

### de Bruijn graph assembler

The main study employs the `velvet` *de novo* genome assembly algorithm to assemble the short sequencing reads into contigs. 

If you want to run the optional study as [detailed below](#optional-experiment), the following describes the simple implementation of a de Bruijn graph assembly.

The breakage probability-biased reads are brute-forced assembled together into long scaffolds using the de Bruijn graph assembly algorithm approach. The k-mers are weighted in the graph. Sources for the weighted k-mer approach are found
[here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4719071/) and
[here](https://www.cs.jhu.edu/~langmea/resources/lecture_notes/assembly_dbg.pdf).

The workflow of the assembler is the following:

1. Obtain weighted k-mers from all the sequencing reads
2. Create a dictionary of prefix and suffix from each k-mer to obtain a de bruijn Graph
3. Count the number of in-degree and out-degree to check if the Eulerian path is balanced
4. Traverse along all branches of the rooted tree and save all de novo assembled genomes

### Calculating k-meric breakage probability score

Each simulated sequencing read is aligned against each assembled sequence through brute force alignment. If a match is found, then the broken genomic position is expanded into an octamer sequence, and map it to the probability of breaking the octamer. The breakage score is obtained through the weighted sum of the count and probability of each octamer. 

### Evaluating closeness of assembled sequences with the reference genome

We evaluated the breakage probability score for each assembled sequence. Here, we aligned each read to the assembled sequence, extracted the octameric breakage probability, and counted the breakage frequency of each octamer. The breakage score was obtained through the weighted sum of the count and probability of each octamer. We calculated the Kolmogorov–Smirnov (KS) statistic between their probability distributions to evaluate the closeness of the assembled solution to the 50 kb-long reference sequence. 

To assess the similarities between their genomic sequences, we used the Levenshtein distance metric from the edlib.

### Optional experiment

The study is based on an existing *de novo* genome assembly algorithm `velvet` and evaluates individual contig outputs by its k-meric breakage scores. 

There are two alternative implementations. First and most ideal case, the k-meric breakage score is built directly into the *de novo* genome assembly algorithm, helping direct the traversal along the de Bruijn graph. However, this implementation takes time and was out of scope for this proof of concept demonstration. 

Second and simpler case, we brute-force *de novo* assemble the contigs to generate long scaffolds and these final solutions are evaluated using the k-meric breakage scores. This implementation was faster to implement, though not efficient. This is an optional study and can be run by uncommenting the last line in the `submit.sh` bash script in the [scripts/](https://github.com/SahakyanLab/GenomeAssembler_dev/tree/master/scripts/) folder.

## Run all setup files

If you wish to run all setups, including all the aforementioned bash scripts, please run the below bash script. 

```bash
bash run_all_setup_files.sh
```

## Run the full genomeassembler_dev study

To run the full study, please run the below bash script. Please ensure you have followed all the aforementioned steps. The full study is expected to take a few days to complete. 

```bash
bash run_genomeassembler_dev.sh
```