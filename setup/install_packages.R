reg_libs <- c(
    # data wrangling
    "data.table", "dplyr", "stringr",

    # plotting
    "ggplot2", "ggsignif", "RColorBrewer",

    # others
    "pbapply", "Rcpp"
)
to_install <- reg_libs[!reg_libs %in% installed.packages()]
for(lib in to_install){
    install.packages(
        lib, 
        dependencies = TRUE,
        repos = 'http://cran.uk.r-project.org'
    )
}

bioscience_libs <- c(
    # data wrangling and more
    "plyranges", "Biostrings",

    # reference genomes
    "BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0"
)
to_install <- bioscience_libs[!bioscience_libs %in% installed.packages()]
for(lib in to_install){
    BiocManager::install(lib, update = TRUE, ask = FALSE)
}