# load dependencies
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(data.table)))
suppressPackageStartupMessages(suppressWarnings(library(Rcpp)))
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))
suppressPackageStartupMessages(suppressWarnings(library(pbapply)))
suppressPackageStartupMessages(suppressWarnings(library(Biostrings)))
suppressPackageStartupMessages(suppressWarnings(library(plyranges)))
suppressPackageStartupMessages(suppressWarnings(library(BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0)))
pbapply::pboptions(char = "=", type = "txt")

args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])

setwd(my.path)

source("../lib/GenerateReads.R")
source("../lib/DeNovoAssembler.R")
Rcpp::sourceCpp("../lib/DeNovoAssembler.cpp")

seed = 1234
break_kmer = 8
total_iters = 200
industry_standard = FALSE
coverage_target = 40
seq_len = 1000

loop_grid <- tibble(
    read_lens = c(12,14,16,18,20,25,40),
    all_dbg_kmers = c(9,9,13,15,15,15,15)
)

for(row in 1:nrow(loop_grid)){
    for(i in 1:total_iters){
        cat(paste0("Industry standard: ", industry_standard, ".\n"))

        set.seed(seed = seed)
        assembler <- DeNovoAssembler$new(
            seq_len = seq_len,
            read_len = loop_grid$read_lens[row],
            coverage_target = coverage_target,
            kmer = break_kmer,
            dbg_kmer = loop_grid$all_dbg_kmers[row],
            seed = seed,
            ind = i,
            only_kmers_from_reads = FALSE,
            save_read_files = TRUE,
            industry_standard = industry_standard,
            plot_results = FALSE
        )
        assembler$run_assembler(total_iters = total_iters)
    }
}

##################################################################################
#' Import all a) real and b) randomly sampled probability experiments.
#' Plot on x-axis a and b, y-axis the recovered contig length.
##################################################################################
res <- pbapply::pblapply(1:nrow(loop_grid), function(x){
    files <- list.files(
        path = "../data/results",
        pattern = paste0(
            "SolutionsTable_SeqLen-", format(seq_len, scientific = FALSE), 
            "_SeqSeed-", seed, 
            "_ReadLen-", loop_grid$read_lens[x], 
            "_DBGKmer-", loop_grid$all_dbg_kmers[x],  
            "_kmer-", break_kmer, 
            "_IndustryModel-", industry_standard,
            ".csv"
        ),
        full.names = TRUE,
        recursive = TRUE
    )
    files <- stringr::str_sort(files, numeric = TRUE)
    files <- files[1:total_iters]

    df <- lapply(1:length(files), function(i){
        dt <- fread(
            files[i], 
            select = c(
                "bp_score_norm_by_len_true",
                "bp_score_norm_by_len_random"
            ),
            showProgress = FALSE
        )
        return(colMeans(dt, na.rm = TRUE))
    })
    df <- do.call(rbind, df)
    df <- as.data.table(df)
    df[, `:=`(
        read_len = loop_grid$read_lens[x],
        dbg_kmer = loop_grid$all_dbg_kmers[x]
    )]
    return(df)
})
df <- rbindlist(res)

df <- as_tibble(df) %>% 
    tidyr::pivot_longer(
        cols = -c(read_len, dbg_kmer),
        names_to = "Key",
        values_to = "Value"
    ) %>% 
    dplyr::mutate(
        random_prob = ifelse(
            grepl("random", Key, ignore.case = TRUE), TRUE, FALSE
        ),
        Key = stringr::str_replace(
            string = Key, 
            pattern = "_[^_]*$", 
            replacement = ""
        )
    )

# save results
base_dir <- paste0("IndustryModel_", industry_standard)
dir.create(paste0("../data/", base_dir), showWarnings = FALSE, recursive = TRUE)
dir.create(paste0("../figures/", base_dir), showWarnings = FALSE, recursive = TRUE)

fwrite(df, paste0("../data/", base_dir, "/results_summary.csv"))

##################################################################################
#' To demonstrate:
#' 1. Show that my solutions are better than random probabilities.
##################################################################################
unique_read_len <- unique(df$read_len)
unique_read_len <- paste0("Read len: ", sort(unique_read_len, decreasing = TRUE))

p3 <- df %>% 
    dplyr::filter(!is.na(Value)) %>% 
    dplyr::mutate(
        solution = ifelse(random_prob, "Random", "Non-random"),
        solution = as.factor(solution),
        read_len = paste0("Read len: ", read_len),
        read_len = factor(read_len, levels = unique_read_len)
    ) %>% 
    ggplot(aes(x = solution, y = Value, fill = solution)) + 
    geom_boxplot(alpha = 0.75) + 
    ggsignif::geom_signif(
        comparisons = list(c("Non-random", "Random")),
        map_signif_level = TRUE,
        test = "t.test",
        textsize = 5,
        y_position = 7.5e-5,
        tip_length = 0
    ) +
    scale_fill_manual(values = c("#2166ac", "#b2182b")) +  
    scale_y_continuous(labels = scales::label_number(scale = 1e5)) +
    facet_wrap(vars(read_len), nrow = 1) + 
    coord_cartesian(ylim = c(0, 8e-5)) + 
    theme_bw() + 
    theme_classic() + 
    theme(text = element_text(size = 15)) + 
    labs(
        x = "",
        y = expression("Breakage score norm. by contig length, x10"^-5*"")
    )

ggsave(
    filename = paste0(
        "../figures/", base_dir, "/",
        "Breakscore_contigs_reference",
        "_IndustryModel-", industry_standard, 
        ".pdf"
    ),
    plot = p3,
    height = 5, width = 16
)

##################################################################################
#' To demonstrate:
#' 1. Breakage score raw and normalised improve with better quality solutions
##################################################################################
all_res <- pbapply::pblapply(1:nrow(loop_grid), function(x){
    files <- list.files(
        path = "../data/results",
        pattern = paste0(
            "SolutionsTable_SeqLen-", format(seq_len, scientific = FALSE), 
            "_SeqSeed-", seed, 
            "_ReadLen-", loop_grid$read_lens[x], 
            "_DBGKmer-", loop_grid$all_dbg_kmers[x],  
            "_kmer-", break_kmer, 
            "_IndustryModel-", industry_standard,
            ".csv"
        ),
        full.names = TRUE,
        recursive = TRUE
    )
    files <- stringr::str_sort(files, numeric = TRUE)
    files <- files[1:total_iters]

    df <- lapply(1:length(files), function(i){
        df <- fread(
            files[i], 
            select = c(
                "sequence_len",
                "kmer_breaks",
                "bp_score_norm_by_break_freqs_true",
                "bp_score_norm_by_len_true",
                "bp_score_true",
                "lev_dist_vs_true",
                "stat_test_KS_true"
            ),
            showProgress = FALSE
        )
        return(df)
    })
    df <- rbindlist(df)
    df[, read_len := loop_grid$read_lens[x]]
    df[, dbg_kmer := loop_grid$all_dbg_kmers[x]]
    return(df)
})
df_all_res <- rbindlist(all_res)
fwrite(df_all_res, paste0("../data/", base_dir, "/results_all.csv"))

#' Show that the top 5% of solutions are better by breakage score
#' compared to all solutions.
top_sol <- as_tibble(df_all_res) %>% 
    dplyr::select(bp_score_norm_by_len_true, read_len) %>% 
    dplyr::rename(bp_score = bp_score_norm_by_len_true) %>% 
    dplyr::group_by(read_len) %>% 
    dplyr::slice_max(order_by = bp_score, prop = 0.05) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(solution = "Top 5%")
all_sol <- as_tibble(df_all_res) %>% 
    dplyr::select(bp_score_norm_by_len_true, read_len) %>% 
    dplyr::rename(bp_score = bp_score_norm_by_len_true) %>% 
    dplyr::slice_min(order_by = bp_score, prop = 0.95) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(solution = "Remaining")
c_sol <- dplyr::bind_rows(top_sol, all_sol)

p4 <- c_sol %>% 
    dplyr::mutate(
        solution = factor(solution, levels = c("Top 5%", "Remaining")),
        read_len = paste0("Read len: ", read_len),
        read_len = factor(read_len, levels = unique_read_len)
    ) %>% 
    ggplot(aes(x = solution, y = bp_score, fill = solution)) + 
    geom_boxplot(alpha = 0.75) + 
    ggsignif::geom_signif(
        comparisons = list(c("Top 5%", "Remaining")),
        map_signif_level = TRUE,
        test = "t.test",
        textsize = 5,
        y_position = 3.2e-4,
        tip_length = 0
    ) +
    scale_fill_manual(values = c("#2166ac", "#b2182b")) +  
    scale_y_continuous(labels = scales::label_number(scale = 1e4)) +
    facet_wrap(vars(read_len), nrow = 1) + 
    coord_cartesian(ylim = c(0, 3.5e-4)) +
    theme_bw() + 
    theme_classic() + 
    theme(text = element_text(size = 15)) + 
    labs(
        x = "",
        y = expression("Breakage score norm. by contig length, x10"^-4*"")
    )

ggsave(
    filename = paste0(
        "../figures/", base_dir, "/",
        "Breakscore_Top-vs-all-solutions",
        "_IndustryModel-", industry_standard, 
        ".png"
    ),
    plot = p4,
    height = 5, width = 16
)

#' Show that the top 5% of solutions are better by normalised breakage score
#' compared to all solutions.
top_sol <- as_tibble(df_all_res) %>% 
    dplyr::select(bp_score_norm_by_len_true, read_len) %>% 
    dplyr::rename(bp_score_norm = bp_score_norm_by_len_true) %>% 
    dplyr::group_by(read_len) %>% 
    dplyr::slice_max(order_by = bp_score_norm, prop = 0.05) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(solution = "Top 5%")
all_sol <- as_tibble(df_all_res) %>% 
    dplyr::select(bp_score_norm_by_len_true, read_len) %>% 
    dplyr::rename(bp_score_norm = bp_score_norm_by_len_true) %>% 
    dplyr::slice_min(order_by = bp_score_norm, prop = 0.95) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(solution = "Remaining")
c_sol <- dplyr::bind_rows(top_sol, all_sol)

p5 <- c_sol %>% 
    dplyr::mutate(
        solution = factor(solution, levels = c("Top 5%", "Remaining")),
        read_len = paste0("Read len: ", read_len),
        read_len = factor(read_len, levels = unique_read_len)
    ) %>% 
    ggplot(aes(x = solution, y = bp_score_norm, fill = solution)) + 
    geom_boxplot(alpha = 0.75) + 
    ggsignif::geom_signif(
        comparisons = list(c("Top 5%", "Remaining")),
        map_signif_level = TRUE,
        test = "t.test",
        textsize = 5,
        y_position = 3.1e-4,
        tip_length = 0
    ) +
    scale_fill_manual(values = c("#2166ac", "#b2182b")) +  
    scale_y_continuous(labels = scales::label_number(scale = 1e4)) +
    facet_wrap(vars(read_len), nrow = 1) + 
    coord_cartesian(ylim = c(0, 3.5e-4)) +
    theme_bw() + 
    theme_classic() + 
    theme(text = element_text(size = 15)) + 
    labs(
        x = "",
        y = expression("Breakage score norm. by contig length, x10"^-4*"")
    )

ggsave(
    filename = paste0(
        "../figures/", base_dir, "/",
        "NormBreakscore_Top-vs-all-solutions",
        "_IndustryModel-", industry_standard, 
        ".png"
    ),
    plot = p5,
    height = 5, width = 16
)

# bin into buckets and plot boxplots of those bins
sol_to_bin <- as_tibble(df_all_res) %>% 
    dplyr::select(
        lev_dist_vs_true, 
        bp_score_norm_by_len_true, 
        bp_score_true,
        bp_score_norm_by_break_freqs_true,
        read_len
    ) %>% 
    dplyr::filter(!is.na(bp_score_norm_by_len_true)) %>% 
    dplyr::mutate(
        read_len = paste0("Read len: ", read_len),
        read_len = factor(read_len, levels = unique_read_len)
    )

nr_bins <- 4
p8 <- sol_to_bin %>% 
    dplyr::filter(
        !is.na(bp_score_norm_by_len_true)
        # lev_dist_vs_true <= 5
    ) %>%
    dplyr::mutate(bins = cut(
        lev_dist_vs_true, 
        breaks = seq(
            from = 0,
            to = max(lev_dist_vs_true), 
            length.out = nr_bins + 1
        ), 
        include.lowest = TRUE
    )) %>% 
    dplyr::ungroup() %>% 
    ggplot(aes(x = bins, y = bp_score_norm_by_len_true, fill = read_len)) + 
    geom_boxplot(alpha = 0.75, outlier.alpha = 0.1) + 
    theme_bw() + 
    theme_classic() + 
    theme(
        text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
    ) + 
    scale_fill_brewer(palette = "Set3") +
    scale_y_continuous(labels = scales::label_number(scale = 1e4)) + 
    coord_cartesian(ylim = c(0, 1e-4)) + 
    facet_wrap(vars(read_len), nrow = 1) + 
    labs(
        x = "Levenshtein distance",
        y = expression("Breakage score norm. by contig length, x10"^-4*"")
    )
    
ggsave(
    filename = paste0(
        "../figures/", base_dir, "/",
        "Binned-Levenshtein-distance_vs_NormBreakscore",
        "_IndustryModel-", industry_standard, 
        ".png"
    ),
    plot = p8,
    height = 5, width = 16
)

p8 <- sol_to_bin %>% 
    dplyr::filter(!is.na(bp_score_true)) %>%
    dplyr::mutate(bins = cut(
        lev_dist_vs_true, 
        breaks = seq(
            from = 0,
            to = max(lev_dist_vs_true), 
            length.out = nr_bins + 1
        ), 
        include.lowest = TRUE
    )) %>% 
    dplyr::ungroup() %>% 
    ggplot(aes(x = bins, y = bp_score_true, fill = read_len)) + 
    geom_boxplot(alpha = 0.75, outlier.alpha = 0.1) + 
    theme_bw() + 
    theme_classic() + 
    theme(
        text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
    ) + 
    scale_fill_brewer(palette = "Set3") +
    # scale_y_continuous(labels = scales::label_number(scale = 1e4)) + 
    # coord_cartesian(ylim = c(0, 0.8e-4)) + 
    facet_wrap(vars(read_len), nrow = 1) + 
    labs(
        x = "Levenshtein distance",
        y = "Breakage score"
    )
    
ggsave(
    filename = paste0(
        "../figures/", base_dir, "/",
        "Binned-Levenshtein-distance_vs_Breakscore",
        "_IndustryModel-", industry_standard, 
        ".png"
    ),
    plot = p8,
    height = 5, width = 16
)

p9 <- as_tibble(df_all_res) %>% 
    dplyr::select(
        lev_dist_vs_true, 
        bp_score_norm_by_len_true, 
        bp_score_true,
        bp_score_norm_by_break_freqs_true,
        read_len
    ) %>% 
    dplyr::filter(
        !is.na(bp_score_norm_by_len_true)
        # read_len >= 16
    ) %>%
    dplyr::mutate(
        read_len = paste0("Read len: ", read_len),
        read_len = factor(read_len, levels = unique_read_len),
        bins = cut(
            lev_dist_vs_true, 
            breaks = seq(
                from = 0,
                to = max(lev_dist_vs_true), 
                length.out = nr_bins + 1
            ), 
            include.lowest = TRUE
        )
    ) %>% 
    dplyr::ungroup() %>% 
    ggplot(aes(x = bins, y = bp_score_norm_by_len_true, fill = read_len)) + 
    geom_boxplot(alpha = 0.75, outlier.alpha = 0.1) + 
    theme_bw() + 
    theme_classic() + 
    theme(
        text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
    ) + 
    scale_fill_brewer(palette = "Set3") + 
    scale_y_continuous(labels = scales::label_number(scale = 1e4)) + 
    coord_cartesian(ylim = c(0, 1e-4)) + 
    facet_wrap(vars(read_len), nrow = 1) + 
    labs(
        x = "",
        y = ""
    )
    
ggsave(
    filename = paste0(
        "../figures/", base_dir, "/",
        "For_Paper_IndustryModel-", industry_standard, 
        ".png"
    ),
    plot = p9,
    height = 7, width = 12
)

p10 <- as_tibble(df_all_res) %>% 
    dplyr::select(
        lev_dist_vs_true, 
        bp_score_norm_by_len_true, 
        bp_score_true,
        bp_score_norm_by_break_freqs_true,
        read_len
    ) %>% 
    dplyr::filter(
        !is.na(bp_score_true)
        # read_len >= 16
    ) %>%
    dplyr::mutate(
        read_len = paste0("Read len: ", read_len),
        read_len = factor(read_len, levels = unique_read_len),
        bins = cut(
            lev_dist_vs_true, 
            breaks = seq(
                from = 0,
                to = max(lev_dist_vs_true), 
                length.out = nr_bins + 1
            ), 
            include.lowest = TRUE
        )
    ) %>% 
    dplyr::ungroup() %>% 
    ggplot(aes(x = bins, y = bp_score_true, fill = read_len)) + 
    geom_boxplot(alpha = 0.75, outlier.alpha = 0.1) + 
    theme_bw() + 
    theme_classic() + 
    theme(
        text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
    ) + 
    scale_fill_brewer(palette = "Set3") + 
    # coord_cartesian(ylim = c(0, 1e-4)) + 
    facet_wrap(vars(read_len), nrow = 1) + 
    labs(
        x = "",
        y = ""
    )

pdf(
    paste0(
        "../figures/", base_dir, "/",
        "For_Paper_IndustryModel-", industry_standard, 
        "_BOTH.pdf"
    ),
    height = 10, width = 12
)
gridExtra::grid.arrange(p10, p9, ncol = 1)
plot.saved <- dev.off()

png(
    paste0(
        "../figures/", base_dir, "/",
        "For_Paper_IndustryModel-", industry_standard, 
        "_BOTH.png"
    ),
    height = 10, width = 12, units = "in", res = 300
)
gridExtra::grid.arrange(p10, p9, ncol = 1)
plot.saved <- dev.off()

# one-way anova test of BPs by each group (read lengths)
df_prep_anova <- as_tibble(df_all_res) %>% 
    dplyr::select(
        lev_dist_vs_true, 
        bp_score_norm_by_len_true, 
        bp_score_true,
        read_len
    ) %>% 
    dplyr::mutate(
        read_len = paste0("Read len: ", read_len),
        read_len = factor(read_len, levels = unique_read_len),
        bins = cut(
            lev_dist_vs_true, 
            breaks = seq(
                from = 0,
                to = max(lev_dist_vs_true), 
                length.out = nr_bins + 1
            ), 
            include.lowest = TRUE
        )
    ) %>% 
    dplyr::ungroup()

anova_test <- aov(lev_dist_vs_true ~ read_len, data = df_prep_anova)
summary(anova_test)

# spearman test
spearman_test <- cor.test(
    x = df_prep_anova$bp_score_true, 
    y = df_prep_anova$lev_dist_vs_true, 
    method = "spearman", 
    alternative = "two.sided", 
    exact = FALSE
)

spearman_rank <- cor(
    x = df_prep_anova$bp_score_true, 
    y = df_prep_anova$lev_dist_vs_true, 
    method = "spearman"
)
abs_spearman_rank <- abs(spearman_rank)