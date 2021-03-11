###########################################################################################
###### This function plots the results obtained from function MetricsAnalysis.R      ###### 
###########################################################################################

# Plot of breakage probability scores versus alignment scores
# Metrics Analysis -> analysis parameter -> "breakpoints"
library("ggplot2")
p <- ggplot(results[[1]], aes(x = prob.table, y = total.align.score)) + geom_point()
p <- p + labs(x = "Breakage Probability Scores")
p

pdf(width=8, height=6, 
    file=paste0("../figures/SCORE-vs-PROBS-", 
                "Runs_",max.runs,"-",
                "Len_",length, "-",
                "Kmer",kmer,"-",
                "ReadLength",read.length,".pdf"))
p
dev.off()

##################################################################
# Plot of breakage probability scores versus all other alignment metrics
# Metrics Analysis -> analysis parameter -> "breakpoints"
library("reshape2")
results.copy <- subset(results[[1]], 
                       select = -c(total.align.score, all.sequences, ref.seq))
dat.m <- melt(results.copy, id.vars = "prob.table")

pdf(width=8, height=6, 
    file=paste0("../figures/ALL_METRICS-vs-PROBS-", 
                "Runs_",max.runs,"-",
                "Len_",length, "-",
                "Kmer",kmer,"-",
                "ReadLength",read.length,".pdf"))
ggplot(dat.m, aes(prob.table, value, colour = variable)) +
  geom_point() + labs(x = "Breakage Probability Scores",
                      y = "Alignment Scores Metrics")
dev.off()

##################################################################
# Barplot of the correlation coefficients of breakage probability scores 
# versus each alignment metric
# Metrics Analysis -> analysis parameter -> "breakpoints"
pdf(width=8, height=6, 
    file=paste0("../figures/Correlation-Coefficients-", 
                "Runs_",max.runs,"-",
                "Len_",length, "-",
                "Kmer",kmer,"-",
                "ReadLength",read.length,".pdf"))
results.rsq <- as.data.frame(t(results[[2]]))
barplot(results.rsq$V1, 
        names.arg = c("local", "global", "global-local", "local-global", "kmer", "levenshtein"),
        beside = TRUE,
        ylim = c(0,1),
        xlab = "Alignment Scores Metrics", ylab = "Correlation Coefficient")
dev.off()

##################################################################
# Barplot showing the percent of correctness of breakpoints versus 
# each alignment metric; we define "correctness" as a breakpoint score that is
# equal to or greater than the alignment score
# Metrics Analysis -> analysis parameter -> "breakpoints.performance"
pdf(width=8, height=6, 
    file=paste0("../figures/Alignment_scores/", 
                "Runs_",max.runs,"-",
                "Len_",length, "-",
                "Kmer",kmer,"-",
                "ReadLength",read.length,".pdf"))
barplot(results,
        names.arg = c("local", "global", "global-local", "local-global", "kmer", "levenshtein"),
        ylim = c(0,100),
        # beside = TRUE, 
        legend.text = TRUE,
        args.legend = list(x = "topright", bty = "y", box.lwd = 1, bg = "white"),
        xlab = "Alignment Scores Metrics", ylab = "Percent of correctness",
        main = "Percent of correctness between breakpoints and alignment scores")
dev.off()

##################################################################
# Plot showing the read length as a function of the levenshtein distance, 
# signifying the breakpoint scoring performance 
# Metrics Analysis -> analysis parameter -> "breakpoints.performance"

# df <- data.frame(levenshtein = 214,
#                  read.length = 30)
results[[2]]
df
df <- rbind(df, results[[2]])
pdf(width=8, height=6, 
    file=paste0("../figures/BREAKPOINTS-PERFORMANCE-", 
                "Runs_",max.runs,"-",
                "Len_",length, "-",
                "Kmer",kmer,"-",
                "ReadLength",read.length,".pdf"))
p <- ggplot(df, aes(x = read.length, y = levenshtein)) + geom_point()
p <- p + labs(x = "Read Length", y = "Levenshtein Distance", 
              title = "Breakpoints Performance")
p
dev.off()

##################################################################
# 3D plotly scatter plot of local versus global versus each of the other
# alignment metric scores
# Metrics Analysis -> analysis parameter -> "r-squared"
library("plotly")
rows.in.df    <- dim(results$alignment.res)[1]
cols          <- 4
max.score.ind <- match(max(results$alignment.res$score), 
                       results$alignment.res$score)

# create data frame in the right format for a 3D interactive scatter plot
normalised.align.res <- data.frame(
  local  = c(rep(results$alignment.res$local, cols)),
  global = c(rep(results$alignment.res$global, cols)),
  all    = c(scale(results$alignment.res[,"kmer"])[,1],
             scale(results$alignment.res[,"levenshtein"])[,1],
             scale(results$alignment.res[,"global.local"])[,1],
             scale(results$alignment.res[,"local.global"])[,1]))

normalised.align.res$type <- c(rep("kmer", rows.in.df),
                               rep("levenshtein", rows.in.df),
                               rep("global.local", rows.in.df),
                               rep("local.global", rows.in.df))

# markers need to be in t he right mode
normalised.align.res$type <- as.factor(normalised.align.res$type)

# 3D scatter plot
fig <- plot_ly(data = normalised.align.res, 
               x = ~local, 
               y = ~global, 
               z = ~all,
               color = ~type)
fig <- fig %>% add_markers()
fig <- fig %>%
  # Highlight the highest score
  add_trace(
    x = c(results$alignment.res[max.score.ind, "local"]), 
    y = c(results$alignment.res[max.score.ind, "global"]), 
    z = c(max(scale(results$alignment.res[,"score"])[,1])),
    marker = list(
      color = 'rgb(255, 0, 0)',
      size = 15,
      line = list(
        color = 'rgb(255, 0, 0)',
        width = 5
      )),
    type = 'scatter3d', 
    mode = 'markers',
    hoverinfo = 'Highest Score',
    name = 'Highest Score',
    showlegend = F)
fig <- fig %>% layout(scene = list(
  xaxis = list(title = 'local score'),
  yaxis = list(title = 'global score'),
  zaxis = list(title = 'normalised scores')))

##################################################################
# Violin plots of multi-regression analysis
# Metrics Analysis -> analysis parameter -> "r-squared"
plots(results = results, plot.type = "violin",
      len = results$sequence.length[1], k = results$kmer[1], 
      max.runs = dim(results)[1])

results.violin <- data.frame(
  freq  = c(results$r.squared.kmer,
            results$r.squared.levenshtein,
            results$r.squared.local.global,
            results$r.squared.global.vs.local),
  metric = c(rep("kmer", dim(results)[1]),
             rep("levenshtein", dim(results)[1]),
             rep("local.global", dim(results)[1]),
             rep("global.vs.local", dim(results)[1])))

p <- ggplot(results.violin, aes(x=metric, y=freq)) + 
  geom_violin() + theme_minimal()
p <- p + geom_violin(aes(fill = factor(metric))) + 
  geom_boxplot(width=0.1, fill="white") +
  labs(x = "Metric", y = "Correlation Coefficient", 
       title = paste("Sequence length = ",length, 
                     "kmer = ", kmer, "runs = ", max.runs))
p

