
#/usr/bin/Rscript

options(warn = - 1)

library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  print("no arguments provided")
}

rankdiffs <- read.csv(args[1], sep="\t", header=FALSE)
outdir <- args[2]

rankdiffs <- as.data.frame(rankdiffs)
colnames(rankdiffs)[1] <- "rankdiff"

ggplot(rankdiffs, aes(rankdiff)) + stat_ecdf(geom="step") + ylab("% genes") + xlab("rank differences")

ggplot(rankdiff, aes(rankdiff)) + geom_histogram(bins=40) + xlab("rank differences") +
  geom_vline(xintercept=median(rankdiff$rankdiff), col="red") +
  geom_vline(xintercept=median(rankdiff$rankdiff)+sd(rankdiffs$rankdiffs), col="orange") +
  annotate("text", x=median(z$z)+sd(z$z)+0.6, y=10, label="-> unclear")

ggsave(paste0(outdir, .Platform$file.sep, "de_scores_ranked_PLOTS.pdf"), width=10, height=10)

