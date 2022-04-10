#/usr/bin/Rscript

options(warn = - 1)

library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
	print("no arguments provided")
}

genes <- genes <- read.csv(args[1], sep="\t", header=TRUE)
outdir <- args[2]

genes <- subset(genes, genes$is_unclear == "false")

ggplot(genes, aes(weighted_score)) + stat_ecdf(geom="step") + ylab("% genes") + xlab("weighted score")
ggsave(paste0(outdir, .Platform$file.sep, "genes_weighted_score_CUMULATIVE.pdf"), width=10, height=10)