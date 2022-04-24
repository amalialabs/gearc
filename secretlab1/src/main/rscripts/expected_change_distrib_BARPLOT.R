
#/usr/bin/Rscript

library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
	print("no arguments provided")
}

genes <- read.csv(args[1], sep="\t", header=TRUE)
outdir <- args[2]

genes <- genes[,c(1,4,7)]
genes <- subset(genes, genes$is_unclear=="false")

genes <- genes[,c(1,2)]


ggplot(genes, aes(geneset)) + geom_bar() + ylab("# genes") + xlab("expected change distribution")
ggsave(paste0(outdir, .Platform$file.sep, "expected_change_distrib_BARPLOT.pdf"), width=10, height=10)


