#/usr/bin/Rscript

library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
	print("no arguments provided")
}

genes <- read.csv(args[1], sep="\t", header=TRUE)
outdir <- args[2]

#groups: sig-diff, unclear, sig-nondiff, (rest)
#df$group <- factor(df$group, levels=c())

genes <- subset(genes, genes$is_unclear == "false")

ggplot(genes, aes(x=geneset, fill=geneset)) + geom_bar() + ylab("# genes") + xlab("")
ggsave(paste0(outdir, .Platform$file.sep, "genes_categories_barplot.pdf"), width=10, height=10)