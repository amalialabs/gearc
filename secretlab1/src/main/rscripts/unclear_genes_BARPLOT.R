#/usr/bin/Rscript

library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
	print("no arguments provided")
}

genes <- read.csv(args[1], sep="\t", header=TRUE)
outdir <- args[2]

genes$unclear <- ifelse(genes$is_unclear=="true", "unclear", "clear")

ggplot(genes, aes(x=unclear, fill=unclear)) + geom_bar() + ylab("# genes") + xlab("")
ggsave(paste0(outdir, .Platform$file.sep, "unclear_genes_BARPLOT.pdf"), width=10, height=10)