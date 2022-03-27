#/usr/bin/Rscript

library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
	print("ERROR: no arguments provided")
}

print(args)

genes <- read.csv(args[1], sep="\t", header=TRUE)
outdir <- args[2]

ggplot(genes, aes(x=log2FC, y=-log10(FDR))) + geom_point() + ylab("-log10(FDR)") + xlab("log2FC")
ggsave(paste0(outdir, .Platform$file.sep, "significant_genes_volcano.pdf"), width=10, height=10)