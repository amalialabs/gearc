#/usr/bin/Rscript

options(warn = - 1)

library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
	print("ERROR: no arguments provided")
}


genes <- read.csv(args[1], sep="\t", header=TRUE)
#fdr <- args[2]
#fc <- args[3]
outdir <- args[2]

genes <- subset(genes, genes$is_unclear == "false")

ggplot(genes, aes(x=log2FC, y=-log10(FDR), col=ifelse(is_sig, "red", "black"))) + geom_jitter() +
ylab("-log10(FDR)") + xlab("log2FC") + scale_color_identity()
ggsave(paste0(outdir, .Platform$file.sep, "significant_genes_VOLCANO.pdf"), width=10, height=10)