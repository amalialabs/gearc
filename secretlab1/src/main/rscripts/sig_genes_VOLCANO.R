#/usr/bin/Rscript

library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
	print("ERROR: no arguments provided")
}

print(args)

genes <- read.csv(args[1], sep="\t", header=TRUE)
#fdr <- args[2]
#fc <- args[3]
outdir <- args[4]

ggplot(genes, aes(x=log2FC, y=-log10(FDR), col=ifelse(is_sig, "red", "black"))) + geom_jitter() +
ylab("-log10(FDR)") + xlab("log2FC") + scale_color_identity()
ggsave(paste0(outdir, .Platform$file.sep, "significant_genes_volcano.pdf"), width=10, height=10)