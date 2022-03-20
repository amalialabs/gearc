#/usr/bin/Rscript

library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
	print("ERROR: no arguments provided")
}

genes <- args[2]
fcs <- args[3]
frds <- args[4]
outdir <- args[5]
df <- as.data.frame(cbind(genes, fcs, fdrs))
ggplot(df, aes(x=fcs, y=-log10(fdrs))) + geom_point() + ylab("-log10(FDR)") + xlab("log2FC")
ggsave(paste0(outdir, .Platform$file.sep, "significant_genes_volcano.pdf"), width=10, height=10)