#/usr/bin/Rscript

library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
	print("ERROR: no arguments provided")
}

nodes <- args[2]
fdrs <- args[3]
outdir <- args[4]
df <- as.data.frame(cbind(rep(genes, each=1000), fdrs))
ggplot(df, aes(x=nodes, y=-log10(fdrs))) + geom_boxplot() + ylab("-log10(FDR)") + xlab("GO node(s)")
ggsave(paste0(outdir, .Platform$file.sep, "selected_nodes_fdrs_BOXPLOT.pdf"), width=10, height=10)