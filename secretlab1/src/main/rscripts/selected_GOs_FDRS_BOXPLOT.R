#/usr/bin/Rscript

library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
	print("ERROR: no arguments provided")
}

nodes <- args[1]
outdir <- args[2]

nodes <- nodes[order(nodes[,502]),]
nodes <- nodes[1:5,]

nodes <- reshape2::melt(nodes, id.vars=c(1,2), variable.name="FDR")

ggplot(nodes, aes(x=V1, y=-log10(FDR))) + geom_boxplot() + ylab("-log10(FDR)") + xlab("GO node(s)")
ggsave(paste0(outdir, .Platform$file.sep, "selected_nodes_fdrs_BOXPLOT.pdf"), width=10, height=10)