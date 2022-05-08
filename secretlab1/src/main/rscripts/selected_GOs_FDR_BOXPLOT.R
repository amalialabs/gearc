#/usr/bin/Rscript

library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
	print("ERROR: no arguments provided")
}

nodes <- read.csv(args[1], sep="\t", header=FALSE)
outdir <- args[2]

num_cols <- ncol(nodes)

nodes <- nodes[order(nodes[,num_cols]),]
nodes <- nodes[1:5,]

nodes <- reshape2::melt(nodes, id.vars=c(1,2), variable.name="FDR")
nodes$value <- as.numeric(as.character(nodes$value))
mn <- max(-log10(nodes$value))

ggplot(nodes, aes(x=V1, y=-log10(value), fill=V1)) + geom_boxplot() + ylab("-log10(FDR)") +
	xlab("GO node(s)") + ylim(0,mn) + theme(legend.position="none")
ggsave(paste0(outdir, .Platform$file.sep, "selected_nodes_fdrs_BOXPLOT.pdf"), width=10, height=10)