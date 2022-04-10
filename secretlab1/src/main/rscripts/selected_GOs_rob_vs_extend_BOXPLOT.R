#/usr/bin/Rscript

library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
	print("ERROR: no arguments provided")
}

nodes <- args[1]
nodes_ext <- args[2]
outdir <- args[3]

nodes <- nodes[order(nodes[,ncol(nodes)]),]
nodes <- nodes[1:2,]
nodes_ext <- nodes_ext[order(nodes_ext[,ncol(nodes_ext)]),]
nodes_ext <- nodes_ext[1:2,]

nodes <- reshape2::melt(nodes, id.vars=c(1,2), variable.name="FDR")
nodes$value <- as.numeric(as.character(nodes$value))
nodes$type <- "robust"
nodes_ext <- reshape2::melt(nodes_ext, id.vars=c(1,2), variable.name="FDR")
nodes_ext$value <- as.numeric(as.character(nodes_ext$value))
nodes_ext$type <- "robust + extended"

n <- rbind(nodes, nodes_ext)

ggplot(n, aes(x=V1, y=-log10(value)), fill=type) + geom_boxplot() + ylab("-log10(FDR)") + xlab("GO node(s)")
ggsave(paste0(outdir, .Platform$file.sep, "selected_nodes_rob_vs_ext_BOXPLOT.pdf"), width=10, height=10)




