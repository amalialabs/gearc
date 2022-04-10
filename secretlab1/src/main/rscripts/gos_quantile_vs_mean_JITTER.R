
#/usr/bin/Rscript

library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
	print("ERROR: no arguments provided")
}

nodes <- args[1]
nodes_ext <- args[2]
outdir <- args[3]

nodes$mean <- lapply(nodes[,3:ncol(nodes)], mean)
nodes_ext$mean <- lapply(nodes_ext[,3:ncol(nodes_ext)], mean)
ncol_95quant <- ncol(nodes)*0.95+2
nodes$quant <- nodes[,ncol_95quant]  #95% quantile
nodes_ext$quant <- nodes[,ncol_95quant]
nodes <- nodes[,c(1,ncol(nodes)-1, ncol(nodes))]   #nur die beiden letzt berechneten Werte
nodes_ext <- nodes_ext[,c(1,ncol(nodes)-1, ncol(nodes))]
nodes$type <- "robust"
nodes_ext$type <- "robust + extended"
n <- rbind(nodes, nodes_ext)


ggplot(n, aes(x=-log10(mean), y=-log10(quant)), fill=type) + geom_boxplot() + ylab("-log10(FDR 95%)") + xlab("-log10(FDR mean)")
ggsave(paste0(outdir, .Platform$file.sep, "gos_mean_vs_quantile_JITTER.pdf"), width=10, height=10)




