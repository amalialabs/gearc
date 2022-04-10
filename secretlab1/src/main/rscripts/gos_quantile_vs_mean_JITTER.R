
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
nodes$quant <- nodes[,477]  #95% quantile
nodes_ext$quant <- nodes[,477]
nodes <- nodes[,c(1,503:504)]   #nur die beiden Werte
nodes_ext <- nodes_ext[,c(1,503:504)]
nodes$type <- "robust"
nodes_ext$type <- "robust + extended"
n <- rbind(nodes, nodes_ext)


ggplot(n, aes(x=-log10(mean), y=-log10(quant)), fill=type) + geom_boxplot() + ylab("-log10(FDR 95%)") + xlab("-log10(FDR mean)")
ggsave(paste0(outdir, .Platform$file.sep, "gos_mean_vs_quantile_JITTER.pdf"), width=10, height=10)




