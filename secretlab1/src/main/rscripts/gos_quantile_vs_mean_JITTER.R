
#/usr/bin/Rscript

library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
	print("ERROR: no arguments provided")
}

nodes <- read.csv(args[1], sep="\t", header=FALSE)
nodes_ext <- read.csv(args[2], sep="\t", header=FALSE)
outdir <- args[3]

num_cols <- ncol(nodes)
num_cols2 <- ncol(nodes_ext)

nodes$mean <- lapply(nodes[,3:num_cols], mean)
nodes_ext$mean <- lapply(nodes_ext[,3:num_cols2], mean)
ncol_95quant <- num_cols*0.95+2
nodes$quant <- nodes[,ncol_95quant]  #95% quantile
nodes_ext$quant <- nodes[,ncol_95quant]
nodes <- nodes[,c(1,num_cols-1, num_cols)]   #nur die beiden letzt berechneten Werte
nodes_ext <- nodes_ext[,c(1,num_cols2-1, num_cols2)]
nodes$type <- "robust"
nodes_ext$type <- "robust + extended"
n <- rbind(nodes, nodes_ext)


ggplot(n, aes(x=-log10(mean), y=-log10(quant)), fill=type) + geom_boxplot() + ylab("-log10(FDR 95%)") + xlab("-log10(FDR mean)")
ggsave(paste0(outdir, .Platform$file.sep, "gos_mean_vs_quantile_JITTER.pdf"), width=10, height=10)




