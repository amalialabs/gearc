
#/usr/bin/Rscript

library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
	print("ERROR: no arguments provided")
}

nodes <- args[1]
nodes_ext <- args[2]
nodes_standard <- args[3]
outdir <- args[4]

fdr <- 0.05

ncol_95quant <- ncol(nodes)*0.95+2
nodes <- nodes[,c(1,ncol_95quant)] #95% quantile
colnames(nodes)[2] <- "FDR"
nodes$type <- "robust"
nodes <- subset(nodes, nodes$FDR <= fdr)

nodes_ext <- ndoes_ext[,c(1,ncol_95quant)]
colnames(nodes_ext)[2] <- "FDR"
nodes_ext$type <- "robust + extended"
nodes_ext <- susbet(nodes_ext, nodes_ext$FDR <= fdr)

colnames(nodes_standard)[2] <- "FDR"
nodes_standard$type <- "standard"
nodes_standard <- subset(nodes_standard, nodes_standard$FDR <= fdr)

n <- rbind(nodes, nodes_ext, nodes_standard)


ggplot(n, aes(x=type, y=-log10(FDR)), fill=type) + geom_boxplot() + ylab("-log10(FDR)") + xlab("")
ggsave(paste0(outdir, .Platform$file.sep, "gos_standard_vs_robust_vs_extended_BARPLOT.pdf"), width=10, height=10)




