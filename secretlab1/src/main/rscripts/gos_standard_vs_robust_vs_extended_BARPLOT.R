
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

ncol_95quant <- ncol(nodes)*0.95+2
nodes <- nodes[,c(1,ncol_95quant)] #95% quantile
colnames(nodes)[2] <- "robust"
nodes_ext <- ndoes_ext[,c(1,ncol_95quant)]
colnames(nodes_ext)[2] <- "robust + extended"
colnames(nodes_standard)[2] <- "standard"

n <- merge(nodes, nodes_ext, by="V1")
n <- merge(n, nodes_standard, by="V1")

n2 <- reshape2::melt(n, id.vars=1, variable.name="type")


ggplot(n2, aes(x=type, y=-log10(value)), fill=type) + geom_boxplot() + ylab("-log10(FDR)") + xlab("")
ggsave(paste0(outdir, .Platform$file.sep, "gos_standard_vs_robust_vs_extended_BARPLOT.pdf"), width=10, height=10)




