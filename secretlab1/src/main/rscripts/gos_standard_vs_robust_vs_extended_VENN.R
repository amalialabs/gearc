
#/usr/bin/Rscript

library(VennDiagram)

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
nodes <- subset(nodes, nodes[,2] <= fdr)

nodes_ext <- nodes_ext[,c(1,ncol_95quant)]
nodes_ext <- subset(nodes_ext, nodes_ext[,2] <= fdr)

nodes_standard <- subset(nodes_standard, nodes_standard[,2] <= fdr)


pdf(paste0(outdir, .Platform$file.sep, "gos_standard_vs_robust_vs_extended_BARPLOT.pdf"), width=10, height=10)
VennDiagram(n1=nrow(nodes), n2=nrow(nodes_ext), n3=nrow(nodes_standard),
                n12=nrow(merge(nodes, nodes_ext, by="V1")),
                n13=nrow(merge(nodes, nodes_standard, by="V1")),
                n23=nrow(merge(nodes_ext, nodes_standard, by="V1")),
                n123=nrow(merge(merge(nodes, nodes_ext, by="V1")), nodes_standard, by="V1"),
                col=c("red", "blue", "green"), labels=c("robust", "robust + extended", "standard"))
dev.off()




