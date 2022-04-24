
#/usr/bin/Rscript

options(warn = - 1)

library(VennDiagram, quietly=TRUE)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
	print("ERROR: no arguments provided")
}

nodes <- read.csv(args[1], sep="\t", header=FALSE)
nodes_ext <- read.csv(args[2], sep="\t", header=FALSE)
nodes_standard <- read.csv(args[3], sep="\t", header=FALSE)
outdir <- args[4]

fdr <- 0.05

ncol_95quant <- round((ncol(nodes)-2)*0.95)+2
nodes <- nodes[,c(1,ncol_95quant)] #95% quantile
nodes <- subset(nodes, nodes[,2] <= fdr)

nodes_ext <- nodes_ext[,c(1,ncol_95quant)]
nodes_ext <- subset(nodes_ext, nodes_ext[,2] <= fdr)

nodes_standard <- subset(nodes_standard, nodes_standard[,2] <= fdr)


pdf(paste0(outdir, .Platform$file.sep, "gos_standard_vs_robust_vs_extended_VENN.pdf"), width=10, height=10)
draw.triple.venn(area1=nrow(nodes), area2=nrow(nodes_ext), area3=nrow(nodes_standard),
                n12=nrow(merge(nodes, nodes_ext, by="V1")),
                n13=nrow(merge(nodes, nodes_standard, by="V1")),
                n23=nrow(merge(nodes_ext, nodes_standard, by="V1")),
                n123=nrow(merge(merge(nodes, nodes_ext, by="V1"), nodes_standard, by="V1")),
                fill=c("red", "blue", "green"), category=c("robust", "robust + extended", "standard"))
dev.off()




