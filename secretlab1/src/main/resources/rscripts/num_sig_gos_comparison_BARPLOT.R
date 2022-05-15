
#/usr/bin/Rscript

library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
	print("no arguments provided")
}

gos <- read.csv(args[1], sep="\t", header=TRUE)
gos_extended <- read.csv(args[2], sep="\t", header=TRUE)
gos_standard <- read.csv(args[3], sep="\t", header=TRUE)
outdir <- args[4]

iters <- ncol(gos)

gos <- gos[, c(1,2,iters)]
gos_extended <- gos_extended[,c(1,2,iters)]

gos$type <- "robust"
gos_extended$type <- "robust + extended"
gos_standard$type <- "standard"
colnames(gos)[3] <- "FDR"
colnames(gos_extended)[3] <- "FDR"
colnames(gos_standard)[3] <- "FDR"
colnames(gos)[1] <- "id"
colnames(gos)[2] <- "name"
colnames(gos_extended)[1] <- "id"
colnames(gos_extended)[2] <- "name"
colnames(gos_standard)[1] <- "id"
colnames(gos_standard)[2] <- "name"

g <- rbind(gos, gos_extended, gos_standard)
g <- subset(g, g$FDR <= 0.05)
g$type <- factor(g$type, levels=c("standard", "robust", "robust + extended"))


ggplot(g, aes(type)) + geom_bar(width=0.5) + ylab("# genes") + xlab("enrichment type")
ggsave(paste0(outdir, .Platform$file.sep, "num_sig_gos_comparison_BARPLOT.pdf"), width=10, height=10)


