
#/usr/bin/Rscript

options(warn = - 1)

library(ggplot2)
library(gridExtra)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
	print("no arguments provided")
}

genes <- read.csv(args[1], sep="\t", header=TRUE)
outdir <- args[2]

genes <- subset(genes, genes$is_unclear == "false")

genes <- genes[,c(1,3,6)]

fccutoff <- 1.0
idx <- nrow(subset(genes, abs(genes$log2FC) >= fccutoff)) / nrow(genes)


num_sig_genes <- nrow(subset(genes, genes$is_sig == "true"))
num_nonsig_genes <- nrow(subset(genes, genes$is_sig == "false"))
num_one_percent_sig <- max(round(0.01 * num_sig_genes), 5)
num_five_percent_nonsig <- round(0.05 * num_nonsig_genes)

num_extend <- min(num_one_percent_sig, num_five_percent_nonsig)

genes$absfc <- abs(genes$log2FC)
genes <- genes[order(genes$absfc, decreasing=TRUE),]
fcext <- genes[num_sig_genes+num_extend,4]

percent_extm <- nrow(subset(genes, genes$log2FC<=-fcext))/nrow(genes)
tmp <- subset(genes, genes$log2FC<=-fcext)
extmfc <- tmp[nrow(tmp), 2]
percent_ext <- nrow(subset(genes, genes$log2FC<=fcext))/nrow(genes)
tmp <- subset(genes, genes$log2FC>=fcext)
extfc <- tmp[nrow(tmp),2]

p1 <- ggplot(genes, aes(log2FC)) + stat_ecdf(geom="step") + ylab("% genes") + xlab("log2FC") +
	geom_vline(xintercept=c(-1,1), col="red") +
	geom_point(aes(extmfc, percent_extm, col="blue")) +
	geom_point(aes(extfc, percent_ext, col="blue")) +
	scale_color_identity()

p2 <- ggplot(genes, aes(log2FC)) + stat_ecdf(geom="step") + ylab("% genes") + xlab("log2FC") +
	geom_vline(xintercept=c(-1,1), col="red") +
	geom_point(aes(extmfc, percent_extm, col="blue")) +
	geom_point(aes(extfc, percent_ext, col="blue")) +
	scale_color_identity() + coord_cartesian(xlim=c(-2,2))


num_extend <- num_one_percent_sig
fcext <- genes[num_sig_genes+num_extend,4]

percent_extm <- nrow(subset(genes, genes$log2FC<=-fcext))/nrow(genes)
tmp <- subset(genes, genes$log2FC<=-fcext)
extmfc <- tmp[nrow(tmp), 2]
percent_ext <- nrow(subset(genes, genes$log2FC<=fcext))/nrow(genes)
tmp <- subset(genes, genes$log2FC>=fcext)
extfc <- tmp[nrow(tmp),2]


p3 <- ggplot(genes, aes(log2FC)) + stat_ecdf(geom="step") + ylab("% genes") + xlab("log2FC") +
	geom_vline(xintercept=c(-1,1), col="red") +
	geom_point(aes(extmfc, percent_extm, col="blue")) +
	geom_point(aes(extfc, percent_ext, col="blue")) +
	scale_color_identity()

p4 <- ggplot(genes, aes(log2FC)) + stat_ecdf(geom="step") + ylab("% genes") + xlab("log2FC") +
	geom_vline(xintercept=c(-1,1), col="red") +
	geom_point(aes(extmfc, percent_extm, col="blue")) +
	geom_point(aes(extfc, percent_ext, col="blue")) +
	scale_color_identity() + coord_cartesian(xlim=c(-2,2))



num_extend <- num_five_percent_nonsig
fcext <- genes[num_sig_genes+num_extend,4]

percent_extm <- nrow(subset(genes, genes$log2FC<=-fcext))/nrow(genes)
tmp <- subset(genes, genes$log2FC<=-fcext)
extmfc <- tmp[nrow(tmp), 2]
percent_ext <- nrow(subset(genes, genes$log2FC<=fcext))/nrow(genes)
tmp <- subset(genes, genes$log2FC>=fcext)
extfc <- tmp[nrow(tmp),2]


p5 <- ggplot(genes, aes(log2FC)) + stat_ecdf(geom="step") + ylab("% genes") + xlab("log2FC") +
	geom_vline(xintercept=c(-1,1), col="red") +
	geom_point(aes(extmfc, percent_extm, col="blue")) +
	geom_point(aes(extfc, percent_ext, col="blue")) +
	scale_color_identity()

p6 <- ggplot(genes, aes(log2FC)) + stat_ecdf(geom="step") + ylab("% genes") + xlab("log2FC") +
	geom_vline(xintercept=c(-1,1), col="red") +
	geom_point(aes(extmfc, percent_extm, col="blue")) +
	geom_point(aes(extfc, percent_ext, col="blue")) +
	scale_color_identity() + coord_cartesian(xlim=c(-2,2))




plots <- list(p1, p2, p3, p4, p5, p6)
ggsave(paste0(outdir, .Platform$file.sep, "fc_cutoff_finding_all_CUMULATIVE.pdf"), marrangeGrob(grobs = plots, nrow=1, ncol=1), device = "pdf")
