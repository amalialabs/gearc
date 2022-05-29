
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

genes <- genes[,c(1,2,6)]

fdrcutoff <- 0.05 #* nrow(genes)
num_sig_genes <- nrow(subset(genes, genes$is_sig == "true"))
num_nonsig_genes <- nrow(subset(genes, genes$is_sig == "false"))
num_one_percent_sig <- max(round(0.01 * num_sig_genes), 5)
num_five_percent_nonsig <- round(0.05 * num_nonsig_genes)

num_extend <- min(num_one_percent_sig, num_five_percent_nonsig)

genes <- genes[order(genes$FDR),]

num_last_gene <- min(num_sig_genes + num_extend, nrow(genes))
difffdr <- abs(genes[num_last_gene,2]-fdrcutoff)
tmp <- subset(genes, genes$FDR<=fdrcutoff+difffdr)
first <- tmp[nrow(tmp),2]
firstp <- nrow(tmp)/nrow(genes)
tmp <- subset(genes, genes$FDR<=fdrcutoff-difffdr)
sec <- tmp[nrow(tmp), 2]
secp <- nrow(tmp)/nrow(genes)


p1 <- ggplot(genes, aes(FDR)) + stat_ecdf(geom="step") + ylab("% genes") + xlab("FDR") +
	geom_vline(xintercept=fdrcutoff, col="red") +
	geom_point(aes(first, firstp, col="blue")) +
	geom_point(aes(sec, secp, col="blue")) +
	scale_color_identity()


num_extend <- num_one_percent_sig

num_last_gene <- min(num_sig_genes + num_extend, nrow(genes))
difffdr <- abs(genes[num_last_gene,2]-fdrcutoff)
tmp <- subset(genes, genes$FDR<=fdrcutoff+difffdr)
first <- tmp[nrow(tmp),2]
firstp <- nrow(tmp)/nrow(genes)
tmp <- subset(genes, genes$FDR<=fdrcutoff-difffdr)
sec <- tmp[nrow(tmp), 2]
secp <- nrow(tmp)/nrow(genes)


p2 <- ggplot(genes, aes(FDR)) + stat_ecdf(geom="step") + ylab("% genes") + xlab("FDR") +
	scale_color_identity() +
	geom_point(aes(first, firstp, col="darkblue")) +
	geom_point(aes(sec, secp, col="darkblue")) +
	geom_vline(xintercept=fdrcutoff, col="red") +
	coord_cartesian(xlim=c(0,0.1))


num_extend <- num_five_percent_nonsig

num_last_gene <- min(num_sig_genes + num_extend, nrow(genes))
difffdr <- abs(genes[num_last_gene,2]-fdrcutoff)
tmp <- subset(genes, genes$FDR<=fdrcutoff+difffdr)
first <- tmp[nrow(tmp),2]
firstp <- nrow(tmp)/nrow(genes)
tmp <- subset(genes, genes$FDR<=fdrcutoff-difffdr)
sec <- tmp[nrow(tmp), 2]
secp <- nrow(tmp)/nrow(genes)


p3 <- ggplot(genes, aes(FDR)) + stat_ecdf(geom="step") + ylab("% genes") + xlab("FDR") +
	scale_color_identity() +
	geom_point(aes(first, firstp, col="darkblue")) +
	geom_point(aes(sec, secp, col="darkblue")) +
	geom_vline(xintercept=fdrcutoff, col="red") +
	coord_cartesian(xlim=c(0,0.1))



plots <- list(p1,p2,p3)
ggsave(paste0(outdir, .Platform$file.sep, "fdr_cutoff_finding_all_CUMULATIVE.pdf"), marrangeGrob(grobs = plots, nrow=1, ncol=1), device = "pdf")


