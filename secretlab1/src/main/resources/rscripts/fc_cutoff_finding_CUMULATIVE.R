
#/usr/bin/Rscript

library(ggplot2)

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

num_last_gene <- min(num_sig_genes + num_extend, nrow(genes))
num_diff_to_center <- num_last_gene - num_sig_genes
num_first_gene <- max(num_sig_genes - num_extend, 1)


ggplot(genes, aes(log2FC)) + stat_ecdf(geom="step") + ylab("% genes") + xlab("log2FC") +
    geom_hline(yintercept=idx, col="red") + geom_point(aes(num_first_gene, genes[num_first_gene, 2], col="blue")) +
    geom_point(aes(num_last_gene, genes[num_last_gene, 2], col="blue"))
ggsave(paste0(outdir, .Platform$file.sep, "fc_cutoff_finding_CUMULATIVE.pdf"), width=10, height=10)


num_last_gene <- min(num_sig_genes + num_one_percent_sig, nrow(genes))
num_first_gene <- max(num_sig_genes - num_one_percent_sig, 1)

ggplot(genes, aes(log2FC)) + stat_ecdf(geom="step") + ylab("% genes") + xlab("log2FC") +
    geom_hline(yintercept=idx, col="red") + geom_point(aes(num_first_gene, genes[num_first_gene, 2], col="blue")) +
    geom_point(aes(num_last_gene, genes[num_last_gene, 2], col="blue"))
ggsave(paste0(outdir, .Platform$file.sep, "fc_cutoff_finding_one_percent_siggenes_ontop_CUMULATIVE.pdf"), width=10, height=10)




num_last_gene <- min(num_sig_genes + num_five_percent_nonsig, nrow(genes))
num_first_gene <- max(num_sig_genes - num_five_percent_nonsig, 1)

ggplot(genes, aes(log2FC)) + stat_ecdf(geom="step") + ylab("% genes") + xlab("log2FC") +
    geom_hline(yintercept=idx, col="red") + geom_point(aes(num_first_gene, genes[num_first_gene, 2], col="blue")) +
    geom_point(aes(num_last_gene, genes[num_last_gene, 2], col="blue")) + scale_color_identity()
ggsave(paste0(outdir, .Platform$file.sep, "fc_cutoff_finding_five_percent_nonsiggenes_ontop_CUMULATIVE.pdf"), width=10, height=10)


