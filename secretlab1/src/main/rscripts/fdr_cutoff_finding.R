
#/usr/bin/Rscript

library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
	print("no arguments provided")
}

genes <- read.csv(args[1], sep="\t", header=TRUE)
outdir <- args[2]

genes <- genes[,c(1,2,6)]

fdrcutoff <- 0.05
num_sig_genes <- nrow(subset(genes, genes$is_sig == "true"))
num_nonsig_genes <- nrow(subset(genes, genes$is_sig == "false"))
num_one_percent_sig <- max(round(0.01 * num_sig_genes), 5)
num_five_percent_nonsig <- round(0.05 * num_nonsig_genes)

num_extend <- min(num_one_percent_sig, num_five_percent_nonsig)

num_last_gene <- num_sig_genes + num_extend
num_diff_to_center <- num_last_gene - num_sig_genes
num_first_gene <- num_sig_genes - num_extend


ggplot(genes, aes(FDR)) + stat_ecdf(geom="step") + ylab("% genes") + xlab("FDR") +
    geom_hline(yintercept=fdrcutoff, col="red") + geom_point(aes(num_first_gene, genes[num_first_gene, 2])) +
    geom_point(aes(num_last_gene, genes[num_last_gene, 2]))
ggsave(paste0(outdir, .Platform$file.sep, "fdr_cutoff_finding_CUMULATIVE.pdf"), width=10, height=10)
