
#/usr/bin/Rscript

library(ggplot2)

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

num_last_gene <- min(num_sig_genes + num_extend, nrow(genes))
num_diff_to_center <- num_last_gene - num_sig_genes
num_first_gene <- max(num_sig_genes - num_extend, 1)


p1 <- ggplot(genes, aes(FDR)) + stat_ecdf(geom="step") + ylab("% genes") + xlab("FDR") +
	geom_vline(xintercept=fdrcutoff, col="red") +
	geom_point(aes(genes[num_first_gene, 2], num_first_gene/nrow(genes), col="blue")) +
	geom_point(aes(genes[num_last_gene, 2], num_last_gene/nrow(genes), col="blue")) +
	scale_color_identity()
#ggsave(paste0(outdir, .Platform$file.sep, "fdr_cutoff_finding_CUMULATIVE.pdf"), width=10, height=10)




p <- ggplot(genes, aes(FDR)) + stat_ecdf(geom="step") + ylab("% genes") + xlab("FDR") +
	geom_vline(xintercept=fdrcutoff, col="red") +
	scale_color_identity() +
	geom_rect(aes(xmin=genes[num_first_gene, 2], xmax=genes[num_last_gene, 2], ymin=0,
				  ymax=1), fill="pink", alpha=0.02) +
	coord_cartesian(xlim=c(0,0.05)) +
	geom_point(aes(genes[num_first_gene, 2], num_first_gene/nrow(genes), col="darkblue")) +
	geom_point(aes(genes[num_last_gene, 2], num_last_gene/nrow(genes), col="darkblue"))

y1 <- layer_data(p)[which(layer_data(p)$x==genes[num_first_gene, 2]),1]
y2 <- layer_data(p)[which(layer_data(p)$x==genes[num_last_gene, 2]),1]

p2 <- ggplot(genes, aes(FDR)) + stat_ecdf(geom="step") + ylab("% genes") + xlab("FDR") +
	scale_color_identity() +
	geom_rect(aes(xmin=genes[num_first_gene, 2], xmax=genes[num_last_gene, 2], ymin=0,
				  ymax=1), fill="pink", alpha=0.02) +
	geom_point(aes(genes[num_first_gene, 2], y1, col="darkblue")) +
	geom_point(aes(genes[num_last_gene, 2], y2, col="darkblue")) +
	geom_vline(xintercept=fdrcutoff, col="red") +
	coord_cartesian(xlim=c(0,0.05))





num_last_gene <- min(num_sig_genes + num_one_percent_sig, nrow(genes))
num_first_gene <- max(num_sig_genes - num_one_percent_sig, 1)

p <- ggplot(genes, aes(FDR)) + stat_ecdf(geom="step") + ylab("% genes") + xlab("FDR") +
	geom_vline(xintercept=fdrcutoff, col="red") +
	geom_point(aes(genes[num_first_gene, 2], num_first_gene/nrow(genes), col="blue")) +
	geom_point(aes(genes[num_last_gene, 2], num_last_gene/nrow(genes), col="blue")) +
	scale_color_identity()

y1 <- layer_data(p)[which(layer_data(p)$x==genes[num_first_gene, 2]),1]
y2 <- layer_data(p)[which(layer_data(p)$x==genes[num_last_gene, 2]),1]


p3 <-ggplot(genes, aes(FDR)) + stat_ecdf(geom="step") + ylab("% genes") + xlab("FDR") +
	geom_vline(xintercept=fdrcutoff, col="red") +
	geom_point(aes(genes[num_first_gene, 2], y1, col="blue")) +
	geom_point(aes(genes[num_last_gene, 2], y2, col="blue")) +
	scale_color_identity() + coord_cartesian(xlim=c(0,0.05))

#p3 <- ggplot(genes, aes(FDR)) + stat_ecdf(geom="step") + ylab("% genes") + xlab("FDR") +
 #   geom_hline(yintercept=fdrcutoff, col="red") + geom_point(aes(num_first_gene, genes[num_first_gene, 2], col="blue")) +
 #   geom_point(aes(num_last_gene, genes[num_last_gene, 2], col="blue")) +
#	scale_color_identity()
#ggsave(paste0(outdir, .Platform$file.sep, "fdr_cutoff_finding_one_percent_siggenes_ontop_CUMULATIVE.pdf"), width=10, height=10)




num_last_gene <- min(num_sig_genes + num_five_percent_nonsig, nrow(genes))
num_first_gene <- max(num_sig_genes - num_five_percent_nonsig, 1)

p <- ggplot(genes, aes(FDR)) + stat_ecdf(geom="step") + ylab("% genes") + xlab("FDR") +
	geom_vline(xintercept=fdrcutoff, col="red") +
	geom_point(aes(genes[num_first_gene, 2], num_first_gene/nrow(genes), col="blue")) +
	geom_point(aes(genes[num_last_gene, 2], num_last_gene/nrow(genes), col="blue")) +
	scale_color_identity()

y1 <- layer_data(p)[which(layer_data(p)$x==genes[num_first_gene, 2]),1]
y2 <- layer_data(p)[which(layer_data(p)$x==genes[num_last_gene, 2]),1]

p4 <- ggplot(genes, aes(FDR)) + stat_ecdf(geom="step") + ylab("% genes") + xlab("FDR") +
	geom_vline(xintercept=fdrcutoff, col="red") +
	geom_point(aes(genes[num_first_gene, 2], num_first_gene/nrow(genes), col="blue")) +
	geom_point(aes(genes[num_last_gene, 2], num_last_gene/nrow(genes), col="blue")) +
	scale_color_identity()

plots <- list(p1,p2,p3,p4)
ggsave(paste0(outdir, .Platform$file.sep, "fdr_cutoff_finding_all_CUMULATIVE.pdf"), marrangeGrob(grobs = plots, nrow=1, ncol=1), device = "pdf")

#ggsave(paste0(outdir, .Platform$file.sep, "fdr_cutoff_finding_five_percent_nonsiggenes_ontop_CUMULATIVE.pdf"), width=10, height=10)

