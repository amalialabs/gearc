
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

num_sig_genes <- nrow(subset(genes, abs(genes$log2FC)>=fccutoff))
num_nonsig_genes <- nrow(subset(genes, abs(genes$log2FC)<fccutoff))


num_one_percent_ontop <- max(round(0.01 * num_sig_genes), 5)
num_five_percent <- round(0.05 * num_nonsig_genes)

genes <- genes[order(abs(genes$log2FC), decreasing=TRUE),]


idx_last_sig_gene <- num_sig_genes-1

idx_extended <- idx_last_sig_gene + min(num_one_percent_ontop , num_five_percent)
diff_to_center <- abs(fccutoff-abs(genes[idx_extended,2]))
upper_bound <- fccutoff+diff_to_center
lower_bound <- abs(genes[idx_extended, 2])


num_extend <- min(num_one_percent_ontop, num_five_percent)

fcext <- lower_bound
percent_extm <- nrow(subset(genes, genes$log2FC<=-fcext))/nrow(genes)
tmp <- subset(genes, genes$log2FC<=-fcext)
extmfc <- tmp[nrow(tmp), 2]
percent_ext <- nrow(subset(genes, genes$log2FC<=fcext))/nrow(genes)
tmp <- subset(genes, genes$log2FC>=fcext)
extfc <- tmp[nrow(tmp),2]

fcext <- upper_bound
percent_extm2 <- nrow(subset(genes, genes$log2FC<=-fcext))/nrow(genes)
tmp <- subset(genes, genes$log2FC<=-fcext)
extmfc2 <- tmp[nrow(tmp), 2]
percent_ext2 <- nrow(subset(genes, genes$log2FC<=fcext))/nrow(genes)
tmp <- subset(genes, genes$log2FC>=fcext)
extfc2 <- tmp[nrow(tmp),2]

p1 <- ggplot(genes, aes(log2FC)) + stat_ecdf(geom="step") + ylab("% genes") + xlab("log2FC") +
	geom_vline(xintercept=c(-1,1), col="red") +
	geom_point(aes(extmfc, percent_extm, col="blue")) +
	geom_point(aes(extfc, percent_ext, col="blue3")) +
	geom_point(aes(extmfc2, percent_extm2, col="blue")) +
	geom_point(aes(extfc2, percent_ext2, col="blue3")) +
	scale_color_identity() +
	annotate("text", x=-5, y=0.5, label=paste0("new interval =  +-[", extfc, ",", extfc2, "]")) +
	ggtitle(paste0("FC cutoff finding (min of 1% sig-genes and 5% nonsig-genes) (+",
				   num_extend, ")"))

p1b <- ggplot(genes, aes(log2FC)) + stat_ecdf(geom="step") + ylab("% genes") + xlab("log2FC") +
	geom_vline(xintercept=c(-1,1), col="red") +
	geom_point(aes(extmfc, percent_extm, col="blue")) +
	geom_point(aes(extfc, percent_ext, col="blue3")) +
	geom_point(aes(extmfc2, percent_extm2, col="blue")) +
	geom_point(aes(extfc2, percent_ext2, col="blue3")) +
	scale_color_identity() +
	#annotate("text", x=-5, y=0.5, label=paste0("new cutoff = +- ", extfc)) +
	ggtitle(paste0("FC cutoff finding (min of 1% sig-genes and 5% nonsig-genes) (+",
				   num_extend, "), zoomed")) +
	coord_cartesian(xlim=c(round(extmfc, 1)-0.1, round(extfc2, 1)+0.1))

pa <- ggplot(genes, aes(log2FC)) + stat_ecdf(geom="step") + ylab("% genes") + xlab("log2FC") +
	geom_vline(xintercept=-1, col="red") +
	geom_point(aes(extmfc, percent_extm, col="blue")) +
	geom_point(aes(extmfc2, percent_extm2, col="blue")) +
	scale_color_identity() +
	annotate("text", x=-1.05, y=0.5, label=paste0("[", extmfc2, ",", extmfc, "]")) +
	ggtitle("left interval zoomed") +
	coord_cartesian(xlim=c(round(extmfc, 1)-0.1, round(extmfc2, 1)+0.1))

pb <- ggplot(genes, aes(log2FC)) + stat_ecdf(geom="step") + ylab("% genes") + xlab("log2FC") +
	geom_vline(xintercept=1, col="red") +
	geom_point(aes(extfc, percent_ext, col="blue3")) +
	geom_point(aes(extfc2, percent_ext2, col="blue3")) +
	scale_color_identity() +
	annotate("text", x=1.05, y=0.5, label=paste0("[", extfc, ",", extfc2, "]")) +
	ggtitle("right interval zoomed") +
	coord_cartesian(xlim=c(round(extfc, 1)-0.1, round(extfc2, 1)+0.1))

p1x <- grid.arrange(pa, pb, ncol=2, top=paste0("FC cutoff finding
            (min of 1% sig-genes and 5% nonsig-genes) (+", num_extend, ")"))




idx_extended <- idx_last_sig_gene + num_one_percent_ontop
diff_to_center <- abs(fccutoff-abs(genes[idx_extended,2]))
upper_bound <- fccutoff+diff_to_center
lower_bound <- abs(genes[idx_extended, 2])


num_extend <- num_one_percent_ontop

fcext <- lower_bound
percent_extm <- nrow(subset(genes, genes$log2FC<=-fcext))/nrow(genes)
tmp <- subset(genes, genes$log2FC<=-fcext)
extmfc <- tmp[nrow(tmp), 2]
percent_ext <- nrow(subset(genes, genes$log2FC<=fcext))/nrow(genes)
tmp <- subset(genes, genes$log2FC>=fcext)
extfc <- tmp[nrow(tmp),2]

fcext <- upper_bound
percent_extm2 <- nrow(subset(genes, genes$log2FC<=-fcext))/nrow(genes)
tmp <- subset(genes, genes$log2FC<=-fcext)
extmfc2 <- tmp[nrow(tmp), 2]
percent_ext2 <- nrow(subset(genes, genes$log2FC<=fcext))/nrow(genes)
tmp <- subset(genes, genes$log2FC>=fcext)
extfc2 <- tmp[nrow(tmp),2]

p2 <- ggplot(genes, aes(log2FC)) + stat_ecdf(geom="step") + ylab("% genes") + xlab("log2FC") +
	geom_vline(xintercept=c(-1,1), col="red") +
	geom_point(aes(extmfc, percent_extm, col="blue")) +
	geom_point(aes(extfc, percent_ext, col="blue3")) +
	geom_point(aes(extmfc2, percent_extm2, col="blue")) +
	geom_point(aes(extfc2, percent_ext2, col="blue3")) +
	scale_color_identity() +
	annotate("text", x=-5, y=0.5, label=paste0("new interval =  +-[", extfc, ",", extfc2, "]")) +
	ggtitle(paste0("FC cutoff finding with 1% sig-genes on top (+",
				   num_extend, ")"))

p2b <- ggplot(genes, aes(log2FC)) + stat_ecdf(geom="step") + ylab("% genes") + xlab("log2FC") +
	geom_vline(xintercept=c(-1,1), col="red") +
	geom_point(aes(extmfc, percent_extm, col="blue")) +
	geom_point(aes(extfc, percent_ext, col="blue3")) +
	geom_point(aes(extmfc2, percent_extm2, col="blue")) +
	geom_point(aes(extfc2, percent_ext2, col="blue3")) +
	scale_color_identity() +
	#annotate("text", x=-5, y=0.5, label=paste0("new cutoff = +- ", extfc)) +
	ggtitle(paste0("FC cutoff finding with 1% sig-genes on top (+",
				   num_extend, "), zoomed")) +
	coord_cartesian(xlim=c(round(extmfc, 1)-0.1, round(extfc2, 1)+0.1))

pc <- ggplot(genes, aes(log2FC)) + stat_ecdf(geom="step") + ylab("% genes") + xlab("log2FC") +
	geom_vline(xintercept=-1, col="red") +
	geom_point(aes(extmfc, percent_extm, col="blue")) +
	geom_point(aes(extmfc2, percent_extm2, col="blue")) +
	scale_color_identity() +
	annotate("text", x=-1.05, y=0.5, label=paste0("[", extmfc2, ",", extmfc, "]")) +
	ggtitle("left interval zoomed") +
	coord_cartesian(xlim=c(round(extmfc, 1)-0.1, round(extmfc2, 1)+0.1))

pd <- ggplot(genes, aes(log2FC)) + stat_ecdf(geom="step") + ylab("% genes") + xlab("log2FC") +
	geom_vline(xintercept=1, col="red") +
	geom_point(aes(extfc, percent_ext, col="blue3")) +
	geom_point(aes(extfc2, percent_ext2, col="blue3")) +
	scale_color_identity() +
	annotate("text", x=1.05, y=0.5, label=paste0("[", extfc, ",", extfc2, "]")) +
	ggtitle("right interval zoomed") +
	coord_cartesian(xlim=c(round(extfc, 1)-0.1, round(extfc2, 1)+0.1))

p2x <- grid.arrange(pc, pd, ncol=2, top=paste0("FC cutoff finding
            with 1% sig-genes on top (+", num_extend, ")"))






idx_extended <- idx_last_sig_gene + num_five_percent
diff_to_center <- abs(fccutoff-abs(genes[idx_extended,2]))
upper_bound <- fccutoff+diff_to_center
lower_bound <- abs(genes[idx_extended, 2])


num_extend <- num_five_percent

fcext <- lower_bound
percent_extm <- nrow(subset(genes, genes$log2FC<=-fcext))/nrow(genes)
tmp <- subset(genes, genes$log2FC<=-fcext)
extmfc <- tmp[nrow(tmp), 2]
percent_ext <- nrow(subset(genes, genes$log2FC<=fcext))/nrow(genes)
tmp <- subset(genes, genes$log2FC>=fcext)
extfc <- tmp[nrow(tmp),2]

fcext <- upper_bound
percent_extm2 <- nrow(subset(genes, genes$log2FC<=-fcext))/nrow(genes)
tmp <- subset(genes, genes$log2FC<=-fcext)
extmfc2 <- tmp[nrow(tmp), 2]
percent_ext2 <- nrow(subset(genes, genes$log2FC<=fcext))/nrow(genes)
tmp <- subset(genes, genes$log2FC>=fcext)
extfc2 <- tmp[nrow(tmp),2]

p3 <- ggplot(genes, aes(log2FC)) + stat_ecdf(geom="step") + ylab("% genes") + xlab("log2FC") +
	geom_vline(xintercept=c(-1,1), col="red") +
	geom_point(aes(extmfc, percent_extm, col="blue")) +
	geom_point(aes(extfc, percent_ext, col="blue3")) +
	geom_point(aes(extmfc2, percent_extm2, col="blue")) +
	geom_point(aes(extfc2, percent_ext2, col="blue3")) +
	scale_color_identity() +
	annotate("text", x=-5, y=0.5, label=paste0("new interval =  +-[", extfc, ",", extfc2, "]")) +
	ggtitle(paste0("FC cutoff finding with 5% nonsig-genes on top (",
				   num_extend, ")"))

p3b <- ggplot(genes, aes(log2FC)) + stat_ecdf(geom="step") + ylab("% genes") + xlab("log2FC") +
	geom_vline(xintercept=c(-1,1), col="red") +
	geom_point(aes(extmfc, percent_extm, col="blue")) +
	geom_point(aes(extfc, percent_ext, col="blue3")) +
	geom_point(aes(extmfc2, percent_extm2, col="blue")) +
	geom_point(aes(extfc2, percent_ext2, col="blue3")) +
	scale_color_identity() +
	#annotate("text", x=-5, y=0.5, label=paste0("new cutoff = +- ", extfc)) +
	ggtitle(paste0("FC cutoff finding with 5% nonsig-genes on top (+",
				   num_extend, "), zoomed")) +
	coord_cartesian(xlim=c(round(extmfc2, 1)-0.1, round(extfc2, 1)+0.1))

pe <- ggplot(genes, aes(log2FC)) + stat_ecdf(geom="step") + ylab("% genes") + xlab("log2FC") +
	geom_vline(xintercept=-1, col="red") +
	geom_point(aes(extmfc, percent_extm, col="blue")) +
	geom_point(aes(extmfc2, percent_extm2, col="blue")) +
	scale_color_identity() +
	annotate("text", x=-1.05, y=0.5, label=paste0("[", extmfc2, ",", extmfc, "]")) +
	ggtitle("left interval zoomed") +
	coord_cartesian(xlim=c(round(extmfc2, 1)-0.1, round(extmfc, 1)+0.1))

pf <- ggplot(genes, aes(log2FC)) + stat_ecdf(geom="step") + ylab("% genes") + xlab("log2FC") +
	geom_vline(xintercept=1, col="red") +
	geom_point(aes(extfc, percent_ext, col="blue3")) +
	geom_point(aes(extfc2, percent_ext2, col="blue3")) +
	scale_color_identity() +
	annotate("text", x=1.05, y=0.5, label=paste0("[", extfc, ",", extfc2, "]")) +
	ggtitle("right interval zoomed") +
	coord_cartesian(xlim=c(round(extfc, 1)-0.1, round(extfc2, 1)+0.1))

p3x <- grid.arrange(pe, pf, ncol=2, top=paste0("FC cutoff finding
            with 5% nonsig-genes on top (+", num_extend, ")"))


plots <- list(p1, p1b, p1x, p2, p2b, p2x, p3, p3b, p3x)
ggsave(paste0(outdir, .Platform$file.sep, "fc_cutoff_finding_all_CUMULATIVE.pdf"),
	   marrangeGrob(grobs = plots, nrow=1, ncol=1), device = "pdf", width=14, height=12)





if (FALSE) {
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
	scale_color_identity() +
	annotate("text", x=-5, y=0.5, label=paste0("new cutoff = +- ", extfc)) +
	ggtitle(paste0("FC cutoff finding (min of 1% sig-genes and 5% nonsig-genes) (", num_extend, ")"))

p2 <- ggplot(genes, aes(log2FC)) + stat_ecdf(geom="step") + ylab("% genes") + xlab("log2FC") +
	geom_vline(xintercept=c(-1,1), col="red") +
	geom_point(aes(extmfc, percent_extm, col="blue")) +
	geom_point(aes(extfc, percent_ext, col="blue")) +
	scale_color_identity() + coord_cartesian(xlim=c(-2,2)) +
	annotate("text", x=0, y=0.1, label=paste0("new cutoff = +- ", extfc)) +
	ggtitle(paste0("FC cutoff finding (min of 1% sig-genes and 5% nonsig-genes) (", num_extend, "), zoomed"))


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
	scale_color_identity() +
	annotate("text", x=-5, y=0.5, label=paste0("new cutoff = +- ", extfc)) +
	ggtitle(paste0("FC cutoff finding with 1% sig-genes ("), num_extend, ")")

p4 <- ggplot(genes, aes(log2FC)) + stat_ecdf(geom="step") + ylab("% genes") + xlab("log2FC") +
	geom_vline(xintercept=c(-1,1), col="red") +
	geom_point(aes(extmfc, percent_extm, col="blue")) +
	geom_point(aes(extfc, percent_ext, col="blue")) +
	scale_color_identity() + coord_cartesian(xlim=c(-2,2)) +
	annotate("text", x=0, y=0.1, label=paste0("new cutoff = +- ", extfc)) +
	ggtitle(paste0("FC cutoff finding with 1% sig-genes ("), num_extend, "), zoomed")



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
	scale_color_identity() +
	annotate("text", x=-5, y=0.5, label=paste0("new cutoff = +- ", extfc)) +
	ggtitle(paste0("FC cutoff finding with 5% nonsig-genes ("), num_extend, ")")

p6 <- ggplot(genes, aes(log2FC)) + stat_ecdf(geom="step") + ylab("% genes") + xlab("log2FC") +
	geom_vline(xintercept=c(-1,1), col="red") +
	geom_point(aes(extmfc, percent_extm, col="blue")) +
	geom_point(aes(extfc, percent_ext, col="blue")) +
	scale_color_identity() + coord_cartesian(xlim=c(-2,2)) +
	annotate("text", x=0, y=0.1, label=paste0("new cutoff = +- ", extfc)) +
	ggtitle(paste0("FC cutoff finding with 5% nonsig-genes ("), num_extend, "), zoomed")

plots <- list(p1, p2, p3, p4, p5, p6)
ggsave(paste0(outdir, .Platform$file.sep, "fc_cutoff_finding_all_CUMULATIVE.pdf"),
	   marrangeGrob(grobs = plots, nrow=1, ncol=1), device = "pdf")
}
