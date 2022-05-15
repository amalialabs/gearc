
#/usr/bin/Rscript

library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
	print("no arguments provided")
}

genes <- read.csv(args[1], sep="\t", header=TRUE)
outdir <- args[2]

genes <- subset(genes, genes$is_unclear == "false")

genes <- genes[,c(1,4,5,7)]
genes$geneset <- factor(genes$geneset, levels=c("SIG_CORE", "FLEX", "SIGNON_CORE"))

idx_last_flex <- nrow(subset(genes, genes$geneset=="FLEX"))

genes <- genes[order(genes$weighted_score, decreasing=TRUE),]

nsig <- nrow(subset(genes, genes$geneset=="SIG_CORE"))
nflex <- 0.2*nrow(subset(genes, genes$geneset=="FLEX"))


#SD <- 3000
SD <- 0.8*nrow(genes)/6
gauss_mean <- nsig

z <- genes[idx_last_flex, 3] + 0.000001
bonus <- 1.0
penalty <- 0.0
idx_current_gene <- idx_last_flex
while (bonus > penalty && idx_current_gene < nrow(genes)) {
    idx_current_gene <- idx_current_gene + 1
    bonus <- dnorm(idx_current_gene, gauss_mean, SD)
    penalty <- abs(genes[idx_current_gene, 3]-z)  * (idx_current_gene-idx_last_flex)
}

idx_profitable_extension <- idx_current_gene - 1

ggplot(data.frame(x=c(0, nrow(genes))), aes(x)) +
  stat_function(fun=dnorm, n=101, args=list(mean=gauss_mean, sd=SD)) +
  xlab("distance") + ylab("bonus") +
  geom_point(aes(nsig, dnorm(nsig, gauss_mean, SD), col="blue")) +
  scale_color_identity() +
  geom_point(aes(nflex+nsig, dnorm(nflex+nsig, gauss_mean, SD), col="yellow")) +
  geom_point(aes(nsig+idx_last_flex, dnorm(nsig+idx_last_flex, gauss_mean, SD), col="orange")) +
  geom_point(aes(idx_profitable_extension, dnorm(idx_profitable_extension, gauss_mean, SD), col="red"))

ggsave(paste0(outdir, .Platform$file.sep, "flexset_extension_CURVE.pdf"), width=10, height=10)


