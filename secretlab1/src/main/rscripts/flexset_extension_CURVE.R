
#/usr/bin/Rscript

library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
	print("no arguments provided")
}

genes <- read.csv(args[1], sep="\t", header=TRUE)
outdir <- args[2]

genes <- genes[,c(1,4,5,7)]
genes$geneset <- factor(genes$geneset, levels=c("SIG_CORE", "FLEX", "SIGNON_CORE"))

idx_last_flex <- nrow(subset(genes, genes$geneset=="SIG_CORE" | genes$geneset=="FLEX"))

genes <- genes[order(genes$weighted_score, decreasing=TRUE),]

z <- genes[idx_last_flex, 3]
bonus <- 1.0
penalty <- 0.0
idx_current_gene <- idx_last_flex
while (bonus > penalty) {
    idx_current_gene <- idx_current_gene + 1
    bonus <- dnorm(idx_current_gene)
    penalty <- abs(genes[idx_current_gene, 3]-z)
}

idx_profitable_extension <- idx_current_gene - 1


ggplot(data.frame(x=c(0, nrow(genes))), aes(x)) +
  stat_function(fun=dnorm, n=101, args=list(mean=0, sd=1)) + ylab("") +
  scale_y_continuous(breaks=NULL) +
  geom_point(aes(idx_profitable_extension, dnorm(idx_profitable_extension), col="red")) +
  xlab("distance") + ylab("bonus") +
  geom_point(aes(0.2*nrow(genes), dnorm(0.2*nrow(genes)), col="blue"))
ggsave(paste0(outdir, .Platform$file.sep, "flexset_extension_CURVE.pdf"), width=10, height=10)

print("plotted file")


