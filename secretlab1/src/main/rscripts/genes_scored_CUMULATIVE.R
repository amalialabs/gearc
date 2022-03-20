#/usr/bin/Rscript

library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
	print("no arguments provided")
}

genes <- args[2]
scores <- args[3]
outdir <- args[4]
df <- as.data.frame(cbind(genes, scores))
ggplot(df, aes(score)) + stat_ecdf(geom="step") + ylab("% genes") + xlab("weighted score")
ggsave(paste0(outdir, .Platform$file.sep, "genes_weighted_score_CUMULATIVE.pdf"), width=10, height=10)