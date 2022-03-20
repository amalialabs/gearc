#/usr/bin/Rscript

library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
	print("ERROR: no arguments provided")
}

print(args)

gos <- args[2]
fdrs_robust <- args[3]  #per go fdr(rob/extend) of x-quantile
fdrs_extend <- args[4]
quantil <- args[5]
outdir <- args[6]
df <- as.data.frame(cbind(gos, fdrs_robust, fdrs_extend))

ggplot(df, aes(x=fdrs_extend, y=fdrs_robust)) + geom_jitter() + ylab("robust GOs") + xlab("extended robust GOs") +
ggtitle(paste0(quantil, " quantile (enrichment)"))
ggsave(paste0(outdir, .Platform$file.sep, "robust_vs_extended_GOs_JITTER.pdf"), width=10, height=10)



