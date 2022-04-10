#/usr/bin/Rscript

library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
	print("ERROR: no arguments provided")
}

print(args)

gos <- read.csv(args[1], sep="\t", header=TRUE)
quantile <- args[2]
c <- paste0("V", 1000*quantile+2)
outdir <- args[3]
gos$mean <- lapply(go[,3:1003], mean)

ggplot(gos, aes(x=-log10(mean), y=-log10(c))) + geom_jitter() +
ylab(paste0("-log10(", quantile , "FDR)")) + xlab("-log10(mean FDR)") + geom_abline(col="red", lty=2)
ggsave(paste0(outdir, .Platform$file.sep, "gos_mean_vs_quantile_fdr_scatter.pdf"), width=10, height=10)