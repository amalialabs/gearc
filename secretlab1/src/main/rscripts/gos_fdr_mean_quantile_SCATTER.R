#/usr/bin/Rscript

library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
	print("ERROR: no arguments provided")
}


gos <- read.csv(args[1], sep="\t", header=FALSE)
quantile <- as.numeric(args[2])
outdir <- args[3]


num_cols <- ncol(gos)
c <- paste0("V", (round(((num_cols-2)*quantile))+2))

gos$mean <- apply(gos[,3:num_cols], 1, mean)

ggplot(gos, aes_string(x="mean", y=c)) + geom_jitter() +
ylab(paste0("-log10(", quantile , " FDR quantile)")) + xlab("-log10(mean FDR)") + geom_abline(col="red", lty=2) +
scale_x_continuous(trans="log10") + scale_y_continuous(trans="log10")
ggsave(paste0(outdir, .Platform$file.sep, "gos_mean_vs_quantile_fdr_SCATTER.pdf"), width=10, height=10)