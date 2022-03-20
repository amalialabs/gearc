#/usr/bin/Rscript

library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
	print("ERROR: no arguments provided")
}

print(args)

gos <- args[2]   #TODO select top 3 enriched gos default if user not defined
fdrs_robust <- args[3]
fdrs_extend <- args[4]
quantil <- args[5]
outdir <- args[6]
df <- as.data.frame(cbind(rep(gos, each=1000), fdrs_robust, fdrs_extend))  #TODO make it variable e.g. if runs=500

for (go in gos) {
    tmp <- subset(df, df$gos==go)
    ggplot(tmp, aes(x=fdrs_extend, y=fdrs_robust)) + geom_jitter() + ylab("robust GOs") + xlab("extended robust GOs")
    ggsave(paste0(outdir, .Platform$file.sep, "robust_vs_extended_", go, "_JITTER.pdf"), width=10, height=10)
}




