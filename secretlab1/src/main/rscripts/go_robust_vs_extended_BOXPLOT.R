#/usr/bin/Rscript

library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
	print("ERROR: no arguments provided")
}

print(args)

gos <- args[2]   #TODO default select top 3 enriched gos if user not defined
fdrs_robust <- args[3]
fdrs_extend <- args[4]
quantil <- args[5]
outdir <- args[6]
tmp <- as.data.frame(cbind(rep(gos, each=1000), fdrs_robust, fdrs_extend))

for (go in gos) {
    df <- subset(tmp, tmp$gos==go)
    df2 <- reshape2::melt(df, id.vars=1, variable.name="type")

    ggplot(df2, aes(x=type, y=value, fill=type)) + geom_boxplot() + ylab("FDRs from robust runs") + xlab("") +
    ggtitle(paste0(quantil, " quantile (enrichment) for GO node ", go))
    ggsave(paste0(outdir, .Platform$file.sep, "robust_vs_extended_GO_", go, "_BOXPLOT.pdf"), width=10, height=10)
}