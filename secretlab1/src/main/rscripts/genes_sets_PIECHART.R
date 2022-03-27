#/usr/bin/Rscript

library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
	print("no arguments provided")
}

genes <- read.csv(args[1], sep="\t", header=TRUE)
outdir <- args[2]

df <- as.data.frame(t(table(genes$geneset)))
colnames(df) <- c("set", "num")

print(df)

ggplot(df, aes(x="", y=set, fill=set)) + geom_bar(stat="identity", width=1) + coord_polar("y", start=0)
 + ylab("") + xlab("")
ggsave(paste0(outdir, .Platform$file.sep, "genes_sets_PIE.pdf"), width=10, height=10)
