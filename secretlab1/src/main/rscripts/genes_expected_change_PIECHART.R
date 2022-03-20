#/usr/bin/Rscript

library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
	print("no arguments provided")
}

genes <- args[2]
expchange <- args[3]
outdir <- args[4]
df <- as.data.frame(cbind(genes, expchange))
df_sum <- as.data.frame(table(expchange))

ggplot(df_sum, aes(x="", y=set, fill=expchange)) + geom_bar(stat="identity", width=1) + coord_polar("y", start=0)
 + ylab("") + xlab("")
ggsave(paste0(outdir, .Platform$file.sep, "genes_expected_change_PIECHART.pdf"), width=10, height=10)
