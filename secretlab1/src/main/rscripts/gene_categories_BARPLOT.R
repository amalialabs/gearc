#/usr/bin/Rscript

library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
	print("no arguments provided")
}

genes <- args[2]
group <- args[3]
outdir <- args[4]
df <- as.data.frame(cbind(genes, group))
#groups: sig-diff, unclear, sig-nondiff, (rest)
#df$group <- factor(df$group, levels=c())
ggplot(df, aes(x=group, fill=group)) + geom_bar() + ylab("# genes") + xlab("")
ggsave(paste0(outdir, .Platform$file.sep, "genes_categories_barplot.pdf"), width=10, height=10)