#/usr/bin/Rscript

library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
	print("no arguments provided")
}

genes <- read.csv(args[1], sep="\t", header=TRUE)
outdir <- args[2]

genes <- subset(genes, genes$is_unclear == false)

genes$geneset <- as.character(genes$geneset)
df <- as.data.frame(table(genes$geneset))
colnames(df) <- c("set", "num")

df$p <- df$num/sum(df$num)

ggplot(df, aes(x="", y=num, fill=set)) + geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0, direction=-1) + xlab("") + ylab("") + theme_void() +
  geom_text(aes(label = scales::percent(p, accuracy = 1)),
            position = position_stack(vjust = 0.5),
            color = "grey20", size = 10)
ggsave(paste0(outdir, .Platform$file.sep, "genes_sets_PIE.pdf"), width=10, height=10)
