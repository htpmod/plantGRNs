library(corrplot)
library(RColorBrewer)
ord <- c('AP1_day2', 'AP1_day4', 'AP1_day8', 'SEP3_day2', 'SEP3_day4', 'SEP3_day8')
gene_table <- read.table('gene_pairwise_overlapping_rate.txt', header=TRUE, row.names=1, check.names=FALSE)
gene_number <- rownames(gene_table)
names(gene_number) <- gsub(',.*', '', gene_number)
gene_table <- gene_table[gene_number[ord], ord]

peak_table <- read.table('peak_pairwise_overlapping_rate.txt', header=TRUE, row.names=1, check.names=FALSE)
peak_number <- rownames(peak_table)
names(peak_number) <- gsub(',.*', '', peak_number)
peak_table <- peak_table[peak_number[ord], ord]

G <- as.matrix(gene_table)
colnames(G) <- rownames(G)
P <- as.matrix(peak_table)
colnames(P) <- rownames(P)
diag(G) <- diag(P) <- 0
GP <- G
GP[lower.tri(GP)] <- -G[lower.tri(G)]
rownames(GP) <- rownames(G)

pdf("corrplot.pdf", width=10, heigh=10, pointsize=10)
corrplot(P, order="original", col=colorRampPalette(brewer.pal(n=8, name="PuOr"))(100), method="pie", tl.col="black", main="Peak")
corrplot(G, order="original", col=colorRampPalette(rev(brewer.pal(n=8, name="PuOr")))(100), method="pie", tl.col="black", main="Gene")
corrplot(GP, method="pie", col=colorRampPalette(brewer.pal(n=8, name="PuOr"))(100), tl.col="black", main="Gene and peak")
dev.off()
