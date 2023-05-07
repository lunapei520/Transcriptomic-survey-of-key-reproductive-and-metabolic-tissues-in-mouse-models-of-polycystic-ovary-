# Same code as Ovary, Hypothalamus, adipose tissue and oocyte
library(DESeq2)
read <- read.table("countsqc.txt", seq="/t", row.names=1)
read <- select(read, contains("Ovary"))
colData <- data.frame(name=colnames(read),condition)
dds <- DESeqDataSetFromMatrix(as.matrix(read),
                              colData = colData,
                              design = ~ condition)
dds2 <- DESeq(dds)
keep <- rowSums(counts(dds2)) >= 5
dds2 <- dds2[keep,]
resultsNames(dds2)
res <- results(dds2)
sum( res$pvalue < 0.05, na.rm=TRUE )
table( is.na(res$pvalue) )
resSig <- res[ which(res$pvalue < 0.05), ]
head(resSig)
write.table(resSig, "DeseqOvary.txt", 
            sep = "\t", row.names = TRUE)

