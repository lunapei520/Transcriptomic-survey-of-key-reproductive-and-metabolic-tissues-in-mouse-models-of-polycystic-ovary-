
# Load WGCNA and flashClust libraries every time you open R
library(WGCNA)
library(flashClust)
library(tidyverse)
library(magrittr) 
library(DESeq2)
setwd("~/Desktop/WGCNA/")
datExpr = read.table("counts_Ovary.txt", header = T, sep = "\t", row.names = 1)
datTraits = read.table("Ovaryfourmodel1.txt", header = T)
condition = as.factor(datTraits$condition)
colData <- data.frame(name=colnames(datExpr),condition)
dds <- DESeqDataSetFromMatrix(round(datExpr),
                              colData = colData,
                              design = ~condition)
dds <- DESeq(dds)
vsd <- varianceStabilizingTransformation(dds)
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
library(genefilter) 
wpn_vsd <- getVarianceStabilizedData(dds)
expr_normalized <- wpn_vsd[ rv_wpn > q75_wpn, ]
expr_normalized[1:5,1:10]
datExpr = t(expr_normalized)
sampleTree = hclust(dist(datExpr, method = "euclidean"), method = "complete")
par(cex =(0.6));
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Dendrogram cluster", sub="", xlab="", cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)
library(factoextra)
fviz_dend(x = sampleTree, cex = 0.7, lwd = 0.7)
traitColors = numbers2colors(datTraits$condition, signed = FALSE);
plotDendroAndColors(sampleTree2, traitColors,groupLabels = names(datTraits),main = "Sample dendrogram and trait heatmap")
powers = c(c(1:10), seq(from = 12, to = 20, by = 2))

sft = pickSoftThreshold(
  datExpr,             # <= Input data
  #blockSize = 30,
  powerVector = powers,
  verbose = 5
)
par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence")
)
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red"
)
abline(h = 0.90, col = "red")
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")

picked_power = 12
temp_cor <- cor       
cor <- WGCNA::cor  

netwk <- blockwiseModules(datExpr,                                          power = picked_power,                
                          networkType = "signed",
                          
                          # == Tree and Block Options ==
                          deepSplit = 2,
                          pamRespectsDendro = F,
                          # detectCutHeight = 0.75,
                          minModuleSize = 30,
                          #maxBlockSize = 4000,
                          
                          # == Module Adjustments ==
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          
                         
                          saveTOMs = T,
                          saveTOMFileBase = "ER",
                          
                          # == Output Options
                          numericLabels = T,
                          verbose = 3)







