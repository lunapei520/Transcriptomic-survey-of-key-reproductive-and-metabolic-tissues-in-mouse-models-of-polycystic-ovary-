library(Seurat)
library(dplyr)
library(plyr)
#QC
data = read.table("rawcounts.txt",header = T, row.names = 1, sep = "\t")
#Metadata
qc_star_log <- read.table(file = "fourmodel.txt", header = T,
                          row.names = 1, sep="\t")
data <- CreateSeuratObject(counts = data, project = "QC", min.cells = 3, min.features = 200)  
data <- AddMetaData(object = data, metadata = qc_star_log)
mito_genes <- rownames(data)[grep("mt-",rownames(data))]
head(mito_genes,10)
total_counts_per_cell = colSums(as.matrix( data@assays$RNA@counts))
data$percent_mito <-colSums( as.matrix( data@assays$RNA@counts[mito_genes,])  ) / total_counts_per_cell
data$total_counts_per_cell = total_counts_per_cell
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), ncol = 3, group.by = "class")

selected_c <- WhichCells(data, expression = nFeature_RNA >  3000 &nCount_RNA >200000 & percent_mito<0.05  )
selected_f <- rownames(data)[ Matrix::rowSums(data) > 1]
data.filt <- subset(test.df, features=selected_f, cells=selected_c)
dim(data.filt)
