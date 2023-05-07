#Same as other models
test.df = read.table("17NFqc.txt",header = T, row.names = 1, sep = "\t")
class=colnames(test.df)
class=ifelse(grepl(pattern = "^MII17B",x = class),'MII17B',class)


class=ifelse(grepl(pattern = "^MII17NF",x = class),'MII17NF',class)

condition <- data.frame(class = class,samples=colnames(test.df),
                        row.names = colnames(test.df))
test.df <- CreateSeuratObject(counts = test.df, project = "MII 17NF", min.cells = 3, min.features = 200)
test.df <- AddMetaData(object = test.df, metadata = condition[colnames(test.df),])
x <- factor(test.df@meta.data$class,c('MII17B','MII17NF'))

test.df@meta.data$class <- x
test.df<- NormalizeData(object = test.df, normalization.method = "LogNormalize", scale.factor = 10000)
test.df <- ScaleData(object = test.df)
test.df<- SetIdent(test.df,value = 'class')

countexp.Seurat<-sc.metabolism.Seurat(obj = test.df, method = "AUCell", imputation = F, ncores = 2, metabolism.type = "KEGG")
DimPlot.metabolism(obj = countexp.Seurat, pathway = "Glycolysis / Gluconeogenesis", dimention.reduction.type = "umap", dimention.reduction.run = F, size = 0.5)

metabolism.matrix <- 
  countexp.Seurat@assays$METABOLISM$score
saveRDS(countexp.Seurat,"countexp.Seurat17NF.rds")
write.table(metabolism.matrix, "metabolism17NFMII.txt", sep = "\t")
gene = as.data.frame(test.df@assays$RNA@data)
NFgene = gene[which(rownames(gene) %in% c("Melk","Acadm","Lonp2","Echs1","Pparg","Auh",
                                          "Diablo","Hspb1","Map3k5",
                                          "Bag5")), ]
write.table(NFgene, "fig5_17NFV.txt", sep = "\t")
