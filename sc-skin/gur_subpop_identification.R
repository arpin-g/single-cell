library(Seurat)
library(dplyr)
library(scSorter)

seurat <- readRDS("rds/gur_fibro_metacells.rds")[[1]]
top <- read.csv2("csv/subpop_markers.csv")


anno <- top[c("cluster", "gene")]  
colnames(anno) <- c("Type", "Marker")
anno <- as.data.frame(anno)
anno <- anno[complete.cases(anno), ]


topgenes <- head(VariableFeatures(seurat), 2000)
expr = GetAssayData(seurat)
topgene_filter = rowSums(as.matrix(expr)[topgenes, ]!=0) > ncol(expr)*.1
topgenes = topgenes[topgene_filter]
picked_genes = unique(c(anno$Marker, topgenes))
expr <- expr[rownames(expr) %in% picked_genes, ]
expr <- as.data.frame(expr)

rts <- scSorter(expr, anno)

annotations <- as.character(rts$Pred_Type)
seurat<- AddMetaData(seurat, metadata = annotations, col.name = "fibro_annotation")

# DimPlot(seurat, group.by= "fibro_annotation", pt.size=1.5, label=T, repel = T)

saveRDS(seurat, "rds/gur_fibro.rds")
