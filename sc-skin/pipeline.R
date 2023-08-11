library(Seurat)

seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^mt-")

VlnPlot(
 seurat,
  features = c("nFeature_RNA", "nCount_RNA","percent.mt" ),
  ncol = 3,
  pt.size = 1
)
seurat <- subset(seurat,
                subset = nFeature_RNA > 200 &
                  nFeature_RNA < 5000 &
                  nCount_RNA < 30000 &
                  percent.mt < 1)


seurat <- NormalizeData(seurat, assay = 'RNA', normalization.method = "LogNormalize")
seurat <- FindVariableFeatures(seurat, selection.method = "vst", assay = 'RNA')
seurat <- ScaleData(seurat, model.use = 'linear')

# plot1 <- VariableFeaturePlot(seurat)
# plot2 <-  LabelPoints(
#   plot = plot1,
#   points =  head(VariableFeatures(seurat), 10),
#   repel = TRUE
# )
# plot1 + plot2

### PCA, UMAP 
seurat <- RunPCA(seurat , npcs = 20)
seurat <- RunUMAP(seurat , assay = 'RNA', reduction = 'pca', dims = 1:20, n.neighbors = 20, n.epochs = 1000,
                   min.dist = .0001, spread = 10)
### clustering
seurat <- FindNeighbors(object = seurat, assay = 'RNA', dims = 1:20)
seurat <- FindClusters(object = seurat, resolution = 0.1)

DimPlot(seurat, reduction = "umap", pt.size = 1.5, label = T) + ggtitle ("Seurat Clusters")

saveRDS(seurat, file = "rds/seurat.rds")




