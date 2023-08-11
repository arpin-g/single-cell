library(Seurat)
library(SeuratData)

###load Seurat object
seurat <- readRDS()
DefaultAssay(seurat) <- "RNA"

### cluster labels celltype X Louvain X group ###
labels <- paste0(seurat$cells_Str.pop, "_", seurat@$seurat_clusters, "_c", seurat$group)
names(labels) <- colnames(seurat)


labels <- labels[ which( labels %in% names( which( table(labels)>=5 ) ) ) ]


### subclustering – determine number of subclusters ###
labels.clusterNo <- ceiling( sort( table(labels) ) / 20 )


o <- order( sapply(strsplit(names(labels.clusterNo), "_" ), function(x) as.numeric(substr(x[2],2,nchar(x[2])))))
labels.clusterNo <- labels.clusterNo[o]


### subclustering - do it! - part 1 ###
metacell.labels <- rep(NA,ncol(seurat)) 
names(metacell.labels) <- colnames(seurat)
metacell.labels[names(labels)] <- labels
seurat[["metacellLabels"]] <- metacell.labels
seurat[["cellInMetacell"]] <- !is.na(metacell.labels)
metacell.data <- matrix(NA,nrow(seurat),0,dimnames=list(rownames(seurat),c()))


### subclustering - do it! – part 2 ###
pb <-txtProgressBar(min = 0, max = length(labels.clusterNo),style=3) 

for( x in names(labels.clusterNo) ) {
  mc.cells <- names(which(metacell.labels==x))
  expr <- data.matrix( seurat@assays$RNA@data[ , mc.cells ] )
  km <- kmeans(t(expr),centers=labels.clusterNo[x])
  lab <- paste0( x, "_k", seq(max(km$cluster)) )
  expr <- t(km$centers)
  colnames(expr) <- lab
  metacell.data <- cbind(metacell.data, expr)
  seurat$metacellLabels[mc.cells] <- paste0( seurat$metacellLabels[mc.cells], " x", km$cluster )
  setTxtProgressBar( pb, pb$getVal()+1 )
}

pb$kill()

saveRDS(list(seurat, metacell.data), file = "/../metacells.rds")
