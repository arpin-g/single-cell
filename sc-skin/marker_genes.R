library(Seurat)
library(dplyr)

# seurat <- SetIdent(seurat, value=seurat$fibro_annotation)

markers <- FindAllMarkers(seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

markers <- markers[markers$p_val < 0.01,]
markers%>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) -> top
write.csv(top, "sc-skin/csv/subpop_markers.csv")
