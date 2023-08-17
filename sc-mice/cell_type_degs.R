library(Seurat)

mice <- readRDS("sc-mice/rds/current.rds")
mice <- SetIdent(mice, value= mice$cell_type)
marker_list <- list()

for(type in unique(mice$cell_type)) 
{ 
marker_list[[type]] <- FindMarkers(mice, ident.1="ko", ident.2 = "wt", group.by = "group", subset.ident = type)
}

data_frame <- do.call(rbind, marker_list)
data_frame$gene <- sub("^[^.]+\\.", "", rownames(data_frame))
data_frame$cell_type <- sub("\\..*", "", rownames(data_frame))  

write.csv(data_frame, "sc-mice/csv/cell_type_degs.csv", row.names= F)
