reindeer <- read_csv("data/marker_genes/f2_reindeer_markers.csv")
fetal <- read_csv("data/marker_genes/f2_all_markers.csv")

joined_df <- join(reindeer, fetal, by="gene", )


zro <- subset(joined_df, cluster == "0", select = gene)
porc <- joined_df[,c("cluster","gene", "Fetal_Fibroblasts")]
colnames(porc) <- c("cluster","gene", "fc")
write.csv(porc, "data/marker_genes/f2_reindeer_markers.csv")


murine$gene <- sapply(murine$gene, toupper)

joined_df <- join(mer, murine, by="gene", type="inner")
joined_df2 <- joined_df[c("cluster", "gene",    "fc"  ,   "p_val" )]
write.csv2(joined_df2, "/storage2/proj/hmskin_scrna_susie_2022/data/marker_genes/f2_reindeer_with_p_val.csv")


saveRDS(joined_df2, "data/marker_genes/x.rds")
colnames(joined_df)

