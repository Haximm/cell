###################################05.Assigning cell type identity to clusters###################################
pbmc <- readRDS("pbmc.rds")
levels(pbmc)
#annotation
new.cluster.ids <- c("Immune cells_1","Immune cells_2","Basal","Spinous","Mitotic","Melanocytes")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
table(Idents(pbmc))
levels(pbmc)
#TSNE 
pbmc <- RunTSNE(object = pbmc, dims = 1:pcSelect)             
TSNEPlot(object = pbmc, pt.size = 2, label = TRUE)    
write.table(pbmc$seurat_clusters,file="05.NTSNECluster.txt",quote=F,sep="\t",col.names=F)
##UMAP
pbmc <- RunUMAP(object = pbmc, dims = 1:pcSelect)             
UMAPPlot(object = pbmc, pt.size = 2, label = TRUE)    
write.table(pbmc$seurat_clusters,file="05.NUMAPCluster.txt",quote=F,sep="\t",col.names=F)
#DoHeatmap
DoHeatmap(object = pbmc, features = top10$gene) + NoLegend()
saveRDS(pbmc, file="Npbmc.rds")