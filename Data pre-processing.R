##########################01.Data pre-processing ###################################
#take 'PU' as an example
#Pakages
library(limma)
library(Seurat)
library(dplyr)
library(patchwork)
library(magrittr)
logFCfilter=1               
adjPvalFilter=0.05          
inputFile="PU.txt"       

#Read file and reorganize
rt=read.table(inputFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
conData=avereps(data)


#Initialize the Seurat object with the raw (non-normalized data).
pbmc=CreateSeuratObject(counts = conData,project = "seurat", min.cells=3, min.features=200, names.delim = "_")
table(Idents(pbmc))
#The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^MT-")
#Cancel science and technology units
options(scipen=200)
#Visualize QC metrics as a violin plot
VlnPlot(object = pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
pbmc=subset(x = pbmc, subset = nFeature_RNA > 200 & percent.mt < 5)   
table(Idents(pbmc))
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt",pt.size=1.5)
plot2 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",pt.size=1.5)
CombinePlots(plots = list(plot1, plot2))
#Identification of highly variable features (feature selection)
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
#extracte Genes with high intercellular coefficient of variation
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 2000)
#cell count
table(Idents(pbmc))
# Identify the 10 most highly variable genes
top10 <- head(x = VariableFeatures(object = pbmc), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(object = pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))


###################################02.PCA Principal Component Analysis###################################
#Scaling the data
pbmc=ScaleData(pbmc)          
#Perform linear dimensional reduction
pbmc=RunPCA(object= pbmc,npcs = 20,pc.genes=VariableFeatures(object = pbmc))     
#visualizing both cells and features that define the PCA
VizDimLoadings(object = pbmc, dims = 1:4, reduction = "pca",nfeatures = 20)
#DimPlot
DimPlot(object = pbmc, reduction = "pca")
#DimHeatmap
DimHeatmap(object = pbmc, dims = 1:4, cells = 381, balanced = TRUE,nfeatures = 30,ncol=2)
#Determine the ‘dimensionality’ of the dataset
pbmc <- JackStraw(object = pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(object = pbmc, dims = 1:20)
JackStrawPlot(object = pbmc, dims = 1:15)
#ElbowPlot
ElbowPlot(pbmc)
ElbowPlot(object = pbmc)
###################################03.Cluster the cells and markers gene###################################
##Run non-linear dimensional reduction (UMAP/tSNE)
pcSelect=15
pbmc <- FindNeighbors(object = pbmc, dims = 1:pcSelect)       

#treeplot，determine 'resolution'
library(clustree)
pbmc <- FindClusters(pbmc, resolution = 0.1)
pbmc <- FindClusters(pbmc, resolution = 0.2)
pbmc <- FindClusters(pbmc, resolution = 0.4)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- FindClusters(pbmc, resolution = 0.6)
pbmc <- FindClusters(pbmc, resolution = 0.7)
pbmc <- FindClusters(pbmc, resolution = 0.8)
pbmc <- FindClusters(pbmc, resolution = 1)
pbmc <- FindClusters(pbmc, resolution = 1.2)
clustree(pbmc)
#
pbmc <- FindClusters(object = pbmc, resolution = 0.8)         
table(Idents(pbmc))
#TSNE
pbmc <- RunTSNE(object = pbmc, dims = 1:pcSelect)            
TSNEPlot(object = pbmc, pt.size = 2, label = TRUE)   
write.table(pbmc$seurat_clusters,file="03.tsneCluster.txt",quote=F,sep="\t",col.names=F)
table(Idents(pbmc))
##UMAP
pbmc <- RunUMAP(object = pbmc, dims = 1:pcSelect)             
UMAPPlot(object = pbmc, pt.size = 2, label = TRUE)    
write.table(pbmc$seurat_clusters,file="03.UMAPCluster.txt",quote=F,sep="\t",col.names=F)
##Finding differentially expressed features (cluster biomarkers)
pbmc.markers <- FindAllMarkers(object = pbmc,
                               only.pos = TRUE,
                               min.pct = 0.25,
                               logfc.threshold = logFCfilter)
sig.markers=pbmc.markers[(abs(as.numeric(as.vector(pbmc.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(pbmc.markers$p_val_adj))<adjPvalFilter),]
write.table(sig.markers,file="03.clusterMarkers.txt",sep="\t",row.names=F,quote=F)

top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.table(top10,file="03.top10markers.txt",sep="\t",row.names=F,quote=F)
#FeaturePlot
top2 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
FeaturePlot(pbmc, features = top2$gene,split.by = pbmc$seurat_clusters)
#save pbmc.markers
save(pbmc.markers,file = 'pbmc.markers.Rdata')
write.table(pbmc.markers,file="03.pbmc.markers.txt",sep="\t",row.names=F,quote=F)
#plotting the top 10 markers for each cluster
DoHeatmap(object = pbmc, features = top10$gene) + NoLegend()
#VlnPlot
VlnPlot(object = pbmc, features = row.names(sig.markers)[1:2])
#The genes that need to be shown can be modified
showGenes=c("KRT5","KRT14","CDH3","KRT1","KRT10","DSG1","CDH1","DSC1","KRT2","IVL","TGM3","DSP","LOR","FLG","SPINK5","TYRP1","PMEL","MLANA","S100B","CD207","CD86","CD74") 
#VlnPlot
VlnPlot(object = pbmc, features = showGenes)
#cluster markers scatterPlot
FeaturePlot(object = pbmc, features = showGenes, cols = c("green", "red"))
#cluster markers bublePlot
cluster10Marker=showGenes
DotPlot(object = pbmc, features = cluster10Marker)
#cluster markers RidgePlot
RidgePlot(pbmc, features =showGenes )
saveRDS(pbmc, file="pbmc.rds")

