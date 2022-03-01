###################################0.data preparation###################################
library(limma)               
ReadFile="rowread.txt"        
NSFile="nsnames.txt"             
AWFile="awnames.txt"          
PUFile="punames.txt"           

#Read the rowread file and reorganize
rt=read.table(ReadFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
#Read sample information
sample1=read.table(NSFile,sep="\t",header=F,check.names=F)
sample2=read.table(AWFile,sep="\t",header=F,check.names=F)
sample3=read.table(PUFile,sep="\t",header=F,check.names=F)
sampleName1=gsub("^ | $", "", as.vector(sample1[,1]))
sampleName2=gsub("^ | $", "", as.vector(sample2[,1]))
sampleName3=gsub("^ | $", "", as.vector(sample3[,1]))
NSData=data[,sampleName1]
AWData=data[,sampleName2]
PUData=data[,sampleName3]
#output NS samples expression levels
NSData=rbind(id=colnames(NSData), NSData)
write.table(NSData, file="NS.txt",sep="\t",quote=F,col.names=F)

#output AW,PU samples expression levels
AWData=rbind(id=colnames(AWData), AWData)
write.table(AWData, file="AW.txt",sep="\t",quote=F,col.names=F)
PUData=rbind(id=colnames(PUData), PUData)
write.table(PUData, file="PU.txt",sep="\t",quote=F,col.names=F)

##########################01.Data pre-processing and correctio###################################
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
###################################03.Cluster the cells和marker基因###################################
##Run non-linear dimensional reduction (UMAP/tSNE)
pcSelect=15
pbmc <- FindNeighbors(object = pbmc, dims = 1:pcSelect)       #计算邻接距离

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

###################################04.cluster Markers gene GO/KEGG analysis################################################
#use compareCluster function in clusterProfiler
load(file = 'pbmc.markers.Rdata')
table(pbmc.markers$cluster)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
ids=bitr(pbmc.markers$gene,'SYMBOL','ENTREZID','org.Hs.eg.db')
write.table(ids,file="03.ids.txt",sep="\t",row.names=F,quote=F)
pbmc.markers=merge(pbmc.markers,ids,by.x='gene',by.y='SYMBOL')
View(pbmc.markers)

gcSample=split(pbmc.markers$ENTREZID, pbmc.markers$cluster)
gcSample 
## KEGG analysis
xx <- compareCluster(gcSample, fun="enrichKEGG",
                     organism="hsa", pvalueCutoff=0.05)
p=dotplot(xx) 
p+ theme(axis.text.x = element_text(angle = 45, 
                                    vjust = 0.5, hjust=0.5))
p
## GO analysis
xx <- compareCluster(gcSample,
                     fun = "enrichGO",
                     OrgDb = "org.Hs.eg.db",
                     ont = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.01,
                     qvalueCutoff = 0.05
)
p <- dotplot(xx)
p + theme(axis.text.x = element_text(
  angle = 45,
  vjust = 0.5, hjust = 0.5
))

# ANALYSE BY cluster seperately
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")

pvalueFilter=0.05       
qvalueFilter=0.05       
inputFile="03.clusterMarkers.txt"      
#read file
rt=read.table(inputFile, header=T, sep="\t", check.names=F)
#Separate data sets by cluster 
cluster0 <- rt[rt$cluster == "0",]
cluster1 <- rt[rt$cluster == "1",]
cluster2 <- rt[rt$cluster == "2",]
cluster3 <- rt[rt$cluster == "3",]
cluster4 <- rt[rt$cluster == "4",]
cluster5 <- rt[rt$cluster == "5",]


#Define Colors
colorSel="qvalue"
if(qvalueFilter>0.05){
  colorSel="pvalue"
}
#cluster0/1/2/3/4/5/....
#Gene name is converted to gene ID
genes=unique(as.vector(cluster0[,7]))#cluster0/1/2/3/4/5/....
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]        
#gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)
##GO analysis
kk=enrichGO(gene=gene,OrgDb=org.Hs.eg.db, pvalueCutoff=1, qvalueCutoff=1, ont="all", readable=T)
GO=as.data.frame(kk)
GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]

write.table(GO, file="04.GOcluster0.txt", sep="\t", quote=F, row.names = F)
#Defines the number of GO to display
showNum=10
if(nrow(GO)<30){
  showNum=nrow(GO)
}
#histogram
bar=barplot(kk, drop=TRUE, showCategory=showNum, split="ONTOLOGY", color=colorSel) + facet_grid(ONTOLOGY~., scale='free_y')
print(bar)
#bubleplot
bub=dotplot(kk, showCategory=showNum, orderBy="GeneRatio", split="ONTOLOGY", color=colorSel) + facet_grid(ONTOLOGY~., scale='free_y')
print(bub)
##kegg analysis
kk=enrichKEGG(gene=gene, organism="hsa", pvalueCutoff=1, qvalueCutoff=1)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$genes[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]
#save results
write.table(KEGG, file="04.KEGGcluster0.txt", sep="\t", quote=F, row.names = F)
#Defines the number of display paths
showNum=30
if(nrow(KEGG)<showNum){
  showNum=nrow(KEGG)
}
#histogram
barplot(kk, drop=TRUE, showCategory=showNum, color=colorSel)
#bubleplot
dotplot(kk, showCategory=showNum, orderBy="GeneRatio", color=colorSel)
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
###################################06.cell communication analysis###############################
library(CellChat)
library(ggplot2)                  
library(patchwork)
options(stringsAsFactors = FALSE)
library(igraph)
library(Seurat)
##Data input & processing and initialization of CellChat object
PU<- readRDS(file = "Npbmc.rds")
data.input <- GetAssayData(PU, assay = "RNA", slot = "data") # normalized data matrix
labels <- Idents(PU)
identity <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels
unique(identity$group) # check the cell labels
#set up cellchat object
cellchat.PU <- createCellChat(object = data.input)
cellchat.PU <- addMeta(cellchat.PU, meta = identity, meta.name = "labels")
cellchat.PU <- setIdent(cellchat.PU, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat.PU@idents) # show factor levels of the cell labels
#
cellchat.PU@data.signaling
cellchat=cellchat.PU
#Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
# use all CellChatDB for cell-cell communication analysis
#CellChatDB.use <- CellChatDB # simply use the default CellChatDB
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
cellchat@DB <- CellChatDB.use
#Preprocessing of expression data for cell communication analysis
cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
future::plan("multiprocess", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
##Part II: Inference of cell-cell communication network
#Compute the communication probability and infer cellular communication network
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
#Extract the inferred cellular communication network as a data frame
#df.net <- subsetCommunication(cellchat)
#Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)
#Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)
# visualize the aggregated cell-cell communication network
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
#Due to the complicated cell-cell communication network, we can examine the signaling sent from each cell group
mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
#Part III: Visualization of cell-cell communication network
#All the signaling pathways showing significant communications can be accessed by cellchat@netP$pathways.
cellchat@netP$pathways
write.table(cellchat@netP$pathways,file="PU.toppathways.txt",quote=F,sep="\t",col.names=F)

#Visualize each signaling pathway using Hierarchy plot, Circle plot or Chord diagram
#Here we take input of one signaling pathway as an example
pathways.show <- c("EGF") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
#Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
#Chord plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
#> Note: The first link end is drawn out of sector 'Inflam. FIB'.
#heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
#> Do heatmap based on a single object
#Compute the contribution of each ligand-receptor pair to the overall signaling pathway and visualize cell-cell communication mediated by a single ligand-receptor pair
netAnalysis_contribution(cellchat, signaling = pathways.show)
#Chord plot
pairLR.EGF <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.EGF[1,] # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
#Circle plot
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
#Chord diagram
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
#Automatically save the plots of the all inferred network for quick exploration
# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat@netP$pathways
pathways.show.all
length(pathways.show.all)
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)
table(cellchat@idents)
vertex.receiver = seq(1,2)
vertex.receiver
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 4, height = 2, units = 'in', dpi = 300)
}
#Visualize cell-cell communication mediated by multiple ligand-receptors or signaling pathways
#doplot
# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_bubble(cellchat, sources.use =1, targets.use = c(1:5), remove.isolate = FALSE)
#> Comparing communications on a single object
# show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_bubble(cellchat, sources.use = 2, targets.use = c(1:5), signaling = c("WNT","EGF"), remove.isolate = FALSE)
#> Comparing communications on a single object
# show all the significant interactions (L-R pairs) based on user's input (defined by `pairLR.use`)
pairLR.use <- extractEnrichedLR(cellchat, signaling = c("EGF"))
netVisual_bubble(cellchat, sources.use = 2, targets.use = c(1:5), pairLR.use = pairLR.use, remove.isolate = TRUE)


#> Comparing communications on a single object
#Chord diagram
# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
# show all the interactions sending from Inflam.FIB

netVisual_chord_gene(cellchat, sources.use = c(1,2), targets.use = c(1:5), lab.cex = 0.5,legend.pos.y = 30)

#> Note: The first link end is drawn out of sector 'MIF'.
# show all the interactions received by Inflam.DC
netVisual_chord_gene(cellchat, sources.use = c(1,2,3), targets.use = 4, legend.pos.x = 15)
# show all the significant interactions (L-R pairs) associated with certain signaling pathways
pathways.show <- c("CALCR")
pdf(file="L-R pairschord by CALCR.pdf",width=9,height=7)
netVisual_chord_gene(cellchat, sources.use = 2, targets.use = c(1:5), signaling = pathways.show,legend.pos.x = 8)
dev.off()

#> Note: The second link end is drawn out of sector 'CXCR4 '.
#> Note: The first link end is drawn out of sector 'CXCL12 '.
# show all the significant signaling pathways from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4), targets.use = c(1:4), slot.name = "netP", legend.pos.x = 10)
#> Note: The second link end is drawn out of sector ' '.
#> Note: The first link end is drawn out of sector 'MIF'.
#> Note: The second link end is drawn out of sector ' '.
#> Note: The first link end is drawn out of sector 'CXCL '.
#Plot the signaling gene expression distribution using violin/dot plot
plotGeneExpression(cellchat, signaling = "EGF")
plotGeneExpression(cellchat, signaling = "WNT", enriched.only = FALSE)
pathways.show <- c("VISFATIN")
plotGeneExpression(cellchat, signaling = pathways.show, enriched.only = FALSE)
#Part IV: Systems analysis of cell-cell communication network
#Identify signaling roles (e.g., dominant senders, receivers) of cell groups as well as the major contributing signaling
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
pathways.show <- c("NT")
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
#Visualize the dominant senders (sources) and receivers (targets) in a 2D space
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("WNT","EGF"))
#> Signaling role analysis on the cell-cell communication network from user's input
gg1 + gg2
#Identify signals contributing most to outgoing or incoming signaling of certain cell groups
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2
# Signaling role analysis on the cell-cell communication networks of interest
ht <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("EGF"))

ht <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("WNT"))


#Identify global communication patterns to explore how multiple cell types and signaling pathways coordinate together
#Identify and visualize outgoing communication pattern of secreting cells
#Loading required packages
library(NMF)
library(ggalluvial)
#Here we run selectK to infer the number of patterns.
pdf(file="selectK_outgoingpattern.pdf",width=11,height=5)
selectK(cellchat, pattern = "outgoing")
dev.off()
#Both Cophenetic and Silhouette values begin to drop suddenly when the number of outgoing patterns is 3.

nPatterns = 3
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)

# river plot

netAnalysis_river(cellchat, pattern = "outgoing")

#> Please make sure you have load `library(ggalluvial)` when running this function
# dot plot

netAnalysis_dot(cellchat, pattern = "outgoing")

#Incoming patterns
#Here we run selectK to infer the number of patterns.

selectK(cellchat, pattern = "incoming")

##Both Cophenetic and Silhouette values begin to drop suddenly when the number of outgoing patterns is 4.
nPatterns = 4
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)

# river plot

netAnalysis_river(cellchat, pattern = "incoming")

#> Please make sure you have load `library(ggalluvial)` when running this function
# dot plot

netAnalysis_dot(cellchat, pattern = "incoming")


saveRDS(cellchat, file = "cellchat.PU.rds")


############################################################################################
##Data input & processing and initialization of CellChat object
PU<- readRDS(file = "Npbmc.rds")
####reannotation of cell types(put the same cell type to a group ,if need)
new.cluster.ids <- c("Immune cells","Immune cells","Keratinocytes","Keratinocytes","Mitotic","Melanocytes")
names(new.cluster.ids) <- levels(PU)
PU <- RenameIdents(PU, new.cluster.ids)
table(Idents(PU))
levels(PU)
##########################Comparison analysis of multiple datasets using CellChat##########################################
#Load CellChat object of each dataset and then merge together
cellchat.AW<- readRDS(file = "cellchat.AW.rds")
cellchat.PU<- readRDS(file = "cellchat.PU.rds")
#unified the order of the cells type in two groups(if need)
group.new = levels(cellchat.AW@idents)
cellchat.PU <- liftCellChat(cellchat.PU, group.new)
#
object.list <- list(AW = cellchat.AW, PU = cellchat.PU)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
cellchat
cellchat@data.signaling

#Part I: Predict general principles of cell-cell communication
#Compare the total number of interactions and interaction strength

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

#Compare the number of interactions and interaction strength among different cell populations
#Differential number of interactions or interaction strength among different cell populations

par(mfrow = c(1,2), xpd=TRUE)
gg1 <-netVisual_diffInteraction(cellchat, weight.scale = T)
gg2 <-netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

#Do heatmap based on a merged object

gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2

#We then can show the number of interactions or interaction strength between any two cell types in each dataset.

weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

###Differential number of interactions or interaction strength among different cell types

weight.max <- getMaxWeight(object.list, slot.name = c("idents","net","net"), attribute = c("idents","count","count"))
weight.max
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= T, edge.weight.max = weight.max[3], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

#

par(mfrow = c(1,2), xpd=TRUE)
gg1 <-netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count", label.edge = T)
gg2 <-netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight", label.edge = T)

###################################################################################
#Compare the major sources and targets in 2D space
pdf(file="signalingRole_scatterAW-pu.pdf",width=6,height=5)
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)
dev.off()
###
#Furthermore, we can identify the specific signaling changes of Inflam.DC and cDC1 between NL and LS. 
## Identify signaling changes associated with one cell group
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Basal", signaling.exclude = "MIF")
#> Visualizing differential outgoing and incoming signaling changes from NL to LS
#> The following `from` values were not present in `x`: 0, -1, 1
#> The following `from` values were not present in `x`: 0, -1
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Mitotic", signaling.exclude = c("MIF"))
gg3 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Spinous", signaling.exclude = c("MIF"))
#> Visualizing differential outgoing and incoming signaling changes from NL to LS
#> The following `from` values were not present in `x`: 0, -1, 1
#> The following `from` values were not present in `x`: 0, -1
patchwork::wrap_plots(plots = list(gg1,gg2,gg3))

#

gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Immune cells", signaling.exclude = "MIF")
patchwork::wrap_plots(plots = list(gg1))


gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Melanocytes", signaling.exclude = "MIF")
patchwork::wrap_plots(plots = list(gg1))

#
saveRDS(cellchat, file = "cellchat.rds")

#########################
#Part II: Identify the conserved and context-specific signaling pathways
#Compare the overall information flow of each signaling pathway

gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2

#Compare outgoing (or incoming) signaling associated with each cell population
library(ComplexHeatmap)
i = 1
# combining all the identified signaling pathways from different datasets 

pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 6, height = 9)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 6, height = 9)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

###
pdf(file="ComplexHeatmap AW-pu incoming.pdf",width=8,height=7)
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 6, height = 9, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 6, height = 9, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

#

pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 6, height = 9, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 6, height = 9, color.heatmap = "OrRd")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

#Part III: Identify the upgulated and down-regulated signaling ligand-receptor pairs
netVisual_bubble(cellchat, sources.use = 5, targets.use = c(1:5),  comparison = c(1, 2), angle.x = 45)
#
gg1 <- netVisual_bubble(cellchat, sources.use = 5, targets.use = c(1:5),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in PU", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = 5, targets.use = c(1:5),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in PU", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2

##

gg1 <- netVisual_bubble(cellchat, sources.use = 2, targets.use = c(1:5),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in PU", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = 2, targets.use = c(1:5),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in PU", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2



# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "PU"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in PU
net.up <- subsetCommunication(cellchat, net = net, datasets = "PU",ligand.logFC = 0.2, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in AW
net.down <- subsetCommunication(cellchat, net = net, datasets = "AW",ligand.logFC = -0.1, receptor.logFC = -0.1)
#
gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)
#

pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = 5, targets.use = c(1:5), comparison = c(1, 2),  angle.x = 45, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = 5, targets.use = c(1:5), comparison = c(1, 2),  angle.x = 45, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
gg1 + gg2

#
# Chord diagram

par(mfrow = c(1,2), xpd=TRUE)
gg1 <-netVisual_chord_gene(object.list[[2]], sources.use = 5, targets.use = c(1:5), slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Note: The first link end is drawn out of sector 'MIF'.
gg2 <-netVisual_chord_gene(object.list[[1]], sources.use = 5, targets.use = c(1:5), slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))

#

par(mfrow = c(1,2), xpd=TRUE)
gg1 <-netVisual_chord_gene(object.list[[2]], sources.use = 2, targets.use = c(1:5), slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Note: The first link end is drawn out of sector 'MIF'.
gg2 <-netVisual_chord_gene(object.list[[1]], sources.use = 2, targets.use = c(1:5), slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))

#Part IV: Visually compare cell-cell communication using Hierarchy plot, Circle plot or Chord diagram
pathways.show <- c("EGF") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

#
pathways.show <- c("EGF") 
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
#> Do heatmap based on a single object 
#> 
#> Do heatmap based on a single object
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
dev.off()
#
pathways.show <- c("EGF") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
}
# Chord diagram
pathways.show <- c("EGF")
par(mfrow = c(1, 2), xpd=TRUE)
# compare all the interactions sending from Inflam.FIB to DC cells
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], sources.use = 2, targets.use = c(1:5), signaling = pathways.show, lab.cex = 0.5, title.name = paste0("EGF Signaling from Basal - ", names(object.list)[i]))
}
#
pathways.show <- c("GRN")
par(mfrow = c(1, 2), xpd=TRUE)
# compare all the interactions sending from Inflam.FIB to DC cells
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], sources.use = 5, targets.use = c(1:5), signaling = pathways.show, lab.cex = 0.5, title.name = paste0("GRN Signaling from Melanocytes - ", names(object.list)[i]))
}
# compare all the interactions sending from fibroblast to inflamatory immune cells
par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], sources.use = c(3), targets.use = c(1,2,4,5),  title.name = paste0("Signaling received by Mel and .Ke - ", names(object.list)[i]), legend.pos.x = 10)
}

#
par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], sources.use = 2, targets.use = c(1,3,4,5),slot.name = "netP", title.name = paste0("Signaling pathways sending from basal - ", names(object.list)[i]), legend.pos.x = 10)
}

#Part V: Compare the signaling gene expression distribution between different datasets
cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("AW", "PU")) # set factor level
plotGeneExpression(cellchat, signaling = "GRN", split.by = "datasets", colors.ggplot = T)
#Save the merged CellChat object
saveRDS(cellchat, file = "cellchat_comparisonAnalysis_AW_vs_PU.rds")

