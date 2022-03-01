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

