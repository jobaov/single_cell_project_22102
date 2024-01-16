library(Seurat)
library(data.table)
library(ggplot2)
library(BiocSingular)
library(scDblFinder)
library(SingleCellExperiment)
library(SingleR)
library(scater)
library(SeuratData)
library(patchwork)
library(dplyr)
library(tidyr)
library(multtest)
library(metap)
library(tibble)
library(purrr)

#Use the integrated and annotated data instead

Teff <- readRDS("/home/projects/Group1/Teff_annotation.rds")

Treg <- readRDS("/home/projects/Group1/treg_annotated.rds")

#Set default Assay, since this is the data we want to work with
DefaultAssay(Teff) <- "integrated"
DefaultAssay(Treg) <- "integrated"

obj.list <- list(Teffector = Teff,
                 Tregulatory = Treg)

#This code was related to batch correction, so don't think we need this, since it is already done
#for(i in 1:length(obj.list)){
#  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
#  obj.list[[i]]<- FindVariableFeatures(object = obj.list[[i]])
#}


#Speed up integration, by running as parallel session
library(future)
plan("multisession", workers = 4)
options(future.globals.maxSize = 8000 * 1024^2)


features <- SelectIntegrationFeatures(object.list = obj.list)

anchors <- FindIntegrationAnchors(object.list = obj.list,
                                  anchor.features = features)

pbmc.integrated <- IntegrateData(anchorset = anchors)


DefaultAssay(pbmc.integrated) <- "integrated"
saveRDS(pbmc.integrated, "/home/projects/Group1/Teff_Treg_integration_annotated.rds")

#################################################################################
plan("sequential")

#Run the standard workflow for visualization and clustering (scaling, pca, find neighbors
#and clusters with different resolutions)

all.genes <- rownames(pbmc.integrated)
pbmc.integrated <- ScaleData(pbmc.integrated, features = all.genes)

pbmc.integrated <- RunPCA(pbmc.integrated, features = VariableFeatures(object = pbmc.integrated))

ElbowPlot(pbmc.integrated)

pbmc.integrated <- FindNeighbors(pbmc.integrated, dims = 1:30)
pbmc.integrated <- FindClusters(pbmc.integrated, resolution = c(0.1, 0.2, 0.3, 0.5, 0.7))

########################################################################################

# Visualize the new clustering with a UMAP split by the batch corrected group.

pbmc.integrated <- RunUMAP(pbmc.integrated, dims = 1:30)

DimPlot(pbmc.integrated, reduction = "umap")

integrated_celltype <- DimPlot(pbmc.integrated, group.by = "celltype", label = TRUE)
integrated_celltype
ggsave("integrated_celltype.png", integrated_celltype, path= "Plots/Teff_Treg_integration")

integrated_cellsubtype <- DimPlot(pbmc.integrated_joinlayer, group.by = "Cell_ann", label = TRUE)
integrated_cellsubtype
ggsave("integrated_cellsubtype.png", integrated_cellsubtype, path= "Plots/Teff_Treg_integration")


saveRDS(pbmc.integrated, "/home/projects/Group1/Teff_Treg_integration_cluster_annotated.rds")


#######################################################################################
#make two cellchat objects one for each condition (NA vs AL)
#Then merge the cellchat objects and follow the vignette for combined CCC

library(CellChat)
library(patchwork)
library(future)
options(stringsAsFactors = FALSE)

#Trying to avoid error
pbmc.integrated_joinlayer <- JoinLayers(pbmc.integrated)

NA_integrated <- subset(pbmc.integrated_joinlayer, subset = diseasegroup == "AS_NA")

AL_integrated <- subset(pbmc.integrated_joinlayer, subset = diseasegroup == "AS_AL")

#create Cell chat object will not accept integrated assay
DefaultAssay(NA_integrated) <- "RNA"
DefaultAssay(AL_integrated) <- "RNA"

#Create a CellChat object
cellchat_NA <- createCellChat(object = NA_integrated,
                           meta = NA_integrated@meta.data,
                           group.by = "Cell_ann")

cellchat_AL <- createCellChat(object = AL_integrated,
                              meta = AL_integrated@meta.data,
                              group.by = "Cell_ann")

#Set interaction-ligand database
CellChatDB <- CellChatDB.human

cellchat_AL@DB <- CellChatDB

cellchat_NA@DB <- CellChatDB

##########################################################################
#ALLERGIC

# subset the expression data of signaling genes for saving computation cost
cellchat_AL <- subsetData(cellchat_AL) # This step is necessary even if using the whole database
plan("multisession", workers = 4) # do parallel

cellchat_AL <- identifyOverExpressedGenes(cellchat_AL)
cellchat_AL <- identifyOverExpressedInteractions(cellchat_AL)

#Part II: Inference of cell-cell communication network
cellchat_AL <- computeCommunProb(cellchat_AL)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_AL <- filterCommunication(cellchat_AL, min.cells = 10)


#Extract the inferred cellular communication network as a data frame
df.net_AL <- subsetCommunication(cellchat_AL)
df.net_AL <- subsetCommunication(cellchat_AL) #returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. Set slot.name = "netP" to access the the inferred communications at the level of signaling pathways


#df.net <- subsetCommunication(cellchat, signaling = c("TNF", "TGFb")) #gives the inferred cell-cell communications mediated by signaling WNT and TGFb.

#Infer the cell-cell communication at a signaling pathway level

cellchat_AL <- computeCommunProbPathway(cellchat_AL)

#Calculate the aggregated cell-cell communication network
cellchat_AL <- aggregateNet(cellchat_AL)
groupSize_AL <- as.numeric(table(cellchat_AL@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat_AL@net$count, vertex.weight = groupSize_AL, weight.scale = T, label.edge= F, title.name = "Number of interactions Allergic")
netVisual_circle(cellchat_AL@net$weight, vertex.weight = groupSize_AL, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength Allergic")

#########################################################

pathways_AL <- cellchat_AL@netP$pathways

#Part III: Visualization of cell-cell communication network

# Circle plot
par(mfrow=c(1,2))
netVisual_aggregate(cellchat_AL, signaling = pathways_AL[1], layout = "circle")
netVisual_aggregate(cellchat_AL, signaling = pathways_AL[2], layout = "circle")
netVisual_aggregate(cellchat_AL, signaling = pathways_AL[3], layout = "circle")
netVisual_aggregate(cellchat_AL, signaling = pathways_AL[4], layout = "circle")
netVisual_aggregate(cellchat_AL, signaling = pathways_AL[5], layout = "circle")
netVisual_aggregate(cellchat_AL, signaling = pathways_AL[6], layout = "circle")

netVisual_heatmap(cellchat_AL, signaling = pathways_AL[1], color.heatmap = "Reds")

saveRDS(cellchat_AL, "/home/projects/Group1/CellChat_AL.rds")

################################################################################
#NON-ALLERGIC

# subset the expression data of signaling genes for saving computation cost
cellchat_NA <- subsetData(cellchat_NA) # This step is necessary even if using the whole database

cellchat_NA <- identifyOverExpressedGenes(cellchat_NA)
cellchat_NA <- identifyOverExpressedInteractions(cellchat_NA)

#Part II: Inference of cell-cell communication network
cellchat_NA <- computeCommunProb(cellchat_NA)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_NA <- filterCommunication(cellchat_NA, min.cells = 10)


#Extract the inferred cellular communication network as a data frame
df.net_NA <- subsetCommunication(cellchat_NA)
df.net_NA <- subsetCommunication(cellchat_NA) #returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. Set slot.name = "netP" to access the the inferred communications at the level of signaling pathways


#df.net <- subsetCommunication(cellchat, signaling = c("TNF", "TGFb")) #gives the inferred cell-cell communications mediated by signaling WNT and TGFb.

#Infer the cell-cell communication at a signaling pathway level

cellchat_NA <- computeCommunProbPathway(cellchat_AL)

#Calculate the aggregated cell-cell communication network
cellchat_NA <- aggregateNet(cellchat_NA)
groupSize_NA <- as.numeric(table(cellchat_NA@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat_NA@net$count, vertex.weight = groupSize_NA, weight.scale = T, label.edge= F, title.name = "Number of interactions Non-allergic")
netVisual_circle(cellchat_NA@net$weight, vertex.weight = groupSize_NA, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength Non-allergic")

#########################################################

pathways_NA <- cellchat_NA@netP$pathways

#Part III: Visualization of cell-cell communication network

# Circle plot
par(mfrow=c(1,2))
netVisual_aggregate(cellchat_NA, signaling = pathways_NA[1], layout = "circle")
netVisual_aggregate(cellchat_NA, signaling = pathways_NA[2], layout = "circle")
netVisual_aggregate(cellchat_NA, signaling = pathways_NA[3], layout = "circle")
netVisual_aggregate(cellchat_NA, signaling = pathways_NA[4], layout = "circle")
netVisual_aggregate(cellchat_NA, signaling = pathways_NA[5], layout = "circle")
netVisual_aggregate(cellchat_NA, signaling = pathways_NA[6], layout = "circle")

par(mfrow=c(1,1))
netVisual_heatmap(cellchat_NA, signaling = pathways_NA[1], color.heatmap = "Reds")

saveRDS(cellchat_NA, "/home/projects/Group1/CellChat_NA.rds")

####################################################################

#THIS WAS CODE FOR THE COMPARISON BETWEEN ALLERGIC AND NON-ALLERGIC, DOES NOT RUN NOW

#object.list <- list(AS_NA = cellchat_NA, AS_AL = cellchat_AL)
#cellchat <- mergeCellChat(object.list, add.names = names(object.list))
#cellchat@DB <- CellChatDB

#Compare the total number of interactions and interaction strength

#CellChat compares the the total number of interactions and interaction strength of the inferred
#cell-cell communication networks from different biological conditions.

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2


#Compare the number of interactions and interaction strength among different cell populations
#The differential number of interactions or interaction strength in the cell-cell communication network between two datasets can be visualized using circle plot, where red
#(or blue ) colored edges represent increased (or decreased) signaling in the second dataset compared to the first one


par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")


#We can also show differential number of interactions or interaction strength in a greater details
#using a heatmap.

gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2










