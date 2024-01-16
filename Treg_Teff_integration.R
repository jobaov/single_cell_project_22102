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
ggsave("integrated_celltype.png", integrated_celltype, path= "Plots/Teff_Treg_integration_ann")

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


object.list <- list(AS_NA = cellchat_NA, AS_AL = cellchat_AL)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

#Set interaction-ligand database
CellChatDB <- CellChatDB.human


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










