#Load libraries

library(Seurat)
library(SingleCellExperiment)
library(tradeSeq)
library(dplyr)
suppressPackageStartupMessages({
  library(slingshot)
})
library(cowplot)

#Load Seurat object
AS_NA <- readRDS("/home/projects/Group1/Teff_annotation.rds" )


#Decide if data should be subsetted and then run standard seurat workflow
#We will probably subset for each condition for the effector T cells

AS_NA <- subset(AS_NA, subset = diseasegroup == "AS_NA")
DefaultAssay(AS_NA) <- "integrated"

#Standard workflow
AS_NA <- FindVariableFeatures(AS_NA, selection.method = "vst", nfeatures = 2000)
AS_NA <- ScaleData(AS_NA, features = VariableFeatures(object = AS_NA))
AS_NA <- RunPCA(AS_NA, features = VariableFeatures(object = AS_NA))
ElbowPlot(AS_NA)

#Change dimensions

#Run FindNeighbors and FindClusters
AS_NA <- FindNeighbors(AS_NA, dims = 1:30)
AS_NA<- FindClusters(AS_NA, resolution = c(0.1,0.2, 0.3, 0.5, 0.7))

DimPlot(AS_NA, group.by = "integrated_snn_res.0.1", label = TRUE, reduction = "umap")

DimPlot(AS_NA, group.by = "integrated_snn_res.0.2", label = TRUE, reduction = "umap")

DimPlot(AS_NA, group.by = "integrated_snn_res.0.3", label = TRUE, reduction = "umap")


#Change resolution (consider clustree tool to decide on an accurate resolution)
library(clustree)
# integrated_snn_res.0.1 to 0.7 correspond to the different clustering resolutions

clustering_tree_res <- clustree(AS_NA, prefix = "integrated_snn_res.")
clustering_tree_res

Idents(AS_NA) <- "integrated_snn_res.0.3"
AS_NA <- RunUMAP(AS_NA, dims = 1:30)

DimPlot(AS_NA, group.by = "integrated_snn_res.0.3", label = TRUE, reduction = "umap")

#############################################################
#Inspect cluster markers
FeaturePlot(AS_NA,
            reduction = "umap",
            features = c("IFNG", "XCL1", "IL2"),
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE) #cluster 3 = TH1 and chemokines

FeaturePlot(AS_NA,
            reduction = "umap",
            features = c("IL5", "IL13"),
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE) # cluster 4 = type 2 cytokine genes TH2 cells

FeaturePlot(AS_NA,
            reduction = "umap",
            features = c("IL17F", "IL22", "IL17A"),
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE) #cluster 0 and cluster 2 = TH17

FeaturePlot(AS_NA,
            reduction = "umap",
            features = c("IFI6", "MX1", "ISG20", "OAS1", "IFIT1", "IFI44L"),
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE) #TH-IFNR subset

FeaturePlot(AS_NA,
            reduction = "umap",
            features = c("IL7R", "S100A4"),
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE) #CD4 memory cells

#cell annotation on subset

#Cluster 0, 1 and 2,3 are biologically uncharacterized activated T-cells
AS_NA <- RenameIdents(object = AS_NA,
                                "0" = "TH-ACT1",
                                "1" = "TH-ACT2_1",
                                "2" = "TH-ACT2_2",
                                "3" = "TH-ACT3",
                                "4" = "TH-1",
                                "5" = "TH-2",
                                "6" = "TH-IFNR"
)

# Create a new column in metadata with cell type annotations from idents
AS_NA@meta.data$Cell_ann2 <- Idents(AS_NA)
######################################################
#Trajectory analysis

#Utilize Slingshot vignette to perform trajectory inference,
#selecting the starting cluster based on your cell annotation.

#DefaultAssay(AS_NA) <- "RNA"
# Save the objects as separate matrices for input in slingshot
dimred <- AS_NA@reductions$umap@cell.embeddings

#Change resolution name

clustering <- AS_NA$integrated_snn_res.0.3
#counts <- as.matrix(AS_NA@assays$RNA@counts[AS_NA@assays$RNA@var.features, ])

# Run default Slingshot lineage identification
set.seed(1)
lineages <- getLineages(data = dimred,
                        clusterLabels = clustering,
                        start.clus = "1") #define where to start the trajectories

lineages <- as.SlingshotDataSet(lineages)
#Look at lineages
lineages


# Define a color pallete to use
pal <- c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8, "Set2"))

#Visualize as curves
curves <- getCurves(lineages, approx_points = 300, thresh = 0.01, stretch = 0.8, allow.breaks = FALSE, shrink = 0.99)
curves <- as.SlingshotDataSet(curves)
curves

plot(dimred, col = pal[clustering], asp = 1, pch = 16)
lines(curves, lwd = 3, col = "black")

NA_annotated <- DimPlot(AS_NA, group.by = "Cell_ann2", label = TRUE, reduction = "umap")
ggsave("NA_annotate.png", NA_annotated, path= "Plots/Trajectory")

saveRDS(AS_NA, "/home/projects/Group1/Teff_trajectory_AS_NA.rds")

#########################################################################
#subset for the Allergic
#Load Seurat object
AS_AL <- readRDS("/home/projects/Group1/Teff_annotation.rds" )


#Decide if data should be subsetted and then run standard seurat workflow
#We will probably subset for each condition for the effector T cells

AS_AL <- subset(AS_AL, subset = diseasegroup == "AS_AL")
DefaultAssay(AS_AL) <- "integrated"

#Standard workflow
AS_AL <- FindVariableFeatures(AS_AL, selection.method = "vst", nfeatures = 2000)
AS_AL <- ScaleData(AS_AL, features = VariableFeatures(object = AS_AL))
AS_AL <- RunPCA(AS_AL, features = VariableFeatures(object = AS_AL))
ElbowPlot(AS_AL)

#Change dimensions

#Run FindNeighbors and FindClusters
AS_AL <- FindNeighbors(AS_AL, dims = 1:30)
AS_AL<- FindClusters(AS_AL, resolution = c(0.1,0.2, 0.3, 0.5, 0.7))

DimPlot(AS_AL, group.by = "integrated_snn_res.0.1", label = TRUE, reduction = "umap")

DimPlot(AS_AL, group.by = "integrated_snn_res.0.2", label = TRUE, reduction = "umap")

DimPlot(AS_AL, group.by = "integrated_snn_res.0.3", label = TRUE, reduction = "umap")


#Change resolution (consider clustree tool to decide on an accurate resolution)
library(clustree)
# integrated_snn_res.0.1 to 0.7 correspond to the different clustering resolutions

clustering_tree_res <- clustree(AS_AL, prefix = "integrated_snn_res.")
clustering_tree_res

Idents(AS_AL) <- "integrated_snn_res.0.3"
AS_AL <- RunUMAP(AS_AL, dims = 1:30)

DimPlot(AS_AL, group.by = "integrated_snn_res.0.3", label = TRUE, reduction = "umap")

#Cluster 0, 1 and 2 are biologically uncharacterized activated T-cells
AS_AL <- RenameIdents(object = AS_AL,
                                "0" = "TH-ACT1",
                                "1" = "TH-ACT2",
                                "2" = "TH-ACT3",
                                "3" = "TH-1",
                                "4" = "TH-2",
                                "5" = "TH-IFNR",
                                "6" = "TH-17"
)

# Create a new column in metadata with cell type annotations from idents
AS_AL@meta.data$Cell_ann_AL <- Idents(AS_AL)


#Trajectory analysis

# Save the objects as separate matrices for input in slingshot
dimred_AL <- AS_AL@reductions$umap@cell.embeddings

#Change resolution name

clustering_AL <- AS_AL$integrated_snn_res.0.3
Cells_AL <- AS_AL$Cell_ann_AL
#counts <- as.matrix(AS_NA@assays$RNA@counts[AS_NA@assays$RNA@var.features, ])

# Run default Slingshot lineage identification
set.seed(1)
lineages_AL <- getLineages(data = dimred_AL,
                        clusterLabels = clustering_AL,
                        start.clus = "0") #define where to start the trajectories

lineages_AL <- as.SlingshotDataSet(lineages_AL)
#Look at lineages
lineages_AL

#Visualize as curves
curves_AL <- getCurves(lineages_AL, approx_points = 300, thresh = 0.01, stretch = 0.8, allow.breaks = FALSE, shrink = 0.99)
curves_AL <- as.SlingshotDataSet(curves_AL)
curves_AL


#AL_trajectory <- recordPlot()
AL_trajectory <- plot(dimred_AL, col = pal[Cells_AL], asp = 1, pch = 16) +
lines(curves_AL, lwd = 3, col = "black")

AL_annotated <- DimPlot(AS_AL, group.by = "Cell_ann_AL", label = TRUE, reduction = "umap")
ggsave("AL_annotate.png", AL_annotated, path= "Plots/Trajectory")

saveRDS(AS_AL, "/home/projects/Group1/Teff_trajectory_AS_AL.rds")

