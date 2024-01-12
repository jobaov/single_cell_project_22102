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
data <- readRDS( )


#Decide if data should be subsetted and then run standard seurat workflow
#We will probably subset for each condition for the effector T cells

AS_NA <- subset(data, subset = diseasegroup == "AS_NA")

#Standard workflow
AS_NA <- FindVariableFeatures(AS_NA, selection.method = "vst", nfeatures = 2000)
AS_NA <- ScaleData(AS_NA, features = VariableFeatures(object = AS_NA))
AS_NA <- RunPCA(AS_NA, features = VariableFeatures(object = AS_NA))
ElbowPlot(AS_NA)

#Change dimensions

#Run FindNeighbors and FindClusters
AS_NA <- FindNeighbors(AS_NA, dims = 1:10)
AS_NA<- FindClusters(AS_NA, resolution = c(0.1,0.2, 0.3))

#Change resolution (consider clustree tool to decide on an accurate resolution)

Idents(AS_NA) <- "RNA_snn_res.0.2"
AS_NA <- RunUMAP(AS_NA, dims = 1:10)

DimPlot(AS_NA, group.by = "RNA_snn_res.0.2", label = TRUE)

#############################################################
#Inspect cluster markers
FeaturePlot(AS_NA,
            reduction = "umap",
            features = c("CCR7", "SELL"),
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE) # 1

FeaturePlot(AS_NA,
            reduction = "umap",
            features = c("CREM", "CD69"),
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE)

######################################################
#Trajectory analysis

#Utilize Slingshot vignette to perform trajectory inference,
#selecting the starting cluster based on your cell annotation.

# Save the objects as separate matrices for input in slingshot
dimred <- AS_NA@reductions$umap@cell.embeddings

#Change resolution name

clustering <- AS_NA$RNA_snn_res.0.2
counts <- as.matrix(AS_NA@assays$RNA@counts[AS_NA@assays$RNA@var.features, ])

# Run default Slingshot lineage identification
set.seed(1)
lineages <- getLineages(data = dimred,
                        clusterLabels = clustering,
                        start.clus = "0") #define where to start the trajectories

lineages <- as.SlingshotDataSet(lineages)
lineages

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












