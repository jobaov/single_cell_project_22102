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
AS_NA <- readRDS("Data/teff_asthm_filtered.rds" )


#Decide if data should be subsetted and then run standard seurat workflow
#We will probably subset for each condition for the effector T cells

AS_NA <- subset(data, subset = diseasegroup == "AS_NA")

#Standard workflow
AS_NA <- NormalizeData(AS_NA, normalization.method = "LogNormalize", scale.factor = 10000)
AS_NA <- FindVariableFeatures(AS_NA, selection.method = "vst", nfeatures = 2000)
AS_NA <- ScaleData(AS_NA, features = VariableFeatures(object = AS_NA))
AS_NA <- RunPCA(AS_NA, features = VariableFeatures(object = AS_NA))
ElbowPlot(AS_NA)

#Change dimensions

#Run FindNeighbors and FindClusters
AS_NA <- FindNeighbors(AS_NA, dims = 1:30)
AS_NA<- FindClusters(AS_NA, resolution = c(0.1,0.2, 0.3))

DimPlot(AS_NA, group.by = "RNA_snn_res.0.1", label = TRUE, reduction = "umap")

DimPlot(AS_NA, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "umap")

DimPlot(AS_NA, group.by = "RNA_snn_res.0.3", label = TRUE, reduction = "umap")


#Change resolution (consider clustree tool to decide on an accurate resolution)

Idents(AS_NA) <- "RNA_snn_res.0.3"
AS_NA <- RunUMAP(AS_NA, dims = 1:30)

DimPlot(AS_NA, group.by = "RNA_snn_res.0.3", label = TRUE)

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
            features = c("ARL6IP5", "ARPC5"),
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE) #Endocytosis and membrane trafficking

FeaturePlot(AS_NA,
            reduction = "umap",
            features = c("IL7R", "S100A4"),
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE) #CD4 memory cells

FeaturePlot(AS_NA,
            reduction = "umap",
            features = c("CD8A"),
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE) #CD8+ T cells

FeaturePlot(AS_NA,
            reduction = "umap",
            features = c("IL9", "IRF4"),
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE) #CD8+ T cells

######################################################
#Trajectory analysis

#Utilize Slingshot vignette to perform trajectory inference,
#selecting the starting cluster based on your cell annotation.

DefaultAssay(AS_NA) <- "RNA"
# Save the objects as separate matrices for input in slingshot
dimred <- AS_NA@reductions$umap@cell.embeddings

#Change resolution name

clustering <- AS_NA$RNA_snn_res.0.3
#counts <- as.matrix(AS_NA@assays$RNA@counts[AS_NA@assays$RNA@var.features, ])

# Run default Slingshot lineage identification
set.seed(1)
lineages <- getLineages(data = dimred,
                        clusterLabels = clustering,
                        start.clus = "0") #define where to start the trajectories

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












