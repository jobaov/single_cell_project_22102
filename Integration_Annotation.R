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


#Load Seurat object
data <- readRDS("")

#Visualize Umap grouped by "batch effect"
DimPlot(data, reduction = "umap", group.by = "")
DimPlot(data, reduction = "umap", split.by = "")

## Sample integration

#Split the dataset into a list of Seurat objects (xx, xx).
#Normalize and identify variable features for each group independently.


#Add group to split by
obj.list <- SplitObject(data, split.by = "")

for(i in 1:length(obj.list)){
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
  obj.list[[i]]<- FindVariableFeatures(object = obj.list[[i]])
}


#Speed up integration, by running as parallel session
library(future)
plan("multisession", workers = 4)
options(future.globals.maxSize = 8000 * 1024^2)


features <- SelectIntegrationFeatures(object.list = obj.list)

anchors <- FindIntegrationAnchors(object.list = obj.list,
                                  anchor.features = features)

pbmc.integrated <- IntegrateData(anchorset = anchors)


DefaultAssay(pbmc.integrated) <- "integrated"
plan("sequential")

#Run the standard workflow for visualization and clustering (scaling, pca, find neighbors
#and clusters with different resolutions)

all.genes <- rownames(pbmc.integrated)
pbmc.integrated <- ScaleData(pbmc.integrated, features = all.genes)

pbmc.integrated <- RunPCA(pbmc.integrated, features = VariableFeatures(object = pbmc.integrated))

ElbowPlot(pbmc.integrated)

pbmc.integrated <- FindNeighbors(pbmc.integrated, dims = 1:13)
pbmc.integrated <- FindClusters(pbmc.integrated, resolution = c(0.3, 0.5, 0.7))

# Visualize the new clustering with a UMAP split by the batch corrected group.

pbmc.integrated <- RunUMAP(pbmc.integrated, dims = 1:30)

DimPlot(pbmc.integrated, reduction = "umap")

DimPlot(pbmc.integrated, group.by = "", label = TRUE)
DimPlot(pbmc.integrated, split.by = "", label = TRUE)


# UMAP of cells in each cluster by sample
DimPlot(pbmc.integrated,
        label = TRUE,
        split.by = "")  + NoLegend()

metrics <-  c("nCount_RNA", "nFeature_RNA", "percent.mt")

FeaturePlot(pbmc.integrated,
            reduction = "umap",
            features = metrics,
            pt.size = 0.4,
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE)

#Determine what resolution is the most optimal using clustree? and set Idents to the best resolution

library(clustree)
# integrated_snn_res.0.1 to 0.7 correspond to the different clustering resolutions

clustree(pbmc.integrated, prefix = "integrated_snn_res.")

# set idents to the best resolution
Idents(pbmc.integrated) <- "RNA_snn_res.0.3"

# 3. Cell Type assignment

#feature plots

# Identify conserved cell type markers to identify the cell types corresponding to the remaining clusters

#Identify conserved cell type markers
DefaultAssay(pbmc.integrated) <- "RNA"

cluster0_conserved_markers <- FindConservedMarkers(pbmc.integrated,
                                                   ident.1 = 0,
                                                   grouping.var = "diseasegroup",
                                                   only.pos = TRUE, min.pct = 0.25,  min.diff.pct = 0.25,
                                                   logfc.threshold = 0.25)


# Rename all identities
pbmc.integrated <- RenameIdents(object = pbmc.integrated,
                                "0" = "T-cells",
                                "1" = "T-cells",
                                "2" = "T-cells",
                                "3" = "T-cells",
                                "4" = "NK cells",
                                "5" = "Monocytes",
                                "6" = "T-cells",
                                "7" = "Monocytes",
                                "8" = "NK cells",
                                "9" = "Monocytes",
                                "10" = "B-cells",
                                "11" = "Monocytes",
                                "12" = "Platelet")

# Create a new column in metadata with cell type annotations from idents
pbmc.integrated@meta.data$Cell_ann <- Idents(pbmc.integrated)

# Now, you can view the updated metadata
View(pbmc.integrated@meta.data)


# visualize data
clusters <- DimPlot(pbmc.integrated, reduction = 'umap', label = TRUE)
xx <- DimPlot(pbmc.integrated, reduction = 'umap', group.by = '')
celltype <- DimPlot(pbmc.integrated, reduction = 'umap', group.by = 'Cell_ann')

xx|clusters
xx|celltype











