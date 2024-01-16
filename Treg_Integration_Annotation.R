# Read libraries

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
library(clustree)

# Load Seurat object
data <- readRDS("Data/treg_asthm_filtered.rds")

# Clustering

# Run FindNeighbors and FindClusters
data <- FindNeighbors(data, dims = 1:30)
data <- FindClusters(data, resolution = c(0.1, 0.2, 0.3, 0.5, 0.7))
ElbowPlot(data)

data <- RunUMAP(data, dims = 1:30)

DimPlot(data, group.by = "RNA_snn_res.0.1", label = TRUE, reduction = "umap")

DimPlot(data, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "umap")

DimPlot(data, group.by = "RNA_snn_res.0.3", label = TRUE, reduction = "umap")

DimPlot(data, group.by = "RNA_snn_res.0.5", label = TRUE, reduction = "umap")

DimPlot(data, group.by = "RNA_snn_res.0.7", label = TRUE, reduction = "umap")

####################################################################################

# Sample integration

# Visualize Umap grouped by "batch effect"
DimPlot(data, reduction = "umap", group.by = "donor")

DimPlot(data, reduction = "umap", split.by = "donor")

# Split the dataset into a list of Seurat objects (xx, xx).
# Normalize and identify variable features for each group independently.

# Add group to split by
obj.list <- SplitObject(data, split.by = "donor")

for(i in 1:length(obj.list)){
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
  obj.list[[i]]<- FindVariableFeatures(object = obj.list[[i]])
}

# Speed up integration, by running as parallel session
library(future)
plan("multisession", workers = 4)
options(future.globals.maxSize = 8000 * 1024^2)

features <- SelectIntegrationFeatures(object.list = obj.list)

anchors <- FindIntegrationAnchors(object.list = obj.list,
                                  anchor.features = features)

data.integrated <- IntegrateData(anchorset = anchors)

DefaultAssay(data.integrated) <- "integrated"
plan("sequential")

saveRDS(data.integrated, "Data/treg_integrated.rds")

# Run the standard workflow for visualization and clustering (scaling, pca, find neighbors
# and clusters with different resolutions)

all.genes <- rownames(data.integrated)
data.integrated <- ScaleData(data.integrated,
                             features = all.genes)

data.integrated <- RunPCA(data.integrated,
                          features = VariableFeatures(object = data.integrated))

ElbowPlot(data.integrated)

data.integrated <- FindNeighbors(data.integrated,
                                 dims = 1:30)
data.integrated <- FindClusters(data.integrated,
                                resolution = c(0.1, 0.2, 0.3, 0.5, 0.7))

########################################################################################

# Visualize the new clustering with a UMAP split by the batch corrected group.

data.integrated <- RunUMAP(data.integrated,
                           dims = 1:30)

DimPlot(data.integrated, reduction = "umap")

DimPlot(data.integrated, group.by = "donor", label = TRUE)
DimPlot(data.integrated, split.by = "donor", label = TRUE)

# UMAP of cells in each cluster by sample
DimPlot(data.integrated,
        label = TRUE,
        split.by = "donor")  + NoLegend()

metrics <-  c("nCount_RNA", "nFeature_RNA", "percent.mt")

FeaturePlot(data.integrated,
            reduction = "umap",
            features = metrics,
            pt.size = 0.4,
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE)

# Determine what resolution is the most optimal using clustree? and set Idents to the best resolution

library(clustree)
# integrated_snn_res.0.1 to 0.7 correspond to the different clustering resolutions

clustree(data.integrated, prefix = "integrated_snn_res.")

# set idents to the best resolution
Idents(data.integrated) <- "integrated_snn_res.0.3"

DimPlot(data.integrated, group.by = "integrated_snn_res.0.3", label = TRUE)

####################################################################################
# Cell Type assignment

#feature plots

#Inspect cluster markers

FeaturePlot(data.integrated,
            reduction = "umap",
            features = c("IFI6", "MX1", "ISG20", "OAS1", "IFIT1", "IFI44L"),
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE) # TH-IFNR subset

FeaturePlot(,
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

FeaturePlot(data.integrated,
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

# Identify conserved cell type markers to identify the cell types corresponding to the remaining clusters

# Identify conserved cell type markers
DefaultAssay(data.integrated) <- "RNA"

data.integrated <- JoinLayers(data.integrated)

cluster0_conserved_markers <- FindConservedMarkers(data.integrated,
                                                   ident.1 = 0,
                                                   grouping.var = "donor",
                                                   only.pos = TRUE,
                                                   min.pct = 0.25,
                                                   min.diff.pct = 0.25,
                                                   logfc.threshold = 0.25)

cluster1_conserved_markers <- FindConservedMarkers(data.integrated,
                                                   ident.1 = 1,
                                                   grouping.var = "donor",
                                                   only.pos = TRUE,
                                                   min.pct = 0.25,
                                                   min.diff.pct = 0.25,
                                                   logfc.threshold = 0.25)

cluster2_conserved_markers <- FindConservedMarkers(data.integrated,
                                                   ident.1 = 2,
                                                   grouping.var = "donor",
                                                   only.pos = TRUE,
                                                   min.pct = 0.25,
                                                   min.diff.pct = 0.25,
                                                   logfc.threshold = 0.25)

annotations <- read.csv("/home/projects/22102_single_cell/day3/annotation.csv")

# Combine markers with gene descriptions
cluster0_ann_markers <- cluster0_conserved_markers %>%
  rownames_to_column(var="gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))

# Combine markers with gene descriptions
cluster1_ann_markers <- cluster1_conserved_markers %>%
  rownames_to_column(var="gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))

# Combine markers with gene descriptions
cluster2_ann_markers <- cluster2_conserved_markers %>%
  rownames_to_column(var="gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))

# Rename all identities
data.integrated <- RenameIdents(object = data.integrated,
                                "2" = "TregIFNR",
                                "0" = "TregACT1",
                                "1" = "TregACT2")

# Create a new column in metadata with cell type annotations from idents
data.integrated@meta.data$Cell_ann <- Idents(data.integrated)

# Now, you can view the updated metadata
View(data.integrated@meta.data)

# Visualize data
clusters <- DimPlot(data.integrated, reduction = 'umap', label = TRUE)
xx <- DimPlot(data.integrated, reduction = 'umap', split.by = 'donor')
celltype <- DimPlot(data.integrated, reduction = 'umap', group.by = 'Cell_ann')

xx|clusters
xx|celltype

saveRDS(data.integrated, "Data/treg_annotated.rds")

a <- table(data.integrated@meta.data$Cell_ann, data.integrated@meta.data$diseasegroup)
a[,1]/sum(a[,1])
a[,2]/sum(a[,2])

