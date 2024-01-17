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
data <- readRDS("Data/teff_asthm_filtered.rds")

#Clustering

#Run FindNeighbors and FindClusters
data <- FindNeighbors(data, dims = 1:30)
data <- FindClusters(data, resolution = c(0.1, 0.2, 0.3, 0.5, 0.7))
ElbowPlot(data)

data <- RunUMAP(data, dims = 1:30)

DimPlot(data, group.by = "RNA_snn_res.0.1", label = TRUE, reduction = "umap")

Teff_nonintegrated_res0.2 <- DimPlot(data, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "umap")
ggsave("Teff_nonintegrated_res0.2.png", Teff_nonintegrated_res0.2, path= "Plots/Teff_integration")

DimPlot(data, group.by = "RNA_snn_res.0.3", label = TRUE, reduction = "umap")

DimPlot(data, group.by = "RNA_snn_res.0.5", label = TRUE, reduction = "umap")

DimPlot(data, group.by = "RNA_snn_res.0.7", label = TRUE, reduction = "umap")


####################################################################################

# Sample integration


#Visualize Umap grouped by "batch effect"
Teff_nonintegrated_donor <- DimPlot(data, reduction = "umap", group.by = "donor")
ggsave("Teff_nonintegrated_donor.png", Teff_nonintegrated_donor, path= "Plots/Teff_integration")

Teff_nonintegrated_donor_clus_res02 <- DimPlot(data, reduction = "umap", split.by = "donor", group.by = "RNA_snn_res.0.2")
ggsave("Teff_nonintegrated_donor_clus_res02.png", Teff_nonintegrated_donor_clus_res02, path= "Plots/Teff_integration")

#Split the dataset into a list of Seurat objects.
#Normalize and identify variable features for each group independently.


#Add group to split by
obj.list <- SplitObject(data, split.by = "donor")

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
saveRDS(pbmc.integrated, "Data/Teff_integration.rds")

#read integrated data
pbmc.integrated <- readRDS("/home/projects/Group1/Teff_integration.rds")

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

umap_integrated_donor <- DimPlot(pbmc.integrated, group.by = "donor", label = TRUE)
umap_integrated_donor
ggsave("umap_integrated_donor.png", umap_integrated_donor, path= "Plots/Teff_integration")

umap_integrated_disease <- DimPlot(pbmc.integrated, split.by = "diseasegroup", label = TRUE, group.by = "integrated_snn_res.0.2")
umap_integrated_disease
ggsave("umap_integrated_disease.png", umap_integrated_disease, path= "Plots/Teff_integration")


# UMAP of cells in each cluster by sample
DimPlot(pbmc.integrated,
        label = TRUE,
        split.by = "diseasegroup")  + NoLegend()

metrics <-  c("nCount_RNA", "nFeature_RNA", "percent_mt")

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

clustering_tree_res <- clustree(pbmc.integrated, prefix = "integrated_snn_res.")
clustering_tree_res
ggsave("clustree_res.png", clustering_tree_res, path= "Plots/Teff_integration")

# set idents to the best resolution
Idents(pbmc.integrated) <- "integrated_snn_res.0.2"

umap_res_0.2 <- DimPlot(pbmc.integrated, group.by = "integrated_snn_res.0.2", label = TRUE, reduction = "umap")
umap_res_0.2
ggsave("umap_res_0.2.png", umap_res_0.2, path= "Plots/Teff_integration")

####################################################################################
# Cell Type assignment

#feature plots
#Inspect cluster markers
Feature_TH1 <- FeaturePlot(pbmc.integrated,
            reduction = "umap",
            features = c("IFNG", "XCL1", "IL2"),
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE) #cluster 3 = TH1 and chemokines
Feature_TH1
ggsave("Feature_TH1.png", Feature_TH1, path= "Plots/Teff_integration")

Feature_TH2 <- FeaturePlot(pbmc.integrated,
            reduction = "umap",
            features = c("IL5", "IL13"),
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE) # cluster 4 = type 2 cytokine genes TH2 cells
Feature_TH2
ggsave("Feature_TH2.png", Feature_TH2, path= "Plots/Teff_integration")

Feature_TH17 <- FeaturePlot(pbmc.integrated,
            reduction = "umap",
            features = c("IL17F", "IL22", "IL17A"),
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE) #cluster 6 = TH17
Feature_TH17
ggsave("Feature_TH17.png", Feature_TH17, path= "Plots/Teff_integration")

Feature_TH_IFNR <- FeaturePlot(pbmc.integrated,
            reduction = "umap",
            features = c("IFI6", "MX1", "ISG20", "OAS1", "IFIT1", "IFI44L"),
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE) #TH-IFNR subset cluster 5
Feature_TH_IFNR
ggsave("Feature_TH_IFNR.png", Feature_TH_IFNR, path= "Plots/Teff_integration")


Feature_memory <- FeaturePlot(pbmc.integrated,
            reduction = "umap",
            features = c("IL7R", "S100A4"),
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE) #CD4 memory cells
Feature_memory
ggsave("Feature_memory.png", Feature_memory, path= "Plots/Teff_integration")

#checking more markers from highly variable features pre-processing plot
FeaturePlot(pbmc.integrated,
                              reduction = "umap",
                              features = c("IL9","IL31"),
                              order = TRUE,
                              min.cutoff = 'q10',
                              label = TRUE)

FeaturePlot(pbmc.integrated,
            reduction = "umap",
            features = c("IL3","IL13", "CSF2" , "IL4", "IL5"),
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE) # cytokine cluster TH-2

FeaturePlot(pbmc.integrated,
            reduction = "umap",
            features = c("CCR6", "CCR4", "CCR10"),
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE) #TH-22 expression (TH-17)



# Identify conserved cell type markers to identify the cell types corresponding to the remaining clusters

#Identify conserved cell type markers
DefaultAssay(pbmc.integrated) <- "RNA"

seurat_joinlayer <- JoinLayers(pbmc.integrated)

cluster6_conserved_markers <- FindConservedMarkers(seurat_joinlayer,
                                                   ident.1 = 6,
                                                   grouping.var = "diseasegroup",
                                                   only.pos = TRUE, min.pct = 0.25,  min.diff.pct = 0.25,
                                                   logfc.threshold = 0.25)

cluster5_conserved_markers <- FindConservedMarkers(seurat_joinlayer,
                                                   ident.1 = 5,
                                                   grouping.var = "diseasegroup",
                                                   only.pos = TRUE, min.pct = 0.25,  min.diff.pct = 0.25,
                                                   logfc.threshold = 0.25)

cluster4_conserved_markers <- FindConservedMarkers(seurat_joinlayer,
                                                   ident.1 = 4,
                                                   grouping.var = "diseasegroup",
                                                   only.pos = TRUE, min.pct = 0.25,  min.diff.pct = 0.25,
                                                   logfc.threshold = 0.25)

cluster3_conserved_markers <- FindConservedMarkers(seurat_joinlayer,
                                                   ident.1 = 3,
                                                   grouping.var = "diseasegroup",
                                                   only.pos = TRUE, min.pct = 0.25,  min.diff.pct = 0.25,
                                                   logfc.threshold = 0.25)

cluster2_conserved_markers <- FindConservedMarkers(seurat_joinlayer,
                                                   ident.1 = 2,
                                                   grouping.var = "diseasegroup",
                                                   only.pos = TRUE, min.pct = 0.25,  min.diff.pct = 0.25,
                                                   logfc.threshold = 0.25)

cluster1_conserved_markers <- FindConservedMarkers(seurat_joinlayer,
                                                   ident.1 = 1,
                                                   grouping.var = "diseasegroup",
                                                   only.pos = TRUE, min.pct = 0.25,  min.diff.pct = 0.25,
                                                   logfc.threshold = 0.25)

cluster0_conserved_markers <- FindConservedMarkers(seurat_joinlayer,
                                                   ident.1 = 0,
                                                   grouping.var = "diseasegroup",
                                                   only.pos = TRUE, min.pct = 0.25,  min.diff.pct = 0.25,
                                                   logfc.threshold = 0.25)



#Load the table with gene descriptions to help you identify to what cell types the genes in the cluster can correspond to
annotations <- read.csv("/home/projects/22102_single_cell/day3/annotation.csv")

# Combine markers with gene descriptions
cluster6_ann_markers <- cluster6_conserved_markers %>%
  rownames_to_column(var="gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))

cluster5_ann_markers <- cluster5_conserved_markers %>%
  rownames_to_column(var="gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))

cluster4_ann_markers <- cluster4_conserved_markers %>%
  rownames_to_column(var="gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))


cluster3_ann_markers <- cluster3_conserved_markers %>%
  rownames_to_column(var="gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))

cluster2_ann_markers <- cluster2_conserved_markers %>%
  rownames_to_column(var="gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))




# Rename all identities

#Cluster 0, 1 and 2 are biologically uncharacterized activated T-cells
pbmc.integrated <- RenameIdents(object = pbmc.integrated,
                                "0" = "TH-ACT1",
                                "1" = "TH-ACT2",
                                "2" = "TH-ACT3",
                                "3" = "TH-1",
                                "4" = "TH-2",
                                "5" = "TH-IFNR",
                                "6" = "TH-17"
                                )

# Create a new column in metadata with cell type annotations from idents
pbmc.integrated@meta.data$Cell_ann <- Idents(pbmc.integrated)

# Now, you can view the updated metadata
View(pbmc.integrated@meta.data)

saveRDS(pbmc.integrated, "Data/Teff_annotation.rds")

# visualize data
clusters <- DimPlot(pbmc.integrated, reduction = 'umap', label = TRUE)
diseasegroup <- DimPlot(pbmc.integrated, reduction = 'umap', group.by = 'diseasegroup')
celltype <- DimPlot(pbmc.integrated, reduction = 'umap', group.by = 'Cell_ann')
donor <- DimPlot(pbmc.integrated, reduction = 'umap', group.by = 'donor')

diseasegroup|
diseasegroup|celltype
donor

clusters

ggsave("umap_cell_annotation.png", clusters, path= "Plots/Teff_integration")

ggsave("umap_integrated_donor.png", donor, path= "Plots/Teff_integration")





