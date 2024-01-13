## 22102 Applied Single Cell Bioinformatics
## Project
## Last modified: 13/1 2024

## PREPROCESSING

###########################################

# Load libraries
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(clustree)
library(scDblFinder)
library(SingleCellExperiment)

# Read Seurat objects
treg_asthm <- readRDS("Data/treg_asthm.rds")
teff_asthm <- readRDS("Data/teff_asthm.rds")
data_asthm <- list(treg_asthm, teff_asthm)

# ---------------------- Quality Control---------------------- #

# Compute commonly used QC metrics.

for (i in 1:2){
  # Calculate novelty score
  data_asthm[[i]]$log10GenesPerUMI <- log10(data_asthm[[i]]$nFeature_RNA) / log10(data_asthm[[i]]$nCount_RNA)

  # Mitochondrial content
  data_asthm[[i]] <- PercentageFeatureSet(data_asthm[[i]],
                                          pattern = "^MT-",
                                          col.name = "percent_mt")

  # Ribosomal content
  data_asthm[[i]] <- PercentageFeatureSet(data_asthm[[i]],
                                          pattern = "^RP[SL]",
                                          col.name = "percent_ribo")

  # Hemoglobin content
  data_asthm[[i]] <- PercentageFeatureSet(data_asthm[[i]],
                                          pattern = "^HB[^(P)]",
                                          col.name = "percent_hb")
}

# Overwrite the old objects to the ones saved in the list.
treg_asthm <- data_asthm[[1]]
teff_asthm <- data_asthm[[2]]

# Plot the results.
plot_vln_treg <- VlnPlot(treg_asthm,
                         features = c("nFeature_RNA", "nCount_RNA", "percent_mt", "percent_ribo", "percent_hb"),
                         cols = "#6699CC",
                         combine = FALSE)
plot_vln_teff <- VlnPlot(teff_asthm,
                         features = c("nFeature_RNA", "nCount_RNA", "percent_mt", "percent_ribo", "percent_hb"),
                         cols = "#6699CC",
                         combine = FALSE)
plot_vln_list <- list(plot_vln_treg, plot_vln_teff)

title <- c("Number of genes", "Number RNA molecules", "Mitochondrial content (%)", "Ribosomal content (%)", "Hemoglobin content (%)")

threshold_treg <- list(c(250, 1500), 500, 5, 1, 1)
threshold_teff <- list(c(250, 1750), 500, 5, 1, 1)
threshold_list <- list(threshold_treg, threshold_teff)

for (i in 1:2){
  plot_vln <- plot_vln_list[[i]]
  threshold <- threshold_list[[i]]
  for (j in 1:5){
    plot_vln[[j]] <- plot_vln[[j]] +
      geom_hline(yintercept = threshold[[j]],
                 col = "#CC3333") +
      ggtitle(title[j]) +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            legend.position = "none")
  }
  plot_vln_list[[i]] <- plot_vln
}

plot_vln_treg <- wrap_plots(plot_vln_list[[1]],
                            ncol = 5)
plot_vln_teff <- wrap_plots(plot_vln_list[[2]],
                            ncol = 5)

plot_dens_list <- list()

for (i in 1:2){
  plot_dens <- list()

  # Visualize the distribution of genes.
  plot_dens[[1]] <- ggplot(data = data_asthm[[i]]@meta.data,
                           aes(x = nFeature_RNA,
                               color = diseasegroup,
                               fill = diseasegroup)) +
    geom_density(alpha = 0.2) +
    scale_x_log10() +
    theme_classic() +
    xlab("RNA molecules") +
    ylab("Cell density") +
    geom_vline(xintercept = 500) +
    scale_fill_discrete(name = "Disease group", labels = c("Asthmatic allergic",
                                                           "Asthmatic non-allergic")) +
    guides(color = FALSE)

  # Visualize the distribution of RNA molecules.
  plot_dens[[2]] <- ggplot(data = data_asthm[[i]]@meta.data,
                           aes(x = nCount_RNA,
                               color = diseasegroup,
                               fill = diseasegroup)) +
    geom_density(alpha = 0.2) +
    scale_x_log10() +
    theme_classic() +
    xlab("RNA molecules") +
    ylab("Cell density") +
    geom_vline(xintercept = 500) +
    scale_fill_discrete(name = "Disease group", labels = c("Asthmatic allergic",
                                                            "Asthmatic non-allergic")) +
    guides(color = FALSE)

  # Visualize the overall complexity by genes detected per UMI (novelty score)
  plot_dens[[3]] <- ggplot(data = data_asthm[[i]]@meta.data,
                           aes(x = log10GenesPerUMI,
                               color = diseasegroup,
                               fill = diseasegroup)) +
    geom_density(alpha = 0.2) +
    scale_x_log10() +
    theme_classic() +
    xlab("RNA molecules") +
    ylab("Cell density") +
    geom_vline(xintercept = 0.8) +
    scale_fill_discrete(name = "Disease group", labels = c("Asthmatic allergic",
                                                           "Asthmatic non-allergic")) +
    guides(color = FALSE)

  # Visualize the distribution of mitochondrial gene expression detected per cell.
  plot_dens[[4]] <- ggplot(data = data_asthm[[i]]@meta.data,
                           aes(x = percent_mt,
                               color = diseasegroup,
                               fill = diseasegroup)) +
    geom_density(alpha = 0.2) +
    scale_x_log10() +
    theme_classic() +
    xlab("RNA molecules") +
    ylab("Cell density") +
    geom_vline(xintercept = 5) +
    scale_fill_discrete(name = "Disease group", labels = c("Asthmatic allergic",
                                                           "Asthmatic non-allergic")) +
    guides(color = FALSE)

  plot_dens_list[[i]] <- plot_dens
}

plot_scatter_list <- list()

for (i in 1:2){
  plot_scatter <- list()
  plot_scatter[[1]] <- FeatureScatter(data_asthm[[i]],
                                      feature1 = "nCount_RNA",
                                      feature2 = "nFeature_RNA")
  plot_scatter[[2]] <- FeatureScatter(data_asthm[[i]],
                                      feature1 = "nCount_RNA",
                                      feature2 = "percent_mt")
  plot_scatter_list[[i]] <- plot_scatter
}

# SAVE RELEVANT FIGURES
#ggsave("UMI_cell.png", UMI_cell, path= "../Plots/QC")
#ggsave("genes_cell.png", genes_cell, path= "../Plots/QC")
#ggsave("novelty.png", novelty, path= "../Plots/QC")
#ggsave("mito.png", mito, path= "../Plots/QC")


# Filter the cells based on your QC metrics.
#asthm <- subset(asthm, subset = nFeature_RNA > 200 & nFeature_RNA < 1500 & percent.mt < 5 & percent_hb < 1 & log10GenesPerUMI > 0.8 & nCount_RNA > 500)

# Gene level filtering
## filter out ribosomal genes
#asthm<- asthm[ ! grepl('^RP[SL]', rownames(asthm)), ]

# Extract counts
counts <- GetAssayData(object = asthm, slot = "counts")

# Output a logical matrix specifying for each gene on whether or not there are more than zero counts per cell
nonzero <- counts > 0

# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 4

# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]

# Reassign to filtered Seurat object
asthm <- CreateSeuratObject(filtered_counts, meta.data = asthm@meta.data)











# Do a log-normalization following Seurat’s standard workflow.
asthm <- NormalizeData(asthm, normalization.method = "LogNormalize", scale.factor = 10000)

# Identifying highly variable features.
asthm <- FindVariableFeatures(asthm, selection.method = "vst", nfeatures = 2000)

# Plot variable features.
VariableFeaturePlot(asthm) +
  scale_y_continuous(limits = c(0, 10))

# Scale the data.
#Scale the data
asthm <- ScaleData(asthm, features = VariableFeatures(object = asthm))

# Perform linear dimensionality reduction.
asthm<- RunPCA(asthm, features = VariableFeatures(object = asthm))

# Elbow plot. The “elbow point” determines the dimensionality of the data.
ElbowPlot(asthm)

# Perform clustering of cells.
asthm <- FindNeighbors(asthm, dims = 1:10)
asthm<- FindClusters(asthm, resolution = c(0.3, 0.5, 0.7))
clustree(asthm)
DimPlot(asthm, group.by = "RNA_snn_res.0.3", label = TRUE)

#  Perform a non-linear reduction(UMAP).
asthm <- RunUMAP(asthm, dims = 1:10)

# Plot UMAP grouped by clusters.
DimPlot(asthm, group.by = "RNA_snn_res.0.3", label = TRUE)

# Determine resolution for clustering.
Idents(asthm) <- "RNA_snn_res.0.3"

# Convert your seurat object to a sce object.
sce <- as.SingleCellExperiment(asthm)

# Check for doublets.
top.var <- VariableFeatures(asthm)
dbl.dens <- computeDoubletDensity(sce, subset.row=top.var, d=ncol(reducedDim(sce)))
sce$DoubletScore <- dbl.dens
dbl.calls <- doubletThresholding(data.frame(score=dbl.dens), method="griffiths", returnType="call")
summary(dbl.calls)

### Prior plots ###

data_asthm

# Visualize the number of cell counts per diseasegroup
cell_counts <- asthm[[]] %>%
  ggplot(aes(x=diseasegroup, fill=diseasegroup)) +
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")

cell_counts

ggsave("Cell_counts.png", cell_counts, path= "../Plots/QC")
