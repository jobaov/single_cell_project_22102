# Load libraries
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(clustree)
library(scDblFinder)
library(SingleCellExperiment)

# Read Seurat object
treg <- readRDS("Data/treg_data")

# Subset and analyze cells from asthmatic non-allergic "AS_NA and asthmatic allergic "AS_AL"
asthm <- subset(treg, subset = diseasegroup %in% c("AS_NA", "AS_AL"))


# --------- Quality Control----------------------#

# Calculate novelty score
asthm$log10GenesPerUMI <- log10(asthm$nFeature_RNA) / log10(asthm$nCount_RNA)

# Mitochondrial content
asthm[["percent.mt"]] <- PercentageFeatureSet(asthm, pattern = "^MT-")

# Ribosomal proportion
asthm <- PercentageFeatureSet(asthm, "^RP[SL]", col.name = "percent_ribo")

# Hemoglobin proportion
asthm <- PercentageFeatureSet(asthm, "^HB[^(P)]", col.name = "percent_hb")

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

#  Visualize the number of transcripts per cell in each diseasegroup.
UMI_cell <- asthm[[]] %>%
  ggplot(aes(color= diseasegroup , x= nCount_RNA, fill=diseasegroup)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)


UMI_cell

ggsave("UMI_cell.png", UMI_cell, path= "../Plots/QC")

# Visualize the distribution of genes detected per cell via histogram
genes_cell <- asthm[[]] %>%
  ggplot(aes(color=diseasegroup, x=nFeature_RNA, fill= diseasegroup)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  scale_x_log10() +
  geom_vline(xintercept = 300)

genes_cell

ggsave("genes_cell.png", genes_cell, path= "../Plots/QC")

# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI (novelty score)
novelty <- asthm[[]] %>%
  ggplot(aes(x=log10GenesPerUMI, color = diseasegroup, fill=diseasegroup)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)

novelty

ggsave("novelty.png", novelty, path= "../Plots/QC")

# Visualize the distribution of mitochondrial gene expression detected per cell
mito <- asthm[[]] %>%
  ggplot(aes(color=diseasegroup, x=percent.mt, fill=diseasegroup)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  geom_vline(xintercept = 5)

mito

ggsave("mito.png", mito, path= "../Plots/QC")

#  Visualize the created QC metrics as a violin plot.
violinQC <- VlnPlot(asthm, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent_hb", "percent_ribo"), ncol = 5)
violinQC

# Visualize count vs feature
count_feature <- FeatureScatter(asthm, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
count_feature

# Visualize count vs mito content
count_mito <- FeatureScatter(asthm, feature1 = "nCount_RNA", feature2 = "percent.mt")
count_mito

# Filter the cells based on your QC metrics.
asthm <- subset(asthm, subset = nFeature_RNA > 200 & nFeature_RNA < 1500 & percent.mt < 5 & percent_hb < 1 & log10GenesPerUMI > 0.8 & nCount_RNA > 500)

# Gene level filtering
## filter out ribosomal genes
asthm<- asthm[ ! grepl('^RP[SL]', rownames(asthm)), ]

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
