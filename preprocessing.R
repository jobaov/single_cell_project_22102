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

##### Compute commonly used QC metrics.

for (i in 1:2){

  # Calculate novelty score
  data_asthm[[i]]$log10GenesPerUMI <- log10(data_asthm[[i]]$nFeature_RNA) / log10(data_asthm[[i]]$nCount_RNA)

  # Update nCount_RNA
  data_asthm[[i]]$nCount_RNA <- colSums(data_asthm[[i]]@assays$RNA$counts)

  # Update nFeatures_RNA
  data_asthm[[i]]$nFeature_RNA <- colSums(data_asthm[[i]]@assays$RNA$counts != 0)

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

##### Plot the results (violin plots).

plot_vln_treg <- VlnPlot(treg_asthm,
                         features = c("nFeature_RNA",
                                      "nCount_RNA",
                                      "percent_mt",
                                      "percent_ribo",
                                      "percent_hb"),
                         cols = "#6699CC",
                         combine = FALSE)
plot_vln_teff <- VlnPlot(teff_asthm,
                         features = c("nFeature_RNA", "nCount_RNA", "percent_mt", "percent_ribo", "percent_hb"),
                         cols = "#6699CC",
                         combine = FALSE)

plot_vln_list <- list(plot_vln_treg, plot_vln_teff)

title <- c("Number of genes",
           "Number RNA molecules",
           "Mitochondrial content (%)",
           "Ribosomal content (%)",
           "Hemoglobin content (%)")

threshold_treg <- list(c(250, 1600), 500, 5, 0, 1)
threshold_teff <- list(c(250, 1600), 500, 5, 0, 1)
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
            legend.position = "none",
            plot.title = element_text(size = 5, face = "bold"))
  }
  plot_vln_list[[i]] <- plot_vln
}

plot_vln_treg <- wrap_plots(plot_vln_list[[1]], # Remove line ribosomal
                            ncol = 5)
plot_vln_teff <- wrap_plots(plot_vln_list[[2]], # Mito to 5
                            ncol = 5)

ggsave("plot_vln_treg.png", plot_vln_treg, path= "Plots/QC")
ggsave("plot_vln_teff.png", plot_vln_teff, path= "Plots/QC")

##### Plot the results (density plots).

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
    xlab("Number of genes") +
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
    xlab("Novelty score") +
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
    xlab("Mitochondrial content") +
    ylab("Cell density") +
    geom_vline(xintercept = 5) +
    scale_fill_discrete(name = "Disease group", labels = c("Asthmatic allergic",
                                                           "Asthmatic non-allergic")) +
    guides(color = FALSE)

  plot_dens_list[[i]] <- plot_dens
}

##### Plot the results (scatter plots).

plot_scatter_list <- list()

for (i in 1:2){

  plot_scatter <- list()

  # Visualize number of RNA molecules vs. number of genes.
  plot_scatter[[1]] <- FeatureScatter(data_asthm[[i]],
                                      feature1 = "nCount_RNA",
                                      feature2 = "nFeature_RNA") +
    xlab("Number of RNA molecules") +
    ylab("Number of genes") +
    theme(legend.position = "none")

  # Visualize number of RNA molecules vs. percent mitochondrial content.
  plot_scatter[[2]] <- FeatureScatter(data_asthm[[i]],
                                      feature1 = "nCount_RNA",
                                      feature2 = "percent_mt") +
    xlab("Number of RNA molecules") +
    ylab("Mitochondrial content")

  plot_scatter_list[[i]] <- plot_scatter
}

# Save figures
ggsave("plot_scatter_treg.png", plot_scatter_list[[1]][[1]], path= "Plots/QC")
ggsave("plot_scatter_teff.png", plot_scatter_list[[2]][[1]], path= "Plots/QC")

# Filter the cells based on your QC metrics.
treg_asthm <- subset(treg_asthm,
                     subset = nFeature_RNA > 250 & nFeature_RNA < 1600 & percent_mt < 5 & percent_hb < 1 & log10GenesPerUMI > 0.8 & nCount_RNA > 500)

# Filter out ribosomal genes.
treg_asthm <- treg_asthm[!grepl('^RP[SL]', rownames(treg_asthm)), ]

# Filter the cells based on your QC metrics.
teff_asthm <- subset(teff_asthm,
                     subset = nFeature_RNA > 250 & nFeature_RNA < 1600 & percent_mt < 5 & percent_hb < 1 & log10GenesPerUMI > 0.8 & nCount_RNA > 500)

# Filter out ribosomal genes.
teff_asthm <- teff_asthm[!grepl('^RP[SL]', rownames(teff_asthm)), ]

# Extract counts
counts_treg <- GetAssayData(object = treg_asthm,
                            slot = "counts")
counts_teff <- GetAssayData(object = teff_asthm,
                            slot = "counts")

# Output a logical matrix specifying for each gene on whether or not there are more than zero counts per cell
nonzero_treg <- counts_treg > 0
nonzero_teff <- counts_teff > 0

# Sums all TRUE values and returns TRUE if more than 4 TRUE values per gene
keep_genes_treg <- Matrix::rowSums(nonzero_treg) >= 4
keep_genes_teff <- Matrix::rowSums(nonzero_teff) >= 4

# Only keeping those genes expressed in more than 4 cells
filtered_counts_treg <- counts_treg[keep_genes_treg, ]
filtered_counts_teff <- counts_teff[keep_genes_teff, ]

# Reassign to filtered Seurat object
treg_asthm <- CreateSeuratObject(filtered_counts_treg,
                                 meta.data = treg_asthm@meta.data)
teff_asthm <- CreateSeuratObject(filtered_counts_teff,
                                 meta.data = teff_asthm@meta.data)

# Standard workflow
# Do a log-normalization following Seurat’s standard workflow.
treg_asthm <- NormalizeData(treg_asthm,
                            normalization.method = "LogNormalize",
                            scale.factor = 10000)
teff_asthm <- NormalizeData(teff_asthm,
                            normalization.method = "LogNormalize",
                            scale.factor = 10000)

# Identifying highly variable features.
treg_asthm <- FindVariableFeatures(treg_asthm,
                                   selection.method = "vst",
                                   nfeatures = 2000)
teff_asthm <- FindVariableFeatures(teff_asthm,
                                   selection.method = "vst",
                                   nfeatures = 2000)
# Plot variable features.
VariableFeaturePlot(treg_asthm) +
  scale_y_continuous(limits = c(0, 10))
top10_treg <- head(VariableFeatures(treg_asthm), 10)

# Plot variable features with labels
plot1_treg <- VariableFeaturePlot(treg_asthm)
plot2_treg <- LabelPoints(plot = plot1_treg, points = top10_treg, repel = TRUE)
plot2_treg

VariableFeaturePlot(teff_asthm) +
  scale_y_continuous(limits = c(0, 10))
top10_teff <- head(VariableFeatures(teff_asthm), 10)

# Plot variable features with labels
plot1_teff <- VariableFeaturePlot(teff_asthm)
plot2_teff <- LabelPoints(plot = plot1_teff, points = top10_teff, repel = TRUE)
plot2_teff

ggsave("plot_varfeature_treg.png", plot2_treg, path= "Plots/QC")
ggsave("plot_varfeature_teff.png", plot2_teff, path= "Plots/QC")

# Scale the data.
treg_asthm <- ScaleData(treg_asthm,
                        features = VariableFeatures(object = treg_asthm))
teff_asthm <- ScaleData(teff_asthm,
                        features = VariableFeatures(object = teff_asthm))

# Perform linear dimensionality reduction.
treg_asthm <- RunPCA(treg_asthm,
                     features = VariableFeatures(object = treg_asthm))
teff_asthm <- RunPCA(teff_asthm,
                     features = VariableFeatures(object = teff_asthm))

treg_asthm <- RunUMAP(treg_asthm, dims = 1:30)
teff_asthm <- RunUMAP(teff_asthm, dims = 1:30)

# Convert your seurat object to a sce object.
treg_sce <- as.SingleCellExperiment(treg_asthm)
teff_sce <- as.SingleCellExperiment(teff_asthm)

# Check for doublets.
treg_top.var <- VariableFeatures(treg_asthm)
teff_top.var <- VariableFeatures(teff_asthm)

treg_dbl.dens <- computeDoubletDensity(treg_sce,
                                       subset.row = treg_top.var,
                                       d = ncol(reducedDim(treg_sce)))
teff_dbl.dens <- computeDoubletDensity(teff_sce,
                                       subset.row = teff_top.var,
                                       d = ncol(reducedDim(teff_sce)))
treg_sce$DoubletScore <- treg_dbl.dens
teff_sce$DoubletScore <- teff_dbl.dens

treg_dbl.calls <- doubletThresholding(data.frame(score = treg_dbl.dens),
                                      method ="griffiths",
                                      returnType ="call")
teff_dbl.calls <- doubletThresholding(data.frame(score = teff_dbl.dens),
                                      method ="griffiths",
                                      returnType ="call")
summary(treg_dbl.calls)
summary(teff_dbl.calls)

plot_dscore_treg <- plotUMAP(treg_sce, colour_by="DoubletScore")
plot_dscore_teff <- plotUMAP(teff_sce, colour_by="DoubletScore")

ggsave("plot_dscore_treg.png", plot_dscore_treg, path= "Plots/QC")
ggsave("plot_dscore_teff.png", plot_dscore_teff, path= "Plots/QC")

# Remove doublets
# Extract singlet cell indices
treg_singlet_indices <- which(treg_dbl.calls == "singlet")
teff_singlet_indices <- which(teff_dbl.calls == "singlet")

# Filter out doublets from the original Seurat object
treg_asthm <- treg_asthm[, treg_singlet_indices]
teff_asthm <- teff_asthm[, teff_singlet_indices]

saveRDS(treg_asthm, "Data/treg_asthm_filtered.rds")
saveRDS(teff_asthm, "Data/teff_asthm_filtered.rds")
