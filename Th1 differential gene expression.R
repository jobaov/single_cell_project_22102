library(DESeq2)
library(Seurat)
library(SingleCellExperiment)
library(DCATS)
library(SeuratData)
library(tidyverse)
library(Seurat)
library(enrichplot)
library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(tibble)
teff_asthm <- readRDS("/net/pupil1/home/projects/Group1/Teff_annotation.rds")
#Preparing the single-cell dataset for pseudobulk analysis

# Create ID column
teff_asthm$diseasegroup[teff_asthm$diseasegroup == 'AS_NA'] <- 'Asthm_non-allergic'
teff_asthm$diseasegroup[teff_asthm$diseasegroup == 'AS_AL'] <- 'Asthm_allergic'

teff_asthm$id <- paste0(teff_asthm$diseasegroup, teff_asthm$donor)

# Aggregate counts to sample level
counts <- AggregateExpression(teff_asthm,
                              group.by = c("Cell_ann", "id"),
                              assays =  "RNA",
                              return.seurat = FALSE)

counts <- counts$RNA

# transpose
counts.t <- t(counts)

# convert to data.frame
counts.t <- as.data.frame(counts.t)

# get values where to split
splitRows <- gsub('_.*', '', rownames(counts.t))

# split data.frame
cts.split <- split.data.frame(counts.t,
                              f = factor(splitRows))
# fix colnames and transpose
cts.split.modified <- lapply(cts.split, function(x){
  rownames(x) <- gsub('.*_(.*)', '\\1', rownames(x))
  t(x)

})


# Step 1: Extract counts for the current cell type
  counts_th1 <- cts.split.modified$`TH-1`



# 2. generate sample level metadata
  colData <- data.frame(samples = colnames(counts_th1))
  colData <- colData %>%
    mutate(diseasegroup = ifelse(grepl('Asthm-non-allergic', samples), 'Non-allergic', 'Allergic')) %>%
    column_to_rownames(var = 'samples')


# Step 3: Create DESeq2 object
  dds <- DESeqDataSetFromMatrix(countData = counts_th1,
                                colData = colData,
                                design = ~ diseasegroup)

# Quality control

  # Add counts to dds object colData
  colData(dds)$num_cells <- colSums(counts(dds) > 0)

  # Transform counts for data visualization
  rld <- rlog(dds, blind=TRUE)

  # Plot PCA
  PCA_th1 <- DESeq2::plotPCA(rld, ntop = 500, intgroup = "diseasegroup")
  PCA_th1
  ggsave("PCA_th1.png", PCA_th1, path= "Plots/DGE/Teff")

  PCA_th1_counts <- DESeq2::plotPCA(rld, ntop = 500, intgroup = "num_cells")
  PCA_th1_counts
  ggsave("PCA_th1_counts.png", PCA_th1_counts, path= "Plots/DGE/Teff")


  rld_mat <- assay(rld)
  rld_cor <- cor(rld_mat)

  # Plot heatmap without clustering rows or columns
  heatmap <- pheatmap(rld_cor,
           annotation_col = colData[, c("diseasegroup"), drop = FALSE],
           cluster_rows = FALSE,
           cluster_cols = FALSE)
  heatmap
  ggsave("heatmap_th1.png", heatmap, path= "Plots/DGE/Teff")

# Filter
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]

  # Run DESeq2
  dds <- DESeq(dds)

  dispersions <- plotDispEsts(dds)
  dispersions
  ggsave("dispersions_th1.png", dispersions, path= "Plots/DGE/Teff")
  # Generate results object
  res <- results(dds, name = "diseasegroup_Non.allergic_vs_Allergic")
  res_tbl <- as.data.frame(res)



##########################################################################################################################
# Data interpretation

# Set thresholds
padj_cutoff <- 0.05

###########################################################################################################################

  res_tbl <- res_tbl %>%
    rownames_to_column(var = "gene") %>%
    as_tibble() %>%
    arrange(padj)

  # Subset the significant results
  sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) %>%
    dplyr::arrange(padj)

  # Order results by padj values and select top 20
  top20_sig_genes <- sig_res %>%
    dplyr::arrange(padj) %>%
    dplyr::pull(gene) %>%
    head(n = 20)

  # Order results by log fold change and select top 20
  top20_sig_genes_change <- sig_res %>%
    dplyr::arrange(log2FoldChange) %>%
    dplyr::pull(gene) %>%
    head(n = 20)

  # Print top 20 genes for each cell type
  print("Top 20 significant genes ordered by padj:")
  print(top20_sig_genes)

  print("Top 20 significant genes ordered by log fold change:")
  print(top20_sig_genes_change)



ggplot(res_tbl, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = ifelse(abs(log2FoldChange) > 1 & padj < 0.05, "red", "black")), size = 2) +
    scale_color_manual(values = c("black", "red")) +
    theme_minimal() +
    labs(title = "Volcano Plot for Differential Gene Expression",
         x = "log2 Fold Change",
         y = "-log10(padj)")
