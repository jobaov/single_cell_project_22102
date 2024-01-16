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


# Create an empty list to store results for each cell type
results_list <- list()
for (Cell_ann in names(cts.split.modified)) {
  # Step 1: Extract counts for the current cell type
  counts_current <- cts.split.modified[[Cell_ann]]



  # 2. generate sample level metadata
  colData <- data.frame(samples = colnames(counts_current))
  colData <- colData %>%
    mutate(diseasegroup = ifelse(grepl('Asthm-non-allergic', samples), 'Non-allergic', 'Allergic')) %>%
    column_to_rownames(var = 'samples')


  # Step 3: Perform DESeq2 analysis
  dds <- DESeqDataSetFromMatrix(countData = counts_current,
                                colData = colData,
                                design = ~ diseasegroup)

  # Filter
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]

  # Run DESeq2
  dds <- DESeq(dds)

  # Generate results object
  res <- results(dds, name = "diseasegroup_Non.allergic_vs_Allergic")
  res_tbl <- as.data.frame(res) ## not sure why yet

  # Save the results for the current cell type
  results_list[[Cell_ann]] <- res_tbl
}
##########################################################################################################################
# Data interpretation

# Set thresholds
padj_cutoff <- 0.05
###########################################################################################################################

# Define function
process_data_frame <- function(res_tbl, padj_cutoff = 0.05) {
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

  # Return the results
  return(list(
    top20_genes_by_padj = top20_sig_genes,
    top20_genes_by_change = top20_sig_genes_change
  ))
}

processed_data_teff <- lapply(results_list, process_data_frame)
saveRDS(processed_data_teff, "Data/DGE/Teff/DGE_processed_data_teff.rds")

for (i in seq_along(processed_data_teff)) {
  cat("Results for object", i, ":\n")
  print("Top 20 significant genes ordered by padj:")
  print(processed_data_teff[[i]]$top20_genes_by_padj)

  print("Top 20 significant genes ordered by log fold change:")
  print(processed_data_teff[[i]]$top20_genes_by_change)
  cat("\n")
}
