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
treg_asthm <- readRDS("Data/treg_asthm_filtered.rds")
#Preparing the single-cell dataset for pseudobulk analysis

# Create ID column
treg_asthm$diseasegroup[treg_asthm$diseasegroup == 'AS_NA'] <- 'Asthm_non-allergic'
treg_asthm$diseasegroup[treg_asthm$diseasegroup == 'AS_AL'] <- 'Asthm_allergic'

treg_asthm$id <- paste0(treg_asthm$diseasegroup, treg_asthm$donor)

# Aggregate counts to sample level
counts <- AggregateExpression(immune.sub,
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
  mutate(diseasegroup = ifelse(grepl('AS_NA', samples), 'Non-allergic', 'Allergic')) %>%
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
res <- results(dds, name = "condition_Non-allergic_vs_Allergic")
res_tbl <- as.data.frame(res) ## not sure why yet

# Save the results for the current cell type
results_list[[Cell_ann]] <- res_tbl
}

# Turn the DESeq2 results object into a tibble for use with tidyverse functions
res_tbl<- res2 %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  arrange(padj)

# Check results output
res_tbl

# Set thresholds
padj_cutoff <- 0.005

# Subset the significant results
sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) %>%
  dplyr::arrange(padj)

# Check significant genes output
sig_res

## Order results by padj values
top20_sig_genes <- sig_res %>%
  dplyr::arrange(padj) %>%
  dplyr::pull(gene) %>%
  head(n=20)

## Order results by log fold change
top20_sig_genes_change <- sig_res %>%
  dplyr::arrange(log2FoldChange) %>%
  dplyr::pull(gene) %>%
  head(n=20)
