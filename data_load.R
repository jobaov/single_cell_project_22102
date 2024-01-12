# Load libraries
library(Seurat)
library(dplyr)

# Read in entire file

counts <- read.delim("data/GSE146170_allergen_TREG_umi.txt.gz",
                     header = TRUE)

# Pull only gene information
tirosh_genes <- tirosh[-1:-3,]

# Duplicate gene names so make names unique (up to you how you want to deal with this part)
gene_list <- tirosh_genes %>%
  pull("Cell") %>%
  make.unique(sep = ".")

# Add back unique rownames
rownames(tirosh_genes) <- gene_list

# Remove Column of gene names
tirosh_genes <- tirosh_genes[, -1]

# Pull meta data columns from original data
tirosh_meta <- tirosh[1:3,]

# Make rownames equal to column 1 values
rownames(tirosh_meta) <- tirosh_meta[, 1]

# Remove column 1
tirosh_meta <- tirosh_meta[, -1]

# Transpose meta data as Seurat expects meta data to have cell names as rows and meta data values as columns
tirosh_meta_transpose <- data.frame(t(tirosh_meta))

# Create Seurat Object
tirosh_seurat <- CreateSeuratObject(counts = tirosh_genes, meta.data = tirosh_meta_transpose)
View(tirosh_seurat@meta.data)
