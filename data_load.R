## 22102 Applied Single Cell Bioinformatics
## Project
## Last modified: 12/1 2024

# Load libraries
library(Seurat)
library(tidyverse)
library(dplyr)

# Read in the file
counts <- read.table(file = "Data/GSE146170_allergen_TREG_umi.txt.gz",
                     header = TRUE)

# Make gene names as row names
rownames(counts) <- counts[ , 1]
counts <- counts[ , -1]

# Read in the meta data
meta_data <- read.table(file = "Data/allergen_TREG_annotation.txt",
                        header = TRUE)

# Make gene names as row names
rownames(meta_data) <- meta_data[ , 1]
meta_data <- meta_data[ , -1]

# Extract relevant information
meta_data <- meta_data[ , c(1:2, 6, 8:10)]

# Create Seurat Object
treg_data <- CreateSeuratObject(counts = counts,
                                meta.data = meta_data)
# Save as RDS file
saveRDS(treg_data, "Data/treg_data")
