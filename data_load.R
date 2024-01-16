## 22102 Applied Single Cell Bioinformatics
## Project
## Last modified: 13/1 2024

## DATA LOAD

###########################################

# Load libraries
library(Seurat)
library(tidyverse)
library(dplyr)

# Read in the files
counts_treg <- read.table(file = "Data/GSE146170_allergen_TREG_umi.txt.gz",
                          header = TRUE)
counts_teff <- read.table(file = "Data/GSE146170_allergen_TEFF_umi.txt.gz",
                          header = TRUE)

# Make gene names as row names
rownames(counts_treg) <- counts_treg[ , 1]
counts_treg <- counts_treg[ , -1]
rownames(counts_teff) <- counts_teff[ , 1]
counts_teff <- counts_teff[ , -1]

# Read in the meta data
meta_treg <- read.table(file = "Data/allergen_TREG_annotation.txt",
                        header = TRUE)
meta_teff <- read.table(file = "Data/allergen_TEFF_annotation.txt",
                        header = TRUE)

# Make cell ID as row names
rownames(meta_treg) <- meta_treg[ , 1]
meta_treg <- meta_treg[ , -1]
rownames(meta_teff) <- meta_teff[ , 1]
meta_teff <- meta_teff[ , -1]

# Extract relevant information
meta_treg <- meta_treg[ , c(1:2, 6:10)]
meta_teff <- meta_teff[ , c(1:2, 6:10)]

# Create nCount and nFeature metadata
nCount_treg <- colSums(counts_treg)
nCount_teff <- colSums(counts_teff)

nFeatures_treg <- colSums(counts_treg != 0)
nFeatures_teff <- colSums(counts_teff != 0)

#meta_treg <- cbind(meta_treg, nCount_treg, nFeatures_treg)
#meta_teff <- cbind(meta_teff, nCount_teff, nFeatures_teff)

# Create Seurat Object
treg_data <- CreateSeuratObject(counts = counts_treg,
                                meta.data = meta_treg)
teff_data <- CreateSeuratObject(counts = counts_teff,
                                meta.data = meta_teff)

# Subset the Seurat object to only analyze cells from
# Control: asthmatic non-allergic (AS_NA)
# Diseased: asthmatic allergic (AS_AL)
treg_asthm <- subset(treg_data, subset = diseasegroup %in% c("AS_NA", "AS_AL"))
teff_asthm <- subset(teff_data, subset = diseasegroup %in% c("AS_NA", "AS_AL"))

# Save all 4 Seurat objects as RDS files
saveRDS(treg_data, "Data/treg_data.rds")
saveRDS(teff_data, "Data/teff_data.rds")
saveRDS(treg_asthm, "Data/treg_asthm.rds")
saveRDS(teff_asthm, "Data/teff_asthm.rds")
