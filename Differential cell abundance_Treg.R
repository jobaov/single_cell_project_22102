treg_asthm <- readRDS("Data/treg_asthm_annotated.rds")
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

## ------ Create misclassification matrix
knn_mat = knn_simMat(treg_asthm@graphs$integrated_snn, treg_asthm$Cell_ann)
print(knn_mat)

## ------ Get count matrix  which contains the numbers of cell for each cell type in each sample
treg_asthm$diseasegroup[treg_asthm$diseasegroup == 'AS_NA'] <- 'Asthm_non-allergic'
treg_asthm$diseasegroup[treg_asthm$diseasegroup == 'AS_AL'] <- 'Asthm_allergic'

treg_asthm$id <- paste0(treg_asthm$diseasegroup, treg_asthm$donor)

count_mat = table(treg_asthm$id, treg_asthm$Cell_ann) ## update cell type annotation
count_mat

## ------ Create design dataframe
condition_vector <- rep(c("Asthm_allergic", "Asthm_non_allergic"), each = 6)

# Create the design dataframe
treg_design <- data.frame(condition = condition_vector)

# Perform differential abundance analysis
results_DA <- dcats_GLM(count_mat, treg_design, knn_mat)

# --- Data interpretation
# Extract log-fold changes, p-values, and cell types
log_fold_changes <- results_DA$ceoffs[, 1]
p_values <- results_DA$LRT_pvals[, 1]
cell_types <- rownames(results_DA$ceoffs)

# Calculate -log10(p-values)
neg_log10_p_values <- -log10(p_values)

# Create a dataframe for plotting
volcano_data <- data.frame(
  log_fold_change = log_fold_changes,
  neg_log10_p_value = neg_log10_p_values,
  cell_type = cell_types
)

# Create a volcano plot with labels
ggplot(volcano_data, aes(x = log_fold_change, y = neg_log10_p_value, label = cell_type)) +
  geom_point(aes(color = factor(results_DA$fdr[, 1] < 0.05)), size = 3) +
  geom_text_repel(aes(label = cell_type), box.padding = 0.5, point.padding = 0.1, size = 2) +
  scale_color_manual(values = c("black", "red"), guide = FALSE) +
  labs(
    title = "Volcano Plot",
    x = "Log-fold Change",
    y = "-log10(p-value)"
  )

count_mat_diseasegroup = table(treg_asthm$diseasegroup, treg_asthm$Cell_ann)
count_mat_diseasegroup
