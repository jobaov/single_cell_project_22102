# Differential cell type abundance - effector T cells
# Last changed: 15/01

# Load libraries
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
library(tidyr)

# Load data
teff_asthm <- readRDS("/net/pupil1/home/projects/Group1/Teff_annotation.rds")

#######################################################################################
# Visual inspection of differential abundances and save plot

dim_diseasegroup <- DimPlot(teff_asthm, split.by = "diseasegroup")
dim_diseasegroup
ggsave("dim_diseasegroup.png", dim_diseasegroup, path= "Plots/DCAT/Teff")

#######################################################################################
# Differential Cell Type abundance using DCAT tool

# Create misclassification matrix
knn_mat = knn_simMat(teff_asthm@graphs$integrated_snn, teff_asthm$Cell_ann)
print(knn_mat)

# Get count matrix  which contains the numbers of cell for each cell type in each sample
teff_asthm$diseasegroup[teff_asthm$diseasegroup == 'AS_NA'] <- 'Asthm_non-allergic'
teff_asthm$diseasegroup[teff_asthm$diseasegroup == 'AS_AL'] <- 'Asthm_allergic'

teff_asthm$id <- paste0(teff_asthm$diseasegroup, teff_asthm$donor)

count_mat = table(teff_asthm$id, teff_asthm$Cell_ann)
count_mat

# Create design dataframe
condition_vector <- rep(c("Asthm_allergic", "Asthm_non_allergic"), each = 6)
teff_design <- data.frame(condition = condition_vector)

# Perform differential abundance analysis
results_DA <- dcats_GLM(count_mat, teff_design, knn_mat)

################################################################################

# Data interpretation and visual inspection

# Extract log-fold changes, p-values, and cell types
log_fold_changes <- results_DA$ceoffs[, 1]
p_values <- results_DA$LRT_pvals[, 1]
cell_types <- rownames(results_DA$ceoffs)

# Calculate -log10(p-values)
neg_log10_p_values <- -log10(p_values)

# Create a dataframe for plotting and save
DCAT_data <- data.frame(
  log_fold_change = log_fold_changes,
  neg_log10_p_value = neg_log10_p_values,
  cell_type = cell_types
)

saveRDS(DCAT_data, "Data/DCAT/Teff/DCAT_result_data.rds")

# Create a volcano plot with labels and save
volcanoplot <- ggplot(DCAT_data, aes(x = log_fold_change, y = neg_log10_p_value, label = cell_type)) +
  geom_point(aes(color = factor(results_DA$fdr[, 1] < 0.05)), size = 3) +
  geom_text_repel(aes(label = cell_type), box.padding = 0.5, point.padding = 0.1, size = 2) +
  scale_color_manual(values = c("black", "red"), guide = FALSE) +
  labs(
    title = "Volcano Plot",
    x = "Log-fold Change",
    y = "-log10(p-value)"
  )
volcanoplot
ggsave("volcanoplot.png", volcanoplot, path= "Plots/DCAT/Teff")

# Create stacked barplots and save
count_mat_diseasegroup <- table(teff_asthm$diseasegroup, teff_asthm$Cell_ann)
count_mat_diseasegroup
saveRDS(count_mat_diseasegroup, "Data/DCAT/Teff/DCAT_stackedplot_data.rds")

# Create a stacked bar plot with each bar representing a condition
# Calculate row wise percentages
percentages <- prop.table(count_mat_diseasegroup, margin = 1) * 100

# Convert data to a data frame for ggplot
df <- data.frame(DiseaseGroup = rep(c("Asthm_allergic", "Asthm_non-allergic"), each = 7), # number has to be adjusted
                 CellType = rep(c("TH-ACT1", "TH-ACT2", "TH-ACT3", "TH-1", "TH-2", "TH-IFNR", "TH-17"), times = 2), ## has to be adjusted!!
                 Percentage = as.vector(t(percentages)))

# Create a stacked barplot using ggplot and save
stacked_bycondition <- ggplot(df, aes(x = DiseaseGroup, y = Percentage, fill = CellType)) +
  geom_bar(stat = "identity") +
  labs(title = "Stacked Barplot of Cell Types by Disease Group (Percentage)",
       x = "Disease Groups",
       y = "Percentage") +
  theme_minimal() +
  scale_fill_manual(values = c("blue", "green", "red","#76EEC6","#CD661D","#FF6EB4","#FFA500")) +
  geom_text(aes(label = paste0(round(Percentage), "%")),
           position = position_stack(vjust = 0.5), color = "white") +
  theme(legend.position = "right", legend.title = element_blank())

stacked_bycondition
ggsave("stacked_bycondition.png", stacked_bycondition, path= "Plots/DCAT/Teff")

# Create a stacked barplot with each bar represeting a celltype
# Calculate column-wise percentages
percentages02 <- prop.table(count_mat_diseasegroup, margin = 2) * 100

# Convert data to a data frame for ggplot
df02 <- data.frame(CellType = rep(c("TH-ACT1", "TH-ACT2", "TH-ACT3", "TH-1", "TH-2", "TH-IFNR", "TH-17"), each = 2),
                   DiseaseGroup = rep(c("Asthm_allergic", "Asthm_non-allergic"), times = 7),
                   Percentage = as.vector(percentages02))

# Create a stacked barplot using ggplot and save
stacked_bycelltype <- ggplot(df02, aes(x = CellType, y = Percentage, fill = DiseaseGroup)) +
  geom_bar(stat = "identity") +
  labs(title = "Stacked Barplot of Cell Types by Condition (Percentage)",
       x = "Cell Types",
       y = "Percentage") +
  theme_minimal() +
  scale_fill_manual(values = c("blue", "green")) +
  geom_text(aes(label = paste0(round(Percentage), "%")),
            position = position_stack(vjust = 0.5), color = "white") +
  theme(legend.position = "top")

stacked_bycelltype
ggsave("stacked_bycelltype.png", stacked_bycelltype, path= "Plots/DCAT/Teff")
