# Differential cell type abundance - regulatory T cells
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
treg_asthm <- readRDS("/net/pupil1/home/projects/Group1/treg_annotated.rds")

#######################################################################################
# Visual inspection of differential abundances and save plot

dim_diseasegroup <- DimPlot(treg_asthm, split.by = "diseasegroup")
dim_diseasegroup
ggsave("dim_diseasegroup.png", dim_diseasegroup, path= "Plots/DCAT/Treg", width = 10, height = 7.5, units = "in", dpi = 300)

#######################################################################################
# Differential Cell Type abundance using DCAT tool

# Create misclassification matrix
knn_mat = knn_simMat(treg_asthm@graphs$integrated_snn, treg_asthm$Cell_ann)
print(knn_mat)

# Get count matrix  which contains the numbers of cell for each cell type in each sample
treg_asthm$diseasegroup[treg_asthm$diseasegroup == 'AS_NA'] <- 'Asthm_non-allergic'
treg_asthm$diseasegroup[treg_asthm$diseasegroup == 'AS_AL'] <- 'Asthm_allergic'

treg_asthm$id <- paste0(treg_asthm$diseasegroup, treg_asthm$donor)

count_mat = table(treg_asthm$id, treg_asthm$Cell_ann)
count_mat

# Create design dataframe
condition_vector <- rep(c("Asthm_allergic", "Asthm_non_allergic"), each = 6)
treg_design <- data.frame(condition = condition_vector)

# Perform differential abundance analysis
results_DA <- dcats_GLM(count_mat, treg_design, knn_mat)

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

saveRDS(DCAT_data, "Data/DCAT/Treg/DCAT_result_data.rds")
write_xlsx(DCAT_data, path = "Data/DCAT/Treg/DCAT_results.xlsx")
# Create a volcano plot with labels and save
volcanoplot <- ggplot(DCAT_data, aes(x = log_fold_change, y = neg_log10_p_value, label = cell_type)) +
  geom_point(aes(color = factor(results_DA$fdr[, 1] < 0.05)), size = 3) +
  geom_text_repel(aes(label = cell_type), box.padding = 0.5, point.padding = 0.1, size = 4) +
  scale_color_manual(values = c("black", "red"), guide = FALSE) +
  labs(
    title = "Allergy vs. Non-allergy",
    x = "Log-fold Change",
    y = "-log10(p-value)"
  ) +
  scale_x_continuous(limits = c(-0.4, 0.4), breaks = seq(-0.4, 0.4, 0.2)) +
  scale_y_continuous(limits = c(0, 4), breaks = seq(0, 3, 1)) +
  theme(title = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 12),)
volcanoplot
ggsave("volcanoplot.png", volcanoplot, path= "Plots/DCAT/Treg")

# Create stacked barplots and save
count_mat_diseasegroup <- table(treg_asthm$diseasegroup, treg_asthm$Cell_ann)
count_mat_diseasegroup
saveRDS(count_mat_diseasegroup, "Data/DCAT/Treg/DCAT_stackedplot_data.rds")

# Create a stacked bar plot with each bar representing a condition
# Calculate row wise percentages
percentages <- prop.table(count_mat_diseasegroup, margin = 1) * 100

# Convert data to a data frame for ggplot
df <- data.frame(DiseaseGroup = rep(c("Asthm_allergic", "Asthm_non-allergic"), each = 3),
                 CellType = rep(c("TregIFNR", "TregACT1", "TregACT2"), times = 2),
                 Percentage = as.vector(t(percentages)))

# Create a stacked barplot using ggplot and save
stacked_bycondition <- ggplot(df, aes(x = DiseaseGroup, y = Percentage, fill = CellType)) +
  geom_bar(stat = "identity") +
  labs(x = NULL,
       y = "Percentage") +
  theme_minimal() +
  scale_fill_manual(values = c("#87CEFF", "#15CD29", "#EE6363")) +
  geom_text(aes(label = paste0(round(Percentage), "%")),
            position = position_stack(vjust = 0.5), color = "white") +
  theme(legend.position = "right",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        axis.text.x = element_text(size = 12))


stacked_bycondition
ggsave("stacked_bycondition.png", stacked_bycondition, path= "Plots/DCAT/Treg",width = 10, height = 7.5, units = "in", dpi = 300)
#############################
# Stacked barplot with total counts

df_tot <- data.frame(DiseaseGroup = rep(c("Asthm_allergic", "Asthm_non-allergic"), each = 3),
                     CellType = rep(c("TregIFNR", "TregACT1", "TregACT2"), times = 2),
                     Counts = as.vector(t(count_mat_diseasegroup)))

# Create a stacked barplot using ggplot and save
stacked_bycondition_tot <- ggplot(df_tot, aes(x = DiseaseGroup, y = Counts, fill = CellType)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(y = "Counts",
       x = NULL) +
  theme_minimal() +
  scale_fill_manual(values = c("#87CEFF", "#15CD29", "#EE6363")) +
  theme(legend.position = "right",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        axis.text.x = element_text(size = 12))

stacked_bycondition_tot
ggsave("stacked_bycondition_tot.png", stacked_bycondition_tot, path= "Plots/DCAT/Treg",width = 10, height = 7.5, units = "in", dpi = 300)
#############################
# Create a stacked barplot with each bar representing a celltype
# Calculate column-wise percentages
percentages02 <- prop.table(count_mat_diseasegroup, margin = 2) * 100

# Convert data to a data frame for ggplot
df02 <- data.frame(CellType = rep(c("TregIFNR", "TregACT1", "TregACT2"), each = 2),
                 DiseaseGroup = rep(c("Asthm_allergic", "Asthm_non-allergic"), times = 3),
                 Percentage = as.vector(percentages02))

# Create a stacked barplot using ggplot and save
stacked_bycelltype <- ggplot(df02, aes(x = CellType, y = Percentage, fill = DiseaseGroup)) +
  geom_bar(stat = "identity") +
  labs(x = NULL,
       y = "Percentage") +
  theme_minimal() +
  scale_fill_manual(values = c("dodgerblue3", "#CD6889")) +
  geom_text(aes(label = paste0(round(Percentage), "%")),
            position = position_stack(vjust = 0.5), color = "white") +
  theme(legend.position = "top",
        legend.title = element_text( size = 12),
        legend.text = element_text(size = 10),
        axis.text.x = element_text(size = 10))

stacked_bycelltype
ggsave("stacked_bycelltype.png", stacked_bycelltype, path= "Plots/DCAT/Treg",width = 10, height = 7.5, units = "in", dpi = 300)

###################################
# With total counts

df02_tot <- data.frame(CellType = rep(c("TregIFNR", "TregACT1", "TregACT2"), each = 2),
                   DiseaseGroup = rep(c("Asthm_allergic", "Asthm_non-allergic"), times = 3),
                   Counts = as.vector(count_mat_diseasegroup))


# Create a stacked barplot using ggplot and save
stacked_bycelltype_tot <- ggplot(df02_tot, aes(x = CellType, y = Counts, fill = DiseaseGroup)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = NULL,
       y = "Counts") +
  theme_minimal() +
  scale_fill_manual(values = c("dodgerblue3", "#CD6889")) +
  theme(legend.position = "top",
        legend.title = element_text( size = 12),
        legend.text = element_text(size = 10),
        axis.text.x = element_text(size = 10))

stacked_bycelltype_tot
ggsave("stacked_bycelltype_tot.png", stacked_bycelltype_tot, path= "Plots/DCAT/Treg",width = 10, height = 7.5, units = "in", dpi = 300)

