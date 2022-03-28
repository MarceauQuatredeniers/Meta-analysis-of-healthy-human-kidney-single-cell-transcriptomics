

# ------------------------------ #
# px2 - Computational analysis of kidney single cell datasets
# The overall goal is to integrate different control samples from various sources to get an average overview of the kidney cell landscape 
# The present script in particular, aims to project cells onto a pseudo-time trajectory using monocle v3 while starting from Seurat objects
# ------------------------------ #


# ------------------------------ #
# INITIALISATION
# Load the necessary packages
library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(Matrix)
library(ggplot2)
library(cowplot)
library(dplyr)
set.seed(1234)
# Then set the path to the relevant folders
RData_dir <- "C:/Users/quatr/Desktop/POSTDOC_U1163/211103_p2_trajectory_inference/data/"
reports_dir <- "C:/Users/quatr/Desktop/POSTDOC_U1163/211103_p2_trajectory_inference/reports/"
data_processed_dir <- "C:/Users/quatr/Desktop/POSTDOC_U1163/211103_p2_trajectory_inference/data/processed/"
figure_dir <- "C:/Users/quatr/Desktop/POSTDOC_U1163/211103_p2_trajectory_inference/reports/figures/"


# ------------------------------ #
# PREPARATION
# ------------------------------ #
sc_ptc <- readRDS(paste(RData_dir, "4A_sc_ptc.rds", sep = "")) # load gsm_merge backup
DefaultAssay(sc_ptc) <- "SCT"
p1 <- DimPlot(sc_ptc, reduction = "umap", group.by = "seurat_clusters", label = T, repel = T, pt.size = 0.1)
p2 <- VlnPlot(sc_ptc, group.by = "seurat_clusters", features = c("CUBN", "ALDOB", "GPX3", "DCXR", "SLC13A3", "SLC34A1", "SLC22A8", "SLC22A7"), stack = T)
p1 + p2
sc_ptc <- subset(sc_ptc, seurat_clusters == "0" | seurat_clusters == "1" | seurat_clusters == "2" | seurat_clusters == "3" | seurat_clusters == "4" | seurat_clusters == "5" | seurat_clusters == "6" | seurat_clusters == "7" | seurat_clusters == "8" | seurat_clusters == "9" | seurat_clusters == "10" | seurat_clusters == "11" | seurat_clusters == "26" | seurat_clusters == "36" | seurat_clusters == "38")
# ------------------------------ #
# TRAJECTORY INFERENCE ANALYSIS
# ------------------------------ #
sc_ptc <- readRDS(paste(RData_dir, "4A_sc_ptc.rds", sep = "")) # load gsm_merge backup
sc_ptc@assays$RNA <- NULL
sc_ptc@assays$SCT <- NULL
DefaultAssay(sc_ptc) <- "seurat.integration"
root <- colnames(subset(sc_ptc, subset = seurat_clusters == "31" | seurat_clusters == "42" | seurat_clusters == "48")) # specify root cells (the starting point of the pseudo-time trajectory)
sc_ptc <- as.cell_data_set(sc_ptc) # transform Seurat S4 object to a CellDataSet object used by monocle3
sc_ptc <- cluster_cells(cds = sc_ptc, reduction_method = "UMAP", cluster_method = "louvain", resolution = 0.0001, num_iter = 3, verbose = T)
sc_ptc <- learn_graph(sc_ptc, use_partition = T, verbose = T)
# Order cells along the trajectory
sc_ptc <- order_cells(sc_ptc, reduction_method = "UMAP", root_cells = root)
# plot trajectories colored by pseudotime
plot_cells(cds = sc_ptc, color_cells_by = "pseudotime", show_trajectory_graph = TRUE)
# ------------------------------ #
# Extract the pseudotime values and add to the Seurat object
sc_ptc_S4 <- readRDS(paste(RData_dir, "4A_sc_ptc.rds", sep = "")) # load gsm_merge backup
sc_ptc_S4 <- AddMetaData(
  object = sc_ptc_S4,
  metadata = sc_ptc@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "sc_ptc"
)
library(viridis)
p2 <- FeaturePlot(sc_ptc_S4, "sc_ptc", pt.size = 0.1, label = T, repel = T) & scale_color_viridis_c() + 
plot_grid(p2, p1, ncol = 2, rel_widths = c(1, 1.05)) %>% ggsave(filename = paste(figure_dir, "6_sc_ptc_SeuratClusters_TI.png", sep = ""), width = 400, height = 180, units = "mm")
VlnPlot(sc_ptc_S4, group.by = "seurat_clusters", features = c("CUBN", "ALDOB", "GPX3", "DCXR", "SLC13A3", "SLC34A1", "SLC22A8", "SLC22A7"), stack = T)
