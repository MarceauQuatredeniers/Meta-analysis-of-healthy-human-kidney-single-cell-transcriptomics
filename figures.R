# ------------------------------ #
# Computational analysis of kidney single cell datasets
# The overall goal is to integrate different control samples from various sources to get an average overview of the kidney cell landscape 
# ------------------------------ #
# INITIALISATION
# Load the necessary packages
library(Rcpp)
library(fgsea)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(Seurat)
library(readr)
library(cowplot)
library(pheatmap)
library(RColorBrewer)
library(ggplotify)
# Then set the path to the relevant folders
RData_dir <- "/data-cbl/mquatre/p2_scrnaseq_renal_landscape/reports/220517_SCTransform_v1_Seurat_v4/data/"
reports_dir <- "/data-cbl/mquatre/p2_scrnaseq_renal_landscape/reports/220517_SCTransform_v1_Seurat_v4/"
data_processed_dir <- "/data-cbl/mquatre/p2_scrnaseq_renal_landscape/reports/220517_SCTransform_v1_Seurat_v4/data/"
figure_dir <- "/data-cbl/mquatre/p2_scrnaseq_renal_landscape/reports/220517_SCTransform_v1_Seurat_v4/figures/"
# Finally set options if needed
options(ggrepel.max.overlaps = Inf)


# --------------------------------------------------
# Figure 2. Healthy human kidney landscape at the single cell level
# --------------------------------------------------
sc_gsm_integrated <- readRDS(paste(reports_dir, "3_sc_gsm_integrated_SeuratAnnot_backup.rds", sep = ""))
# Panels A & B
p1 <- DimPlot(sc_gsm_integrated, reduction = "umap", group.by = "seurat_clusters", label = T, repel = T, pt.size = 0.1) + ggtitle("Louvain clustering")
p1 <- p1 + theme(axis.title.x = element_text("arial", "plain", "grey5", 10)) + theme(axis.title.y = element_text("arial", "plain", "grey5", 10)) + theme(plot.title = element_text("arial", "bold", "grey5", 12)) + theme(legend.text = element_text(size = 10)) + guides(color = guide_legend(override.aes = list(size = 2), ncol = 7)) + theme(legend.position = "bottom")
p1[[1]]$layers[[1]]$aes_params$alpha = 0.8
colors <- c("grey90", "darkseagreen2", "chartreuse3", "cornflowerblue", "cyan3", "cyan4", "grey90", "darkslategray2", "indianred1", "indianred2", "indianred3", "firebrick3", "burlywood", "gold2", "darkorange", "skyblue2", "grey90", "lightpink2", "lightpink3", "lightpink4", "grey90", "yellow2", "olivedrab1", "olivedrab3", "khaki3", "grey90", "plum3", "plum4", "grey90")
p2 <- DimPlot(sc_gsm_integrated, reduction = "umap", group.by = "annot_clusters", label = T, repel = T, pt.size = 0.1, cols = colors) + ggtitle("Kidney cell landscape")
p2 <- p2 + theme(axis.title.x = element_text("arial", "plain", "grey5", 10)) + theme(axis.title.y = element_text("arial", "plain", "grey5", 10)) + theme(plot.title = element_text("arial", "bold", "grey5", 12)) + theme(legend.text = element_text(size = 10)) + guides(color = guide_legend(override.aes = list(size = 2), ncol = 4)) + theme(legend.position = "bottom")
p2[[1]]$layers[[1]]$aes_params$alpha = 0.8
pAB <- ggarrange(p1, p2, labels = c("A", "B"), align = "hv", nrow = 1, ncol = 2, widths = c(1, 1), font.label = list(size = 12, color = "grey5", face = "bold", family = "arial"))
# Panel C
sc_annot_markers <- read.table(paste(reports_dir, "3_sc_gsm_integrated_AnnotClusters_CellTypeMarkers_SCT.csv", sep=""), sep=",")
sc_top_annot_markers <- sc_annot_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
p1 <- DotPlot(sc_gsm_integrated, features = make.unique(sc_top_annot_markers$gene), assay = "SCT", cols = c("green", "red")) + RotatedAxis() + theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10))
pC <- p1 + xlab(NULL) + ylab(NULL) + theme(plot.title = element_text("arial", "bold", "grey5", 11)) + NoLegend()
# Figure
fig <- plot_grid(pAB, pC, labels = c("", "C"), align = "none", nrow = 2, ncol = 1, rel_heights = c(1, 0.7), label_fontfamily = "arial", label_size = 12)
# export | res. w1000 x h1100 px.


# --------------------------------------------------
# Supplementary Figure 1. Evaluation of batch effects correction for scRNA-seq dataset
# --------------------------------------------------
sc_gsm_merge <- readRDS(paste(RData_dir, "2_sc_gsm_integrated_harmony_backup.rds", sep = ""))
# Panel A
p1 <- DimPlot(sc_gsm_merge, reduction = "pca", group.by = "orig.ident", pt.size = 0.1) + ggtitle("Batch effects")
p1 <- p1 + theme(axis.title.x = element_text("arial", "plain", "grey5", 10)) + theme(axis.title.y = element_text("arial", "plain", "grey5", 10)) + theme(plot.title = element_text("arial", "bold", "grey5", 12)) + theme(legend.text = element_text(size = 10)) + guides(color = guide_legend(override.aes = list(size = 2)))
p1[[1]]$layers[[1]]$aes_params$alpha = 0.8
p2 <- DimPlot(sc_gsm_merge, reduction = "harmony", group.by = "orig.ident", pt.size = 0.1) + ggtitle("Harmony correction")
p2 <- p2 + theme(axis.title.x = element_text("arial", "plain", "grey5", 10)) + theme(axis.title.y = element_text("arial", "plain", "grey5", 10)) + theme(plot.title = element_text("arial", "bold", "grey5", 12)) + theme(legend.text = element_text(size = 10)) + guides(color = guide_legend(override.aes = list(size = 2)))
p2[[1]]$layers[[1]]$aes_params$alpha = 0.8
p3 <- DimPlot(sc_gsm_integrated, reduction = "pca", group.by = "orig.ident", pt.size = 0.1) + ggtitle("Seurat correction")
p3 <- p3 + theme(axis.title.x = element_text("arial", "plain", "grey5", 10)) + theme(axis.title.y = element_text("arial", "plain", "grey5", 10)) + theme(plot.title = element_text("arial", "bold", "grey5", 12)) + theme(legend.text = element_text(size = 10)) + guides(color = guide_legend(override.aes = list(size = 2)))
p3[[1]]$layers[[1]]$aes_params$alpha = 0.8
pA <- ggarrange(p1, p2, p3, nrow = 1, ncol = 3, common.legend = TRUE, legend = "bottom")
# Panels B, C, D
p1 <- DimPlot(sc_gsm_integrated, reduction = "umap", pt.size = 0.1, group.by = "orig.ident", label = F) + ggtitle("Sample") + NoLegend()
p1 <- p1 + theme(axis.title.x = element_text("arial", "plain", "grey5", 10)) + theme(axis.title.y = element_text("arial", "plain", "grey5", 10)) + theme(plot.title = element_text("arial", "bold", "grey5", 12)) + theme(legend.text = element_text(size = 0)) + guides(color = guide_legend(override.aes = list(size = 0), ncol = 11))
p1[[1]]$layers[[1]]$aes_params$alpha = 0.8
p2 <- DimPlot(sc_gsm_integrated, reduction = "umap", pt.size = 0.1, group.by = "batch.ident", label = F) + ggtitle("Batch")
p2 <- p2 + theme(axis.title.x = element_text("arial", "plain", "grey5", 10)) + theme(axis.title.y = element_text("arial", "plain", "grey5", 10)) + theme(plot.title = element_text("arial", "bold", "grey5", 12)) + theme(legend.text = element_text(size = 10)) + guides(color = guide_legend(override.aes = list(size = 2), ncol = 2))
p2[[1]]$layers[[1]]$aes_params$alpha = 0.8
sc_gsm_integrated$sex.ident <- factor(x = sc_gsm_integrated$sex.ident, levels = c("M", "F", "?"))
p3 <- DimPlot(sc_gsm_integrated, reduction = "umap", pt.size = 0.1, group.by = "sex.ident", label = F, cols = c("lightcoral", "lightcyan4", "grey90")) + ggtitle("Gender")
p3 <- p3 + theme(axis.title.x = element_text("arial", "plain", "grey5", 10)) + theme(axis.title.y = element_text("arial", "plain", "grey5", 10)) + theme(plot.title = element_text("arial", "bold", "grey5", 12)) + theme(legend.text = element_text(size = 10)) + guides(color = guide_legend(override.aes = list(size = 2)))
p3[[1]]$layers[[1]]$aes_params$alpha = 0.8
pBCD <- ggarrange(p1, p2, p3, nrow = 1, ncol = 3, align = "hv", legend = "bottom", labels = c("B", "C", "D"), font.label = list(size = 12, color = "grey5", face = "bold", family = "arial"))
# Figure
fig <- plot_grid(pA, pBCD, labels = c("A", ""), nrow = 2, ncol = 1, rel_heights = c(1.2, 1), label_fontfamily = "arial", label_size = 12)
# export | res. w1000 x h950 px.


# --------------------------------------------------
# Supplementary Figure 2. Cell type allocation in scRNA-seq dataset
# --------------------------------------------------
# Panel A
cluster_count <- as.data.frame(table(subset(sc_gsm_integrated, seurat_clusters == as.character(0))@meta.data$orig.ident)) # set the structure
names(cluster_count)[1] <- "GSM_ID"
names(cluster_count)[2] <- "ToRemove"
cluster_count$ToRemove <- 0
for (i in levels(sc_gsm_integrated$seurat_clusters)) {
  print(paste("cluster ", i, sep = ""))
  tmp_df <- as.data.frame(table(subset(sc_gsm_integrated, seurat_clusters == as.character(i))@meta.data$orig.ident))
  names(tmp_df)[1] <- "GSM_ID"
  names(tmp_df)[2] <- as.character(paste("Cluster_", i, sep =""))
  cluster_count <- merge(cluster_count, tmp_df, all.x = T, all.y = T)
}
rownames(cluster_count) <- as.character(cluster_count$GSM_ID)
cluster_count <- cluster_count[, 3:length(cluster_count)]
for (i in 1:length(cluster_count)) {
  for (j in 1:length(cluster_count[[i]])) {
    if (is.na(cluster_count[j, i]) == T) {
      cluster_count[j, i] <- 0
    }
  }
}
tmp_df <- data.frame(table(sc_gsm_integrated$orig.ident, sc_gsm_integrated$batch.ident))
lines_to_manage <- c()
for (i in 1:length(rownames(tmp_df))) {
  if (tmp_df$Freq[i] == 0) {
    lines_to_manage <- c(lines_to_manage, i)
  }
} 
tmp_df <- tmp_df[-lines_to_manage,]
rownames(tmp_df) <- tmp_df[[1]]
tmp_df <- tmp_df[, -1]
colnames(tmp_df) <- c("Batch", "Cell.Count")
pA <- as.ggplot(pheatmap(cluster_count, scale = "column", cluster_rows = F, cluster_cols = F, display_numbers = F, color = colorRampPalette(brewer.pal(n = 11, name = "PRGn"))(100), cellheight = 10, cellwidth = 10, annotation_row = tmp_df))
# Panel B
p1 <- VlnPlot(sc_gsm_integrated, group.by = "annot_clusters", assay = "SCT", features = c("S100A8", "S100A9", "CD68", "FCN1", "LILRA5", "CLEC10A", "FCER1A", "MS4A1", "CD79A", "IL7R", "CD3D", "CD8A", "GZMA", "GNLY", "NKG7", "PLVAP", "ENG", "EMCN", "PLAT", "EHD3", "SOX17", "CAV1", "ACTA2", "TAGLN", "PDGFRB", "PLK2", "PLK3", "CFH", "CTGF", "PODXL", "NPHS2", "ALDOB", "APOE", "MIOX", "CRYAB", "VCAM1", "CLDN4", "CLDN10", "SLC12A1", "UMOD", "KNG1", "SLC12A3", "SLC8A1", "CALB1", "AQP2", "AQP3", "FXYD4", "ATP6V1G3", "FOXI1", "SLC4A1", "DMRT2", "SLC26A4", "INSRR"), stack = T, fill.by = "ident", cols = colors)
p1 <- p1 + theme(legend.text = element_text(size = 10)) + guides(color = guide_legend(override.aes = list(size = 2)))
pB <- p1 + xlab(NULL) + ylab(NULL) + theme(plot.title = element_text("arial", "bold", "grey5", 12), axis.title.x = element_text("arial", "plain", "grey5", 10), axis.title.y = element_text("arial", "plain", "grey5", 10))
# Figure
fig <- plot_grid(pA, pB, labels = c("A", "B"), nrow = 2, ncol = 1, rel_heights = c(1, 1), rel_widths = c(1, 1), label_fontfamily = "arial", label_size = 12)
# export | res. w1000 x h1100 px.


# --------------------------------------------------
# Figure 3. Healthy human kidney landscape at the single nucleus level
# --------------------------------------------------
sn_gsm_integrated <- readRDS(paste(reports_dir, "3_sn_gsm_integrated_SeuratAnnot_backup.rds", sep = ""))
# Panels A & B
p1 <- DimPlot(sn_gsm_integrated, reduction = "umap", group.by = "seurat_clusters", label = T, repel = T, pt.size = 0.1) + ggtitle("Louvain clustering")
p1 <- p1 + theme(axis.title.x = element_text("arial", "plain", "grey5", 10)) + theme(axis.title.y = element_text("arial", "plain", "grey5", 10)) + theme(plot.title = element_text("arial", "bold", "grey5", 12)) + theme(legend.text = element_text(size = 10)) + guides(color = guide_legend(override.aes = list(size = 2), ncol = 7)) + theme(legend.position = "bottom")
p1[[1]]$layers[[1]]$aes_params$alpha = 0.8
colors <- c("grey90", "indianred1", "indianred2", "indianred3", "firebrick3", "firebrick1", "burlywood", "gold2", "darkorange", "skyblue2", "grey90", "lightpink2", "lightpink3", "lightpink4", "grey90", "yellow2", "olivedrab1", "olivedrab3", "khaki3", "grey90", "plum3", "plum4")
p2 <- DimPlot(sn_gsm_integrated, reduction = "umap", group.by = "annot_clusters", label = T, repel = T, pt.size = 0.1, cols = colors) + ggtitle("Kidney cell landscape")
p2 <- p2 + theme(axis.title.x = element_text("arial", "plain", "grey5", 10)) + theme(axis.title.y = element_text("arial", "plain", "grey5", 10)) + theme(plot.title = element_text("arial", "bold", "grey5", 12)) + theme(legend.text = element_text(size = 10)) + guides(color = guide_legend(override.aes = list(size = 2), ncol = 4)) + theme(legend.position = "bottom")
p2[[1]]$layers[[1]]$aes_params$alpha = 0.8
pAB <- ggarrange(p1, p2, labels = c("A", "B"), align = "hv", nrow = 1, ncol = 2, widths = c(1, 1), font.label = list(size = 12, color = "grey5", face = "bold", family = "arial"))
# Panel C
sn_annot_markers <- read.table(paste(reports_dir, "3_sn_gsm_integrated_AnnotClusters_CellTypeMarkers_SCT.csv", sep=""), sep=",")
sn_top_annot_markers <- sn_annot_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
p1 <- DotPlot(sn_gsm_integrated, features = make.unique(sn_top_annot_markers$gene), assay = "SCT", cols = c("green", "red")) + RotatedAxis() + theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10))
pC <- p1 + xlab(NULL) + ylab(NULL) + theme(plot.title = element_text("arial", "bold", "grey5", 11)) + NoLegend()
# Figure
fig <- plot_grid(pAB, pC, labels = c("", "C"), align = "none", nrow = 2, ncol = 1, rel_heights = c(1, 0.7), label_fontfamily = "arial", label_size = 12)
# export | res. w1000 x h1100 px.


# --------------------------------------------------
# Supplementary Figure 3. Evaluation of batch effects correction for snRNA-seq dataset
# --------------------------------------------------
sn_gsm_merge <- readRDS(paste(RData_dir, "2_sn_gsm_integrated_harmony_backup.rds", sep = ""))
# Panel A
p1 <- DimPlot(sn_gsm_merge, reduction = "pca", group.by = "orig.ident", pt.size = 0.1) + ggtitle("Batch effects")
p1 <- p1 + theme(axis.title.x = element_text("arial", "plain", "grey5", 10)) + theme(axis.title.y = element_text("arial", "plain", "grey5", 10)) + theme(plot.title = element_text("arial", "bold", "grey5", 12)) + theme(legend.text = element_text(size = 10)) + guides(color = guide_legend(override.aes = list(size = 2)))
p1[[1]]$layers[[1]]$aes_params$alpha = 0.8
p2 <- DimPlot(sn_gsm_merge, reduction = "harmony", group.by = "orig.ident", pt.size = 0.1) + ggtitle("Harmony correction")
p2 <- p2 + theme(axis.title.x = element_text("arial", "plain", "grey5", 10)) + theme(axis.title.y = element_text("arial", "plain", "grey5", 10)) + theme(plot.title = element_text("arial", "bold", "grey5", 12)) + theme(legend.text = element_text(size = 10)) + guides(color = guide_legend(override.aes = list(size = 2)))
p2[[1]]$layers[[1]]$aes_params$alpha = 0.8
p3 <- DimPlot(sn_gsm_integrated, reduction = "pca", group.by = "orig.ident", pt.size = 0.1) + ggtitle("Seurat correction")
p3 <- p3 + theme(axis.title.x = element_text("arial", "plain", "grey5", 10)) + theme(axis.title.y = element_text("arial", "plain", "grey5", 10)) + theme(plot.title = element_text("arial", "bold", "grey5", 12)) + theme(legend.text = element_text(size = 10)) + guides(color = guide_legend(override.aes = list(size = 2)))
p3[[1]]$layers[[1]]$aes_params$alpha = 0.8
pA <- ggarrange(p1, p2, p3, nrow = 1, ncol = 3, common.legend = TRUE, legend = "bottom")
# Panels B, C, D
p1 <- DimPlot(sn_gsm_integrated, reduction = "umap", pt.size = 0.1, group.by = "orig.ident", label = F) + ggtitle("Sample") + NoLegend()
p1 <- p1 + theme(axis.title.x = element_text("arial", "plain", "grey5", 10)) + theme(axis.title.y = element_text("arial", "plain", "grey5", 10)) + theme(plot.title = element_text("arial", "bold", "grey5", 12)) + theme(legend.text = element_text(size = 0)) + guides(color = guide_legend(override.aes = list(size = 0), ncol = 3))
p1[[1]]$layers[[1]]$aes_params$alpha = 0.8
p2 <- DimPlot(sn_gsm_integrated, reduction = "umap", pt.size = 0.1, group.by = "batch.ident", label = F) + ggtitle("Batch")
p2 <- p2 + theme(axis.title.x = element_text("arial", "plain", "grey5", 10)) + theme(axis.title.y = element_text("arial", "plain", "grey5", 10)) + theme(plot.title = element_text("arial", "bold", "grey5", 12)) + theme(legend.text = element_text(size = 10)) + guides(color = guide_legend(override.aes = list(size = 2), ncol = 2))
p2[[1]]$layers[[1]]$aes_params$alpha = 0.8
sn_gsm_integrated$sex.ident <- factor(x = sn_gsm_integrated$sex.ident, levels = c("M", "F", "?"))
p3 <- DimPlot(sn_gsm_integrated, reduction = "umap", pt.size = 0.1, group.by = "sex.ident", label = F, cols = c("lightcoral", "lightcyan4", "grey90")) + ggtitle("Gender")
p3 <- p3 + theme(axis.title.x = element_text("arial", "plain", "grey5", 10)) + theme(axis.title.y = element_text("arial", "plain", "grey5", 10)) + theme(plot.title = element_text("arial", "bold", "grey5", 12)) + theme(legend.text = element_text(size = 10)) + guides(color = guide_legend(override.aes = list(size = 2)))
p3[[1]]$layers[[1]]$aes_params$alpha = 0.8
pBCD <- ggarrange(p1, p2, p3, nrow = 1, ncol = 3, align = "hv", legend = "bottom", labels = c("B", "C", "D"), font.label = list(size = 12, color = "grey5", face = "bold", family = "arial"))
# Figure
fig <- plot_grid(pA, pBCD, labels = c("A", ""), nrow = 2, ncol = 1, rel_heights = c(1, 1), label_fontfamily = "arial", label_size = 12)
# export | res. w1000 x h850 px.


# --------------------------------------------------
# Supplementary Figure 4. Cell type allocation in snRNA-seq dataset
# --------------------------------------------------
# Panel A
cluster_count <- as.data.frame(table(subset(sn_gsm_integrated, seurat_clusters == as.character(0))@meta.data$orig.ident)) # set the structure
names(cluster_count)[1] <- "GSM_ID"
names(cluster_count)[2] <- "ToRemove"
cluster_count$ToRemove <- 0
for (i in levels(sn_gsm_integrated$seurat_clusters)) {
  print(paste("cluster ", i, sep = ""))
  tmp_df <- as.data.frame(table(subset(sn_gsm_integrated, seurat_clusters == as.character(i))@meta.data$orig.ident))
  names(tmp_df)[1] <- "GSM_ID"
  names(tmp_df)[2] <- as.character(paste("Cluster_", i, sep =""))
  cluster_count <- merge(cluster_count, tmp_df, all.x = T, all.y = T)
}
rownames(cluster_count) <- as.character(cluster_count$GSM_ID)
cluster_count <- cluster_count[, 3:length(cluster_count)]
for (i in 1:length(cluster_count)) {
  for (j in 1:length(cluster_count[[i]])) {
    if (is.na(cluster_count[j, i]) == T) {
      cluster_count[j, i] <- 0
    }
  }
}
tmp_df <- data.frame(table(sn_gsm_integrated$orig.ident, sn_gsm_integrated$batch.ident))
lines_to_manage <- c()
for (i in 1:length(rownames(tmp_df))) {
  if (tmp_df$Freq[i] == 0) {
    lines_to_manage <- c(lines_to_manage, i)
  }
} 
tmp_df <- tmp_df[-lines_to_manage,]
rownames(tmp_df) <- tmp_df[[1]]
tmp_df <- tmp_df[, -1]
colnames(tmp_df) <- c("Batch", "Cell.Count")
pA <- as.ggplot(pheatmap(cluster_count, scale = "column", cluster_rows = F, cluster_cols = F, display_numbers = F, color = colorRampPalette(brewer.pal(n = 11, name = "PRGn"))(100), cellheight = 10, cellwidth = 10, annotation_row = tmp_df))
# Panel B
p1 <- VlnPlot(sn_gsm_integrated, group.by = "annot_clusters", assay = "SCT", features = c("PTPRC", "IL7R", "EMCN", "ENG", "PLVAP", "KDR", "EHD3", "CD34", "VEGFC", "ACTA2", "PDGFRB", "ITGA8", "COL12A1", "COL6A2", "CTGF", "CFH", "WT1", "NPHS2", "CUBN", "ALDOB", "MIOX", "GPX3", "CRYAB", "VCAM1", "CLDN4", "CLDN10", "SLC12A1", "UMOD", "KNG1", "SLC12A3", "SLC8A1", "CALB1", "AQP2", "AQP3", "FXYD4", "ATP6V1G3", "FOXI1", "SLC4A1", "DMRT2", "SLC26A4", "INSRR"), stack = T, fill.by = "ident", cols = colors)
p1 <- p1 + theme(legend.text = element_text(size = 10), axis.title.x = element_text("arial", "plain", "grey5", 10), axis.title.y = element_text("arial", "plain", "grey5", 10)) + guides(color = guide_legend(override.aes = list(size = 2)))
pB <- p1 + xlab(NULL) + ylab(NULL) + theme(plot.title = element_text("arial", "bold", "grey5", 11))
# Figure
fig <- plot_grid(pA, pB, labels = c("A", "B"), nrow = 2, ncol = 1, rel_heights = c(1, 2), label_fontfamily = "arial", label_size = 12)
# export | res. w1000 x h700 px.


# --------------------------------------------------
# Figure 4. Joint analysis of scRNA-seq and snRNA-seq datasets
# --------------------------------------------------
scsn_integrated <- readRDS(paste(RData_dir, "5_scsn_integrated_backup3.rds", sep = ""))
scsn_gsm_merge <- readRDS(paste(RData_dir, "5_scsn_gsm_integrated_harmony_backup.rds", sep = ""))
# Panels A et B
p1 <- DimPlot(scsn_integrated, reduction = "umap", pt.size = 0.1, group.by = "batch.ident", label = F, raster = F) + ggtitle("Sample") + theme(legend.position = "right")
p1 <- p1 + theme(axis.title.x = element_text("arial", "plain", "grey5", 10)) + theme(axis.title.y = element_text("arial", "plain", "grey5", 10)) + theme(plot.title = element_text("arial", "bold", "grey5", 12)) + theme(legend.text = element_text(size = 10)) + guides(color = guide_legend(override.aes = list(size = 2), ncol = 1))
p1[[1]]$layers[[1]]$aes_params$alpha = 0.6
p2 <- DimPlot(scsn_integrated, reduction = "umap", pt.size = 0.1, group.by = "techno.ident", label = F, raster = F) + ggtitle("Batch") + theme(legend.position = "right")
p2 <- p2 + theme(axis.title.x = element_text("arial", "plain", "grey5", 10)) + theme(axis.title.y = element_text("arial", "plain", "grey5", 10)) + theme(plot.title = element_text("arial", "bold", "grey5", 12)) + theme(legend.text = element_text(size = 10)) + guides(color = guide_legend(override.aes = list(size = 2), ncol = 1))
p2[[1]]$layers[[1]]$aes_params$alpha = 0.6
pAB <- plot_grid(p1, p2, labels = c("A", "B"), nrow = 1, ncol = 2, rel_widths = c(1, 1), align = "hv", label_fontfamily = "arial", label_size = 12)
# Panel C
colors <- c("grey90", "darkseagreen2", "chartreuse3", "cornflowerblue", "cyan3", "cyan4", "grey90", "darkslategray2", "indianred1", "indianred2", "indianred3", "firebrick3", "firebrick1", "burlywood", "gold2", "darkorange", "skyblue2", "grey90", "lightpink2", "lightpink3", "lightpink4", "grey90", "yellow2", "olivedrab1", "olivedrab3", "khaki3", "grey90", "plum3", "plum4", "grey90")
p1 <- DimPlot(scsn_integrated, group.by = "annot_clusters", split.by = "techno.ident", reduction = "umap", pt.size = 0.1, label = T, repel = T, raster = F, cols = colors) + ggtitle("Kidney single-cell and single-nucleus landscape")
p1[[1]]$layers[[1]]$aes_params$alpha = 0.6
pC <- p1 + theme(axis.title.x = element_text("arial", "plain", "grey5", 10)) + theme(axis.title.y = element_text("arial", "plain", "grey5", 10)) + theme(plot.title = element_text("arial", "bold", "grey5", 12), legend.position = "right") + theme(legend.text = element_text(size = 10)) + guides(color = guide_legend(override.aes = list(size = 2), ncol = 2))
# Figure
fig <- plot_grid(pAB, pC, labels = c("", "C"), nrow = 2, ncol = 1, rel_heights = c(1, 1.25), label_fontfamily = "arial", label_size = 12)
# export | res. w1000 x h800 px.


# --------------------------------------------------
# Supplementary Figure 5. Evaluation of batch effects correction for the integrated single cell and single nucleus dataset
# --------------------------------------------------
# Panel A
p1 <- DimPlot(scsn_gsm_merge, reduction = "pca", group.by = "orig.ident", pt.size = 0.1, raster = F) + ggtitle("Batch effects")
p1 <- p1 + theme(axis.title.x = element_text("arial", "plain", "grey5", 10)) + theme(axis.title.y = element_text("arial", "plain", "grey5", 10)) + theme(plot.title = element_text("arial", "bold", "grey5", 12)) + theme(legend.text = element_text(size = 10)) + guides(color = guide_legend(override.aes = list(size = 2)))
p1[[1]]$layers[[1]]$aes_params$alpha = 0.6
p2 <- DimPlot(scsn_gsm_merge, reduction = "harmony", group.by = "orig.ident", pt.size = 0.1, raster = F) + ggtitle("Harmony correction")
p2 <- p2 + theme(axis.title.x = element_text("arial", "plain", "grey5", 10)) + theme(axis.title.y = element_text("arial", "plain", "grey5", 10)) + theme(plot.title = element_text("arial", "bold", "grey5", 12)) + theme(legend.text = element_text(size = 10)) + guides(color = guide_legend(override.aes = list(size = 2)))
p2[[1]]$layers[[1]]$aes_params$alpha = 0.6
p3 <- DimPlot(scsn_integrated, reduction = "pca", group.by = "orig.ident", pt.size = 0.1, raster = F) + ggtitle("Seurat correction")
p3 <- p3 + theme(axis.title.x = element_text("arial", "plain", "grey5", 10)) + theme(axis.title.y = element_text("arial", "plain", "grey5", 10)) + theme(plot.title = element_text("arial", "bold", "grey5", 12)) + theme(legend.text = element_text(size = 10)) + guides(color = guide_legend(override.aes = list(size = 2)))
p3[[1]]$layers[[1]]$aes_params$alpha = 0.6
# Figure 
fig <- ggarrange(p1, p2, p3, nrow = 1, ncol = 3, common.legend = TRUE, legend = "bottom")
# export | res. w1000 x h500 px.


# --------------------------------------------------
# Figure 5. Enrichments of consensus signatures reveal kidney cell identities within a single-cell RNA-seq dataset
# --------------------------------------------------
zenodo4059315_integrated <- readRDS(paste(RData_dir, "6_sc_zenodo4059315_PredHGT_integrated_backup.rds", sep = ""))
# Panel AB
colors = c("darkseagreen", "darkseagreen2", "limegreen", "chartreuse3", "cornflowerblue", "cyan", "darkslategray2", "indianred1", "indianred2", "indianred3", "lightsalmon1", "lightsalmon3", "grey90", "firebrick3", "darkred", "burlywood", "burlywood3", "darkorange", "mediumorchid", "skyblue2", "grey90", "lightpink2", "lightpink4", "yellow2", "olivedrab1", "khaki3", "plum3", "plum4", "grey90", "forestgreen", "mediumvioletred")
p1 <- DimPlot(zenodo4059315_integrated, reduction = "mcumap", group.by = "annot_clusters", pt.size = 0.1, cols = colors) + ggtitle("Labelling from Kuppe C., et al.")
p1[[1]]$layers[[1]]$aes_params$alpha = 0.6
p1 <- p1 + theme(axis.title.x = element_text("arial", "plain", "grey5", 10)) + theme(axis.title.y = element_text("arial", "plain", "grey5", 10)) + theme(plot.title = element_text("arial", "bold", "grey5", 12)) + theme(legend.text = element_text(size = 10), legend.position = "bottom") + guides(color = guide_legend(override.aes = list(size = 2), ncol = 5))
colors <- c("darkseagreen2", "chartreuse3", "cornflowerblue", "cyan3", "cyan4", "darkslategray2", "indianred1", "indianred2", "indianred3", "firebrick3", "burlywood", "gold2", "darkorange", "skyblue2", "lightpink2", "lightpink3", "lightpink4", "yellow2", "olivedrab1", "olivedrab1", "khaki3", "plum3", "plum4", "grey90")
p2 <- DimPlot(zenodo4059315_integrated, reduction = "mcumap", group.by = "PredHGT.p2", pt.size = 0.1, cols = colors) + ggtitle("Consensus signatures")
p2[[1]]$layers[[1]]$aes_params$alpha = 0.6
p2 <- p2 + theme(axis.title.x = element_text("arial", "plain", "grey5", 10)) + theme(axis.title.y = element_text("arial", "plain", "grey5", 10)) + theme(plot.title = element_text("arial", "bold", "grey5", 12)) + theme(legend.text = element_text(size = 10), legend.position = "bottom") + guides(color = guide_legend(override.aes = list(size = 2), ncol = 5))
pAB <- ggarrange(p1, p2, labels = c("A", "B"), align = "hv", legend = "bottom", common.legend = F, nrow = 1, ncol = 2, widths = c(1, 1), font.label = list(size = 12, color = "grey5", face = "bold", family = "arial"))
# Panel C
comp1 <- as.data.frame(table(zenodo4059315_integrated$annot_clusters))
names(comp1)[1] <- "Cell.type"
names(comp1)[2] <- "Count"
comp1$Annot. <- "Authors.labelling"
comp2 <- as.data.frame(table(zenodo4059315_integrated$PredHGT.p2))
names(comp2)[1] <- "Cell.type"
names(comp2)[2] <- "Count"
comp2$Annot. <- "Prediction"
comp <- merge(comp1, comp2, by = "Cell.type", all.x = T, all.y = T)
for (i in 1:length(comp$Cell.type)) {
  if (is.na(comp$Count.x[i]) == T) {
    comp$Count.x[i] <- 0
    comp$Annot..x <- "Authors.labelling"
  } else if (is.na(comp$Count.y[i]) == T) {
    comp$Count.y[i] <- 0
    comp$Annot..y <- "Prediction"
  }
}
comp1 <- comp[, c(1, 2, 3)]
names(comp1)[2] <- "Count"
names(comp1)[3] <- "Annot."
comp2 <- comp[, c(1, 4, 5)]
names(comp2)[2] <- "Count"
names(comp2)[3] <- "Annot."
comp <- rbind(comp1, comp2)
pC <- ggplot(comp, aes(x = Cell.type, y = Count, fill = Annot.)) + 
  geom_bar(stat = "identity", position = "dodge", color = "grey5") +
  facet_wrap(~ Cell.type, scales = "free") +
  scale_color_manual(values = c("lightsalmon", "lightskyblue")) +
  scale_fill_manual(values = c("lightsalmon", "lightskyblue")) +
  theme(legend.text = element_text(size = 10), axis.title.x = element_text("arial", "plain", "grey5", 0), axis.title.y = element_text("arial", "plain", "grey5", 0)) + 
  theme(strip.background = element_blank(), strip.text = element_text(size = 8, face = "bold")) +
  theme(axis.text.x = element_text(color = "grey5", size = 0)) +
  guides(color = guide_legend(override.aes = list(size = 2)), fill = guide_legend(title = ""))
# Figure
fig <- plot_grid(pAB, pC, labels = c("", "C"), align = "none", nrow = 2, ncol = 1, rel_heights = c(1.25, 1), label_fontfamily = "arial", label_size = 12)
# export | res. w1000 x h1100 px.


# --------------------------------------------------
# Figure 6. Enrichments of consensus signatures reveal kidney cell identities within a single-nucleus RNA-seq dataset
# --------------------------------------------------
gse121862_integrated <- readRDS(paste(RData_dir, "6_sn_gse121862_integrated_PredHGT_backup.rds", sep = ""))
# Panel AB
colors = c("darkseagreen2", "indianred1", "indianred2", "grey90", "firebrick3", "burlywood", "firebrick1", "darkorange", "grey90", "skyblue2", "grey90", "lightpink2", "lightpink3", "lightpink4", "yellow2", "olivedrab1", "khaki3", "grey90", "plum3", "plum4")
p1 <- DimPlot(gse121862_integrated, reduction = "mcumap", group.by = "annot_clusters", pt.size = 0.1, cols = colors) + ggtitle("Labelling from Kuppe C., et al.")
p1[[1]]$layers[[1]]$aes_params$alpha = 0.6
p1 <- p1 + theme(axis.title.x = element_text("arial", "plain", "grey5", 10)) + theme(axis.title.y = element_text("arial", "plain", "grey5", 10)) + theme(plot.title = element_text("arial", "bold", "grey5", 12)) + theme(legend.text = element_text(size = 10), legend.position = "bottom") + guides(color = guide_legend(override.aes = list(size = 2), ncol = 7))
colors <- c("indianred1", "indianred2", "indianred3", "firebrick3", "firebrick1", "burlywood", "gold2", "darkorange", "skyblue2", "lightpink2", "lightpink3", "lightpink4", "yellow2", "olivedrab1", "olivedrab3", "khaki3", "plum3", "plum4", "grey90")
p2 <- DimPlot(gse121862_integrated, reduction = "mcumap", group.by = "PredHGT.p2", pt.size = 0.1, cols = colors) + ggtitle("Consensus signatures")
p2[[1]]$layers[[1]]$aes_params$alpha = 0.6
p2 <- p2 + theme(axis.title.x = element_text("arial", "plain", "grey5", 10)) + theme(axis.title.y = element_text("arial", "plain", "grey5", 10)) + theme(plot.title = element_text("arial", "bold", "grey5", 12)) + theme(legend.text = element_text(size = 10), legend.position = "bottom") + guides(color = guide_legend(override.aes = list(size = 2), ncol = 7))
pAB <- ggarrange(p1, p2, labels = c("A", "B"), align = "hv", legend = "bottom", common.legend = T, nrow = 1, ncol = 2, widths = c(1, 1), font.label = list(size = 12, color = "grey5", face = "bold", family = "arial"))
# Panel C
comp1 <- as.data.frame(table(gse121862_integrated$annot_clusters))
names(comp1)[1] <- "Cell.type"
names(comp1)[2] <- "Count"
comp1$Annot. <- "Authors.labelling"
comp2 <- as.data.frame(table(gse121862_integrated$PredHGT.p2))
names(comp2)[1] <- "Cell.type"
names(comp2)[2] <- "Count"
comp2$Annot. <- "Prediction"
comp <- merge(comp1, comp2, by = "Cell.type", all.x = T, all.y = T)
for (i in 1:length(comp$Cell.type)) {
  if (is.na(comp$Count.x[i]) == T) {
    comp$Count.x[i] <- 0
    comp$Annot..x <- "Authors.labelling"
  } else if (is.na(comp$Count.y[i]) == T) {
    comp$Count.y[i] <- 0
    comp$Annot..y <- "Prediction"
  }
}
comp1 <- comp[, c(1, 2, 3)]
names(comp1)[2] <- "Count"
names(comp1)[3] <- "Annot."
comp2 <- comp[, c(1, 4, 5)]
names(comp2)[2] <- "Count"
names(comp2)[3] <- "Annot."
comp <- rbind(comp1, comp2)
pC <- ggplot(comp, aes(x = Cell.type, y = Count, fill = Annot.)) + 
  geom_bar(stat = "identity", position = "dodge", color = "grey5") +
  facet_wrap(~ Cell.type, scales = "free", ncol = 7) +
  scale_color_manual(values = c("lightsalmon", "lightskyblue")) +
  scale_fill_manual(values = c("lightsalmon", "lightskyblue")) +
  theme(legend.text = element_text(size = 10), axis.title.x = element_text("arial", "plain", "grey5", 0), axis.title.y = element_text("arial", "plain", "grey5", 0)) + 
  theme(strip.background = element_blank(), strip.text = element_text(size = 8, face = "bold")) +
  theme(axis.text.x = element_text(color = "grey5", size = 0)) +
  guides(color = guide_legend(override.aes = list(size = 2)), fill = guide_legend(title = ""))
# Figure
fig <- plot_grid(pAB, pC, labels = c("", "C"), align = "none", nrow = 2, ncol = 1, rel_heights = c(1.4, 1), label_fontfamily = "arial", label_size = 12)
# export | res. w1000 x h850 px.


# --------------------------------------------------
# Supplementary Figure 6. Sample distribution across single-cell and single-nucleus RNA-seq datasets
# --------------------------------------------------
# Panel A
p1 <- DimPlot(zenodo4059315_integrated, reduction = "mcumap", group.by = "Patient.ID", pt.size = 0.1) + ggtitle("Dataset from Kuppe C., et al.")
p1[[1]]$layers[[1]]$aes_params$alpha = 0.6
pA <- p1 + theme(axis.title.x = element_text("arial", "plain", "grey5", 10)) + theme(axis.title.y = element_text("arial", "plain", "grey5", 10)) + theme(plot.title = element_text("arial", "bold", "grey5", 12)) + theme(legend.text = element_text(size = 10), legend.position = "bottom") + guides(color = guide_legend(override.aes = list(size = 2), ncol = 3))
# Panel B
p1 <- DimPlot(gse121862_integrated, reduction = "mcumap", group.by = "orig.ident", pt.size = 0.1) + ggtitle("Dataset from Lake BB., et al.")
p1[[1]]$layers[[1]]$aes_params$alpha = 0.6
pB <- p1 + theme(axis.title.x = element_text("arial", "plain", "grey5", 10)) + theme(axis.title.y = element_text("arial", "plain", "grey5", 10)) + theme(plot.title = element_text("arial", "bold", "grey5", 12)) + theme(legend.text = element_text(size = 10), legend.position = "bottom") + guides(color = guide_legend(override.aes = list(size = 2), ncol = 5))
# Figure
fig <- plot_grid(pA, pB, labels = c("A", "B"), align = "hv", nrow = 1, ncol = 2, rel_widths = c(1, 1), label_fontfamily = "arial", label_size = 12)
# export | res. w1000 x h700 px.


# --------------------------------------------------
# Supplementary Figure 7. Annotation accuracy in scRNA-seq dataset
# --------------------------------------------------
# Panel A
colors = c("darkseagreen", "darkseagreen2", "limegreen", "chartreuse3", "cornflowerblue", "cyan", "darkslategray2", "indianred1", "indianred2", "indianred3", "lightsalmon1", "lightsalmon3", "grey90", "firebrick3", "darkred", "burlywood", "burlywood3", "darkorange", "mediumorchid", "skyblue2", "grey90", "lightpink2", "lightpink4", "yellow2", "olivedrab1", "khaki3", "plum3", "plum4", "grey90", "forestgreen", "mediumvioletred")
p1 <- VlnPlot(zenodo4059315_integrated, group.by = "annot_clusters", assay = "SCT", features = c("S100A8", "S100A9", "CD68", "FCN1", "LILRA5", "CLEC10A", "FCER1A", "MS4A1", "CD79A", "IL7R", "CD3D", "CD8A", "GZMA", "GNLY", "NKG7", "PLVAP", "ENG", "EMCN", "PLAT", "EHD3", "SOX17", "CAV1", "ACTA2", "TAGLN", "PDGFRB", "PLK2", "PLK3", "CFH", "CTGF", "PODXL", "NPHS2", "ALDOB", "APOE", "MIOX", "CRYAB", "VCAM1", "CLDN4", "CLDN10", "SLC12A1", "UMOD", "KNG1", "SLC12A3", "SLC8A1", "CALB1", "AQP2", "AQP3", "FXYD4", "ATP6V1G3", "FOXI1", "SLC4A1", "DMRT2", "SLC26A4", "INSRR", "UPK3A"), stack = T, fill.by = "ident", cols = colors)
p1 <- p1 + theme(legend.text = element_text(size = 10), legend.position = "none", axis.title.x = element_text("arial", "plain", "grey5", 9), axis.title.y = element_text("arial", "plain", "grey5", 9)) + guides(color = guide_legend(override.aes = list(size = 2)))
p1 <- p1 + xlab(NULL) + ylab(NULL) + theme(plot.title = element_text("arial", "bold", "grey5", 11))
# Panel B
colors <- c("darkseagreen2", "chartreuse3", "cornflowerblue", "cyan3", "cyan4", "darkslategray2", "indianred1", "indianred2", "indianred3", "firebrick3", "burlywood", "gold2", "darkorange", "skyblue2", "lightpink2", "lightpink3", "lightpink4", "yellow2", "olivedrab1", "olivedrab1", "khaki3", "plum3", "plum4", "grey90")
p2 <- VlnPlot(zenodo4059315_integrated, group.by = "PredHGT.p2", assay = "SCT", features = c("S100A8", "S100A9", "CD68", "FCN1", "LILRA5", "CLEC10A", "FCER1A", "MS4A1", "CD79A", "IL7R", "CD3D", "CD8A", "GZMA", "GNLY", "NKG7", "PLVAP", "ENG", "EMCN", "PLAT", "EHD3", "SOX17", "CAV1", "ACTA2", "TAGLN", "PDGFRB", "PLK2", "PLK3", "CFH", "CTGF", "PODXL", "NPHS2", "ALDOB", "APOE", "MIOX", "CRYAB", "VCAM1", "CLDN4", "CLDN10", "SLC12A1", "UMOD", "KNG1", "SLC12A3", "SLC8A1", "CALB1", "AQP2", "AQP3", "FXYD4", "ATP6V1G3", "FOXI1", "SLC4A1", "DMRT2", "SLC26A4", "INSRR", "UPK3A"), stack = T, fill.by = "ident", cols = colors)
p2 <- p2 + theme(legend.text = element_text(size = 10), legend.position = "none", axis.title.x = element_text("arial", "plain", "grey5", 9), axis.title.y = element_text("arial", "plain", "grey5", 9)) + guides(color = guide_legend(override.aes = list(size = 2)))
p2 <- p2 + xlab(NULL) + ylab(NULL) + theme(plot.title = element_text("arial", "bold", "grey5", 11))
# Figure 
fig <- plot_grid(p1, p2, labels = c("A", "B"), align = "hv", nrow = 2, ncol = 1, rel_heights = c(1.1, 1), label_fontfamily = "arial", label_size = 12)
# export | res. w1000 x h1100 px.


# --------------------------------------------------
# Supplementary Figure 8. Annotation accuracy in snRNA-seq dataset
# --------------------------------------------------
# Panel A
colors = c("darkseagreen2", "indianred1", "indianred2", "grey90", "firebrick3", "burlywood", "firebrick1", "darkorange", "grey90", "skyblue2", "grey90", "lightpink2", "lightpink3", "lightpink4", "yellow2", "olivedrab1", "khaki3", "grey90", "plum3", "plum4")
p1 <- VlnPlot(gse121862_integrated, group.by = "annot_clusters", assay = "SCT", features = c("PTPRC", "IL7R", "EMCN", "ENG", "PLVAP", "KDR", "EHD3", "CD34", "VEGFC", "ACTA2", "PDGFRB", "ITGA8", "COL12A1", "COL6A2", "CTGF", "CFH", "WT1", "NPHS2", "CUBN", "ALDOB", "MIOX", "GPX3", "CRYAB", "VCAM1", "CLDN4", "CLDN10", "SLC12A1", "UMOD", "KNG1", "SLC12A3", "SLC8A1", "CALB1", "AQP2", "AQP3", "FXYD4", "ATP6V1G3", "FOXI1", "SLC4A1", "DMRT2", "SLC26A4"), stack = T, fill.by = "ident", cols = colors)
p1 <- p1 + theme(legend.text = element_text(size = 10), legend.position = "none", axis.title.x = element_text("arial", "plain", "grey5", 9), axis.title.y = element_text("arial", "plain", "grey5", 9)) + guides(color = guide_legend(override.aes = list(size = 2)))
p1 <- p1 + xlab(NULL) + ylab(NULL) + theme(plot.title = element_text("arial", "bold", "grey5", 11))
# Panel B
colors <- c("indianred1", "indianred2", "indianred3", "firebrick3", "firebrick1", "burlywood", "gold2", "darkorange", "skyblue2", "lightpink2", "lightpink3", "lightpink4", "yellow2", "olivedrab1", "olivedrab3", "khaki3", "plum3", "plum4", "grey90")
p2 <- VlnPlot(gse121862_integrated, group.by = "PredHGT.p2", assay = "SCT", features = c("PTPRC", "IL7R", "EMCN", "ENG", "PLVAP", "KDR", "EHD3", "CD34", "VEGFC", "ACTA2", "PDGFRB", "ITGA8", "COL12A1", "COL6A2", "CTGF", "CFH", "WT1", "NPHS2", "CUBN", "ALDOB", "MIOX", "GPX3", "CRYAB", "VCAM1", "CLDN4", "CLDN10", "SLC12A1", "UMOD", "KNG1", "SLC12A3", "SLC8A1", "CALB1", "AQP2", "AQP3", "FXYD4", "ATP6V1G3", "FOXI1", "SLC4A1", "DMRT2", "SLC26A4"), stack = T, fill.by = "ident", cols = colors)
p2 <- p2 + theme(legend.text = element_text(size = 10), legend.position = "none", axis.title.x = element_text("arial", "plain", "grey5", 9), axis.title.y = element_text("arial", "plain", "grey5", 9)) + guides(color = guide_legend(override.aes = list(size = 2)))
p2 <- p2 + xlab(NULL) + ylab(NULL) + theme(plot.title = element_text("arial", "bold", "grey5", 11))
# Figure 
fig <- plot_grid(p1, p2, labels = c("A", "B"), align = "hv", nrow = 2, ncol = 1, rel_heights = c(1, 1), label_fontfamily = "arial", label_size = 12)
# export | res. w900 x h900 px.



