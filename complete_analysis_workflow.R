# ------------------------------ #
# Complete analysis workflow of kidney scRNA-seq and snRNA-seq datasets
# publication: Quatredeniers M, et al. Sci Data. 2022. doi: XXX
# ------------------------------ #

# ------------------------------ #
# INITIALISATION
# ------------------------------ #
# Load the necessary packages
library(Rcpp)
library(fgsea)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(Seurat)
library(readr)
library(cowplot)
# Set paths to the relevant folders
RData_dir <- "~/wd/meta_analysis/reports/"
reports_dir <- "~/wd/meta_analysis/reports/"
data_processed_dir <- "~/wd/meta_analysis/data/processed/"
figure_dir <- "~/wd/meta_analysis/reports/figures/"
# Set options if needed
options(ggrepel.max.overlaps = Inf)


# ------------------------------ #
# DATA LOADING
# ------------------------------ #
# (sn) Wu H, et al. Comparative Analysis and Refinement of Human PSC-Derived Kidney Organoid Differentiation with Single-Cell Transcriptomics. Cell Stem Cell 2018 Dec 6;23(6):869-881.e8. | PMID 30449713 | GSE118184
gsm3320197 <- read.table("/data-cbl/mquatre/p2_scrnaseq_renal_landscape/data/raw/GSE118184/human.dge.txt", sep = "\t", row.names = 1, header = T)
gsm3320197 <- CreateSeuratObject(gsm3320197, project = "GSM3320197", min.cells = 10, min.features = 200) # and GSM3320198, there are 2 replicates of the same individual in the matrix
gsm3320197@meta.data$batch.ident <- "GSE118184"
gsm3320197@meta.data$techno.ident <- "snRNAseq"
gsm3320197@meta.data$orig.ident <- "GSM3320197-8"
gsm3320197@meta.data$sex.ident <- "M"
# ------------------------------ #
# (sn) Wilson PC, et al. The single-cell transcriptomic landscape of early human diabetic nephropathy. Proc Natl Acad Sci U S A. 2019 Sep 24;116(39):19619-19625. | PMID 31506348 | GSE131882 (GSM3823939/GSM4572192, GSM3823940/GSM4572193, GSM3823941/GSM4572194)
gsm4572192 <- Read10X_h5("/data-cbl/mquatre/p2_scrnaseq_renal_landscape/data/raw/GSE131882/GSM4572192_Control1_filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
gsm4572193 <- Read10X_h5("/data-cbl/mquatre/p2_scrnaseq_renal_landscape/data/raw/GSE131882/GSM4572193_Control2_filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
gsm4572194 <- Read10X_h5("/data-cbl/mquatre/p2_scrnaseq_renal_landscape/data/raw/GSE131882/GSM4572194_Control3_filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
gsm4572192 <- CreateSeuratObject(gsm4572192, project = "GSM4572192", min.cells = 10, min.features = 200)
gsm4572193 <- CreateSeuratObject(gsm4572193, project = "GSM4572193", min.cells = 10, min.features = 200)
gsm4572194 <- CreateSeuratObject(gsm4572194, project = "GSM4572194", min.cells = 10, min.features = 200)
gse131882 <- c(gsm4572192, gsm4572193, gsm4572194)
for (i in 1:length(gse131882)) {
  gse131882[[i]]@meta.data$batch.ident <- "GSE131882"
  gse131882[[i]]@meta.data$techno.ident <- "snRNAseq"
}
gse131882[[1]]@meta.data$sex.ident <- "M"
gse131882[[2]]@meta.data$sex.ident <- "M"
gse131882[[3]]@meta.data$sex.ident <- "F"

# ------------------------------ #
# (sc) Liao J, et al. Single-cell RNA sequencing of human kidney. Sci Data. 2020 Jan 2;7(1):4. | PMID 31896769 | GSE131685 (GSM4145204, GSM4145205, GSM4145206)
gsm4145204 <- Read10X("/data-cbl/mquatre/p2_scrnaseq_renal_landscape/data/raw/GSE131685/GSM4145204/")
gsm4145205 <- Read10X("/data-cbl/mquatre/p2_scrnaseq_renal_landscape/data/raw/GSE131685/GSM4145205/")
gsm4145206 <- Read10X("/data-cbl/mquatre/p2_scrnaseq_renal_landscape/data/raw/GSE131685/GSM4145206/")
gsm4145204 <- CreateSeuratObject(gsm4145204, project = "GSM4145204", min.cells = 10, min.features = 200)
gsm4145205 <- CreateSeuratObject(gsm4145205, project = "GSM4145205", min.cells = 10, min.features = 200)
gsm4145206 <- CreateSeuratObject(gsm4145206, project = "GSM4145206", min.cells = 10, min.features = 200)
gse131685 <- c(gsm4145204, gsm4145205, gsm4145206)
for (i in 1:length(gse131685)) {
  gse131685[[i]]@meta.data$batch.ident <- "GSE131685"
  gse131685[[i]]@meta.data$techno.ident <- "scRNAseq"
}
gse131685[[1]]@meta.data$sex.ident <- "M"
gse131685[[2]]@meta.data$sex.ident <- "F"
gse131685[[3]]@meta.data$sex.ident <- "M"
# ------------------------------ #
# (sn) Wu H, et al. Single-Cell Transcriptomics of a Human Kidney Allograft Biopsy Specimen Defines a Diverse Inflammatory Response. J Am Soc Nephrol. 2018 Aug;29(8):2069-2080. | PMID 29980650 | GSE114156
gsm3135714 <- read.table("/data-cbl/mquatre/p2_scrnaseq_renal_landscape/data/raw/GSE114156/MTS.kidney.dge.txt", sep = "\t", row.names = 1, header = T)
gsm3135714 <- CreateSeuratObject(gsm3135714, project = "GSM3135714", min.cells = 10, min.features = 200)
gsm3135714@meta.data$batch.ident <- "GSE114156"
gsm3135714@meta.data$techno.ident <- "snRNAseq"
gsm3135714@meta.data$sex.ident <- "M"
# ------------------------------ #
# (sn) Muto Y, et al. Single cell transcriptional and chromatin accessibility profiling redefine cellular heterogeneity in the adult human kidney. Nat Commun. 2021 Apr 13;12(1):2190. | PMID 33850129 | GSE151302
gsm4572195 <- Read10X_h5("/data-cbl/mquatre/p2_scrnaseq_renal_landscape/data/raw/GSE151302/GSM4572195_Control4_filtered_feature_bc_matrix.h5")
gsm4572196 <- Read10X_h5("/data-cbl/mquatre/p2_scrnaseq_renal_landscape/data/raw/GSE151302/GSM4572196_Control5_filtered_feature_bc_matrix.h5")
gsm4572195 <- CreateSeuratObject(gsm4572195, project = "GSM4572195", min.cells = 10, min.features = 200)
gsm4572196 <- CreateSeuratObject(gsm4572196, project = "GSM4572196", min.cells = 10, min.features = 200)
gse151302 <- c(gsm4572195, gsm4572196)
for (i in 1:length(gse151302)) {
  gse151302[[i]]@meta.data$batch.ident <- "GSE151302"
  gse151302[[i]]@meta.data$techno.ident <- "snRNAseq"
}
gse151302[[1]]@meta.data$sex.ident <- "M"
gse151302[[2]]@meta.data$sex.ident <- "F"
# ------------------------------ #
# (sc) Zhang Y, et al. Single-cell analyses of renal cell cancers reveal insights into tumor microenvironment, cell of origin, and therapy response. Proc Natl Acad Sci U S A 2021 Jun 15;118(24). | PMID 34099557 | GSE159115 (GSM4819727, GSM4819729, GSM4819730, GSM4819731, GSM4819734, GSM4819736)
gsm4819726 <- Read10X_h5("/data-cbl/mquatre/p2_scrnaseq_renal_landscape/data/raw/GSE159115/GSM4819726_SI_18856_filtered_gene_bc_matrices_h5.h5")
gsm4819728 <- Read10X_h5("/data-cbl/mquatre/p2_scrnaseq_renal_landscape/data/raw/GSE159115/GSM4819728_SI_19704_filtered_gene_bc_matrices_h5.h5")
gsm4819730 <- Read10X_h5("/data-cbl/mquatre/p2_scrnaseq_renal_landscape/data/raw/GSE159115/GSM4819730_SI_21255_filtered_gene_bc_matrices_h5.h5")
gsm4819731 <- Read10X_h5("/data-cbl/mquatre/p2_scrnaseq_renal_landscape/data/raw/GSE159115/GSM4819731_SI_21256_filtered_gene_bc_matrices_h5.h5")
gsm4819733 <- Read10X_h5("/data-cbl/mquatre/p2_scrnaseq_renal_landscape/data/raw/GSE159115/GSM4819733_SI_22369_filtered_gene_bc_matrices_h5.h5")
gsm4819735 <- Read10X_h5("/data-cbl/mquatre/p2_scrnaseq_renal_landscape/data/raw/GSE159115/GSM4819735_SI_22605_filtered_gene_bc_matrices_h5.h5")
gsm4819726 <- CreateSeuratObject(gsm4819726, project = "GSM4819726", min.cells = 10, min.features = 200)
gsm4819728 <- CreateSeuratObject(gsm4819728, project = "GSM4819728", min.cells = 10, min.features = 200)
gsm4819730 <- CreateSeuratObject(gsm4819730, project = "GSM4819730", min.cells = 10, min.features = 200)
gsm4819731 <- CreateSeuratObject(gsm4819731, project = "GSM4819731", min.cells = 10, min.features = 200)
gsm4819730 <- merge(gsm4819730, gsm4819731)
gsm4819730@meta.data$orig.ident <- "GSM4819730-1"
gsm4819733 <- CreateSeuratObject(gsm4819733, project = "GSM4819733", min.cells = 10, min.features = 200)
gsm4819735 <- CreateSeuratObject(gsm4819735, project = "GSM4819735", min.cells = 10, min.features = 200)
gse159115 <- c(gsm4819726, gsm4819728, gsm4819730, gsm4819733, gsm4819735)
for (i in 1:length(gse159115)) {
  gse159115[[i]]@meta.data$batch.ident <- "GSE159115"
  gse159115[[i]]@meta.data$techno.ident <- "scRNAseq"
}
gse159115[[1]]@meta.data$sex.ident <- "M"
gse159115[[2]]@meta.data$sex.ident <- "M"
gse159115[[3]]@meta.data$sex.ident <- "F"
gse159115[[4]]@meta.data$sex.ident <- "M"
gse159115[[5]]@meta.data$sex.ident <- "M"
# ------------------------------ #
# (sc) Menon R, et al. Single cell transcriptomics identifies focal segmental glomerulosclerosis remission endothelial biomarker. JCI Insight. 2020 Mar 26;5(6):e133267. | PMID 32107344 | GSE140989 (24 healthy adult kidney samples)
gsm4191941 <- Read10X("/data-cbl/mquatre/p2_scrnaseq_renal_landscape/data/raw/GSE140989/GSM4191941/")
gsm4191942 <- Read10X("/data-cbl/mquatre/p2_scrnaseq_renal_landscape/data/raw/GSE140989/GSM4191942/")
gsm4191943 <- Read10X("/data-cbl/mquatre/p2_scrnaseq_renal_landscape/data/raw/GSE140989/GSM4191943/")
gsm4191944 <- Read10X("/data-cbl/mquatre/p2_scrnaseq_renal_landscape/data/raw/GSE140989/GSM4191944/")
gsm4191945 <- Read10X("/data-cbl/mquatre/p2_scrnaseq_renal_landscape/data/raw/GSE140989/GSM4191945/")
gsm4191946 <- Read10X("/data-cbl/mquatre/p2_scrnaseq_renal_landscape/data/raw/GSE140989/GSM4191946/")
gsm4191947 <- Read10X("/data-cbl/mquatre/p2_scrnaseq_renal_landscape/data/raw/GSE140989/GSM4191947/")
gsm4191948 <- Read10X("/data-cbl/mquatre/p2_scrnaseq_renal_landscape/data/raw/GSE140989/GSM4191948/")
gsm4191949 <- Read10X("/data-cbl/mquatre/p2_scrnaseq_renal_landscape/data/raw/GSE140989/GSM4191949/")
gsm4191950 <- Read10X("/data-cbl/mquatre/p2_scrnaseq_renal_landscape/data/raw/GSE140989/GSM4191950/")
gsm4191951 <- Read10X("/data-cbl/mquatre/p2_scrnaseq_renal_landscape/data/raw/GSE140989/GSM4191951/")
gsm4191952 <- Read10X("/data-cbl/mquatre/p2_scrnaseq_renal_landscape/data/raw/GSE140989/GSM4191952/")
gsm4191953 <- Read10X("/data-cbl/mquatre/p2_scrnaseq_renal_landscape/data/raw/GSE140989/GSM4191953/")
gsm4191954 <- Read10X("/data-cbl/mquatre/p2_scrnaseq_renal_landscape/data/raw/GSE140989/GSM4191954/")
gsm4191955 <- Read10X("/data-cbl/mquatre/p2_scrnaseq_renal_landscape/data/raw/GSE140989/GSM4191955/")
gsm4191956 <- Read10X("/data-cbl/mquatre/p2_scrnaseq_renal_landscape/data/raw/GSE140989/GSM4191956/")
gsm4191957 <- Read10X("/data-cbl/mquatre/p2_scrnaseq_renal_landscape/data/raw/GSE140989/GSM4191957/")
gsm4191958 <- Read10X("/data-cbl/mquatre/p2_scrnaseq_renal_landscape/data/raw/GSE140989/GSM4191958/")
gsm4191959 <- Read10X("/data-cbl/mquatre/p2_scrnaseq_renal_landscape/data/raw/GSE140989/GSM4191959/")
gsm4191960 <- Read10X("/data-cbl/mquatre/p2_scrnaseq_renal_landscape/data/raw/GSE140989/GSM4191960/")
gsm4191961 <- Read10X("/data-cbl/mquatre/p2_scrnaseq_renal_landscape/data/raw/GSE140989/GSM4191961/")
gsm4191962 <- Read10X("/data-cbl/mquatre/p2_scrnaseq_renal_landscape/data/raw/GSE140989/GSM4191962/")
gsm4191963 <- Read10X("/data-cbl/mquatre/p2_scrnaseq_renal_landscape/data/raw/GSE140989/GSM4191963/")
gsm4191964 <- Read10X("/data-cbl/mquatre/p2_scrnaseq_renal_landscape/data/raw/GSE140989/GSM4191964/")
gsm4191941 <- CreateSeuratObject(gsm4191941, project = "GSM4191941", min.cells = 10, min.features = 200)
gsm4191942 <- CreateSeuratObject(gsm4191942, project = "GSM4191942", min.cells = 10, min.features = 200)
gsm4191943 <- CreateSeuratObject(gsm4191943, project = "GSM4191943", min.cells = 10, min.features = 200)
gsm4191944 <- CreateSeuratObject(gsm4191944, project = "GSM4191944", min.cells = 10, min.features = 200)
gsm4191945 <- CreateSeuratObject(gsm4191945, project = "GSM4191945", min.cells = 10, min.features = 200)
gsm4191946 <- CreateSeuratObject(gsm4191946, project = "GSM4191946", min.cells = 10, min.features = 200)
gsm4191947 <- CreateSeuratObject(gsm4191947, project = "GSM4191947", min.cells = 10, min.features = 200)
gsm4191948 <- CreateSeuratObject(gsm4191948, project = "GSM4191948", min.cells = 10, min.features = 200)
gsm4191949 <- CreateSeuratObject(gsm4191949, project = "GSM4191949", min.cells = 10, min.features = 200)
gsm4191950 <- CreateSeuratObject(gsm4191950, project = "GSM4191950", min.cells = 10, min.features = 200)
gsm4191951 <- CreateSeuratObject(gsm4191951, project = "GSM4191951", min.cells = 10, min.features = 200)
gsm4191952 <- CreateSeuratObject(gsm4191952, project = "GSM4191952", min.cells = 10, min.features = 200)
gsm4191953 <- CreateSeuratObject(gsm4191953, project = "GSM4191953", min.cells = 10, min.features = 200)
gsm4191954 <- CreateSeuratObject(gsm4191954, project = "GSM4191954", min.cells = 10, min.features = 200)
gsm4191955 <- CreateSeuratObject(gsm4191955, project = "GSM4191955", min.cells = 10, min.features = 200)
gsm4191956 <- CreateSeuratObject(gsm4191956, project = "GSM4191956", min.cells = 10, min.features = 200)
gsm4191957 <- CreateSeuratObject(gsm4191957, project = "GSM4191957", min.cells = 10, min.features = 200)
gsm4191958 <- CreateSeuratObject(gsm4191958, project = "GSM4191958", min.cells = 10, min.features = 200)
gsm4191959 <- CreateSeuratObject(gsm4191959, project = "GSM4191959", min.cells = 10, min.features = 200)
gsm4191960 <- CreateSeuratObject(gsm4191960, project = "GSM4191960", min.cells = 10, min.features = 200)
gsm4191961 <- CreateSeuratObject(gsm4191961, project = "GSM4191961", min.cells = 10, min.features = 200)
gsm4191962 <- CreateSeuratObject(gsm4191962, project = "GSM4191962", min.cells = 10, min.features = 200)
gsm4191963 <- CreateSeuratObject(gsm4191963, project = "GSM4191963", min.cells = 10, min.features = 200)
gsm4191964 <- CreateSeuratObject(gsm4191964, project = "GSM4191964", min.cells = 10, min.features = 200)
gse140989 <- c(gsm4191941, gsm4191942, gsm4191943, gsm4191944, gsm4191945, gsm4191946, gsm4191947, gsm4191948, gsm4191949, gsm4191950, gsm4191951, gsm4191952, gsm4191953, gsm4191954, 
               gsm4191955, gsm4191956, gsm4191957, gsm4191958, gsm4191959, gsm4191960, gsm4191961, gsm4191962, gsm4191963, gsm4191964)
for (i in 1:length(gse140989)) {
  gse140989[[i]]@meta.data$batch.ident <- "GSE140989"
  gse140989[[i]]@meta.data$techno.ident <- "scRNAseq"
  gse140989[[i]]@meta.data$sex.ident <- "?"
}
# We need 2 different pipelines to process data from sn- and scRNA-seq studies
# (the overall strategy is to integrate scRNA-seq samples in one hand, integrate snRNA-seq studies in another hand, then integrate the 2 datasets together)
# Create a list of scRNA-seq S4 objects to integrate
sc_gsm_list <- c(gse131685, gse159115, gse140989)
sc_gsm_names <- c("GSM4145204", "GSM4145205", "GSM4145206", "GSM4819726", "GSM4819728", "GSM4819730-1", "GSM4819733", "GSM4819735", "GSM4191941", "GSM4191942", "GSM4191943", "GSM4191944", "GSM4191945", "GSM4191946", "GSM4191947", "GSM4191948", "GSM4191949", "GSM4191950", "GSM4191951", "GSM4191952", "GSM4191953", "GSM4191954", "GSM4191955", "GSM4191956", "GSM4191957", "GSM4191958", "GSM4191959", "GSM4191960", "GSM4191961", "GSM4191962", "GSM4191963", "GSM4191964")
# Then create a list of snRNA-seq S4 objects to integrate
sn_gsm_list <- c(gsm3320197, gse131882, gsm3135714, gse151302)
sn_gsm_names <- c("GSM3823939", "GSM3823940", "GSM3823941", "GSM3320197-8", "GSM3135714", "GSM4572195", "GSM4572196")


# ------------------------------ #
# PRE-PROCESSING
# ------------------------------ #
# A. For scRNA-seq dataset
# ------------------------------ #
print('QC metrics, before any data manipulation :')
for (i in 1:length(sc_gsm_list)) {
  sc_gsm_list[[i]][["percent.mt"]] <- PercentageFeatureSet(sc_gsm_list[[i]], pattern = "^MT-")
  sc_gsm_list[[i]][["percent.rpl"]] <- PercentageFeatureSet(sc_gsm_list[[i]], pattern = c("^RPL"))
  sc_gsm_list[[i]][["percent.rps"]] <- PercentageFeatureSet(sc_gsm_list[[i]], pattern = c("^RPS"))
  sc_gsm_list[[i]][["percent.ribo"]] <- PercentageFeatureSet(sc_gsm_list[[i]], pattern = c("^RPL", "^RPS"))
  print(paste(sc_gsm_names[[i]], ': ',  'nCells = ', as.numeric(table(sc_gsm_list[[i]]$orig.ident)), ' | nFeatures = ', format(round(mean(sc_gsm_list[[i]]$nFeature_RNA), digits = 2), nsmall = 2), ' | nCount = ', format(round(mean(sc_gsm_list[[i]]$nCount_RNA), digits = 2), nsmall = 2), ' | %MT = ', format(round(mean(sc_gsm_list[[i]]$percent.mt), digits = 2), nsmall = 2), sep = ''))
  p1 <- VlnPlot(sc_gsm_list[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rpl", "percent.rps"), pt.size = 0.1, ncol = 5)
  p1 > ggsave(filename = paste(figure_dir, "1_", sc_gsm_names[[i]], "_sc_QC_metrics.png", sep=""), width = 250, height = 150, units = "mm")
}
# Store and save these values in a dataframe
tmp_df <- data.frame("techno.ident" = c(0), "batch.ident" = c(0), "orig.ident" = c(0), "n.cells" = c(0), "n.features" = c(0), "n.count" = c(0), "percent.mt" = c(0), "percent.ribo" = c(0))
for (i in 1:length(sc_gsm_list)) {
  tmp_df[i, "techno.ident"] <- as.character(as.data.frame(table(sc_gsm_list[[i]]$techno.ident))[[1]])
  tmp_df[i, "batch.ident"] <- as.character(as.data.frame(table(sc_gsm_list[[i]]$batch.ident))[[1]])
  tmp_df[i, "orig.ident"] <- sc_gsm_names[[i]]
  tmp_df[i, "n.cells"] <- as.numeric(table(sc_gsm_list[[i]]$orig.ident))
  tmp_df[i, "n.features"] <- format(round(mean(sc_gsm_list[[i]]$nFeature_RNA), digits = 2), nsmall = 2)
  tmp_df[i, "n.count"] <- format(round(mean(sc_gsm_list[[i]]$nCount_RNA), digits = 2), nsmall = 2)
  tmp_df[i, "percent.mt"] <- format(round(mean(sc_gsm_list[[i]]$percent.mt), digits = 2), nsmall = 2)
  tmp_df[i, "percent.ribo"] <- format(round(mean(sc_gsm_list[[i]]$percent.ribo), digits = 2), nsmall = 2)
}
write.table(tmp_df, paste(reports_dir, "1_sc_gsm_list_SampleMetrics_BeforeFiltering.csv", sep=""), sep=",") # export sample characteristics as a table
# Thresholds settings for scRNA-seq samples: percent.mt < 30% 
tmp_df <- data.frame("techno.ident" = c(0), "batch.ident" = c(0), "orig.ident" = c(0), "n.cells" = c(0), "n.features" = c(0), "n.count" = c(0), "percent.mt" = c(0), "percent.ribo" = c(0))
for (i in 1:length(sc_gsm_list)) {
  sc_gsm_list[[i]] <- subset(sc_gsm_list[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < 3500 & percent.mt < 30)
  print(paste(sc_gsm_names[[i]], ': ',  'nCells = ', as.numeric(table(sc_gsm_list[[i]]$orig.ident)), ' | nFeatures = ', format(round(mean(sc_gsm_list[[i]]$nFeature_RNA), digits = 2), nsmall = 2), ' | nCount = ', format(round(mean(sc_gsm_list[[i]]$nCount_RNA), digits = 2), nsmall = 2), ' | %MT = ', format(round(mean(sc_gsm_list[[i]]$percent.mt), digits = 2), nsmall = 2), ' | %RPL = ', format(round(mean(sc_gsm_list[[i]]$percent.rpl), digits = 2), nsmall = 2), ' | %RPS = ', format(round(mean(sc_gsm_list[[i]]$percent.rps), digits = 2), nsmall = 2), sep = ''))
  tmp_df[i, "techno.ident"] <- as.character(as.data.frame(table(sc_gsm_list[[i]]$techno.ident))[[1]])
  tmp_df[i, "batch.ident"] <- as.character(as.data.frame(table(sc_gsm_list[[i]]$batch.ident))[[1]])
  tmp_df[i, "orig.ident"] <- sc_gsm_names[[i]]
  tmp_df[i, "n.cells"] <- as.numeric(table(sc_gsm_list[[i]]$orig.ident))
  tmp_df[i, "n.features"] <- format(round(mean(sc_gsm_list[[i]]$nFeature_RNA), digits = 2), nsmall = 2)
  tmp_df[i, "n.count"] <- format(round(mean(sc_gsm_list[[i]]$nCount_RNA), digits = 2), nsmall = 2)
  tmp_df[i, "percent.mt"] <- format(round(mean(sc_gsm_list[[i]]$percent.mt), digits = 2), nsmall = 2)
  tmp_df[i, "percent.ribo"] <- format(round(mean(sc_gsm_list[[i]]$percent.ribo), digits = 2), nsmall = 2)
} 
write.table(tmp_df, paste(reports_dir, "1_sc_gsm_list_SampleMetrics_AfterFiltering.csv", sep=""), sep=",") # export sample characteristics as a table
sc_gsm_merge <- merge(sc_gsm_list[[1]], c(sc_gsm_list[2:length(sc_gsm_list)]))
saveRDS(sc_gsm_merge, paste(RData_dir, "1_sc_gsm_merge_backup.rds", sep = "")) # save sc_gsm_merge in a .rds file
# sc_gsm_merge <- readRDS(paste(RData_dir, "1_sc_gsm_merge_backup.rds", sep = "")) # load sc_gsm_merge backup
saveRDS(sc_gsm_list, paste(RData_dir, "1_sc_gsm_list_backup.rds", sep = "")) # save sc_gsm_list in a .rds file
# sc_gsm_list <- readRDS(paste(RData_dir, "1_sc_gsm_list_backup.rds", sep = "")) # load sc_gsm_list backup
# ------------------------------ #
# B. For snRNA-seq dataset
# ------------------------------ #
print('QC metrics, before any data manipulation :')
for (i in 1:length(sn_gsm_list)) {
  sn_gsm_list[[i]][["percent.mt"]] <- PercentageFeatureSet(sn_gsm_list[[i]], pattern = "^MT-")
  sn_gsm_list[[i]][["percent.rpl"]] <- PercentageFeatureSet(sn_gsm_list[[i]], pattern = c("^RPL"))
  sn_gsm_list[[i]][["percent.rps"]] <- PercentageFeatureSet(sn_gsm_list[[i]], pattern = c("^RPS"))
  sn_gsm_list[[i]][["percent.ribo"]] <- PercentageFeatureSet(sn_gsm_list[[i]], pattern = c("^RPL", "^RPS"))
  print(paste(sn_gsm_names[[i]], ': ',  'nCells = ', as.numeric(table(sn_gsm_list[[i]]$orig.ident)), ' | nFeatures = ', format(round(mean(sn_gsm_list[[i]]$nFeature_RNA), digits = 2), nsmall = 2), ' | nCount = ', format(round(mean(sn_gsm_list[[i]]$nCount_RNA), digits = 2), nsmall = 2), ' | %MT = ', format(round(mean(sn_gsm_list[[i]]$percent.mt), digits = 2), nsmall = 2), sep = ''))
  p1 <- VlnPlot(sn_gsm_list[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rpl", "percent.rps"), pt.size = 0.1, ncol = 5)
  p1 > ggsave(filename = paste(figure_dir, "1_", sn_gsm_names[[i]], "_sn_QC_metrics.png", sep=""), width = 250, height = 150, units = "mm")
}
# Store and save these values in a dataframe
tmp_df <- data.frame("techno.ident" = c(0), "batch.ident" = c(0), "orig.ident" = c(0), "n.cells" = c(0), "n.features" = c(0), "n.count" = c(0), "percent.mt" = c(0), "percent.ribo" = c(0))
for (i in 1:length(sn_gsm_list)) {
  tmp_df[i, "techno.ident"] <- as.character(as.data.frame(table(sn_gsm_list[[i]]$techno.ident))[[1]])
  tmp_df[i, "batch.ident"] <- as.character(as.data.frame(table(sn_gsm_list[[i]]$batch.ident))[[1]])
  tmp_df[i, "orig.ident"] <- sn_gsm_names[[i]]
  tmp_df[i, "n.cells"] <- as.numeric(table(sn_gsm_list[[i]]$orig.ident))
  tmp_df[i, "n.features"] <- format(round(mean(sn_gsm_list[[i]]$nFeature_RNA), digits = 2), nsmall = 2)
  tmp_df[i, "n.count"] <- format(round(mean(sn_gsm_list[[i]]$nCount_RNA), digits = 2), nsmall = 2)
  tmp_df[i, "percent.mt"] <- format(round(mean(sn_gsm_list[[i]]$percent.mt), digits = 2), nsmall = 2)
  tmp_df[i, "percent.ribo"] <- format(round(mean(sn_gsm_list[[i]]$percent.ribo), digits = 2), nsmall = 2)
}
write.table(tmp_df, paste(reports_dir, "1_sn_gsm_list_SampleMetrics_BeforeFiltering.csv", sep=""), sep=",") # export sample characteristics as a table
# Thresholds settings for snRNA-seq samples: percent.mt < 5% 
tmp_df <- data.frame("techno.ident" = c(0), "batch.ident" = c(0), "orig.ident" = c(0), "n.cells" = c(0), "n.features" = c(0), "n.count" = c(0), "percent.mt" = c(0), "percent.ribo" = c(0))
for (i in 1:length(sn_gsm_list)) {
  sn_gsm_list[[i]] <- subset(sn_gsm_list[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < 3500 & percent.mt < 30)
  print(paste(sn_gsm_names[[i]], ': ',  'nCells = ', as.numeric(table(sn_gsm_list[[i]]$orig.ident)), ' | nFeatures = ', format(round(mean(sn_gsm_list[[i]]$nFeature_RNA), digits = 2), nsmall = 2), ' | nCount = ', format(round(mean(sn_gsm_list[[i]]$nCount_RNA), digits = 2), nsmall = 2), ' | %MT = ', format(round(mean(sn_gsm_list[[i]]$percent.mt), digits = 2), nsmall = 2), ' | %RPL = ', format(round(mean(sn_gsm_list[[i]]$percent.rpl), digits = 2), nsmall = 2), ' | %RPS = ', format(round(mean(sn_gsm_list[[i]]$percent.rps), digits = 2), nsmall = 2), sep = ''))
  tmp_df[i, "techno.ident"] <- as.character(as.data.frame(table(sn_gsm_list[[i]]$techno.ident))[[1]])
  tmp_df[i, "batch.ident"] <- as.character(as.data.frame(table(sn_gsm_list[[i]]$batch.ident))[[1]])
  tmp_df[i, "orig.ident"] <- sn_gsm_names[[i]]
  tmp_df[i, "n.cells"] <- as.numeric(table(sn_gsm_list[[i]]$orig.ident))
  tmp_df[i, "n.features"] <- format(round(mean(sn_gsm_list[[i]]$nFeature_RNA), digits = 2), nsmall = 2)
  tmp_df[i, "n.count"] <- format(round(mean(sn_gsm_list[[i]]$nCount_RNA), digits = 2), nsmall = 2)
  tmp_df[i, "percent.mt"] <- format(round(mean(sn_gsm_list[[i]]$percent.mt), digits = 2), nsmall = 2)
  tmp_df[i, "percent.ribo"] <- format(round(mean(sn_gsm_list[[i]]$percent.ribo), digits = 2), nsmall = 2)
} 
write.table(tmp_df, paste(reports_dir, "1_sn_gsm_list_SampleMetrics_AfterFiltering.csv", sep=""), sep=",") # export sample characteristics as a table
sn_gsm_merge <- merge(sn_gsm_list[[1]], c(sn_gsm_list[2:length(sn_gsm_list)]))
saveRDS(sn_gsm_merge, paste(RData_dir, "1_sn_gsm_merge_backup.rds", sep = "")) # save sn_gsm_merge in a .rds file
# sn_gsm_merge <- readRDS(paste(RData_dir, "1_sn_gsm_merge_backup.rds", sep = "")) # load sn_gsm_merge backup
saveRDS(sn_gsm_list, paste(RData_dir, "1_sn_gsm_list_backup.rds", sep = "")) # save sn_gsm_list in a .rds file
# sn_gsm_list <- readRDS(paste(RData_dir, "1_sn_gsm_list_backup.rds", sep = "")) # load sn_gsm_list backup
# Remove useless variables to free up memory
rm(gsm3320197, gsm4572192, gsm4572193, gsm4572194, gsm4145204, gsm4145205, gsm4145206, gsm3135714, gsm4819726, gsm4819728, gsm4819730, gsm4819731, gsm4819733, gsm4819735, gsm4191941, 
   gsm4191942, gsm4191943, gsm4191944, gsm4191945, gsm4191946, gsm4191947, gsm4191948, gsm4191949, gsm4191950, gsm4191951, gsm4191952, gsm4191953, gsm4191954, gsm4191955, 
   gsm4191956, gsm4191957, gsm4191958, gsm4191959, gsm4191960, gsm4191961, gsm4191962, gsm4191963, gsm4191964, gse131685, gse131882, gse140989, gse159115, gse151302, gsm4572195, gsm4572196, sc_gsm_names, sn_gsm_names, tmp_df)
rm(sn_gsm_merge, sn_gsm_list)


# ------------------------------ #
# II. INTEGRATION USING Seurat v4
# ------------------------------ #
# A. For scRNA-seq dataset
# ------------------------------ #
# First, check the batch effect (thus without integration)
sc_gsm_merge <- SCTransform(sc_gsm_merge, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = T)
DefaultAssay(sc_gsm_merge) <- "SCT"
# sc_gsm_merge <- NormalizeData(sc_gsm_merge, assay = "RNA", verbose = T)
# sc_gsm_merge <- FindVariableFeatures(sc_gsm_merge)
# sc_gsm_merge <- ScaleData(sc_gsm_merge, vars.to.regress = "percent.mt")
sc_gsm_merge <- RunPCA(sc_gsm_merge, nmcs = 50, features = VariableFeatures(sc_gsm_merge), verbose = T)
p1 <- DimPlot(sc_gsm_merge, reduction = "pca", group.by = "orig.ident", pt.size = 0.1)
p2 <- VlnPlot(sc_gsm_merge, features = "PC_1", group.by = "orig.ident", pt.size = 0.1)
plot_grid(p1, p2, ncol = 2, rel_widths = c(1, 1.5)) %>% ggsave(filename = paste(figure_dir, "2_sc_gsm_merge_check_PCA.png", sep=""), width = 500, height = 150, units = "mm")
saveRDS(sc_gsm_merge, paste(RData_dir, "2_sc_gsm_merge_backup.rds", sep = "")) # save sc_gsm_merge in a .rds file
# sc_gsm_merge <- readRDS(paste(RData_dir, "2_sc_gsm_merge_backup.rds", sep = "")) # load sc_gsm_merge backup
rm(sc_gsm_merge)
# ------------------------------ #
# Preparation
for (i in 1:length(sc_gsm_list)) {
  sc_gsm_list[[i]] <- SCTransform(sc_gsm_list[[i]], method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = T)
  DefaultAssay(sc_gsm_list[[i]]) <- "SCT"
  sc_gsm_list[[i]] <- RunPCA(sc_gsm_list[[i]], features = VariableFeatures(sc_gsm_list[[i]]), assay = "SCT", npcs = 50, verbose = T)
}
integration_features <- SelectIntegrationFeatures(object.list = sc_gsm_list, nfeatures = 2500)
sc_gsm_list <- PrepSCTIntegration(object.list = sc_gsm_list, anchor.features = integration_features)
sc_gsm_anchors <- FindIntegrationAnchors(object.list = sc_gsm_list, normalization.method = "SCT", anchor.features = integration_features, reduction = "rpca")
saveRDS(sc_gsm_anchors, paste(RData_dir, "2_sc_gsm_anchors_backup.rds", sep = ""))
saveRDS(sc_gsm_list, paste(RData_dir, "2_sc_gsm_list_backup.rds", sep = "")) 
# sc_gsm_list <- readRDS(paste(RData_dir, "2_sc_gsm_list_backup.rds", sep = "")) 
rm(sc_gsm_list)
sc_gsm_integrated <- IntegrateData(anchorset = sc_gsm_anchors, normalization.method = "SCT", new.assay.name = "seurat.integration")
DefaultAssay(sc_gsm_integrated) <- "seurat.integration"
rm(sc_gsm_anchors)
saveRDS(sc_gsm_integrated, paste(RData_dir, "2_sc_gsm_integrated_backup.rds", sep = "")) 
# sc_gsm_integrated <- readRDS(paste(RData_dir, "2_sc_gsm_integrated_backup.rds", sep = ""))
# Run dimensional reductions
sc_gsm_integrated <- RunPCA(sc_gsm_integrated, features = VariableFeatures(sc_gsm_integrated@assays$seurat.integration), assay = "seurat.integration", npcs = 50, verbose = T)
ElbowPlot(sc_gsm_integrated, ndims = 30, reduction = "pca") %>% ggsave(filename = paste(figure_dir, "2_sc_gsm_integrated_check_PCA_dims_ElbowPlot.png", sep=""), width = 150, height = 150, units = "mm")
VizDimLoadings(sc_gsm_integrated, dims = 1, balanced = T) %>% ggsave(filename = paste(figure_dir, "2_sc_gsm_integrated_check_PC1_associated_genes.png", sep=""), width = 150, height = 150, units = "mm") # Visualize the top genes associated with PC1 (pseudotime)
sc_gsm_integrated <- RunUMAP(sc_gsm_integrated, reduction = "pca", dims = 1:30, verbose = T)
sc_gsm_integrated <- RunTSNE(sc_gsm_integrated, reduction = "pca", dims = 1:30, verbose = T)
p1 <- DimPlot(sc_gsm_integrated, reduction = "umap", group.by = "orig.ident", label = F, pt.size = 0.1) # display the profile of the integrated dataset
p2 <- DimPlot(sc_gsm_integrated, reduction = "umap", group.by = "batch.ident", label = F, pt.size = 0.1) # check also batch integration
plot_grid(p1, p2, ncol = 2, rel_widths = c(1.2, 1)) %>% ggsave(filename = paste(figure_dir, "2_sc_gsm_integrated_batch_vs_sample_ident.png", sep=""), width = 450, height = 150, units = "mm")
p2 <- DimPlot(sc_gsm_integrated, reduction = "umap", group.by = "techno.ident", label = F, pt.size = 0.1) # check also batch integration
plot_grid(p1, p2, ncol = 2, rel_widths = c(1.2, 1)) %>% ggsave(filename = paste(figure_dir, "2_sc_gsm_integrated_techno_vs_sample_ident.png", sep=""), width = 450, height = 150, units = "mm")
# Finally, check the correction of the batch effect (thus after integration)
p1 <- DimPlot(sc_gsm_integrated, reduction = "pca", group.by = "orig.ident", pt.size = 0.1)
p2 <- VlnPlot(sc_gsm_integrated, features = "PC_1", group.by = "orig.ident", pt.size = 0.1)
plot_grid(p1, p2, ncol = 2, rel_widths = c(1, 1.5)) %>% ggsave(filename = paste(figure_dir, "2_sc_gsm_integrated_check_PCA.png", sep=""), width = 500, height = 150, units = "mm")
# Save sc_gsm_integrated before clustering
saveRDS(sc_gsm_integrated, paste(RData_dir, "2_sc_gsm_integrated_seuratv4_backup2.rds", sep = ""))
# sc_gsm_integrated <- readRDS(paste(RData_dir, "2_sc_gsm_integrated_seuratv4_backup2.rds", sep = "")) # load sc_gsm_integrated backup
# ------------------------------ #
# Integrate with Harmony for the sake of comparison
library(harmony)
sc_gsm_merge <- readRDS(paste(RData_dir, "2_sc_gsm_merge_backup.rds", sep = "")) # load sc_gsm_merge backup
# Integrating using Harmony
DefaultAssay(sc_gsm_merge) <- "SCT"
sc_gsm_merge <- RunHarmony(sc_gsm_merge, group.by.vars = "orig.ident", project.dim = F, plot_convergence = T, verbose = T, assay.use = "SCT")
# harmony_embeddings <- Embeddings(sc_gsm_merge, "harmony")
p1 <- DimPlot(sc_gsm_merge, reduction = "harmony", group.by = "orig.ident", pt.size = 0.1)
p2 <- VlnPlot(sc_gsm_merge, features = "harmony_1", group.by = "orig.ident", pt.size = 0.1)
plot_grid(p1, p2, ncol = 2, rel_widths = c(1, 1.5)) %>% ggsave(filename = paste(figure_dir, "2_sc_gsm_integrated_check_harmony.png", sep=""), width = 500, height = 150, units = "mm")
ElbowPlot(sc_gsm_merge, ndims = 30, reduction = "harmony") %>% ggsave(filename = paste(figure_dir, "2_sc_gsm_integrated_check_harmony_dims_ElbowPlot.png", sep=""), width = 150, height = 150, units = "mm")
sc_gsm_merge <- RunUMAP(sc_gsm_merge, reduction = "harmony", umap.method = "umap-learn", dims = 1:20)
p1 <- DimPlot(sc_gsm_merge, reduction = "umap", group.by = "orig.ident", label = F, pt.size = 0.1)
p1[[1]]$layers[[1]]$aes_params$alpha = 0.2 
p1 %>% ggsave(filename = paste(figure_dir, "2_sc_gsm_integrated_check_UMAP_harmony.png", sep=""), width = 150, height = 130, units = "mm") # display the profile of the integrated dataset
# Save pheno.ident, stage.ident and sex.ident
p1 <- DimPlot(sc_gsm_merge, reduction = "umap", group.by = "orig.ident", label = F, pt.size = 0.1)
p1[[1]]$layers[[1]]$aes_params$alpha = 0.2
p1 %>% ggsave(filename = paste(figure_dir, "2_sc_gsm_integrated_check_harmony_OrigIdent.png", sep=""), width = 150, height = 130, units = "mm")
p1 <- DimPlot(sc_gsm_merge, reduction = "umap", group.by = "batch.ident", label = F, pt.size = 0.1)
p1[[1]]$layers[[1]]$aes_params$alpha = 0.2
p1 %>% ggsave(filename = paste(figure_dir, "2_sc_gsm_integrated_check_harmony_BatchIdent.png", sep=""), width = 150, height = 130, units = "mm")
# save sc_gsm_merge before clustering
saveRDS(sc_gsm_merge, paste(RData_dir, "2_sc_gsm_integrated_harmony_backup.rds", sep = ""))
# sc_gsm_merge <- readRDS(paste(RData_dir, "2_sc_gsm_integrated_harmony_backup.rds", sep = "")) # load sc_gsm_merge backup
# ------------------------------ #
# B. For snRNA-seq dataset
# ------------------------------ #
# Reload the datasets
sn_gsm_merge <- readRDS(paste(RData_dir, "1_sn_gsm_merge_backup.rds", sep = "")) # load sn_gsm_merge backup
sn_gsm_list <- readRDS(paste(RData_dir, "1_sn_gsm_list_backup.rds", sep = "")) # load sn_gsm_list backup
# First, check the batch effect (thus without integration)
sn_gsm_merge <- SCTransform(sn_gsm_merge, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = T)
DefaultAssay(sn_gsm_merge) <- "SCT"
# sn_gsm_merge <- NormalizeData(sn_gsm_merge, assay = "RNA", verbose = T)
# sn_gsm_merge <- FindVariableFeatures(sn_gsm_merge)
# sn_gsm_merge <- ScaleData(sn_gsm_merge, vars.to.regress = "percent.mt")
sn_gsm_merge <- RunPCA(sn_gsm_merge, nmcs = 50, features = VariableFeatures(sn_gsm_merge), verbose = T)
p1 <- DimPlot(sn_gsm_merge, reduction = "pca", group.by = "orig.ident", pt.size = 0.1)
p2 <- VlnPlot(sn_gsm_merge, features = "PC_1", group.by = "orig.ident", pt.size = 0.1)
plot_grid(p1, p2, ncol = 2, rel_widths = c(1, 1.5)) %>% ggsave(filename = paste(figure_dir, "2_sn_gsm_merge_check_PCA.png", sep=""), width = 500, height = 150, units = "mm")
saveRDS(sn_gsm_merge, paste(RData_dir, "2_sn_gsm_merge_backup.rds", sep = "")) # save sn_gsm_merge in a .rds file
# sn_gsm_merge <- readRDS(paste(RData_dir, "2_sn_gsm_merge_backup.rds", sep = "")) # load sn_gsm_merge backup
rm(sn_gsm_merge)
# ------------------------------ #
# Preparation
for (i in 1:length(sn_gsm_list)) {
  sn_gsm_list[[i]] <- SCTransform(sn_gsm_list[[i]], method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = T)
  DefaultAssay(sn_gsm_list[[i]]) <- "SCT"
  sn_gsm_list[[i]] <- RunPCA(sn_gsm_list[[i]], features = VariableFeatures(sn_gsm_list[[i]]), assay = "SCT", npcs = 50, verbose = T)
}
integration_features <- SelectIntegrationFeatures(object.list = sn_gsm_list, nfeatures = 2500)
sn_gsm_list <- PrepSCTIntegration(object.list = sn_gsm_list, anchor.features = integration_features)
sn_gsm_anchors <- FindIntegrationAnchors(object.list = sn_gsm_list, normalization.method = "SCT", anchor.features = integration_features, reduction = "rpca")
saveRDS(sn_gsm_anchors, paste(RData_dir, "2_sn_gsm_anchors_backup.rds", sep = ""))
saveRDS(sn_gsm_list, paste(RData_dir, "2_sn_gsm_list_backup.rds", sep = "")) 
# sn_gsm_list <- readRDS(paste(RData_dir, "2_sn_gsm_list_backup.rds", sep = "")) 
rm(sn_gsm_list)
sn_gsm_integrated <- IntegrateData(anchorset = sn_gsm_anchors, normalization.method = "SCT", new.assay.name = "seurat.integration")
DefaultAssay(sn_gsm_integrated) <- "seurat.integration"
rm(sn_gsm_anchors)
saveRDS(sn_gsm_integrated, paste(RData_dir, "2_sn_gsm_integrated_backup.rds", sep = "")) 
# sn_gsm_integrated <- readRDS(paste(RData_dir, "2_sn_gsm_integrated_backup.rds", sep = ""))
# Run dimensional reductions
sn_gsm_integrated <- RunPCA(sn_gsm_integrated, features = VariableFeatures(sn_gsm_integrated@assays$seurat.integration), assay = "seurat.integration", npcs = 50, verbose = T)
ElbowPlot(sn_gsm_integrated, ndims = 30, reduction = "pca") %>% ggsave(filename = paste(figure_dir, "2_sn_gsm_integrated_check_PCA_dims_ElbowPlot.png", sep=""), width = 150, height = 150, units = "mm")
VizDimLoadings(sn_gsm_integrated, dims = 1, balanced = T) %>% ggsave(filename = paste(figure_dir, "2_sn_gsm_integrated_check_PC1_associated_genes.png", sep=""), width = 150, height = 150, units = "mm") # Visualize the top genes associated with PC1 (pseudotime)
sn_gsm_integrated <- RunUMAP(sn_gsm_integrated, reduction = "pca", dims = 1:30, verbose = T)
sn_gsm_integrated <- RunTSNE(sn_gsm_integrated, reduction = "pca", dims = 1:30, verbose = T)
p1 <- DimPlot(sn_gsm_integrated, reduction = "umap", group.by = "orig.ident", label = F, pt.size = 0.1) # display the profile of the integrated dataset
p2 <- DimPlot(sn_gsm_integrated, reduction = "umap", group.by = "batch.ident", label = F, pt.size = 0.1) # check also batch integration
plot_grid(p1, p2, ncol = 2, rel_widths = c(1.2, 1)) %>% ggsave(filename = paste(figure_dir, "2_sn_gsm_integrated_batch_vs_sample_ident.png", sep=""), width = 450, height = 150, units = "mm")
p2 <- DimPlot(sn_gsm_integrated, reduction = "umap", group.by = "techno.ident", label = F, pt.size = 0.1) # check also batch integration
plot_grid(p1, p2, ncol = 2, rel_widths = c(1.2, 1)) %>% ggsave(filename = paste(figure_dir, "2_sn_gsm_integrated_techno_vs_sample_ident.png", sep=""), width = 450, height = 150, units = "mm")
# Finally, check the correction of the batch effect (thus after integration)
p1 <- DimPlot(sn_gsm_integrated, reduction = "pca", group.by = "orig.ident", pt.size = 0.1)
p2 <- VlnPlot(sn_gsm_integrated, features = "PC_1", group.by = "orig.ident", pt.size = 0.1)
plot_grid(p1, p2, ncol = 2, rel_widths = c(1, 1.5)) %>% ggsave(filename = paste(figure_dir, "2_sn_gsm_integrated_check_PCA.png", sep=""), width = 500, height = 150, units = "mm")
# Save sn_gsm_integrated before clustering
saveRDS(sn_gsm_integrated, paste(RData_dir, "2_sn_gsm_integrated_seuratv4_backup2.rds", sep = ""))
# sn_gsm_integrated <- readRDS(paste(RData_dir, "2_sn_gsm_integrated_seuratv4_backup2.rds", sep = "")) # load sn_gsm_integrated backup
# ------------------------------ #
# Integrate with Harmony for the sake of comparison
library(harmony)
sn_gsm_merge <- readRDS(paste(RData_dir, "2_sn_gsm_merge_backup.rds", sep = "")) # load sn_gsm_merge backup
# Integrating using Harmony
DefaultAssay(sn_gsm_merge) <- "SCT"
sn_gsm_merge <- RunHarmony(sn_gsm_merge, group.by.vars = "orig.ident", project.dim = F, plot_convergence = T, verbose = T, assay.use = "SCT")
# harmony_embeddings <- Embeddings(gsm_merge, "harmony")
p1 <- DimPlot(sn_gsm_merge, reduction = "harmony", group.by = "orig.ident", pt.size = 0.1)
p2 <- VlnPlot(sn_gsm_merge, features = "harmony_1", group.by = "orig.ident", pt.size = 0.1)
plot_grid(p1, p2, ncol = 2, rel_widths = c(1, 1.5)) %>% ggsave(filename = paste(figure_dir, "2_sn_gsm_integrated_check_harmony.png", sep=""), width = 500, height = 150, units = "mm")
ElbowPlot(sn_gsm_merge, ndims = 30, reduction = "harmony") %>% ggsave(filename = paste(figure_dir, "2_sn_gsm_integrated_check_harmony_dims_ElbowPlot.png", sep=""), width = 150, height = 150, units = "mm")
sn_gsm_merge <- RunUMAP(sn_gsm_merge, reduction = "harmony", umap.method = "umap-learn", dims = 1:20)
p1 <- DimPlot(sn_gsm_merge, reduction = "umap", group.by = "orig.ident", label = F, pt.size = 0.1)
p1[[1]]$layers[[1]]$aes_params$alpha = 0.2 
p1 %>% ggsave(filename = paste(figure_dir, "2_sn_gsm_integrated_check_UMAP_harmony.png", sep=""), width = 150, height = 130, units = "mm") # display the profile of the integrated dataset
# Save pheno.ident, stage.ident and sex.ident
p1 <- DimPlot(sn_gsm_merge, reduction = "umap", group.by = "orig.ident", label = F, pt.size = 0.1)
p1[[1]]$layers[[1]]$aes_params$alpha = 0.2
p1 %>% ggsave(filename = paste(figure_dir, "2_sn_gsm_integrated_check_harmony_OrigIdent.png", sep=""), width = 150, height = 130, units = "mm")
p1 <- DimPlot(sn_gsm_merge, reduction = "umap", group.by = "batch.ident", label = F, pt.size = 0.1)
p1[[1]]$layers[[1]]$aes_params$alpha = 0.2
p1 %>% ggsave(filename = paste(figure_dir, "2_sn_gsm_integrated_check_harmony_BatchIdent.png", sep=""), width = 150, height = 130, units = "mm")
# save sn_gsm_merge before clustering
saveRDS(sn_gsm_merge, paste(RData_dir, "2_sn_gsm_integrated_harmony_backup.rds", sep = ""))
# sn_gsm_merge <- readRDS(paste(RData_dir, "2_sn_gsm_integrated_harmony_backup.rds", sep = "")) # load sn_gsm_merge backup



# ------------------------------ #
# III. CLUSTERING
# ------------------------------ #
# A. For scRNA-seq dataset
# ------------------------------ #
sc_gsm_integrated <- FindNeighbors(sc_gsm_integrated, reduction = "pca", dims = 1:30)
sc_gsm_integrated <- FindClusters(sc_gsm_integrated, resolution = 3, verbose = T)
p1 <- DimPlot(sc_gsm_integrated, reduction = "umap", group.by = "seurat_clusters", label = T, repel = T, pt.size = 0.1)
p2 <- DimPlot(sc_gsm_integrated, reduction = "umap", group.by = "orig.ident", pt.size = 0.1)
plot_grid(p1, p2, ncol = 2, rel_widths = c(1, 1.15)) %>% ggsave(filename = paste(figure_dir, "3_sc_gsm_integrated_SeuratClusters_res3.png", sep=""), width = 450, height = 150, units = "mm")
saveRDS(sc_gsm_integrated, paste(RData_dir, "3_sc_gsm_integrated_seuratv4_backup.rds", sep = ""))
# sc_gsm_integrated <- readRDS(paste(RData_dir, "3_sc_gsm_integrated_seuratv4_backup.rds", sep = "")) # load sc_gsm_integrated backup
# ------------------------------ #
# CHECK CLUSTER COMPOSITION
# Check how much samples are represented in each of the clusters
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
# Draw a heatmap presenting the representation of each sample in each cluster
library(pheatmap)
library(RColorBrewer)
# Create the annotations needed for the heatmap (which needs dataframe structures for annotations)
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
# Then draw the pheatmap
p1 <- pheatmap(cluster_count, scale = "column", cluster_rows = F, cluster_cols = F, display_numbers = F, color = colorRampPalette(brewer.pal(n = 11, name = "PRGn"))(100), cellheight = 10, cellwidth = 10, annotation_row = tmp_df)
# You should save the heatmap from RStudio | 1000 x 535
write.table(cluster_count, paste(reports_dir, "3_sc_gsm_integrated_SeuratClusters_ClusterCount.csv", sep=""), sep=",") # export "cluster_count" as a table  
# ------------------------------ #
# COMPUTE MARKER GENES
DefaultAssay(sc_gsm_integrated) <- "SCT"
sc_cluster_markers <- FindAllMarkers(sc_gsm_integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(sc_cluster_markers, paste(reports_dir, "3_sc_gsm_integrated_SeuratClusters_CellTypeMarkers_SCT.csv", sep=""), sep=",")
sc_top_cluster_markers <- sc_cluster_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
p1 <- DoHeatmap(subset(sc_gsm_integrated, downsample = 200), features = sc_top_cluster_markers$gene, group.by = "seurat_clusters", group.bar = T, slot = "data") + NoLegend()
p1 %>% ggsave(filename = paste(figure_dir, "3_sc_gsm_integrated_SeuratClusters_Heatmap_SCT.png", sep = ""), width = 500, height = 500, units = "mm")
p1 <- DotPlot(sc_gsm_integrated, features = make.unique(sc_top_cluster_markers$gene), assay = "SCT", cols = c("green", "red")) + RotatedAxis()
p1 %>% ggsave(filename = paste(figure_dir, "3_sc_gsm_integrated_SeuratClusters_DotPlot_SCT.png", sep = ""), width = 800, height = 450, units = "mm")
saveRDS(sc_gsm_integrated, paste(reports_dir, "3_sc_gsm_integrated_backup.rds", sep = ""))
# sc_gsm_integrated <- readRDS(paste(RData_dir, "3_sc_gsm_integrated_backup.rds", sep = "")) # load gsm_merge backup
# ------------------------------ #
# B. For snRNA-seq dataset
# ------------------------------ #
sn_gsm_integrated <- FindNeighbors(sn_gsm_integrated, reduction = "pca", dims = 1:30)
sn_gsm_integrated <- FindClusters(sn_gsm_integrated, resolution = 3, verbose = T)
p1 <- DimPlot(sn_gsm_integrated, reduction = "umap", group.by = "seurat_clusters", label = T, repel = T, pt.size = 0.1)
p2 <- DimPlot(sn_gsm_integrated, reduction = "umap", group.by = "orig.ident", pt.size = 0.1)
plot_grid(p1, p2, ncol = 2, rel_widths = c(1.05, 1)) %>% ggsave(filename = paste(figure_dir, "3_sn_gsm_integrated_SeuratClusters_res3.png", sep=""), width = 450, height = 150, units = "mm")
saveRDS(sn_gsm_integrated, paste(RData_dir, "3_sn_gsm_integrated_seuratv4_backup.rds", sep = ""))
# sn_gsm_integrated <- readRDS(paste(RData_dir, "3_sn_gsm_integrated_seuratv4_backup.rds", sep = "")) # load sn_gsm_integrated backup
# ------------------------------ #
# CHECK CLUSTER COMPOSITION
# Check how much samples are represented in each of the clusters
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
# Draw a heatmap presenting the representation of each sample in each cluster
library(pheatmap)
library(RColorBrewer)
# Create the annotations needed for the heatmap (which needs dataframe structures for annotations)
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
colnames(tmp_df) <- c("Batch", "Nucleus.Count")
# Then draw the pheatmap
p1 <- pheatmap(cluster_count, scale = "column", cluster_rows = F, cluster_cols = F, display_numbers = F, color = colorRampPalette(brewer.pal(n = 11, name = "PRGn"))(100), cellheight = 10, cellwidth = 10, annotation_row = tmp_df)
# You should save the heatmap from RStudio | 1000 x 225
write.table(cluster_count, paste(reports_dir, "3_sn_gsm_integrated_SeuratClusters_ClusterCount.csv", sep=""), sep=",") # export "cluster_count" as a table  
# ------------------------------ #
# COMPUTE MARKER GENES
DefaultAssay(sn_gsm_integrated) <- "SCT"
sn_cluster_markers <- FindAllMarkers(sn_gsm_integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(sn_cluster_markers, paste(reports_dir, "3_sn_gsm_integrated_SeuratClusters_CellTypeMarkers_SCT.csv", sep=""), sep=",")
sn_top_cluster_markers <- sn_cluster_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
p1 <- DoHeatmap(subset(sn_gsm_integrated, downsample = 200), features = sn_top_cluster_markers$gene, group.by = "seurat_clusters", group.bar = T, slot = "data") + NoLegend()
p1 %>% ggsave(filename = paste(figure_dir, "3_sn_gsm_integrated_SeuratClusters_Heatmap_SCT.png", sep = ""), width = 500, height = 500, units = "mm")
p1 <- DotPlot(sn_gsm_integrated, features = make.unique(sn_top_cluster_markers$gene), assay = "SCT", cols = c("green", "red")) + RotatedAxis()
p1 %>% ggsave(filename = paste(figure_dir, "3_sn_gsm_integrated_SeuratClusters_DotPlot_SCT.png", sep = ""), width = 800, height = 450, units = "mm")
saveRDS(sn_gsm_integrated, paste(reports_dir, "3_sn_gsm_integrated_backup.rds", sep = ""))
# sn_gsm_integrated <- readRDS(paste(RData_dir, "3_sn_gsm_integrated_backup.rds", sep = "")) # load gsm_merge backup
rm(sn_gsm_integrated)

# ------------------------------ #
# IV. CELL TYPE ASSIGNMENT
# ------------------------------ #
# A. For scRNA-seq dataset
# ------------------------------ #
sc_gsm_integrated <- readRDS(paste(RData_dir, "3_sc_gsm_integrated_backup.rds", sep = "")) # load sc_gsm_integrated backup
DefaultAssay(sc_gsm_integrated) <- "SCT"
# Immune cells 
VlnPlot(sc_gsm_integrated, group.by = "seurat_clusters", assay = "SCT", features = c("CD3D", "GZMA", "GNLY", "NKG7", "IL7R", "CXCR4", "CD79B", "CD19", "MS4A1", "FCER1A", "CLEC10A", "CD68", "LILRA5", "S100A8", "S100A9"), stack = T, fill.by = "ident")
# Vascular cells
VlnPlot(sc_gsm_integrated, group.by = "seurat_clusters", assay = "SCT", features = c("PLVAP", "EMCN", "ENG", "CAV1", "MYH11", "ACTA2", "PDGFRB", "ITGA8", "SMTN", "PODXL", "NPHS2", "WT1"), stack = T, fill.by = "ident")
# Nephron epithelial cells
VlnPlot(sc_gsm_integrated, group.by = "seurat_clusters", assay = "SCT", features = c("ATP6V1G3", "FOXI1", "SLC26A4", "INSRR", "SLC4A1", "DMRT2", "FXYD4", "AQP2", "AQP3", "CALB1", "SLC8A1", "SLC12A3", "TRPM6", "KNG1", "UMOD", "SLC12A1", "CLCNKA", "TIMP3", "S100A6", "CUBN", "SLC5A12", "SLC22A6", "ALDOB", "GPX3", "LRP2", "ITGB8", "PIGR", "PROM1"), stack = T, fill.by = "ident")
# All cells
p1 <- VlnPlot(sc_gsm_integrated, group.by = "seurat_clusters", assay = "SCT", features = c("CD3D", "GZMA", "GNLY", "NKG7", "IL7R", "CXCR4", "CD79B", "CD19", "FCER1A", "CLEC10A", "CD68", "LILRA5", "S100A8", "S100A9", "PLVAP", "EMCN", "ENG", "CAV1", "MYH11", "ACTA2", "PDGFRB", "ITGA8", "SMTN", "PODXL", "NPHS2", "WT1", "ATP6V1G3", "FOXI1", "SLC26A4", "INSRR", "SLC4A1", "DMRT2", "FXYD4", "AQP2", "AQP3", "CALB1", "SLC8A1", "SLC12A3", "TRPM6", "KNG1", "UMOD", "SLC12A1", "CUBN", "SLC5A12", "SLC22A6", "ALDOB", "GPX3", "LRP2", "ITGB8", "PIGR", "PROM1"), stack = T, fill.by = "ident")
p1 %>% ggsave(filename = paste(figure_dir, "3_sc_gsm_integrated_SeuratClusters_VlnPLot.png", sep=""), width = 550, height = 200, units = "mm")
# Assign cell types to clusters
old_cluster_names <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "39", "40", "41", "42", "43", "44", "45", "46", "47", "48", "49", "50", "51", "52")
new_cluster_names <- c("PTC.na", "PTC.na", "PTC.S2", "PTC.S3", "PTC.S1", "PTC.S2", "PTC.na", "PTC.S3", "PTC.na", "PTC.na", "PTC.S1", "PTC.na", "CD4.T.cells", "EC.vei", "PTC.na", "LoH.TAL", "LoH.DTL", "CD4.T.cells", "NK.cells", "Neutro.", "EC.glom", "LoH.TAL", "IC.A", "CD8.T.cells", "LoH.DTL", "DC", "PTC.na", "vSMC", "DCT", "LoH.na", "EC.art", "LoH.DTL", "CNT", "PC.IMCD", "LoH.na", "Trans.cells", "PTC.na", "LoH.ATL", "PTC.na", "EC.na", "PC.OMCD", "CD4.T.cells", "LoH.DTL", "PC.CCD", "B.cells", "IC.B", "Podo.", "Mono.", "PEC", "IC.A", "Macro.", "PC.na", "EC.na")
sc_gsm_integrated@meta.data$annot_clusters <- plyr::mapvalues(x = sc_gsm_integrated@meta.data$seurat_clusters, from = old_cluster_names, to = new_cluster_names) # replacing cluster names by cell origins
sorted_cell_types <- c("Mono.", "Macro.", "Neutro.", "DC", "B.cells", "CD4.T.cells", "CD8.T.cells", "NK.cells", "EC.vei", "EC.glom", "EC.art", "EC.na", "vSMC", "Podo.", "PEC", "PTC.S1", "PTC.S2", "PTC.S3", "PTC.na", "LoH.DTL", "LoH.ATL", "LoH.TAL", "LoH.na", "DCT", "CNT", "PC.CCD", "PC.OMCD", "PC.IMCD", "PC.na", "IC.A", "IC.B", "Trans.cells")
sc_gsm_integrated$annot_clusters <- factor(x = sc_gsm_integrated$annot_clusters, levels = sorted_cell_types)
p1 <- DimPlot(sc_gsm_integrated, reduction = "umap", pt.size = 0.1, group.by = "annot_clusters", label = T, repel = T)
p1 %>% ggsave(filename = paste(figure_dir, "3_sc_gsm_integrated_AnnotClusters_DimPlot.png", sep=""), width = 200, height = 150, units = "mm")
p1 <- VlnPlot(sc_gsm_integrated, group.by = "annot_clusters", assay = "SCT", features = c("CD3D", "GZMA", "GNLY", "NKG7", "IL7R", "CXCR4", "CD79B", "CD19", "FCER1A", "CLEC10A", "CD68", "LILRA5", "S100A8", "S100A9", "PLVAP", "EMCN", "ENG", "CAV1", "MYH11", "ACTA2", "PDGFRB", "ITGA8", "SMTN", "PODXL", "NPHS2", "WT1", "ATP6V1G3", "FOXI1", "SLC26A4", "INSRR", "SLC4A1", "DMRT2", "FXYD4", "AQP2", "AQP3", "CALB1", "SLC8A1", "SLC12A3", "TRPM6", "KNG1", "UMOD", "SLC12A1", "CUBN", "SLC5A12", "SLC22A6", "ALDOB", "GPX3", "LRP2", "ITGB8", "PIGR", "PROM1"), stack = T, fill.by = "ident")
p1 %>% ggsave(filename = paste(figure_dir, "3_sc_gsm_integrated_AnnotClusters_VlnPLot.png", sep=""), width = 550, height = 200, units = "mm")
# Now with other colors and an alpha parameter set
# colors <- c("brown", "brown1", "coral2", "cornflowerblue", "darkseagreen1", "darkseagreen3", "cyan4", "chocolate1", "firebrick1", "firebrick3", "violetred", "darkmagenta", "goldenrod1", "aquamarine3", "skyblue1", "skyblue3", "grey90", "hotpink3", "hotpink1", "lightpink", "grey90", "grey90", "turquoise3", "yellow3", "olivedrab3", "khaki3", "lemonchiffon4", "palegreen3", "grey90", "plum3", "plum4", "plum1")
colors <- c("brown1", "brown3", "violetred", "darkgoldenrod", "cornflowerblue", "cyan2", "cyan3", "cyan4", "firebrick3", "firebrick2", "firebrick4", "grey90", "coral", "darkorange", "blueviolet", "darkseagreen2", "darkseagreen3", "darkseagreen4", "grey90", "deepskyblue2", "deepskyblue3", "deepskyblue4", "grey90", "darksalmon", "sandybrown", "rosybrown2", "rosybrown3", "rosybrown4", "grey90", "plum3", "plum4", "plum1")
p1 <- DimPlot(sc_gsm_integrated, reduction = "umap", pt.size = 0.1, group.by = "annot_clusters", label = T, repel = T, cols = colors)
p1 %>% ggsave(filename = paste(figure_dir, "3_sc_gsm_integrated_AnnotClusters_DimPlot_Colors_scale.png", sep=""), width = 200, height = 150, units = "mm")
p1[[1]]$layers[[1]]$aes_params$alpha = 0.2
p1 %>% ggsave(filename = paste(figure_dir, "3_sc_gsm_integrated_AnnotClusters_DimPlot_Colors.png", sep=""), width = 200, height = 150, units = "mm")
p1 <- VlnPlot(sc_gsm_integrated, group.by = "annot_clusters", assay = "SCT", features = c("PTPRC", "CD14", "LYZ", "S100A12", "FCN1", "FCGR3A", "CSF3R", "CD68", "FCGR3B", "IL7R", "CLEC10A", "FCER1A", "MS4A1", "CD79A", "CXCR4", "CD3D", "CD69", "CCR7", "GZMA", "GNLY", "CD8A", "TRAC", "NKG7", "EMCN", "ENG", "PLVAP", "EHD3", "PLAT", "TSPAN7", "CAV1", "SOX17", "ACTA2", "PDGFRB", "MYH11", "TAGLN"), stack = T, fill.by = "ident", cols = colors)
p1 %>% ggsave(filename = paste(figure_dir, "3_sc_gsm_integrated_AnnotClusters_VlnPLot_Immune_Others.png", sep=""), width = 550, height = 200, units = "mm")
p1 <- VlnPlot(sc_gsm_integrated, group.by = "annot_clusters", assay = "SCT", features = c("NPHS2", "PODXL", "CLIC5", "WT1", "CTGF", "PTGDS", "ITGB8", "OCIAD2", "APOE", "SLC5A2", "SLC5A12", "SLC22A6", "SLC22A8", "ACSM2A", "MIOX", "AQP1", "CRYAB", "SPP1", "SLC14A2", "CLDN1", "S100A6", "CLDN3", "CLCNKA", "PROX1", "TACSTD2", "UMOD", "SLC12A1", "DEFB1", "KNG1", "SLC8A1", "SLC12A3", "CALB1", "AQP2", "AQP3", "FXYD4", "ATP6V1G3", "FOXI1", "SLC4A1", "DMRT2", "SLC26A4", "INSRR"), stack = T, fill.by = "ident", cols = colors)
p1 %>% ggsave(filename = paste(figure_dir, "3_sc_gsm_integrated_AnnotClusters_VlnPLot_Nephron.png", sep=""), width = 550, height = 200, units = "mm")
p1 <- VlnPlot(sc_gsm_integrated, group.by = "annot_clusters", assay = "SCT", features = c("PTPRC", "CD14", "LYZ", "S100A12", "FCN1", "FCGR3A", "CSF3R", "CD68", "FCGR3B", "IL7R", "CLEC10A", "FCER1A", "MS4A1", "CD79A", "CXCR4", "CD3D", "CD69", "CCR7", "GZMA", "GNLY", "CD8A", "TRAC", "NKG7", "EMCN", "ENG", "PLVAP", "EHD3", "PLAT", "TSPAN7", "CAV1", "SOX17", "ACTA2", "PDGFRB", "MYH11", "TAGLN", "NPHS2", "PODXL", "CLIC5", "WT1", "CTGF", "PTGDS", "ITGB8", "OCIAD2", "APOE", "SLC5A2", "SLC5A12", "SLC22A6", "SLC22A8", "ACSM2A", "MIOX", "AQP1", "CRYAB", "SPP1", "SLC14A2", "CLDN1", "S100A6", "CLDN3", "CLCNKA", "PROX1", "TACSTD2", "UMOD", "SLC12A1", "DEFB1", "KNG1", "SLC8A1", "SLC12A3", "CALB1", "AQP2", "AQP3", "FXYD4", "ATP6V1G3", "FOXI1", "SLC4A1", "DMRT2", "SLC26A4", "INSRR"), stack = T, fill.by = "ident", cols = colors)
p1 %>% ggsave(filename = paste(figure_dir, "3_sc_gsm_integrated_AnnotClusters_VlnPLot_AllSubsets.png", sep=""), width = 600, height = 200, units = "mm")
rm(colors, cols_to_keep, indx, indx005, integration_features, new_cluster_names, old_cluster_names, sorted_cell_types) # free up memory
# Compute consensus cell type signatures
Idents(sc_gsm_integrated) <- "annot_clusters"
sc_annot_markers <- FindAllMarkers(sc_gsm_integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.6)
write.table(sc_annot_markers, paste(reports_dir, "3_sc_gsm_integrated_AnnotClusters_CellTypeMarkers_SCT.csv", sep=""), sep=",")
sc_top_annot_markers <- sc_annot_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
p1 <- DoHeatmap(subset(sc_gsm_integrated, downsample = 200), features = sc_top_annot_markers$gene, group.by = "annot_clusters", group.bar = T, slot = "data") + NoLegend()
p1 %>% ggsave(filename = paste(figure_dir, "3_sc_gsm_integrated_AnnotClusters_Heatmap_SCT.png", sep = ""), width = 500, height = 500, units = "mm")
p1 <- DotPlot(sc_gsm_integrated, features = make.unique(sc_top_annot_markers$gene), assay = "SCT", cols = c("green", "red")) + RotatedAxis()
p1 %>% ggsave(filename = paste(figure_dir, "3_sc_gsm_integrated_AnnotClusters_DotPlot_SCT.png", sep = ""), width = 800, height = 450, units = "mm")
# Save the annotated dataset
saveRDS(sc_gsm_integrated, paste(reports_dir, "3_sc_gsm_integrated_SeuratAnnot_backup.rds", sep = ""))
# sc_gsm_integrated <- readRDS(paste(RData_dir, "3_sc_gsm_integrated_SeuratAnnot_backup.rds", sep = "")) # load sc_gsm_integrated backup


# ------------------------------ #
# B. For snRNA-seq dataset
# ------------------------------ #
sn_gsm_integrated <- readRDS(paste(RData_dir, "3_sn_gsm_integrated_backup.rds", sep = "")) # load sn_gsm_integrated backup
DefaultAssay(sn_gsm_integrated) <- "SCT"
# Immune cells 
VlnPlot(sn_gsm_integrated, group.by = "seurat_clusters", assay = "SCT", features = c("CD3D", "GZMA", "GNLY", "NKG7", "IL7R", "CXCR4", "CD79B", "CD19", "FCER1A", "CLEC10A", "CD68", "LILRA5", "S100A8", "S100A9"), stack = T, fill.by = "ident")
# Vascular cells
VlnPlot(sn_gsm_integrated, group.by = "seurat_clusters", assay = "SCT", features = c("PLVAP", "EMCN", "ENG", "CAV1", "MYH11", "ACTA2", "PDGFRB", "ITGA8", "SMTN", "PODXL", "NPHS2", "WT1"), stack = T, fill.by = "ident")
# Nephron epithelial cells
VlnPlot(sn_gsm_integrated, group.by = "seurat_clusters", assay = "SCT", features = c("ATP6V1G3", "FOXI1", "SLC26A4", "INSRR", "SLC4A1", "DMRT2", "FXYD4", "AQP2", "AQP3", "CALB1", "SLC8A1", "SLC12A3", "TRPM6", "KNG1", "UMOD", "SLC12A1", "CLCNKA", "TIMP3", "S100A6", "CUBN", "SLC5A12", "SLC22A6", "ALDOB", "GPX3", "LRP2", "ITGB8", "PIGR", "PROM1"), stack = T, fill.by = "ident")
# All cells
p1 <- VlnPlot(sn_gsm_integrated, group.by = "seurat_clusters", assay = "SCT", features = c("CD3D", "GZMA", "GNLY", "NKG7", "IL7R", "CXCR4", "CD79B", "CD19", "FCER1A", "CLEC10A", "CD68", "LILRA5", "S100A8", "S100A9", "PLVAP", "EMCN", "ENG", "CAV1", "MYH11", "ACTA2", "PDGFRB", "ITGA8", "SMTN", "PODXL", "NPHS2", "WT1", "ATP6V1G3", "FOXI1", "SLC26A4", "INSRR", "SLC4A1", "DMRT2", "FXYD4", "AQP2", "AQP3", "CALB1", "SLC8A1", "SLC12A3", "TRPM6", "KNG1", "UMOD", "SLC12A1", "CUBN", "SLC5A12", "SLC22A6", "ALDOB", "GPX3", "LRP2", "ITGB8", "PIGR", "PROM1"), stack = T, fill.by = "ident")
p1 %>% ggsave(filename = paste(figure_dir, "3_sn_gsm_integrated_SeuratClusters_VlnPLot.png", sep=""), width = 550, height = 200, units = "mm")
# Assign cell types to clusters
old_cluster_names <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "39", "40", "41", "42", "43", "44", "45", "46", "47", "48", "49", "50")
# new_cluster_names <- c("CNT", "LoH.TAL", "DCT", "PTC.diff", "PTC.S3", "PTC.S1", "IC.A", "PTC.S2", "PTC.diff", "DCT", "PC.CNT", "LoH.TAL", "PTC.S2", "PTC.S1", "LoH.TAL", "PC", "DCT", "PTC.S1", "LoH.TAL", "PTC.S3", "PTC.S1", "PTC.diff", "PTC.S3", "rMSC.glom", "IC.A", "DCT", "LoH.TAL", "LoH.TAL", "PC", "LoH.TAL", "PTC.S2", "LoH.diff", "PTC.S1", "DCT", "EC.vei", "Podo.", "EC", "PC", "EC.glom", "PTC.S2", "IC.B", "Podo.", "Fibro.", "IC.B", "LoH.diff", "EC.art", "vSMC", "Pericytes", "Macro.", "LoH.DTL1", "EC")
# new_cluster_names <- c("CNT", "LoH.TAL", "DCT", "PTC.na", "PTC.na", "PTC.na", "IC.A", "PTC.S3", "PTC.na", "DCT", "CNT", "LoH.TAL", "PTC.S3", "PTC.S2", "LoH.TAL", "PC.IMCD", "DCT", "PTC.S2", "LoH.TAL", "PTC.na", "PTC.S1", "PTC.na", "PTC.na", "EPC", "IC.A", "DCT", "LoH.TAL", "LoH.TAL", "PC.OMCD", "LoH.TAL", "PTC.na", "LoH.ATL", "PTC.S1", "DCT", "EC.vei", "Podo.", "EC", "PC.CCD", "EC.glom", "PTC.na", "IC.B", "Podo.", "Fibro.", "IC.B", "LoH.DTL", "EC.art", "vSMC", "Pericytes", "Macro.", "DCT.na", "EC")
new_cluster_names <- c("CNT", "LoH.TAL", "DCT", "LoH.na", "PTC.na", "PTC.S2", "IC.A", "PTC.S3", "LoH.DTL", "DCT", "CNT", "LoH.TAL", "PTC.S3", "PTC.S2", "LoH.TAL", "PC.CCD", "DCT", "PTC.S1", "LoH.TAL", "PTC.na", "PTC.S1", "LoH.DTL", "PTC.na", "PEC", "IC.A", "DCT", "LoH.TAL", "LoH.TAL", "PC.OMCD", "LoH.TAL", "PTC.na", "LoH.DTL", "PTC.S1", "DCT", "EC.vei", "Podo.", "EC.vei", "PC.IMCD", "EC.glom", "PTC.na", "IC.B", "Podo.", "Fibro.", "IC.B", "LoH.ATL", "EC.art", "vSMC", "Pericytes", "Macro.", "LoH.na", "EC.na")
sn_gsm_integrated@meta.data$annot_clusters <- plyr::mapvalues(x = sn_gsm_integrated@meta.data$seurat_clusters, from = old_cluster_names, to = new_cluster_names) # replacing cluster names by cell types
sorted_cell_types <- c("Macro.", "EC.vei", "EC.glom", "EC.art", "EC.na", "vSMC", "Pericytes", "Fibro.", "Podo.", "PEC", "PTC.S1", "PTC.S2", "PTC.S3", "PTC.na", "LoH.DTL", "LoH.ATL", "LoH.TAL", "LoH.na", "DCT", "CNT", "PC.CCD", "PC.OMCD", "PC.IMCD", "IC.A", "IC.B")
sn_gsm_integrated$annot_clusters <- factor(x = sn_gsm_integrated$annot_clusters, levels = sorted_cell_types)
p1 <- DimPlot(sn_gsm_integrated, reduction = "umap", pt.size = 0.1, group.by = "annot_clusters", label = T, repel = T)
p1 %>% ggsave(filename = paste(figure_dir, "3_sn_gsm_integrated_AnnotClusters_DimPlot.png", sep=""), width = 200, height = 150, units = "mm")
p1 <- VlnPlot(sn_gsm_integrated, group.by = "annot_clusters", assay = "SCT", features = c("CD3D", "GZMA", "GNLY", "NKG7", "IL7R", "CXCR4", "CD79B", "CD19", "FCER1A", "CLEC10A", "CD68", "LILRA5", "S100A8", "S100A9", "PLVAP", "EMCN", "ENG", "CAV1", "MYH11", "ACTA2", "PDGFRB", "ITGA8", "SMTN", "PODXL", "NPHS2", "WT1", "ATP6V1G3", "FOXI1", "SLC26A4", "INSRR", "SLC4A1", "DMRT2", "FXYD4", "AQP2", "AQP3", "CALB1", "SLC8A1", "SLC12A3", "TRPM6", "KNG1", "UMOD", "SLC12A1", "CUBN", "SLC5A12", "SLC22A6", "ALDOB", "GPX3", "LRP2", "ITGB8", "PIGR", "PROM1"), stack = T, fill.by = "ident")
p1 %>% ggsave(filename = paste(figure_dir, "3_sn_gsm_integrated_AnnotClusters_VlnPLot.png", sep=""), width = 550, height = 200, units = "mm")
# Now with other colors and an alpha parameter set
colors <- c("brown3", "firebrick3", "firebrick2", "firebrick4", "grey90", "coral", "coral2", "coral3", "darkorange", "blueviolet", "darkseagreen2", "darkseagreen3", "darkseagreen4", "grey90", "deepskyblue2", "deepskyblue3", "deepskyblue4", "grey90", "darksalmon", "sandybrown", "rosybrown2", "rosybrown3", "rosybrown4", "plum3", "plum4")
p1 <- DimPlot(sn_gsm_integrated, reduction = "umap", pt.size = 0.1, group.by = "annot_clusters", label = T, repel = T, cols = colors)
p1 %>% ggsave(filename = paste(figure_dir, "3_sn_gsm_integrated_AnnotClusters_DimPlot_Colors_scale.png", sep=""), width = 200, height = 150, units = "mm")
p1[[1]]$layers[[1]]$aes_params$alpha = 0.2
p1 %>% ggsave(filename = paste(figure_dir, "3_sn_gsm_integrated_AnnotClusters_DimPlot_Colors.png", sep=""), width = 200, height = 150, units = "mm")
p1 <- VlnPlot(sn_gsm_integrated, group.by = "annot_clusters", assay = "SCT", features = c("ATP6V1G3", "FOXI1", "SLC26A4", "INSRR", "SLC4A1", "DMRT2", "FXYD4", "AQP2", "AQP3", "CALB1", "SLC8A1", "SLC12A3", "TRPM6", "KNG1", "MGP", "CTGF", "UMOD", "SLC12A1", "SLC14A2", "S100A6", "CLCNKA", "CLDN3", "SPP1", "CRYAB", "AQP1", "CUBN", "SLC5A12", "SLC22A6", "ALDOB", "GPX3", "LRP2", "ITGB8", "PIGR", "CD24", "PROM1", "PLVAP", "EMCN", "ENG", "CAV1", "MYH11", "ACTA2", "PDGFRB", "ITGA8", "SMTN", "PODXL", "NPHS2", "WT1", "CD3D", "GZMA", "GNLY", "NKG7", "IL7R", "CXCR4", "CD79B", "CD19", "FCER1A", "CLEC10A", "CD68", "LILRA5", "S100A8", "S100A9"), stack = T, fill.by = "ident", cols = colors)
p1 %>% ggsave(filename = paste(figure_dir, "3_sn_gsm_integrated_AnnotClusters_VlnPLot_Colors.png", sep=""), width = 620, height = 200, units = "mm")
p1 <- VlnPlot(sn_gsm_integrated, group.by = "annot_clusters", assay = "SCT", features = c("PTPRC", "CD14", "CD68", "ITGAM", "PECAM1", "FLT1", "EMCN", "ENG", "PLVAP", "EHD3", "PLAT", "TSPAN7", "CAV1", "ACTA2", "PDGFRB", "MYH11", "TAGLN", "EDNRA"), stack = T, fill.by = "ident", cols = colors)
p1 %>% ggsave(filename = paste(figure_dir, "3_sn_gsm_integrated_AnnotClusters_VlnPLot_Immune_Others.png", sep=""), width = 550, height = 200, units = "mm")
p1 <- VlnPlot(sn_gsm_integrated, group.by = "annot_clusters", assay = "SCT", features = c("NPHS2", "PODXL", "CLIC5", "WT1", "CTGF", "PTGDS", "ITGB8", "CUBN", "GPX3", "SLC5A2", "SLC5A12", "SLC22A6", "SLC22A8", "HNF4A", "MIOX", "CRYAB", "SPP1", "S100A2", "S100A6", "UMOD", "SLC12A1", "DEFB1", "KNG1", "SLC8A1", "SLC12A3", "CALB1", "AQP2", "AQP3", "FXYD4", "ATP6V1G3", "FOXI1", "SLC4A1", "DMRT2", "SLC26A4", "INSRR"), stack = T, fill.by = "ident", cols = colors)
p1 %>% ggsave(filename = paste(figure_dir, "3_sn_gsm_integrated_AnnotClusters_VlnPLot_Nephron.png", sep=""), width = 550, height = 200, units = "mm")
p1 <- VlnPlot(sn_gsm_integrated, group.by = "annot_clusters", assay = "SCT", features = c("PTPRC", "CD14", "CD68", "ITGAM", "PECAM1", "FLT1", "EMCN", "ENG", "PLVAP", "EHD3", "PLAT", "TSPAN7", "CAV1", "ACTA2", "PDGFRB", "MYH11", "TAGLN", "EDNRA", "NPHS2", "PODXL", "CLIC5", "WT1", "CTGF", "PTGDS", "ITGB8", "CUBN", "GPX3", "SLC5A2", "SLC5A12", "SLC22A6", "SLC22A8", "HNF4A", "MIOX", "CRYAB", "SPP1", "S100A2", "S100A6", "UMOD", "SLC12A1", "DEFB1", "KNG1", "SLC8A1", "SLC12A3", "CALB1", "AQP2", "AQP3", "FXYD4", "ATP6V1G3", "FOXI1", "SLC4A1", "DMRT2", "SLC26A4", "INSRR"), stack = T, fill.by = "ident", cols = colors)
p1 %>% ggsave(filename = paste(figure_dir, "3_sn_gsm_integrated_AnnotClusters_VlnPLot_AllSubsets.png", sep=""), width = 500, height = 200, units = "mm")
# Compute consensus cell type signatures
Idents(sn_gsm_integrated) <- "annot_clusters"
sn_annot_markers <- FindAllMarkers(sn_gsm_integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.6)
write.table(sn_annot_markers, paste(reports_dir, "3_sn_gsm_integrated_AnnotClusters_CellTypeMarkers_SCT.csv", sep=""), sep=",")
sn_top_annot_markers <- sn_annot_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
p1 <- DoHeatmap(subset(sn_gsm_integrated, downsample = 200), features = sn_top_annot_markers$gene, group.by = "annot_clusters", group.bar = T, slot = "data") + NoLegend()
p1 %>% ggsave(filename = paste(figure_dir, "3_sn_gsm_integrated_AnnotClusters_Heatmap_SCT.png", sep = ""), width = 500, height = 500, units = "mm")
p1 <- DotPlot(sn_gsm_integrated, features = make.unique(sn_top_annot_markers$gene), assay = "SCT", cols = c("green", "red")) + RotatedAxis()
p1 %>% ggsave(filename = paste(figure_dir, "3_sn_gsm_integrated_AnnotClusters_DotPlot_SCT.png", sep = ""), width = 800, height = 450, units = "mm")
# Save the annotated dataset
saveRDS(sn_gsm_integrated, paste(reports_dir, "3_sn_gsm_integrated_SeuratAnnot_backup.rds", sep = ""))
# sn_gsm_integrated <- readRDS(paste(RData_dir, "3_sn_gsm_integrated_SeuratAnnot_backup.rds", sep = "")) # load sn_gsm_integrated backup


# ------------------------------ #
# V. INTEGRATION OF sc- AND sn-RNAseq
# ------------------------------ #
# First, check the batch effect (thus without integration)
scsn_gsm_merge <- merge(sc_gsm_integrated, sn_gsm_integrated)
rm(sc_gsm_integrated, sn_gsm_integrated) # free up RAM
DefaultAssay(scsn_gsm_merge) <- "RNA"
scsn_gsm_merge <- SCTransform(scsn_gsm_merge, method = "glmGamPoi", vars.to.regress = "percent.mt", assay = "RNA", verbose = T)
DefaultAssay(scsn_gsm_merge) <- "SCT"
# scsn_gsm_merge <- NormalizeData(scsn_gsm_merge, assay = "RNA", verbose = T)
# scsn_gsm_merge <- FindVariableFeatures(sc_gsm_merge)
# scsn_gsm_merge <- ScaleData(sc_gsm_merge, vars.to.regress = "percent.mt")
scsn_gsm_merge <- RunPCA(scsn_gsm_merge, nmcs = 50, features = VariableFeatures(scsn_gsm_merge), verbose = T)
p1 <- DimPlot(scsn_gsm_merge, reduction = "pca", group.by = "orig.ident", pt.size = 0.1, raster = F)
p2 <- VlnPlot(scsn_gsm_merge, features = "PC_1", group.by = "orig.ident", pt.size = 0.1)
plot_grid(p1, p2, ncol = 2, rel_widths = c(1, 1.5)) %>% ggsave(filename = paste(figure_dir, "5_scsn_merge_check_PCA.png", sep=""), width = 500, height = 150, units = "mm")
saveRDS(scsn_gsm_merge, paste(RData_dir, "5_scsn_gsm_merge_backup.rds", sep = "")) # save scsn_gsm_merge in a .rds file
# scsn_gsm_merge <- readRDS(paste(RData_dir, "5_scsn_gsm_merge_backup.rds", sep = "")) # load scsn_gsm_merge backup
# ------------------------------ #
# Integrate with Harmony for the sake of comparison
library(harmony)
# Integrating using Harmony
DefaultAssay(scsn_gsm_merge) <- "SCT"
scsn_gsm_merge <- RunHarmony(scsn_gsm_merge, group.by.vars = "orig.ident", project.dim = F, plot_convergence = T, verbose = T, assay.use = "SCT")
p1 <- DimPlot(scsn_gsm_merge, reduction = "harmony", group.by = "orig.ident", pt.size = 0.1, raster = F)
p2 <- VlnPlot(scsn_gsm_merge, features = "harmony_1", group.by = "orig.ident", pt.size = 0.1)
plot_grid(p1, p2, ncol = 2, rel_widths = c(1, 1.5)) %>% ggsave(filename = paste(figure_dir, "5_scsn_gsm_integrated_check_harmony.png", sep=""), width = 500, height = 150, units = "mm")
ElbowPlot(scsn_gsm_merge, ndims = 30, reduction = "harmony") %>% ggsave(filename = paste(figure_dir, "5_scsn_gsm_integrated_check_harmony_dims_ElbowPlot.png", sep=""), width = 150, height = 150, units = "mm")
scsn_gsm_merge <- RunUMAP(scsn_gsm_merge, reduction = "harmony", umap.method = "umap-learn", dims = 1:20)
p1 <- DimPlot(scsn_gsm_merge, reduction = "umap", group.by = "orig.ident", label = F, pt.size = 0.1, raster = F)
p1[[1]]$layers[[1]]$aes_params$alpha = 0.2 
p1 %>% ggsave(filename = paste(figure_dir, "5_scsn_gsm_integrated_check_UMAP_harmony.png", sep=""), width = 150, height = 130, units = "mm") # display the profile of the integrated dataset
# Save pheno.ident, stage.ident and sex.ident
p1 <- DimPlot(scsn_gsm_merge, reduction = "umap", group.by = "orig.ident", label = F, pt.size = 0.1, raster = F)
p1[[1]]$layers[[1]]$aes_params$alpha = 0.2
p1 %>% ggsave(filename = paste(figure_dir, "5_scsn_gsm_integrated_check_harmony_OrigIdent.png", sep=""), width = 150, height = 130, units = "mm")
p1 <- DimPlot(scsn_gsm_merge, reduction = "umap", group.by = "batch.ident", label = F, pt.size = 0.1, raster = F)
p1[[1]]$layers[[1]]$aes_params$alpha = 0.2
p1 %>% ggsave(filename = paste(figure_dir, "5_scsn_gsm_integrated_check_harmony_BatchIdent.png", sep=""), width = 150, height = 130, units = "mm")
# save scsn_gsm_merge before clustering
saveRDS(scsn_gsm_merge, paste(RData_dir, "5_scsn_gsm_integrated_harmony_backup.rds", sep = ""))
# scsn_gsm_merge <- readRDS(paste(RData_dir, "5_scsn_gsm_integrated_harmony_backup.rds", sep = "")) # load scsn_gsm_merge backup
rm(scsn_gsm_merge)
# Preparation
sc_gsm_integrated <- readRDS(paste(RData_dir, "3_sc_gsm_integrated_SeuratAnnot_backup.rds", sep = ""))
sn_gsm_integrated <- readRDS(paste(RData_dir, "3_sn_gsm_integrated_SeuratAnnot_backup.rds", sep = ""))
scsn_list <- merge(sc_gsm_integrated, sn_gsm_integrated)
rm(sc_gsm_integrated, sn_gsm_integrated) # free up RAM
scsn_list <- SplitObject(scsn_list, split.by = "orig.ident")
for (i in 1:length(scsn_list)) {
  DefaultAssay(scsn_list[[i]]) <- "RNA"
  scsn_list[[i]] <- SCTransform(scsn_list[[i]], method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = T)
  DefaultAssay(scsn_list[[i]]) <- "SCT"
  scsn_list[[i]] <- RunPCA(scsn_list[[i]], features = VariableFeatures(scsn_list[[i]]), assay = "SCT", npcs = 50, verbose = T)
}
integration_features <- SelectIntegrationFeatures(object.list = scsn_list, nfeatures = 2500)
scsn_integrated <- PrepSCTIntegration(object.list = scsn_list, anchor.features = integration_features)
rm(scsn_list)
scsn_anchors <- FindIntegrationAnchors(object.list = scsn_integrated, normalization.method = "SCT", anchor.features = integration_features, reduction = "rpca")
scsn_integrated <- IntegrateData(anchorset = scsn_anchors, normalization.method = "SCT", new.assay.name = "scsn.integration", verbose = T)
DefaultAssay(scsn_integrated) <- "scsn.integration"
saveRDS(scsn_integrated, paste(RData_dir, "5_scsn_integrated_backup.rds", sep = ""))
# scsn_integrated <- readRDS(paste(RData_dir, "5_scsn_integrated_backup.rds", sep = ""))
rm(scsn_anchors)
# Run dimensional reductions
scsn_integrated <- RunPCA(scsn_integrated, features = VariableFeatures(scsn_integrated@assays$scsn.integration), assay = "scsn.integration", npcs = 50, verbose = T)
ElbowPlot(scsn_integrated, ndims = 30, reduction = "pca") %>% ggsave(filename = paste(figure_dir, "5_scsn_integrated_check_PCA_dims_ElbowPlot.png", sep=""), width = 150, height = 150, units = "mm")
VizDimLoadings(scsn_integrated, dims = 1, balanced = T) %>% ggsave(filename = paste(figure_dir, "5_scsn_integrated_check_PC1_associated_genes.png", sep=""), width = 150, height = 150, units = "mm") # Visualize the top genes associated with PC1 (pseudotime)
scsn_integrated <- RunUMAP(scsn_integrated, reduction = "pca", dims = 1:30, verbose = T)
scsn_integrated <- RunTSNE(scsn_integrated, reduction = "pca", dims = 1:30, verbose = T)
p1 <- DimPlot(scsn_integrated, reduction = "umap", group.by = "orig.ident", label = F, pt.size = 0.1, raster = F) # display the profile of the integrated dataset
p2 <- DimPlot(scsn_integrated, reduction = "umap", group.by = "batch.ident", label = F, pt.size = 0.1, raster = F) # check also batch integration
plot_grid(p1, p2, ncol = 2, rel_widths = c(1.2, 1)) %>% ggsave(filename = paste(figure_dir, "5_scsn_integrated_batch_vs_sample_ident.png", sep=""), width = 450, height = 150, units = "mm")
p2 <- DimPlot(scsn_integrated, reduction = "umap", group.by = "techno.ident", label = F, pt.size = 0.1, raster = F) # check also batch integration
plot_grid(p1, p2, ncol = 2, rel_widths = c(1.2, 1)) %>% ggsave(filename = paste(figure_dir, "5_scsn_integrated_techno_vs_sample_ident.png", sep=""), width = 450, height = 150, units = "mm")
# Finally, check the correction of the batch effect (thus after integration)
p1 <- DimPlot(scsn_integrated, reduction = "pca", group.by = "orig.ident", pt.size = 0.1, raster = F)
p2 <- VlnPlot(scsn_integrated, features = "PC_1", group.by = "orig.ident", pt.size = 0.1)
plot_grid(p1, p2, ncol = 2, rel_widths = c(1, 1.5)) %>% ggsave(filename = paste(figure_dir, "5_scsn_integrated_check_PCA.png", sep=""), width = 500, height = 150, units = "mm")
# Save scsn_integrated before clustering
saveRDS(scsn_integrated, paste(RData_dir, "5_scsn_integrated_backup2.rds", sep = ""))
# scsn_integrated <- readRDS(paste(RData_dir, "5_scsn_integrated_backup2.rds", sep = ""))
# Finally, plot "annot_clusters" depending on techno.ident
sorted_cell_types <- c("Mono.", "Macro.", "Neutro.", "DC", "B.cells", "CD4.T.cells", "CD8.T.cells", "NK.cells", "EC.vei", "EC.glom", "EC.art", "EC.na", "vSMC", "Pericytes", "Fibro.", "Podo.", "PEC", "PTC.S1", "PTC.S2", "PTC.S3", "PTC.na", "LoH.DTL", "LoH.ATL", "LoH.TAL", "LoH.na", "DCT", "CNT", "PC.CCD", "PC.OMCD", "PC.IMCD", "PC.na", "IC.A", "IC.B", "Trans.cells")
scsn_integrated$annot_clusters <- factor(x = scsn_integrated$annot_clusters, levels = sorted_cell_types)
colors <- c("brown1", "brown3", "violetred", "darkgoldenrod", "cornflowerblue", "cyan2", "cyan3", "cyan4", "firebrick3", "firebrick2", "firebrick4", "grey90", "coral", "coral2", "coral3", "darkorange", "blueviolet", "darkseagreen2", "darkseagreen3", "darkseagreen4", "grey90", "deepskyblue2", "deepskyblue3", "deepskyblue4", "grey90", "darksalmon", "sandybrown", "rosybrown2", "rosybrown3", "rosybrown4", "grey90", "plum3", "plum4", "plum1")
p1 <- DimPlot(scsn_integrated, group.by = "annot_clusters", reduction = "umap", pt.size = 0.1, label = T, repel = T, raster = F, cols = colors)
p1 %>% ggsave(filename = paste(figure_dir, "5_scsn_integrated_AnnotClusters_DimPlot_colors.png", sep=""), width = 200, height = 150, units = "mm")
p1 <- DimPlot(scsn_integrated, group.by = "annot_clusters", reduction = "umap", pt.size = 0.1, label = T, repel = T, raster = F)
p1[[1]]$layers[[1]]$aes_params$alpha = 0.2
p2 <- DimPlot(scsn_integrated, reduction = "umap", group.by = "batch.ident", label = F, pt.size = 0.1, raster = F)
p2[[1]]$layers[[1]]$aes_params$alpha = 0.2
plot_grid(p1, p2, ncol = 2, rel_widths = c(1.1, 1)) %>% ggsave(filename = paste(figure_dir, "5_scsn_integrated_AnnotClusters_vs_SampleIdent.png", sep=""), width = 400, height = 150, units = "mm")
p2 <- DimPlot(scsn_integrated, group.by = "techno.ident", reduction = "umap", pt.size = 0.1, raster = F)
p2[[1]]$layers[[1]]$aes_params$alpha = 0.2
plot_grid(p1, p2, ncol = 2, rel_widths = c(1.1, 1)) %>% ggsave(filename = paste(figure_dir, "5_scsn_integrated_AnnotClusters_vs_TechnoIdent.png", sep=""), width = 400, height = 150, units = "mm")
p1 <- DimPlot(scsn_integrated, group.by = "annot_clusters", split.by = "techno.ident", reduction = "umap", pt.size = 0.1, label = T, repel = T, raster = F, cols = colors)
p1 %>% ggsave(filename = paste(figure_dir, "5_scsn_integrated_AnnotClusters_SplitByTechnoIdent_scale.png", sep=""), width = 300, height = 150, units = "mm")
p1[[1]]$layers[[1]]$aes_params$alpha = 0.2
p1 %>% ggsave(filename = paste(figure_dir, "5_scsn_integrated_AnnotClusters_SplitByTechnoIdent_colors.png", sep=""), width = 300, height = 150, units = "mm")
saveRDS(scsn_integrated, paste(RData_dir, "5_scsn_integrated_backup3.rds", sep = ""))


# ------------------------------ #
# VI. TESTING THE SIGNATURES IN CELL TYPE ASSIGNMENT
# ------------------------------ #
library(CelliD)
set.seed(1)
# ------------------------------ #
# A. For scRNA-seq dataset
# ------------------------------ #
# I. Using a non-annotated dataset 
# ------------------------------ #
# First, load the query dataset (that must contain metadata with the labels from authors)
gsm4630027 <- Read10X("/data-cbl/mquatre/p2_scrnaseq_renal_landscape/data/raw/GSE152936/GSM4630027_pRCC/")
gsm4630028 <- Read10X("/data-cbl/mquatre/p2_scrnaseq_renal_landscape/data/raw/GSE152936/GSM4630028_ccRCC1/")
gsm4630029 <- Read10X("/data-cbl/mquatre/p2_scrnaseq_renal_landscape/data/raw/GSE152936/GSM4630029_ccRCC2/")
gsm4630030 <- Read10X("/data-cbl/mquatre/p2_scrnaseq_renal_landscape/data/raw/GSE152936/GSM4630030_chRCC/")
gsm4630031 <- Read10X("/data-cbl/mquatre/p2_scrnaseq_renal_landscape/data/raw/GSE152936/GSM4630031_Normal/")
gsm4630027 <- CreateSeuratObject(gsm4630027, project = "gsm4630027", min.cells = 10, min.features = 200)
gsm4630028 <- CreateSeuratObject(gsm4630028, project = "gsm4630028", min.cells = 10, min.features = 200)
gsm4630029 <- CreateSeuratObject(gsm4630029, project = "gsm4630029", min.cells = 10, min.features = 200)
gsm4630030 <- CreateSeuratObject(gsm4630030, project = "gsm4630030", min.cells = 10, min.features = 200)
gsm4630031 <- CreateSeuratObject(gsm4630031, project = "gsm4630031", min.cells = 10, min.features = 200)
# QC filtering
gse152938 <- merge(gsm4630031, y = c(gsm4630027, gsm4630028, gsm4630029, gsm4630030))
gse152938[["percent.mt"]] <- PercentageFeatureSet(gse152938, pattern = "^MT-")
gse152938 <- subset(gse152938, subset = nFeature_RNA > 200 & nFeature_RNA < 3500 & percent.mt < 30)
# Then process the data with normalization, scaling (using SCTransform) and dimension reductions
DefaultAssay(gse152938) <- "RNA"
gse152938 <- SplitObject(gse152938, split.by = "orig.ident")
for (i in 1:length(gse152938)) {
  gse152938[[i]] <- SCTransform(gse152938[[i]], method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = T)
  DefaultAssay(gse152938[[i]]) <- "SCT"
}
integration_features <- SelectIntegrationFeatures(object.list = gse152938, nfeatures = 2500)
gse152938 <- PrepSCTIntegration(object.list = gse152938, anchor.features = integration_features)
gse152938_anchors <- FindIntegrationAnchors(object.list = gse152938, normalization.method = "SCT", anchor.features = integration_features) # reduction = "rpca"
saveRDS(gse152938_anchors, paste(RData_dir, "6_sc_gse152938_anchors_backup.rds", sep = ""))
# gse152938_anchors <- readRDS(paste(RData_dir, "6_sn_gse152938_anchors_backup.rds", sep = ""))
rm(gse152938, gsm4630027, gsm4630028, gsm4630029, gsm4630030, gsm4630031)
gse152938_integrated <- IntegrateData(anchorset = gse152938_anchors, normalization.method = "SCT", new.assay.name = "seurat.integration")
rm(gse152938_anchors)
saveRDS(gse152938_integrated, paste(RData_dir, "6_sc_gse152938_integrated_backup.rds", sep = "")) 
# gse152938_integrated <- readRDS(paste(RData_dir, "6_sc_gse152938_integrated_backup.rds", sep = ""))
gse152938_integrated <- RunMCA(gse152938_integrated, nmcs = 50, features = VariableFeatures(object = gse152938_integrated), assay = "SCT", verbose = T) # assay = "RNA"
# Switch back to the "seurat.integration" assay as default, to compute every next step with this assay
DefaultAssay(gse152938_integrated) <- "seurat.integration"
gse152938_integrated <- RunPCA(gse152938_integrated, features = VariableFeatures(gse152938_integrated@assays$seurat.integration), assay = "seurat.integration", npcs = 50, verbose = T)
gse152938_integrated <- RunUMAP(gse152938_integrated, reduction = "pca", dims = 1:30, verbose = T)
DimPlot(gse152938_integrated, reduction = "umap", group.by = "orig.ident", label = T, repel = T, pt.size = 0.1)
# Save again gse152938_integrated once everything is ok
saveRDS(gse152938_integrated, paste(RData_dir, "6_sc_gse152938_integrated_backup.rds", sep = "")) 
# gse152938_integrated <- readRDS(paste(RData_dir, "6_sc_gse152938_integrated_backup.rds", sep = ""))
# Now try to identify cell types, based on the consensus signatures
library(CelliD)
set.seed(1)
# Setting the signature of cell types as a list; signatures are extracted from Stewart BJ., et al. Science. 2019 --- PMID 31604275
signatures <- as.list(read.table("/data-cbl/mquatre/p2_scrnaseq_renal_landscape/data/external/220302_signatures_p2_scrnaseq.csv", header = T, sep = ";", fill = TRUE))
DefaultAssay(gse152938_integrated) <- "SCT"
# Perform HGT-based enrichment to identify cell types
log10_pval_matrix_is <- RunCellHGT(gse152938_integrated, reduction = "mca", pathways = signatures, dims = 1:50, log.trans = T, n.features = 500, minSize = 10)
# Predict cell type based on signatures
prediction <- rownames(log10_pval_matrix_is)[apply(log10_pval_matrix_is, 2, FUN = which.max)]
# Remove low p-value predictions > 0.01
prediction <- ifelse(apply(log10_pval_matrix_is, 2, FUN = max) > 2, prediction, "unassigned") # this is the step where it's possible to play with p-val threshold, even if it's not adviced
# Store the prediction as metadata of the Seurat object
gse152938_integrated@meta.data$PredHGT.p2 <- prediction
gse152938_integrated@assays[["PredHGT.p2"]] <- CreateAssayObject(data = log10_pval_matrix_is)
# Then plot predicted cell types
p1 <- DimPlot(gse152938_integrated, reduction = "umap", group.by = "PredHGT.p2", label = T, repel = T, pt.size = 0.1)
saveRDS(gse152938_integrated, paste(RData_dir, "6_sc_gse152938_integrated_PredHGT_backup.rds", sep = "")) 


# ------------------------------ #
# I. Using an annotated dataset 
# ------------------------------ #
library(CelliD)
# First, load the query dataset (that must contain metadata with the labels from authors)
# Kuppe C, et al. Decoding myofibroblast origins in human kidney fibrosis. Nature. 2021 Jan;589(7841):281-286. | PMID  33176333 | zenodo 4059315
cd10_neg <- ReadMtx(mtx = "/data-cbl/mquatre/p2_scrnaseq_renal_landscape/data/raw/zenodo4059315/CD10_neg/kidneyMap_UMI_counts.mtx", 
                    features = "/data-cbl/mquatre/p2_scrnaseq_renal_landscape/data/raw/zenodo4059315/CD10_neg/cd10neg_Map_UMI_counts_rowData.txt", 
                    cells = "/data-cbl/mquatre/p2_scrnaseq_renal_landscape/data/raw/zenodo4059315/CD10_neg/cd10neg_Map_UMI_counts_colData.txt",
                    cell.column = 1, skip.cell = 1,
                    feature.column = 1, skip.feature = 1)
metadata <- read.csv("/data-cbl/mquatre/p2_scrnaseq_renal_landscape/data/raw/zenodo4059315/CD10_neg/cd10neg_metadata.csv", sep = ";")
cd10_neg <- CreateSeuratObject(cd10_neg, project = "CD10-", min.cells = 10, min.features = 200, meta.data = metadata)
cd10_pos <- ReadMtx(mtx = "/data-cbl/mquatre/p2_scrnaseq_renal_landscape/data/raw/zenodo4059315/CD10_pos/PTmap_UMI_counts.mtx", 
                    features = "/data-cbl/mquatre/p2_scrnaseq_renal_landscape/data/raw/zenodo4059315/CD10_pos/cd10pos_Map_UMI_counts_rowData.txt", 
                    cells = "/data-cbl/mquatre/p2_scrnaseq_renal_landscape/data/raw/zenodo4059315/CD10_pos/cd10pos_Map_UMI_counts_colData.txt",
                    cell.column = 1, skip.cell = 1,
                    feature.column = 1, skip.feature = 1)
metadata <- read.csv("/data-cbl/mquatre/p2_scrnaseq_renal_landscape/data/raw/zenodo4059315/CD10_pos/cd10pos_metadata.csv", sep = ";")
cd10_pos <- CreateSeuratObject(cd10_pos, project = "CD10+", min.cells = 10, min.features = 200, meta.data = metadata) 
zenodo4059315 <- merge(cd10_neg, cd10_pos) # kidney full dataset is reconstructed
rm(cd10_neg, cd10_pos)
# ------------------------------ #
# PROCESS DATA
# First add some metadata
zenodo4059315[["percent.mt"]] <- PercentageFeatureSet(zenodo4059315, pattern = "^MT-")
zenodo4059315 <- subset(zenodo4059315, subset = nFeature_RNA > 200 & nFeature_RNA < 3500 & percent.mt < 30)
# Then process the data with normalization, scaling (using SCTransform) and dimension reductions
DefaultAssay(zenodo4059315) <- "RNA"
zenodo4059315 <- SplitObject(zenodo4059315, split.by = "Patient.ID")
for (i in 1:length(zenodo4059315)) {
  zenodo4059315[[i]] <- SCTransform(zenodo4059315[[i]], method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = T)
  DefaultAssay(zenodo4059315[[i]]) <- "SCT"
}
# Filter out samples that contain less than 100 cells
# tmp_ls <- c()
# for (i in 1:length(zenodo4059315)) {
#   if (length(colnames(zenodo4059315[[i]])) < 100) {
#     tmp_ls <- c(tmp_ls, i)
#   }
# }
# zenodo4059315 <- zenodo4059315[-tmp_ls]
# Now we can run the PCA on every remaining sample (samples with > 100 cells)
# for (i in 1:length(zenodo4059315)) {
#   zenodo4059315[[i]] <- RunPCA(zenodo4059315[[i]], features = VariableFeatures(zenodo4059315[[i]]), assay = "SCT", npcs = 50, verbose = T)
# }
integration_features <- SelectIntegrationFeatures(object.list = zenodo4059315, nfeatures = 2500)
zenodo4059315 <- PrepSCTIntegration(object.list = zenodo4059315, anchor.features = integration_features)
zenodo4059315_anchors <- FindIntegrationAnchors(object.list = zenodo4059315, normalization.method = "SCT", anchor.features = integration_features) # reduction = "rpca"
saveRDS(zenodo4059315_anchors, paste(RData_dir, "6_sc_zenodo4059315_anchors_backup.rds", sep = ""))
# zenodo4059315_anchors <- readRDS(paste(RData_dir, "6_sc_zenodo4059315_anchors_backup.rds", sep = ""))
rm(zenodo4059315)
zenodo4059315_integrated <- IntegrateData(anchorset = zenodo4059315_anchors, normalization.method = "SCT", new.assay.name = "seurat.integration")
rm(zenodo4059315_anchors)
saveRDS(zenodo4059315_integrated, paste(RData_dir, "6_sc_zenodo4059315_integrated_backup.rds", sep = "")) 
# zenodo4059315_integrated <- readRDS(paste(RData_dir, "6_sc_zenodo4059315_integrated_backup.rds", sep = ""))
zenodo4059315_integrated <- RunMCA(zenodo4059315_integrated, nmcs = 50, features = VariableFeatures(object = zenodo4059315_integrated), assay = "SCT", verbose = T) # assay = "RNA"
# Switch back to the "seurat.integration" assay as default, to compute every next step with this assay
DefaultAssay(zenodo4059315_integrated) <- "seurat.integration"
zenodo4059315_integrated <- RunPCA(zenodo4059315_integrated, features = VariableFeatures(zenodo4059315_integrated@assays$seurat.integration), assay = "seurat.integration", npcs = 50, verbose = T)
zenodo4059315_integrated <- RunUMAP(zenodo4059315_integrated, reduction = "pca", dims = 1:30, verbose = T)
zenodo4059315_integrated <- RunTSNE(zenodo4059315_integrated, reduction = "pca", dims = 1:30, verbose = T)
zenodo4059315_integrated <- RunMCUMAP(zenodo4059315_integrated, reduction = "mca", dims = 1:30, verbose = T)
DimPlot(zenodo4059315_integrated, reduction = "mcumap", group.by = "Annotation.3", label = T, repel = T, pt.size = 0.1)
# Retrieve labeled clusters from the Kuppe C, et al.
old_cluster_names <- c("Arteriolar Endothelium", "B Cells", "Collecting Duct Principal Cells", "Connecting Tubule", "Dendritic Cells", "Descending Thin Limb", "Distal Convoluted Tubule", "Fibroblast 2", "Fibroblast 4", "Fibroblast 6", "Glomerular Capillaries", "Injured Endothelial Cells", "Injured Proximal tubule", "Intercalated Cells 3", "Intercalated Cells 4", "Intercalated Cells 5", "Intercalated Cells 6", "Intercalated Cells 7", "Intercalated Cells 8", "Intercalated Cells A", "Intercalated Cells B", "Lymph Endothelium", "Macrophages 1", "Macrophages 2", "Macrophages 3", "Macula Densa Cells", "Mast Cells", "Monocytes", "Myofibroblast 1a", "Myofibroblast 1b", "Natural Killer Cells", "Pericytes 1", "Pericytes 2", "Plasma Cells", "Podocytes", "Proximal Tubule", "S1", "S1/2 1", "S1/2 2", "S1/2 3", "S3 1", "S3 2", "S3 3", "Schwann Cells", "T Cells", "Thick Ascending Limb 2", "Thick Ascending Limb 3", "Thick Ascending Limb 4", "Uroethlial Cells", "Vasa Recta 1", "Vasa Recta 2", "Vasa Recta 3", "Vasa Recta 4", "Vasa Recta 5", "Vasa Recta 6", "Vascular Smooth Muscle Cells", "Venular Endothelium")
new_cluster_names <- c("EC.art", "B.cells", "PC.CCD", "CNT", "DC", "LoH.DTL", "DCT", "Fibro.", "Fibro.", "Fibro.", "EC.glom", "EC.na", "PTC.na", "IC.na", "IC.na", "IC.na", "IC.na", "IC.na", "IC.na", "IC.A", "IC.B", "EC.lym", "Macro.", "Macro.", "Macro.", "MD.cells", "Mast.cells", "Mono.", "Myofibro.", "Myofibro.", "NK.cells", "Pericytes", "Pericytes", "Plasma.cells", "Podo.", "PTC.na", "PTC.S1", "PTC.S2", "PTC.S2", "PTC.S2", "PTC.S3", "PTC.S3", "PTC.S3", "Schwann.cells", "T.cells", "LoH.TAL", "LoH.TAL", "LoH.TAL", "Uro.", "EC.vasa.recta", "EC.vasa.recta", "EC.vasa.recta", "EC.vasa.recta", "EC.vasa.recta", "EC.vasa.recta", "vSMC", "EC.vei")
zenodo4059315_integrated@meta.data$annot_clusters <- plyr::mapvalues(x = zenodo4059315_integrated@meta.data$Annotation.3, from = old_cluster_names, to = new_cluster_names) # replacing cluster names by cell types
sorted_cell_types <- c("Mono.", "Macro.", "Mast.cells", "DC", "B.cells", "T.cells", "NK.cells", "Plasma.cells", "EC.vei", "EC.glom", "EC.art", "EC.lym", "EC.vasa.recta", "EC.na", "vSMC", "Pericytes", "Fibro.", "Myofibro.", "Podo.", "MD.cells", "PTC.S1", "PTC.S2", "PTC.S3", "PTC.na", "LoH.DTL", "LoH.TAL", "DCT", "CNT", "PC.CCD", "IC.A", "IC.B", "IC.na", "Schwann.cells", "Uro.")
zenodo4059315_integrated$annot_clusters <- factor(x = zenodo4059315_integrated$annot_clusters, levels = sorted_cell_types)
DimPlot(zenodo4059315_integrated, reduction = "mcumap", group.by = "annot_clusters", label = T, repel = T, pt.size = 0.1)
colors = c("brown1", "brown3", "green", "darkgoldenrod", "cornflowerblue", "cyan2", "cyan4", "cornflowerblue", "firebrick3", "firebrick2", "firebrick4", "yellow", "firebrick1", "grey90", "coral", "coral2", "coral3", "coral4", "darkorange", "pink", "darkseagreen2", "darkseagreen3", "darkseagreen4", "grey90", "deepskyblue2", "deepskyblue4", "darksalmon", "sandybrown", "rosybrown2", "plum3", "plum4", "grey90", "orange", "magenta")
p1 <- DimPlot(zenodo4059315_integrated, reduction = "mcumap", group.by = "annot_clusters", pt.size = 0.1, label = T, repel = T, cols = colors) + ggtitle("Kuppe C, et al. labeling") # display the profile of the integrated dataset
p1 %>% ggsave(filename = paste(figure_dir, "6_sc_zenodo4059315_integrated_AnnotClusters_DimPlot_colors.png", sep=""), width = 200, height = 150, units = "mm")
p1[[1]]$layers[[1]]$aes_params$alpha = 0.5
p2 <- DimPlot(zenodo4059315_integrated, reduction = "mcumap", group.by = "Patient.ID", label = F, pt.size = 0.1) + ggtitle("Samples") # check also batch integration
plot_grid(p1, p2, ncol = 2, rel_widths = c(1, 1.1)) %>% ggsave(filename = paste(figure_dir, "6_sc_zenodo4059315_integrated_AnnotClusters_vs_OrigIdent.png", sep=""), width = 450, height = 150, units = "mm")
# Save again zenodo4059315_integrated once everything is ok
saveRDS(zenodo4059315_integrated, paste(RData_dir, "6_sc_zenodo4059315_integrated_backup.rds", sep = "")) 
# zenodo4059315_integrated <- readRDS(paste(RData_dir, "6_sc_zenodo4059315_integrated_backup.rds", sep = ""))
# Now try to identify cell types, based on the consensus signatures
library(CelliD)
set.seed(1) 
# Setting the signature of cell types as a list
signatures <- as.list(read.table("/data-cbl/mquatre/p2_scrnaseq_renal_landscape/data/external/220317_p2_signatures_scrnaseq.csv", header = T, sep = ";", fill = TRUE))
DefaultAssay(zenodo4059315_integrated) <- "SCT"
# Perform HGT-based enrichment to identify cell types
log10_pval_matrix_is <- RunCellHGT(zenodo4059315_integrated, reduction = "mca", pathways = signatures, dims = 1:50, log.trans = T, n.features = 500, minSize = 10)
# Predict cell type based on signatures
prediction <- rownames(log10_pval_matrix_is)[apply(log10_pval_matrix_is, 2, FUN = which.max)]
# Remove low p-value predictions > 0.01
prediction <- ifelse(apply(log10_pval_matrix_is, 2, FUN = max) > 2, prediction, "unassigned") # this is the step where it's possible to play with p-val threshold, even if it's not adviced
# Store the prediction as metadata of the Seurat object
zenodo4059315_integrated@meta.data$PredHGT.p2 <- prediction
zenodo4059315_integrated@assays[["PredHGT.p2"]] <- CreateAssayObject(data = log10_pval_matrix_is)
# Then plot predicted cell types
DimPlot(zenodo4059315_integrated, reduction = "mcumap", group.by = "PredHGT.p2", label = T, repel = T, pt.size = 0.1)
sorted_cell_types <- c("Mono.", "Macro.", "Neutro.", "DC", "B.cells", "CD4.T.cells", "CD8.T.cells", "NK.cells", "EC.vei", "EC.glom", "EC.art", "vSMC", "Podo.", "PEC", "PTC.S1", "PTC.S2", "PTC.S3", "LoH.DTL", "LoH.ATL", "LoH.TAL", "DCT", "CNT", "PC.CCD", "PC.IMCD", "IC.A", "IC.B", "unassigned")
zenodo4059315_integrated$PredHGT.p2 <- factor(x = zenodo4059315_integrated$PredHGT.p2, levels = sorted_cell_types)
colors <- c("brown1", "brown3", "violetred", "darkgoldenrod", "cornflowerblue", "cyan2", "cyan3", "cyan4", "firebrick3", "firebrick2", "firebrick4", "coral", "darkorange", "blueviolet", "darkseagreen2", "darkseagreen3", "darkseagreen4", "deepskyblue2", "deepskyblue3", "deepskyblue4", "darksalmon", "sandybrown", "rosybrown2", "rosybrown4", "plum3", "plum4", "grey90")
p1 <- DimPlot(zenodo4059315_integrated, reduction = "mcumap", group.by = "PredHGT.p2", pt.size = 0.1, label = T, repel = T, cols = colors) + ggtitle("Consensus signatures") # display the profile of the integrated dataset
p1[[1]]$layers[[1]]$aes_params$alpha = 0.5
p1 %>% ggsave(filename = paste(figure_dir, "6_sc_zenodo4059315_integrated_AnnotClusters_PredHGT_DimPlot_colors.png", sep=""), width = 200, height = 150, units = "mm")
saveRDS(zenodo4059315_integrated, paste(RData_dir, "6_sc_zenodo4059315_PredHGT_integrated_backup.rds", sep = "")) 


# ------------------------------ #
# B. For snRNA-seq dataset
# ------------------------------ #
# First, load the query dataset (that must contain metadata with the labels from authors)
# Lake BB, et al. A single-nucleus RNA-sequencing pipeline to decipher the molecular anatomy and pathophysiology of human kidneys. | PMID 31249312 | GSE121852
gse121862 <- read.table("/data-cbl/mquatre/p2_scrnaseq_renal_landscape/data/raw/GSE121862/GSE121862_UCSD-WU_Single_Nuclei_Cluster_Annotated_Raw_UMI_Matrix.tsv", sep = "\t", row.names = 1, header = T)
gse121862 <- CreateSeuratObject(gse121862, project = "GSE121862", min.cells = 10, min.features = 200)
# Retrieve labeled clusters from the Lake BB, et al.
gse121862@meta.data$seurat_clusters <- gse121862@meta.data$orig.ident
old_cluster_names <- c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10", "C11", "C12", "C13", "C14", "C15", "C16", "C17", "C18", "C19", "C20", "C21", "C22", "C23", "C24", "C25", "C26", "C27", "C28", "C29", "C30")
new_cluster_names <- c("Epi.na", "Podo.", "PTC.S1", "PTC.S2", "PTC.na", "PTC.na", "PTC.S3", "LoH.DTL", "LoH.ATL", "LoH.ATL", "LoH.ATL", "LoH.TAL", "LoH.TAL", "DCT", "CNT", "PC.CCD", "PC.na", "PC.IMCD", "IC.A", "IC.A", "IC.B", "EC.glom", "EC.vei", "EC.vei", "EC.na", "Mes.", "vSMC", "Fibro.", "PTC.na", "Macro.")
gse121862@meta.data$annot_clusters <- plyr::mapvalues(x = gse121862@meta.data$seurat_clusters, from = old_cluster_names, to = new_cluster_names) # replacing cluster names by cell types
# ... and retrieve the sample origin of every single nucleus
tmp_ls <- c()
for (i in 1:length(colnames(gse121862))) {
  tmp_tx <- as.list(str_split(colnames(gse121862)[[i]], "_", n = Inf, simplify = T))
  tmp_ls <- c(tmp_ls, tmp_tx[[2]])
}
gse121862@meta.data$orig.ident <- tmp_ls
rm(tmp_tx, tmp_ls)
# ... finally, keep only control samples (donor kidney from MAT, see GSE121862_family.soft file)
# gse121862 <- subset(gse121862, subset = orig.ident == "NK37" | orig.ident == "NK38" | orig.ident == "NK41" | orig.ident == "NK42" | orig.ident == "NK44" | orig.ident == "NK45" | orig.ident == "NK46" | orig.ident == "NK93" | orig.ident == "NK94")
# ------------------------------ #
# PROCESS DATA
# First add some metadata
gse121862[["percent.mt"]] <- PercentageFeatureSet(gse121862, pattern = "^MT-")
gse121862 <- subset(gse121862, subset = nFeature_RNA > 200 & nFeature_RNA < 3500 & percent.mt < 5)
# Then process the data with normalization, scaling (using SCTransform) and dimension reductions
DefaultAssay(gse121862) <- "RNA"
gse121862 <- SplitObject(gse121862, split.by = "orig.ident")
for (i in 1:length(gse121862)) {
  gse121862[[i]] <- SCTransform(gse121862[[i]], method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = T)
  DefaultAssay(gse121862[[i]]) <- "SCT"
}
# Filter out samples that contain less than 100 cells
tmp_ls <- c()
for (i in 1:length(gse121862)) {
  if (length(colnames(gse121862[[i]])) < 100) {
    tmp_ls <- c(tmp_ls, i)
  }
}
gse121862 <- gse121862[-tmp_ls]
# Now we can run the PCA on every remaining sample (samples with > 100 cells)
# for (i in 1:length(gse121862)) {
#   gse121862[[i]] <- RunPCA(gse121862[[i]], features = VariableFeatures(gse121862[[i]]), assay = "SCT", npcs = 50, verbose = T)
# }
integration_features <- SelectIntegrationFeatures(object.list = gse121862, nfeatures = 2500)
gse121862 <- PrepSCTIntegration(object.list = gse121862, anchor.features = integration_features)
gse121862_anchors <- FindIntegrationAnchors(object.list = gse121862, normalization.method = "SCT", anchor.features = integration_features) # reduction = "rpca"
saveRDS(gse121862_anchors, paste(RData_dir, "6_sn_gse121862_anchors_backup.rds", sep = ""))
# gse121862_anchors <- readRDS(paste(RData_dir, "6_sn_gse121862_anchors_backup.rds", sep = ""))
rm(gse121862)
gse121862_integrated <- IntegrateData(anchorset = gse121862_anchors, normalization.method = "SCT", new.assay.name = "seurat.integration")
rm(gse121862_anchors)
saveRDS(gse121862_integrated, paste(RData_dir, "6_sn_gse121862_integrated_backup.rds", sep = "")) 
# gse121862_integrated <- readRDS(paste(RData_dir, "6_sn_gse121862_integrated_backup.rds", sep = ""))
gse121862_integrated <- RunMCA(gse121862_integrated, nmcs = 50, features = VariableFeatures(object = gse121862_integrated), assay = "SCT", verbose = T) # assay = "RNA"
# Switch back to the "seurat.integration" assay as default, to compute every next step with this assay
DefaultAssay(gse121862_integrated) <- "seurat.integration"
gse121862_integrated <- RunPCA(gse121862_integrated, features = VariableFeatures(gse121862_integrated@assays$seurat.integration), assay = "seurat.integration", npcs = 50, verbose = T)
gse121862_integrated <- RunUMAP(gse121862_integrated, reduction = "pca", dims = 1:30, verbose = T)
gse121862_integrated <- RunTSNE(gse121862_integrated, reduction = "pca", dims = 1:30, verbose = T)
gse121862_integrated <- RunMCUMAP(gse121862_integrated, reduction = "mca", dims = 1:30, verbose = T)
DimPlot(gse121862_integrated, reduction = "mcumap", group.by = "annot_clusters", label = T, repel = T, pt.size = 0.1)
sorted_cell_types <- c("Macro.", "EC.vei", "EC.glom", "EC.na", "vSMC", "Fibro.", "Mes.", "Podo.", "Epi.na", "PTC.S1", "PTC.S2", "PTC.S3", "PTC.na", "LoH.DTL", "LoH.ATL", "LoH.TAL", "DCT", "CNT", "PC.CCD", "PC.IMCD", "PC.na", "IC.A", "IC.B")
gse121862_integrated$annot_clusters <- factor(x = gse121862_integrated$annot_clusters, levels = sorted_cell_types)
colors = c("brown3", "firebrick3", "firebrick2", "grey90", "coral", "coral3", "blueviolet", "darkorange", "grey90", "darkseagreen2", "darkseagreen3", "darkseagreen4", "grey90", "deepskyblue2", "deepskyblue3", "deepskyblue4", "darksalmon", "sandybrown", "rosybrown2", "rosybrown4", "grey90", "plum3", "plum4")
p1 <- DimPlot(gse121862_integrated, reduction = "mcumap", group.by = "annot_clusters", pt.size = 0.1, label = T, repel = T, cols = colors) + ggtitle("Lake BB, et al. labeling") # display the profile of the integrated dataset
p1 %>% ggsave(filename = paste(figure_dir, "6_sn_gse121862_integrated_AnnotClusters_DimPlot_colors.png", sep=""), width = 200, height = 150, units = "mm")
p1[[1]]$layers[[1]]$aes_params$alpha = 0.5
p2 <- DimPlot(gse121862_integrated, reduction = "mcumap", group.by = "orig.ident", label = F, pt.size = 0.1) + ggtitle("Samples") # check also batch integration
plot_grid(p1, p2, ncol = 2, rel_widths = c(1, 1.1)) %>% ggsave(filename = paste(figure_dir, "6_sn_gse121862_integrated_AnnotClusters_vs_OrigIdent.png", sep=""), width = 450, height = 150, units = "mm")
# Save again gse121862_integrated once everything is ok
saveRDS(gse121862_integrated, paste(RData_dir, "6_sn_gse121862_integrated_backup.rds", sep = "")) 
# gse121862_integrated <- readRDS(paste(RData_dir, "6_sn_gse121862_integrated_backup.rds", sep = ""))
# Now try to identify cell types, based on the consensus signatures
library(CelliD)
set.seed(1)
# Setting the signature of cell types as a list
signatures <- as.list(read.table("/data-cbl/mquatre/p2_scrnaseq_renal_landscape/data/external/220317_p2_signatures_snrnaseq.csv", header = T, sep = ";", fill = TRUE))
# signatures <- signatures[-c(13, 17)] # remove ".na" signatures
DefaultAssay(gse121862_integrated) <- "SCT"
# Perform HGT-based enrichment to identify cell types
log10_pval_matrix_is <- RunCellHGT(gse121862_integrated, reduction = "mca", pathways = signatures, dims = 1:50, log.trans = T, n.features = 500, minSize = 10)
# Predict cell type based on signatures
prediction <- rownames(log10_pval_matrix_is)[apply(log10_pval_matrix_is, 2, FUN = which.max)]
# Remove low p-value predictions > 0.01
prediction <- ifelse(apply(log10_pval_matrix_is, 2, FUN = max) > 2, prediction, "unassigned") # this is the step where it's possible to play with p-val threshold, even if it's not adviced
# Store the prediction as metadata of the Seurat object
gse121862_integrated@meta.data$PredHGT.p2 <- prediction
gse121862_integrated@assays[["PredHGT.p2"]] <- CreateAssayObject(data = log10_pval_matrix_is)
# Then plot predicted cell types
DimPlot(gse121862_integrated, reduction = "mcumap", group.by = "PredHGT.p2", label = T, repel = T, pt.size = 0.1)
sorted_cell_types <- c("Macro.", "EC.vei", "EC.glom", "EC.art", "vSMC", "Pericytes", "Fibro.", "Podo.", "PEC", "PTC.S1", "PTC.S2", "PTC.S3", "LoH.DTL", "LoH.ATL", "LoH.TAL", "DCT", "CNT", "PC.CCD", "PC.OMCD", "PC.IMCD", "IC.A", "IC.B", "unassigned")
gse121862_integrated$PredHGT.p2 <- factor(x = gse121862_integrated$PredHGT.p2, levels = sorted_cell_types)
colors <- c("brown3", "firebrick3", "firebrick2", "firebrick4", "coral", "coral2", "coral3", "darkorange", "blueviolet", "darkseagreen2", "darkseagreen3", "darkseagreen4", "deepskyblue2", "deepskyblue3", "deepskyblue4", "darksalmon", "sandybrown", "rosybrown2", "rosybrown3", "rosybrown4", "plum3", "plum4", "grey90")
p1 <- DimPlot(gse121862_integrated, reduction = "mcumap", group.by = "PredHGT.p2", pt.size = 0.1, label = T, repel = T, cols = colors) + ggtitle("Consensus signatures") # display the profile of the integrated dataset
p1[[1]]$layers[[1]]$aes_params$alpha = 0.5
p1 %>% ggsave(filename = paste(figure_dir, "6_sn_gse121862_integrated_AnnotClusters_PredHGT_DimPlot_colors.png", sep=""), width = 200, height = 150, units = "mm")
saveRDS(gse121862_integrated, paste(RData_dir, "6_sn_gse121862_integrated_PredHGT_backup.rds", sep = ""))


# --------------------------------------------------
# Figures pour l'article p2 - Quatredeniers, et al.
# --------------------------------------------------
# 1. scRNA-seq
# --------------------------------------------------
sc_gsm_integrated <- readRDS(paste(RData_dir, "3_sc_gsm_integrated_SeuratAnnot_backup.rds", sep = ""))
sc_gsm_merge <- readRDS(paste(RData_dir, "2_sc_gsm_integrated_harmony_backup.rds", sep = ""))
# --------------------------------------------------
# Figure 2
p1 <- DimPlot(sc_gsm_integrated, reduction = "umap", pt.size = 0.1, group.by = "orig.ident", label = F) + NoLegend()
p2 <- DimPlot(sc_gsm_integrated, reduction = "umap", pt.size = 0.1, group.by = "batch.ident", label = F) + NoLegend()
sc_gsm_integrated$sex.ident <- factor(x = sc_gsm_integrated$sex.ident, levels = c("M", "F", "?"))
p3 <- DimPlot(sc_gsm_integrated, reduction = "umap", pt.size = 0.1, group.by = "sex.ident", label = F, cols = c("lightcoral", "lightcyan4", "grey90")) + NoLegend()
p4 <- DimPlot(sc_gsm_integrated, reduction = "umap", pt.size = 0.1, group.by = "seurat_clusters", label = T, repel = T) + NoLegend()
colors <- c("brown1", "brown3", "violetred", "darkgoldenrod", "cornflowerblue", "cyan2", "cyan3", "cyan4", "firebrick3", "firebrick2", "firebrick4", "grey90", "coral", "darkorange", "blueviolet", "darkseagreen2", "darkseagreen3", "darkseagreen4", "grey90", "deepskyblue2", "deepskyblue3", "deepskyblue4", "grey90", "darksalmon", "sandybrown", "rosybrown2", "rosybrown3", "rosybrown4", "grey90", "plum3", "plum4", "plum1")
p5 <- DimPlot(sc_gsm_integrated, reduction = "umap", pt.size = 0.1, group.by = "annot_clusters", label = T, repel = T, cols = colors) + NoLegend()
p1[[1]]$layers[[1]]$aes_params$alpha = 0.5
p2[[1]]$layers[[1]]$aes_params$alpha = 0.5
p3[[1]]$layers[[1]]$aes_params$alpha = 0.5
plot_grid(p1, p2, p3, ncol = 3, rel_widths = c(1, 1, 1)) %>% ggsave(filename = paste(figure_dir, "FigS2BCD.png", sep=""), width = 300, height = 100, units = "mm")
p5[[1]]$layers[[1]]$aes_params$alpha = 0.5
p1 %>% ggsave(filename = paste(figure_dir, "FigS2B.png", sep=""), width = 150, height = 150, units = "mm")
p2 %>% ggsave(filename = paste(figure_dir, "FigS2C.png", sep=""), width = 150, height = 150, units = "mm")
p3 %>% ggsave(filename = paste(figure_dir, "FigS2D.png", sep=""), width = 150, height = 150, units = "mm")
p4 %>% ggsave(filename = paste(figure_dir, "Fig2A.png", sep=""), width = 150, height = 150, units = "mm")
p5 %>% ggsave(filename = paste(figure_dir, "Fig2B.png", sep=""), width = 150, height = 150, units = "mm")
# Legends
p1 <- DimPlot(sc_gsm_integrated, reduction = "umap", pt.size = 0.1, group.by = "orig.ident", label = F) + theme(legend.position = "bottom")
p2 <- DimPlot(sc_gsm_integrated, reduction = "umap", pt.size = 0.1, group.by = "batch.ident", label = F) + theme(legend.position = "bottom")
p3 <- DimPlot(sc_gsm_integrated, reduction = "umap", pt.size = 0.1, group.by = "sex.ident", label = F, cols = c("lightcoral", "lightcyan4", "grey90")) + theme(legend.position = "bottom")
p4 <- DimPlot(sc_gsm_integrated, reduction = "umap", pt.size = 0.1, group.by = "seurat_clusters", label = T, repel = T)
p5 <- DimPlot(sc_gsm_integrated, reduction = "umap", pt.size = 0.1, group.by = "annot_clusters", label = T, repel = T, cols = colors)
p1 %>% ggsave(filename = paste(figure_dir, "FigS2B_legend.png", sep=""), width = 220, height = 150, units = "mm")
p2 %>% ggsave(filename = paste(figure_dir, "FigS2C_legend.png", sep=""), width = 150, height = 150, units = "mm")
p3 %>% ggsave(filename = paste(figure_dir, "FigS2D_legend.png", sep=""), width = 150, height = 150, units = "mm")
p4 %>% ggsave(filename = paste(figure_dir, "Fig2A_legend.png", sep=""), width = 150, height = 150, units = "mm")
p5 %>% ggsave(filename = paste(figure_dir, "Fig2B_legend.png", sep=""), width = 150, height = 150, units = "mm")
# ViolinPlot & DotPlot
p1 <- VlnPlot(sc_gsm_integrated, group.by = "annot_clusters", assay = "SCT", features = c("PTPRC", "CD14", "LYZ", "S100A12", "FCN1", "FCGR3A", "CSF3R", "CD68", "FCGR3B", "IL7R", "CLEC10A", "FCER1A", "MS4A1", "CD79A", "CXCR4", "CD3D", "CD69", "CCR7", "GZMA", "GNLY", "CD8A", "TRAC", "NKG7", "EMCN", "ENG", "PLVAP", "EHD3", "PLAT", "TSPAN7", "SOX17", "CAV1", "ACTA2", "PDGFRB", "MYH11", "TAGLN", "NPHS2", "PODXL", "CLIC5", "NPHS1", "WT1", "CTGF", "PTGDS", "PAX2", "CLDN1", "PROM1", "OCIAD2", "APOE", "ALDOB", "SLC5A2", "SLC5A12", "SLC22A6", "SLC22A8", "MIOX", "RBP4", "CRYAB", "SPP1", "SLC14A2", "CLDN3", "TACSTD2", "UMOD", "SLC12A1", "DEFB1", "KNG1", "SLC8A1", "SLC12A3", "CALB1", "AQP2", "AQP3", "FXYD4", "ATP6V1G3", "FOXI1", "SLC4A1", "DMRT2", "SLC26A4", "INSRR"), stack = T, fill.by = "ident", cols = colors) + NoLegend()
p1 %>% ggsave(filename = paste(figure_dir, "Fig2C_VlnPlot.png", sep=""), width = 405, height = 190, units = "mm")
sc_annot_markers <- read.table(paste(reports_dir, "3_sc_gsm_integrated_AnnotClusters_CellTypeMarkers_SCT.csv", sep=""), sep=",")
sc_top_annot_markers <- sc_annot_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
p1 <- DotPlot(sc_gsm_integrated, features = make.unique(sc_top_annot_markers$gene), assay = "SCT", cols = c("green", "red")) + RotatedAxis()
p1 %>% ggsave(filename = paste(figure_dir, "Fig2D_DotPlot.png", sep = ""), width = 450, height = 200, units = "mm")
# Supp. Figure 2: PCA / Harmony / Seurat v4
p1 <- DimPlot(sc_gsm_merge, reduction = "pca", group.by = "orig.ident", pt.size = 0.1, raster = F) + theme(legend.position = "bottom")
p1 %>% ggsave(filename = paste(figure_dir, "FigS2A_PCA_Harmony_SeuratV4_legend.png", sep=""), width = 300, height = 100, units = "mm")
p1 <- DimPlot(sc_gsm_merge, reduction = "pca", group.by = "orig.ident", pt.size = 0.1, raster = F) + NoLegend() + ggtitle("Batch effects") 
p2 <- DimPlot(sc_gsm_merge, reduction = "harmony", group.by = "orig.ident", pt.size = 0.1, raster = F) + NoLegend() + ggtitle("Harmony correction")
p3 <- DimPlot(sc_gsm_integrated, reduction = "pca", group.by = "orig.ident", pt.size = 0.1, raster = F) + NoLegend() + ggtitle("Seurat.v4 correction")
p1[[1]]$layers[[1]]$aes_params$alpha = 0.5
p2[[1]]$layers[[1]]$aes_params$alpha = 0.5
p3[[1]]$layers[[1]]$aes_params$alpha = 0.5
p4 <- plot_grid(p1, p2, p3, ncol = 3, rel_widths = c(1, 1, 1)) 
p4 %>% ggsave(filename = paste(figure_dir, "FigS2A_PCA_Harmony_SeuratV4.png", sep=""), width = 300, height = 100, units = "mm")
rm(sc_gsm_integrated, sc_gsm_merge)
# --------------------------------------------------
# 2. snRNA-seq
# --------------------------------------------------
sn_gsm_integrated <- readRDS(paste(RData_dir, "3_sn_gsm_integrated_SeuratAnnot_backup.rds", sep = ""))
sn_gsm_merge <- readRDS(paste(RData_dir, "2_sn_gsm_integrated_harmony_backup.rds", sep = ""))
# --------------------------------------------------
# Figure 2
p1 <- DimPlot(sn_gsm_integrated, reduction = "umap", pt.size = 0.1, group.by = "orig.ident", label = F) + NoLegend()
p2 <- DimPlot(sn_gsm_integrated, reduction = "umap", pt.size = 0.1, group.by = "batch.ident", label = F) + NoLegend()
sn_gsm_integrated$sex.ident <- factor(x = sn_gsm_integrated$sex.ident, levels = c("M", "F"))
p3 <- DimPlot(sn_gsm_integrated, reduction = "umap", pt.size = 0.1, group.by = "sex.ident", label = F, cols = c("lightcoral", "lightcyan4")) + NoLegend()
p4 <- DimPlot(sn_gsm_integrated, reduction = "umap", pt.size = 0.1, group.by = "seurat_clusters", label = T, repel = T) + NoLegend()
colors <- c("brown3", "firebrick3", "firebrick2", "firebrick4", "grey90", "coral", "coral2", "coral3", "darkorange", "blueviolet", "darkseagreen2", "darkseagreen3", "darkseagreen4", "grey90", "deepskyblue2", "deepskyblue3", "deepskyblue4", "grey90", "darksalmon", "sandybrown", "rosybrown2", "rosybrown3", "rosybrown4", "plum3", "plum4")
p5 <- DimPlot(sn_gsm_integrated, reduction = "umap", pt.size = 0.1, group.by = "annot_clusters", label = T, repel = T, cols = colors) + NoLegend()
p1[[1]]$layers[[1]]$aes_params$alpha = 0.5
p2[[1]]$layers[[1]]$aes_params$alpha = 0.5
p3[[1]]$layers[[1]]$aes_params$alpha = 0.5
plot_grid(p1, p2, p3, ncol = 3, rel_widths = c(1, 1, 1)) %>% ggsave(filename = paste(figure_dir, "FigS3BCD.png", sep=""), width = 300, height = 100, units = "mm")
p5[[1]]$layers[[1]]$aes_params$alpha = 0.5
p1 %>% ggsave(filename = paste(figure_dir, "FigS3B.png", sep=""), width = 150, height = 150, units = "mm")
p2 %>% ggsave(filename = paste(figure_dir, "FigS3C.png", sep=""), width = 150, height = 150, units = "mm")
p3 %>% ggsave(filename = paste(figure_dir, "FigS3D.png", sep=""), width = 150, height = 150, units = "mm")
p4 %>% ggsave(filename = paste(figure_dir, "Fig3A.png", sep=""), width = 150, height = 150, units = "mm")
p5 %>% ggsave(filename = paste(figure_dir, "Fig3B.png", sep=""), width = 150, height = 150, units = "mm")
# Legends
p1 <- DimPlot(sn_gsm_integrated, reduction = "umap", pt.size = 0.1, group.by = "orig.ident", label = F) + theme(legend.position = "bottom")
p2 <- DimPlot(sn_gsm_integrated, reduction = "umap", pt.size = 0.1, group.by = "batch.ident", label = F) + theme(legend.position = "bottom")
p3 <- DimPlot(sn_gsm_integrated, reduction = "umap", pt.size = 0.1, group.by = "sex.ident", label = F, cols = c("lightcoral", "lightcyan4")) + theme(legend.position = "bottom")
p4 <- DimPlot(sn_gsm_integrated, reduction = "umap", pt.size = 0.1, group.by = "seurat_clusters", label = T, repel = T)
p5 <- DimPlot(sn_gsm_integrated, reduction = "umap", pt.size = 0.1, group.by = "annot_clusters", label = T, repel = T, cols = colors)
p1 %>% ggsave(filename = paste(figure_dir, "FigS3B_legend.png", sep=""), width = 220, height = 150, units = "mm")
p2 %>% ggsave(filename = paste(figure_dir, "FigS3C_legend.png", sep=""), width = 220, height = 150, units = "mm")
p3 %>% ggsave(filename = paste(figure_dir, "FigS3D_legend.png", sep=""), width = 150, height = 150, units = "mm")
p4 %>% ggsave(filename = paste(figure_dir, "Fig3A_legend.png", sep=""), width = 150, height = 150, units = "mm")
p5 %>% ggsave(filename = paste(figure_dir, "Fig3B_legend.png", sep=""), width = 150, height = 150, units = "mm")
# ViolinPlot & DotPlot
p1 <- VlnPlot(sn_gsm_integrated, group.by = "annot_clusters", assay = "SCT", features = c("PTPRC", "CD14", "CD68", "ITGAM", "PECAM1", "FLT1", "EMCN", "ENG", "PLVAP", "EHD3", "PLAT", "TSPAN7", "CAV1", "PDGFRB", "EDNRA", "ACTA2", "NPHS2", "PODXL", "CLIC5", "WT1", "NPHS1", "CTGF", "CLDN1", "PROM1", "CUBN", "GPX3", "HNF4A", "ALDOB", "SLC5A12", "SLC22A6", "SLC22A8", "MIOX", "CRYAB", "SPP1", "S100A2", "S100A6", "UMOD", "SLC12A1", "DEFB1", "KNG1", "SLC12A3", "SLC8A1", "CALB1", "AQP2", "AQP3", "FXYD4", "ATP6V1G3", "FOXI1", "SLC4A1", "DMRT2", "SLC26A4", "INSRR"), stack = T, fill.by = "ident", cols = colors) + NoLegend()
p1 %>% ggsave(filename = paste(figure_dir, "Fig3C_VlnPlot.png", sep=""), width = 405, height = 155, units = "mm")
sc_annot_markers <- read.table(paste(reports_dir, "3_sn_gsm_integrated_AnnotClusters_CellTypeMarkers_SCT.csv", sep=""), sep=",")
sc_top_annot_markers <- sc_annot_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
p1 <- DotPlot(sn_gsm_integrated, features = make.unique(sc_top_annot_markers$gene), assay = "SCT", cols = c("green", "red")) + RotatedAxis()
p1 %>% ggsave(filename = paste(figure_dir, "Fig3D_DotPlot.png", sep = ""), width = 450, height = 160, units = "mm")
# Supp. Figure 2: PCA / Harmony / Seurat v4
p1 <- DimPlot(sn_gsm_merge, reduction = "pca", group.by = "orig.ident", pt.size = 0.1, raster = F) + theme(legend.position = "bottom")
p1 %>% ggsave(filename = paste(figure_dir, "FigS3A_PCA_Harmony_SeuratV4_legend.png", sep=""), width = 300, height = 100, units = "mm")
p1 <- DimPlot(sn_gsm_merge, reduction = "pca", group.by = "orig.ident", pt.size = 0.1, raster = F) + NoLegend() + ggtitle("Batch effects") 
p2 <- DimPlot(sn_gsm_merge, reduction = "harmony", group.by = "orig.ident", pt.size = 0.1, raster = F) + NoLegend() + ggtitle("Harmony correction")
p3 <- DimPlot(sn_gsm_integrated, reduction = "pca", group.by = "orig.ident", pt.size = 0.1, raster = F) + NoLegend() + ggtitle("Seurat.v4 correction")
p1[[1]]$layers[[1]]$aes_params$alpha = 0.5
p2[[1]]$layers[[1]]$aes_params$alpha = 0.5
p3[[1]]$layers[[1]]$aes_params$alpha = 0.5
p4 <- plot_grid(p1, p2, p3, ncol = 3, rel_widths = c(1, 1, 1)) 
p4 %>% ggsave(filename = paste(figure_dir, "FigS3A_PCA_Harmony_SeuratV4.png", sep=""), width = 300, height = 100, units = "mm")
rm(sn_gsm_integrated, sn_gsm_merge)
# --------------------------------------------------
# 3. scsn_integrated
# --------------------------------------------------
scsn_integrated <- readRDS(paste(RData_dir, "5_scsn_integrated_backup3.rds", sep = ""))
scsn_gsm_merge <- readRDS(paste(RData_dir, "5_scsn_gsm_integrated_harmony_backup.rds", sep = ""))
# Figures
colors <- c("brown1", "brown3", "violetred", "darkgoldenrod", "cornflowerblue", "cyan2", "cyan3", "cyan4", "firebrick3", "firebrick2", "firebrick4", "grey90", "coral", "coral2", "coral3", "darkorange", "blueviolet", "darkseagreen2", "darkseagreen3", "darkseagreen4", "grey90", "deepskyblue2", "deepskyblue3", "deepskyblue4", "grey90", "darksalmon", "sandybrown", "rosybrown2", "rosybrown3", "rosybrown4", "grey90", "plum3", "plum4", "plum1")
p1 <- DimPlot(scsn_integrated, reduction = "umap", pt.size = 0.1, group.by = "orig.ident", label = F, raster = F) + NoLegend()
p2 <- DimPlot(scsn_integrated, reduction = "umap", pt.size = 0.1, group.by = "techno.ident", label = F, raster = F) + NoLegend()
p3 <- DimPlot(scsn_integrated, reduction = "umap", pt.size = 0.1, group.by = "sex.ident", label = F, raster = F, cols = c("lightcoral", "lightcyan4", "grey90")) + NoLegend()
p4 <- DimPlot(scsn_integrated, reduction = "umap", pt.size = 0.1, group.by = "annot_clusters", split.by = "techno.ident", label = T, repel = T, raster = F, cols = colors) + NoLegend()
p1[[1]]$layers[[1]]$aes_params$alpha = 0.5
p2[[1]]$layers[[1]]$aes_params$alpha = 0.5
p3[[1]]$layers[[1]]$aes_params$alpha = 0.5
p4[[1]]$layers[[1]]$aes_params$alpha = 0.5
p1 %>% ggsave(filename = paste(figure_dir, "Fig4B.png", sep=""), width = 150, height = 150, units = "mm")
p2 %>% ggsave(filename = paste(figure_dir, "Fig4C.png", sep=""), width = 150, height = 150, units = "mm")
p3 %>% ggsave(filename = paste(figure_dir, "Fig4D.png", sep=""), width = 150, height = 150, units = "mm")
p4 %>% ggsave(filename = paste(figure_dir, "Fig4E.png", sep=""), width = 280, height = 150, units = "mm")
# Legends
p1 <- DimPlot(scsn_integrated, reduction = "umap", pt.size = 0.1, group.by = "orig.ident", label = F, raster = F)
p2 <- DimPlot(scsn_integrated, reduction = "umap", pt.size = 0.1, group.by = "techno.ident", label = F, raster = F)
p3 <- DimPlot(scsn_integrated, reduction = "umap", pt.size = 0.1, group.by = "sex.ident", label = F, raster = F, cols = c("lightcoral", "lightcyan4", "grey90")) + theme(legend.position = "bottom")
p4 <- DimPlot(scsn_integrated, reduction = "umap", pt.size = 0.1, group.by = "annot_clusters", split.by = "techno.ident", label = T, repel = T, raster = F, cols = colors)
p1 %>% ggsave(filename = paste(figure_dir, "Fig4B_legend.png", sep=""), width = 150, height = 150, units = "mm")
p2 %>% ggsave(filename = paste(figure_dir, "Fig4C_legend.png", sep=""), width = 150, height = 150, units = "mm")
p3 %>% ggsave(filename = paste(figure_dir, "Fig4D_legend.png", sep=""), width = 150, height = 150, units = "mm")
p4 %>% ggsave(filename = paste(figure_dir, "Fig4E_legend.png", sep=""), width = 280, height = 150, units = "mm")
# Figure 5A: PCA / Harmony / Seurat v4
p1 <- DimPlot(scsn_gsm_merge, reduction = "pca", group.by = "orig.ident", pt.size = 0.1, raster = F) + theme(legend.position = "bottom")
p1 %>% ggsave(filename = paste(figure_dir, "Fig4A_PCA_Harmony_SeuratV4_legend.png", sep=""), width = 300, height = 100, units = "mm")
p1 <- DimPlot(scsn_gsm_merge, reduction = "pca", group.by = "orig.ident", pt.size = 0.1, raster = F) + NoLegend() + ggtitle("Batch effects") 
p2 <- DimPlot(scsn_gsm_merge, reduction = "harmony", group.by = "orig.ident", pt.size = 0.1, raster = F) + NoLegend() + ggtitle("Harmony correction")
p3 <- DimPlot(scsn_integrated, reduction = "pca", group.by = "orig.ident", pt.size = 0.1, raster = F) + NoLegend() + ggtitle("Seurat.v4 correction")
p1[[1]]$layers[[1]]$aes_params$alpha = 0.5
p2[[1]]$layers[[1]]$aes_params$alpha = 0.5
p3[[1]]$layers[[1]]$aes_params$alpha = 0.5
p4 <- plot_grid(p1, p2, p3, ncol = 3, rel_widths = c(1, 1, 1)) 
p4 %>% ggsave(filename = paste(figure_dir, "Fig4A_PCA_Harmony_SeuratV4.png", sep=""), width = 300, height = 100, units = "mm")
rm(scsn_integrated, scsn_gsm_merge)
# --------------------------------------------------
# 4. zenodo4059315_integrated (sc)
# --------------------------------------------------
zenodo4059315_integrated <- readRDS(paste(RData_dir, "6_sc_zenodo4059315_PredHGT_integrated_backup.rds", sep = ""))
p1 <- DimPlot(zenodo4059315_integrated, reduction = "mcumap", pt.size = 0.1, group.by = "Patient.ID", label = F) # orig.ident
colors = c("brown1", "brown3", "green", "darkgoldenrod", "cornflowerblue", "cyan2", "cyan4", "cornflowerblue", "firebrick3", "firebrick2", "firebrick4", "yellow", "firebrick1", "grey90", "coral", "coral2", "coral3", "coral4", "darkorange", "pink", "darkseagreen2", "darkseagreen3", "darkseagreen4", "grey90", "deepskyblue2", "deepskyblue4", "darksalmon", "sandybrown", "rosybrown2", "plum3", "plum4", "grey90", "orange", "magenta")
p2 <- DimPlot(zenodo4059315_integrated, reduction = "mcumap", pt.size = 0.1, group.by = "annot_clusters", label = F, repel = F, cols = colors) + ggtitle("Labelling from Kuppe C, et al.") # Annotation.3
colors <- c("brown1", "brown3", "violetred", "darkgoldenrod", "cornflowerblue", "cyan2", "cyan3", "cyan4", "firebrick3", "firebrick2", "firebrick4", "coral", "darkorange", "blueviolet", "darkseagreen2", "darkseagreen3", "darkseagreen4", "deepskyblue2", "deepskyblue3", "deepskyblue4", "darksalmon", "sandybrown", "rosybrown2", "rosybrown4", "plum3", "plum4", "grey90")
p3 <- DimPlot(zenodo4059315_integrated, reduction = "mcumap", pt.size = 0.1, group.by = "PredHGT.p2", label = F, repel = F, cols = colors) + ggtitle("Enrichment of consensus signatures")
# Legends
p1 %>% ggsave(filename = paste(figure_dir, "Fig5A_legends.png", sep=""), width = 150, height = 150, units = "mm")
p2 %>% ggsave(filename = paste(figure_dir, "Fig5B_legends.png", sep=""), width = 150, height = 150, units = "mm")
p3 %>% ggsave(filename = paste(figure_dir, "Fig5C_legends.png", sep=""), width = 150, height = 150, units = "mm")
# Figures
p1[[1]]$layers[[1]]$aes_params$alpha = 0.5
p2[[1]]$layers[[1]]$aes_params$alpha = 0.5
p3[[1]]$layers[[1]]$aes_params$alpha = 0.5
p1 <- p1 + NoLegend() 
p1 %>% ggsave(filename = paste(figure_dir, "Fig5A.png", sep=""), width = 150, height = 150, units = "mm")
p2 <- p2 + NoLegend() 
p2 %>% ggsave(filename = paste(figure_dir, "Fig5B.png", sep=""), width = 150, height = 150, units = "mm")
p3 <- p3 + NoLegend() 
p3 %>% ggsave(filename = paste(figure_dir, "Fig5C.png", sep=""), width = 150, height = 150, units = "mm")
rm(zenodo4059315_integrated)
# --------------------------------------------------
# 5. gse121862_integrated (sn)
# --------------------------------------------------
gse121862_integrated <- readRDS(paste(RData_dir, "6_sn_gse121862_integrated_PredHGT_backup.rds", sep = ""))
p1 <- DimPlot(gse121862_integrated, reduction = "mcumap", pt.size = 0.1, group.by = "orig.ident", label = F)
colors = c("brown3", "firebrick3", "firebrick2", "grey90", "coral", "coral3", "blueviolet", "darkorange", "grey90", "darkseagreen2", "darkseagreen3", "darkseagreen4", "grey90", "deepskyblue2", "deepskyblue3", "deepskyblue4", "darksalmon", "sandybrown", "rosybrown2", "rosybrown4", "grey90", "plum3", "plum4")
p2 <- DimPlot(gse121862_integrated, reduction = "mcumap", pt.size = 0.1, group.by = "annot_clusters", label = F, repel = F, cols = colors) + ggtitle("Labelling from Lake BB, et al.")
colors <- c("brown3", "firebrick3", "firebrick2", "firebrick4", "coral", "coral2", "coral3", "darkorange", "blueviolet", "darkseagreen2", "darkseagreen3", "darkseagreen4", "deepskyblue2", "deepskyblue3", "deepskyblue4", "darksalmon", "sandybrown", "rosybrown2", "rosybrown3", "rosybrown4", "plum3", "plum4", "grey90")
p3 <- DimPlot(gse121862_integrated, reduction = "mcumap", pt.size = 0.1, group.by = "PredHGT.p2", label = F, repel = F, cols = colors) + ggtitle("Enrichment of consensus signatures")
# Legends
p1 %>% ggsave(filename = paste(figure_dir, "Fig5D_legends.png", sep=""), width = 150, height = 150, units = "mm")
p2 %>% ggsave(filename = paste(figure_dir, "Fig5E_legends.png", sep=""), width = 150, height = 150, units = "mm")
p3 %>% ggsave(filename = paste(figure_dir, "Fig5F_legends.png", sep=""), width = 150, height = 150, units = "mm")
# Figures
p1[[1]]$layers[[1]]$aes_params$alpha = 0.5
p2[[1]]$layers[[1]]$aes_params$alpha = 0.5
p3[[1]]$layers[[1]]$aes_params$alpha = 0.5
p1 <- p1 + NoLegend() 
p1 %>% ggsave(filename = paste(figure_dir, "Fig5D.png", sep=""), width = 150, height = 150, units = "mm")
p2 <- p2 + NoLegend() 
p2 %>% ggsave(filename = paste(figure_dir, "Fig5E.png", sep=""), width = 150, height = 150, units = "mm")
p3 <- p3 + NoLegend() 
p3 %>% ggsave(filename = paste(figure_dir, "Fig5F.png", sep=""), width = 150, height = 150, units = "mm")
rm(gse121862_integrated)




# --------------------------------------------------
# 19/01/2022 - Figures pour concours Amandine
cytokine_list <- c("BTC", "CCL2", "CCL5", "CCL8", "CCL17", "CCL19", "CX3CL1", "CXCL1", "CXCL9", "CXCL10", "CXCL12", "CXCL14", "CXCL16", "IL1RN", "IL33", "IL34", "LCN2", "LGALS9", "LIF", "PDGFD", "PTN", "TIMP1", "TNFSF8")
receptor_list <- c("CCR5", "CCR1", "CCR2", "CCR7", "CX3CR1", "CXCR2", "CXCR3", "CXCR4", "CXCR6", "GPR35", "IL1R1", "IL1RL1", "CSF1R", "HAVCR2")
p1 <- DotPlot(sc_gsm_integrated, features = c(cytokine_list, receptor_list), group.by = "annot_clusters", assay = "SCT", cols = c("blue", "red")) + RotatedAxis()
p1 %>% ggsave(filename = paste(figure_dir, "x_dotplot_concours_av.png", sep = ""), width = 400, height = 225, units = "mm")
p1 <- DoHeatmap(subset(sc_gsm_integrated, downsample = 200), features = c(cytokine_list, receptor_list), group.by = "annot_clusters", group.bar = T, assay = "SCT", slot = "data") + NoLegend()
p1 %>% ggsave(filename = paste(figure_dir, "x_heatmap_concours_av.png", sep=""), width = 450, height = 225, units = "mm")
