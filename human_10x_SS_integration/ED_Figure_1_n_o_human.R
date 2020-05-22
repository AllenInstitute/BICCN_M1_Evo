library(Seurat)
library(Matrix)
library(matrixStats)
library(dplyr)
library(feather)
library(scrattch.hicat)
library(scrattch.io)

#load raw UMIs, create seurat object, and load with metadata
data <- readRDS("~/NeMO_analysis_folder/Transcriptomics/human_10x_SS_integration/Figure_ED1_source_data/human_data.RDS")
anno <- read_feather("~/NeMO_analysis_folder/Transcriptomics/human_10x_SS_integration/Figure_ED1_source_data/anno.feather")

data <- data[ , which(colnames(data) %in% anno$sample_id)]
anno <- anno[order(match(anno$sample_id, colnames(data))), ]

data <- CreateSeuratObject(counts = data)
data$cluster_label <- anno$cluster_label
data$cluster_id <- anno$cluster_id
data$cluster_color <- anno$cluster_color
data$class_label <- anno$class_label
data$class_id <- anno$class_id
data$class_color <- anno$class_color
data$donor <- anno$donor_id
data$sample_id <- anno$sample_id

#SCT-normalize
de_genes <- read.csv("~/NeMO_analysis_folder/Transcriptomics/human_10x_SS_integration/Figure_ED1_source_data/all_subclass_DE_genes.csv")
de_genes <- de_genes[which(de_genes$human_de), ]                     
de_genes$genes <- as.character(de_genes$genes)
de_genes <- unique(de_genes$genes)
de_genes <- de_genes[which(de_genes %in% rownames(data))]
VariableFeatures(data) <- de_genes

data <- SCTransform(data, verbose = TRUE, return.only.var.genes = FALSE) 
data <- RunPCA(data, features = de_genes)
data <- RunTSNE(data, dims = 1:50)                     

library(ggplot2)                     
Idents(data) <- data$cluster_label
order.of.clusters <- levels(data@active.ident)
tsne_colors <- data.frame(clusters = order.of.clusters)
tsne_colors$colors <- data$cluster_color[match(tsne_colors$clusters, data$cluster_label)]
DimPlot(data, reduction = "tsne", pt.size = 0.7, label = TRUE, label.size = 2) + 
  ggtitle("Clusters") + 
  NoLegend() +
  scale_color_manual(values = tsne_colors$colors)

DimPlot(data, reduction = "tsne", group.by = "donor", pt.size = 1.2)          + ggtitle("Species") 
