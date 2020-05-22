library(Seurat)
library(Matrix)
library(matrixStats)
library(dplyr)
library(ggplot2)
library(feather)
library(scrattch.hicat)

#Load pre-made seurat objects
sample.1.data <- readRDS("~/NeMO_analysis_folder/Transcriptomics/human_m1_mtg_integration/source_data_m1/m1_seurat_object.RDS")
sample.2.data <- readRDS("~/NeMO_analysis_folder/Transcriptomics/human_m1_mtg_integration/source_data_mtg/mtg_seurat_object.RDS")

#select only class of interest for integration
Idents(sample.1.data) <- sample.1.data$class
Idents(sample.2.data) <- sample.2.data$class
sample.1.data <- subset(sample.1.data, idents = "exc") #change to inh for glia to integrate different class
sample.2.data <- subset(sample.2.data, idents = "exc") #change to inh for glia to integrate different class


#add metadata to know which region cells came from
sample.1.data$orig.ident <- "M1"
sample.2.data$orig.ident <- "MTG"

sample.1.data <- NormalizeData(sample.1.data)
sample.2.data <- NormalizeData(sample.2.data)

#downsample to max of 300 cells per cluster
Idents(sample.1.data) <- sample.1.data$cluster_label
Idents(sample.2.data) <- sample.2.data$cluster_label
sample.1.data.ds <- subset(sample.1.data, downsample = 300)
sample.2.data.ds <- subset(sample.2.data, downsample = 300)

#Find DE genes to guide alignment
Var.genes.sample.1         <- select_markers(sample.1.data.ds@assays$RNA@data, sample.1.data.ds$cluster_id, n.markers = 100)
Var.genes.sample.1.markers <- Var.genes.sample.1$markers 
Var.genes.sample.2         <- select_markers(sample.2.data.ds@assays$RNA@data, sample.2.data.ds$cluster_id, n.markers = 100)
Var.genes.sample.2.markers <- Var.genes.sample.2$markers 
total.Var.genes <- combine(Var.genes.sample.1.markers, Var.genes.sample.2.markers)
total.Var.genes <- unique(total.Var.genes)

#Integrate 
all.data <- merge(x = sample.1.data, y = sample.2.data, add.cell.ids = c("M1", "MTG"))
combined.list <- SplitObject(all.data, split.by = "orig.ident")

VariableFeatures(combined.list$M1) <- total.Var.genes
VariableFeatures(combined.list$MTG) <- total.Var.genes

combined.anchors <- FindIntegrationAnchors(object.list = combined.list, dims = 1:100, verbose = TRUE)
sample.combined <- IntegrateData(anchorset = combined.anchors, features = total.Var.genes, dims = 1:100, verbose = TRUE)
sample.combined <- ScaleData(sample.combined)
sample.combined <- RunPCA(sample.combined, npcs = 100)
sample.combined <- FindNeighbors(sample.combined, reduction = "pca", dims = 1:100, nn.eps = 0)
sample.combined <- FindClusters(sample.combined, resolution = 6, n.start = 100) #change resolution to 3 for glia, 6 for neuron classes
sample.combined$seurat_clusters.new <- as.integer(sample.combined$seurat_clusters)
sample.combined <- RunTSNE(sample.combined, dims = 1:100)
