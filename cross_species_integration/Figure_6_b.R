library(mgcv)
library(Seurat)
library(cowplot)
library(feather)
library(dplyr)
library(Matrix)
library(Hmisc)
library(tidyr)
library(gplots)
library(gdata)
library(stringr)
library(rlang)
library(scrattch.hicat)
library(scrattch.io)
library(matrixStats)
library(ggplot2)
library(igraph)
library(visNetwork)
library(riverplot)
library(scales)
library(dendextend)
library(reshape)
library(openxlsx)
library(rhdf5)
options(stringsAsFactors = FALSE)

#load integrated inhibitory data
sample.combined <- readRDS("~/NeMO_analysis_folder/Transcriptomics/cross_species_integration/additional_intermediate_files/sample.combined_glia_integration.RDS")

#subset PVALB cells for chandelier and basket
Idents(sample.combined) <- sample.combined$subclass_label
sample.combined <- subset(sample.combined, idents = "Pvalb")

#label basket and chandelier cells
chandelier_clusters <- c(unique(sample.combined$cluster_label[grep("FAM19A4", sample.combined$cluster_label)]),
                         unique(sample.combined$cluster_label[grep("COL15A1", sample.combined$cluster_label)]),
                         unique(sample.combined$cluster_label[grep("Vipr2", sample.combined$cluster_label)]))

sample.combined$PVALB_label <- "Basket"
sample.combined$PVALB_label[which(sample.combined$cluster_label %in% chandelier_clusters)] <- "Chandelier"

#subset by species
Idents(sample.combined) <- sample.combined$orig.ident
human <- subset(sample.combined, idents = "human")
marmoset <- subset(sample.combined, idents = "marmoset")
mouse <- subset(sample.combined, idents = "mouse")

#Find Chandelier specific genes for each species
Idents(human) <- human$PVALB_label
Idents(marmoset) <- marmoset$PVALB_label
Idents(mouse) <- mouse$PVALB_label

human_genes <- FindAllMarkers(human, assay = "SCT", slot = "counts", test.use = "roc")
marmoset_genes <- FindAllMarkers(marmoset, assay = "SCT", slot = "counts", test.use = "roc")
mouse_genes <- FindAllMarkers(mouse, assay = "SCT", slot = "counts", test.use = "roc")

############################################################################################################################
########################## Find DE genes               #####################################################################
############################################################################################################################

human_genes <- human_genes[which(human_genes$avg_diff > 0), ]
marmoset_genes <- marmoset_genes[which(marmoset_genes$avg_diff > 0), ]
mouse_genes <- mouse_genes[which(mouse_genes$avg_diff > 0), ]

human_genes <- human_genes[which(human_genes$cluster == "Chandelier"), ]
marmoset_genes <- marmoset_genes[which(marmoset_genes$cluster == "Chandelier"), ]
mouse_genes <- mouse_genes[which(mouse_genes$cluster == "Chandelier"), ]

#venn Diagrams
all.genes <- data.frame(genes = unique(c(human_genes$gene, marmoset_genes$gene, mouse_genes$gene)))
all.genes$Human <- as.character(match(all.genes$genes, human_genes$gene))
all.genes$Marmoset <- as.character(match(all.genes$genes, marmoset_genes$gene))
all.genes$Mouse <- as.character(match(all.genes$genes, mouse_genes$gene))
all.genes$Human[which(is.na(all.genes$Human))] <- FALSE 
all.genes$Human[which(all.genes$Human != FALSE)] <- TRUE
all.genes$Marmoset[which(is.na(all.genes$Marmoset))] <- FALSE 
all.genes$Marmoset[which(all.genes$Marmoset != FALSE)] <- TRUE
all.genes$Mouse[which(is.na(all.genes$Mouse))] <- FALSE 
all.genes$Mouse[which(all.genes$Mouse != FALSE)] <- TRUE
all.genes$Human <- as.logical(all.genes$Human)
all.genes$Marmoset <- as.logical(all.genes$Marmoset)
all.genes$Mouse <- as.logical(all.genes$Mouse)

library(eulerr)
plot(euler(
  all.genes[ ,2:4]),
  quantities = list(cex = 3),
  labels = NULL,
  main = paste0("Chandelier vs. Basket"),
  fills = c("royalblue1", "maroon4", "sienna2")
)

genes_use <- all.genes[which(all.genes$Human == T & all.genes$Marmoset == T & all.genes$Mouse == T), ]

#conserved genes heatmap

sample.combined$new_subclass <- sample.combined$PVALB_label
sample.combined$new_subclass <- paste(sample.combined$orig.ident, sample.combined$new_subclass, sep = "_")
Idents(sample.combined) <- sample.combined$new_subclass

to_plot <- subset(sample.combined, downsample = 100)
to_plot@assays$SCT@data <- to_plot@assays$SCT@counts
to_plot@assays$SCT@data@x <- log2(to_plot@assays$SCT@data@x + 1)
to_plot <- ScaleData(to_plot, assay = "SCT")
levels(to_plot) <- c("human_Chandelier","human_Basket", "marmoset_Chandelier","marmoset_Basket", "mouse_Chandelier",  "mouse_Basket")
genes_use <- sort(genes_use$genes)
DoHeatmap(to_plot, assay = "SCT", slot = "scale.data", features = genes_use) +
  scale_fill_viridis(option = "magma")
