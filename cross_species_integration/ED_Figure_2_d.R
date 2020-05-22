library(Seurat)
library(dplyr)
library(Matrix)
library(matrixStats)
library(gplots)
library(ggplot2)
library(feather)
library(viridis)

#load data and format
anno <- read_feather("~/NeMO_analysis_folder/Transcriptomics/cross_species_integration/additional_intermediate_files/3_species_anno.feather")
sample.combined <- readRDS("~/NeMO_analysis_folder/Transcriptomics/cross_species_integration/sample.combined_inh_integration.RDS")

anno$cross_species_cluster_label <- sub("human_", "", anno$cross_species_cluster_label)
anno$cross_species_cluster_label <- sub("marmoset_", "", anno$cross_species_cluster_label)
anno$cross_species_cluster_label <- sub("mouse_", "", anno$cross_species_cluster_label)

sample.combined$new_sample_id <- paste(sample.combined$orig.ident, sample.combined$sample_id, sep = "_")
sample.combined$cross_species_cluster_label <- anno$cross_species_cluster_label[match(sample.combined$new_sample_id, anno$sample_id)]

#exclude nuclei without an integrated cluster association
sample.combined$keep <- "YES"
sample.combined$keep[which(sample.combined$cross_species_cluster_label == "exclude")] <- "NO"
Idents(sample.combined) <- sample.combined$keep
sample.combined <- subset(sample.combined, idents = "YES")

#downsample to 20 nuclei per cluster
data.subset <- sample.combined
Idents(data.subset) <- data.subset$cluster_label
data.subset <- subset(data.subset, downsample = 20)

#subset each species and find DE genes for each species
Idents(data.subset) <- data.subset$orig.ident
data.human <- subset(data.subset, idents = "human")
data.marmoset <- subset(data.subset, idents = "marmoset")
data.mouse <- subset(data.subset, idents = "mouse")
Idents(data.human) <- data.human$cross_species_cluster_label
Idents(data.marmoset) <- data.marmoset$cross_species_cluster_label
Idents(data.mouse) <- data.mouse$cross_species_cluster_label

human.markers <- FindAllMarkers(data.human, assay = "SCT", slot = "counts", only.pos = TRUE)
marmoset.markers <- FindAllMarkers(data.marmoset, assay = "SCT", slot = "counts", only.pos = TRUE)
mouse.markers <- FindAllMarkers(data.mouse, assay = "SCT", slot = "counts", only.pos = TRUE)

#find genes that are DE in same in cl across all 3 species
int_cl <- as.character(unique(human.markers$cluster))
genes_use <- data.frame(matrix(NA, nrow = 1, ncol = ncol(human.markers)))
colnames(genes_use) <- colnames(human.markers)

for(i in 1:length(int_cl)){
  tmp1 <- human.markers[which(human.markers$cluster == int_cl[i]), ]
  tmp2 <- marmoset.markers[which(marmoset.markers$cluster == int_cl[i]), ]
  tmp3 <- mouse.markers[which(mouse.markers$cluster == int_cl[i]), ]
  
  tmp1 <- tmp1[which(tmp1$gene %in% tmp2$gene), ]
  tmp1 <- tmp1[which(tmp1$gene %in% tmp3$gene), ]
  genes_use <- rbind(tmp1, genes_use)
}

#order genes and prepare data order for heatmap
order <- c("Lamp5_1", "Lamp5_2", "Lamp5_3", "Lamp5_4", "Sncg_1", "Sncg_2", "Sncg_3",
           "Lamp5_5", "Sncg_4", "Vip_1", "Vip_2", "Vip_3", "Vip_4", "Sst Chodl", 
           "Sst_1", "Sst_2","Sst_3","Sst_4","Sst_5","Sst_6","Sst_7", "Pvalb_1", "Pvalb_2", "Chandelier")
genes_use_all <- genes_use[order(match(genes_use$cluster, order)), ]

################################################################################
#########  find markers for subclasses - Lamp5       ##########################
################################################################################

cl_of_interest <- c("Lamp5_1", "Lamp5_2", "Lamp5_3", "Lamp5_4", "Lamp5_5")
Idents(data.human) <- data.human$cross_species_cluster_label
Idents(data.marmoset) <- data.marmoset$cross_species_cluster_label
Idents(data.mouse) <- data.mouse$cross_species_cluster_label


human.sub <- subset(data.human, idents = cl_of_interest)
marmoset.sub <- subset(data.marmoset, idents = cl_of_interest)
mouse.sub <- subset(data.mouse, idents = cl_of_interest)

human.sub.markers <- FindAllMarkers(human.sub, assay = "SCT", slot = "counts", only.pos = TRUE)
marmoset.sub.markers <- FindAllMarkers(marmoset.sub, assay = "SCT", slot = "counts", only.pos = TRUE)
mouse.sub.markers <- FindAllMarkers(mouse.sub, assay = "SCT", slot = "counts", only.pos = TRUE)

int_cl <- as.character(unique(human.sub.markers$cluster))
genes_use <- data.frame(matrix(NA, nrow = 1, ncol = ncol(human.sub.markers)))
colnames(genes_use) <- colnames(human.sub.markers)

for(i in 1:length(int_cl)){
  tmp1 <- human.sub.markers[which(human.sub.markers$cluster == int_cl[i]), ]
  tmp2 <- marmoset.sub.markers[which(marmoset.sub.markers$cluster == int_cl[i]), ]
  tmp3 <- mouse.sub.markers[which(mouse.sub.markers$cluster == int_cl[i]), ]
  
  tmp1 <- tmp1[which(tmp1$gene %in% tmp2$gene), ]
  tmp1 <- tmp1[which(tmp1$gene %in% tmp3$gene), ]
  genes_use <- rbind(tmp1, genes_use)
}

genes_use_LAMP5 <- genes_use[order(match(genes_use$cluster, order)), ]

################################################################################
#########  find markers for subclasses - Sst        ##########################
################################################################################

cl_of_interest <- c("Sst_1", "Sst_2", "Sst_3", "Sst_4", "Sst_5", "Sst_6", "Sst_7")
Idents(data.human) <- data.human$cross_species_cluster_label
Idents(data.marmoset) <- data.marmoset$cross_species_cluster_label
Idents(data.mouse) <- data.mouse$cross_species_cluster_label


human.sub <- subset(data.human, idents = cl_of_interest)
marmoset.sub <- subset(data.marmoset, idents = cl_of_interest)
mouse.sub <- subset(data.mouse, idents = cl_of_interest)

human.sub.markers <- FindAllMarkers(human.sub, assay = "SCT", slot = "counts", only.pos = TRUE)
marmoset.sub.markers <- FindAllMarkers(marmoset.sub, assay = "SCT", slot = "counts", only.pos = TRUE)
mouse.sub.markers <- FindAllMarkers(mouse.sub, assay = "SCT", slot = "counts", only.pos = TRUE)

int_cl <- as.character(unique(human.sub.markers$cluster))
genes_use <- data.frame(matrix(NA, nrow = 1, ncol = ncol(human.sub.markers)))
colnames(genes_use) <- colnames(human.sub.markers)

for(i in 1:length(int_cl)){
  tmp1 <- human.sub.markers[which(human.sub.markers$cluster == int_cl[i]), ]
  tmp2 <- marmoset.sub.markers[which(marmoset.sub.markers$cluster == int_cl[i]), ]
  tmp3 <- mouse.sub.markers[which(mouse.sub.markers$cluster == int_cl[i]), ]
  
  tmp1 <- tmp1[which(tmp1$gene %in% tmp2$gene), ]
  tmp1 <- tmp1[which(tmp1$gene %in% tmp3$gene), ]
  genes_use <- rbind(tmp1, genes_use)
}

genes_use_Sst <- genes_use[order(match(genes_use$cluster, order)), ]

################################################################################
#########  plot genes                                 ##########################
################################################################################

genes_use_all <- genes_use_all %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
genes_use_all <- as.data.frame(genes_use_all)
genes_use_Sst <- genes_use_Sst %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
genes_use_Sst <- as.data.frame(genes_use_Sst)
genes_use_LAMP5 <- genes_use_LAMP5 %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
genes_use_LAMP5 <- as.data.frame(genes_use_LAMP5)


genes_use <- rbind(genes_use_all, genes_use_Sst, genes_use_LAMP5)
genes_use <- genes_use[order(match(genes_use$cluster, order)), ]

data_use <- data.human
Idents(data_use) <- data_use$cross_species_cluster_label
data_use <- subset(data_use, downsample = 20)
data_use <- ScaleData(data_use, assay = "SCT")
levels(data_use) <- order

DoHeatmap(data_use, assay = "SCT", slot = "scale.data", features = unique(genes_use$gene)) + NoLegend() +
  scale_fill_gradientn(colors = c("white", "white", "#ae5a8c")) +
  theme(axis.text.y = element_text(size = 4))
