library(Seurat)
library(dplyr)
library(matrixStats)
library(Matrix)
library(ggplot2)
library(scrattch.hicat)

sample.combined <- readRDS("~/NeMO_analysis_folder/Transcriptomics/cross_species_integration/sample.combined_glia_integration.RDS")

#exclude LRP4 (QC'd as doublet)
Idents(sample.combined) <- sample.combined$cluster_label
exclude <- unique(sample.combined$cluster_label[grep("LRP4", sample.combined$cluster_label)])
keep <- unique(sample.combined$cluster_label[which(sample.combined$cluster_label %in% exclude == FALSE)])
sample.combined <- subset(sample.combined, idents = keep)
sample.combined <- RunUMAP(sample.combined, dims = 1:100)

#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
############################### Visualization                                        ################################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################

#UMAP cluster
data_to_plot <- data.frame(sample.combined@reductions$umap@cell.embeddings)
sample.combined$new_id <- colnames(sample.combined)
data_to_plot$cluster <- sample.combined$cluster_label[match(rownames(data_to_plot), sample.combined$new_id)]
data_to_plot$species <- sample.combined$orig.ident[match(rownames(data_to_plot), sample.combined$new_id)]
data_to_plot$plot_order <- 0
for(i in 1:nrow(data_to_plot)){
  data_to_plot$plot_order[i] <- sample(1:100000, 1)
  print(i)
}
gc()
data_to_plot <- arrange(data_to_plot, plot_order)
data_to_plot <- data_to_plot[which(data_to_plot$species == "marmoset"), ]
data_to_plot$color <- sample.combined$cluster_color[match(data_to_plot$cluster, sample.combined$cluster_label)]

ggplot(data_to_plot,
       aes(x = UMAP_1, y = UMAP_2, color = color)) +
  geom_point(size = 2.5, color = data_to_plot$color)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 

#UMAP subclass
data_to_plot <- data.frame(sample.combined@reductions$umap@cell.embeddings)
sample.combined$new_id <- colnames(sample.combined)
data_to_plot$cluster <- sample.combined$subclass_label[match(rownames(data_to_plot), sample.combined$new_id)]
data_to_plot$species <- sample.combined$orig.ident[match(rownames(data_to_plot), sample.combined$new_id)]
data_to_plot$plot_order <- 0
for(i in 1:nrow(data_to_plot)){
  data_to_plot$plot_order[i] <- sample(1:100000, 1)
  print(i)
}
gc()
data_to_plot <- arrange(data_to_plot, plot_order)
data_to_plot$color <- sample.combined$subclass_color[match(data_to_plot$cluster, sample.combined$subclass_label)]

ggplot(data_to_plot,
       aes(x = UMAP_1, y = UMAP_2, color = color)) +
  geom_point(size = 2.5, color = data_to_plot$color)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 

DimPlot(sample.combined, group.by = "subclass_label", reduction = "umap")
DimPlot(sample.combined, group.by = "cluster_label", reduction = "umap", split.by = "orig.ident", label = TRUE, label.size = 3) + NoLegend()

#UMAP species
data_to_plot <- data.frame(sample.combined@reductions$umap@cell.embeddings)
sample.combined$new_id <- colnames(sample.combined)
data_to_plot$species <- sample.combined$orig.ident[match(rownames(data_to_plot), sample.combined$new_id)]
data_to_plot$plot_order <- 0
for(i in 1:nrow(data_to_plot)){
  data_to_plot$plot_order[i] <- sample(1:100000, 1)
  print(i)
}
data_to_plot <- arrange(data_to_plot, plot_order)

ggplot(data_to_plot,
       aes(x = UMAP_1, y = UMAP_2, color = species)) +
  geom_point(size = 2.5)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_color_manual(values = c("royalblue1", "maroon4", "sienna2")) 

#plot mouse oligo maturation genes
Idents(sample.combined) <- sample.combined$orig.ident
mouse <- subset(sample.combined, idents = "mouse")
Idents(mouse) <- mouse$cluster_label

oligo_mouse <- c("Oligo", "OPC")
Idents(mouse) <- mouse$subclass_label
oligo_mouse <- subset(mouse, idents = oligo_mouse)

FeaturePlot(oligo_mouse,  features = c("PDGFRA"), reduction = "umap", min.cutoff = 0, pt.size = 1, cols = c("cyan", "red"))
FeaturePlot(oligo_mouse,  features = c("LHFPL3"), reduction = "umap", min.cutoff = 0, pt.size = 1, cols = c("cyan", "red"))
FeaturePlot(oligo_mouse,  features = c("TMEM2"), reduction = "umap", min.cutoff = 0, pt.size = 1, cols = c("cyan", "red"))
FeaturePlot(oligo_mouse,  features = c("ARAP2"), reduction = "umap", min.cutoff = 0, pt.size = 1, cols = c("cyan", "red"))
FeaturePlot(oligo_mouse,  features = c("OPALIN"), reduction = "umap", min.cutoff = 0, pt.size = 1, cols = c("cyan", "red"))
FeaturePlot(oligo_mouse,  features = c("MBP"), reduction = "umap", min.cutoff = 0, pt.size = 1, cols = c("cyan", "red"))



