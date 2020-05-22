library(Seurat)
library(dplyr)
library(matrixStats)
library(Matrix)
library(ggplot2)
library(scrattch.hicat)

sample.combined <- readRDS("/NeMO_analysis_folder/Transcriptomics/cross_species_integration/sample.combined_inh_integration.RDS")

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
data_to_plot <- data_to_plot[which(data_to_plot$species == "mouse"), ]
data_to_plot$color <- sample.combined$cluster_color[match(data_to_plot$cluster, sample.combined$cluster_label)]

ggplot(data_to_plot,
       aes(x = UMAP_1, y = UMAP_2, color = color)) +
  geom_point(size = 1.5, color = data_to_plot$color)+
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
  geom_point(size = 1.5, color = data_to_plot$color)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 

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
  geom_point(size = 1.5)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_color_manual(values = c("royalblue1", "maroon4", "sienna2")) 



##############################################################################################################################################
##############################################################################################################################################
############################################# DE genes plots     #############################################################################
##############################################################################################################################################
##############################################################################################################################################

Idents(sample.combined) <- sample.combined$orig.ident
human_data <- subset(sample.combined, idents = "human")
marmoset_data <- subset(sample.combined, idents = "marmoset")
mouse_data <- subset(sample.combined, idents = "mouse")

Idents(human_data) <- human_data$subclass_label
Idents(marmoset_data) <- marmoset_data$subclass_label
Idents(mouse_data) <- mouse_data$subclass_label

#Get rid of Meis2 for comparison... not in human
subclasses <- unique(human_data$subclass_label)
marmoset_data <- subset(marmoset_data, idents = subclasses)
mouse_data <- subset(mouse_data, idents = subclasses)

#subset to run quicker
Idents(human_data) <- human_data$cluster_label
Idents(marmoset_data) <- marmoset_data$cluster_label
Idents(mouse_data) <- mouse_data$cluster_label

human_data <- subset(human_data, downsample = 200)
marmoset_data <- subset(marmoset_data, downsample = 200)
mouse_data <- subset(mouse_data, downsample = 200)

Idents(human_data) <- human_data$subclass_label
Idents(marmoset_data) <- marmoset_data$subclass_label
Idents(mouse_data) <- mouse_data$subclass_label

#Find DE genes
human_cells_markers <- FindAllMarkers(human_data, assay = "SCT", slot = "data", test.use = "roc")
marmoset_cells_markers <- FindAllMarkers(marmoset_data, assay = "SCT", slot = "data", test.use = "roc")
mouse_cells_markers <- FindAllMarkers(mouse_data, assay = "SCT", slot = "data", test.use = "roc")

subclasses
tmp <- 1 #change this value to generate different subclass Venn diagrams
human_genes <- human_cells_markers[grep(subclasses[tmp], human_cells_markers$cluster), ]
marmoset_genes <- marmoset_cells_markers[grep(subclasses[tmp], marmoset_cells_markers$cluster), ]
mouse_genes <- mouse_cells_markers[grep(subclasses[tmp], mouse_cells_markers$cluster), ]
human_genes <- human_genes[which(human_genes$avg_diff > 0), ]
marmoset_genes <- marmoset_genes[which(marmoset_genes$avg_diff > 0), ]
mouse_genes <- mouse_genes[which(mouse_genes$avg_diff > 0), ]

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
  main = paste0(subclasses[tmp], " vs. All Inh"),
  fills = c("royalblue1", "maroon4", "sienna2")
)

subclasses <- c("Lamp5", "Sncg", "Vip", "Sst Chodl", "Sst", "Pvalb")

#heatmap - conserved genes
genes_to_plot <- NA
for(i in 1:length(subclasses)){
  tmp <- i
  human_genes <- human_cells_markers[grep(subclasses[tmp], human_cells_markers$cluster), ]
  marmoset_genes <- marmoset_cells_markers[grep(subclasses[tmp], marmoset_cells_markers$cluster), ]
  mouse_genes <- mouse_cells_markers[grep(subclasses[tmp], mouse_cells_markers$cluster), ]
  human_genes <- human_genes[which(human_genes$avg_diff > 0), ]
  marmoset_genes <- marmoset_genes[which(marmoset_genes$avg_diff > 0), ]
  mouse_genes <- mouse_genes[which(mouse_genes$avg_diff > 0), ]
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
  all.genes <- all.genes[which(all.genes[,2] == TRUE), ]
  all.genes <- all.genes[which(all.genes[,3] == TRUE), ]
  all.genes <- all.genes[which(all.genes[,4] == TRUE), ]
  gc()
  genes_to_plot <- c(genes_to_plot, as.character(all.genes$genes))
}
genes_to_plot1 <- genes_to_plot[-1]


data_to_plot <- subset(mouse_data, downsample = 50)
data_to_plot <- ScaleData(data_to_plot, assay = "SCT")
levels(data_to_plot) <- subclasses
library(viridis)
DoHeatmap(data_to_plot, assay = "SCT", slot = "scale.data", features = unique(genes_to_plot)) +
  scale_fill_viridis(option = "viridis")

#heatmap - human_specific genes
genes_to_plot <- NA
for(i in 1:length(subclasses)){
  tmp <- i
  human_genes <- human_cells_markers[grep(subclasses[tmp], human_cells_markers$cluster), ]
  marmoset_genes <- marmoset_cells_markers[grep(subclasses[tmp], marmoset_cells_markers$cluster), ]
  mouse_genes <- mouse_cells_markers[grep(subclasses[tmp], mouse_cells_markers$cluster), ]
  human_genes <- human_genes[which(human_genes$avg_diff > 0), ]
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
  all.genes <- all.genes[which(all.genes[,2] == TRUE), ]
  all.genes <- all.genes[which(all.genes[,3] == FALSE), ]
  all.genes <- all.genes[which(all.genes[,4] == FALSE), ]
  gc()
  genes_to_plot <- c(genes_to_plot, as.character(all.genes$genes))
}
genes_to_plot2 <- genes_to_plot[-1]

#heatmap - marmoset_specific genes
genes_to_plot <- NA
for(i in 1:length(subclasses)){
  tmp <- i
  human_genes <- human_cells_markers[grep(subclasses[tmp], human_cells_markers$cluster), ]
  marmoset_genes <- marmoset_cells_markers[grep(subclasses[tmp], marmoset_cells_markers$cluster), ]
  mouse_genes <- mouse_cells_markers[grep(subclasses[tmp], mouse_cells_markers$cluster), ]
  marmoset_genes <- marmoset_genes[which(marmoset_genes$avg_diff > 0), ]
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
  all.genes <- all.genes[which(all.genes[,2] == FALSE), ]
  all.genes <- all.genes[which(all.genes[,3] == TRUE), ]
  all.genes <- all.genes[which(all.genes[,4] == FALSE), ]
  gc()
  genes_to_plot <- c(genes_to_plot, as.character(all.genes$genes))
}
genes_to_plot3 <- genes_to_plot[-1]

#heatmap - mouse_specific genes
genes_to_plot <- NA
for(i in 1:length(subclasses)){
  tmp <- i
  human_genes <- human_cells_markers[grep(subclasses[tmp], human_cells_markers$cluster), ]
  marmoset_genes <- marmoset_cells_markers[grep(subclasses[tmp], marmoset_cells_markers$cluster), ]
  mouse_genes <- mouse_cells_markers[grep(subclasses[tmp], mouse_cells_markers$cluster), ]
  mouse_genes <- mouse_genes[which(mouse_genes$avg_diff > 0), ]
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
  all.genes <- all.genes[which(all.genes[,2] == FALSE), ]
  all.genes <- all.genes[which(all.genes[,3] == FALSE), ]
  all.genes <- all.genes[which(all.genes[,4] == TRUE), ]
  gc()
  genes_to_plot <- c(genes_to_plot, as.character(all.genes$genes))
}
genes_to_plot4 <- genes_to_plot[-1]

genes_to_plot <- c(genes_to_plot1, genes_to_plot2, genes_to_plot3, genes_to_plot4)

data_to_plot <- subset(marmoset_data, downsample = 50)
data_to_plot <- ScaleData(data_to_plot, assay = "SCT")
levels(data_to_plot) <- subclasses

DoHeatmap(data_to_plot, assay = "SCT", slot = "scale.data", features = unique(genes_to_plot)) +
  scale_fill_gradientn(colors = c("white", "white", "#ae5a8c")) +
  theme(axis.text.y = element_text(size = 2))

#"royalblue1", "maroon4", "sienna2", #ae5a8c
