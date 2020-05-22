library(Seurat)
library(dplyr)
library(Matrix)
library(matrixStats)
library(gplots)
library(ggplot2)
library(feather)


sample.combined <- readRDS("~/NeMO_analysis_folder/Transcriptomics/cross_species_integration/sample.combined_exc_4_species_integration.RDS")

Idents(sample.combined) <- sample.combined$orig.ident
human_data <- subset(sample.combined, idents = "human")
marmoset_data <- subset(sample.combined, idents = "marmoset")
mouse_data <- subset(sample.combined, idents = "mouse")
macaque_data <- subset(sample.combined, idents = "macaque")

Idents(human_data) <- human_data$subclass_label
Idents(marmoset_data) <- marmoset_data$subclass_label
Idents(mouse_data) <- mouse_data$subclass_label
Idents(macaque_data) <- macaque_data$subclass_label

subclasses <- c("L5 ET", "L5 IT")
marmoset_data <- subset(marmoset_data, idents = subclasses)
macaque_data <- subset(macaque_data, idents = subclasses)
mouse_data <- subset(mouse_data, idents = subclasses)
human_data <- subset(human_data, idents = subclasses)

#subset to run quicker
Idents(human_data) <- human_data$cluster_label
Idents(marmoset_data) <- marmoset_data$cluster_label
Idents(mouse_data) <- mouse_data$cluster_label
Idents(macaque_data) <- macaque_data$cluster_label


human_data <- subset(human_data, downsample = 200)
marmoset_data <- subset(marmoset_data, downsample = 200)
mouse_data <- subset(mouse_data, downsample = 200)
macaque_data <- subset(macaque_data, downsample = 200)


Idents(human_data) <- human_data$subclass_label
Idents(marmoset_data) <- marmoset_data$subclass_label
Idents(mouse_data) <- mouse_data$subclass_label
Idents(macaque_data) <- macaque_data$subclass_label


#Find DE genes
human_cells_markers <- FindAllMarkers(human_data, assay = "SCT", slot = "data", test.use = "roc")
marmoset_cells_markers <- FindAllMarkers(marmoset_data, assay = "SCT", slot = "data", test.use = "roc")
mouse_cells_markers <- FindAllMarkers(mouse_data, assay = "SCT", slot = "data", test.use = "roc")
macaque_cells_markers <- FindAllMarkers(macaque_data, assay = "SCT", slot = "data", test.use = "roc")

human <- human_cells_markers[which(human_cells_markers$avg_diff >0 ), ]
marmoset <- marmoset_cells_markers[which(marmoset_cells_markers$avg_diff >0 ), ]
macaque <- macaque_cells_markers[which(macaque_cells_markers$avg_diff >0 ), ]
mouse <- mouse_cells_markers[which(mouse_cells_markers$avg_diff >0 ), ]

subclass_of_interest <- "L5 ET"
human <- human[which(human$cluster == subclass_of_interest), ]
macaque <- macaque[which(macaque$cluster == subclass_of_interest), ]
marmoset <- marmoset[which(marmoset$cluster == subclass_of_interest), ]
mouse <- mouse[which(mouse$cluster == subclass_of_interest), ]

unique_genes <- unique(c(human$gene, macaque$gene, marmoset$gene, mouse$gene))
all_genes <- data.frame(ID = seq(1, length(unique_genes)))
all_genes$gene <- unique_genes
all_genes$Human <- 0
all_genes$Macaque <- 0
all_genes$Marmoset <- 0
all_genes$Mouse <- 0

for(i in 1:nrow(all_genes)){
  if(all_genes$gene[i] %in% human$gene){
    all_genes$Human[i] <- 1
  }
  if(all_genes$gene[i] %in% macaque$gene){
    all_genes$Macaque[i] <- 1
  }
  if(all_genes$gene[i] %in% marmoset$gene){
    all_genes$Marmoset[i] <- 1
  }
  if(all_genes$gene[i] %in% mouse$gene){
    all_genes$Mouse[i] <- 1
  }
}
library(UpSetR)
upset(all_genes, 
      order.by = "freq", sets = c("Mouse", "Marmoset", "Macaque", "Human"),
      keep.order= TRUE, text.scale = 1.5, show.numbers = FALSE, 
      query.legend = "none", number.angles = 30, point.size = 3.5, line.size = 2,
      mainbar.y.label = "Shared DE Genes", sets.x.label = "Total DE genes \n (L5 ET vs. other Exc)",
      queries = list(
        list(
          query = intersects,
          params = list("Human"),
          active = T,
          color = "royalblue1"
        ),
        list(
          query = intersects,
          params = list("Marmoset"),
          active = T,
          color = "maroon4"
        ),
        list(
          query = intersects,
          params = list("Mouse"),
          active = T,
          color = "sienna2"
        ),
        list(
          query = intersects,
          params = list("Macaque"),
          active = T,
          color = "green4"
        )
      )
)

