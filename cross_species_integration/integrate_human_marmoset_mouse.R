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

#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
############################### Load in premade raw UMI Seurat objects with metadata    #############################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################

human_data <- readRDS("~/NeMO_analysis_folder/Transcriptomics/cross_species_integration/source_data_human/human_seurat_trimmed_genes.RDS")
marmoset_data <- readRDS("~/NeMO_analysis_folder/Transcriptomics/cross_species_integration/source_data_marmoset/marmoset_seurat_trimmed_genes.RDS")
mouse_data <- readRDS("~/NeMO_analysis_folder/Transcriptomics/cross_species_integration/source_data_mouse/mouse_seurat_trimmed_genes.RDS")
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
############################### Filter only cell class of interest to integrate         #############################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################

Idents(human_data) <- human_data$class_label
Idents(marmoset_data) <- marmoset_data$class_label
Idents(mouse_data) <- mouse_data$class_label
human_data <- subset(human_data, idents = "Glutamatergic")         #change to GABAergic or Non-Neuronal to integrate different class
marmoset_data <- subset(marmoset_data, idents = "Glutamatergic")   #change to GABAergic or Non-Neuronal to integrate different class
mouse_data <- subset(mouse_data, idents = "Glutamatergic")         #change to GABAergic or Non-Neuronal to integrate different class

#weighted average downsample... 200 cells max per pre-integrated cluster for species with most subclasses
human_subclass <- table(human_data$subclass_label, human_data$cluster_label)
marmoset_subclass <- table(marmoset_data$subclass_label, marmoset_data$cluster_label)
mouse_subclass <- table(mouse_data$subclass_label, mouse_data$cluster_label)

all_subclasses <- data.frame(subclass = unique(c(rownames(human_subclass), rownames(marmoset_subclass), rownames(mouse_subclass))))
all_subclasses$human <- 0
all_subclasses$marmoset <- 0
all_subclasses$mouse <- 0

for(i in 1:nrow(all_subclasses)){ #calculate how many clusters per subclass
  all_subclasses$human[i] <- length(which(human_subclass[which(rownames(human_subclass) == all_subclasses$subclass[i]), ] > 0))
  all_subclasses$marmoset[i] <- length(which(marmoset_subclass[which(rownames(marmoset_subclass) == all_subclasses$subclass[i]), ] > 0))
  all_subclasses$mouse[i] <- length(which(mouse_subclass[which(rownames(mouse_subclass) == all_subclasses$subclass[i]), ] > 0))
}
all_subclasses

all_subclasses$max_cells <- 0
for(i in 1:nrow(all_subclasses)){ #calculate max cells per subclass
  all_subclasses$max_cells[i] <- all_subclasses[i, 1 + which.max(all_subclasses[i, 2:4])] * 200
}
all_subclasses

all_subclasses$human_cells_per_cl <- 0
all_subclasses$marmoset_cells_per_cl <- 0
all_subclasses$mouse_cells_per_cl <- 0
for(i in 1:nrow(all_subclasses)){
  all_subclasses$human_cells_per_cl[i] <- round(all_subclasses$max_cells[i] / all_subclasses$human[i])
  all_subclasses$marmoset_cells_per_cl[i] <- round(all_subclasses$max_cells[i] / all_subclasses$marmoset[i])
  all_subclasses$mouse_cells_per_cl[i] <- round(all_subclasses$max_cells[i] / all_subclasses$mouse[i])
}
all_subclasses

#subset human data
cells_to_keep <- NA
for(i in 1:nrow(all_subclasses)){
  clusters_to_subset <- unique(human_data$cluster_label[which(human_data$subclass_label %in% all_subclasses$subclass[i])])
  for(p in 1:length(clusters_to_subset)){
    if(length(which(human_data$cluster_label == clusters_to_subset[p])) > all_subclasses$human_cells_per_cl[i]){
      cells_to_keep <- c(cells_to_keep, sample(names(which(human_data$cluster_label == clusters_to_subset[p])), all_subclasses$human_cells_per_cl[i]))
    }
    if(length(which(human_data$cluster_label == clusters_to_subset[p])) <= all_subclasses$human_cells_per_cl[i]){
      cells_to_keep <- c(cells_to_keep, names(which(human_data$cluster_label == clusters_to_subset[p])))
    }
  }
}
cells_to_keep <- cells_to_keep[-1]
human_data$sample_id <- colnames(human_data)
Idents(human_data) <- human_data$sample_id
human_data <- subset(human_data, idents = cells_to_keep)

#subset marmoset data
cells_to_keep <- NA
for(i in 1:nrow(all_subclasses)){
  clusters_to_subset <- unique(marmoset_data$cluster_label[which(marmoset_data$subclass_label %in% all_subclasses$subclass[i])])
  for(p in 1:length(clusters_to_subset)){
    if(length(which(marmoset_data$cluster_label == clusters_to_subset[p])) > all_subclasses$marmoset_cells_per_cl[i]){
      cells_to_keep <- c(cells_to_keep, sample(names(which(marmoset_data$cluster_label == clusters_to_subset[p])), all_subclasses$marmoset_cells_per_cl[i]))
    }
    if(length(which(marmoset_data$cluster_label == clusters_to_subset[p])) <= all_subclasses$marmoset_cells_per_cl[i]){
      cells_to_keep <- c(cells_to_keep, names(which(marmoset_data$cluster_label == clusters_to_subset[p])))
    }
  }
}
cells_to_keep <- cells_to_keep[-1]
marmoset_data$sample_id <- colnames(marmoset_data)
Idents(marmoset_data) <- marmoset_data$sample_id
marmoset_data <- subset(marmoset_data, idents = cells_to_keep)

#subset mouse data
cells_to_keep <- NA
for(i in 1:nrow(all_subclasses)){
  clusters_to_subset <- unique(mouse_data$cluster_label[which(mouse_data$subclass_label %in% all_subclasses$subclass[i])])
  for(p in 1:length(clusters_to_subset)){
    if(length(which(mouse_data$cluster_label == clusters_to_subset[p])) > all_subclasses$mouse_cells_per_cl[i]){
      cells_to_keep <- c(cells_to_keep, sample(names(which(mouse_data$cluster_label == clusters_to_subset[p])), all_subclasses$mouse_cells_per_cl[i]))
    }
    if(length(which(mouse_data$cluster_label == clusters_to_subset[p])) <= all_subclasses$mouse_cells_per_cl[i]){
      cells_to_keep <- c(cells_to_keep, names(which(mouse_data$cluster_label == clusters_to_subset[p])))
    }
  }
}
cells_to_keep <- cells_to_keep[-1]
mouse_data$sample_id <- colnames(mouse_data)
Idents(mouse_data) <- mouse_data$sample_id
mouse_data <- subset(mouse_data, idents = cells_to_keep)
table(mouse_data$cluster_label)


#####################################################################################################################
###############################       SCT norm and integration                        ###############################
#####################################################################################################################

human_data$orig.ident <- "human"
marmoset_data$orig.ident <- "marmoset"
mouse_data$orig.ident <- "mouse"

all.data <- merge(x = human_data, y = list(marmoset_data, mouse_data), add.cell.ids = c("human", "marmoset", "mouse"))

combined.list <- SplitObject(all.data, split.by = "orig.ident")

for (i in 1:length(combined.list)) {
  combined.list[[i]] <- SCTransform(combined.list[[i]],  verbose = TRUE)
}

#Find DE genes to guide alignment

Var.genes.human         <- select_markers(combined.list$human@assays$SCT@counts, combined.list$human$cluster_id, n.markers = 100)
Var.genes.human.markers <- Var.genes.human$markers
Var.genes.marmoset         <- select_markers(combined.list$marmoset@assays$SCT@counts, combined.list$marmoset$cluster_id, n.markers = 100)
Var.genes.marmoset.markers <- Var.genes.marmoset$markers
Var.genes.mouse         <- select_markers(combined.list$mouse@assays$SCT@counts, combined.list$mouse$cluster_id, n.markers = 100)
Var.genes.mouse.markers <- Var.genes.mouse$markers
total.Var.genes <- combine(Var.genes.human.markers, Var.genes.marmoset.markers, Var.genes.mouse.markers)
total.Var.genes <- unique(total.Var.genes$data)

total.Var.genes <- total.Var.genes[which(total.Var.genes %in% rownames(combined.list$human@assays$SCT@counts))]
total.Var.genes <- total.Var.genes[which(total.Var.genes %in% rownames(combined.list$marmoset@assays$SCT@counts))]
total.Var.genes <- total.Var.genes[which(total.Var.genes %in% rownames(combined.list$mouse@assays$SCT@counts))]


#Need to increase memory for server usage
library(future)
options(future.globals.maxSize = 10000 * 1024^2)

combined.list <- PrepSCTIntegration(object.list = combined.list, anchor.features = total.Var.genes,
                                    verbose = FALSE)

combined.anchors <- FindIntegrationAnchors(object.list = combined.list, normalization.method = "SCT",
                                           anchor.features = total.Var.genes, verbose = TRUE)

sample.combined <- IntegrateData(anchorset = combined.anchors, normalization.method = "SCT", 
                                 verbose = TRUE)

sample.combined <- RunPCA(sample.combined, features = total.Var.genes, npcs = 100)
ElbowPlot(sample.combined, ndims = 100)
sample.combined <- FindNeighbors(sample.combined, reduction = "pca", dims = 1:100, nn.eps = 0)
sample.combined <- FindClusters(sample.combined, resolution = 6, n.start = 100) 
sample.combined$seurat_clusters.new <- as.integer(sample.combined$seurat_clusters)
sample.combined <- RunTSNE(sample.combined, dims = 1:100)


