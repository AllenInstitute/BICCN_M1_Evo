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

human_data <- readRDS("~/NeMO_analysis_folder/Transcriptomics/cross_species_integration/source_data_human/human_loaded_seurat.RDS")
marmoset_data <- readRDS("~/NeMO_analysis_folder/Transcriptomics/cross_species_integration/source_data_marmoset/marmoset_loaded_seurat.RDS")
macaque_data <- readRDS("~/NeMO_analysis_folder/Transcriptomics/cross_species_integration/source_data_macaque/macaque_loaded_seurat.RDS")
mouse_data <- readRDS("~/NeMO_analysis_folder/Transcriptomics/cross_species_integration/source_data_mouse/mouse_loaded_seurat.RDS")

tmp <- unique(macaque$cluster_label)
macaque$subclass_label <- NA

macaque$subclass_label[which(macaque$cluster_label == tmp[1])] <- "L5 ET"
macaque$subclass_label[which(macaque$cluster_label == tmp[2])] <- "L5 IT"
macaque$subclass_label[which(macaque$cluster_label == tmp[3])] <- "L6b"
macaque$subclass_label[which(macaque$cluster_label == tmp[4])] <- "L6 IT Car3"
macaque$subclass_label[which(macaque$cluster_label == tmp[5])] <- "L5 IT"
macaque$subclass_label[which(macaque$cluster_label == tmp[6])] <- "L5/6 NP"
macaque$subclass_label[which(macaque$cluster_label == tmp[7])] <- "L6 CT"
macaque$subclass_label[which(macaque$cluster_label == tmp[8])] <- "L6 IT"
macaque$subclass_label[which(macaque$cluster_label == tmp[9])] <- "L5 IT"
macaque$subclass_label[which(macaque$cluster_label == tmp[10])] <- "L5 ET"
macaque$subclass_label[which(macaque$cluster_label == tmp[11])] <- "L5/6 NP"
macaque$subclass_label[which(macaque$cluster_label == tmp[12])] <- "L6b"
macaque$subclass_label[which(macaque$cluster_label == tmp[13])] <- "L5/6 NP"

#weighted average downsample... 200 cells max per pre-integrated cluster for species with most subclasses
human_subclass <- table(human_data$subclass_label, human_data$cluster_label)
marmoset_subclass <- table(marmoset_data$subclass_label, marmoset_data$cluster_label)
macaque_subclass <- table(macaque_data$subclass_label, macaque_data$cluster_label)
mouse_subclass <- table(mouse_data$subclass_label, mouse_data$cluster_label)

all_subclasses <- data.frame(subclass = unique(c(rownames(human_subclass), rownames(marmoset_subclass), rownames(mouse_subclass), rownames(macaque_subclass))))
all_subclasses$human <- 0
all_subclasses$marmoset <- 0
all_subclasses$mouse <- 0
all_subclasses$macaque <- 0

for(i in 1:nrow(all_subclasses)){ #calculate how many clusters per subclass
  all_subclasses$human[i] <- length(which(human_subclass[which(rownames(human_subclass) == all_subclasses$subclass[i]), ] > 0))
  all_subclasses$marmoset[i] <- length(which(marmoset_subclass[which(rownames(marmoset_subclass) == all_subclasses$subclass[i]), ] > 0))
  all_subclasses$mouse[i] <- length(which(mouse_subclass[which(rownames(mouse_subclass) == all_subclasses$subclass[i]), ] > 0))
  all_subclasses$macaque[i] <- length(which(macaque_subclass[which(rownames(macaque_subclass) == all_subclasses$subclass[i]), ] > 0))
}
all_subclasses

all_subclasses$max_cells <- 0
for(i in 1:nrow(all_subclasses)){ #calculate max cells per subclass
  all_subclasses$max_cells[i] <- all_subclasses[i, 1 + which.max(all_subclasses[i, 2:5])] * 200
}
all_subclasses

all_subclasses$human_cells_per_cl <- 0
all_subclasses$marmoset_cells_per_cl <- 0
all_subclasses$mouse_cells_per_cl <- 0
all_subclasses$macaque_cells_per_cl <- 0
for(i in 1:nrow(all_subclasses)){
  all_subclasses$human_cells_per_cl[i] <- round(all_subclasses$max_cells[i] / all_subclasses$human[i])
  all_subclasses$marmoset_cells_per_cl[i] <- round(all_subclasses$max_cells[i] / all_subclasses$marmoset[i])
  all_subclasses$mouse_cells_per_cl[i] <- round(all_subclasses$max_cells[i] / all_subclasses$mouse[i])
  all_subclasses$macaque_cells_per_cl[i] <- round(all_subclasses$max_cells[i] / all_subclasses$macaque[i])
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
  clusters_to_subset <- unique(macaque_data$cluster_label[which(macaque_data$subclass_label %in% all_subclasses$subclass[i])])
  for(p in 1:length(clusters_to_subset)){
    if(length(which(macaque_data$cluster_label == clusters_to_subset[p])) > all_subclasses$macaque_cells_per_cl[i]){
      cells_to_keep <- c(cells_to_keep, sample(names(which(macaque_data$cluster_label == clusters_to_subset[p])), all_subclasses$macaque_cells_per_cl[i]))
    }
    if(length(which(macaque_data$cluster_label == clusters_to_subset[p])) <= all_subclasses$macaque_cells_per_cl[i]){
      cells_to_keep <- c(cells_to_keep, names(which(macaque_data$cluster_label == clusters_to_subset[p])))
    }
  }
}
cells_to_keep <- cells_to_keep[-1]
macaque_data$sample_id <- colnames(macaque_data)
Idents(macaque_data) <- macaque_data$sample_id
macaque_data <- subset(macaque_data, idents = cells_to_keep)

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
macaque_data$orig.ident <- "macaque"
marmoset_data$orig.ident <- "marmoset"
mouse_data$orig.ident <- "mouse"

all.data <- merge(x = human_data, y = list(macaque_data, marmoset_data, mouse_data), add.cell.ids = c("human","macaque", "marmoset", "mouse"))

combined.list <- SplitObject(all.data, split.by = "orig.ident")

for (i in 1:length(combined.list)) {
  combined.list[[i]] <- SCTransform(combined.list[[i]],  verbose = TRUE)
}

#Find DE genes to guide alignment

Var.genes.human         <- select_markers(combined.list$human@assays$SCT@counts, combined.list$human$cluster_id, n.markers = 100)
Var.genes.human.markers <- Var.genes.human$markers
Var.genes.macaque         <- select_markers(combined.list$macaque@assays$SCT@counts, combined.list$macaque$cluster_id, n.markers = 100)
Var.genes.macaque.markers <- Var.genes.macaque$markers
Var.genes.marmoset         <- select_markers(combined.list$marmoset@assays$SCT@counts, combined.list$marmoset$cluster_id, n.markers = 100)
Var.genes.marmoset.markers <- Var.genes.marmoset$markers
Var.genes.mouse         <- select_markers(combined.list$mouse@assays$SCT@counts, combined.list$mouse$cluster_id, n.markers = 100)
Var.genes.mouse.markers <- Var.genes.mouse$markers
total.Var.genes <- combine(Var.genes.human.markers, Var.genes.macaque.markers, Var.genes.marmoset.markers, Var.genes.mouse.markers)
total.Var.genes <- unique(total.Var.genes$data)

total.Var.genes <- total.Var.genes[which(total.Var.genes %in% rownames(combined.list$human@assays$SCT@counts))]
total.Var.genes <- total.Var.genes[which(total.Var.genes %in% rownames(combined.list$macaque@assays$SCT@counts))]
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
sample.combined <- FindClusters(sample.combined, resolution = 2, n.start = 100) 
sample.combined$seurat_clusters.new <- as.integer(sample.combined$seurat_clusters)
sample.combined <- RunTSNE(sample.combined, dims = 1:100)
