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
############################### Load in premade raw UMI Seurat objects with metadata    #############################
#####################################################################################################################

human_data <- readRDS("http://data.nemoarchive.org/biccn/lab/lein/2020_M1_study_analysis/Transcriptomics/cross_species_integration/source_data_human/human_seurat_trimmed_genes.RDS")
marmoset_data <- readRDS("http://data.nemoarchive.org/biccn/lab/lein/2020_M1_study_analysis/Transcriptomics/cross_species_integration/source_data_marmoset/marmoset_seurat_trimmed_genes.RDS")

#####################################################################################################################
############################### Filter only cell class of interest to integrate         #############################
#####################################################################################################################

Idents(human_data) <- human_data$class_label
Idents(marmoset_data) <- marmoset_data$class_label
human_data <- subset(human_data, idents = "Glutamatergic")
marmoset_data <- subset(marmoset_data, idents = "Glutamatergic")


#weighted average downsample... 
human_subclass <- table(human_data$subclass_label, human_data$cluster_label)
marmoset_subclass <- table(marmoset_data$subclass_label, marmoset_data$cluster_label)


all_subclasses <- data.frame(subclass = unique(c(rownames(human_subclass), rownames(marmoset_subclass))))
all_subclasses$human <- 0
all_subclasses$marmoset <- 0


for(i in 1:nrow(all_subclasses)){ #calculate how many clusters per subclass
  all_subclasses$human[i] <- length(which(human_subclass[which(rownames(human_subclass) == all_subclasses$subclass[i]), ] > 0))
  all_subclasses$marmoset[i] <- length(which(marmoset_subclass[which(rownames(marmoset_subclass) == all_subclasses$subclass[i]), ] > 0))
}
all_subclasses

all_subclasses$max_cells <- 0
for(i in 1:nrow(all_subclasses)){ #calculate max cells per subclass
  all_subclasses$max_cells[i] <- all_subclasses[i, 1 + which.max(all_subclasses[i, 2:3])] * 200
}
all_subclasses

all_subclasses$human_cells_per_cl <- 0
all_subclasses$marmoset_cells_per_cl <- 0

for(i in 1:nrow(all_subclasses)){
  all_subclasses$human_cells_per_cl[i] <- round(all_subclasses$max_cells[i] / all_subclasses$human[i])
  all_subclasses$marmoset_cells_per_cl[i] <- round(all_subclasses$max_cells[i] / all_subclasses$marmoset[i])
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

#####################################################################################################################
###############################       SCT norm and integration                        ###############################
#####################################################################################################################

human_data$orig.ident <- "human"
marmoset_data$orig.ident <- "marmoset"

all.data <- merge(x = human_data, y = list(marmoset_data), add.cell.ids = c("human", "marmoset"))

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

##############################################################################################################################################
############################################# Confusion heatmap  #############################################################################
##############################################################################################################################################

sample.combined$label_for_heatmap <- paste(sample.combined$orig.ident, sample.combined$cluster_label, sep = "_")

reorder_matrix <- function(matrix1, by.rows = TRUE) {
  if (by.rows == TRUE) {
    conf.order <- order(apply(matrix1, 1, which.max))
    matrix1.reordered <- matrix1[conf.order, ]
  } else {
    conf.order <- order(apply(matrix1, 2, which.max))
    matrix1.reordered <- matrix1[, conf.order]
  }
}

compare_cl <- function(cl, ref.cl,
                       plot.title = NA, plot.silent = TRUE,
                       heat.colors = colorRampPalette(c("white", "grey70", "black"))(100),
                       row.cl.num = min(length(unique(cl)),
                                        length(unique(ref.cl)))) {
  library(grid)
  library(pheatmap)
  
  conf1 <- table(cl, ref.cl)
  conf1 <- sweep(conf1, 1, rowSums(conf1), "/")
  conf2 <- reorder_matrix(conf1)
  
  # Cluster co-occurence
  cl.prop.cocl <- apply(conf1, 2, function(x) {
    grid1 <- expand.grid(x, x)
    min.prop <- apply(grid1, 1, min)
  })
  
  cl.prop.cocl.total <- apply(cl.prop.cocl, 1, sum)
  cl.prop.cocl.m <- matrix(cl.prop.cocl.total, nrow(conf1), nrow(conf1),
                           dimnames = list(rownames(conf1), rownames(conf1)))
  
  ph1 <- pheatmap(conf2, cutree_rows = row.cl.num, clustering_method = "ward.D2",
                  color = heat.colors, fontsize = 6,
                  main = plot.title, silent = plot.silent)
  return(list(conf = conf2, cocl = cl.prop.cocl.m, ph = ph1))
}

# Heatmap palette
heat.colors <- colorRampPalette(c("grey99", "orange", "red"))(100)
#heat.colors <- colorRampPalette(c("white", "grey70", "black"))(100)

# Compare clustering
ref.cl <- sample.combined$label_for_heatmap
cca.cl <- sample.combined$seurat_clusters.new
compare.species <- unique(sample.combined$orig.ident)

cl.conf <- compare_cl(ref.cl, cca.cl)
cocl <- cl.conf$cocl

cocl.subset <- cocl[grepl(compare.species[1], row.names(cocl)),
                    grepl(compare.species[2], row.names(cocl))]

row.names(cocl.subset) <- sub(paste0(compare.species[1], "_"), "",
                              row.names(cocl.subset))

colnames(cocl.subset) <- sub(paste0(compare.species[2], "_"), "",
                             colnames(cocl.subset))

dend1 <- readRDS("http://data.nemoarchive.org/biccn/lab/lein/2020_M1_study_analysis/Transcriptomics/human_10x_SS_integration/dend.RData")
dend.order <- labels(dend1)

cl.order <- intersect(dend.order, row.names(cocl.subset))

cocl.subset2 <- reorder_matrix(cocl.subset[cl.order, ], by.rows = FALSE)

grid.draw(cl.conf[["ph"]]$gtable)


library(fpc)
clus.method <- "single"

clus.num <- pamk(cocl, 1:(min(nrow(cocl), ncol(cocl)) - 1))$nc

ph1 <- pheatmap(cocl, clustering_method = clus.method,
                cutree_cols = clus.num, cutree_rows = clus.num,
                color = heat.colors,
                fontsize = 6)

pheatmap(cocl.subset2, cluster_rows = FALSE, cluster_cols = FALSE,
         color = heat.colors, 
         fontsize = 6)





