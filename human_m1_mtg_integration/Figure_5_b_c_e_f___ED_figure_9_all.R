library(Seurat)
library(Matrix)
library(matrixStats)
library(dplyr)
library(ggplot2)
library(feather)
library(scrattch.hicat)

#load integrated seurat object by class
sample.combined <- readRDS("~/NeMO_analysis_folder/Transcriptomics/human_10x_SS_integration/sample.combined_exc_integration.rds")

##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
############################################# tSNEs                 ##########################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################

#TSNE
DimPlot(sample.combined, reduction = "tsne", group.by = "orig.ident", pt.size = .2)          + ggtitle("Region") 

Idents(sample.combined) <- sample.combined$cluster_label
order.of.clusters <- levels(sample.combined@active.ident)
tsne_colors <- data.frame(clusters = order.of.clusters)
tsne_colors$colors <- sample.combined$cluster_color[match(tsne_colors$clusters, sample.combined$cluster_label)]
DimPlot(sample.combined, reduction = "tsne", pt.size = 0.7, label = TRUE, label.size = 2, split.by = "orig.ident") + 
  ggtitle("Clusters") + 
  NoLegend() +
  scale_color_manual(values = tsne_colors$colors)


Idents(sample.combined) <- sample.combined$seurat_clusters.new
DimPlot(sample.combined,label = TRUE, label.size = 3,  reduction = "tsne", pt.size = 0.7) + ggtitle("Integrated Clusters") + NoLegend()

sample.1.data$pial_depth <- 0
sample.2.data$pial_depth <- 0
sample.1.data$pial_depth[which(sample.1.data$layer == "L1")] <- 0.02
sample.1.data$pial_depth[which(sample.1.data$layer == "L2")] <- 0.07
sample.1.data$pial_depth[which(sample.1.data$layer == "L3")] <- 0.26
sample.1.data$pial_depth[which(sample.1.data$layer == "L5")] <- 0.55
sample.1.data$pial_depth[which(sample.1.data$layer == "L6")] <- 0.83
sample.2.data$pial_depth[which(sample.2.data$layer == "L1")] <- 0.03
sample.2.data$pial_depth[which(sample.2.data$layer == "L2")] <- 0.09
sample.2.data$pial_depth[which(sample.2.data$layer == "L3")] <- 0.22
sample.2.data$pial_depth[which(sample.2.data$layer == "L4")] <- 0.35
sample.2.data$pial_depth[which(sample.2.data$layer == "L5")] <- 0.50
sample.2.data$pial_depth[which(sample.2.data$layer == "L6")] <- 0.81

sample.combined$pial_depth <- 0
sample.combined$pial_depth[which(sample.combined$sample_id %in% names(sample.1.data$pial_depth))] <- sample.1.data$pial_depth
sample.combined$pial_depth[which(sample.combined$sample_id %in% names(sample.2.data$pial_depth))] <- sample.2.data$pial_depth

Idents(sample.combined) <- sample.combined$pial_depth
order.of.cells <- sort(unique(sample.combined$pial_depth))
levels(sample.combined) <- order.of.cells
tsne_colors <- data.frame(clusters = levels(sample.combined))
color_lookup <- data.frame(relative_depth = seq(0,1, by = 0.01))
colfunc <- colorRampPalette(c("#C4D4E5", "dodgerblue4"))
color_lookup$colors <- colfunc(101)
tsne_colors$color <- color_lookup$colors[match(tsne_colors$clusters, color_lookup$relative_depth)]

DimPlot(sample.combined, reduction = "tsne", pt.size = 0.7, group.by = "pial_depth", split.by = "orig.ident") + 
  ggtitle("Layers") + 
  scale_color_manual(values = tsne_colors$color) 

##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
############################################# Cluster overlap plot  ##########################################################################
##############################################################################################################################################
##############################################################################################################################################
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
                  # annotation_row = ref.cl.anno[, -grep("cluster_label", colnames(ref.cl.anno))],
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
                    grepl(compare.species[2], row.names(cocl))] #change 2 to 1 to plot M1 vs M1 

row.names(cocl.subset) <- sub(paste0(compare.species[1], "_"), "",
                              row.names(cocl.subset))

colnames(cocl.subset) <- sub(paste0(compare.species[2], "_"), "", #change 2 to 1 to plot M1 vs M1 
                             colnames(cocl.subset))

dend1 <- readRDS("~/NeMO_analysis_folder/Transcriptomics/human_10x_SS_integration/dend.RData")
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