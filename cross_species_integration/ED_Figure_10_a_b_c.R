library(Seurat)
library(dplyr)
library(matrixStats)
library(Matrix)
library(ggplot2)
library(scrattch.hicat)

sample.combined <- readRDS("~/NeMO_analysis_folder/Transcriptomics/cross_species_integration/sample.combined_exc_4_species_integration.RDS")

sample.combined <- RunUMAP(sample.combined, dims = 1:100)

#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
############################### Visualization                                        ################################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################

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
  scale_color_manual(values = c("royalblue1", "green4", "maroon4", "sienna2")) 

##############################################################################################################################################
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
                  # annotation_row = ref.cl.anno[, -grep("cluster_label", colnames(ref.cl.anno))],
                  color = heat.colors, fontsize = 6,
                  main = plot.title, silent = plot.silent)
  return(list(conf = conf2, cocl = cl.prop.cocl.m, ph = ph1))
}

# Heatmap palette
#heat.colors <- colorRampPalette(c("grey99", "orange", "red"))(100)
heat.colors <- colorRampPalette(c("white", "grey70", "black"))(100)

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

####################
#customise with subclass 
annotation_row = data.frame(Subclass = sample.combined$subclass_label[match(rownames(cocl.subset2), sample.combined$cluster_label)])
annotation_col = data.frame(Subclass = sample.combined$subclass_label[match(colnames(cocl.subset2), sample.combined$cluster_label)])
rownames(annotation_row) <- rownames(cocl.subset2)
rownames(annotation_col) <- colnames(cocl.subset2)

tmp <- data.frame(subclass = unique(sample.combined$subclass_label))
tmp$color <- sample.combined$subclass_color[match(tmp$subclass, sample.combined$subclass_label)]
Subclass <- tmp$color
names(Subclass) <- tmp$subclass

ann_colors = list(
  Subclass = Subclass
)


pheatmap(cocl.subset2, cluster_rows = FALSE, cluster_cols = FALSE,
         color = heat.colors,  annotation_col = annotation_col, annotation_row = annotation_row,
         fontsize = 6, annotation_colors = ann_colors, cellwidth = 5, cellheight = 5)
