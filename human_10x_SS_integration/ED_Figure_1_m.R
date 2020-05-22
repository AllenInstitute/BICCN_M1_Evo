library(Seurat)
library(dplyr)
library(Matrix)
library(matrixStats)
library(feather)
library(readxl)
options(stringsAsFactors = FALSE)

inh <- readRDS("~/NeMO_analysis_folder/Transcriptomics/human_10x_SS_integration/sample.combined_glia_integration.rds")
unique(inh$orig.ident)
table(inh$cluster_label, inh$orig.ident)
table(inh$merged.cl, inh$orig.ident)
length(unique(inh$merged.cl))


final_cl <- read.csv("~/NeMO_analysis_folder/Transcriptomics/human_10x_SS_integration/Figure_ED1_source_data/final_human_anno.csv")
tmp <- data.frame(sample_id = inh$sample_id)
clusters <- data.frame(sample_id = tmp[which(tmp$sample_id %in% final_cl$sample_id), ])
clusters$cluster_label <- final_cl$RNA_cluster_label[match(clusters$sample_id, final_cl$sample_id)]
clusters$seurat_cl <- inh$merged.cl[match(clusters$sample_id, inh$sample_id)]
clusters$cluster_color <- final_cl$RNA_cluster_color[match(clusters$sample_id, final_cl$sample_id)]
clusters$subclass_label <- final_cl$subclass_label[match(clusters$sample_id, final_cl$sample_id)]
clusters$subclass_color <- final_cl$subclass_color[match(clusters$sample_id, final_cl$sample_id)]

inh$final_cl <- clusters$cluster_label[match(inh$merged.cl, clusters$seurat_cl)]
inh$final_cl_color <- clusters$cluster_color[match(inh$merged.cl, clusters$seurat_cl)]
inh$subclass <- clusters$subclass_label[match(inh$merged.cl, clusters$seurat_cl)]
inh$subclass_color <- clusters$subclass_color[match(inh$merged.cl, clusters$seurat_cl)]

##############################################################################################################################################
##############################################################################################################################################
############################################# tSNE               #############################################################################
##############################################################################################################################################

DimPlot(inh, reduction = "tsne", group.by = "orig.ident", order = c("sm", "tenx"), pt.size = 1)

#UMAP cluster
data_to_plot <- data.frame(inh@reductions$tsne@cell.embeddings)
inh$new_id <- colnames(inh)
data_to_plot$cluster <- inh$final_cl[match(rownames(data_to_plot), inh$new_id)]
data_to_plot$species <- inh$orig.ident[match(rownames(data_to_plot), inh$new_id)]
data_to_plot$plot_order <- 0
for(i in 1:nrow(data_to_plot)){
  data_to_plot$plot_order[i] <- sample(1:100000, 1)
  print(i)
}
gc()

data_to_plot <- arrange(data_to_plot, plot_order)
data_to_plot$color <- inh$final_cl_color[match(data_to_plot$cluster, inh$final_cl)]

ggplot(data_to_plot,
       aes(x = tSNE_1, y = tSNE_2, color = color)) +
  geom_point(size = 1.5, color = data_to_plot$color)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 

#UMAP subclass
data_to_plot <- data.frame(inh@reductions$tsne@cell.embeddings)
inh$new_id <- colnames(inh)
data_to_plot$cluster <- inh$subclass[match(rownames(data_to_plot), inh$new_id)]
data_to_plot$species <- inh$orig.ident[match(rownames(data_to_plot), inh$new_id)]
data_to_plot$plot_order <- 0
for(i in 1:nrow(data_to_plot)){
  data_to_plot$plot_order[i] <- sample(1:100000, 1)
  print(i)
}
gc()

data_to_plot <- arrange(data_to_plot, plot_order)
data_to_plot$color <- inh$subclass_color[match(data_to_plot$cluster, inh$subclass)]

ggplot(data_to_plot,
       aes(x = tSNE_1, y = tSNE_2, color = color)) +
  geom_point(size = 1.5, color = data_to_plot$color)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 

##############################################################################################################################################
##############################################################################################################################################
############################################# Confusion Plot     #############################################################################
##############################################################################################################################################
Idents(inh) <- inh$orig.ident

smarter <- subset(inh, idents = "sm")
tenx <- subset(inh, idents = "tenx")


iter.cl <- setNames(as.factor(smarter$cluster_id), colnames(smarter))
iter.cl.df <- data.frame(cluster_id = 1:nlevels(iter.cl),
                         cluster_label = levels(iter.cl),
                         stringsAsFactors = FALSE)
cl <- setNames(as.factor(smarter$final_cl), colnames(smarter))
row.names(iter.cl.df) <- iter.cl.df$cluster_label
compare.result <- compare_annotate(cl, ref.cl =  iter.cl, ref.cl.df =  iter.cl.df, rename = FALSE, reorder = TRUE)
compare.result$g +  theme(axis.text.y = element_text(size = 5), axis.text.x = element_text(angle = 90)) + scale_size(range = c(0, 5))

iter.cl <- setNames(as.factor(tenx$cluster_id), colnames(tenx))
iter.cl.df <- data.frame(cluster_id = 1:nlevels(iter.cl),
                         cluster_label = levels(iter.cl),
                         stringsAsFactors = FALSE)
cl <- setNames(as.factor(tenx$final_cl), colnames(tenx))
row.names(iter.cl.df) <- iter.cl.df$cluster_label
compare.result <- compare_annotate(cl, ref.cl =  iter.cl, ref.cl.df =  iter.cl.df, rename = FALSE, reorder = TRUE)
compare.result$g +  theme(axis.text.y = element_text(size = 5), axis.text.x = element_text(angle = 90)) + scale_size(range = c(0, 5))

