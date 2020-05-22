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
options(stringsAsFactors = FALSE)

#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
############################### Create Comb.dat object for merging steps later       ################################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################


dat.list=list()
cl.list = list()
cl.df.list=list()
anno.list = list()
shiny.dir = c(tenx="~/NeMO_analysis_folder/Transcriptomics/human_10x_SS_integration/source_data_10x/", sm="~/NeMO_analysis_folder/Transcriptomics/human_10x_SS_integration/source_data_10x/")


for(d in names(shiny.dir)){
  anno = read_feather(file.path(shiny.dir[d], "anno.feather"))
  cl.df = as.data.frame(unique(anno[,c("cluster_id","cluster_label","cluster_color","cell_class_label")]))
  row.names(cl.df) = cl.df$cluster_label
  cl = setNames(anno$cluster_label,anno$sample_id)
  
  select.cells <- anno %>%
    filter(cell_class_label == "exc") %>% #inh for inhibitory neurons
    select(sample_id) %>%
    unlist()
  
  counts = as.data.frame(read_feather(file.path(shiny.dir[d], "data_t.feather")))
  rownames(counts) <- counts$gene
  counts <- counts[ , -1] 
  norm.dat <- as.matrix(counts[ , select.cells])
  
  norm.dat <- Matrix(norm.dat, sparse = TRUE)
  norm.dat@x <- log2(norm.dat@x + 1)
  
  dat.list[[d]] = norm.dat
  cl.list[[d]] = cl
  cl.df.list[[d]] = cl.df
  anno.list[[d]] = anno
}


sets = names(dat.list)
comb.dat = prepare_unify(dat.list = dat.list, cl.list = cl.list, cl.df.list = cl.df.list)
de.param.list = list(de_param(q1.th=.4, q.diff.th=0.7, de.score.th=150,min.cells=10, min.genes = 4),  de_param(q1.th=0.5, q.diff.th=0.7, de.score.th=100,min.cells=2, min.genes = 1))
names(de.param.list) = sets
comb.dat$de.param.list = de.param.list


##create joint meta data.
for(d in sets){
  anno.list[[d]]$sample_id = paste(d, anno.list[[d]]$sample_id, sep=".")
}

common.cols = Reduce("intersect", sapply(anno.list, colnames))
tmp.df = do.call("rbind",lapply(anno.list, function(x)x[,common.cols]))
row.names(tmp.df) = tmp.df$sample_id

meta.df = comb.dat$meta.df
meta.df = cbind(meta.df, tmp.df[row.names(meta.df),])
comb.dat$meta.df = meta.df


#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
############################### Load in datasets / Create Seurat objects             ################################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################

sample.1_path <- "~/NeMO_analysis_folder/Transcriptomics/human_10x_SS_integration/source_data_10x/"
sample.1.anno.file <- paste0(sample.1_path, "/anno.feather")
sample.1.anno <- as.data.frame(read_feather(sample.1.anno.file))
sample.1.data_t.file <- paste0(sample.1_path, "/data_t.feather")
sample.1.data_t.df <- as.data.frame(read_feather(sample.1.data_t.file))
rownames(sample.1.data_t.df) = toupper(sample.1.data_t.df$gene)
sample.1.data_t.df = subset(sample.1.data_t.df, select = -gene)
norm.dat.sample.1 <- as.matrix(sample.1.data_t.df)
norm.dat.sample.1 <- Matrix(norm.dat.sample.1, sparse = TRUE)
norm.dat.sample.1@x <- log2(norm.dat.sample.1@x + 1)

sample.2_path <- "~/NeMO_analysis_folder/Transcriptomics/human_10x_SS_integration/source_data_10x/"
sample.2.anno.file <- paste0(sample.2_path, "/anno.feather")
sample.2.anno <- as.data.frame(read_feather(sample.2.anno.file))
sample.2.data_t.file <- paste0(sample.2_path, "/data_t.feather")
sample.2.data_t.df <- as.data.frame(read_feather(sample.2.data_t.file))
rownames(sample.2.data_t.df) = toupper(sample.2.data_t.df$gene)
sample.2.data_t.df = subset(sample.2.data_t.df, select = -gene)

#Breaks matrix into 10,000 column chunks and log2 normalizes them
num.breaks <- round(ncol(sample.2.data_t.df) / 10000)
sample.2.data_t.df_tmp <- list()
d <- 1
while(d < num.breaks + 2){
  last.break <- d * 10000
  first.break <- last.break - 9999
  if(last.break < ncol(sample.2.data_t.df)){
    matrix.to.add <- as.matrix(sample.2.data_t.df[ , first.break:last.break])
    matrix.to.add <- Matrix(matrix.to.add, sparse = TRUE)
    matrix.to.add@x <- log2(matrix.to.add@x + 1)
    sample.2.data_t.df_tmp[[d]] <- matrix.to.add   
  } 
  if(last.break > ncol(sample.2.data_t.df)){
    matrix.to.add <- as.matrix(sample.2.data_t.df[ , first.break:ncol(sample.2.data_t.df)])
    matrix.to.add <- Matrix(matrix.to.add, sparse = TRUE)
    matrix.to.add@x <- log2(matrix.to.add@x + 1)
    sample.2.data_t.df_tmp[[d]] <- matrix.to.add   
  }
  d <- d + 1
}

#Stitches subdivided log2 normalized chunks back together into single matrix
norm.dat.sample.2 <- sample.2.data_t.df_tmp[[1]]
d <- 2
while(d < num.breaks + 2){
  norm.dat.sample.2 <- cbind(norm.dat.sample.2, sample.2.data_t.df_tmp[[d]])
  d <- d + 1
}

sample.1.combined <- CreateSeuratObject(counts = norm.dat.sample.1, project = "sm") 
sample.2.combined <- CreateSeuratObject(counts = norm.dat.sample.2, project = "tenx") 

#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
############################### Pulling in metadata from scrattch .anno file and     ################################
############################### save to Seurat objects.                              ################################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################

metadata.to.add.sample.1                  <- data.frame(sample.1.anno$sample_id)          
metadata.to.add.sample.1$cluster_label    <- sample.1.anno$cluster_label          
metadata.to.add.sample.1$cell_class_label <- sample.1.anno$cell_class_label
metadata.to.add.sample.1$sample_id        <- sample.1.anno$sample_id
metadata.to.add.sample.1$cluster_id       <- sample.1.anno$cluster_id
metadata.to.add.sample.1$cluster_color    <- sample.1.anno$cluster_color
metadata.to.add.sample.1$layer_label      <- sample.1.anno$layer_label

sample.1.combined <- AddMetaData(sample.1.combined, metadata = metadata.to.add.sample.1$cluster_label,    col.name = "cluster_label") #Writes metadata to "Class" column in Seurat object
sample.1.combined <- AddMetaData(sample.1.combined, metadata = metadata.to.add.sample.1$cell_class_label, col.name = "cell_class_label")
sample.1.combined <- AddMetaData(sample.1.combined, metadata = metadata.to.add.sample.1$sample_id,        col.name = "sample_id")
sample.1.combined <- AddMetaData(sample.1.combined, metadata = metadata.to.add.sample.1$cluster_id,       col.name = "cluster_id")
sample.1.combined <- AddMetaData(sample.1.combined, metadata = metadata.to.add.sample.1$cluster_color,    col.name = "cluster_color")
sample.1.combined <- AddMetaData(sample.1.combined, metadata = metadata.to.add.sample.1$layer_label,      col.name = "layer_label")


metadata.to.add.sample.2                  <- data.frame(sample.2.anno$sample_id)          
metadata.to.add.sample.2$cluster_label    <- sample.2.anno$cluster_label          
metadata.to.add.sample.2$cell_class_label <- sample.2.anno$cell_class_label
metadata.to.add.sample.2$sample_id        <- sample.2.anno$sample_id
metadata.to.add.sample.2$cluster_id       <- sample.2.anno$cluster_id
metadata.to.add.sample.2$cluster_color    <- sample.2.anno$cluster_color

sample.2.combined <- AddMetaData(sample.2.combined, metadata = metadata.to.add.sample.2$cluster_label,    col.name = "cluster_label") #Writes metadata to "Class" column in Seurat object
sample.2.combined <- AddMetaData(sample.2.combined, metadata = metadata.to.add.sample.2$cell_class_label, col.name = "cell_class_label")
sample.2.combined <- AddMetaData(sample.2.combined, metadata = metadata.to.add.sample.2$sample_id,        col.name = "sample_id")
sample.2.combined <- AddMetaData(sample.2.combined, metadata = metadata.to.add.sample.2$cluster_id,       col.name = "cluster_id")
sample.2.combined <- AddMetaData(sample.2.combined, metadata = metadata.to.add.sample.2$cluster_color,    col.name = "cluster_color")

Idents(sample.1.combined) <- sample.1.combined@meta.data$cell_class_label
Idents(sample.2.combined) <- sample.2.combined@meta.data$cell_class_label

sample.1.combined <- subset(sample.1.combined, idents = c("exc")) #inh for inhibitory neurons
sample.2.combined <- subset(sample.2.combined, idents = c("exc")) #inh for inhibitory neurons

#####################################################################################################################
############################### Find marker genes for each cluster                    ###############################
#####################################################################################################################

Idents(sample.1.combined) <- sample.1.combined@meta.data$cluster_id
Idents(sample.2.combined) <- sample.2.combined@meta.data$cluster_id

sample.1.combined.ds <- subset(sample.1.combined, downsample = 300)
sample.2.combined.ds <- subset(sample.2.combined, downsample = 300)

Var.genes.sample.1 <- select_markers(sample.1.combined.ds@assays$RNA@counts, sample.1.combined.ds$cluster_id, n.markers = 20)
Var.genes.sample.1.markers <- Var.genes.sample.1$markers 
Var.genes.sample.2 <- select_markers(sample.2.combined.ds@assays$RNA@counts, sample.2.combined.ds$cluster_id, n.markers = 20)
Var.genes.sample.2.markers <- Var.genes.sample.2$markers 


#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
############################### Integrate and cluster                                 ###############################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################

VariableFeatures(sample.1.combined) <- Var.genes.sample.1.markers
VariableFeatures(sample.2.combined) <- Var.genes.sample.2.markers
total.Var.genes <- combine(Var.genes.sample.1.markers, Var.genes.sample.2.markers)
total.Var.genes <- unique(total.Var.genes$data)

sample.anchors <- FindIntegrationAnchors(object.list = list(sample.1.combined, sample.2.combined), anchor.features = total.Var.genes,  dims = 1:100)


sample.combined <- IntegrateData(anchorset = sample.anchors, dims = 1:100)
DefaultAssay(sample.combined) <- "integrated"

sample.combined <- ScaleData(sample.combined)
sample.combined <- RunPCA(sample.combined, npcs = 100)
ElbowPlot(sample.combined, ndims = 100)
sample.combined <- FindNeighbors(sample.combined, reduction = "pca", dims = 1:100, nn.eps = 0)
length(unique(sample.combined$cluster_label))
sample.combined <- FindClusters(sample.combined, resolution = 30, n.start = 100) 
sample.combined$seurat_clusters.new <- as.integer(sample.combined$seurat_clusters)
sample.combined <- RunTSNE(sample.combined, dims = 1:100)


#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
###############################              merging                                 ################################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
Idents(sample.combined) <- sample.combined$seurat_clusters.new
sample.combined.ds <- subset(sample.combined, downsample = 100)
test.genes <- select_markers(sample.combined.ds@assays$RNA@counts, sample.combined.ds$seurat_clusters.new, n.markers = 20)
merge.genes <- test.genes$markers 
genes.to.remove <- which(is.na(match(merge.genes, rownames(comb.dat$dat.list$tenx))) == TRUE)
merge.genes <- merge.genes[-genes.to.remove]

sets = names(dat.list)
de.param.list = list(de_param(q1.th=0.4, q.diff.th=0.7, de.score.th=80,min.cells=10, min.genes = 4),  de_param(q1.th=0.5, q.diff.th=0.7, de.score.th=20,min.cells=2, min.genes = 1))
names(de.param.list) = sets
comb.dat$de.param.list = de.param.list

cl.to.test <- sample.combined$seurat_clusters.new
names(cl.to.test) <- paste(sample.combined$orig.ident, ".", sample.combined$sample_id, sep = "")

merge.cl =  merge_cl_multiple(comb.dat, comb.dat$dat.list, cl=cl.to.test, anchor.genes=merge.genes, verbose=TRUE, merge.type = "directional")


length(unique(merge.cl))
merge.results <- as.numeric(merge.cl)
names(merge.results) <- gsub("^.*\\.","", names(merge.cl))


sample.combined <- AddMetaData(sample.combined, metadata = merge.results, col.name = "merged.cl")
