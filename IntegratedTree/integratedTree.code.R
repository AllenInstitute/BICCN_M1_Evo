# build integrated tree
# Qiwen Hu - 2020

library(conos)
library(pagoda2)
library(speciesTree)
library(clues)
library(parallel)
library(Seurat)
library(ggplot2)

buildTree <- function(d, cell_clusters, upperlevelinfo, stability=TRUE){
  dendr <- hclust(as.dist(d), method='ward.D2')
  dend <- as.dendrogram(dendr)
  dendr <- TransferDend(dend, renameCluster=TRUE, cls.groups = cell_clusters)
  cls.groups <- dendr$new.groups
  dend <- dendr$dendrogram
  leafcontent <- dendr$leafcontent
  
  # raise tree
  dend <- raise.dendrogram(dend, 0.15)
  dend <- AddTreeAttribute(dend, species, leafcontent)
  dend <- dendSetWidthBysize(dend, scale=8)
  dend <- UpperLevelInfo(dend, cellannot=upperlevelinfo, leafcontent, propCutoff = 0.1)
  upperLevelnodes <- getUpperLevelNode(dend, cutoff=0.65)
  
  #normalize tree
  dend <- NormTree(dend, upperLevelnodes, upperlevelinfo, species)
  dend <- dendSetColorByNormMixing(dend)
  
  return(list(dend=dend, upperlevelannot=upperLevelnodes))
}

############################
# building inhibitory tree
# reading expression data
expression.matrix <- readRDS("integrated_matrix_all_cells.RDS")
cell_clusters <- readRDS("all_cells_cluster_membership.RDS")

# loading annotation
human.annot <- readRDS("human_annotation_file.RDS")
human.annot$sample_id <- paste("human", human.annot$sample_id, sep="_")
mouse.annot <- readRDS("mouse_annotation_file.RDS")
mouse.annot$sample_id <- paste("mouse", mouse.annot$sample_id, sep="_")
marmoset.annot <- readRDS("marmoset_annotation_file.RDS")
marmoset.annot$sample_id <- paste("marmoset", marmoset.annot$sample_id, sep="_")
upperlevelinfo <- c(human.annot$subclass_label, mouse.annot$subclass_label, marmoset.annot$subclass_label)
names(upperlevelinfo) <- c(human.annot$sample_id, mouse.annot$sample_id, marmoset.annot$sample_id)

human.cellannot <- setNames(human.annot$within_species_cluster_label, human.annot$sample_id)
marmoset.cellannot <- setNames(marmoset.annot$within_species_cluster_label, marmoset.annot$sample_id)
mouse.cellannot <- setNames(mouse.annot$within_species_cluster_label, mouse.annot$sample_id)

species <- sapply(strsplit(names(cell_clusters), "_"), `[`, 1)
names(species) <- names(cell_clusters)
inh.upperlevelinfo <- upperlevelinfo[names(upperlevelinfo) %in% names(cell_clusters)]

# building tree
d <- cluster.matrix.expression.distances(t(expression.matrix), groups=cell_clusters, dist="cor", 
                                         useVariablegenes=FALSE,  use.scaled.data=TRUE)

inh.dend <- buildTree(d, cell_clusters, inh.upperlevelinfo, stability=FALSE)

# plot inhibitory neuron tree
plot(inh.dend$dend, leaflab = "none")
text(inh.dend$upperlevelannot$xy, labels=inh.dend$upperlevelannot$upperlabel, adj=c(0.2,0.2), cex=1, col="red")

## add stability score
files <- dir()
files <- files[grep("cluster_membership_iteration", files)]
subsampledClusters <- lapply(files, function(f){readRDS(f)})

subsampled.dend <- subSampleTree(t(expression.matrix), subsample.groups=subsampledClusters)
stability.measurements <- TreeStabilityDend(inh.dend$dend, cls.groups, subsampled.dend, n.cores=10)
inh.dend$dend <- stability.measurements$dendrogram
saveRDS(inh.dend, "inh.tree.obj.rds")

# prune tree
inh.dend$dend <- TreeEntropy(inh.dend$dend, entropy.cutoff = 2.9)
dend_pruned <- pruneTreeEntropy(inh.dend$dend, cutoff=2.9)
saveRDS(dend_pruned, "inh.pruned.dend.rds")

######
# building excitatory neuron tree
expression.matrix <- readRDS("integrated_matrix_all_cells.RDS")
cell_clusters <- readRDS("all_cells_cluster_membership.RDS")
cell_clusters <- cell_clusters[names(cell_clusters) %in% colnames(expression.matrix)]
upperlevelinfo <- readRDS("cross_species.upperleve.info.rds")
exc.upperlevelinfo <- upperlevelinfo[names(upperlevelinfo) %in% names(cell_clusters)]
species <- sapply(strsplit(names(cell_clusters), "_"), `[`, 1)
names(species) <- names(cell_clusters)

d <- cluster.matrix.expression.distances(t(expression.matrix), groups=cell_clusters, dist="cor", 
                                         useVariablegenes=FALSE, use.scaled.data=TRUE)

exc.dend <- buildTree(d, cell_clusters, exc.upperlevelinfo)

# plot excitatory neuron tree
plot(exc.dend$dend, leaflab = "none")
text(exc.dend$upperlevelannot$xy, labels=inh.dend$upperlevelannot$upperlabel, adj=c(0.2,0.2), cex=1, col="red")

# add stability
files <- dir()
files <- files[grep("cluster_membership_iteration", files)]
subsampledClusters <- lapply(files, function(f){readRDS(f)})

subsampled.dend <- subSampleTree(t(expression.matrix), subsample.groups=subsampledClusters)
stability.measurements <- TreeStabilityDend(exc.dend$dend, cell_clusters, subsampled.dend, n.cores=10)
exc.dend$dend <- stability.measurements$dendrogram
saveRDS(inh.dend, "exc.tree.obj.rds")

#prune tree
dend <- TreeEntropy(exc.dend$dend, entropy.cutoff=2.75)
dend_pruned <- pruneTreeEntropy(dend, cutoff=2.75)
dend_pruned <- removeLowStabilityLeaf(dend_pruned)

upperLevelnodes <- getUpperLevelNode(dend_pruned, cutoff=0.65)
plot(dend_pruned, leaflab="none", main="Excitatory neurons")
text(upperLevelnodes$xy, labels=upperLevelnodes$upperlabel, adj=c(0.5, -0.5), cex=1, col="red")
saveRDS(dend_pruned, "exc.pruned.tree.rds")
