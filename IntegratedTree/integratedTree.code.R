# Build integrated tree
# This script is used to get robust homologous clusters across species based on tree method
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
# read related file for inhibitory neurons
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
stability.measurements <- TreeStabilityDend(inh.dend$dend, cell_clusters, subsampled.dend, n.cores=10)
inh.dend$dend <- stability.measurements$dendrogram
saveRDS(inh.dend, "inh.tree.obj.rds")

# prune tree
inh.dend$dend <- TreeEntropy(inh.dend$dend, entropy.cutoff = 2.9)
dend_pruned <- pruneTreeEntropy(inh.dend$dend, cutoff=2.9)
saveRDS(dend_pruned, "inh.pruned.dend.rds")

# get homologous clusters
inh.homo.clusters <- getClusters(dend_pruned, plotTree=TRUE)

# homologous clusters when considering within-species composition
# add within-species cell composition and pruned tree based on new heurististic criteria
dend_pruned <- addWithinSpeciesComp(dend_pruned, humanAnnot = human.cellannot, 
                                    mouseAnnot = mouse.cellannot, marmoAnnot = marmoset.cellannot)
dend_pruned <- removeLowStabilityLeaf(dend_pruned)
upperLevelnodes <- getUpperLevelNode(dend_pruned, cutoff=0.65)

dend_pruned <- remove_nodes_nodePar(dend_pruned)
# get homologous clusters
inh.clusters.new <- getClusters(dend_pruned, upperlevelannot=upperLevelnodes, 
                                plotleafLabel=FALSE, furtherStop=TRUE, 
                                withinSpeciesStop=TRUE,plotClusterTree = TRUE)

#################################
# building excitatory neuron tree
# read related files for excitatory neurons
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
text(exc.dend$upperlevelannot$xy, labels=exc.dend$upperlevelannot$upperlabel, adj=c(0.2,0.2), cex=1, col="red")

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

# get homologous cluster for excitatory neurons
exc.homo.clusters <- getClusters(dend_pruned, plotTree=TRUE)

# add within-species cell composition and pruned tree based on new heurististic criteria
dend_pruned <- addWithinSpeciesComp(dend_pruned, humanAnnot = human.cellannot, 
                                    mouseAnnot = mouse.cellannot, marmoAnnot = marmoset.cellannot)
dend_pruned <- removeLowStabilityLeaf(dend_pruned)
upperLevelnodes <- getUpperLevelNode(dend_pruned, cutoff=0.65)

dend_pruned <- remove_nodes_nodePar(dend_pruned)
# get homologous clusters
exc.clusters.new <- getClusters(dend_pruned, upperlevelannot=upperLevelnodes, 
                                plotleafLabel=FALSE, furtherStop=TRUE, 
                                withinSpeciesStop=TRUE,plotClusterTree = TRUE)


######################################
# glia tree
# reading glia files
setwd("/home/qiwenhu/BICCN/dataset/allcells/integration/human_mouse_marmaset_integration/datasets/allen_data/Allen_tree/expression.matrix/glia_new_data")
expression.matrix <- readRDS("integrated_matrix_all_cells.RDS")
cell_clusters <- readRDS("all_cells_cluster_membership.RDS")

# exclude Oligo L3-6 OPALIN LRP4-AS1
exclude.cells <- names(human.cellannot[human.cellannot == "Oligo L3-6 OPALIN LRP4-AS1"])
expression.matrix <- expression.matrix[, -1*which(colnames(expression.matrix) %in% exclude.cells)]

cell_clusters <- cell_clusters[names(cell_clusters) %in% colnames(expression.matrix)]
upperlevelinfo <- readRDS("cross_species.upperleve.info.rds")
glia.upperlevelinfo <- upperlevelinfo[names(upperlevelinfo) %in% names(cell_clusters)]
species <- sapply(strsplit(names(cell_clusters), "_"), `[`, 1)
names(species) <- names(cell_clusters)

# building tree
d <- cluster.matrix.expression.distances(t(expression.matrix), 
                                         groups=cell_clusters, dist="cor", 
                                         useVariablegenes=FALSE, use.scaled.data=TRUE)
glia.dend <- buildTree(d, cell_clusters, glia.upperlevelinfo)
# plot excitatory neuron tree
plot(glia.dend$dend, leaflab = "none")
text(glia.dend$upperlevelannot$xy, labels=glia.dend$upperlevelannot$upperlabel, adj=c(0.2,0.2), cex=1, col="red")

dend <- TreeEntropy(glia.dend$dend, entropy.cutoff=2.6)
nnodes <- length(get_nodes_attr(dend, "height"))
dend <- assign_values_to_nodes(dend, "stability", rep(0.6, nnodes))

dend_pruned <- pruneTreeEntropy(dend, cutoff=2.6)
dend_pruned <- removeLowStabilityLeaf(dend_pruned, sizeCutoff=15)
dend_pruned <- speciesTree:::addWithinSpeciesComp(dend_pruned, humanAnnot = human.cellannot, 
                                    mouseAnnot = mouse.cellannot, marmoAnnot = marmoset.cellannot)

# get homologous clusters
glia.clusters <- getClusters(dend_pruned, plotleafLabel=FALSE, furtherStop=TRUE, 
                             withinSpeciesStop=TRUE)



