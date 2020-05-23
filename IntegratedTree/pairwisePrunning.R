# pairwise pruning tree - figure 2g, Extended figure 2b, 3h
# Qiwen Hu - 2020

library(conos)
library(pagoda2)
library(clues)
library(parallel)
library(speciesTree)
source("tree.method.R")

###################################
# ihibitory neurons
### pariwise prunning and similarity
# human v.s. marmoset

# dend - ihibitory neuron origial tree
dend <- remove_nodes_nodePar(dend)
human.marmoset.dend <- pairwiseNormTree(dend, facPrefix=c("human", "marmoset"))
human.marmoset.dend <- speciesTree:::removeLowStabilityLeaf(human.marmoset.dend, sizeOnly = TRUE, sizeCutoff = 5)
dend.human.marmoset.pruned <- pruneTree(human.marmoset.dend, minCutoff=0.23, maxCutoff=0.77)
dend.human.marmoset.pruned <- speciesTree:::removeLowStabilityLeaf(dend.human.marmoset.pruned)
dend.human.marmoset.pruned.adj <- speciesTree:::adjustTreeheight(dend.human.marmoset.pruned, scale = 1)

# human v.s. mouse
human.mouse.dend <- pairwiseNormTree(dend, facPrefix=c("human", "mouse"))
dend.human.mouse.pruned <- pruneTree(human.mouse.dend, minCutoff=0.23, maxCutoff=0.77)
dend.human.mouse.pruned <- speciesTree:::removeLowStabilityLeaf(dend.human.mouse.pruned)
dend.human.mouse.pruned.adj <- speciesTree:::adjustTreeheight(dend.human.mouse.pruned)

# marmoset v.s. mouse
marmoset.mouse.dend <- pairwiseNormTree(dend, facPrefix=c("marmoset", "mouse"))
dend.marmoset.mouse.pruned <- pruneTree(marmoset.mouse.dend, minCutoff=0.23, maxCutoff=0.77)
dend.marmoset.mouse.pruned <- speciesTree:::removeLowStabilityLeaf(dend.marmoset.mouse.pruned)
dend.marmoset.mouse.pruned.adj <- speciesTree:::adjustTreeheight(dend.marmoset.mouse.pruned)

# plot tree panels
speciesTree:::alignTreePlot(list(dend.human.marmoset.pruned.adj, dend.human.mouse.pruned.adj, 
                                 dend.marmoset.mouse.pruned.adj))

## complexity plot
# load annotation
upperlevelinfo <- readRDS("cross_species.upperleve.info.rds")

# load subsampled trees
expression.matrix <- readRDS("integrated_matrix_all_cells.RDS")
files <- dir()
files <- files[grep("cluster_membership_iteration", files)]
subsampledClusters <- lapply(files, function(f){readRDS(f)})
subsampled.dend <- subSampleTree(t(expression.matrix), subsample.groups=subsampledClusters)

facPrefix <- c("human", "marmoset")
human.marmoset.leafs <- leaf_stability_sum(subsampled.dend, facPrefix, upperlevelinfo, subsampledFrac=10)
facPrefix <- c("human", "mouse")
human.mouse.leafs <- leaf_stability_sum(subsampled.dend, facPrefix, upperlevelinfo, subsampledFrac=10)
facPrefix <- c("marmoset", "mouse")
marmoset.mouse.leafs <- leaf_stability_sum(subsampled.dend, facPrefix, upperlevelinfo, subsampledFrac=10)
leafs.similarity <- rbind(human.marmoset.leafs, human.mouse.leafs, marmoset.mouse.leafs)
leafs.similarity <- leafs.similarity[-1*which(leafs.similarity$celltype=="Sst Chodl"), ]

# estimate standard deviation
human.marmoset.sd <- dplyr::bind_rows(lapply(unique(human.marmoset.leafs$celltype), function(r){
  nleafs <- human.marmoset.leafs[human.marmoset.leafs$celltype == r, ]
  data.frame(celltype=r, sd=sd(nleafs$nleaves))}))
human.mouse.sd <- dplyr::bind_rows(lapply(unique(human.mouse.leafs$celltype), function(r){
  nleafs <- human.mouse.leafs[human.mouse.leafs$celltype == r, ]
  data.frame(celltype=r, sd=sd(nleafs$nleaves))}))
marmoset.mouse.sd <- dplyr::bind_rows(lapply(unique(marmoset.mouse.leafs$celltype), function(r){
  nleafs <- marmoset.mouse.leafs[marmoset.mouse.leafs$celltype == r, ]
  data.frame(celltype=r, sd=sd(nleafs$nleaves))}))

# plotting
human.marmoset <- leavesSubtree(dend.human.marmoset.pruned.adj)[[1]]
human.marmoset$comp <- "human.marmo"
human.mouse <- leavesSubtree(dend.human.mouse.pruned.adj)[[1]]
human.mouse$comp <- "human.mouse"
mouse.marmoset <- leavesSubtree(dend.marmoset.mouse.pruned.adj)[[1]]
mouse.marmoset$comp <- "mouse.marmoset"

human.marmoset <- merge(human.marmoset, human.marmoset.sd, by=c(1))
human.mouse <- merge(human.mouse, human.mouse.sd, by=c(1))
mouse.marmoset <- merge(mouse.marmoset, marmoset.mouse.sd, by=c(1))
spComp.complexity <- rbind(human.marmoset, human.mouse, mouse.marmoset)

ggplot(data=spComp.complexity, aes(x=celltype, y=nleaves, fill=comp)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  geom_errorbar(aes(ymin=nleaves-sd, ymax=nleaves+sd),width=.2,
                position=position_dodge(.9)) + 
  theme_classic() + scale_x_discrete(limits=c("Lamp5", "Sncg", "Vip", "Sst", "Pvalb"))

############################
# excitatory neurons
# read related files
expression.matrix <- readRDS("integrated_matrix_all_cells.RDS")
cell_clusters <- readRDS("all_cells_cluster_membership.RDS")
cell_clusters <- cell_clusters[names(cell_clusters) %in% colnames(expression.matrix)]
upperlevelinfo <- readRDS("cross_species.upperleve.info.rds")
exc.upperlevelinfo <- upperlevelinfo[names(upperlevelinfo) %in% names(cell_clusters)]
species <- sapply(strsplit(names(cell_clusters), "_"), `[`, 1)
names(species) <- names(cell_clusters)

# original dendrogram for excitatory neurons
dend <- remove_branches_edgePar(dend)
dend <- remove_nodes_nodePar(dend)
dend <- dendSetWidthBysize(dend, scale=8)

# human vs marmoset
facPrefix <- c("human", "marmoset")
colorpallete <- colorRampPalette(c("cyan", "grey", "grey", "darkorchid1"))(101)
dend.exc.human.marmoset.pruned.adj <- generate_pairwise_tree(dend, facPrefix, colorpallete,
                                                             minCutoff = 0.28, maxCutoff = 0.72)[[1]]
dend.exc.human.marmoset.pruned <- generate_pairwise_tree(dend, facPrefix, colorpallete,
                                                         minCutoff = 0.28, maxCutoff = 0.72)[[2]]

# human vs mouse
colorpallete <- colorRampPalette(c("cyan", "grey", "grey","yellow"))(101) 
facPrefix <- c("human", "mouse")
dend.exc.human.mouse.pruned.adj <- generate_pairwise_tree(dend, facPrefix, colorpallete)[[1]]
dend.exc.human.mouse.pruned <- generate_pairwise_tree(dend, facPrefix, colorpallete)[[2]]

# marmose vs mouse
facPrefix <- c("marmoset", "mouse")
colorpallete <- colorRampPalette(c("darkorchid1", "grey", "grey", "yellow"))(101)
dend.exc.mouse.marmoset.pruned.adj <- generate_pairwise_tree(dend, facPrefix, colorpallete)[[1]]
dend.exc.mouse.marmoset.pruned <- generate_pairwise_tree(dend, facPrefix, colorpallete)[[2]]

speciesTree:::alignTreePlot(list(dend.exc.human.marmoset.pruned.adj, dend.exc.human.mouse.pruned.adj, 
                   dend.exc.marmoset.mouse.pruned.adj))

# complexity plot
# subsample dendrograms
files <- dir()
files <- files[grep("cluster_membership_iteration", files)]
subsampledClusters <- lapply(files, function(f){readRDS(f)})
subsampled.dend <- subSampleTree(t(expression.matrix), subsample.groups=subsampledClusters)

facPrefix <- c("human", "marmoset")
human.marmoset.leafs <- leaf_stability_sum(subsampled.dend, facPrefix, exc.upperlevelinfo, subsampledFrac=10)
facPrefix <- c("human", "mouse")
human.mouse.leafs <- leaf_stability_sum(subsampled.dend, facPrefix, exc.upperlevelinfo, subsampledFrac=10)
facPrefix <- c("marmoset", "mouse")
marmoset.mouse.leafs <- leaf_stability_sum(subsampled.dend, facPrefix, exc.upperlevelinfo, subsampledFrac=10)
leafs.similarity <- rbind(human.marmoset.leafs, human.mouse.leafs, marmoset.mouse.leafs)

# estimate standard deviation
human.marmoset.sd <- dplyr::bind_rows(lapply(unique(human.marmoset.leafs$celltype), function(r){
  nleafs <- human.marmoset.leafs[human.marmoset.leafs$celltype == r, ]
  data.frame(celltype=r, sd=sd(nleafs$nleaves))}))
human.mouse.sd <- dplyr::bind_rows(lapply(unique(human.mouse.leafs$celltype), function(r){
  nleafs <- human.mouse.leafs[human.mouse.leafs$celltype == r, ]
  data.frame(celltype=r, sd=sd(nleafs$nleaves))}))
marmoset.mouse.sd <- dplyr::bind_rows(lapply(unique(marmoset.mouse.leafs$celltype), function(r){
  nleafs <- marmoset.mouse.leafs[marmoset.mouse.leafs$celltype == r, ]
  data.frame(celltype=r, sd=sd(nleafs$nleaves))}))

human.marmoset <- leavesSubtree(dend.exc.human.marmoset.pruned.adj)[[1]]
human.marmoset$comp <- "human.marmo"
human.mouse <- leavesSubtree(dend.exc.human.mouse.pruned.adj)[[1]]
human.mouse$comp <- "human.mouse"
mouse.marmoset <- leavesSubtree(dend.exc.marmoset.mouse.pruned.adj)[[1]]
mouse.marmoset$comp <- "mouse.marmoset"

human.marmoset <- merge(human.marmoset, human.marmoset.sd, by=c(1))
human.mouse <- merge(human.mouse, human.mouse.sd, by=c(1))
human.mouse <- rbind(human.mouse, data.frame(celltype="L5/6 NP", nleaves=0, comp="human.mouse", sd=0))
mouse.marmoset <- merge(mouse.marmoset, marmoset.mouse.sd, by=c(1))
mouse.marmoset <- rbind(mouse.marmoset, data.frame(celltype="L5/6 NP", nleaves=0, comp="mouse.marmoset", sd=0))
mouse.marmoset <- rbind(mouse.marmoset, data.frame(celltype="L6b", nleaves=0, comp=mouse.marmoset[1,]$comp, sd=0))

spComp.complexity <- rbind(human.marmoset, human.mouse, mouse.marmoset)

ggplot(data=spComp.complexity, aes(x=celltype, y=nleaves, fill=comp)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  geom_errorbar(aes(ymin=nleaves-sd, ymax=nleaves+sd),width=.2,
                position=position_dodge(.9)) + theme_classic() +
  scale_x_discrete(limits=c("L2/3`` IT", "L5 IT", "L6 IT", "L5 ET", "L5/6 NP", "L6 CT",
                            "L6b"))


