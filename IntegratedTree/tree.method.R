# tree related functions

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

# generate pariwise pruned tree, 3-species
generate_pairwise_tree <- function(dend, facPrefix, colorpallete, minCutoff=0.28, maxCutoff=0.72){
  dend <- remove_branches_edgePar(dend)
  dend <- remove_nodes_nodePar(dend)
  dend <- dendSetWidthBysize(dend, scale=8)
  
  pairwise.dend <- dendSetColor2factorNormMixing(dend, reNormalize=TRUE, 
                                                 facPrefix=facPrefix, addValueToNode=TRUE, colorpallete)
  pairwise.dend <- speciesTree:::removeLowStabilityLeaf(pairwise.dend, sizeOnly=TRUE, sizeCutoff = 5)
  pairwise.dend.pruned <- pruneTree(pairwise.dend, minCutoff=minCutoff, maxCutoff=maxCutoff, heuristicStop = TRUE)
  pairwise.dend.pruned <- speciesTree:::removeLowStabilityLeaf(pairwise.dend.pruned, sizeCutoff = 1)
  
  upperLevelnodes <- getUpperLevelNode(pairwise.dend.pruned, cutoff=0.65)
  plot(pairwise.dend.pruned, leaflab="none")
  text(upperLevelnodes$xy, labels=upperLevelnodes$upperlabel, adj=c(0.5, -0.5), cex=1, col="red")
  
  pairwise.dend.pruned.adj <- speciesTree:::adjustTreeheight(pairwise.dend.pruned, scale=1)
  return(list(adj=pairwise.dend.pruned.adj, pruned=pairwise.dend.pruned))
}

# 2 species pairwise prunning
pairwise_prunning_tree <- function(expression.matrix, facPrefix, colorpallete, upperlevelinfo){
  ematrix <- expression.matrix[, grep(paste(facPrefix[1], facPrefix[2], sep="|"), colnames(expression.matrix))]
  pcell_clusters <- cell_clusters[names(cell_clusters) %in% colnames(ematrix)]
  species <- sapply(strsplit(names(pcell_clusters), "_"), `[`, 1)
  names(species) <- names(pcell_clusters)
  exc.upperlevelinfo <- upperlevelinfo[names(upperlevelinfo) %in% names(pcell_clusters)]
  
  d <- cluster.matrix.expression.distances(t(ematrix), groups=pcell_clusters, dist="cor", 
                                           useVariablegenes=FALSE,  use.scaled.data=TRUE)
  # manually building up tree
  dendr <- hclust(as.dist(d), method='ward.D2')
  dend <- as.dendrogram(dendr)
  
  dendr <- TransferDend(dend, renameCluster=TRUE, cls.groups = pcell_clusters)
  cls.groups <- dendr$new.groups
  
  dend <- dendr$dendrogram
  leafcontent <- dendr$leafcontent
  
  # raise tree
  dend <- raise.dendrogram(dend, 0.15)
  
  dend <- AddTreeAttribute(dend, species, leafcontent)
  dend <- dendSetWidthBysize(dend, scale=8)
  
  dend <- UpperLevelInfo(dend, cellannot=exc.upperlevelinfo, leafcontent, propCutoff = 0.1)
  upperLevelnodes <- getUpperLevelNode(dend, cutoff=0.65)
  
  dend <- NormTree(dend, upperLevelnodes, exc.upperlevelinfo, species)
  dend.pairwise <- dendSetColor2factorNormMixing(dend, reNormalize=TRUE, facPrefix=facPrefix, 
                                                 addValueToNode=TRUE, colorpallete)
  
  ## prunning
  dend.pruned.pairwise <- pruneTree(dend.pairwise, minCutoff=0.3,
                                    maxCutoff=0.65, heuristicStop = TRUE)
  dend.pruned.pairwise <- speciesTree:::removeLowStabilityLeaf(dend.pruned.pairwise, 
                                                 sizeOnly=TRUE, sizeCutoff = 1)
  upperLevelnodes <- getUpperLevelNode(dend.pruned.pairwise, cutoff=0.65)
  
  plot(dend.pruned.pairwise, leaflab="none")
  text(upperLevelnodes$xy, labels=upperLevelnodes$upperlabel, adj=c(0.5, -0.5), cex=1, col="red")
  
  dend.pruned.pairwise.adj <- speciesTree:::adjustTreeheight(dend.pruned.pairwise, scale=1)
  return(dend.pruned.pairwise.adj)
}

leaf_stability_sum <- function(sampledDend, facPrefix, upperlevelinfo, subsampledFrac=10){
  leafSum <- list()
  sampledDend <- sampledDend[1:subsampledFrac]
  for(i in 4:length(sampledDend)){
    dendr <- TransferDend(sampledDend[[i]]$dend, renameCluster=TRUE, cls.groups = sampledDend[[i]]$clusters)
    cls.groups <- dendr$new.groups
    
    species <- sapply(strsplit(names(cls.groups), "_"), `[`, 1)
    names(species) <- names(cls.groups)
    
    upperlevelannot <- upperlevelinfo[names(upperlevelinfo) %in% names(cls.groups)]
    
    dend <- dendr$dendrogram
    leafcontent <- dendr$leafcontent
    
    stability.measurements <- TreeStabilityDend(dend, cls.groups=cls.groups, sampledDend[-i], n.cores=10)
    dend <- stability.measurements$dendrogram
    
    # add cluster attribute to dendrogram
    dend <- AddTreeAttribute(dend, species, leafcontent)
    dend <- dendSetWidthBysize(dend, scale=8)
    
    dend <- UpperLevelInfo(dend, cellannot=upperlevelannot, leafcontent, propCutoff = 0.1)
    upperLevelnodes <- getUpperLevelNode(dend, cutoff=0.65)
    
    # normalize Tree
    dend <- NormTree(dend, upperLevelnodes, upperlevelannot, species)
    
    colorpallete <- colorRampPalette(c("cyan", "grey", "grey", "darkorchid1"))(101)
    dend.pairwise <- dendSetColor2factorNormMixing(dend, reNormalize=TRUE, facPrefix=facPrefix, 
                                                   addValueToNode=TRUE, colorpallete)
    dend.pruned <- pruneTree(dend.pairwise, minCutoff=0.15, maxCutoff=0.85, heuristicStop = TRUE)
    dend.pruned <- speciesTree:::removeLowStabilityLeaf(dend.pruned)
    #upperLevelnodes <- getUpperLevelNode(dend.pruned, cutoff=0.65)
    leafSum[[i]] <- leavesSubtree(dend.pruned)[[1]]
    leafSum[[i]]$comp <- paste(facPrefix[1], "vs", facPrefix[2], sep=" ")
    leafSum[[i]]$rep <- i
  }
  leafSum <- dplyr::bind_rows(leafSum)
  return(leafSum)
}

leafboxplot <- function(d){
  root_height <- attr(d, "height")
  leaflabels <- get_nodes_attr(d, "label")
  height <- get_nodes_attr(d, "height_original")
  depth <- round((root_height - height), 2)
  leaf.flag <- get_nodes_attr(d, "leaf")
  leaf.info <- data.frame(celltype=leaflabels, depth=depth, leaf.flag=leaf.flag)
  leaf.info <- leaf.info[which(leaf.info$leaf.flag==TRUE), ]
  if(length(which(is.na(leaf.info$depth))) > 0){
    leaf.info[is.na(leaf.info$depth), ]$depth <- root_height
  }
  leaf.info$celltype <- gsub("_.*", "", leaf.info$celltype)
  return(leaf.info)
}
