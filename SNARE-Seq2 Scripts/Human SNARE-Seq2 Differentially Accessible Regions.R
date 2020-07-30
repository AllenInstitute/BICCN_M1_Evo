#Code Underlying SNARE-Seq2 M1 Data Analyses 
#Blue B. Lake
#b1lake@eng.ucsd.edu

# Human SNARE-Seq2 Differentially Accessible Regions ---------------------------------------------
library(Signac)
library(Seurat)
library(dplyr)
library(Matrix)
library(methods)
library(swne)
library(cicero)
library(chromfunks)
library(ggplot2)
library(viridis)


####
##Sample Depth Function
####

#' sample cells based on the distribution of depth in peak matrix
#' @param pmatrix peak matrix, rows are peak positions and columns are cells
#' @param cls.groups factor contains cell annotation for each celltype (cluster)
#' @param posCluster factor contains positive cluster
#' @param sampleSize number of sampled cells for negative class (background) - defaut equal to posCluster
#' @return vector of cells from positive and negative class
sampleCellDepth <- function(pmatrix, cls.groups, posCluster, sampleSize=length(posCluster)){
  depth <- Matrix::colSums(pmatrix)
  depthPos <- depth[names(depth) %in% names(posCluster)]
  depthNeg <- depth[-1*which(names(depth) %in% names(posCluster))]
  negCells <- names(cls.groups)[-1*which(names(cls.groups) %in% names(posCluster))]
  # fit distribution
  densityEst <- density(depthPos, kernel = "gaussian", bw = 1)
  weights <- approx(densityEst$x, densityEst$y, xout=depthNeg,
                    yright = 0.00001,
                    yleft = 0.00001)$y
  sampledCells <- negCells[sample(seq(length(negCells)), size = sampleSize,
                                  prob = weights, replace=FALSE)]
  totalCells <- c(names(posCluster), sampledCells)
  return(totalCells)
}




####
##AC-Cluster DARs
####

load("~/NeMO_analysis_folder/SNARE/Analysis/Human/Zhang_BICCN-H_20190523-20190611_huMOp_Seurat.rda")
DefaultAssay(MOp.atac) <- "ATAC"
Idents(object = MOp.atac) <- "AC_cluster_tree_order"
Idents(object = MOp.atac) <- factor(Idents(object = MOp.atac), levels = 1:43)

count.matrix <- GetAssayData(object = MOp.atac, slot = "counts")
count.matrix@x[count.matrix@x > 0] <- 1 #binarize counts
clusters.use <- levels(Idents(MOp.atac))

DAR.list <- lapply(clusters.use, function(cl) {
  print(paste("Running for Cluster:", cl))
  
  pmatrix = count.matrix
  cls.groups = Idents(object = MOp.atac)
  cls.x <- WhichCells(MOp.atac, idents = cl, downsample = 10000)
  posCluster = droplevels(cls.groups[cls.x])
  sampleSize = 10000
  cls.x.samp <- sampleCellDepth(pmatrix,cls.groups,posCluster,sampleSize)
  
  clusters <- cls.groups[cls.x.samp]
  current.cluster.ids <- levels(clusters)
  new.cluster.ids <- current.cluster.ids
  new.cluster.ids[!new.cluster.ids %in% levels(posCluster)] <- "other"
  clusters <- plyr::mapvalues(clusters, from = current.cluster.ids, to = new.cluster.ids)
  
  cls.x.DARs <- CalcDiffAccess(count.matrix[,cls.x.samp], clusters, min.p = 1e-100)
  
  cls.x.DARs[[cl]]
  
})

names(DAR.list) <- clusters.use
Top_acDARs <- lapply(DAR.list, function(df) subset(df, qval < 0.01 & logfc > 1))
Top_acDARs.df <- do.call("rbind", lapply(Top_acDARs, as.data.frame)) 
cl<-unlist(lapply(rownames(Top_acDARs.df),function(x) unlist(strsplit(x,"[.]"))[1]))
loc <- unlist(lapply(rownames(Top_acDARs.df),function(x) unlist(strsplit(x,"[.]"))[2]))
Top_acDARs.df$cluster <- cl
Top_acDARs.df$site <- loc
rownames(Top_acDARs.df) <- gsub(" ", "_", rownames(Top_acDARs.df))
write.table(Top_acDARs.df, file="Zhang_BICCN-H_20190523-20190611_huMOp_AC-CL_All_DARs_Depth-Corrected.txt", sep = "\t", row.names=TRUE, col.names=TRUE)

#Plot Extended Data Fig. 5h
cl.mark <- Top_acDARs.df
cl.mark <- distinct(cl.mark, site, .keep_all = TRUE)
cl.mark %>% group_by(cluster) %>% top_n(100, logfc) -> top100
ave.DAR <- AverageExpression(MOp.atac, assays = "ATAC", features = top100$site, slot = "data" )

scaled <- t(scale(t(ave.DAR$ATAC)))
scaled <- scale(scaled)
range(scaled)
scaled[scaled < 0] <- 0
scaled[scaled > 5] <- 5
write.table(scaled, file="huMOp_AC-CL_DARs_Heatmap_Table.txt", sep = "\t", row.names=TRUE, col.names=TRUE)

ggHeat((scaled), rescaling = "none", clustering = "none", x.lab.size = 8, y.lab.size = 11,
) + scale_fill_viridis()




####
##Subclass DARs
####

load("~/NeMO_analysis_folder/SNARE/Analysis/Human/Zhang_BICCN-H_20190523-20190611_huMOp_Seurat.rda")
DefaultAssay(MOp.atac) <- "ATAC"
Idents(object = MOp.atac) <- "subclass"
order <- c("LAMP5", "SNCG", "VIP", "SST CHODL", "SST", "PVALB", 
           "L2-3 IT", "L5 IT", "L6 IT", "L6 IT Car3", "L5 ET", "L5-6 NP",  "L6 CT", "L6b",                                          
           "OPC", "Astro", "Oligo", "Micro-PVM", "VLMC", "Endo")  

Idents(object = MOp.atac) <- factor(Idents(object = MOp.atac), levels = order)
clusters.use <- levels(Idents(MOp.atac))

DAR.list <- lapply(clusters.use, function(cl) {
  print(paste("Running for Cluster:", cl))
  
  pmatrix = count.matrix
  cls.groups = Idents(object = MOp.atac)
  cls.x <- WhichCells(MOp.atac, idents = cl, downsample = 10000)
  posCluster = droplevels(cls.groups[cls.x])
  sampleSize = 10000
  cls.x.samp <- sampleCellDepth(pmatrix,cls.groups,posCluster,sampleSize)
  
  clusters <- cls.groups[cls.x.samp]
  current.cluster.ids <- levels(clusters)
  new.cluster.ids <- current.cluster.ids
  new.cluster.ids[!new.cluster.ids %in% levels(posCluster)] <- "other"
  clusters <- plyr::mapvalues(clusters, from = current.cluster.ids, to = new.cluster.ids)
  
  cls.x.DARs <- CalcDiffAccess(count.matrix[,cls.x.samp], clusters, min.p = 1e-100)
  
  cls.x.DARs[[cl]]
  
})

names(DAR.list) <- clusters.use
Top_acDARs <- lapply(DAR.list, function(df) subset(df, qval < 0.001 & logfc > 1))
Top_acDARs.df <- do.call("rbind", lapply(Top_acDARs, as.data.frame)) 
cl<-unlist(lapply(rownames(Top_acDARs.df),function(x) unlist(strsplit(x,"[.]"))[1]))
loc <- unlist(lapply(rownames(Top_acDARs.df),function(x) unlist(strsplit(x,"[.]"))[2]))
Top_acDARs.df$cluster <- cl
Top_acDARs.df$site <- loc
rownames(Top_acDARs.df) <- gsub(" ", "_", rownames(Top_acDARs.df))
write.table(Top_acDARs.df, file="Zhang_BICCN-H_20190523-20190611_huMOp_subclass-CL_all_DARs_Depth-Corrected.txt", sep = "\t", row.names=TRUE, col.names=TRUE)



####
##RNA-Level Cluster DARs
####

load("~/NeMO_analysis_folder/SNARE/Analysis/Human/Zhang_BICCN-H_20190523-20190611_huMOp_Seurat.rda")
DefaultAssay(MOp.atac) <- "ATAC"
Idents(object = MOp.atac) <- "RNA_cluster_tree_order"
Idents(object = MOp.atac) <- factor(Idents(object = MOp.atac), levels = 1:128)

count.matrix <- GetAssayData(object = MOp.atac, slot = "counts")
count.matrix@x[count.matrix@x > 0] <- 1 #binarize counts
clusters.use <- levels(Idents(MOp.atac))

DAR.list <- lapply(clusters.use, function(cl) {
  print(paste("Running for Cluster:", cl))
  
  pmatrix = count.matrix
  cls.groups = Idents(object = MOp.atac)
  cls.x <- WhichCells(MOp.atac, idents = cl, downsample = 10000)
  posCluster = droplevels(cls.groups[cls.x])
  sampleSize = 10000
  cls.x.samp <- sampleCellDepth(pmatrix,cls.groups,posCluster,sampleSize)
  
  clusters <- cls.groups[cls.x.samp]
  current.cluster.ids <- levels(clusters)
  new.cluster.ids <- current.cluster.ids
  new.cluster.ids[!new.cluster.ids %in% levels(posCluster)] <- "other"
  clusters <- plyr::mapvalues(clusters, from = current.cluster.ids, to = new.cluster.ids)
  
  cls.x.DARs <- CalcDiffAccess(count.matrix[,cls.x.samp], clusters, min.p = 1e-100)
  
  cls.x.DARs[[cl]]
  
})

names(DAR.list) <- clusters.use
Top_acDARs <- lapply(DAR.list, function(df) subset(df, qval < 0.05 & logfc > 1))
Top_acDARs.df <- do.call("rbind", lapply(Top_acDARs, as.data.frame)) 
cl<-unlist(lapply(rownames(Top_acDARs.df),function(x) unlist(strsplit(x,"[.]"))[1]))
loc <- unlist(lapply(rownames(Top_acDARs.df),function(x) unlist(strsplit(x,"[.]"))[2]))
Top_acDARs.df$cluster <- cl
Top_acDARs.df$site <- loc
rownames(Top_acDARs.df) <- gsub(" ", "_", rownames(Top_acDARs.df))
write.table(Top_acDARs.df, file="Zhang_BICCN-H_20190523-20190611_huMOp_Consensus-CL_All_DARs_Depth-Corrected.txt", sep = "\t", row.names=TRUE, col.names=TRUE)

#Plot Extended Data Fig. 5h
cl.mark <- Top_acDARs.df
cl.mark <- distinct(cl.mark, site, .keep_all = TRUE)
cl.mark %>% group_by(cluster) %>% top_n(100, logfc) -> top100
ave.DAR <- AverageExpression(MOp.atac, assays = "ATAC", features = na.omit(top100$site), slot = "data" )

scaled <- t(scale(t(ave.DAR$ATAC)))
scaled <- scale(scaled)
range(scaled)
scaled[scaled < 0] <- 0
scaled[scaled > 5] <- 5
write.table(scaled, file="huMOp_Consensus-CL_DARs_Heatmap_Table.txt", sep = "\t", row.names=TRUE, col.names=TRUE)

ggHeat((scaled), rescaling = "none", clustering = "none", x.lab.size = 8, y.lab.size = 11,
) + scale_fill_viridis()




####
##Subclass DARs (sub-sampled)
####

load("~/NeMO_analysis_folder/SNARE/Analysis/Human/Zhang_BICCN-H_20190523-20190611_huMOp_Seurat.rda")
DefaultAssay(MOp.atac) <- "ATAC"
Idents(object = MOp.atac) <- "subclass"
order <- c("LAMP5", "SNCG", "VIP", "SST CHODL", "SST", "PVALB", 
           "L2-3 IT", "L5 IT", "L6 IT", "L6 IT Car3", "L5 ET", "L5-6 NP",  "L6 CT", "L6b",                                          
           "OPC", "Astro", "Oligo", "Micro-PVM", "VLMC", "Endo")  

Idents(object = MOp.atac) <- factor(Idents(object = MOp.atac), levels = order)

#Subset to neuronal
MOp.atac <- subset(MOp.atac, idents = c("LAMP5", "SNCG", "VIP", "SST CHODL", "SST", "PVALB", 
                                        "L2-3 IT", "L5 IT", "L6 IT", "L6 IT Car3", "L5 ET", 
                                        "L5-6 NP",  "L6 CT", "L6b"))

count.matrix <- GetAssayData(object = MOp.atac, slot = "counts")
count.matrix@x[count.matrix@x > 0] <- 1 #binarize counts
clusters.use <- levels(Idents(MOp.atac))
clusterSize = 500
backgroundSize = 10000

DAR.list <- lapply(clusters.use, function(cl) {
  print(paste("Running for Cluster:", cl))
  
  pmatrix = count.matrix
  cls.groups = Idents(object = MOp.atac)
  cls.x <- WhichCells(MOp.atac, idents = cl, downsample = clusterSize)
  posCluster = droplevels(cls.groups[cls.x])
  sampleSize = backgroundSize
  cls.x.samp <- sampleCellDepth(pmatrix,cls.groups,posCluster,sampleSize)
  
  clusters <- cls.groups[cls.x.samp]
  current.cluster.ids <- levels(clusters)
  new.cluster.ids <- current.cluster.ids
  new.cluster.ids[!new.cluster.ids %in% levels(posCluster)] <- "other"
  clusters <- plyr::mapvalues(clusters, from = current.cluster.ids, to = new.cluster.ids)
  
  cls.x.DARs <- CalcDiffAccess(count.matrix[,cls.x.samp], clusters, min.p = 1e-100)
  
  cls.x.DARs[[cl]]
  
})

names(DAR.list) <- clusters.use
Top_acDARs <- lapply(DAR.list, function(df) subset(df, qval < 0.05 & logfc > 1))
Top_acDARs.df <- do.call("rbind", lapply(Top_acDARs, as.data.frame)) 
cl<-unlist(lapply(rownames(Top_acDARs.df),function(x) unlist(strsplit(x,"[.]"))[1]))
loc <- unlist(lapply(rownames(Top_acDARs.df),function(x) unlist(strsplit(x,"[.]"))[2]))
Top_acDARs.df$cluster <- cl
Top_acDARs.df$site <- loc
rownames(Top_acDARs.df) <- gsub(" ", "_", rownames(Top_acDARs.df))
write.table(Top_acDARs.df, file="Zhang_BICCN-H_20190523-20190611_huMOp_subclass-CL_all_DARs_Depth-Corrected_sub-sampled.txt", sep = "\t", row.names=TRUE, col.names=TRUE)



#Plot Figure 4f
AUC.dars = read.table("~/NeMO_analysis_folder/SNARE/Analysis/Human/Zhang_BICCN-H_20190523-20190611_huMOp_subclass-CL_all_DARs_Depth-Corrected_sub-sampled.AUC.txt", header = TRUE, sep = "\t")
cuttoff <- quantile(AUC.dars$AUC, 0.25)
#0.5017065

cl.mark <- AUC.dars[AUC.dars$AUC > cuttoff,]
cl.mark <- cl.mark[cl.mark$qval < 0.005,]
cl.mark <- cl.mark[cl.mark$logfc > 1,]
write.table(cl.mark, file="Zhang_BICCN-H_20190523-20190611_huMOp_subclass-CL_all_DARs_Depth-Corrected_sub-sampled_AUC.txt", sep = "\t", row.names=TRUE, col.names=TRUE)
cl.mark <- distinct(cl.mark, site, .keep_all = TRUE)
ave.DAR <- AverageExpression(MOp.atac, assays = "ATAC", features = as.character(cl.mark$site), slot = "data" )

scaled <- t(scale(t(ave.DAR$ATAC)))
scaled <- scale(scaled)
scaled[scaled < 0] <- 0
scaled[scaled > 5] <- 5
write.table(scaled, file="huMOp_Subclass_DARs_Heatmap_Table_sub-sampled_AUC.txt", sep = "\t", row.names=TRUE, col.names=TRUE)

ggHeat((scaled), rescaling = "none", clustering = "none", x.lab.size = 8, y.lab.size = 11,
) + scale_fill_viridis()



#Barplot of DAR proportions
order <- c("LAMP5", "SNCG", "VIP", "SST CHODL", "SST", "PVALB", 
           "L2-3 IT", "L5 IT", "L6 IT", "L6 IT Car3", "L5 ET", 
           "L5-6 NP",  "L6 CT", "L6b")
hu.dars <- as.data.frame(table(cl.mark$cluster))
colnames(hu.dars) <- c("subclass", "dars")
rownames(hu.dars) <- hu.dars$subclass
hu.dars <- hu.dars[order,]
write.table(hu.dars, file="huMOp_DAR_Totals_Subclass-Level_sub-sampled_AUC.txt", sep = "\t", row.names=TRUE, col.names=TRUE)
hu.dars <- mutate(data.frame(hu.dars), 
                  dars = (dars / sum(dars)) * 100)

barplot(hu.dars[,2], col = c("royalblue1"),cex.names = 0.3, names.arg = hu.dars$subclass, horiz = FALSE)



