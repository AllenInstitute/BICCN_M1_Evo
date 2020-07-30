#Code Underlying SNARE-Seq2 M1 Data Analyses 
#Blue B. Lake
#b1lake@eng.ucsd.edu

# Marmoset SNARE-Seq2 Differentially Accessible Regions ---------
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

load("~/NeMO_analysis_folder/SNARE/Analysis/Marmoset/Zhang_BICCN-H_20190730_20190903_marMOp_Seurat.rda")
DefaultAssay(marMOp.atac) <- 'ATAC'
Idents(object = marMOp.atac) <- "AC_cluster"
Idents(object = marMOp.atac) <- factor(Idents(object = marMOp.atac), levels = 1:24)

clusters.use <- levels(Idents(marMOp.atac))

DAR.list <- lapply(clusters.use, function(cl) {
  print(paste("Running for Cluster:", cl))
  
  pmatrix = count.matrix
  cls.groups = Idents(object = marMOp.atac)
  cls.x <- WhichCells(marMOp.atac, idents = cl, downsample = 2000)
  posCluster = droplevels(cls.groups[cls.x])
  sampleSize = 2000
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
write.table(Top_acDARs.df, file="Zhang_BICCN-H_20190730_20190903_marMOp_AC-CL_All_DARs_Depth-Corrected.txt", sep = "\t", row.names=TRUE, col.names=TRUE)





####
##Subclass DARs
####

load("~/NeMO_analysis_folder/SNARE/Analysis/Marmoset/Zhang_BICCN-H_20190730_20190903_marMOp_Seurat.rda")
DefaultAssay(marMOp.atac) <- "ATAC"
Idents(object = marMOp.atac) <- "subclass"
order <- c("Lamp5", "Sncg", "Meis2", "Vip", "Sst Chodl", "Sst", "Pvalb", 
           "L2-3 IT", "L5 IT", "L6 IT", "L6 IT Car3", "L5 ET", "L5-6 NP", "L6 CT", "L6b",                                          
           "OPC", "Astro", "Oligo", "Micro-PVM", "VLMC","Peri", "Endo")

Idents(object = marMOp.atac) <- factor(Idents(object = marMOp.atac), levels = order)

#remove  clusters that are too small
marMOp.atac <- subset(marMOp.atac, idents = c("Meis2", "Sst Chodl"), invert = TRUE)

count.matrix <- GetAssayData(object = marMOp.atac, slot = "counts")
count.matrix@x[count.matrix@x > 0] <- 1 #binarize counts
clusters.use <- levels(Idents(marMOp.atac))

DAR.list <- lapply(clusters.use, function(cl) {
  print(paste("Running for Cluster:", cl))
  
  pmatrix = count.matrix
  cls.groups = Idents(object = marMOp.atac)
  cls.x <- WhichCells(marMOp.atac, idents = cl, downsample = 2000)
  posCluster = droplevels(cls.groups[cls.x])
  sampleSize = 2000
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
write.table(Top_acDARs.df, file="Zhang_BICCN-H_20190730_20190903_marMOp_subclass-CL_all_DARs_Depth-Corrected.txt", sep = "\t", row.names=TRUE, col.names=TRUE)





####
##Subclass DARs (Sub-sampled)
####

load("~/NeMO_analysis_folder/SNARE/Analysis/Marmoset/Zhang_BICCN-H_20190730_20190903_marMOp_Seurat.rda")
DefaultAssay(marMOp.atac) <- "ATAC"
Idents(object = marMOp.atac) <- "subclass"
order <- c("Lamp5", "Sncg", "Meis2", "Vip", "Sst Chodl", "Sst", "Pvalb", 
           "L2-3 IT", "L5 IT", "L6 IT", "L6 IT Car3", "L5 ET", "L5-6 NP", "L6 CT", "L6b",                                          
           "OPC", "Astro", "Oligo", "Micro-PVM", "VLMC","Peri", "Endo")

Idents(object = marMOp.atac) <- factor(Idents(object = marMOp.atac), levels = order)
table(Idents(object = marMOp.atac))

#Subset to neuronal
#remove too small clusters
marMOp.atac <- subset(marMOp.atac, idents = c("Lamp5", "Sncg", "Vip", "Sst", "Pvalb", 
                                              "L2-3 IT", "L5 IT", "L6 IT", "L6 IT Car3",
                                              "L5 ET", "L5-6 NP", "L6 CT", "L6b"))

count.matrix <- GetAssayData(object = marMOp.atac, slot = "counts")
count.matrix@x[count.matrix@x > 0] <- 1 #binarize counts
clusters.use <- levels(Idents(marMOp.atac))
table(Idents(marMOp.atac))
clusterSize = 200
backgroundSize = 2000

DAR.list <- lapply(clusters.use, function(cl) {
  print(paste("Running for Cluster:", cl))
  
  pmatrix = count.matrix
  cls.groups = Idents(object = marMOp.atac)
  cls.x <- WhichCells(marMOp.atac, idents = cl, downsample = clusterSize)
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
write.table(Top_acDARs.df, file="Zhang_BICCN-H_20190730_20190903_marMOp_subclass-CL_all_DARs_Depth-Corrected_sub-sampled-200.txt", sep = "\t", row.names=TRUE, col.names=TRUE)



#Plot Figure 4f
load("~/NeMO_analysis_folder/SNARE/Analysis/Marmoset/Zhang_BICCN-H_20190730_20190903_marMOp_Seurat.rda")
DefaultAssay(marMOp.atac) <- "ATAC"
Idents(object = marMOp.atac) <- "subclass"
order <- c("Lamp5", "Sncg", "Meis2", "Vip", "Sst Chodl", "Sst", "Pvalb", 
           "L2-3 IT", "L5 IT", "L6 IT", "L6 IT Car3", "L5 ET", "L5-6 NP", "L6 CT", "L6b",                                          
           "OPC", "Astro", "Oligo", "Micro-PVM", "VLMC","Peri", "Endo")

Idents(object = marMOp.atac) <- factor(Idents(object = marMOp.atac), levels = order)
neu.sub <- c("Lamp5", "Sncg", "Vip", "Sst Chodl", "Sst", "Pvalb", 
             "L2-3 IT", "L5 IT", "L6 IT", "L6 IT Car3", "L5 ET", "L5-6 NP", "L6 CT", "L6b")  
marMOp.atac <- subset(marMOp.atac, idents = neu.sub)

AUC.dars = read.table("~/NeMO_analysis_folder/SNARE/Analysis/Marmoset/Zhang_BICCN-H_20190730_20190903_marMOp_subclass-CL_all_DARs_Depth-Corrected_sub-sampled-200.AUC.txt", header = TRUE, sep = "\t")
cuttoff <- quantile(AUC.dars$AUC, 0.25)
#0.5064047

cl.mark <- AUC.dars[AUC.dars$AUC > cuttoff,]
cl.mark <- cl.mark[cl.mark$qval < 0.05,]
cl.mark <- cl.mark[cl.mark$logfc > 1,]
write.table(cl.mark, file="Zhang_BICCN-H_20190730_20190903_marMOp_subclass-CL_all_DARs_Depth-Corrected_sub-sampled_AUC.txt", sep = "\t", row.names=TRUE, col.names=TRUE)
cl.mark <- distinct(cl.mark, site, .keep_all = TRUE)
ave.DAR <- AverageExpression(marMOp.atac, assays = "ATAC", features = as.character(cl.mark$site), slot = "data" )

scaled <- t(scale(t(ave.DAR$ATAC)))
scaled <- scale(scaled)
scaled[scaled < 0] <- 0
scaled[scaled > 5] <- 5
write.table(scaled, file="marMOp_Subclass_DARs_Heatmap_Table_sub-sampled_AUC.txt", sep = "\t", row.names=TRUE, col.names=TRUE)

ggHeat((scaled), rescaling = "none", clustering = "none", x.lab.size = 8, y.lab.size = 11,
) + scale_fill_viridis()



#Barplot of subclass DAR proportions
order <- c("Lamp5", "Sncg", "Vip", "Sst", "Pvalb", 
           "L2-3 IT", "L5 IT", "L6 IT", "L6 IT Car3", "L5 ET", "L5-6 NP", "L6 CT", "L6b")
mar.dars <- as.data.frame(table(cl.mark$cluster))
colnames(mar.dars) <- c("subclass", "dars")
rownames(mar.dars) <- mar.dars$subclass
mar.dars <- mar.dars[order,]
write.table(mar.dars, file="marMOp_DAR_Totals_Subclass-Level_sub-sampled_AUC.txt", sep = "\t", row.names=TRUE, col.names=TRUE)
mar.dars <- mutate(data.frame(mar.dars), 
                   dars = (dars / sum(dars)) * 100)
rownames(mar.dars) <- mar.dars$subclass
neu.sub <- c("Lamp5", "Sncg", "Vip", "Sst Chodl", "Sst", "Pvalb", 
             "L2-3 IT", "L5 IT", "L6 IT", "L6 IT Car3", "L5 ET", "L5-6 NP", "L6 CT", "L6b")
mar.dars <- mar.dars[neu.sub,]
rownames(mar.dars)[rownames(mar.dars) == "NA"] <- "Sst Chodl"
mar.dars[is.na(mar.dars)] <- 0

barplot(mar.dars[,2], col = c("royalblue1"),cex.names = 0.3, names.arg = rownames(mar.dars), horiz = FALSE)







