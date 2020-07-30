#Code Underlying SNARE-Seq2 M1 Data Analyses 
#Blue B. Lake
#b1lake@eng.ucsd.edu

# Marmoset SNARE-Seq2 TFBS Activities --------
library(Signac)
library(Seurat)
library(dplyr)
library(Matrix)
library(methods)
library(swne)
library(chromfunks)
library(ggplot2)
library(viridis)


####
##TFBS Activities by subclasses
####

load("~/NeMO_analysis_folder/SNARE/Analysis/Marmoset/Zhang_BICCN-H_20190730_20190903_marMOp_Seurat.rda")
DefaultAssay(marMOp.atac) <- "chromvar"
Idents(object = marMOp.atac) <- "subclass"
order <- c("Lamp5", "Sncg", "Meis2", "Vip", "Sst Chodl", "Sst", "Pvalb", 
           "L2-3 IT", "L5 IT", "L5 ET", "L5-6 NP", "L6 IT", "L6 IT Car3", "L6 CT", "L6b",                                          
           "OPC", "Astro", "Oligo", "Micro-PVM", "VLMC","Peri", "Endo")
Idents(object = marMOp.atac) <- factor(Idents(object = marMOp.atac), levels = order)

#remove clusters that are too small 
marMOp.atac <- subset(marMOp.atac, idents = c("Meis2", "Sst Chodl"), invert = TRUE)

tf.markers <- FindAllMarkers(
  object = marMOp.atac,
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'nCount_ATAC'
)
write.table(tf.markers, file="Zhang_Marmoset_MOp_Subclass_TF_Activity_Marker_table.txt", sep = "\t", row.names=TRUE, col.names=TRUE)



####
##TFBS Activities by AC-level
####

load("~/NeMO_analysis_folder/SNARE/Analysis/Marmoset/Zhang_BICCN-H_20190730_20190903_marMOp_Seurat.rda")
DefaultAssay(marMOp.atac) <- "chromvar"
Idents(object = marMOp.atac) <- "AC_cluster"
Idents(object = marMOp.atac) <- factor(Idents(object = marMOp.atac), levels = 1:24)

tf.markers <- FindAllMarkers(
  object = marMOp.atac,
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'nCount_ATAC'
)
write.table(tf.markers, file="Zhang_Marmoset_MOp_AC-Level_TF_Activity_Marker_table.txt", sep = "\t", row.names=TRUE, col.names=TRUE)




####
##TFBS Activities ChC vs Basket
####

load("~/NeMO_analysis_folder/SNARE/Analysis/Marmoset/Zhang_BICCN-H_20190730_20190903_marMOp_Seurat.rda")
DefaultAssay(marMOp.atac) <- "chromvar"
Idents(object = marMOp.atac) <- "AC_cluster"
Idents(object = marMOp.atac) <- factor(Idents(object = marMOp.atac), levels = 1:24)

#subset to pvalb subclass
marMOp.atac.pvalb <- subset(marMOp.atac, idents = 5:6)
marMOp.atac.pvalb <- RenameIdents(marMOp.atac.pvalb, '5' = "BC", '6' = "ChC")
tf.markers <- FindAllMarkers(
  object = marMOp.atac.pvalb,
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'nCount_ATAC'
)
write.table(tf.markers, file="Zhang_Marmoset_MOp_AC-Level_PVALB_TF_Activity_Marker_table.txt", sep = "\t", row.names=TRUE, col.names=TRUE)

#Plot Figure 6e
DotPlot(marMOp.atac.pvalb, features = c("NFIB", "RORA")) + RotatedAxis()


