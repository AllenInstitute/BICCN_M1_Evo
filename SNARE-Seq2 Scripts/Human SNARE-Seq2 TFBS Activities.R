#Code Underlying SNARE-Seq2 M1 Data Analyses 
#Blue B. Lake
#b1lake@eng.ucsd.edu

# Human SNARE-Seq2 TFBS Activities --------
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
##TFBS activities by subclass
####

load("~/NeMO_analysis_folder/SNARE/Analysis/Human/Zhang_BICCN-H_20190523-20190611_huMOp_Seurat.rda")
DefaultAssay(MOp.atac) <- "chromvar"
Idents(object = MOp.atac) <- "subclass"
order <- c("LAMP5", "SNCG", "VIP", "SST CHODL", "SST", "PVALB", 
           "L2-3 IT", "L5 IT", "L5 ET", "L5-6 NP", "L6 IT", "L6 IT Car3", "L6 CT", "L6b",                                          
           "OPC", "Astro", "Oligo", "Micro-PVM", "VLMC", "Endo")  
Idents(object = MOp.atac) <- factor(Idents(object = MOp.atac), levels = order)

tf.markers <- FindAllMarkers(
  object = MOp.atac,
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'nCount_ATAC'
)
write.table(tf.markers, file="Zhang_huMOp_Subclass_TF_Activity_Marker_table.txt", sep = "\t", row.names=TRUE, col.names=TRUE)



####
##TFBS activities by AC-level
####

DefaultAssay(MOp.atac) <- "chromvar"
Idents(object = MOp.atac) <- "AC_cluster_tree_order"
Idents(object = MOp.atac) <- factor(Idents(object = MOp.atac), levels = 1:43)

tf.markers <- FindAllMarkers(
  object = MOp.atac,
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'nCount_ATAC'
)
write.table(tf.markers, file="Zhang_huMOp_AC-level_TF_Activity_Marker_table.txt", sep = "\t", row.names=TRUE, col.names=TRUE)



####
##TFBS activities by RNA-level clusters - Glutamatergic
####

Idents(object = MOp.atac) <- "level1"
MOp.atac.glu <- subset(MOp.atac, idents = "Glutamatergic")

Idents(object = MOp.atac.glu) <- "RNA_cluster"
order <- read.delim(file="huMOp_Consensus_Cluster_Order.txt")
order <- order$cluster_label
Idents(object = MOp.atac.glu) <- factor(Idents(object = MOp.atac.glu), levels = order)

DefaultAssay(MOp.atac.glu) <- "chromvar"
tf.markers <- FindAllMarkers(
  object = MOp.atac.glu,
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'nCount_ATAC'
)
write.table(tf.markers, file="Zhang_huMOp_Consensus_Glutamatergic_TF_Activity_Marker_table.txt", sep = "\t", row.names=TRUE, col.names=TRUE)




####
##TFBS activities by RNA-level clusters - GABAergic
####

Idents(object = MOp.atac) <- "level1"
MOp.atac.gaba <- subset(MOp.atac, idents = "GABAergic")

Idents(object = MOp.atac.gaba) <- "RNA_cluster"
order <- read.delim(file="huMOp_Consensus_Cluster_Order.txt")
order <- order$cluster_label
Idents(object = MOp.atac.gaba) <- factor(Idents(object = MOp.atac.gaba), levels = order)

DefaultAssay(MOp.atac.gaba) <- "chromvar"
tf.markers <- FindAllMarkers(
  object = MOp.atac.gaba,
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'nCount_ATAC'
)
write.table(tf.markers, file="Zhang_huMOp_Consensus_GABAergic_TF_Activity_Marker_table.txt", sep = "\t", row.names=TRUE, col.names=TRUE)




####
##TFBS activities by RNA-level clusters - Non-neuronal
####

Idents(object = MOp.atac) <- "level1"
MOp.atac.nn <- subset(MOp.atac, idents = "Non-Neuronal")

Idents(object = MOp.atac.nn) <- "RNA_cluster"
order <- read.delim(file="huMOp_Consensus_Cluster_Order.txt")
order <- order$cluster_label
Idents(object = MOp.atac.nn) <- factor(Idents(object = MOp.atac.nn), levels = order)

DefaultAssay(MOp.atac.nn) <- "chromvar"
tf.markers <- FindAllMarkers(
  object = MOp.atac.nn,
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'nCount_ATAC'
)
write.table(tf.markers, file="Zhang_huMOp_Consensus_Non-Neuronal_TF_Activity_Marker_table.txt", sep = "\t", row.names=TRUE, col.names=TRUE)




####
##TFBS Activities ChC vs Basket
####

DefaultAssay(MOp.atac) <- "chromvar"
Idents(object = MOp.atac) <- "AC_cluster_tree_order"
Idents(object = MOp.atac) <- factor(Idents(object = MOp.atac), levels = 1:43)

#subset to pvalb subclass
MOp.atac.pvalb <- subset(MOp.atac, idents = 14:18)
MOp.atac.pvalb <- RenameIdents(MOp.atac.pvalb, '14' = "BC", '15' = "BC",
                               '16' = "BC", '17' = "BC", '18' = "ChC")
tf.markers <- FindAllMarkers(
  object = MOp.atac.pvalb,
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'nCount_ATAC'
)
write.table(tf.markers, file="Zhang_Human_MOp_AC-Level_PVALB_TF_Activity_Marker_table.txt", sep = "\t", row.names=TRUE, col.names=TRUE)

#Plot Figure 6e
DotPlot(MOp.atac.pvalb, features = c("NFIB", "RORA")) + RotatedAxis()





