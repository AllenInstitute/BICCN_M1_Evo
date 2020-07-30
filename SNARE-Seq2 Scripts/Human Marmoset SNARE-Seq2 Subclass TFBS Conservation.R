#Code Underlying SNARE-Seq2 M1 Data Analyses 
#Blue B. Lake
#b1lake@eng.ucsd.edu

# Human Marmoset SNARE-Seq2 Subclass TFBS Conservation ------------------
library(Signac)
library(Seurat)
library(dplyr)
library(Matrix)
library(methods)
library(swne)
library(chromfunks)
library(ggplot2)
library(viridis)
set.seed(1234)


###Analyses underlying Fig. 4i


####
##load seurat objects 
####

#Marmoset
load("~/NeMO_analysis_folder/SNARE/Analysis/Marmoset/Zhang_BICCN-H_20190730_20190903_marMOp_Seurat.rda")
DefaultAssay(marMOp.atac) <- "chromvar"
Idents(object = marMOp.atac) <- "subclass"
order <- c("Lamp5", "Sncg", "Meis2", "Vip", "Sst Chodl", "Sst", "Pvalb", 
           "L2-3 IT", "L5 IT", "L5 ET", "L5-6 NP", "L6 IT", "L6 IT Car3", "L6 CT", "L6b",                                          
           "OPC", "Astro", "Oligo", "Micro-PVM", "VLMC","Peri", "Endo")
Idents(object = marMOp.atac) <- factor(Idents(object = marMOp.atac), levels = order)
#remove clusters not present in human, too sparse, non-neuronal
marMOp.atac <- subset(marMOp.atac, idents = c("Meis2", "Sst Chodl","OPC", "Astro", "Oligo", "Micro-PVM", "VLMC","Peri", "Endo"), invert = TRUE)


#Human
load("~/NeMO_analysis_folder/SNARE/Analysis/Human/Zhang_BICCN-H_20190523-20190611_huMOp_Seurat.rda")
DefaultAssay(MOp.atac) <- "chromvar"
Idents(object = MOp.atac) <- "subclass"
order <- c("LAMP5", "SNCG", "VIP", "SST CHODL", "SST", "PVALB", 
           "L2-3 IT", "L5 IT", "L5 ET", "L5-6 NP", "L6 IT", "L6 IT Car3", "L6 CT", "L6b",                                          
           "OPC", "Astro", "Oligo", "Micro-PVM", "VLMC", "Endo")  
Idents(object = MOp.atac) <- factor(Idents(object = MOp.atac), levels = order)
#remove clusters not present in marmoset
MOp.atac <- subset(MOp.atac, idents = c("SST CHODL","OPC", "Astro", "Oligo", "Micro-PVM", "VLMC", "Endo"), invert = TRUE)



####
##Correlation on differentially active TFBSs
####

DefaultAssay(MOp.atac) <- "chromvar"
tf.markers <- FindAllMarkers(
  object = MOp.atac,
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'nCount_ATAC'
)

hu.tf.mark <- tf.markers[tf.markers$p_val < 0.05,]$gene

DefaultAssay(marMOp.atac) <- "chromvar"
tf.markers <- FindAllMarkers(
  object = marMOp.atac,
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'nCount_ATAC'
)

mar.tf.mark <- tf.markers[tf.markers$p_val < 0.05,]$gene

common.TFBS <- unique(c(hu.tf.mark,mar.tf.mark))
length(common.TFBS) #579

DefaultAssay(MOp.atac) <- "chromvar"
ave.tf.hs <- AverageExpression(MOp.atac, assays = "chromvar", features = common.TFBS, slot = "data" )
scaled <- t(scale(t(ave.tf.hs$chromvar)))
scaled <- scale(scaled)
tf.hs <- scaled


DefaultAssay(marMOp.atac) <- "chromvar"
ave.tf.mar <- AverageExpression(marMOp.atac, assays = "chromvar", features = common.TFBS, slot = "data" )
scaled <- t(scale(t(ave.tf.mar$chromvar)))
scaled <- scale(scaled)
tf.mar <- scaled

hu.order <- c("SNCG","VIP","LAMP5","SST","PVALB","L2-3 IT","L5 IT","L6 IT","L6 IT Car3","L5-6 NP",
              "L5 ET","L6 CT","L6b")
mar.order <- c("Sncg","Vip","Lamp5","Sst","Pvalb","L2-3 IT","L5 IT","L6 IT","L6 IT Car3","L5-6 NP",
               "L5 ET","L6 CT","L6b")

library("corrplot")
ave.tf.cor<-cor(tf.hs[,hu.order], tf.mar[,mar.order], use = "complete.obs")
write.table(ave.tf.cor, file="HuMOp_marMOp_scaled_average_TFBS_Activities_cor_matrix.txt", sep = "\t", row.names=TRUE, col.names=TRUE)
col2 <- colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#FFFFFF","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))
corrplot(ave.tf.cor, col = col2(500), order = "original", hclust.method = "ward.D2", cl.lim=c(-1,1))




####
##Corr on expression of TFs with diff active TFBSs
####

all.diff.tf <- unique(c(sapply(strsplit(common.TFBS, "::"), "[", 1),
                        sapply(strsplit(common.TFBS, "::"), "[", 2),
                        sapply(strsplit(common.TFBS, "::"), "[", 3)
))
all.diff.tf <- as.character(all.diff.tf)
all.diff.tf <- gsub("\\(var.2)","",all.diff.tf)
all.diff.tf <- gsub("\\(var.3)","",all.diff.tf)

DefaultAssay(MOp.atac) <- "RNA"
DefaultAssay(marMOp.atac) <- "RNA"
common.tf <- intersect(all.diff.tf, intersect(rownames(MOp.atac),rownames(marMOp.atac)))
#418

ave.tf.hs <- AverageExpression(MOp.atac, assays = "RNA", features = common.tf, slot = "scale.data" )
tf.hs <- ave.tf.hs$RNA

ave.tf.mar <- AverageExpression(marMOp.atac, assays = "RNA", features = common.tf, slot = "scale.data" )
tf.mar <- ave.tf.mar$RNA
common.genes <- intersect(rownames(tf.hs), rownames(tf.mar))
hu.order <- c("SNCG","VIP","LAMP5","SST","PVALB","L2-3 IT","L5 IT","L6 IT","L6 IT Car3","L5-6 NP",
              "L5 ET","L6 CT","L6b")
mar.order <- c("Sncg","Vip","Lamp5","Sst","Pvalb","L2-3 IT","L5 IT","L6 IT","L6 IT Car3","L5-6 NP",
               "L5 ET","L6 CT","L6b")
library("corrplot")
ave.tf.cor<-cor(tf.hs[common.genes,hu.order], tf.mar[common.genes,mar.order])
write.table(ave.tf.cor, file="HuMOp_marMOp_scaled_average_TF_Expression_values_cor_matrix.txt", sep = "\t", row.names=TRUE, col.names=TRUE)
col2 <- colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#FFFFFF","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))
corrplot(ave.tf.cor, col = col2(500), order = "original", hclust.method = "ward.D2", cl.lim=c(-1,1))





####
##Corr on SNARE variant gene expression
####

DefaultAssay(MOp.atac) <- "RNA"
DefaultAssay(marMOp.atac) <- "RNA"
MOp.atac <- FindVariableFeatures(object = MOp.atac, selection.method = "vst", nfeatures = 3000)
marMOp.atac <- FindVariableFeatures(object = marMOp.atac, selection.method = "vst", nfeatures = 3000)
head(VariableFeatures(MOp.atac), 10)

var.genes <- c(VariableFeatures(object = MOp.atac),VariableFeatures(object = marMOp.atac))
common.var.genes <- intersect(var.genes, intersect(rownames(MOp.atac),rownames(marMOp.atac)))
#2729 genes

ave.tf.hs <- AverageExpression(MOp.atac, assays = "RNA", features = common.var.genes, slot = "scale.data" )
tf.hs <- ave.tf.hs$RNA

ave.tf.mar <- AverageExpression(marMOp.atac, assays = "RNA", features = common.var.genes, slot = "scale.data" )
tf.mar <- ave.tf.mar$RNA
common.genes <- intersect(rownames(tf.hs), rownames(tf.mar))
hu.order <- c("SNCG","VIP","LAMP5","SST","PVALB","L2-3 IT","L5 IT","L6 IT","L6 IT Car3","L5-6 NP",
              "L5 ET","L6 CT","L6b")
mar.order <- c("Sncg","Vip","Lamp5","Sst","Pvalb","L2-3 IT","L5 IT","L6 IT","L6 IT Car3","L5-6 NP",
               "L5 ET","L6 CT","L6b")
library("corrplot")
ave.tf.cor<-cor(tf.hs[common.genes,hu.order], tf.mar[common.genes,mar.order])
write.table(ave.tf.cor, file="HuMOp_marMOp_scaled_average_Varient_Gene_Expression_values_cor_matrix.txt", sep = "\t", row.names=TRUE, col.names=TRUE)
col2 <- colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#FFFFFF","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))
corrplot(ave.tf.cor, col = col2(500), order = "original", hclust.method = "ward.D2", cl.lim=c(-1,1))












