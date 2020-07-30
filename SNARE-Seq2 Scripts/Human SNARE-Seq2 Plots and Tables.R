#Code Underlying SNARE-Seq2 M1 Data Analyses 
#Blue B. Lake
#b1lake@eng.ucsd.edu

# Human SNARE-Seq2 Plots and Tables --------------------------------------------------
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
set.seed(1234)


#Load Seurat object (MOp.atac)
load("~/NeMO_analysis_folder/SNARE/Analysis/Human/Zhang_BICCN-H_20190523-20190611_huMOp_Seurat.rda")

#Prepare count files
DefaultAssay(MOp.atac) <- "RNA"
rna.counts <- GetAssayData(object = MOp.atac, slot = "counts")
saveRDS(rna.counts, file="Zhang_BICCN-H_20190523-20190611_huMOp_Final_RNA_Counts.RDS")

DefaultAssay(MOp.atac) <- "ATAC"
atac.counts <- GetAssayData(object = MOp.atac, slot = "counts")
saveRDS(atac.counts, file="Zhang_BICCN-H_20190523-20190611_huMOp_Final_AC_Peaks.RDS")

DefaultAssay(MOp.atac) <- "ACTIVITY"
activity.counts <- GetAssayData(object = MOp.atac, slot = "counts")
saveRDS(activity.counts, file="Zhang_BICCN-H_20190523-20190611_huMOp_Final_Cicero_Activities.RDS")

write.table(MOp.atac@meta.data, file="Zhang_BICCN-H_20190523-20190611_huMOp_Final_Sample_Metadata.txt", sep = "\t", row.names=TRUE, col.names=TRUE)





####
##UMAP Plots
####

fac2col <- function(x,s=0.6,v=0.7,shuffle=TRUE,min.group.size=1,return.level.colors=F,unclassified.cell.color='gray80',col=NULL) {
  x <- as.factor(x);
  if(min.group.size>1) {
    x <- factor(x,exclude=levels(x)[unlist(tapply(rep(1,length(x)),x,length))<min.group.size])
  }
  if(!is.null(col)) { 
    if(length(col)<length(levels(x))) stop("fewer colors supplied by col parameter then levels in the factor")
  } else {
    col <- rainbow(length(levels(x)),s=s,v=v);
  }
  if(shuffle) col <- sample(col);
  if(return.level.colors) { names(col) <- levels(x); return(col); }
  y <- col[as.integer(x)]; names(y) <- names(x);
  y[is.na(y)] <- unclassified.cell.color;
  y
}


#AC-level Clusters - Figure 4a, Extended Data Fig. 4k, Extended Data Fig. 5a
Idents(object = MOp.atac) <- "AC_cluster_label"
atac.int.cols<-as.character(unique(meta$AC_cluster_color))
names(atac.int.cols)<-unique(meta$AC_cluster_label)
Idents(object = MOp.atac) <- factor(Idents(object = MOp.atac), levels = unique(meta$AC_cluster_label))

DimPlot(MOp.atac, reduction = "umap", label = TRUE, 
        label.size = 2.5, repel = TRUE) + NoLegend() + ggtitle("AC Clusters"
        ) + scale_color_manual(values = alpha(atac.int.cols, 0.4), name = "AC Clusters"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1)))

DimPlot(MOp.atac, reduction = "r.umap", label = TRUE, 
        label.size = 2.5, repel = TRUE) + NoLegend() + ggtitle("AC Clusters"
        ) + scale_color_manual(values = alpha(atac.int.cols, 0.4), name = "AC Clusters"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1)))

#RNA-level Clusters - Extended Data Fig. 4f
Idents(object = MOp.atac) <- "RNA_cluster"
meta <- MOp.atac@meta.data
cons.cols<-as.character(unique(meta$RNA_cluster_color))
names(cons.cols)<-unique(meta$RNA_cluster)
Idents(object = MOp.atac) <- factor(Idents(object = MOp.atac), levels = unique(meta$RNA_cluster))


DimPlot(MOp.atac, reduction = "umap", label = TRUE, 
        label.size = 1.6, repel = TRUE) + NoLegend() + ggtitle("RNA-level Clusters"
        ) + scale_color_manual(values = alpha(cons.cols, 0.4), name = "RNA-level Clusters"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1)))

DimPlot(MOp.atac, reduction = "r.umap", label = TRUE, 
        label.size = 1.6, repel = TRUE) + NoLegend() + ggtitle("RNA-level Clusters"
        ) + scale_color_manual(values = alpha(cons.cols, 0.4), name = "RNA-level Clusters"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1)))


#Subclass Clusters
Idents(object = MOp.atac) <- "subclass"
subclass.cols<-as.character(unique(MOp.atac$subclass_color))
names(subclass.cols)<-unique(MOp.atac$subclass)
Idents(object = MOp.atac) <- factor(Idents(object = MOp.atac), levels = unique(MOp.atac$subclass))

DimPlot(MOp.atac, reduction = "umap", label = TRUE, 
        label.size = 4, repel = TRUE) + NoLegend() + ggtitle("Subclasss"
        ) + scale_color_manual(values = alpha(subclass.cols, 0.4), name = "Subclass"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1)))

DimPlot(MOp.atac, reduction = "r.umap", label = TRUE, 
        label.size = 4, repel = TRUE) + NoLegend() + ggtitle("Subclass"
        ) + scale_color_manual(values = alpha(subclass.cols, 0.4), name = "Subclass"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1)))



#Patient - Extended Data Fig. 5a
Idents(object = MOp.atac) <- "patient"
patient.cols<-fac2col(levels(Idents(object = MOp.atac)))
names(patient.cols)<-levels(Idents(object = MOp.atac))

DimPlot(MOp.atac, reduction = "umap", label = FALSE, 
        label.size = 4, repel = TRUE) + ggtitle("Patient"
        ) + scale_color_manual(values = alpha(patient.cols, 0.4), name = "Patient"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1)))

DimPlot(MOp.atac, reduction = "r.umap", label = FALSE, 
        label.size = 4, repel = TRUE) + ggtitle("Patient"
        ) + scale_color_manual(values = alpha(patient.cols, 0.4), name = "Patient"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1)))



#Batch
Idents(object = MOp.atac) <- "library"
batch.cols<-fac2col(levels(Idents(object = MOp.atac)))
names(batch.cols)<-levels(Idents(object = MOp.atac))

DimPlot(MOp.atac, reduction = "umap", label = FALSE, 
        label.size = 4, repel = TRUE) + NoLegend() + ggtitle("Batch"
        ) + scale_color_manual(values = alpha(batch.cols, 0.4), name = "Batch"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1)))

DimPlot(MOp.atac, reduction = "r.umap", label = FALSE, 
        label.size = 4, repel = TRUE) + NoLegend() + ggtitle("Batch"
        ) + scale_color_manual(values = alpha(batch.cols, 0.4), name = "Batch"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1)))





#Experiment
Idents(object = MOp.atac) <- "experiment_short"
exp.cols<-fac2col(levels(Idents(object = MOp.atac)))
names(exp.cols)<-levels(Idents(object = MOp.atac))

DimPlot(MOp.atac, reduction = "umap", label = FALSE, 
        label.size = 4, repel = TRUE) + ggtitle("Experiment"
        ) + scale_color_manual(values = alpha(exp.cols, 0.4), name = "Experiment"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1)))

DimPlot(MOp.atac, reduction = "r.umap", label = FALSE, 
        label.size = 4, repel = TRUE) + ggtitle("Experiment"
        ) + scale_color_manual(values = alpha(exp.cols, 0.4), name = "Experiment"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1)))





#Plot QC stats by AC level cluster - Extended Data Fig. 5b
layout(matrix(c(1,1,2,3,4), nrow = 5, ncol = 1, byrow = TRUE))
barplot(prop.table(table(MOp.atac$experiment_short, MOp.atac$AC_cluster_tree_order), margin = 2), 
        main = "Experiment Proportions", cex.names = 0.5, col = exp.cols)
barplot(prop.table(table(MOp.atac$patient, MOp.atac$AC_cluster_tree_order), margin = 2), 
        main = "Patient Proportions", cex.names = 0.5, col = patient.cols)
batch.entropy<-table(MOp.atac$experiment_short, MOp.atac$AC_cluster_tree_order)
entropy <- function(x) {
  value = entropy::entropy(x ,method='MM',unit='log2')/log2(length(x))
  return(value)
}
batch.entropy<-apply(batch.entropy,2,entropy)
barplot(batch.entropy, col = "gray", main = "Experiment Entropy", cex.names = 0.5, ylim=c(0,1))
barplot(table(MOp.atac$AC_cluster_tree_order), col = "gray", main = "Cluster Size", cex.names = 0.5)


#Violin plots - stats
Idents(object = MOp.atac) <- "AC_cluster_tree_order"
Idents(object = MOp.atac) <- factor(Idents(object = MOp.atac), levels = 1:43)
meta <- MOp.atac@meta.data
meta <- meta[order(meta$AC_cluster_tree_order),]
atac.int.cols<-as.character(unique(meta$AC_cluster_color))
names(atac.int.cols)<-unique(meta$AC_cluster_tree_order)

VlnPlot(object = MOp.atac, features = c("RNA_nUMI", "RNA_nGenes","nCount_ATAC","nFeature_ACTIVITY"), 
        ncol = 1, pt.size = -1, group.by = "AC_cluster_tree_order",cols = atac.int.cols)





####
##Correlation of SNARE-Seq2 AC data with SNARE-Seq2 RNA data 
####

DefaultAssay(MOp.atac) <- "ACTIVITY"
all.activity <- rownames(MOp.atac)

#Correlation of RNA-level clusters - Extended Data Fig. 5e
Idents(object = MOp.atac) <- "RNA_cluster_tree_order"
Idents(object = MOp.atac) <- factor(Idents(object = MOp.atac), levels = 1:128)
DefaultAssay(MOp.atac) <- "RNA"
all.genes <- rownames(MOp.atac)
common.genes <- intersect(all.activity,all.genes)
ave.exp<-AverageExpression(MOp.atac, features = common.genes, 
                           slot = "scale.data", assays = c("ACTIVITY","RNA"))

#marker genes
abi.markers <- read.delim("~/NeMO_analysis_folder/SNARE/Analysis/Human/Lein_Smart-seq_Markers_betaprop1_0-4.txt", header = FALSE)
abi.markers <- intersect(abi.markers$V1, common.genes)

library("corrplot")
ave.cor<-cor(ave.exp$RNA[abi.markers,],ave.exp$ACTIVITY[abi.markers,])
write.table(ave.cor, file="huMOp_Correlation_SNARE-R_SNARE-AC_Consensus-Clusters.txt", sep = "\t", row.names=TRUE, col.names=TRUE)

col2 <- colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#FFFFFF","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))
corrplot(ave.cor, col = col2(500),  method="color",order = "original",
         hclust.method = "ward.D", cl.lim=c(-1,1))



#Correlation of AC-level clusters - Extended Data Fig. 5f
Idents(object = MOp.atac) <- "AC_cluster_tree_order"
Idents(object = MOp.atac) <- factor(Idents(object = MOp.atac), levels = 1:43)
ave.exp<-AverageExpression(MOp.atac, features = common.genes, 
                           slot = "scale.data", assays = c("ACTIVITY","RNA"))

ave.cor<-cor(ave.exp$RNA[abi.markers,],ave.exp$ACTIVITY[abi.markers,])
write.table(ave.cor, file="huMOp_Correlation_SNARE-R_SNARE-AC_AC-Clusters.txt", sep = "\t", row.names=TRUE, col.names=TRUE)
col2 <- colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#FFFFFF","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))
corrplot(ave.cor, col = col2(500),  method="color",order = "original",
         hclust.method = "ward.D", cl.lim=c(-1,1))




####
##Coembedding plots 
####

load("Zhang_BICCN-H_20190523-20190611_MOp_RNA-AC_Coembed_Seurat_Final.rda")
fac2col <- function(x,s=0.7,v=0.8,shuffle=TRUE,min.group.size=1,return.level.colors=F,unclassified.cell.color='gray80',col=NULL) {
  x <- as.factor(x);
  if(min.group.size>1) {
    x <- factor(x,exclude=levels(x)[unlist(tapply(rep(1,length(x)),x,length))<min.group.size])
  }
  if(!is.null(col)) { 
    if(length(col)<length(levels(x))) stop("fewer colors supplied by col parameter then levels in the factor")
  } else {
    col <- rainbow(length(levels(x)),s=s,v=v);
  }
  if(shuffle) col <- sample(col);
  if(return.level.colors) { names(col) <- levels(x); return(col); }
  y <- col[as.integer(x)]; names(y) <- names(x);
  y[is.na(y)] <- unclassified.cell.color;
  y
}

res4<-coembed$RNA-AC_snn_res.4.merged
names(res4) <- rownames(coembed@meta.data)
res4.cols<-fac2col(levels(factor(res4)))
names(res4.cols)<-levels(factor(res4))

#Extended Data Fig. 4i
DimPlot(coembed, reduction = "umap", group.by = "RNA-AC_snn_res.4.merged", label = TRUE, 
        label.size = 4, repel = TRUE) + NoLegend() + ggtitle("SNARE-AC"
        ) + scale_color_manual(values = alpha(res4.cols, 0.4), name = "Res4 Clusters") 

tech<-coembed$tech
names(tech) <- rownames(coembed@meta.data)
tech.cols<-fac2col(levels(factor(tech)))
names(tech.cols)<-levels(factor(tech))

#Extended Data Fig. 4h
DimPlot(coembed, reduction = "umap", group.by = "tech", label = FALSE, 
        label.size = 4, repel = TRUE) + ggtitle("SNARE-AC"
        ) + scale_color_manual(values = alpha(tech.cols, 0.4), name = "Technology"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1)))

FeaturePlot(coembed, features = "prediction.score.max")


prediction.table <- coembed@meta.data[coembed@meta.data$tech == "atac", 
                                      c("subclass.predicted.id","subclass.prediction.score.max","cons-cl.predicted.id",
                                        "cons-cl.prediction.score.max")]

write.table(prediction.table, file="huMOp_RNA-AC_Coembed_Prediction_Table.txt", sep = "\t", row.names=TRUE, col.names=TRUE)


#Extended Data Fig. 4g - Consensus cluster predictions
hist(prediction.table$`cons-cl.prediction.score.max`)
abline(v = 0.5, col = "red")

#Extended Data Fig. 4g - Subclass predictions
hist(prediction.table$subclass.prediction.score.max)
abline(v = 0.5, col = "red")


###Extended Data Fig. 4j
#Compare RNA-AC joint embedding clusters vs RNA-level clusters
meta <- coembed@meta.data
meta.AC <- meta[meta$tech == "atac",]
final<-meta.AC$RNA_cluster_ID
names(final)<-rownames(meta.AC)
final <- factor(final)
final <- na.omit(final)
final <- factor(final)
final.df <- data.frame(cluster_id = unique(final),
                       cluster_label = paste0("RNA-level_",unique(final)),
                       cluster_color = rainbow(length(unique(final))))
rownames(final.df) <- final.df$cluster_id
final.df.order<-final.df[order(final.df$cluster_id),]

prop <- meta.AC$`RNA-AC_snn_res.4.merged`
names(prop) <- rownames(meta.AC)
prop <- factor(prop)

library(scrattch.hicat)
compare.result <- compare_annotate(prop, final, final.df.order, reorder = TRUE)

compare.result$g + scale_x_discrete(name ="Co-Embedded Clusters", labels=compare.result$cl.id.map$old)





