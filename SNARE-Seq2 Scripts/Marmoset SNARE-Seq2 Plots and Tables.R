#Code Underlying SNARE-Seq2 M1 Data Analyses 
#Blue B. Lake
#b1lake@eng.ucsd.edu

# Marmoset SNARE-Seq2 Plots and Tables --------------------------------------------------
library(Signac)
library(Seurat)
library(dplyr)
library(Matrix)
library(methods)
library(swne)
library(cicero)
library(ggplot2)
set.seed(1234)

load("~/NeMO_analysis_folder/SNARE/Analysis/Marmoset/Zhang_BICCN-H_20190730_20190903_marMOp_Seurat.rda")

DefaultAssay(marMOp.atac) <- "RNA"
rna.counts <- GetAssayData(object = marMOp.atac, slot = "counts")
saveRDS(rna.counts, file="Zhang_BICCN-H_20190730_20190903_marMOp_Final_RNA_Counts.RDS")

DefaultAssay(marMOp.atac) <- "ATAC"
atac.counts <- GetAssayData(object = marMOp.atac, slot = "counts")
saveRDS(atac.counts, file="Zhang_BICCN-H_20190730_20190903_marMOp_Final_AC_Peaks.RDS")

write.table(marMOp.atac@meta.data, file="Zhang_BICCN-H_20190730_20190903_marMOp_Final_Sample_Metadata.txt", sep = "\t", row.names=TRUE, col.names=TRUE)



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


#AC-level Clusters - Extended Data Fig. 5c
Idents(object = marMOp.atac) <- "AC_cluster_label"
atac.int.cols<-unique(marMOp.atac$AC_cluster_color)
names(atac.int.cols)<-unique(marMOp.atac$AC_cluster_label)
Idents(object = marMOp.atac) <- factor(Idents(object = marMOp.atac), levels = unique(marMOp.atac$AC_cluster_label))

DimPlot(marMOp.atac, reduction = "umap", label = TRUE, pt.size = 0.4,
        label.size = 2.5, repel = TRUE) + NoLegend() + ggtitle("AC Clusters"
        ) + scale_color_manual(values = alpha(atac.int.cols, 0.5), name = "AC Clusters"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1)))

DimPlot(marMOp.atac, reduction = "r.umap", label = TRUE, pt.size = 0.4, 
        label.size = 2.5, repel = TRUE) + NoLegend() + ggtitle("AC Clusters"
        ) + scale_color_manual(values = alpha(atac.int.cols, 0.5), name = "AC Clusters"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1)))


#Subclass Clusters
Idents(object = marMOp.atac) <- "subclass"
order <- c("Lamp5", "Sncg", "Meis2", "Vip", "Sst Chodl", "Sst", "Pvalb", 
           "L2-3 IT", "L5 IT", "L5 ET", "L5-6 NP", "L6 IT", "L6 IT Car3", "L6 CT", "L6b",                                          
           "OPC", "Astro", "Oligo", "Micro-PVM", "VLMC","Peri", "Endo")
marMOp.atac <- subset(marMOp.atac, idents = "Meis2", invert = TRUE)

subclass.cols<-unique(marMOp.atac$subclass_color)
names(subclass.cols)<-unique(marMOp.atac$subclass)
Idents(object = marMOp.atac) <- factor(Idents(object = marMOp.atac), levels = unique(marMOp.atac$subclass))

DimPlot(marMOp.atac, reduction = "umap", label = TRUE, pt.size = 0.4, 
        label.size = 4, repel = TRUE) + NoLegend() + ggtitle("Subclasss"
        ) + scale_color_manual(values = alpha(subclass.cols, 0.5), name = "Subclass"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1)))

DimPlot(marMOp.atac, reduction = "r.umap", label = TRUE, pt.size = 0.4, 
        label.size = 4, repel = TRUE) + NoLegend() + ggtitle("Subclass"
        ) + scale_color_manual(values = alpha(subclass.cols, 0.5), name = "Subclass"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1)))


#Patient - Extended Data Fig. 5c
Idents(object = marMOp.atac) <- "patient"
patient.cols<-fac2col(levels(Idents(object = marMOp.atac)))
names(patient.cols)<-levels(Idents(object = marMOp.atac))

DimPlot(marMOp.atac, reduction = "umap", label = FALSE, pt.size = 0.4, 
        label.size = 4, repel = TRUE) + ggtitle("Patient"
        ) + scale_color_manual(values = alpha(patient.cols, 0.5), name = "Patient"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1)))

DimPlot(marMOp.atac, reduction = "r.umap", label = FALSE, pt.size = 0.4, 
        label.size = 4, repel = TRUE) + ggtitle("Patient"
        ) + scale_color_manual(values = alpha(patient.cols, 0.5), name = "Patient"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1)))


#Batch
Idents(object = marMOp.atac) <- "library"
batch.cols<-fac2col(levels(Idents(object = marMOp.atac)))
names(batch.cols)<-levels(Idents(object = marMOp.atac))

DimPlot(marMOp.atac, reduction = "umap", label = FALSE, pt.size = 0.4, 
        label.size = 4, repel = TRUE) + NoLegend() + ggtitle("Batch"
        ) + scale_color_manual(values = alpha(batch.cols, 0.5), name = "Batch"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1)))

DimPlot(marMOp.atac, reduction = "r.umap", label = FALSE, pt.size = 0.4, 
        label.size = 4, repel = TRUE) + NoLegend() + ggtitle("Batch"
        ) + scale_color_manual(values = alpha(batch.cols, 0.5), name = "Batch"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1)))





####
##Plot QC statistics by cluster
####

#Extended Data Fig. 5d
layout(matrix(c(1,1,2,3,4), nrow = 5, ncol = 1, byrow = TRUE))
barplot(prop.table(table(marMOp.atac$library, marMOp.atac$AC_cluster), margin = 2), 
        main = "Library Proportions", cex.names = 0.5, col = batch.cols)
barplot(prop.table(table(marMOp.atac$patient, marMOp.atac$AC_cluster), margin = 2), 
        main = "Patient Proportions", cex.names = 0.5, col = patient.cols)
batch.entropy<-table(marMOp.atac$library, marMOp.atac$AC_cluster)
entropy <- function(x) {
  value = entropy::entropy(x ,method='MM',unit='log2')/log2(length(x))
  return(value)
}
batch.entropy<-apply(batch.entropy,2,entropy)
barplot(batch.entropy, col = "gray", main = "Library Entropy", cex.names = 0.5, ylim=c(0,1))
barplot(table(marMOp.atac$AC_cluster), col = "gray", main = "Cluster Size", cex.names = 0.5)


#violin plot of stats
DefaultAssay(marMOp.atac) <- 'ATAC'
Idents(object = marMOp.atac) <- "AC_cluster"
Idents(object = marMOp.atac) <- factor(Idents(object = marMOp.atac), levels = 1:24)
meta <- marMOp.atac@meta.data
meta <- meta[order(meta$AC_cluster),]
atac.int.cols<-as.character(unique(meta$AC_cluster_color))
names(atac.int.cols)<-unique(meta$AC_cluster)

VlnPlot(object = marMOp.atac, features = c("nCount_RNA", "nFeature_RNA","nCount_ATAC","nFeature_ACTIVITY"), 
        ncol = 1, pt.size = -1, group.by = "AC_cluster",cols = atac.int.cols)




