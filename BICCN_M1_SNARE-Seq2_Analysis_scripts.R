#Code Underlying SNARE-Seq2 M1 Data Analyses 
#Blue B. Lake
#b1lake@eng.ucsd.edu

# Required Packages and Citations -----------------------------------------
library(DropletUtils)
#Lun ATL, Riesenfeld S, Andrews T, Dao T, Gomes T, participants in the 1st Human Cell
#Atlas Jamboree, Marioni JC (2019). "EmptyDrops: distinguishing cells from empty droplets
#in droplet-based single-cell RNA sequencing data." _Genome Biol._, *20*, 63. doi:
#10.1186/s13059-019-1662-y (URL: https://doi.org/10.1186/s13059-019-1662-y).

#Griffiths JA, Richard AC, Bach K, Lun ATL, Marioni JC (2018). "Detection and removal of
#barcode swapping in single-cell RNA-seq data." _Nat. Commun._, *9*(1), 2667. doi:
#10.1038/s41467-018-05083-x (URL: https://doi.org/10.1038/s41467-018-05083-x).

library(Matrix)
#Douglas Bates and Martin Maechler (2019). Matrix: Sparse and Dense Matrix Classes and
#Methods. R package version 1.2-17. https://CRAN.R-project.org/package=Matrix

library(Seurat)
#Butler, A., Hoffman, P., Smibert, P. et al. Integrating single-cell transcriptomic data 
#across different conditions, technologies, and species. Nat Biotechnol 36, 411–420 (2018).
#https://doi.org/10.1038/nbt.4096

#Tim Stuart, Andrew Butler, Paul Hoffman, Christoph Hafemeister, Efthymia Papalexi, William M.
#Mauck, Yuhan Hao, Marlon Stoeckius, Peter Smibert, Rahul Satija. Comprehensive Integration of
#Single-Cell Data. Cell, Volume 177, Issue 7, 13 June 2019, Pages 1888-1902.e21
#https://satijalab.org/seurat/
  
library(Signac)
#Tim Stuart (2020). Signac: Tools for the Analysis of Single-Cell Chromatin Data.
#https://github.com/timoast/signac, https://satijalab.org/signac.

library(dplyr)
#Hadley Wickham, Romain François, Lionel Henry and Kirill Müller (2019). dplyr: A
#Grammar of Data Manipulation. R package version 0.8.3.
#https://CRAN.R-project.org/package=dplyr

library(pagoda2)
#Peter Kharchenko and Nikolas Barkas (2019). pagoda2: Single Cell Analysis and
#Differential Expression. R package version 0.1.0.
#https://github.com/hms-dbmi/pagoda2

require(parallel)
#R Core Team (2019). R: A language and environment for statistical computing. R
#Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.

library(methods)
#R Core Team (2019). R: A language and environment for statistical computing. R
#Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.

library(stringr)
#Hadley Wickham (2019). stringr: Simple, Consistent Wrappers for Common String
#Operations. R package version 1.4.0. https://CRAN.R-project.org/package=stringr

library(SummarizedExperiment)
#Martin Morgan, Valerie Obenchain, Jim Hester and Hervé Pagès (2019).
#SummarizedExperiment: SummarizedExperiment container. R package version 1.14.1.

library("corrplot")
#Taiyun Wei and Viliam Simko (2017). R package "corrplot": Visualization of a
#Correlation Matrix (Version 0.84). Available from https://github.com/taiyun/corrplot

library(SnapATAC)
#Rongxin Fang (2019). SnapATAC: Single Nucleus Analysis Package for ATAC-Seq. R package
#version 1.0.0. https://github.com/r3fang/SnapATAC

library(GenomicRanges)
#Lawrence M, Huber W, Pag\`es H, Aboyoun P, Carlson M, et al. (2013) Software for
#Computing and Annotating Genomic Ranges. PLoS Comput Biol 9(8): e1003118.
#doi:10.1371/journal.pcbi.1003118

library(swne)
#Yan Wu (2019). swne: Similarity Weighted Nonnegative Embedding: A method for
#visualizing high dimensional datasets. R package version 0.5.7.

library(ggplot2)
#H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York,
#2016.

library(viridis)
#Simon Garnier (2018). viridis: Default Color Maps from 'matplotlib'. R package version
#0.5.1. https://CRAN.R-project.org/package=viridis

library(cicero)
#Hannah A. Pliner, Jay Shendure & Cole Trapnell et. al. (2018). Cicero Predicts
#cis-Regulatory DNA Interactions from Single-Cell Chromatin Accessibility Data.
#Molecular Cell, 71, 858-871.e8.
#https://cole-trapnell-lab.github.io/cicero-release/docs/

library(ChIPseeker)
#Guangchuang Yu, Li-Gen Wang, and Qing-Yu He. ChIPseeker: an R/Bioconductor package for
#ChIP peak annotation, comparison and visualization. Bioinformatics 2015,
#31(14):2382-2383

library(GenomeInfoDb)
#Sonali Arora, Martin Morgan, Marc Carlson and H. Pagès (2019). GenomeInfoDb: Utilities
#for manipulating chromosome names, including modifying them to follow a particular
#naming style. R package version 1.20.0.

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#Bioconductor Core Team and Bioconductor Package Maintainer (2019).
#TxDb.Hsapiens.UCSC.hg38.knownGene: Annotation package for TxDb object(s). R package
#version 3.4.6.

library(BSgenome.Hsapiens.UCSC.hg38)
#The Bioconductor Dev Team (2015). BSgenome.Hsapiens.UCSC.hg38: Full genome sequences
#for Homo sapiens (UCSC version hg38). R package version 1.4.1.

library(EnsDb.Hsapiens.v75)
#Johannes Rainer (2017). EnsDb.Hsapiens.v75: Ensembl based annotation package. R package
#version 2.99.0.

library('org.Hs.eg.db')
#Marc Carlson (2019). org.Hs.eg.db: Genome wide annotation for Human. R package version
#3.8.2.

library(BSgenome.Cjacchus.UCSC.calJac3)
#The Bioconductor Dev Team (2019). BSgenome.Cjacchus.UCSC.calJac3: Full genome sequences
#for Callithrix jacchus (UCSC version calJac3). R package version 1.4.2.

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
#Bioconductor Core Team and Bioconductor Package Maintainer (2019).
#TxDb.Mmusculus.UCSC.mm10.knownGene: Annotation package for TxDb object(s). R package
#version 3.4.7.

library(BSgenome.Mmusculus.UCSC.mm10)
#The Bioconductor Dev Team (2014). BSgenome.Mmusculus.UCSC.mm10: Full genome sequences
#for Mus musculus (UCSC version mm10). R package version 1.4.0.

library("scrattch.hicat")
#Zizhen Yao, Lucas Graybuck and Trygve Bakken (2019). scrattch.hicat: Hierarchical
#Iterative Clustering Analysis for Transcriptomic data. R package version 0.0.16.

library(Gviz)
#Hahne F, Ivanek R. Visualizing Genomic Data Using Gviz and Bioconductor. Methods Mol
#Biol. 1418:335-51 (2016).

library(JASPAR2018)
#Ge Tan (2017). JASPAR2018: Data package for JASPAR 2018. R package version 1.1.1.
#http://jaspar.genereg.net/
  
library(TFBSTools)
#Tan, G., and Lenhard, B. (2016). TFBSTools: an R/bioconductor package for transcription
#factor binding site analysis. Bioinformatics 32, 1555-1556.

library(chromfunks)
#https://github.com/yanwu2014/chromfunks/


# Section 1: Human SNARE-Seq2 RNA - Preparation of Count Matrices and QC Filtering ---------------------------------

library(DropletUtils)
library(Matrix)

#Combine counts by sample and apply library ID's
merge.sparse = function(listMatrixes) {
  # takes a list of sparse matrixes with different columns and adds them row wise
  
  allRownames <- sort(unique(unlist(lapply(listMatrixes,rownames))))
  allColnames <- as.character(unlist(lapply(listMatrixes,colnames)))
  for (currentMatrix in listMatrixes) {
    newRowLocations <- match(rownames(currentMatrix),allRownames)
    indexes <- which(currentMatrix>0, arr.ind = T)
    newRows <- newRowLocations[indexes[,1]]
    columns <- indexes[,2]
    newMatrix <- sparseMatrix(i=newRows,j=columns, x=currentMatrix@x,
                              dims=c(length(allRownames),max(columns)))
    if (!exists("matrixToReturn")) {
      matrixToReturn <- newMatrix
    }
    else {
      matrixToReturn <- cbind2(matrixToReturn,newMatrix)
    }
  }
  rownames(matrixToReturn) <- allRownames
  colnames(matrixToReturn) <- allColnames
  matrixToReturn  
}
dir = "dir/to/RNA_RDS_Files/"
dtn6.file = "dir/to/R1_dTN6_pairs.txt"
merge.dtn6 = function(rds, library = "MOp1", dtn6.file = dtn6.file){
  dt.barcodes <- read.table(dtn6.file, header = F, stringsAsFactors = F)[,1]
  n6.barcodes <- read.table(dtn6.file, header = F, stringsAsFactors = F)[,2]
  
  n6.dt.map <- dt.barcodes
  names(n6.dt.map) <- n6.barcodes
  is.dt <- sapply(colnames(rds$cm_raw), function(x) substr(x, 17, 24) %in% dt.barcodes)
  
  cell.barcodes <- sapply(colnames(rds$cm_raw), function(x) substr(x, 1, 16))
  dt.barcodes <- sapply(colnames(rds$cm_raw), function(x) {
    x <- substr(x, 17, 24)
    if (x %in% names(n6.dt.map)) {
      return(n6.dt.map[[x]])
    } else {
      return(x)
    }
  })
  
  new.barcodes <- paste0(cell.barcodes, dt.barcodes)
  names(new.barcodes) <- colnames(rds$cm_raw)
  unique.barcodes <- unique(new.barcodes)
  
  combined.matrix <- matrix(0, nrow(rds$cm_raw), length(unique.barcodes))
  rownames(combined.matrix) <- rownames(rds$cm_raw)
  colnames(combined.matrix) <- unique.barcodes
  
  combined.matrix[,new.barcodes[colnames(rds$cm_raw)[is.dt]]] <- as.matrix(rds$cm_raw[,is.dt])
  combined.matrix[,new.barcodes[colnames(rds$cm_raw)[!is.dt]]] <- combined.matrix[,new.barcodes[colnames(rds$cm_raw)[!is.dt]]] + as.matrix(rds$cm_raw[,!is.dt])
  combined.matrix <- as(combined.matrix, "dgCMatrix")
  
  colnames(combined.matrix) <- paste(library, colnames(combined.matrix), sep="_")
  combined.matrix
  
}

###Sample 1 - H18-30-001_sorted
MOp1 <- readRDS(paste0(dir,"BICCN-H_20190523A-SPL.hBICCN_20190523_SPL-R_P1_N717_S1.rds"))
MOp1 <- merge.dtn6(MOp1, library = "MOp1", dtn6.file = dtn6.file)

MOp2 <- readRDS(paste0(dir,"BICCN-H_20190523A-SPL.hBICCN_20190523_SPL-R_P2_N718_S2.rds"))
MOp2 <- merge.dtn6(MOp2, library = "MOp2", dtn6.file = dtn6.file)

MOp3 <- readRDS(paste0(dir,"BICCN-H_20190523A-SPL.hBICCN_20190523_SPL-R_P3_N719_S3.rds"))
MOp3 <- merge.dtn6(MOp3, library = "MOp3", dtn6.file = dtn6.file)

MOp4 <- readRDS(paste0(dir,"BICCN-H_20190523A-SPL.hBICCN_20190523_SPL-R_P4_N720_S4.rds"))
MOp4 <- merge.dtn6(MOp4, library = "MOp4", dtn6.file = dtn6.file)

MOp5 <- readRDS(paste0(dir,"BICCN-H_20190523A-SPL.hBICCN_20190523_SPL-R_P5_N701_S9.rds"))
MOp5 <- merge.dtn6(MOp5, library = "MOp5", dtn6.file = dtn6.file)

MOp6 <- readRDS(paste0(dir,"BICCN-H_20190523A-SPL.hBICCN_20190523_SPL-R_P6_N702_S10.rds"))
MOp6 <- merge.dtn6(MOp6, library = "MOp6", dtn6.file = dtn6.file)

MOp7 <- readRDS(paste0(dir,"BICCN-H_20190523A-SPL.hBICCN_20190523_SPL-R_P7_N703_S11.rds"))
MOp7 <- merge.dtn6(MOp7, library = "MOp7", dtn6.file = dtn6.file)

MOp8 <- readRDS(paste0(dir,"BICCN-H_20190523A-SPL.hBICCN_20190523_SPL-R_P8_N704_S12.rds"))
MOp8 <- merge.dtn6(MOp8, library = "MOp8", dtn6.file = dtn6.file)

MOp.list = mget(ls(pattern = "\\MOp"))
H18.30.001.sorted <- merge.sparse(MOp.list)

saveRDS(H18.30.001.sorted, file = paste0(dir,"BICCN-h_20190523A_H18-30-001_sorted_SPL-R_UMI_counts_dTn6.rds"))

#Run DropletUtils to identify empty barcodes
br.out <- barcodeRanks(H18.30.001.sorted)
plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")
o <- order(br.out$rank)
lines(br.out$rank[o], br.out$fitted[o], col="red")
abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"), 
       legend=c("knee", "inflection"))
set.seed(100)
e.out <- emptyDrops(H18.30.001.sorted)
e.out
is.cell <- e.out$FDR <= 0.01
sum(is.cell, na.rm=TRUE)
plot(e.out$Total, -e.out$LogProb, col=ifelse(is.cell, "red", "black"),
     xlab="Total UMI count", ylab="-Log Probability")

e.out$is.cell <- e.out$FDR <= 0.01
is.cell <- rownames(e.out[e.out$is.cell %in% "TRUE",])

#Empty Barcode Filter
H18.30.001.sorted <- H18.30.001.sorted[rowSums(H18.30.001.sorted)>0, ]
H18.30.001.sorted <- H18.30.001.sorted[, is.cell]

#Identify and remove mitochodrial genes 
mt.genes <- grep("^MT-", rownames(H18.30.001.sorted), value = T)
H18.30.001.sorted.noMT <- H18.30.001.sorted[!(rownames(H18.30.001.sorted) %in% mt.genes),]

saveRDS(H18.30.001.sorted, file = paste0(dir,"BICCN-h_20190523A_H18-30-001_sorted_SPL-R_UMI_counts_dTn6_EmptyBC_Filter.rds"))
saveRDS(H18.30.001.sorted.noMT, file = paste0(dir,"BICCN-h_20190523A_H18-30-001_sorted_SPL-R_UMI_counts_dTn6_EmptyBC_MT-Filter.rds"))
writeMM(H18.30.001.sorted.noMT, file= paste0(dir,"BICCN-h_20190523A_H18-30-001_sorted_SPL-R_UMI_counts_dTn6_EmptyBC_MT-Filter.mtx"))

#process remaining libraries as above


#Doublet Filters
dir1 = "/dir/to/RNA_RDS_Files/"
dir2 = "/dir/to/DoubletDetection/"
dir3 = "/dir/to/MOp_Analyses/"

#Sample 1
countMatrix <- readRDS(paste0(dir1, "BICCN-h_20190523A_H18-30-001_sorted_SPL-R_UMI_counts_dTn6_EmptyBC_MT-Filter.rds"))

#Run DoubletDetection (Python) to identify doublets
doublet.det<-read.delim(paste0(dir2, "BICCN-h_20190523A_H18-30-001_sorted_SPL-R_dTn6_DoubletDetection_doublets.txt"),sep="\t",header=FALSE)
doublet.det<-as.character(doublet.det$V1)
names(doublet.det)<-colnames(countMatrix)

doublets<-names(doublet.det)[doublet.det == "1"]
singlets<-names(doublet.det)[doublet.det == "0"]
countMatrix.dd<-countMatrix[,!colnames(countMatrix) %in% doublets]

saveRDS(countMatrix.dd, file = paste0(dir3,"BICCN-h_20190523A_H18-30-001_sorted_SPL-R_UMI_dTn6_EmptyBC_MT_DD_Filter.rds"))

#continue with remaining samples

#Combine all Samples
merge.sparse = function(listMatrixes) {
  # takes a list of sparse matrixes with different columns and adds them row wise
  
  allRownames <- sort(unique(unlist(lapply(listMatrixes,rownames))))
  allColnames <- as.character(unlist(lapply(listMatrixes,colnames)))
  for (currentMatrix in listMatrixes) {
    newRowLocations <- match(rownames(currentMatrix),allRownames)
    indexes <- which(currentMatrix>0, arr.ind = T)
    newRows <- newRowLocations[indexes[,1]]
    columns <- indexes[,2]
    newMatrix <- sparseMatrix(i=newRows,j=columns, x=currentMatrix@x,
                              dims=c(length(allRownames),max(columns)))
    if (!exists("matrixToReturn")) {
      matrixToReturn <- newMatrix
    }
    else {
      matrixToReturn <- cbind2(matrixToReturn,newMatrix)
    }
  }
  rownames(matrixToReturn) <- allRownames
  colnames(matrixToReturn) <- allColnames
  matrixToReturn  
}

MOp.s1 <- readRDS(paste0(dir3, "BICCN-h_20190523A_H18-30-001_sorted_SPL-R_UMI_dTn6_EmptyBC_MT_DD_Filter.rds"))
MOp.s2 <- readRDS(paste0(dir3, "BICCN-h_20190523B_H18-30-001_unsorted_SPL-R_UMI_dTn6_EmptyBC_MT_DD_Filter.rds"))
MOp.s3 <- readRDS(paste0(dir3, "BICCN-h_20190523C_H18-30-002_sorted_SPL-R_UMI_dTn6_EmptyBC_MT_DD_Filter.rds"))
MOp.s4 <- readRDS(paste0(dir3, "BICCN-h_20190523D_H18-30-002_unsorted_SPL-R_UMI_dTn6_EmptyBC_MT_DD_Filter.rds"))
MOp.s5 <- readRDS(paste0(dir3, "BICCN-h_20190611A_H18-30-001_sorted_SPL-R_UMI_dTn6_EmptyBC_MT_DD_Filter.rds"))
MOp.s6 <- readRDS(paste0(dir3, "BICCN-h_20190611B_H18-30-002_sorted_SPL-R_UMI_dTn6_EmptyBC_MT_DD_Filter.rds"))


MOp.list = mget(ls(pattern = "\\MOp"))
MOp <- merge.sparse(MOp.list)
saveRDS(MOp, file = paste0(dir3, "BICCN-h_hMOp_SPL-R_UMI_dTn6_EmptyBC_MT_DD_Filter_07252019.rds"))

#Create Seurat Object and Apply Gene Molecule filters
library(Seurat)
library(dplyr)
library(Matrix)
library(pagoda2)
require(parallel)

MOp <- CreateSeuratObject(counts = MOp, project = "Human M1", min.cells = 3, min.features = 200)

#Gene Filter >200 and <7500 genes
MOp <- subset(x = MOp, subset = nFeature_RNA > 200 & nFeature_RNA < 7500)

#Gene/UMI Ratio Filter to remove low qualty barcodes
countMatrix <- GetAssayData(object = MOp, slot = "counts")
MOp.gmcf <- gene.vs.molecule.cell.filter(countMatrix,min.cell.size=200)
MOp <- subset(MOp, cells = colnames(MOp.gmcf))

saveRDS(MOp, file = "BICCN-h_hMOp_SPL-R_UMI_dTn6_EmptyBC_MT_DD_Gene_Filter_Seurat_v3.RDS")


# Section 2: Human SNARE-Seq2 RNA - Seurat Setup and Pagoda Clustering  ------------------------------------------------------------

#Setup Seurat Object
library(Seurat)
library(dplyr)
library(Matrix)

MOp <- readRDS("BICCN-h_hMOp_SPL-R_UMI_dTn6_EmptyBC_MT_DD_Gene_Filter_Seurat_v3.RDS")

#Normalizing the data
MOp <- NormalizeData(object = MOp, normalization.method = "LogNormalize", 
                     scale.factor = 10000)

#Detection of variable genes across the single cells
MOp <- FindVariableFeatures(object = MOp, selection.method = "vst", nfeatures = 3000)

#Scaling the data and removing unwanted sources of variation
all.genes <- rownames(x = MOp)
MOp <- ScaleData(object = MOp, features = all.genes, vars.to.regress = c("orig.ident"))

#Perform linear dimensional reduction
MOp <- RunPCA(object = MOp, features = VariableFeatures(object = MOp), npcs = 150)
MOp <- FindNeighbors(object = MOp, dims = 1:75)
MOp <- FindClusters(object = MOp, resolution = 1)
MOp <- RunUMAP(object = MOp, dims = 1:75)


#Pagoda2 clustering
library(pagoda2)
require(parallel)

#Down sample Oligodendroctyes (identified as cluster 41)
all.oli <- WhichCells(MOp, idents = 41)
oli.keep <- WhichCells(MOp, idents = 41, downsample = 5000)
oli.remove <- all.oli[!all.oli %in% oli.keep]
MOp <- subset(x = MOp, cells = oli.remove, invert = TRUE)

#Filtered Count Matrix from Seurat
countMatrix <- GetAssayData(object = MOp, slot = "counts")

## Batch annotations
countMatrix.batch <- as.factor(unlist(lapply(colnames(countMatrix),function(x) unlist(strsplit(x,"_"))[1])))
names(countMatrix.batch) <- colnames(countMatrix)

# Generate a new pagoda2 object
p2 <- Pagoda2$new(x = countMatrix, n.cores = 12, trim=10, batch=countMatrix.batch)

# Adjust variance
p2$adjustVariance(plot = T, gam.k = 10)

# Calculate a PCA reduction with the number of PCs specified by nPCs
# and using only the n.odgenes overdispersed genes 
p2$calculatePcaReduction(nPcs = 75, n.odgenes = 6000, maxit=1000)

# Generate K-nearest neighbour graph
p2$makeKnnGraph(k = 100, type = 'PCA', center = T, weight.type = 'none', n.cores = 12, distance = 'cosine')

# Generate an embedding with tSNE on the basis of the PCA reduction
p2$getEmbedding(type = 'PCA', embeddingType = 'tSNE', verbose=TRUE, perplexity=30)

# Identify clusters using the infomap.community method
# on the basis of the reduction called 'PCA' (generated above)
# Save the resulting clustering as 'infomap'
p2$getKnnClusters(method = igraph::infomap.community, type = 'PCA', name = 'infomap')

# Plot
par(mfrow=c(1,1))
p2$plotEmbedding(type = 'PCA',
                 clusterType = 'infomap',
                 embeddingType = 'tSNE',
                 mark.clusters = T,
                 main="infomap",
                 min.group.size=30)

p2$plotEmbedding(type = 'PCA',
                 clusterType = 'walktrap',
                 embeddingType = 'tSNE',
                 mark.clusters = T,
                 main="walktrap",
                 min.group.size=30)


p2$plotEmbedding(type = 'PCA',
                 groups=countMatrix.batch,
                 embeddingType = 'tSNE',
                 mark.clusters = F,
                 main="Batch",
                 mark.cluster.cex = 1,
                 alpha=0.5,show.legend=T) 


p2$plotEmbedding(type='PCA',embeddingType='tSNE',
                 colors=p2$depth,shuffle.colors=F,
                 mark.cluster.cex=1,alpha=0.1,main='depth')

#Repeat makeKnnGraph and getKnnClusters for different k values

k50infomap <- p2$clusters$PCA$infomap
k100infomap <- p2$clusters$PCA$infomap
k200infomap <- p2$clusters$PCA$infomap
k500infomap <- p2$clusters$PCA$infomap


#Transfer Pagoda2 Clusters and PCA values (k = 100) to Seurat
p2 #k = 100
cell.embeddings <- p2$reductions$PCA[rownames(Embeddings(MOp, reduction = "pca")),]
MOp[["pca"]] <- CreateDimReducObject(embeddings = cell.embeddings, key = "pca_", assay = DefaultAssay(MOp))

MOp <- RunUMAP(object = MOp, reduction = "pca", dims = 1:75, n.neighbors = 15L,
               min.dist = 0.4)

MOp[["pagoda_k50_infomap"]] <- k50infomap[rownames(MOp@meta.data)]
MOp[["pagoda_k100_infomap"]] <- k100infomap[rownames(MOp@meta.data)]
MOp[["pagoda_k200_infomap"]] <- k200infomap[rownames(MOp@meta.data)]
MOp[["pagoda_k500_infomap"]] <- k500infomap[rownames(MOp@meta.data)]


#Remove low quality (very low UMI) Exc cluster (cluster 35, k = 200)
Idents(object = MOp) <- "pagoda_k200_infomap"
VlnPlot(MOp, features = "nCount_RNA", pt.size = -1)
MOp <- subset(MOp, idents = 35, invert = TRUE)


##Final annotations 
#Use k = 50 for inhibitory neurons and mic/end clusters
Idents(object = MOp) <- "pagoda_k50_infomap"

all.inh <- rownames(MOp@meta.data[MOp@meta.data$cell_type == "Inh",])
Inh1 <- WhichCells(object = MOp, idents = 23) 
Inh2 <- WhichCells(object = MOp, idents = 84) 
Inh3 <- WhichCells(object = MOp, idents = 18) 
Inh4 <- WhichCells(object = MOp, idents = 55) 
Inh5 <- WhichCells(object = MOp, idents = 36) 
Inh6 <- WhichCells(object = MOp, idents = 62) 
Inh7 <- WhichCells(object = MOp, idents = 57) 
Inh8 <- WhichCells(object = MOp, idents = 61) 
Inh9 <- WhichCells(object = MOp, idents = 59) 
Inh10 <- WhichCells(object = MOp, idents = 79) 
Inh11 <- WhichCells(object = MOp, idents = 69) 
Inh12 <- WhichCells(object = MOp, idents = 67) 
Inh13 <- WhichCells(object = MOp, idents = 78) 
Inh14 <- WhichCells(object = MOp, idents = 58) 
Inh15 <- WhichCells(object = MOp, idents = 45) 
Inh16 <- WhichCells(object = MOp, idents = 42) 
Inh17 <- WhichCells(object = MOp, idents = 63) 
Inh18 <- WhichCells(object = MOp, idents = 37) 
Inh19 <- WhichCells(object = MOp, idents = 22) 

Inh20 <- WhichCells(object = MOp, idents = 81) 
Inh21 <- WhichCells(object = MOp, idents = 15) 
Inh22 <- WhichCells(object = MOp, idents = 20) 
Inh23 <- WhichCells(object = MOp, idents = 40) 
Inh24 <- WhichCells(object = MOp, idents = 60) 
Inh25 <- WhichCells(object = MOp, idents = 68) 
Inh26 <- WhichCells(object = MOp, idents = 70) 
Inh27 <- WhichCells(object = MOp, idents = 46) 
Inh28 <- WhichCells(object = MOp, idents = 56) 
Inh29 <- WhichCells(object = MOp, idents = 89) 
Inh30 <- WhichCells(object = MOp, idents = 64) 
Inh31 <- WhichCells(object = MOp, idents = 35) 
Inh32 <- WhichCells(object = MOp, idents = 27) 
Inh33 <- WhichCells(object = MOp, idents = 48) 
Inh34 <- WhichCells(object = MOp, idents = 16) 
Inh35 <- WhichCells(object = MOp, idents = 80) 
Inh36 <- WhichCells(object = MOp, idents = 87) 
Inh37 <- WhichCells(object = MOp, idents = 28) 
Inh38 <- WhichCells(object = MOp, idents = 24) 

Mic <- WhichCells(object = MOp, idents = 5)

END1 <- WhichCells(object = MOp, idents = 92)
END2 <- WhichCells(object = MOp, idents = 86)
END3 <- WhichCells(object = MOp, idents = 82)
END4 <- WhichCells(object = MOp, idents = 77)
END5 <- WhichCells(object = MOp, idents = 85)
END6 <- WhichCells(object = MOp, idents = 88)


#Use k = 100 for OPC, Oli, Ast
Idents(object = MOp) <- "pagoda_k100_infomap"

OPC <- WhichCells(object = MOp, idents = 2)                    

Ast1 <- WhichCells(object = MOp, idents = 1)                    
Ast2 <- WhichCells(object = MOp, idents = 37)                    
Ast3 <- WhichCells(object = MOp, idents = 44)                    
Ast4 <- WhichCells(object = MOp, idents = 11)                    

Oli1 <- WhichCells(object = MOp, idents = 3)
Oli2 <- WhichCells(object = MOp, idents = 34)


#Use k = 500 for exc L2/3
Idents(object = MOp) <- "pagoda_k500_infomap"

Exc1 <- WhichCells(object = MOp, idents = 7)
Exc2 <- WhichCells(object = MOp, idents = 5)
Exc3 <- WhichCells(object = MOp, idents = 4)
Exc4 <- WhichCells(object = MOp, idents = 12)
Exc5 <- WhichCells(object = MOp, idents = 6)
Exc6 <- WhichCells(object = MOp, idents = 16)
Exc7 <- WhichCells(object = MOp, idents = 26)
Exc8 <- WhichCells(object = MOp, idents = 15)
Exc9 <- WhichCells(object = MOp, idents = 17)


#Use k = 100 for exc L5/6
Idents(object = MOp) <- "pagoda_k100_infomap"

Exc10 <- WhichCells(object = MOp, idents = 32)
Exc11 <- WhichCells(object = MOp, idents = 27)
Exc12 <- WhichCells(object = MOp, idents = 49)
Exc13 <- WhichCells(object = MOp, idents = 15)
Exc14 <- WhichCells(object = MOp, idents = 46)
Exc15 <- WhichCells(object = MOp, idents = 47)
Exc16 <- WhichCells(object = MOp, idents = 48)
Exc17 <- WhichCells(object = MOp, idents = 40)
Exc18 <- WhichCells(object = MOp, idents = 9)
Exc19 <- WhichCells(object = MOp, idents = 41)
Exc20 <- WhichCells(object = MOp, idents = 50)
Exc21 <- WhichCells(object = MOp, idents = 28)

#Rename Clusters
MOp <- SetIdent(MOp, cells = Inh1, value =  1)
MOp <- SetIdent(MOp, cells = Inh2, value =  2)
MOp <- SetIdent(MOp, cells = Inh3, value =  3)
MOp <- SetIdent(MOp, cells = Inh4, value =  4)
MOp <- SetIdent(MOp, cells = Inh5, value =  5)
MOp <- SetIdent(MOp, cells = Inh6, value =  6)
MOp <- SetIdent(MOp, cells = Inh7, value =  7)
MOp <- SetIdent(MOp, cells = Inh8, value =  8)
MOp <- SetIdent(MOp, cells = Inh9, value =  9)
MOp <- SetIdent(MOp, cells = Inh10, value =  10)
MOp <- SetIdent(MOp, cells = Inh11, value =  11)
MOp <- SetIdent(MOp, cells = Inh12, value =  12)
MOp <- SetIdent(MOp, cells = Inh13, value =  13)
MOp <- SetIdent(MOp, cells = Inh14, value =  14)
MOp <- SetIdent(MOp, cells = Inh15, value =  15)
MOp <- SetIdent(MOp, cells = Inh16, value =  16)
MOp <- SetIdent(MOp, cells = Inh17, value =  17)
MOp <- SetIdent(MOp, cells = Inh18, value =  18)
MOp <- SetIdent(MOp, cells = Inh19, value =  19)
MOp <- SetIdent(MOp, cells = Inh20, value =  20)
MOp <- SetIdent(MOp, cells = Inh21, value =  21)
MOp <- SetIdent(MOp, cells = Inh22, value =  22)
MOp <- SetIdent(MOp, cells = Inh23, value =  23)
MOp <- SetIdent(MOp, cells = Inh24, value =  24)
MOp <- SetIdent(MOp, cells = Inh25, value =  25)
MOp <- SetIdent(MOp, cells = Inh26, value =  26)
MOp <- SetIdent(MOp, cells = Inh27, value =  27)
MOp <- SetIdent(MOp, cells = Inh28, value =  28)
MOp <- SetIdent(MOp, cells = Inh29, value =  29)
MOp <- SetIdent(MOp, cells = Inh30, value =  30)
MOp <- SetIdent(MOp, cells = Inh31, value =  31)
MOp <- SetIdent(MOp, cells = Inh32, value =  32)
MOp <- SetIdent(MOp, cells = Inh33, value =  33)
MOp <- SetIdent(MOp, cells = Inh34, value =  34)
MOp <- SetIdent(MOp, cells = Inh35, value =  35)
MOp <- SetIdent(MOp, cells = Inh36, value =  36)
MOp <- SetIdent(MOp, cells = Inh37, value =  37)
MOp <- SetIdent(MOp, cells = Inh38, value =  38)
MOp <- SetIdent(MOp, cells = Exc1, value =  39)
MOp <- SetIdent(MOp, cells = Exc2, value =  40)
MOp <- SetIdent(MOp, cells = Exc3, value =  41)
MOp <- SetIdent(MOp, cells = Exc4, value =  42)
MOp <- SetIdent(MOp, cells = Exc5, value =  43)
MOp <- SetIdent(MOp, cells = Exc6, value =  44)
MOp <- SetIdent(MOp, cells = Exc7, value =  45)
MOp <- SetIdent(MOp, cells = Exc8, value =  46)
MOp <- SetIdent(MOp, cells = Exc9, value =  47)
MOp <- SetIdent(MOp, cells = Exc10, value =  48)
MOp <- SetIdent(MOp, cells = Exc11, value =  49)
MOp <- SetIdent(MOp, cells = Exc12, value =  50)
MOp <- SetIdent(MOp, cells = Exc13, value =  51)
MOp <- SetIdent(MOp, cells = Exc14, value =  52)
MOp <- SetIdent(MOp, cells = Exc15, value =  53)
MOp <- SetIdent(MOp, cells = Exc16, value =  54)
MOp <- SetIdent(MOp, cells = Exc17, value =  55)
MOp <- SetIdent(MOp, cells = Exc18, value =  56)
MOp <- SetIdent(MOp, cells = Exc19, value =  57)
MOp <- SetIdent(MOp, cells = Exc20, value =  58)
MOp <- SetIdent(MOp, cells = Exc21, value =  59)
MOp <- SetIdent(MOp, cells = OPC, value =  60)
MOp <- SetIdent(MOp, cells = AST1, value =  61)
MOp <- SetIdent(MOp, cells = AST2, value =  62)
MOp <- SetIdent(MOp, cells = AST3, value =  63)
MOp <- SetIdent(MOp, cells = AST4, value =  64)
MOp <- SetIdent(MOp, cells = Oli1, value =  65)
MOp <- SetIdent(MOp, cells = Oli2, value =  66)
MOp <- SetIdent(MOp, cells = END1, value =  67)
MOp <- SetIdent(MOp, cells = END2, value =  68)
MOp <- SetIdent(MOp, cells = END3, value =  69)
MOp <- SetIdent(MOp, cells = END4, value =  70)
MOp <- SetIdent(MOp, cells = END5, value =  71)
MOp <- SetIdent(MOp, cells = END6, value =  72)
MOp <- SetIdent(MOp, cells = Mic, value =  73)


to.keep <- c(Inh1,Inh2,Inh3,Inh4,Inh5,Inh6,Inh7,Inh8,Inh9,Inh10,
             Inh11,Inh12,Inh13,Inh14,Inh15,Inh16,Inh17,Inh18,Inh19,Inh20,
             Inh21,Inh22,Inh23,Inh24,Inh25,Inh26,Inh27,Inh28,Inh29,Inh30,
             Inh31,Inh32,Inh33,Inh34,Inh35,Inh36,Inh37,Inh38,
             Exc1,Exc2,Exc3,Exc4,Exc5,Exc6,Exc7,Exc8,Exc9,Exc10,
             Exc11,Exc12,Exc13,Exc14,Exc15,Exc16,Exc17,Exc18,Exc19,Exc20,
             Exc21,OPC,AST1,AST2,AST3,AST4,Oli1,Oli2,END1,END2,
             END3,END4,END5,END6,Mic)

#Remove any cells that were classified to more than one cluster
to.keep <- to.keep[!duplicated(to.keep)]

MOp <- subset(MOp, cells = to.keep)

#check for marker genes
p2 #k = 100
var.info <- p2$misc$varinfo; var.info <- var.info[order(var.info$lp, decreasing = F),];
sn.od.genes <- rownames(var.info)[1:6673]

MOp.markers <- FindAllMarkers(MOp, only.pos = TRUE, features = sn.od.genes,
                              logfc.threshold = 0.25, min.pct = 0.25)




####Compare to SMARTer Clusters

#Prepare Seurat Object for SMARTer Data
genes <- read.csv("SMARTer_nuclei_MOp/20190118_human_MOp_rows-genes.csv",
                  header=TRUE)
abi.data.introns <- read.csv("SMARTer_nuclei_MOp/20190118_human_MOp_intron.counts.csv",
                             header=TRUE)
abi.data.exons <- read.csv("SMARTer_nuclei_MOp/20190118_human_MOp_exon.counts.csv",
                           header=TRUE)
abi.data <- (abi.data.introns + abi.data.exons)
rownames(abi.data) <- genes$gene

abi.anno <- read.csv("SMARTer_nuclei_MOp/cluster.annotation.csv")
abi.clust <- read.csv("SMARTer_nuclei_MOp/cluster.membership.csv")
rownames(abi.clust) <- abi.clust$sample_id
rownames(abi.clust) <- gsub("-",".",rownames(abi.clust))

abi.data <- abi.data[,rownames(abi.clust)]

abi <- CreateSeuratObject(counts = abi.data, project = "ABI Smarter Data", min.cells = 3, min.features = 200)

abi[["cluster"]] <- abi.clust[rownames(abi@meta.data),]$cluster_id

library<-unlist(lapply(rownames(abi@meta.data),function(x) unlist(strsplit(x,"_"))[1])) 
abi[["orig.ident"]]<-library

#Normalize and scale the data
abi <- NormalizeData(abi, normalization.method = "LogNormalize", scale.factor = 10000)
abi <- FindVariableFeatures(abi, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(abi)
abi <- ScaleData(abi, features = all.genes)
abi <- RunPCA(abi, features = VariableFeatures(object = abi), npcs = 50)
Idents(object = abi) <- "cluster"


###Compare SMARTer Clusters with SNARE-Seq2 Pagoda2 Clusters - marker genes
abi.markers <- read.delim("Lein_Smart-seq_Markers_betaprop1_0-4.txt", header = FALSE)
abi.markers <- abi.markers$V1

p2 #k = 100
var.info <- p2$misc$varinfo; var.info <- var.info[order(var.info$lp, decreasing = F),];
sn.od.genes <- rownames(var.info)[1:2000]

common.genes<-c(sn.od.genes, abi.markers)  
common.genes<-common.genes[common.genes %in% rownames(x = abi)]
common.genes<-common.genes[common.genes %in% rownames(x = MOp)]  #5284

ave.MOp<-AverageExpression(MOp, features = common.genes, slot = "scale.data")
ave.abi<-AverageExpression(abi, features = common.genes, slot = "scale.data")
rownames(abi.anno) <- abi.anno$cluster_id
ave.abi$RNA<-ave.abi$RNA[,as.character(1:86)]
colnames(ave.abi$RNA) <- abi.anno$cluster_short_name

library(corrplot)
ave.mop.cor<-cor(cbind(ave.MOp$RNA,ave.abi$RNA))
order.hc2 <- corrMatOrder(ave.mop.cor, order = "hclust", hclust.method = "ward.D")
ave.mop.cor<-ave.mop.cor[order.hc2,order.hc2]
x <- ave.mop.cor[rownames(ave.mop.cor) %in% colnames(ave.MOp$RNA),
                 colnames(ave.mop.cor) %in% colnames(ave.abi$RNA)[1:77]]

col2 <- colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#FFFFFF","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))
corrplot(x, col = col2(500),  method="color",order = "original",
         hclust.method = "ward.D", cl.lim=c(-1,1))






# Section 3: Human SNARE-Seq2 AC - SnapATAC QC Filtering  ---------------------------------
library(SnapATAC)
dir <- "path/to/snap"

file.list = list.files(path = dir, pattern = ".snap", all.files = TRUE,
                       full.names = TRUE)

sample.list = c("MOp1","MOp10","MOp11","MOp12","MOp13","MOp14","MOp15","MOp16","MOp17","MOp18","MOp19",
                "MOp2","MOp20","MOp21","MOp22","MOp23","MOp24","MOp25","MOp26","MOp27","MOp28","MOp29",
                "MOp3","MOp30","MOp31","MOp32","MOp33","MOp34","MOp35","MOp36","MOp37","MOp38","MOp39",
                "MOp4","MOp40","MOp41","MOp42","MOp43","MOp44","MOp45","MOp46","MOp47","MOp48","MOp49",
                "MOp5","MOp50","MOp51","MOp52","MOp53","MOp54","MOp55","MOp56","MOp57","MOp58","MOp59",
                "MOp6","MOp60","MOp61","MOp62","MOp63","MOp64","MOp65","MOp66","MOp7","MOp8","MOp9")

x.sp = createSnap(file=file.list, sample=sample.list)
summarySnap(x.sp)


#Filter by RNA Post-QC cell barcodes
rna.table #RNA Metadata table
rownames(rna.table) <- gsub("MOp","MOP",rownames(rna.table))
rna.matched <- which(x.sp@barcode %in% rownames(rna.table))
x.sp <- x.sp[rna.matched,]
summarySnap(x.sp)

meta <- rna.table
meta <- meta[x.sp@barcode,]

#Add RNA clusters
clust <- meta$cluster
names(clust)<- rownames(meta)
clust <- factor(clust)
x.sp@cluster <- clust


####add in bin matrix
#Barcode selection
x.sp = filterCells(
  obj=x.sp, 
  subset.names=c("fragment.num", "UMI"),
  low.thresholds=c(1000,500),
  high.thresholds=c(Inf, Inf)
)

### bin size selection 

x.sp = addBmatToSnap(
  obj=x.sp, 
  bin.size=5000, 
  num.cores=8
)
calBmatCor(x.sp)


### fragments in promoter ratio and filtering outliers
library(GenomicRanges);
library(ggplot2)
library(dplyr)

promoter.df = read.table("hg38.promoters.bed"); 
promoter.gr = GRanges(promoter.df[,1], IRanges(promoter.df[,2], promoter.df[,3]));
ov = findOverlaps(x.sp@feature, promoter.gr);
idy = queryHits(ov);
promoter_ratio = SnapATAC::rowSums(x.sp[,idy, mat="bmat"], mat="bmat") / SnapATAC::rowSums(x.sp, mat="bmat");

plot(
  x=log(SnapATAC::rowSums(x.sp, mat="bmat") + 1,10), 
  y=promoter_ratio, 
  cex=0.5, 
  col="grey", 
  xlab="log(count)", 
  ylab="FIP Ratio (promoter)",
  ylim=c(0, 1)
)

idx = which(promoter_ratio > 0.1 & promoter_ratio < 0.80)
read_counts = SnapATAC::rowSums(x.sp, mat="bmat")
samples.df <- data.frame(name = x.sp@barcode, 
                         sample = x.sp@sample, 
                         ratio = promoter_ratio, 
                         unique_reads = log(read_counts,10)
)

samples.df <- samples.df %>% 
  group_by(sample) %>% 
  mutate(n = n()) %>% 
  mutate(label = paste0(sample,', N = ',n))

p <- ggplot(samples.df, aes(x=label, y=ratio, fill=sample)) + geom_boxplot() + ggtitle("Promoters - FIP") + ylim(c(0,0.75)) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ylab("")
p

#Filter by FIP
x.sp.filtered <- x.sp[idx,]
x.sp.binary = makeBinary(x.sp.filtered, mat="bmat", outlier.filter=1e-3)
peaks_ov = SnapATAC::rowSums(x.sp.binary[,idy, mat="bmat"], mat="bmat")/SnapATAC::rowSums(x.sp.binary, mat="bmat")

summary(peaks_ov)

plot(
  x=log(SnapATAC::rowSums(x.sp.filtered, mat="bmat") + 1,10), 
  y=peaks_ov, 
  cex=0.5, 
  col="grey", 
  xlab="log(count)", 
  ylab="FIP Ratio (promoter)",
  ylim=c(0, 1)
)

read_counts = SnapATAC::rowSums(x.sp.filtered, mat="bmat")
samples.df <- data.frame(name = x.sp.filtered@barcode, 
                         sample = x.sp.filtered@sample, 
                         ratio = peaks_ov, 
                         unique_reads = log(read_counts,10)
)

samples.df <- samples.df %>% 
  group_by(sample) %>% 
  mutate(n = n()) %>% 
  mutate(label = paste0(sample,', N = ',n))

p <- ggplot(samples.df, aes(x=label, y=ratio, fill=sample)) + geom_boxplot() + ggtitle("Promoters - FIP") + ylim(c(0,0.75)) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ylab("")
p


###Bin filtration -  filter out any bins overlapping with the ENCODE blacklist and bins belonging to chrM

system("wget http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/hg38-human/hg38.blacklist.bed.gz")
black_list <- read.table("hg38.blacklist.bed.gz")
black_list.gr = GRanges(
  black_list[,1], 
  IRanges(black_list[,2], black_list[,3])
)

idy1 = queryHits(findOverlaps(x.sp.binary@feature, black_list.gr))
idy2 = grep("chrM|random", x.sp.binary@feature)
idy = unique(c(idy1, idy2))
x.sp.rmsk = x.sp.binary[,-idy, mat="bmat"];
plotBinCoverage(
  obj=x.sp.binary,
  col="grey",
  border="grey",
  breaks=20,
  xlim=c(-6, 6)
)

x.sp = filterBins(
  x.sp.rmsk,
  low.threshold=-2,
  high.threshold=2,
  mat="bmat"
)

summarySnap(x.sp)


















# Section 4: Human SNARE-Seq2 AC - SnapATAC Peak Calling - round 1 --------------------------------------------

library(SnapATAC)
library(parallel)
library(GenomicRanges);
library(ggplot2)
library(dplyr)


meta #AC metadata table

###Consensus Clusters
clust <- meta$consensus_cluster
names(clust)<- rownames(meta)
clust <- factor(clust)
x.sp@cluster <- clust

#Peak calling
clusters.sel = names(table(x.sp@cluster))[which(table(x.sp@cluster) > 100)];
peaks.ls = mclapply(seq(clusters.sel), function(i){
  print(clusters.sel[i]);
  peaks = runMACS(
    obj=x.sp[which(x.sp@cluster==clusters.sel[i]),], 
    output.prefix=paste0("Cons_clusters.", gsub(" ", "_", clusters.sel)[i]),
    path.to.snaptools="/media/Home_Raid1/b1lake/anaconda2/bin/snaptools",
    path.to.macs="/media/Home_Raid1/b1lake/anaconda2/bin/macs2",
    gsize="hs", # mm, hs, etc
    buffer.size=500, 
    num.cores=6,
    macs.options="--nomodel --shift 100 --ext 200 --qval 5e-2 -B --SPMR",
    tmp.folder=tempdir()
  );
  peaks
}, mc.cores=5);



###Sublass Clusters
clust <- meta$subclass
names(clust)<- rownames(meta)
clust <- factor(clust)
x.sp@cluster <- clust

#Peak calling
clusters.sel = names(table(x.sp@cluster))[which(table(x.sp@cluster) > 50)];
peaks.ls = mclapply(seq(clusters.sel), function(i){
  print(clusters.sel[i]);
  peaks = runMACS(
    obj=x.sp[which(x.sp@cluster==clusters.sel[i]),], 
    output.prefix=paste0("Subclass.", gsub(" ", "_", clusters.sel)[i]),
    path.to.snaptools="/media/Home_Raid1/b1lake/anaconda2/bin/snaptools",
    path.to.macs="/media/Home_Raid1/b1lake/anaconda2/bin/macs2",
    gsize="hs", # mm, hs, etc
    buffer.size=500, 
    num.cores=6,
    macs.options="--nomodel --shift 100 --ext 200 --qval 5e-2 -B --SPMR",
    tmp.folder=tempdir()
  );
  peaks
}, mc.cores=5);


###Class Clusters
clust <- meta$class
names(clust)<- rownames(meta)
clust <- factor(clust)
x.sp@cluster <- clust

#Peak calling
clusters.sel = names(table(x.sp@cluster))[which(table(x.sp@cluster) > 100)];
peaks.ls = mclapply(seq(clusters.sel), function(i){
  print(clusters.sel[i]);
  peaks = runMACS(
    obj=x.sp[which(x.sp@cluster==clusters.sel[i]),], 
    output.prefix=paste0("Class.", gsub(" ", "_", clusters.sel)[i]),
    path.to.snaptools="/media/Home_Raid1/b1lake/anaconda2/bin/snaptools",
    path.to.macs="/media/Home_Raid1/b1lake/anaconda2/bin/macs2",
    gsize="hs", # mm, hs, etc
    buffer.size=500, 
    num.cores=6,
    macs.options="--nomodel --shift 100 --ext 200 --qval 5e-2 -B --SPMR",
    tmp.folder=tempdir()
  );
  peaks
}, mc.cores=5);


###Generate master peak count matrix
peaks.names = system("ls | grep narrowPeak", intern=TRUE);

peak.gr.ls = lapply(peaks.names, function(x){
  peak.df = read.table(x)
  GRanges(peak.df[,1], IRanges(peak.df[,2], peak.df[,3]))
})
peak.gr = reduce(Reduce(c, peak.gr.ls));

x.sp = createPmat(
  obj=x.sp,
  peaks=peak.gr,
  ncell.chunk=20,
  do.par=TRUE,
  num.cores=10
)

# extract peaks x cells matrix from snapATAC
peak.names = paste0(seqnames(peak.gr), ":", start(peak.gr), "-", end(peak.gr))
length(peak.names)

counts <- x.sp@pmat
rownames(counts) <- x.sp@barcode
colnames(counts) <- peak.names

counts = t(counts)
colnames(counts) <- gsub("@", "", colnames(counts))
colnames(counts) <- gsub("#", "_", colnames(counts))



























# Section 5: Human SNARE-Seq2 RNA/AC - Integrative Analysis --------------------------
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)

###make atac seurat using peak file
#load peak by cell barcode counts file
counts

# keep only peaks from chr 1-22 or chrX or chrY
chrom <- sapply(rownames(counts), function(x) strsplit(x, split = ":")[[1]][[1]])
chrom <- sapply(chrom, function(x) strsplit(x, split = "_")[[1]][[1]])
loc <- sapply(rownames(counts), function(x) strsplit(x, split = ":")[[1]][[2]])
rownames(counts) <- paste0(chrom, ":", loc)
chrom.selected <- as.factor(chrom)
chrom.remove <- names(chrom.selected[chrom.selected %in% c("chrUn", "chrM")])
chrom.keep <- setdiff(names(chrom.selected), chrom.remove)
filtered.counts<- counts[which(rownames(counts) %in% chrom.keep),]
counts <- filtered.counts

#create and pre-process atac seurat (MOp.atac)
MOp.atac <- CreateSeuratObject(counts = counts, assay = "ATAC", project = "MOp_ATAC")

meta #metadata tables

MOp.atac <- AddMetaData(MOp.atac, metadata = meta)
MOp.atac$tech <- "atac"

#pre-processing
DefaultAssay(MOp.atac) <- "ATAC"
VariableFeatures(MOp.atac) <- names(which(Matrix::rowSums(MOp.atac) > 100))
MOp.atac <- RunLSI(MOp.atac, n = 50, scale.max = NULL) #latent semantic indexing
MOp.atac <- RunUMAP(MOp.atac, reduction = "lsi", dims = 1:50)

#load SNARE_Seq2 RNA Seurat object (MOp)
MOp
MOp$tech <- "rna"



###Generate gene activity matrix using cisero
library(cicero)

#Cicero analysis: generate predicted gene activity matrix from chrom data
Idents(object = MOp.atac) <- "consensus_cluster"

## Run cicero
hg38.chr.lengths <- read.table("hg38.chr.lengths.txt", header = F, sep = "\t") ## A tab separated text file with chromosome lengths
hg38.chr.lengths[[2]] <- as.numeric(hg38.chr.lengths[[2]])

clusters <- Idents(object = MOp.atac)
count.matrix <- GetAssayData(object = MOp.atac, slot = "counts")
count.matrix@x[count.matrix@x > 0] <- 1 #binarize counts
umap.coordinates <- Embeddings(object = MOp.atac, reduction = "umap")

#update counts matrix
rownames <- rownames(MOp.atac)
chrom <- sapply(rownames, function(x) strsplit(x, split = ":")[[1]][[1]])
loc <- sapply(rownames, function(x) strsplit(x, split = ":")[[1]][[2]])
loc_start <- sapply(loc, function(x) strsplit(x, split = "-")[[1]][[1]])
loc_end <- sapply(loc, function(x) strsplit(x, split = "-")[[1]][[2]])

pData <- data.frame(clusters)
fData <- data.frame(site_name = rownames(count.matrix), chromosome = chrom, bp1 = loc_start, bp2 = loc_end)
input_cds <- newCellDataSet(count.matrix, phenoData = new("AnnotatedDataFrame", data = pData),
                            featureData = new("AnnotatedDataFrame", data = fData),
                            expressionFamily = VGAM::binomialff())

cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = umap.coordinates, k = 50)
conns <- run_cicero(cicero_cds, hg38.chr.lengths) # Takes a few minutes to run

## Assign peaks to modules
ccan.assigns <- generate_ccans(conns)

# Format peak annotation dataframe
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
peaks.file <- makeGRangesFromDataFrame(fData,
                                       keep.extra.columns=TRUE,
                                       ignore.strand=TRUE,
                                       seqinfo=NULL,
                                       seqnames.field= "chromosome",
                                       start.field="bp1",
                                       end.field="bp2",
                                       strand.field="strand",
)
peak_anno <- annotatePeak(peaks.file, tssRegion = c(-10000, 10000),
                          TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,annoDb ="org.Hs.eg.db")

peak_anno_df <- as.data.frame(peak_anno)
rownames(peak_anno_df) <- paste0(peak_anno_df$seqnames, ":", peak_anno_df$start, "-", peak_anno_df$end)
peak_anno_df[!grepl("Promoter|Exon|Intron", peak_anno_df$annotation), "SYMBOL"] <- NA

peak_anno_df <- peak_anno_df[,c("seqnames", "start", "end", "SYMBOL")]
colnames(peak_anno_df) <- c("chromosome", "start", "end", "gene")

peak_anno_df <- peak_anno_df[intersect(rownames(input_cds), rownames(peak_anno_df)),]

## Annotate sites by gene
input_cds <- input_cds[intersect(rownames(input_cds), rownames(peak_anno_df)),]
input_cds <- annotate_cds_by_site(input_cds, peak_anno_df)

## Generate unnormalized gene activity matrix
unnorm_ga <- build_gene_activity_matrix(input_cds, conns, coaccess_cutoff = 0.25)

# remove any rows/columns with all zeroes
unnorm_ga <- unnorm_ga[!Matrix::rowSums(unnorm_ga) == 0, !Matrix::colSums(unnorm_ga) == 0]

#Add activity matrix to new Seurat assay slot 
MOp.atac[["ACTIVITY"]] <- CreateAssayObject(counts = unnorm_ga[,names(clusters)])
DefaultAssay(MOp.atac) <- "ACTIVITY"
MOp.atac <- NormalizeData(object = MOp.atac)
MOp.atac <- ScaleData(object = MOp.atac)

#Find Transfer Anchors using top markers from SMARTer analyses and SNARE-Seq2 RNA
p2 # k = 100
var.info <- p2$misc$varinfo; var.info <- var.info[order(var.info$lp, decreasing = F),];
variable.genes <- rownames(var.info)[1:2000]

DefaultAssay(MOp.atac) <- "ACTIVITY"
abi.markers <- read.delim("Lein_Smart-seq_Markers_betaprop1_0-4.txt", header = FALSE)
abi.markers <- abi.markers$V1
abi.markers <- abi.markers[abi.markers %in% rownames(MOp)]
length(abi.markers) #2492 genes

DefaultAssay(MOp.atac) <- "ATAC"

genes.use <- unique(c(variable.genes,abi.markers)) #3885 genes

transfer.anchors <- FindTransferAnchors(reference = MOp, query = MOp.atac, features = genes.use, 
                                        reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")


#Predict subclass identities in AC data from RNA data
celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = MOp$subclass, 
                                     weight.reduction = MOp.atac[["lsi"]])

MOp.atac <- AddMetaData(MOp.atac, metadata = celltype.predictions[,c(1,22)])

#Plot Predictions for subclass
hist(MOp.atac$prediction.score.max)
abline(v = 0.5, col = "red")


#Co-embedding
refdata <- GetAssayData(
  object = MOp, 
  assay = "RNA", 
  slot = "data"
)

imputation <- TransferData(
  anchorset = transfer.anchors, 
  refdata = refdata, 
  weight.reduction = MOp.atac[["lsi"]]
)

MOp.atac[["RNA"]] <- imputation

coembed <- merge(x = MOp, y = MOp.atac)
variable.genes <- genes.use #3885
coembed <- ScaleData(coembed, features = variable.genes, do.scale = FALSE)
coembed <- RunPCA(coembed, features = variable.genes, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 1:50)

DimPlot(coembed, group.by = "tech")
DimPlot(coembed, group.by = "subclass", label = TRUE, repel = TRUE)
FeaturePlot(coembed, features = "prediction.score.max")

coembed <- FindNeighbors(object = coembed, dims = 1:50)
coembed <- FindClusters(object = coembed, resolution = c(1,2,4))


#Plot predictions for RNA consensus clusters
cl.predictions <- TransferData(anchorset = transfer.anchors, refdata = as.character(MOp$consensus_cluster), 
                               weight.reduction = MOp.atac[["lsi"]])
pred.meta <- cl.predictions[,c(1,130)]
colnames(pred.meta) <- c("cons-cl.predicted.id","cons-cl.prediction.score.max")
MOp.atac <- AddMetaData(MOp.atac, metadata = pred.meta)
hist(MOp.atac$`cons-cl.prediction.score.max`)
abline(v = 0.5, col = "red")




###Compare cluster resolutions
library("scrattch.hicat")
Idents(object = coembed) <- "RNA_snn_res.4"

#Coembed clusters vs consensus clusters (1n clusters)
meta <- coembed@meta.data
meta.AC <- meta[meta$tech == "atac",]
final<-meta.AC$consensus_cluster_number
names(final)<-rownames(meta.AC)
final <- factor(final)
final <- na.omit(final)
final <- factor(final)
final.df <- data.frame(cluster_id = unique(final),
                       cluster_label = paste0("1n_",unique(final)),
                       cluster_color = rainbow(length(unique(final))))
rownames(final.df) <- final.df$cluster_id
final.df.order<-final.df[order(final.df$cluster_id),]

prop<-meta.AC$RNA_snn_res.4
names(prop)<-rownames(meta.AC)
prop <- factor(prop)

compare.result <- compare_annotate(prop, final, final.df.order, reorder = TRUE)
compare.result$g + scale_x_discrete(name ="RNA-AC Integrated Cluster", labels=compare.result$cl.id.map$old)


#Remove low quality clusters and t-cells
Idents(object = coembed) <- "RNA_snn_res.4"

coembed <- subset(coembed, idents = c(48,52,84,62) , invert = TRUE) #clusters showing contribution across cell types
to.remove <- c(na.omit(names(coembed$local_cluster[coembed$local_cluster == 68])),
               na.omit(names(coembed$RNA_cluster[coembed$RNA_cluster == 68])))
coembed <- subset(coembed, cells = to.remove , invert = TRUE) #T-cells

#Merge coembed (res 4) clusters showing no further separation into different RNA clusters
Idents(object = coembed) <- "RNA_snn_res.4"
current.cluster.ids <- c("0","1","10","11","12","13","14","15","16","17","18","19","2","20","21","22","23","24","25","26","27",
                         "28","29","3","30","31","32","33","34","35","36","37","38","39","4","40","41","42","43","44","45","46",
                         "47","49","5","50","51","53","54","55","56","57","58","59","6","60","61","63","64","65","66","67","68",
                         "69","7","70","71","72","73","74","75","76","77","78","79","8","80","81","82","83","85","86","87","88",
                         "89","9","90","91","92","93")
new.cluster.ids <- c("0","1","8","11","12","13","14","15","16","17","8","19","2","20","21","22","23","24","25","20","27",
                     "28","27","3","0","11","32","33","34","35","36","12","38","39","4","40","41","20","43","44","12","4",
                     "47","49","20","50","51","53","54","55","56","57","53","0","6","44","61","63","64","65","66","67","68",
                     "69","7","70","71","47","73","67","75","20","77","78","77","8","80","81","6","83","85","86","87","88",
                     "89","9","90","91","75","0")
Idents(object = coembed) <- plyr::mapvalues(Idents(object = coembed), from = current.cluster.ids, to = new.cluster.ids)
coembed$res4_merged <- Idents(object = coembed)

#Re-compare clusters vs consensus clusters
meta <- coembed@meta.data
meta.AC <- meta[meta$tech == "atac",]
final<-meta.AC$consensus_cluster_number
names(final)<-rownames(meta.AC)
final <- factor(final)
final <- na.omit(final)
final <- factor(final)
final.df <- data.frame(cluster_id = unique(final),
                       cluster_label = paste0("1n_",unique(final)),
                       cluster_color = rainbow(length(unique(final))))
rownames(final.df) <- final.df$cluster_id
final.df.order<-final.df[order(final.df$cluster_id),]

prop <- meta.AC$res4_merged
names(prop) <- rownames(meta.AC)
prop <- factor(prop)

compare.result <- compare_annotate(prop, final, final.df.order, reorder = TRUE)
compare.result$g + scale_x_discrete(name ="RNA-AC Integrated Cluster", labels=compare.result$cl.id.map$old)

#Match RNA/AC integrated clusters with consensus clusters by jaccard similarity and frequency 
comp.1n <- compare.result$tb.df
comp.1n$cl <- as.character(comp.1n$cl)
comp.1n.id <- compare.result$cl.id.map
comp.1n.id$new <- as.character(comp.1n.id$new)
comp.1n <- left_join(comp.1n, comp.1n.id, 
                     by = c("cl" = "new"))
comp.1n$cl <- comp.1n$old
comp.1n %>% group_by(ref.cl) %>% top_n(2, jaccard) -> top2.jacc #select top two consensus clusters by jaccard similarity
top2.jacc %>% group_by(ref.cl) %>% top_n(1, Freq) -> top1 #select final top cluster by frequency
#Using these top matched clusters as a guide, consensus clusters were merged
#taking into account their subclasses, to best match the coembed-level clusters  
#and generate a new AC-level cluster annotation. These were then added to the metadata tables


###Generate Co-embedded Plots

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

res4<-coembed$res4_merged
names(res4) <- rownames(coembed@meta.data)
res4.cols<-fac2col(levels(factor(res4)))
names(res4.cols)<-levels(factor(res4))

DimPlot(coembed, reduction = "umap", group.by = "res4_merged", label = TRUE, 
        label.size = 4, repel = TRUE) + NoLegend() + ggtitle("SNARE-AC"
        ) + scale_color_manual(values = alpha(res4.cols, 0.4), name = "Res4 Clusters") 

tech<-coembed$tech
names(tech) <- rownames(coembed@meta.data)
tech.cols<-fac2col(levels(factor(tech)))
names(tech.cols)<-levels(factor(tech))

DimPlot(coembed, reduction = "umap", group.by = "tech", label = FALSE, 
        label.size = 4, repel = TRUE) + ggtitle("SNARE-AC"
        ) + scale_color_manual(values = alpha(tech.cols, 0.4), name = "Technology"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1)))

FeaturePlot(coembed, features = "prediction.score.max")









# Section 6: Human SNARE-Seq2 AC - SnapATAC Peak Calling - round2 --------------------------------------------

library(SnapATAC)
library(parallel)
library(GenomicRanges);
library(ggplot2)
library(dplyr)


meta #AC metadata table
meta <- meta[x.sp@barcode,]

###AC CLuster Peak calling
clust <- meta$AC_tree_order #AC Level clusters
names(clust)<- rownames(meta)
clust <- factor(clust)
x.sp@cluster <- clust
clusters.sel = names(table(x.sp@cluster))[which(table(x.sp@cluster) > 100)];
peaks.ls = mclapply(seq(clusters.sel), function(i){
  print(clusters.sel[i]);
  peaks = runMACS(
    obj=x.sp[which(x.sp@cluster==clusters.sel[i]),], 
    output.prefix=paste0("AC_clusters.", gsub(" ", "_", clusters.sel)[i]),
    path.to.snaptools="/media/Home_Raid1/b1lake/anaconda2/bin/snaptools",
    path.to.macs="/media/Home_Raid1/b1lake/anaconda2/bin/macs2",
    gsize="hs", # mm, hs, etc
    buffer.size=500, 
    num.cores=6,
    macs.options="--nomodel --shift 100 --ext 200 --qval 5e-2 -B --SPMR",
    tmp.folder=tempdir()
  );
  peaks
}, mc.cores=5);


###Consensus Cluster Peak calling
clust <- meta$consensus_cluster
names(clust)<- rownames(meta)
clust <- factor(clust)
x.sp@cluster <- clust
clusters.sel = names(table(x.sp@cluster))[which(table(x.sp@cluster) > 100)];
peaks.ls = mclapply(seq(clusters.sel), function(i){
  print(clusters.sel[i]);
  peaks = runMACS(
    obj=x.sp[which(x.sp@cluster==clusters.sel[i]),], 
    output.prefix=paste0("Cons_clusters.", gsub(" ", "_", clusters.sel)[i]),
    path.to.snaptools="/media/Home_Raid1/b1lake/anaconda2/bin/snaptools",
    path.to.macs="/media/Home_Raid1/b1lake/anaconda2/bin/macs2",
    gsize="hs", # mm, hs, etc
    buffer.size=500, 
    num.cores=6,
    macs.options="--nomodel --shift 100 --ext 200 --qval 5e-2 -B --SPMR",
    tmp.folder=tempdir()
  );
  peaks
}, mc.cores=5);


###Sublass Cluster Peak calling
clust <- meta$subclass
names(clust)<- rownames(meta)
clust <- factor(clust)
x.sp@cluster <- clust
clusters.sel = names(table(x.sp@cluster))[which(table(x.sp@cluster) > 50)];
peaks.ls = mclapply(seq(clusters.sel), function(i){
  print(clusters.sel[i]);
  peaks = runMACS(
    obj=x.sp[which(x.sp@cluster==clusters.sel[i]),], 
    output.prefix=paste0("Subclass.", gsub(" ", "_", clusters.sel)[i]),
    path.to.snaptools="/media/Home_Raid1/b1lake/anaconda2/bin/snaptools",
    path.to.macs="/media/Home_Raid1/b1lake/anaconda2/bin/macs2",
    gsize="hs", # mm, hs, etc
    buffer.size=500, 
    num.cores=6,
    macs.options="--nomodel --shift 100 --ext 200 --qval 5e-2 -B --SPMR",
    tmp.folder=tempdir()
  );
  peaks
}, mc.cores=5);


###Class Cluster Peak calling
clust <- meta$class
names(clust)<- rownames(meta)
clust <- factor(clust)
x.sp@cluster <- clust
clusters.sel = names(table(x.sp@cluster))[which(table(x.sp@cluster) > 100)];
peaks.ls = mclapply(seq(clusters.sel), function(i){
  print(clusters.sel[i]);
  peaks = runMACS(
    obj=x.sp[which(x.sp@cluster==clusters.sel[i]),], 
    output.prefix=paste0("Class.", gsub(" ", "_", clusters.sel)[i]),
    path.to.snaptools="/media/Home_Raid1/b1lake/anaconda2/bin/snaptools",
    path.to.macs="/media/Home_Raid1/b1lake/anaconda2/bin/macs2",
    gsize="hs", # mm, hs, etc
    buffer.size=500, 
    num.cores=6,
    macs.options="--nomodel --shift 100 --ext 200 --qval 5e-2 -B --SPMR",
    tmp.folder=tempdir()
  );
  peaks
}, mc.cores=5);


###Generate master peak count matrix
peaks.names = system("ls | grep narrowPeak", intern=TRUE);

peak.gr.ls = lapply(peaks.names, function(x){
  peak.df = read.table(x)
  GRanges(peak.df[,1], IRanges(peak.df[,2], peak.df[,3]))
})
peak.gr = reduce(Reduce(c, peak.gr.ls));

x.sp = createPmat(
  obj=x.sp,
  peaks=peak.gr,
  ncell.chunk=20,
  do.par=TRUE,
  num.cores=10
)

# extract peaks x cells matrix from snapATAC
peak.names = paste0(seqnames(peak.gr), ":", start(peak.gr), "-", end(peak.gr))
length(peak.names)

counts <- x.sp@pmat
rownames(counts) <- x.sp@barcode
colnames(counts) <- peak.names

counts = t(counts)
colnames(counts) <- gsub("@", "", colnames(counts))
colnames(counts) <- gsub("#", "_", colnames(counts))


# Section 7: Human SNARE-Seq2 RNA/AC - Generate Final Combined Seurat object ----------------------------------------------------------------

library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)
library(dplyr)
library(Matrix)
library(ggplot2)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(methods)
library(swne)
library(cicero)
set.seed(1234)

#load final combined peak by cell barcode count matrix
counts


###Clean up count matrices
chrom <- sapply(rownames(counts), function(x) strsplit(x, split = ":")[[1]][[1]])
chrom <- sapply(chrom, function(x) strsplit(x, split = "_")[[1]][[1]])
loc <- sapply(rownames(counts), function(x) strsplit(x, split = ":")[[1]][[2]])
rownames(counts) <- paste0(chrom, ":", loc)
chrom.selected <- as.factor(chrom)
chrom.remove <- names(chrom.selected[chrom.selected %in% c("chrUn", "chrM")])
chrom.keep <- setdiff(names(chrom.selected), chrom.remove)
filtered.counts<- counts[which(rownames(counts) %in% chrom.keep),]

counts <- filtered.counts


###re-create and pre-process atac seurat (MOp.atac)
MOp.atac <- CreateSeuratObject(counts = counts, assay = "ATAC", project = "MOp_ATAC")

meta #AC metadata table
meta <- meta[colnames(MOp.atac), ]
MOp.atac <- AddMetaData(MOp.atac, metadata = meta)
MOp.atac$tech <- "atac"

#pre-processing
DefaultAssay(MOp.atac) <- "ATAC"
VariableFeatures(MOp.atac) <- names(which(Matrix::rowSums(MOp.atac) > 100))
MOp.atac <- RunLSI(MOp.atac, n = 50, scale.max = NULL) #latent semantic indexing
MOp.atac <- RunUMAP(MOp.atac, reduction = "lsi", dims = 1:50)


###Add in motif matrix
file = "JASPAR2020_CORE_vertebrates_non-redundant_pfms_jaspar.txt"
pfm <- readJASPARMatrix(file, matrixClass="PFM")

# Scan the DNA sequence of each peak for the presence of each motif
motif.matrix <- CreateMotifMatrix(
  features = StringToGRanges(rownames(MOp.atac), sep = c(":", "-")),
  pwm = pfm,
  genome = 'hg38',
  sep = c(":", "-")
)

# Create a new Mofif object to store the results
motif <- CreateMotifObject(
  data = motif.matrix,
  pwm = pfm
)

# Add the Motif object to the assay
MOp.atac[['ATAC']] <- AddMotifObject(
  object = MOp.atac[['ATAC']],
  motif.object = motif
)

MOp.atac <- RegionStats(
  object = MOp.atac,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  sep = c(":", "-")
)


###Add in Chromvar TF activites
MOp.atac <- RunChromVAR(
  object = MOp.atac,
  genome = BSgenome.Hsapiens.UCSC.hg38
)


###Add in RNA expression data
DefaultAssay(MOp.atac) <- "ATAC"
cells.use = colnames(MOp.atac)
cells.use <- gsub("MOP","MOp", cells.use)
#load RNA seurat object (MOp)
MOp <- subset(MOp, cells = cells.use)
rna.counts <- GetAssayData(object = MOp, slot = "counts")
colnames(rna.counts) <- gsub("MOp","MOP", colnames(rna.counts))
MOp.atac[["RNA"]] <- CreateAssayObject(counts = rna.counts)
DefaultAssay(MOp.atac) <- "RNA"
umap.coordinates <- Embeddings(object = MOp, reduction = "umap")
rownames(umap.coordinates) <- gsub("MOp","MOP", rownames(umap.coordinates))
MOp.atac[["r.umap"]] <- CreateDimReducObject(embeddings = umap.coordinates, key = "r_umap_", assay = DefaultAssay(MOp.atac))
norm.data <- GetAssayData(object = MOp, slot = "data")
colnames(norm.data) <- gsub("MOp","MOP", colnames(norm.data))
MOp.atac <- SetAssayData(object = MOp.atac, slot = "data", new.data = norm.data)
scale.data <- GetAssayData(object = MOp, slot = "scale.data")
colnames(scale.data) <- gsub("MOp","MOP", colnames(scale.data))
MOp.atac <- SetAssayData(object = MOp.atac, slot = "scale.data", new.data = scale.data)


###Add cicero gene activity scores
#generate predicted gene activity matrix from chrom data
Idents(object = MOp.atac) <- "AC_tree_order"
Idents(object = MOp.atac) <- factor(Idents(object = MOp.atac), levels = 1:43)

## Run cicero
hg38.chr.lengths <- read.table("hg38.chr.lengths.txt", header = F, sep = "\t") ## A tab separated text file with chromosome lengths
hg38.chr.lengths[[2]] <- as.numeric(hg38.chr.lengths[[2]])

clusters <- Idents(object = MOp.atac)
count.matrix <- GetAssayData(object = MOp.atac, slot = "counts")
count.matrix@x[count.matrix@x > 0] <- 1 #binarize counts
umap.coordinates <- Embeddings(object = MOp.atac, reduction = "umap")

#update counts matrix
rownames <- rownames(MOp.atac)
chrom <- sapply(rownames, function(x) strsplit(x, split = ":")[[1]][[1]])
loc <- sapply(rownames, function(x) strsplit(x, split = ":")[[1]][[2]])
loc_start <- sapply(loc, function(x) strsplit(x, split = "-")[[1]][[1]])
loc_end <- sapply(loc, function(x) strsplit(x, split = "-")[[1]][[2]])

pData <- data.frame(clusters)
fData <- data.frame(site_name = rownames(count.matrix), chromosome = chrom, bp1 = loc_start, bp2 = loc_end)
input_cds <- newCellDataSet(count.matrix, phenoData = new("AnnotatedDataFrame", data = pData),
                            featureData = new("AnnotatedDataFrame", data = fData),
                            expressionFamily = VGAM::binomialff())

cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = umap.coordinates, k = 50)
conns <- run_cicero(cicero_cds, hg38.chr.lengths) # Takes a few minutes to run

# Assign peaks to modules
ccan.assigns <- generate_ccans(conns)

# Format peak annotation dataframe
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
peaks.file <- makeGRangesFromDataFrame(fData,
                                       keep.extra.columns=TRUE,
                                       ignore.strand=TRUE,
                                       seqinfo=NULL,
                                       seqnames.field= "chromosome",
                                       start.field="bp1",
                                       end.field="bp2",
                                       strand.field="strand",
)
peak_anno <- annotatePeak(peaks.file, tssRegion = c(-10000, 10000),
                          TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,annoDb ="org.Hs.eg.db")

peak_anno_df <- as.data.frame(peak_anno)
rownames(peak_anno_df) <- paste0(peak_anno_df$seqnames, ":", peak_anno_df$start, "-", peak_anno_df$end)
peak_anno_df[!grepl("Promoter|Exon|Intron", peak_anno_df$annotation), "SYMBOL"] <- NA

peak_anno_df <- peak_anno_df[,c("seqnames", "start", "end", "SYMBOL")]
colnames(peak_anno_df) <- c("chromosome", "start", "end", "gene")

peak_anno_df <- peak_anno_df[intersect(rownames(input_cds), rownames(peak_anno_df)),]

# Annotate sites by gene
input_cds <- input_cds[intersect(rownames(input_cds), rownames(peak_anno_df)),]
input_cds <- annotate_cds_by_site(input_cds, peak_anno_df)

# Generate unnormalized gene activity matrix
unnorm_ga <- build_gene_activity_matrix(input_cds, conns, coaccess_cutoff = 0.25)

MOp.atac[["ACTIVITY"]] <- CreateAssayObject(counts = unnorm_ga[,names(clusters)])
DefaultAssay(MOp.atac) <- "ACTIVITY"
MOp.atac <- NormalizeData(object = MOp.atac)
MOp.atac <- ScaleData(MOp.atac)



###Compare SNARE-Seq2 AC with SNARE-Seq2 RNA data
library("corrplot")
#consensus clusters
Idents(object = MOp.atac) <- "consensus_tree_order"
Idents(object = MOp.atac) <- factor(Idents(object = MOp.atac), levels = 1:128)
DefaultAssay(MOp.atac) <- "RNA"
all.genes <- rownames(MOp.atac)
common.genes <- intersect(all.activity,all.genes)
ave.exp<-AverageExpression(MOp.atac, features = common.genes, 
                           slot = "scale.data", assays = c("ACTIVITY","RNA"))

#marker genes
abi.markers <- read.delim("Lein_Smart-seq_Markers_betaprop1_0-4.txt", header = FALSE)
abi.markers <- intersect(abi.markers$V1, common.genes)

x <- cor(ave.exp$RNA[abi.markers,],ave.exp$ACTIVITY[abi.markers,])
col2 <- colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#FFFFFF","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))
corrplot(x, col = col2(500),  method="color",order = "original",
         hclust.method = "ward.D", cl.lim=c(-1,1))


#AC clusters
#marker genes
Idents(object = MOp.atac) <- "AC_tree_order"
Idents(object = MOp.atac) <- factor(Idents(object = MOp.atac), levels = 1:43)
ave.exp<-AverageExpression(MOp.atac, features = common.genes, 
                           slot = "scale.data", assays = c("ACTIVITY","RNA"))

x <- cor(ave.exp$RNA[abi.markers,],ave.exp$ACTIVITY[abi.markers,])
col2 <- colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#FFFFFF","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))
corrplot(x, col = col2(500),  method="color",order = "original",
         hclust.method = "ward.D", cl.lim=c(-1,1))





# Section 8: Human SNARE-Seq2 AC - Identify and visualize DARs ---------------------------------------------
library(Signac)
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(chromfunks)


#DARs for ATAC clusters
Idents(object = MOp.atac) <- "AC_tree_order"
Idents(object = MOp.atac) <- factor(Idents(object = MOp.atac), levels = 1:43)

atac.cl.DARs <- find_all_diff(counts, Idents(object = MOp.atac), min.p = 1e-100)

Top_acDARs <- lapply(atac.cl.DARs, function(df) subset(df, qval < 0.001 & logfc > 1))

Top_acDARs.df <- do.call("rbind", lapply(Top_acDARs, as.data.frame)) 
cl<-unlist(lapply(rownames(Top_acDARs.df),function(x) unlist(strsplit(x,"[.]"))[1]))
loc <- unlist(lapply(rownames(Top_acDARs.df),function(x) unlist(strsplit(x,"[.]"))[2]))
Top_acDARs.df$Cluster <- cl
Top_acDARs.df$Location <- loc
rownames(Top_acDARs.df) <- gsub(" ", "_", rownames(Top_acDARs.df))
Top_acDARs.df

#DARs for subclass clusters
Idents(object = MOp.atac) <- "subclass"
levels(Idents(object = MOp.atac))
order <- c("LAMP5", "SNCG", "VIP", "SST CHODL", "SST", "PVALB", 
           "L2-3 IT", "L5 IT", "L5 ET", "L5-6 NP", "L6 IT", "L6 IT Car3", "L6 CT", "L6b",                                          
           "OPC", "Astro", "Oligo", "Micro-PVM", "VLMC", "Endo")  

Idents(object = MOp.atac) <- factor(Idents(object = MOp.atac), levels = order)

subclass.cl.DARs <- find_all_diff(counts, Idents(object = MOp.atac), min.p = 1e-100)

Top_scDARs <- lapply(subclass.cl.DARs, function(df) subset(df, qval < 0.001 & logfc > 1))

Top_scDARs.df <- do.call("rbind", lapply(Top_scDARs, as.data.frame)) 
cl<-unlist(lapply(rownames(Top_scDARs.df),function(x) unlist(strsplit(x,"[.]"))[1]))
loc <- unlist(lapply(rownames(Top_scDARs.df),function(x) unlist(strsplit(x,"[.]"))[2]))
Top_scDARs.df$Cluster <- cl
Top_scDARs.df$Location <- loc
rownames(Top_scDARs.df) <- gsub(" ", "_", rownames(Top_scDARs.df))
Top_scDARs.df


#DAR Visualizations - AC clusters
Idents(object = MOp.atac) <- "AC_tree_order"
Idents(object = MOp.atac) <- factor(Idents(object = MOp.atac), levels = 1:43)

cl.mark <- Top_acDARs.df[Top_acDARs.df$qval < 0.05,]
cl.mark <- cl.mark[cl.mark$logfc > 1,]
cl.mark <- distinct(cl.mark, Location, .keep_all = TRUE)
cl.mark %>% group_by(Cluster) %>% top_n(5, logfc) -> top5

DotPlot(MOp.atac, features = top5$Location) + RotatedAxis()



# Section 9: Human SNARE-Seq2 RNA/AC - Plots and Tables --------------------------------------------------
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


#UMAP Plots

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


#AC-level Clusters
DefaultAssay(MOp.atac) <- 'ATAC'
Idents(object = MOp.atac) <- "AC_tree_order"
Idents(object = MOp.atac) <- factor(Idents(object = MOp.atac), levels = 1:43)
meta <- MOp.atac@meta.data
meta <- meta[order(meta$AC_tree_order),]
atac.int.cols<-as.character(unique(meta$AC_consensus_color))
names(atac.int.cols)<-unique(meta$AC_tree_order)

DimPlot(MOp.atac, reduction = "umap", label = TRUE, 
        label.size = 4, repel = TRUE) + NoLegend() + ggtitle("AC Clusters"
        ) + scale_color_manual(values = alpha(atac.int.cols, 0.4), name = "AC Clusters"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1)))

DimPlot(MOp.atac, reduction = "r.umap", label = TRUE, 
        label.size = 4, repel = TRUE) + NoLegend() + ggtitle("AC Clusters"
        ) + scale_color_manual(values = alpha(atac.int.cols, 0.4), name = "AC Clusters"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1)))


Idents(object = MOp.atac) <- "AC_consensus_label"
atac.int.cols<-as.character(unique(meta$AC_consensus_color))
names(atac.int.cols)<-unique(meta$AC_consensus_label)
Idents(object = MOp.atac) <- factor(Idents(object = MOp.atac), levels = unique(meta$AC_consensus_label))

DimPlot(MOp.atac, reduction = "umap", label = TRUE, 
        label.size = 2.5, repel = TRUE) + NoLegend() + ggtitle("AC Clusters"
        ) + scale_color_manual(values = alpha(atac.int.cols, 0.4), name = "AC Clusters"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1)))

DimPlot(MOp.atac, reduction = "r.umap", label = TRUE, 
        label.size = 2.5, repel = TRUE) + NoLegend() + ggtitle("AC Clusters"
        ) + scale_color_manual(values = alpha(atac.int.cols, 0.4), name = "AC Clusters"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1)))

#Consensus Clusters
Idents(object = MOp.atac) <- "consensus_cluster"
meta <- MOp.atac@meta.data
cons.cols<-as.character(unique(meta$consensus_color))
names(cons.cols)<-unique(meta$consensus_cluster)
Idents(object = MOp.atac) <- factor(Idents(object = MOp.atac), levels = unique(meta$consensus_cluster))


DimPlot(MOp.atac, reduction = "umap", label = TRUE, 
        label.size = 1.6, repel = TRUE) + NoLegend() + ggtitle("Consensus Clusters"
        ) + scale_color_manual(values = alpha(cons.cols, 0.4), name = "Consensus Clusters"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1)))

DimPlot(MOp.atac, reduction = "r.umap", label = TRUE, 
        label.size = 1.6, repel = TRUE) + NoLegend() + ggtitle("Consensus Clusters"
        ) + scale_color_manual(values = alpha(cons.cols, 0.4), name = "Consensus Clusters"
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



#Patient
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





#Plot QC stats by AC level cluster
layout(matrix(c(1,1,2,3,4), nrow = 5, ncol = 1, byrow = TRUE))
barplot(prop.table(table(MOp.atac$experiment_short, MOp.atac$AC_tree_order), margin = 2), 
        main = "Experiment Proportions", cex.names = 0.5, col = exp.cols)
barplot(prop.table(table(MOp.atac$patient, MOp.atac$AC_tree_order), margin = 2), 
        main = "Patient Proportions", cex.names = 0.5, col = patient.cols)
batch.entropy<-table(MOp.atac$experiment_short, MOp.atac$AC_tree_order)
entropy <- function(x) {
  value = entropy::entropy(x ,method='MM',unit='log2')/log2(length(x))
  return(value)
}
batch.entropy<-apply(batch.entropy,2,entropy)
barplot(batch.entropy, col = "gray", main = "Experiment Entropy", cex.names = 0.5, ylim=c(0,1))
barplot(table(MOp.atac$AC_tree_order), col = "gray", main = "Cluster Size", cex.names = 0.5)


#Violin plots - stats
Idents(object = MOp.atac) <- "AC_tree_order"
Idents(object = MOp.atac) <- factor(Idents(object = MOp.atac), levels = 1:43)
meta <- MOp.atac@meta.data
meta <- meta[order(meta$AC_tree_order),]
atac.int.cols<-as.character(unique(meta$AC_consensus_color))
names(atac.int.cols)<-unique(meta$AC_tree_order)

VlnPlot(object = MOp.atac, features = c("RNA_nUMI", "RNA_nGenes","nCount_ATAC","nFeature_ACTIVITY"), 
        ncol = 1, pt.size = -1, group.by = "AC_tree_order",cols = atac.int.cols)




# Section 10: Human SNARE-Seq2 RNA/AC - Motif CCAN analyses - Subclasses ------------------------------------
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

#load MOp.atac Seurat object

##Identify and plot RNA marker genes by subclass
#RNA Marker genes for Subclass-level clusters
DefaultAssay(MOp.atac) <- "RNA"
Idents(object = MOp.atac) <- "subclass"
order <- c("LAMP5", "SNCG", "VIP", "SST CHODL", "SST", "PVALB", 
           "L2-3 IT", "L5 IT", "L5 ET", "L5-6 NP", "L6 IT", "L6 IT Car3", "L6 CT", "L6b",                                          
           "OPC", "Astro", "Oligo", "Micro-PVM", "VLMC", "Endo")  

Idents(object = MOp.atac) <- factor(Idents(object = MOp.atac), levels = order)
genes.use <- rownames(MOp.atac)
SC.level.markers <- FindAllMarkers(MOp.atac, features = genes.use,
                                   only.pos = TRUE, logfc.threshold = 0.25)
RNA.mark <- SC.level.markers
RNA.mark <- RNA.mark[RNA.mark$p_val_adj < 0.05, ]
RNA.mark <- RNA.mark[RNA.mark$avg_logFC > 1, ]
RNA.mark <- distinct(RNA.mark, gene, .keep_all = TRUE)
RNA.mark %>% group_by(cluster) %>% top_n(10, avg_logFC) -> top10 

ave.RNA <- AverageExpression(MOp.atac, assays = "RNA", features = top10$gene, slot = "data" )

scaled <- t(scale(t(ave.RNA$RNA)))
scaled <- scale(scaled)
range(scaled)
scaled[scaled < 0] <- 0
scaled[scaled > 4] <- 4


ggHeat((scaled), rescaling = "none", clustering = "none", x.lab.size = 8, y.lab.size = 11,
) + scale_fill_viridis()


#Plot associated cicero gene activites
DefaultAssay(MOp.atac) <- "ACTIVITY"
gene.use <- top10$gene[top10$gene %in% rownames(MOp.atac)]
ave.act <- AverageExpression(MOp.atac, assays = "ACTIVITY", features = gene.use, slot = "data" )

scaled <- t(scale(t(ave.act$ACTIVITY)))
scaled <- scale(scaled)
range(scaled)
scaled[scaled < 0] <- 0
scaled[scaled > 4] <- 4

ggHeat((scaled), rescaling = "none", clustering = "none", x.lab.size = 8, y.lab.size = 11,
) + scale_fill_viridis()


#TF motif analysis in subclass DEGs
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




##Marker gene TFBS enrichments

#load conns, ccan.assigns and peak_anno_df objects from Cicero analyses
#load SC.level.markers (subclass level RNA markers)

peak_anno_df <- peak_anno_df[rownames(ccan.assigns),]
peak_anno_df$CCAN <- ccan.assigns$CCAN

#Associate CCANs with subclass Markers
DefaultAssay(MOp.atac) <- "ATAC"
Idents(object = MOp.atac) <- "subclass"
order <- c("LAMP5", "SNCG", "VIP", "SST CHODL", "SST", "PVALB", 
           "L2-3 IT", "L5 IT", "L5 ET", "L5-6 NP", "L6 IT", "L6 IT Car3", "L6 CT", "L6b",                                          
           "OPC", "Astro", "Oligo", "Micro-PVM", "VLMC", "Endo")  
Idents(object = MOp.atac) <- factor(Idents(object = MOp.atac), levels = order)

RNA.mark <- SC.level.markers
RNA.mark <- RNA.mark[RNA.mark$p_val_adj < 0.05, ]
RNA.mark <- RNA.mark[RNA.mark$avg_logFC > 1, ]
RNA.mark <- distinct(RNA.mark, gene, .keep_all = TRUE)

peak_anno_df <- peak_anno_df[peak_anno_df$gene %in% RNA.mark$gene,]
peak_anno_df_sc <- left_join(peak_anno_df, RNA.mark, by = "gene")
rownames(peak_anno_df_sc) <- rownames(peak_anno_df)
peak_anno_df_sc <-na.omit(peak_anno_df_sc)
peak_anno_df_sc$site <- rownames(peak_anno_df_sc)

#Re-order subclasses
peak_anno_df_sc$cluster_order <- factor(peak_anno_df_sc$cluster)
peak_anno_df_sc$cluster_order <- plyr::mapvalues(peak_anno_df_sc$cluster_order,
                                 from = c("LAMP5", "SNCG", "VIP", "SST CHODL", "SST", "PVALB", 
                                 "L2-3 IT", "L5 IT", "L5 ET", "L5-6 NP", "L6 IT Car3", "L6 CT", "L6b",                                          
                                 "OPC", "Astro", "Oligo", "Micro-PVM", "VLMC", "Endo"), 
                                 to = c(1:19))
peak.ccam.mark <-arrange(peak_anno_df_sc, cluster_order, desc(avg_logFC))
rownames(peak.ccam.mark) <- peak.ccam.mark$site

DefaultAssay(MOp.atac) <- "ATAC"

#find gene motif enrichments (genes having >1 site)
tab <- table(peak.ccam.mark$gene)
unique.genes <- names(tab[tab > 1])

gene.tfbs.list <- lapply(unique.genes, function(cl) {
  print(paste("Running for gene:", cl))
  
  cl.sites <- rownames(peak_anno_df[peak_anno_df$gene %in% cl,])
  EM <- FindMotifs(object = MOp.atac,
                   features = cl.sites)
  EM <- EM[EM$pvalue < 0.05,]
  EM
})
names(gene.tfbs.list) <- unique.genes

gene.tfbs <- do.call("rbind", lapply(gene.tfbs.list, as.data.frame)) 

#Add associated cluster (predicted by DEG analysis) for each marker
marker<-unlist(lapply(rownames(gene.tfbs),function(x) unlist(strsplit(x,"[.]"))[1]))
gene.tfbs$gene <- marker
cl.order <- peak.ccam.mark$cluster
names(cl.order) <- peak.ccam.mark$gene
cl.order <- factor(cl.order) 
cl.order <- cl.order[gene.tfbs$gene]
gene.tfbs$cluster <- cl.order


#find active TFBS 
DefaultAssay(MOp.atac) <- "chromvar"
tf.markers <- FindAllMarkers(
  object = MOp.atac,
  features = unique(gene.tfbs$motif),
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'nCount_ATAC'
)
colnames(tf.markers)[colnames(tf.markers) == "gene"] <- "motif"
tf.markers <- left_join(gene.tfbs, tf.markers, by = c("motif","cluster"))
rownames(tf.markers) <- rownames(gene.tfbs)
hu.tf.markers <-na.omit(tf.markers)

#re-order by cluster
hu.tf.markers$rownames <- rownames(hu.tf.markers)
hu.tf.markers$cluster_order <- hu.tf.markers$cluster
hu.tf.markers$cluster_order <- plyr::mapvalues(hu.tf.markers$cluster_order,
                               from = c("LAMP5", "SNCG", "VIP", "SST CHODL", "SST", "PVALB", 
                               "L2-3 IT", "L5 IT", "L5 ET", "L5-6 NP", "L6 IT Car3", "L6 CT", "L6b",                                          
                               "OPC", "Astro", "Oligo", "Micro-PVM", "VLMC", "Endo"), 
                               to = 1:19)
hu.tf.markers$cluster_order <- as.numeric(hu.tf.markers$cluster_order)
hu.tf.mark <-arrange(hu.tf.markers, cluster_order, desc(avg_logFC))
rownames(hu.tf.mark) <- hu.tf.mark$rownames

tf.mark <- hu.tf.mark[hu.tf.mark$p_val_adj < 0.05,]
tf.mark <- tf.mark[tf.mark$avg_logFC > 0.5,]
tf.mark <- distinct(tf.mark, motif, .keep_all = TRUE)
tf.mark %>% group_by(cluster) %>% top_n(10, avg_logFC) -> top10 

sites.use <- top10$motif
ave.sites <- AverageExpression(MOp.atac, assays = "chromvar", features = sites.use, slot = "data" )

names.sites <- top10$motif.name
ave.sites <- ave.sites$chromvar
rownames(ave.sites) <- names.sites

scaled <- t(scale(t(ave.sites)))
scaled <- scale(scaled)
range(scaled)
scaled[scaled < 0] <- 0
scaled[scaled > 4] <- 4


#Plot marker-enriched TFBSs that are active in the subclasses
ggHeat(scaled, rescaling = "none", clustering = "none", x.lab.size = 8, y.lab.size = 11,
) + scale_fill_viridis()




# Section 11: Human SNARE-Seq2 RNA/AC - ChC versus Basket Cells - Motif CCAN analyses ------------------------------------
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

#load MOp.atac object


###Run Cicero only on PVALB neurons
Idents(object = MOp.atac) <- "AC_tree_order"
Idents(object = MOp.atac) <- factor(Idents(object = MOp.atac), levels = 1:43)

MOp.atac.pvalb <- subset(MOp.atac, idents = 14:18)
DefaultAssay(MOp.atac.pvalb) <- "ATAC"

#Cicero analysis: generate predicted gene activity matrix from chrom data
hg38.chr.lengths <- read.table("hg38.chr.lengths.txt", header = F, sep = "\t") ## A tab separated text file with chromosome lengths
hg38.chr.lengths[[2]] <- as.numeric(hg38.chr.lengths[[2]])

clusters <- Idents(object = MOp.atac.pvalb)
count.matrix <- GetAssayData(object = MOp.atac.pvalb, slot = "counts")
count.matrix@x[count.matrix@x > 0] <- 1 #binarize counts
umap.coordinates <- Embeddings(object = MOp.atac.pvalb, reduction = "umap")

#update counts matrix
rownames <- rownames(MOp.atac.pvalb)
chrom <- sapply(rownames, function(x) strsplit(x, split = ":")[[1]][[1]])
loc <- sapply(rownames, function(x) strsplit(x, split = ":")[[1]][[2]])
loc_start <- sapply(loc, function(x) strsplit(x, split = "-")[[1]][[1]])
loc_end <- sapply(loc, function(x) strsplit(x, split = "-")[[1]][[2]])

pData <- data.frame(clusters)
fData <- data.frame(site_name = rownames(count.matrix), chromosome = chrom, bp1 = loc_start, bp2 = loc_end)
input_cds <- newCellDataSet(count.matrix, phenoData = new("AnnotatedDataFrame", data = pData),
                            featureData = new("AnnotatedDataFrame", data = fData),
                            expressionFamily = VGAM::binomialff())

cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = umap.coordinates, k = 50)
conns <- run_cicero(cicero_cds, hg38.chr.lengths) # Takes a few minutes to run

## Assign peaks to modules
ccan.assigns <- generate_ccans(conns)
# "Coaccessibility cutoff used: 0.08"


#annotate ccans with nearby genes
rownames <- rownames(ccan.assigns)
chrom <- sapply(rownames, function(x) strsplit(x, split = ":")[[1]][[1]])
loc <- sapply(rownames, function(x) strsplit(x, split = ":")[[1]][[2]])
loc_start <- sapply(loc, function(x) strsplit(x, split = "-")[[1]][[1]])
loc_end <- sapply(loc, function(x) strsplit(x, split = "-")[[1]][[2]])

ccan.data <- data.frame(site_name = rownames, chromosome = chrom, bp1 = loc_start, bp2 = loc_end)
ccan.assigns 

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
peaks.file <- makeGRangesFromDataFrame(ccan.data,
                                       keep.extra.columns=TRUE,
                                       ignore.strand=TRUE,
                                       seqinfo=NULL,
                                       seqnames.field= "chromosome",
                                       start.field="bp1",
                                       end.field="bp2",
                                       strand.field="strand",
)
peak_anno <- annotatePeak(peaks.file, tssRegion = c(-10000, 10000),
                          TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,annoDb ="org.Hs.eg.db")

peak_anno_df <- as.data.frame(peak_anno)
rownames(peak_anno_df) <- paste0(peak_anno_df$seqnames, ":", peak_anno_df$start, "-", peak_anno_df$end)
peak_anno_df <- peak_anno_df[,c("seqnames", "start", "end", "SYMBOL")]
colnames(peak_anno_df) <- c("chromosome", "start", "end", "gene")

peak_anno_df <- peak_anno_df[rownames(ccan.assigns),]
peak_anno_df$CCAN <- ccan.assigns$CCAN
head(peak_anno_df)

# Annotate sites by gene
input_cds <- input_cds[intersect(rownames(input_cds), rownames(peak_anno_df)),]
input_cds <- annotate_cds_by_site(input_cds, peak_anno_df)
unnorm_ga <- build_gene_activity_matrix(input_cds, conns, coaccess_cutoff = 0.25)

MOp.atac.pvalb[["ACTIVITY"]] <- CreateAssayObject(counts = unnorm_ga[,names(clusters)])
DefaultAssay(MOp.atac.pvalb) <- "ACTIVITY"
MOp.atac.pvalb <- NormalizeData(object = MOp.atac.pvalb)
MOp.atac.pvalb <- ScaleData(MOp.atac.pvalb)

#Associate CCANs with ChC Markers
DefaultAssay(MOp.atac.pvalb) <- "RNA"

#Find all ChC markers
genes.use <- rownames(MOp.atac.pvalb)
ChC.markers <- FindMarkers(MOp.atac.pvalb, 
                           ident.1 = 18,
                           ident.2 = c(14:17),
                           features = genes.use,
                           only.pos = TRUE, logfc.threshold = 0.25)
RNA.mark <- ChC.markers
RNA.mark <- RNA.mark[RNA.mark$p_val_adj < 0.05, ]
RNA.mark <- RNA.mark[RNA.mark$avg_logFC > 0.5, ]

peak_anno_df <- peak_anno_df[peak_anno_df$gene %in% rownames(RNA.mark),]
RNA.mark$gene <- rownames(RNA.mark)
peak_anno_df_pv <- left_join(peak_anno_df, RNA.mark, by = "gene")
rownames(peak_anno_df_pv) <- rownames(peak_anno_df)
peak_anno_df_pv <-na.omit(peak_anno_df_pv)

peak_anno_df_pv$site <- rownames(peak_anno_df_pv)
peak_anno_df_pv

MOp.atac.chand <- subset(MOp.atac.pvalb, idents = 18)
DefaultAssay(MOp.atac.chand) <- "ATAC"

#find gene motif enrichments (genes having >1 site)
tab <- table(peak_anno_df_pv$gene)
unique.genes <- names(tab[tab > 1])

gene.tfbs.list <- lapply(unique.genes, function(cl) {
  print(paste("Running for gene:", cl))
  
  cl.sites <- rownames(peak_anno_df[peak_anno_df$gene %in% cl,])
  EM <- FindMotifs(object = MOp.atac.chand,
                   features = cl.sites)
  EM <- EM[EM$pvalue < 0.05,]
  EM
})
names(gene.tfbs.list) <- unique.genes
gene.tfbs <- do.call("rbind", lapply(gene.tfbs.list, as.data.frame)) 
marker<-unlist(lapply(rownames(gene.tfbs),function(x) unlist(strsplit(x,"[.]"))[1]))
gene.tfbs$gene <- marker

#find active TFBS in ChC cells over remaining PVALB
DefaultAssay(MOp.atac.pvalb) <- "chromvar"
tf.markers <- FindMarkers(
  object = MOp.atac.pvalb, 
  ident.1 = 18,
  ident.2 = c(14:17),
  features = unique(gene.tfbs$motif),
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'nCount_ATAC'
)
tf.markers$motif <- rownames(tf.markers)
tf.markers <- left_join(gene.tfbs, tf.markers, by = c("motif"))
rownames(tf.markers) <- rownames(gene.tfbs)
hu.tf.markers <-na.omit(tf.markers)
hu.tf.markers

tf.mark <- hu.tf.markers[hu.tf.markers$p_val_adj < 0.05,]
tf.mark <- distinct(tf.mark, motif, .keep_all = TRUE)

sites.use <- tf.mark$motif
ave.sites <- AverageExpression(MOp.atac.pvalb, assays = "chromvar", features = sites.use, slot = "data" )

names.sites <- tf.mark$motif.name
ave.sites <- ave.sites$chromvar
rownames(ave.sites) <- names.sites

#Plot ChC marker enriched TFBSs that are active in ChC
ggHeat(ave.sites, rescaling = "row", clustering = "row", x.lab.size = 8, y.lab.size = 11,
) + scale_fill_viridis()





#Feature plots 
DefaultAssay(MOp.atac) <- "RNA"
FeaturePlot(
  object = MOp.atac,
  features = "RORA",
  #min.cutoff = 'q10',
  #max.cutoff = 'q90',
  pt.size = 0.1, reduction = "r.umap",
  cols = c("lightgrey","#586f87","#c43131")
) 

FeaturePlot(
  object = MOp.atac,
  features = "GSG1L",
  #min.cutoff = 'q10',
  #max.cutoff = 'q90',
  pt.size = 0.1, reduction = "r.umap",
  cols = c("lightgrey","#586f87","#c43131")
) 

DefaultAssay(MOp.atac) <- "chromvar"
FeaturePlot(
  object = MOp.atac,
  features = "MA0072.1",
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1, reduction = "r.umap",
  cols = c("lightgrey","#586f87","#c43131")
)

DefaultAssay(MOp.atac) <- 'ATAC'
MotifPlot(
  object = MOp.atac,
  motifs = "MA0072.1",
  assay = 'ATAC'
)




# Section 12: Human SNARE-Seq2 RNA/AC - Connection and Feature Plots -------------------------------------------------------
library(Signac)
library(Seurat)
library(dplyr)
library(Matrix)
library(methods)
library(ggplot2)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library('org.Hs.eg.db')
library(ChIPseeker)
library(cicero)
library(SummarizedExperiment)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
set.seed(1234)

#load MOp.atac seurat object

## Filter connections
filter_conns <- function(conns, min.coaccess) {
  conns$coaccess <- abs(conns$coaccess)
  conns <- subset(conns, coaccess > min.coaccess)
  conns$Peak1 <- as.character(conns$Peak1)
  conns$Peak2 <- as.character(conns$Peak2)
  conns
}

## Load conns file

conns.file <- "path/to/conns"
load(conns.file)
conns <- filter_conns(conns, min.coaccess = 0.1)

## Chromosomes to visualization
#Exc L2/3 CUX2: chr12:111030000-111250000
viz.chr <- "chr12" ## Chromosome
viz.min <- 111030000 ## Start
viz.max <- 111250000 ## End


## Make plot
plot_list <- plot_connections(conns, viz.chr, viz.min, viz.max,
                              coaccess_cutoff = 0,
                              alpha_by_coaccess = T,
                              include_axis_track = T,
                              return_as_list = T)

## Load accessibility data bedGraph tracks and add them to the plot list
dir <- "path/to/bedGraph"

dTrack1 <- DataTrack(range = paste0(dir, "/Subclass.LAMP5_treat_pileup.bedGraph"),
                     genome = "hg38", type = "gradient", chromosome = viz.chr,
                     name = "bedGraph")
plot_list[[4]] <- dTrack1

dTrack2 <- DataTrack(range = paste0(dir, "/Subclass.SNCG_treat_pileup.bedGraph"),
                     genome = "hg38", type = "gradient", chromosome = viz.chr,
                     name = "bedGraph")
plot_list[[5]] <- dTrack2

dTrack3 <- DataTrack(range = paste0(dir, "/Subclass.VIP_treat_pileup.bedGraph"),
                     genome = "hg38", type = "gradient", chromosome = viz.chr,
                     name = "bedGraph")
plot_list[[6]] <- dTrack3

dTrack4 <- DataTrack(range = paste0(dir, "/Subclass.SST_CHODL_treat_pileup.bedGraph"),
                     genome = "hg38", type = "gradient", chromosome = viz.chr,
                     name = "bedGraph")
plot_list[[7]] <- dTrack4

dTrack5 <- DataTrack(range = paste0(dir, "/Subclass.SST_treat_pileup.bedGraph"),
                     genome = "hg38", type = "gradient", chromosome = viz.chr,
                     name = "bedGraph")
plot_list[[8]] <- dTrack5

dTrack6 <- DataTrack(range = paste0(dir, "/Subclass.PVALB_treat_pileup.bedGraph"),
                     genome = "hg38", type = "gradient", chromosome = viz.chr,
                     name = "bedGraph")
plot_list[[9]] <- dTrack6

dTrack7 <- DataTrack(range = paste0(dir, "/Subclass.L2-3_IT_treat_pileup.bedGraph"),
                     genome = "hg38", type = "gradient", chromosome = viz.chr,
                     name = "bedGraph")
plot_list[[10]] <- dTrack7

dTrack8 <- DataTrack(range = paste0(dir, "/Subclass.L5_IT_treat_pileup.bedGraph"),
                     genome = "hg38", type = "gradient", chromosome = viz.chr,
                     name = "bedGraph")
plot_list[[11]] <- dTrack8

dTrack9 <- DataTrack(range = paste0(dir, "/Subclass.L5_ET_treat_pileup.bedGraph"),
                     genome = "hg38", type = "gradient", chromosome = viz.chr,
                     name = "bedGraph")
plot_list[[12]] <- dTrack9

dTrack10 <- DataTrack(range = paste0(dir, "/Subclass.L5-6_NP_treat_pileup.bedGraph"),
                      genome = "hg38", type = "gradient", chromosome = viz.chr,
                      name = "bedGraph")
plot_list[[13]] <- dTrack10

dTrack11 <- DataTrack(range = paste0(dir, "/Subclass.L6_IT_treat_pileup.bedGraph"),
                      genome = "hg38", type = "gradient", chromosome = viz.chr,
                      name = "bedGraph")
plot_list[[14]] <- dTrack11

dTrack12 <- DataTrack(range = paste0(dir, "/Subclass.L6_IT_Car3_treat_pileup.bedGraph"),
                      genome = "hg38", type = "gradient", chromosome = viz.chr,
                      name = "bedGraph")
plot_list[[15]] <- dTrack12

dTrack13 <- DataTrack(range = paste0(dir, "/Subclass.L6_CT_treat_pileup.bedGraph"),
                      genome = "hg38", type = "gradient", chromosome = viz.chr,
                      name = "bedGraph")
plot_list[[16]] <- dTrack13

dTrack14 <- DataTrack(range = paste0(dir, "/Subclass.L6b_treat_pileup.bedGraph"),
                      genome = "hg38", type = "gradient", chromosome = viz.chr,
                      name = "bedGraph")
plot_list[[17]] <- dTrack14

dTrack15 <- DataTrack(range = paste0(dir, "/Subclass.OPC_treat_pileup.bedGraph"),
                      genome = "hg38", type = "gradient", chromosome = viz.chr,
                      name = "bedGraph")
plot_list[[18]] <- dTrack15

dTrack16 <- DataTrack(range = paste0(dir, "/Subclass.Astro_treat_pileup.bedGraph"),
                      genome = "hg38", type = "gradient", chromosome = viz.chr,
                      name = "bedGraph")
plot_list[[19]] <- dTrack16

dTrack17 <- DataTrack(range = paste0(dir, "/Subclass.Oligo_treat_pileup.bedGraph"),
                      genome = "hg38", type = "gradient", chromosome = viz.chr,
                      name = "bedGraph")
plot_list[[20]] <- dTrack17

dTrack18 <- DataTrack(range = paste0(dir, "/Subclass.Micro-PVM_treat_pileup.bedGraph"),
                      genome = "hg38", type = "gradient", chromosome = viz.chr,
                      name = "bedGraph")
plot_list[[21]] <- dTrack18

dTrack19 <- DataTrack(range = paste0(dir, "/Subclass.VLMC_treat_pileup.bedGraph"),
                      genome = "hg38", type = "gradient", chromosome = viz.chr,
                      name = "bedGraph")
plot_list[[22]] <- dTrack19

dTrack20 <- DataTrack(range = paste0(dir, "/Subclass.Endo_treat_pileup.bedGraph"),
                      genome = "hg38", type = "gradient", chromosome = viz.chr,
                      name = "bedGraph")
plot_list[[23]] <- dTrack20


## Replace gene annotation track with gene symbols
plot_list[[24]] <- GeneRegionTrack(txdb, chromosome = viz.chr, start = viz.min, end = viz.max,
                                   collapseTranscripts = "meta", geneSymbols = T,
                                   transcriptAnnotation = "symbol")
symbol(plot_list[[24]]) <- mapIds(org.Hs.eg.db, gene(plot_list[[24]]), 'SYMBOL', 'ENTREZID')

#add chromosomal location track
plot_list[[25]] <- IdeogramTrack(genome="hg38", chromosome=viz.chr)

## Make plots
Gviz::plotTracks(plot_list,
                 sizes = c(1.5, 0.1, 0.3, rep(0.25, 22)),
                 from = viz.min, to = viz.max, chromosome = viz.chr,
                 transcriptAnnotation = "symbol",
                 col.axis = "black",
                 fontsize.group = 10,
                 fontcolor.legend = "black",
                 lwd = 0.3,
                 title.width = 0.75,
                 background.title = "transparent",
                 col.border.title = "transparent",
                 cex.main = 4)

#Plot RNA expression
DefaultAssay(MOp.atac) <- "RNA"
Idents(object = MOp.atac) <- "subclass"
order <- c("LAMP5", "SNCG", "VIP", "SST CHODL", "SST", "PVALB", 
           "L2-3 IT", "L5 IT", "L5 ET", "L5-6 NP", "L6 IT", "L6 IT Car3", "L6 CT", "L6b",                                          
           "OPC", "Astro", "Oligo", "Micro-PVM", "VLMC", "Endo")  
Idents(object = MOp.atac) <- factor(Idents(object = MOp.atac), levels = order)
DotPlot(MOp.atac, features = "CUX2") + RotatedAxis()



#Feature Plots
#Make feature plots for CUX2
DefaultAssay(MOp.atac) <- "RNA"
FeaturePlot(
  object = MOp.atac,
  features = "CUX2",
  #min.cutoff = 'q10',
  #max.cutoff = 'q90',
  pt.size = 0.1, reduction = "r.umap",
  cols = c("lightgrey","#c43131")
) 

DefaultAssay(MOp.atac) <- "chromvar"
FeaturePlot(
  object = MOp.atac,
  features = "MA0755.1",
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1, reduction = "r.umap",
  cols = c("lightgrey","#c43131")
)


DefaultAssay(MOp.atac) <- "RNA"
FeaturePlot(
  object = MOp.atac,
  features = "EGR3",
  #min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1, reduction = "r.umap",
  cols = c("lightgrey","#c43131")
) 

DefaultAssay(MOp.atac) <- "chromvar"
FeaturePlot(
  object = MOp.atac,
  features = "MA0732.1",
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1, reduction = "r.umap",
  cols = c("lightgrey","#c43131")
)

DefaultAssay(MOp.atac) <- "ACTIVITY"
FeaturePlot(
  object = MOp.atac,
  features = c("CUX2"),
  #min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1, reduction = "r.umap",
  cols = c("lightgrey","#c43131") 
) 

DefaultAssay(MOp.atac) <- "RNA"
FeaturePlot(
  object = MOp.atac,
  features = "EGR3",
  #min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1, reduction = "r.umap",
  cols = c("lightgrey", "#012345")
) 

DefaultAssay(MOp.atac) <- "chromvar"
FeaturePlot(
  object = MOp.atac,
  features = "MA0732.1",
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1, reduction = "r.umap",
  cols = c("lightgrey", "#012345")
)

#Plot EGR3 Motif
MotifPlot(
  object = MOp.atac,
  motifs = "MA0732.1",
  assay = 'ATAC'
)



##


# Section 13: Marmoset SNARE-Seq2 RNA/AC - Generate Final Combined Seurat object --------------------------
library(Seurat)
library(Signac)
library(dplyr)
library(Matrix)
library(ggplot2)
library(cicero)
library(GenomeInfoDb)
library(TFBSTools)
library(BSgenome.Cjacchus.UCSC.calJac3)
set.seed(1234)


#Initial analyses to generate RNA count matrices, RNA Seurat objects,  
#RNA/AC coembedding, AC level consensus clusters and final combined 
#peak counts were performed according to human sections with any 
#exceptions outlined in the manuscript methods

#load RNA Seurat object (marMOp)
#load final AC peak by cell barcode count matrix 
counts

# keep only peaks from chr 1-22 or chrX or chrY
chrom <- sapply(rownames(counts), function(x) strsplit(x, split = ":")[[1]][[1]])
chrom <- sapply(chrom, function(x) strsplit(x, split = "_")[[1]][[1]])
loc <- sapply(rownames(counts), function(x) strsplit(x, split = ":")[[1]][[2]])
rownames(counts) <- paste0(chrom, ":", loc)
chrom.selected <- as.factor(chrom)
chrom.remove <- names(chrom.selected[chrom.selected %in% c("chrUn", "chrM")])
chrom.keep <- setdiff(names(chrom.selected), chrom.remove)
filtered.counts<- counts[which(rownames(counts) %in% chrom.keep),]

counts <- filtered.counts

###create and pre-process atac seurat (marMOp.atac)
marMOp.atac <- CreateSeuratObject(counts = counts, assay = "ATAC", project = "MOp_ATAC")

meta #metadata table
meta <- meta[colnames(marMOp.atac), ]
marMOp.atac <- AddMetaData(marMOp.atac, metadata = meta)

#pre-processing
DefaultAssay(marMOp.atac) <- "ATAC"
VariableFeatures(marMOp.atac) <- names(which(Matrix::rowSums(marMOp.atac) > 20))
marMOp.atac <- RunLSI(marMOp.atac, n = 50, scale.max = NULL) #latent semantic indexing
marMOp.atac <- RunUMAP(marMOp.atac, reduction = "lsi", dims = 1:40,
                       n.neighbors = 50L, min.dist = 0.5)


###Add in RNA data
#load RNA seurat object = marMOp
DefaultAssay(marMOp.atac) <- "ATAC"

cells.use = rownames(marMOp.atac)
cells.use <- gsub("MARMOP","marMOp", cells.use)
cells.use <- cells.use[cells.use %in% colnames(marMOp)]
marMOp <- subset(marMOp, cells = cells.use)
rna.counts <- GetAssayData(object = marMOp, slot = "counts")
colnames(rna.counts) <- gsub("marMOp","MARMOP", colnames(rna.counts))
marMOp.atac[["RNA"]] <- CreateAssayObject(counts = rna.counts)
DefaultAssay(marMOp.atac) <- "RNA"
umap.coordinates <- Embeddings(object = marMOp, reduction = "umap")
rownames(umap.coordinates) <- gsub("marMOp","MARMOP", rownames(umap.coordinates))
marMOp.atac[["r.umap"]] <- CreateDimReducObject(embeddings = umap.coordinates, key = "r_umap_", assay = DefaultAssay(marMOp.atac))
norm.data <- GetAssayData(object = marMOp, slot = "data")
colnames(norm.data) <- gsub("marMOp","MARMOP", colnames(norm.data))
marMOp.atac <- SetAssayData(object = marMOp.atac, slot = "data", new.data = norm.data)
scale.data <- GetAssayData(object = marMOp, slot = "scale.data")
colnames(scale.data) <- gsub("marMOp","MARMOP", colnames(scale.data))
marMOp.atac <- SetAssayData(object = marMOp.atac, slot = "scale.data", new.data = scale.data)

DefaultAssay(marMOp.atac) <- "ATAC"


###Add in motif matrix
file = "JASPAR2020_CORE_vertebrates_non-redundant_pfms_jaspar.txt"
pfm <- readJASPARMatrix(file, matrixClass="PFM")

# Scan the DNA sequence of each peak for the presence of each motif
motif.matrix <- CreateMotifMatrix(
  features = StringToGRanges(rownames(marMOp.atac), sep = c(":", "-")),
  pwm = pfm,
  genome = 'BSgenome.Cjacchus.UCSC.calJac3',
  sep = c(":", "-")
)

# Create a new Mofif object to store the results
motif <- CreateMotifObject(
  data = motif.matrix,
  pwm = pfm
)

# Add the Motif object to the assay
marMOp.atac[['ATAC']] <- AddMotifObject(
  object = marMOp.atac[['ATAC']],
  motif.object = motif
)

marMOp.atac <- RegionStats(
  object = marMOp.atac,
  genome = BSgenome.Cjacchus.UCSC.calJac3,
  sep = c(":", "-")
)

###Add in Chromvar TF activites
library(chromVAR)
library(motifmatchr)
library(GenomicRanges)
library(SummarizedExperiment)
library(chromVARmotifs)
library(TFBSTools)
source("/media/Scratch_SSD_Voyager/Blue/chrom.R")

counts <- GetAssayData(object = marMOp.atac, slot = "counts")
counts@x[counts@x > 0] <- 1 #binarize counts

## Function for splitting peak names into bed file dataframes
peak2df <- function(peak.names, keep.colnames = F, metadata.df = NULL) {
  chrom <- sapply(peak.names, function(x) strsplit(x, split = ":")[[1]][[1]])
  bp.range <- sapply(peak.names, function(x) strsplit(x, split = ":")[[1]][[2]])
  bp1 <- as.integer(sapply(bp.range, function(x) strsplit(x, split = "-")[[1]][[1]]))
  bp2 <- as.integer(sapply(bp.range, function(x) strsplit(x, split = "-")[[1]][[2]]))
  
  bed.df <- data.frame("chr" = chrom, "start" = bp1, "end" = bp2)
  if (!is.null(metadata.df)) bed.df <- cbind(bed.df, metadata.df)
  
  if (!keep.colnames) {
    colnames(bed.df) <- NULL
  }
  bed.df
}

# Make peaks GRanges
peaks.df <- peak2df(rownames(counts)); colnames(peaks.df) <- c("chr", "start", "end");
peaks <- makeGRangesFromDataFrame(peaks.df, ignore.strand = T)

# Make SummarizedExperiment and add GC Bias
fragment_counts <- SummarizedExperiment(assays = list(counts = counts), rowRanges = peaks)
fragment_counts <- addGCBias(fragment_counts, genome = BSgenome.Cjacchus.UCSC.calJac3)
dim(fragment_counts)

# addGCBias introduces NAs in the data, remove those rows from filtered.counts
# and recreate the objects
filtered.counts <- counts[!is.na(fragment_counts@rowRanges@elementMetadata@listData$bias), ]
filtered.counts <- filtered.counts[rowSums(filtered.counts) > 10,]
peaks.df <- peak2df(rownames(filtered.counts))
colnames(peaks.df) <- c("chr", "start", "end")
peaks.gr <- makeGRangesFromDataFrame(peaks.df, ignore.strand = T)
fragment_counts <- SummarizedExperiment(assays = list(counts = filtered.counts), rowRanges = peaks.gr)
fragment_counts <- addGCBias(fragment_counts, genome = BSgenome.Cjacchus.UCSC.calJac3)

## Get motifs
motif_ix <- matchMotifs(pfm, fragment_counts, genome = BSgenome.Cjacchus.UCSC.calJac3)

## Computing deviations
BiocParallel::register(BiocParallel::SerialParam())
dev <- computeDeviations(object = fragment_counts, annotations = motif_ix)

## Compute variability
variability <- computeVariability(dev)
variability <- variability[order(variability$p_value_adj),]
dev.mat <- as(dev@assays$data$z, "dgCMatrix")

marMOp.atac[["chromvar"]] <- CreateAssayObject(counts = dev.mat)


###Cicero analysis: generate predicted gene activity matrix from chrom data
Idents(object = marMOp.atac) <- "AC_cluster"
Idents(object = marMOp.atac) <- factor(Idents(object = marMOp.atac), levels = 1:24)

## Run cicero
calJac3.chr.lengths <- read.table("calJac3.chrom.sizes", header = F, sep = "\t") ## A tab separated text file with chromosome lengths
calJac3.chr.lengths[[2]] <- as.numeric(calJac3.chr.lengths[[2]])

clusters <- Idents(object = marMOp.atac)
count.matrix <- GetAssayData(object = marMOp.atac, slot = "counts")
count.matrix@x[count.matrix@x > 0] <- 1 #binarize counts
umap.coordinates <- Embeddings(object = marMOp.atac, reduction = "umap")

#update counts matrix
rownames <- rownames(marMOp.atac)
chrom <- sapply(rownames, function(x) strsplit(x, split = ":")[[1]][[1]])
loc <- sapply(rownames, function(x) strsplit(x, split = ":")[[1]][[2]])
loc_start <- sapply(loc, function(x) strsplit(x, split = "-")[[1]][[1]])
loc_end <- sapply(loc, function(x) strsplit(x, split = "-")[[1]][[2]])

pData <- data.frame(clusters)
fData <- data.frame(site_name = rownames(count.matrix), chromosome = chrom, bp1 = loc_start, bp2 = loc_end)
input_cds <- newCellDataSet(count.matrix, phenoData = new("AnnotatedDataFrame", data = pData),
                            featureData = new("AnnotatedDataFrame", data = fData),
                            expressionFamily = VGAM::binomialff())

cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = umap.coordinates, k = 50)
conns <- run_cicero(cicero_cds, calJac3.chr.lengths) # Takes a few minutes to run

## Assign peaks to modules
ccan.assigns <- generate_ccans(conns)
#"Coaccessibility cutoff used: 0.07"

#generate annotation file
gtf <- rtracklayer::import('refdata-cellranger-marmoset-3.2-GCF_000004665.1/Callithrix_jacchus.GCF_000004665.1/genes/genes.gtf')
ch.conv <- read.delim("refdata-cellranger-marmoset-3.2-GCF_000004665.1/calJac3.2_ucscToRefSeq.txt.gz")
rownames(ch.conv) <- ch.conv$name
ch.conv <- ch.conv[seqlevels(gtf),]
seqlevels(gtf) <- as.character(ch.conv$X.chrom)
gene.coords <- gtf[gtf$type == 'gene']
seqlevelsStyle(gene.coords) <- 'UCSC'
gene.coords <- keepStandardChromosomes(gene.coords, pruning.mode = 'coarse')
genebodyandpromoter.coords <- Extend(x = gene.coords, upstream = 10000, downstream = 10000)

peaks.file <- makeGRangesFromDataFrame(fData,
                                       keep.extra.columns=TRUE,
                                       ignore.strand=TRUE,
                                       seqinfo=NULL,
                                       seqnames.field= "chromosome",
                                       start.field="bp1",
                                       end.field="bp2",
                                       strand.field="strand",
)
peak_anno_df <- as.data.frame(peaks.file)

peak_anno <- ClosestFeature(
  regions = rownames(peak_anno_df),
  annotation = genebodyandpromoter.coords,
  sep = c(':', '-')
)

rownames(peak_anno) <- paste0(peak_anno_df$seqnames, ":", peak_anno_df$start, "-", peak_anno_df$end)
peak_anno_df$gene <- peak_anno$gene
peak_anno_df <- peak_anno_df[,c("seqnames", "start", "end", "gene")]
colnames(peak_anno_df) <- c("chromosome", "start", "end", "gene")
peak_anno_df <- peak_anno_df[intersect(rownames(input_cds), rownames(peak_anno_df)),]

## Annotate sites by gene
input_cds <- input_cds[intersect(rownames(input_cds), rownames(peak_anno_df)),]
input_cds <- annotate_cds_by_site(input_cds, peak_anno_df)

## Generate unnormalized gene activity matrix
unnorm_ga <- build_gene_activity_matrix(input_cds, conns, coaccess_cutoff = 0.25)

marMOp.atac[["ACTIVITY"]] <- CreateAssayObject(counts = unnorm_ga[,names(clusters)])
DefaultAssay(marMOp.atac) <- "ACTIVITY"
marMOp.atac <- NormalizeData(object = marMOp.atac)
marMOp.atac <- ScaleData(marMOp.atac)



# Section 14: Marmoset SNARE-Seq2 RNA/AC - Plots and Tables --------------------------------------------------
library(Signac)
library(Seurat)
library(dplyr)
library(Matrix)
library(methods)
library(swne)
library(cicero)
library(ggplot2)
set.seed(1234)

#load marMOp.atac

DefaultAssay(marMOp.atac) <- "RNA"
rna.counts <- GetAssayData(object = marMOp.atac, slot = "counts")
saveRDS(rna.counts, file="Zhang_BICCN-H_20190730_20190903_marMOp_Final_RNA_Counts.RDS")

DefaultAssay(marMOp.atac) <- "ATAC"
atac.counts <- GetAssayData(object = marMOp.atac, slot = "counts")
saveRDS(atac.counts, file="Zhang_BICCN-H_20190730_20190903_marMOp_Final_AC_Peaks.RDS")

write.table(marMOp.atac@meta.data, file="Zhang_BICCN-H_20190730_20190903_marMOp_Final_Sample_Metadata.txt", sep = "\t", row.names=TRUE, col.names=TRUE)


#UMAP Plots

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


#AC-level Clusters
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


#Patient
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










##Plot QC statistics by cluster

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




# Section 15: Marmoset SNARE-Seq2 RNA/AC - Motif CCAN analyses - Subclasses ------------------------------------
library(Signac)
library(Seurat)
library(dplyr)
library(Matrix)
library(methods)
library(swne)
library(cicero)
library(ggplot2)
library(viridis)

#load marMOp.atac Seurat object

##Identify and plot RNA marker genes by subclass
#RNA Marker genes for Subclass-level clusters
DefaultAssay(marMOp.atac) <- "RNA"
Idents(object = marMOp.atac) <- "subclass"
order <- c("Lamp5", "Sncg", "Meis2", "Vip", "Sst Chodl", "Sst", "Pvalb", 
           "L2-3 IT", "L5 IT", "L5 ET", "L5-6 NP", "L6 IT", "L6 IT Car3", "L6 CT", "L6b",                                          
           "OPC", "Astro", "Oligo", "Micro-PVM", "VLMC","Peri", "Endo")

Idents(object = marMOp.atac) <- factor(Idents(object = marMOp.atac), levels = order)
table(Idents(object = marMOp.atac))
#remove clusters that are too small
marMOp.atac <- subset(marMOp.atac, idents = c("Meis2", "Sst Chodl"), invert = TRUE)

genes.use <- rownames(marMOp.atac)
SC.level.markers <- FindAllMarkers(marMOp.atac, features = genes.use,
                                   only.pos = TRUE, logfc.threshold = 0.25)
RNA.mark <- SC.level.markers
RNA.mark <- RNA.mark[RNA.mark$p_val_adj < 0.05, ]
RNA.mark <- RNA.mark[RNA.mark$avg_logFC > 1, ]
RNA.mark <- distinct(RNA.mark, gene, .keep_all = TRUE)
RNA.mark %>% group_by(cluster) %>% top_n(10, avg_logFC) -> top10 

ave.RNA <- AverageExpression(marMOp.atac, assays = "RNA", features = top10$gene, slot = "data" )

scaled <- t(scale(t(ave.RNA$RNA)))
scaled <- scale(scaled)
range(scaled)
scaled[scaled < 0] <- 0
scaled[scaled > 4] <- 4

ggHeat((scaled), rescaling = "none", clustering = "none", x.lab.size = 8, y.lab.size = 11,
) + scale_fill_viridis()

#plot associated activites
DefaultAssay(marMOp.atac) <- "ACTIVITY"
gene.use <- top10$gene[top10$gene %in% rownames(marMOp.atac)]
ave.act <- AverageExpression(marMOp.atac, assays = "ACTIVITY", features = gene.use, slot = "data" )

scaled <- t(scale(t(ave.act$ACTIVITY)))
scaled <- scale(scaled)
range(scaled)
scaled[scaled < 0] <- 0
scaled[scaled > 4] <- 4

ggHeat((scaled), rescaling = "none", clustering = "none", x.lab.size = 8, y.lab.size = 11,
) + scale_fill_viridis()


#TF motif analysis in subclass DEGs
DefaultAssay(marMOp.atac) <- "chromvar"
Idents(object = marMOp.atac) <- "subclass"
order <- c("Lamp5", "Sncg", "Vip", "Sst", "Pvalb", 
           "L2-3 IT", "L5 IT", "L5 ET", "L5-6 NP", "L6 IT", "L6 IT Car3", "L6 CT", "L6b",                                          
           "OPC", "Astro", "Oligo", "Micro-PVM", "VLMC","Peri", "Endo")
Idents(object = marMOp.atac) <- factor(Idents(object = marMOp.atac), levels = order)

tf.markers <- FindAllMarkers(
  object = marMOp.atac,
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'nCount_ATAC'
)





##Marker gene TFBS enrichments
#load conns, ccan.assigns and peak_anno_df objects from Cicero analyses
#load SC.level.markers (subclass level RNA markers)

peak_anno_df <- peak_anno_df[rownames(ccan.assigns),]
peak_anno_df$CCAN <- ccan.assigns$CCAN
head(peak_anno_df)

#Associate CCANs with subclass Markers
DefaultAssay(marMOp.atac) <- "ATAC"
Idents(object = marMOp.atac) <- "subclass"
order <- c("Lamp5", "Sncg", "Meis2", "Vip", "Sst Chodl", "Sst", "Pvalb", 
           "L2-3 IT", "L5 IT", "L5 ET", "L5-6 NP", "L6 IT", "L6 IT Car3", "L6 CT", "L6b",                                          
           "OPC", "Astro", "Oligo", "Micro-PVM", "VLMC","Peri", "Endo")

Idents(object = marMOp.atac) <- factor(Idents(object = marMOp.atac), levels = order)
table(Idents(object = marMOp.atac))
#remove clusters that are too small
marMOp.atac <- subset(marMOp.atac, idents = c("Meis2", "Sst Chodl"), invert = TRUE)

RNA.mark <- SC.level.markers
RNA.mark <- RNA.mark[RNA.mark$p_val_adj < 0.05, ]
RNA.mark <- RNA.mark[RNA.mark$avg_logFC > 1, ]
RNA.mark <- distinct(RNA.mark, gene, .keep_all = TRUE)

peak_anno_df <- peak_anno_df[peak_anno_df$gene %in% RNA.mark$gene,]
peak_anno_df_sc <- left_join(peak_anno_df, RNA.mark, by = "gene")
rownames(peak_anno_df_sc) <- rownames(peak_anno_df)
peak_anno_df_sc <-na.omit(peak_anno_df_sc)

#Re-order subclasses
peak_anno_df_sc$site <- rownames(peak_anno_df_sc)
peak_anno_df_sc$cluster_order <- factor(peak_anno_df_sc$cluster)
peak_anno_df_sc$cluster_order <- plyr::mapvalues(peak_anno_df_sc$cluster_order,
                                 from = c("Lamp5", "Sncg", "Vip", "Sst", "Pvalb", 
                                 "L2-3 IT", "L5 IT", "L5 ET", "L5-6 NP", "L6 IT Car3", "L6 CT", "L6b",                                          
                                 "OPC", "Astro", "Oligo", "Micro-PVM", "VLMC","Peri", "Endo"), 
                                 to = c(1:19))
peak.ccam.mark <-arrange(peak_anno_df_sc, cluster_order, desc(avg_logFC))
rownames(peak.ccam.mark) <- peak.ccam.mark$site
peak.ccam.mark


DefaultAssay(marMOp.atac) <- "ATAC"

#find gene motif enrichments (genes having >1 site)
tab <- table(peak.ccam.mark$gene)
unique.genes <- names(tab[tab > 1])

gene.tfbs.list <- lapply(unique.genes, function(cl) {
  print(paste("Running for gene:", cl))
  
  cl.sites <- rownames(peak_anno_df[peak_anno_df$gene %in% cl,])
  EM <- FindMotifs(object = marMOp.atac,
                   features = cl.sites)
  EM <- EM[EM$pvalue < 0.05,]
  EM
})
names(gene.tfbs.list) <- unique.genes
gene.tfbs <- do.call("rbind", lapply(gene.tfbs.list, as.data.frame)) 

#Add associated cluster (predicted by DEG analysis) for each marker
marker<-unlist(lapply(rownames(gene.tfbs),function(x) unlist(strsplit(x,"[.]"))[1]))
gene.tfbs$gene <- marker
cl.order <- peak.ccam.mark$cluster
names(cl.order) <- peak.ccam.mark$gene
cl.order <- factor(cl.order) 
cl.order <- cl.order[gene.tfbs$gene]
gene.tfbs$cluster <- cl.order


#find active TFBS 
DefaultAssay(marMOp.atac) <- "chromvar"
tf.markers <- FindAllMarkers(
  object = marMOp.atac,
  features = unique(gene.tfbs$motif),
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'nCount_ATAC'
)
colnames(tf.markers)[colnames(tf.markers) == "gene"] <- "motif"
tf.markers <- left_join(gene.tfbs, tf.markers, by = c("motif","cluster"))
rownames(tf.markers) <- rownames(gene.tfbs)
mar.tf.markers <-na.omit(tf.markers)
mar.tf.markers


#re-order by cluster
mar.tf.markers$rownames <- rownames(mar.tf.markers)
mar.tf.markers$cluster_order <- mar.tf.markers$cluster
mar.tf.markers$cluster_order <- plyr::mapvalues(mar.tf.markers$cluster_order,
                                from = c("Lamp5", "Sncg", "Vip", "Sst", "Pvalb", 
                                "L2-3 IT", "L5 IT", "L5 ET", "L5-6 NP", "L6 IT Car3", "L6 CT", "L6b",                                          
                                "OPC", "Astro", "Oligo", "Micro-PVM", "VLMC","Peri", "Endo"), 
                                to = c(1:19))
mar.tf.markers$cluster_order <- as.numeric(mar.tf.markers$cluster_order)
mar.tf.mark <-arrange(mar.tf.markers, cluster_order, desc(avg_logFC))
rownames(mar.tf.mark) <- mar.tf.mark$rownames
mar.tf.mark


tf.mark <- mar.tf.mark[mar.tf.mark$p_val_adj < 0.05,]
tf.mark <- tf.mark[tf.mark$avg_logFC > 0.5,]
tf.mark <- distinct(tf.mark, motif, .keep_all = TRUE)

tf.mark %>% group_by(cluster) %>% top_n(5, avg_logFC) -> top5 
tf.mark %>% group_by(cluster) %>% top_n(10, avg_logFC) -> top10 

sites.use <- top10$motif
ave.sites <- AverageExpression(marMOp.atac, assays = "chromvar", features = sites.use, slot = "data" )

names.sites <- top10$motif.name
ave.sites <- ave.sites$chromvar
rownames(ave.sites) <- names.sites

scaled <- t(scale(t(ave.sites)))
scaled <- scale(scaled)
range(scaled)
scaled[scaled < 0] <- 0
scaled[scaled > 4] <- 4


#Plot marker-enriched TFBSs that are active in the subclasses
ggHeat(scaled, rescaling = "none", clustering = "none", x.lab.size = 8, y.lab.size = 11,
) + scale_fill_viridis()



# Section 16: Marmoset SNARE-Seq2 RNA/AC - ChC versus Basket Cells - Motif CCAN analyses -----------------------------------------------------------
library(Signac)
library(Seurat)
library(dplyr)
library(Matrix)
library(methods)
library(cicero)
library(swne)
library(ggplot2)
library(viridis)

#load marMOp.atac

###Run Cicero only on PVALB neurons
Idents(object = marMOp.atac) <- "AC_cluster"
Idents(object = marMOp.atac) <- factor(Idents(object = marMOp.atac), levels = 1:24)

marMOp.atac.pvalb <- subset(marMOp.atac, idents = 5:6)
DefaultAssay(marMOp.atac.pvalb) <- "ATAC"

#Cicero analysis: generate predicted gene activity matrix from chrom data
calJac3.chr.lengths <- read.table("calJac3.chrom.sizes", header = F, sep = "\t") ## A tab separated text file with chromosome lengths
calJac3.chr.lengths[[2]] <- as.numeric(calJac3.chr.lengths[[2]])

clusters <- Idents(object = marMOp.atac.pvalb)
count.matrix <- GetAssayData(object = marMOp.atac.pvalb, slot = "counts")
count.matrix@x[count.matrix@x > 0] <- 1 #binarize counts
umap.coordinates <- Embeddings(object = marMOp.atac.pvalb, reduction = "umap")

#update counts matrix
rownames <- rownames(marMOp.atac.pvalb)
chrom <- sapply(rownames, function(x) strsplit(x, split = ":")[[1]][[1]])
loc <- sapply(rownames, function(x) strsplit(x, split = ":")[[1]][[2]])
loc_start <- sapply(loc, function(x) strsplit(x, split = "-")[[1]][[1]])
loc_end <- sapply(loc, function(x) strsplit(x, split = "-")[[1]][[2]])

pData <- data.frame(clusters)
fData <- data.frame(site_name = rownames(count.matrix), chromosome = chrom, bp1 = loc_start, bp2 = loc_end)
input_cds <- newCellDataSet(count.matrix, phenoData = new("AnnotatedDataFrame", data = pData),
                            featureData = new("AnnotatedDataFrame", data = fData),
                            expressionFamily = VGAM::binomialff())

cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = umap.coordinates, k = 50)
conns <- run_cicero(cicero_cds, calJac3.chr.lengths) # Takes a few minutes to run

# Assign peaks to modules
ccan.assigns <- generate_ccans(conns, coaccess_cutoff_override = 0.07)
# "Coaccessibility cutoff used: 0.19"

#generate annotation file
gtf <- rtracklayer::import('refdata-cellranger-marmoset-3.2-GCF_000004665.1/Callithrix_jacchus.GCF_000004665.1/genes/genes.gtf')
ch.conv <- read.delim("refdata-cellranger-marmoset-3.2-GCF_000004665.1/calJac3.2_ucscToRefSeq.txt.gz")
rownames(ch.conv) <- ch.conv$name
ch.conv <- ch.conv[seqlevels(gtf),]
seqlevels(gtf) <- as.character(ch.conv$X.chrom)
gene.coords <- gtf[gtf$type == 'gene']
seqlevelsStyle(gene.coords) <- 'UCSC'
gene.coords <- keepStandardChromosomes(gene.coords, pruning.mode = 'coarse')
genebodyandpromoter.coords <- Extend(x = gene.coords, upstream = 10000, downstream = 10000)

peaks.file <- makeGRangesFromDataFrame(fData,
                                       keep.extra.columns=TRUE,
                                       ignore.strand=TRUE,
                                       seqinfo=NULL,
                                       seqnames.field= "chromosome",
                                       start.field="bp1",
                                       end.field="bp2",
                                       strand.field="strand",
)
peak_anno_df <- as.data.frame(peaks.file)


peak_anno <- ClosestFeature(
  regions = rownames(peak_anno_df),
  annotation = genebodyandpromoter.coords,
  sep = c(':', '-')
)

rownames(peak_anno) <- paste0(peak_anno_df$seqnames, ":", peak_anno_df$start, "-", peak_anno_df$end)
peak_anno_df$gene <- peak_anno$gene
peak_anno_df <- peak_anno_df[,c("seqnames", "start", "end", "gene")]
colnames(peak_anno_df) <- c("chromosome", "start", "end", "gene")

peak_anno_df <- peak_anno_df[rownames(ccan.assigns),]
peak_anno_df$CCAN <- ccan.assigns$CCAN
head(peak_anno_df)

#re-make activity matrix
## Annotate sites by gene
input_cds <- input_cds[intersect(rownames(input_cds), rownames(peak_anno_df)),]
input_cds <- annotate_cds_by_site(input_cds, peak_anno_df)

## Generate unnormalized gene activity matrix
unnorm_ga <- build_gene_activity_matrix(input_cds, conns, coaccess_cutoff = 0.25)
marMOp.atac.pvalb[["ACTIVITY"]] <- CreateAssayObject(counts = unnorm_ga[,names(clusters)])
DefaultAssay(marMOp.atac.pvalb) <- "ACTIVITY"
marMOp.atac.pvalb <- NormalizeData(object = marMOp.atac.pvalb)
marMOp.atac.pvalb <- ScaleData(marMOp.atac.pvalb)

#Associate CCANs with ChC Markers
DefaultAssay(marMOp.atac.pvalb) <- "RNA"

#Find all ChC markers
genes.use <- rownames(marMOp.atac.pvalb)
ChC.markers <- FindMarkers(marMOp.atac.pvalb, 
                           ident.1 = 6,
                           ident.2 = 5,
                           features = genes.use,
                           only.pos = TRUE, logfc.threshold = 0.25)
RNA.mark <- ChC.markers
RNA.mark <- RNA.mark[RNA.mark$p_val_adj < 0.05, ]
RNA.mark <- RNA.mark[RNA.mark$avg_logFC > 0.5, ]

ave.RNA <- AverageExpression(marMOp.atac, assays = "RNA", features = rownames(RNA.mark), slot = "data" )

ggHeat((ave.RNA$RNA[1:6]), rescaling = "row", clustering = "row", x.lab.size = 8, y.lab.size = 11,
) + scale_fill_viridis()

#Co-accessible peaks associated with marker genes
peak_anno_df <- peak_anno_df[peak_anno_df$gene %in% rownames(RNA.mark),]
RNA.mark$gene <- rownames(RNA.mark)
peak_anno_df_pv <- left_join(peak_anno_df, RNA.mark, by = "gene")
rownames(peak_anno_df_pv) <- rownames(peak_anno_df)
peak_anno_df_pv <-na.omit(peak_anno_df_pv)
peak_anno_df_pv$site <- rownames(peak_anno_df_pv)
write.table(peak_anno_df_pv, file="Zhang_marMOp_ChC_CCAN_peak_gene_annotation_cluster_DEGs.txt", sep = "\t", row.names=TRUE, col.names=TRUE)

length(unique(peak_anno_df_pv$gene))
#122

marMOp.atac.chand <- subset(marMOp.atac.pvalb, idents = 6)
DefaultAssay(marMOp.atac.chand) <- "ATAC"

#find gene motif enrichments (genes having >1 site)
tab <- table(peak_anno_df_pv$gene)
unique.genes <- names(tab[tab > 1])

gene.tfbs.list <- lapply(unique.genes, function(cl) {
  print(paste("Running for gene:", cl))
  
  cl.sites <- rownames(peak_anno_df[peak_anno_df$gene %in% cl,])
  EM <- FindMotifs(object = marMOp.atac.chand,
                   features = cl.sites)
  EM <- EM[EM$pvalue < 0.05,]
  EM
})
names(gene.tfbs.list) <- unique.genes
gene.tfbs <- do.call("rbind", lapply(gene.tfbs.list, as.data.frame)) 
marker<-unlist(lapply(rownames(gene.tfbs),function(x) unlist(strsplit(x,"[.]"))[1]))
gene.tfbs$gene <- marker

#find active TFBS in ChC cells over Basket Cells
DefaultAssay(marMOp.atac.pvalb) <- "chromvar"
tf.markers <- FindMarkers(
  object = marMOp.atac.pvalb, 
  ident.1 = 6,
  #ident.2 = 5,
  features = unique(gene.tfbs$motif),
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'nCount_ATAC'
)
tf.markers$motif <- rownames(tf.markers)
tf.markers <- left_join(gene.tfbs, tf.markers, by = c("motif"))
rownames(tf.markers) <- rownames(gene.tfbs)
mar.tf.markers <-na.omit(tf.markers)
mar.tf.markers

tf.mark <- mar.tf.markers[mar.tf.markers$p_val < 0.05,]
tf.mark <- distinct(tf.mark, motif, .keep_all = TRUE)

sites.use <- tf.mark$motif
ave.sites <- AverageExpression(marMOp.atac, assays = "chromvar", features = sites.use, slot = "data" )

names.sites <- tf.mark$motif.name
ave.sites <- ave.sites$chromvar
rownames(ave.sites) <- names.sites

#Plot TFBSs that are enriched in ChC markers and active in the ChC cluster
ggHeat(ave.sites[1:6], rescaling = "row", clustering = "row", x.lab.size = 8, y.lab.size = 11,
) + scale_fill_viridis()




# Section 17: Mouse scATAC-seq - ChC versus Basket Cells - Motif CCAN analyses ------------------------------------
library(Signac)
library(Seurat)
library(dplyr)
library(Matrix)
library(methods)
library(swne)
library(cicero)
library(ggplot2)

#load mouse peak by cell barcode count matrix
counts

msMOp.atac <- CreateSeuratObject(counts = counts, assay = "ATAC", project = "MOp_ATAC")
labels #table of cell cluster annotations
msMOp.atac <- AddMetaData(msMOp.atac, metadata = labels)

#pre-processing
DefaultAssay(msMOp.atac) <- "ATAC"
VariableFeatures(msMOp.atac) <- names(which(Matrix::rowSums(msMOp.atac) > 100))
msMOp.atac <- RunLSI(msMOp.atac, n = 50, scale.max = NULL) #latent semantic indexing
msMOp.atac <- RunUMAP(msMOp.atac, reduction = "lsi", dims = 1:50)

###Run Cicero only on PVALB neurons
Idents(object = msMOp.atac) <- "labels"
msMOp.atac <- subset(msMOp.atac, idents = c("Pv_Ntf3_Trim63","Pv_Tac1","Pvalb_Vipr2"))

#Cicero analysis: generate predicted gene activity matrix from chrom data
"wget http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes"
mm10.chr.lengths <- read.table("mm10.chrom.sizes", header = F, sep = "\t") ## A tab separated text file with chromosome lengths
mm10.chr.lengths[[2]] <- as.numeric(mm10.chr.lengths[[2]])

clusters <- Idents(object = msMOp.atac)
count.matrix <- GetAssayData(object = msMOp.atac, slot = "counts")
count.matrix@x[count.matrix@x > 0] <- 1 #binarize counts
umap.coordinates <- Embeddings(object = msMOp.atac, reduction = "umap")

#update counts matrix
rownames <- rownames(msMOp.atac)
chrom <- sapply(rownames, function(x) strsplit(x, split = ":")[[1]][[1]])
loc <- sapply(rownames, function(x) strsplit(x, split = ":")[[1]][[2]])
loc_start <- sapply(loc, function(x) strsplit(x, split = "-")[[1]][[1]])
loc_end <- sapply(loc, function(x) strsplit(x, split = "-")[[1]][[2]])

pData <- data.frame(clusters)
fData <- data.frame(site_name = rownames(count.matrix), chromosome = chrom, bp1 = loc_start, bp2 = loc_end)
input_cds <- newCellDataSet(count.matrix, phenoData = new("AnnotatedDataFrame", data = pData),
                            featureData = new("AnnotatedDataFrame", data = fData),
                            expressionFamily = VGAM::binomialff())

cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = umap.coordinates, k = 50)
conns <- run_cicero(cicero_cds, mm10.chr.lengths) # Takes a few minutes to run

## Assign peaks to modules
ccan.assigns <- generate_ccans(conns, coaccess_cutoff = 0.08)
# "Coaccessibility cutoff used: 0.08"

library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
peaks.file <- makeGRangesFromDataFrame(ccan.data,
                                       keep.extra.columns=TRUE,
                                       ignore.strand=TRUE,
                                       seqinfo=NULL,
                                       seqnames.field= "chromosome",
                                       start.field="bp1",
                                       end.field="bp2",
                                       strand.field="strand",
)
peak_anno <- annotatePeak(peaks.file, tssRegion = c(-10000, 10000),
                          TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene,annoDb ="org.Mm.eg.db")

peak_anno_df <- as.data.frame(peak_anno)
rownames(peak_anno_df) <- paste0(peak_anno_df$seqnames, ":", peak_anno_df$start, "-", peak_anno_df$end)
peak_anno_df <- peak_anno_df[,c("seqnames", "start", "end", "SYMBOL")]
colnames(peak_anno_df) <- c("chromosome", "start", "end", "gene")

peak_anno_df <- peak_anno_df[rownames(ccan.assigns),]
peak_anno_df$CCAN <- ccan.assigns$CCAN

###add activity matrix
## Format peak annotation dataframe
peak_anno_df <- as.data.frame(peak_anno)
rownames(peak_anno_df) <- paste0(peak_anno_df$seqnames, ":", peak_anno_df$start, "-", peak_anno_df$end)
peak_anno_df <- peak_anno_df[,c("seqnames", "start", "end", "SYMBOL")]
colnames(peak_anno_df) <- c("chromosome", "start", "end", "gene")

peak_anno_df <- peak_anno_df[intersect(rownames(input_cds), rownames(peak_anno_df)),]

## Annotate sites by gene
input_cds <- input_cds[intersect(rownames(input_cds), rownames(peak_anno_df)),]
input_cds <- annotate_cds_by_site(input_cds, peak_anno_df)

## Generate unnormalized gene activity matrix
unnorm_ga <- build_gene_activity_matrix(input_cds, conns, coaccess_cutoff = 0.25)

msMOp.atac[["ACTIVITY"]] <- CreateAssayObject(counts = unnorm_ga[,names(clusters)])
DefaultAssay(msMOp.atac) <- "ACTIVITY"
msMOp.atac <- NormalizeData(object = msMOp.atac)
msMOp.atac <- ScaleData(msMOp.atac)


#Add in motif matrix
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
DefaultAssay(msMOp.atac) <- "ATAC"

file = "/media/Scratch_SSD_Voyager/Blue/JASPAR2020_CORE_vertebrates_non-redundant_pfms_jaspar.txt"
pfm <- readJASPARMatrix(file, matrixClass="PFM")

# Scan the DNA sequence of each peak for the presence of each motif
motif.matrix <- CreateMotifMatrix(
  features = StringToGRanges(rownames(msMOp.atac), sep = c(":", "-")),
  pwm = pfm,
  genome = 'mm10',
  sep = c(":", "-")
)

# Create a new Mofif object to store the results
motif <- CreateMotifObject(
  data = motif.matrix,
  pwm = pfm
)

# Add the Motif object to the assay
msMOp.atac[['ATAC']] <- AddMotifObject(
  object = msMOp.atac[['ATAC']],
  motif.object = motif
)

msMOp.atac <- RegionStats(
  object = msMOp.atac,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  sep = c(":", "-")
)


#Add in Chromvar TF activites
msMOp.atac <- RunChromVAR(
  object = msMOp.atac,
  genome = BSgenome.Mmusculus.UCSC.mm10
)



###Find TFBS enrichments within marker gene-associated CCANs

##use markers enriched in ChC clusters based on RNA data
ChC.markers <- read.delim("Mouse_Pvalb_Vipr2_ChC_markers_pvaladj_0.05.txt", header = FALSE)

peak_anno_df <- peak_anno_df[peak_anno_df$gene %in% ChC.markers$V1,]
peak_anno_df$site <- rownames(peak_anno_df)

DefaultAssay(msMOp.atac) <- "ATAC"

#find gene motif enrichments (genes having >1 site)
tab <- table(peak_anno_df$gene)
unique.genes <- names(tab[tab > 1])

gene.tfbs.list <- lapply(unique.genes, function(cl) {
  print(paste("Running for gene:", cl))
  
  cl.sites <- rownames(peak_anno_df[peak_anno_df$gene %in% cl,])
  EM <- FindMotifs(object = msMOp.atac,
                   features = cl.sites)
  EM <- EM[EM$pvalue < 0.05,]
  EM
})
names(gene.tfbs.list) <- unique.genes
gene.tfbs <- do.call("rbind", lapply(gene.tfbs.list, as.data.frame)) 
marker<-unlist(lapply(rownames(gene.tfbs),function(x) unlist(strsplit(x,"[.]"))[1]))
gene.tfbs$gene <- marker

#find active TFBS in ChC cells over Basket Cells
DefaultAssay(msMOp.atac) <- "chromvar"
tf.markers <- FindMarkers(
  object = msMOp.atac, 
  ident.1 = "Pvalb_Vipr2",
  ident.2 = c("Pv_Ntf3_Trim63","Pv_Tac1"),
  features = unique(gene.tfbs$motif),
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'nCount_ATAC'
)
tf.markers$motif <- rownames(tf.markers)
tf.markers <- left_join(gene.tfbs, tf.markers, by = c("motif"))
rownames(tf.markers) <- rownames(gene.tfbs)
ms.tf.markers <-na.omit(tf.markers)
ms.tf.markers

tf.mark <- ms.tf.markers[ms.tf.markers$p_val_adj < 0.05,]
tf.mark <- distinct(tf.mark, motif, .keep_all = TRUE)
tf.mark$cluster <- "1"
tf.mark %>% group_by(cluster) %>% top_n(20, avg_logFC) -> top20 

sites.use <- top20$motif
ave.sites <- AverageExpression(msMOp.atac, assays = "chromvar", features = sites.use, slot = "data" )

names.sites <- top20$motif.name
ave.sites <- ave.sites$chromvar
rownames(ave.sites) <- names.sites

#Plot TFBSs that are enriched in ChC marker genes and active in the ChC cluster
ggHeat(ave.sites, rescaling = "row", clustering = "row", x.lab.size = 8, y.lab.size = 11,
) + scale_fill_viridis()





# Section 18: Human Marmoset SNARE-Seq2 AC - Subclass TFBS Conservation ------------------
library(Signac)
library(Seurat)
library(dplyr)
library(Matrix)
library(methods)
library(swne)
library(ggplot2)
library(viridis)
set.seed(1234)

#load marmoset data 
#load marMOp.atac

DefaultAssay(marMOp.atac) <- "chromvar"
Idents(object = marMOp.atac) <- "subclass"
order <- c("Lamp5", "Sncg", "Meis2", "Vip", "Sst Chodl", "Sst", "Pvalb", 
           "L2-3 IT", "L5 IT", "L5 ET", "L5-6 NP", "L6 IT", "L6 IT Car3", "L6 CT", "L6b",                                          
           "OPC", "Astro", "Oligo", "Micro-PVM", "VLMC","Peri", "Endo")

Idents(object = marMOp.atac) <- factor(Idents(object = marMOp.atac), levels = order)
#remove clusters not present in human, too sparse or non-neuronal
marMOp.atac <- subset(marMOp.atac, idents = c("Meis2", "Sst Chodl","OPC", "Astro", "Oligo", "Micro-PVM", "VLMC","Peri", "Endo"), invert = TRUE)


#load Human data
#load MOp.atac
DefaultAssay(MOp.atac) <- "chromvar"
Idents(object = MOp.atac) <- "subclass"
order <- c("LAMP5", "SNCG", "VIP", "SST CHODL", "SST", "PVALB", 
           "L2-3 IT", "L5 IT", "L5 ET", "L5-6 NP", "L6 IT", "L6 IT Car3", "L6 CT", "L6b",                                          
           "OPC", "Astro", "Oligo", "Micro-PVM", "VLMC", "Endo")  

Idents(object = MOp.atac) <- factor(Idents(object = MOp.atac), levels = order)
#remove clusters not present in marmoset
MOp.atac <- subset(MOp.atac, idents = c("SST CHODL","OPC", "Astro", "Oligo", "Micro-PVM", "VLMC", "Endo"), invert = TRUE)


#Get average scaled chromvar values - Human
ave.tf.hs <- AverageExpression(MOp.atac, assays = "chromvar", features = rownames(MOp.atac), slot = "data" )

PFMatrixList <- MOp.atac@assays$ATAC@misc$motif@motif.names
PFMatrix <- do.call("rbind", lapply(PFMatrixList, as.data.frame)) 
colnames(PFMatrix) <- "motif.names"
PFMatrix <- as.character(PFMatrix[sites.use,])
rownames(ave.tf$chromvar) <- PFMatrix

scaled <- t(scale(t(ave.tf.hs$chromvar)))
scaled <- scale(scaled)
range(scaled)
scaled[scaled < 0] <- 0
scaled[scaled > 4] <- 4

tf.hs <- scaled

#Get average scaled chromvar values - Marmoset
ave.tf.mar <- AverageExpression(marMOp.atac, assays = "chromvar", features = rownames(marMOp.atac), slot = "data" )

PFMatrixList <- MOp.atac@assays$ATAC@misc$motif@motif.names
PFMatrix <- do.call("rbind", lapply(PFMatrixList, as.data.frame)) 
colnames(PFMatrix) <- "motif.names"
PFMatrix <- as.character(PFMatrix[sites.use,])
rownames(ave.tf$chromvar) <- PFMatrix

scaled <- t(scale(t(ave.tf.mar$chromvar)))
scaled <- scale(scaled)
range(scaled)
scaled[scaled < 0] <- 0
scaled[scaled > 4] <- 4

tf.mar <- scaled

library("corrplot")
ave.tf.cor<-cor(tf.hs, tf.mar, use = "complete.obs")
col2 <- colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#FFFFFF","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))
corrplot(ave.tf.cor, col = col2(500), order = "hclust", hclust.method = "ward.D2", cl.lim=c(-1,1))


#Correlation on Gene activities
int.genes <- read.csv("marm_hum_integration_genes_inh_exc.csv")
DefaultAssay(MOp.atac) <- "ACTIVITY"
DefaultAssay(marMOp.atac) <- "ACTIVITY"
common.int.genes <- intersect(int.genes$x, intersect(rownames(MOp.atac),rownames(marMOp.atac)))

ave.tf.hs <- AverageExpression(MOp.atac, assays = "ACTIVITY", features = common.int.genes, slot = "scale.data" )
tf.hs <- ave.tf.hs$ACTIVITY

ave.tf.mar <- AverageExpression(marMOp.atac, assays = "ACTIVITY", features = common.int.genes, slot = "scale.data" )
tf.mar <- ave.tf.mar$ACTIVITY

common.genes <- intersect(rownames(tf.hs), rownames(tf.mar))
hu.order <- c("VIP","LAMP5","SNCG","SST","PVALB","L6 IT Car3","L2-3 IT","L6 IT","L5 ET",
              "L6b","L5-6 NP","L5 IT","L6 CT")
mar.order <- c("Vip","Lamp5","Sncg","Sst","Pvalb","L6 IT Car3","L2-3 IT","L6 IT","L5 ET",
               "L6b","L5-6 NP","L5 IT","L6 CT")
library("corrplot")
ave.tf.cor<-cor(tf.hs[common.genes,hu.order], tf.mar[common.genes,mar.order])
col2 <- colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#FFFFFF","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))
corrplot(ave.tf.cor, col = col2(500), order = "hclust", hclust.method = "ward.D2", cl.lim=c(-1,1))



#Correlation on gene expression
int.genes <- read.csv("marm_hum_integration_genes_inh_exc.csv") #gene set used for inh exc cross species integration
DefaultAssay(MOp.atac) <- "RNA"
DefaultAssay(marMOp.atac) <- "RNA"
common.int.genes <- intersect(int.genes$x, intersect(rownames(MOp.atac),rownames(marMOp.atac)))

ave.tf.hs <- AverageExpression(MOp.atac, assays = "RNA", features = common.int.genes, slot = "scale.data" )
tf.hs <- ave.tf.hs$RNA

ave.tf.mar <- AverageExpression(marMOp.atac, assays = "RNA", features = common.int.genes, slot = "scale.data" )
tf.mar <- ave.tf.mar$RNA

common.genes <- intersect(rownames(tf.hs), rownames(tf.mar))
hu.order <- c("VIP","LAMP5","SNCG","SST","PVALB","L6 IT Car3","L2-3 IT","L6 IT","L5 ET",
              "L6b","L5-6 NP","L5 IT","L6 CT")
mar.order <- c("Vip","Lamp5","Sncg","Sst","Pvalb","L6 IT Car3","L2-3 IT","L6 IT","L5 ET",
               "L6b","L5-6 NP","L5 IT","L6 CT")
library("corrplot")
ave.tf.cor<-cor(tf.hs[common.genes,hu.order], tf.mar[common.genes,mar.order])
col2 <- colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#FFFFFF","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))
corrplot(ave.tf.cor, col = col2(500), order = "hclust", hclust.method = "ward.D2", cl.lim=c(-1,1))











# Section 19: Human Marmoset SNARE-Seq2 AC - TF Motif Families -----------

##Plots for TF Motif Families
library(stringr)
motif.clust <- read.delim("clusters_motif_names.tab.txt",sep="\t",header=FALSE) #available from http://jaspar.genereg.net/downloads/

cl.motif.list <- lapply(clusters, function(cl) {
  print(paste("Running for cluster:", cl))
  
  cluster <- str_split(motif.clust[motif.clust$V1 == cl,2],",", simplify = FALSE)
  cluster <- as.data.frame(cluster); colnames(cluster) <- "motif.name"
  cluster$cluster <- as.character(rep(cl, length(cluster$motif.name)))
  cluster
  
})
names(cl.motif.list) <- clusters
cl.motif <- do.call("rbind", lapply(cl.motif.list, as.data.frame)) 
cl.motif$motif.name <- gsub("_var_2_", "(var.2)",cl.motif$motif.name)
cl.motif$motif.name <- gsub("_var_3_", "(var.3)",cl.motif$motif.name)


#prepare human seurat object
DefaultAssay(MOp.atac) <- "chromvar"
Idents(object = MOp.atac) <- "subclass"
order <- c("LAMP5", "SNCG", "VIP", "SST CHODL", "SST", "PVALB", 
           "L2-3 IT", "L5 IT", "L5 ET", "L5-6 NP", "L6 IT", "L6 IT Car3", "L6 CT", "L6b",                                          
           "OPC", "Astro", "Oligo", "Micro-PVM", "VLMC", "Endo")  
Idents(object = MOp.atac) <- factor(Idents(object = MOp.atac), levels = order)
Idents(MOp.atac)


#prepare marmoset seurat object
DefaultAssay(marMOp.atac) <- "chromvar"
Idents(object = marMOp.atac) <- "subclass"
order <- c("Lamp5", "Sncg", "Meis2", "Vip", "Sst Chodl", "Sst", "Pvalb", 
           "L2-3 IT", "L5 IT", "L5 ET", "L5-6 NP", "L6 IT", "L6 IT Car3", "L6 CT", "L6b",                                          
           "OPC", "Astro", "Oligo", "Micro-PVM", "VLMC","Peri", "Endo")
Idents(object = marMOp.atac) <- factor(Idents(object = marMOp.atac), levels = order)



#Human Average Scaled by Cluster
hu.cl.motif <- cl.motif[cl.motif$motif.name %in% rownames(MOp.atac),]
rownames(hu.cl.motif) <- hu.cl.motif$motif.name
ave.tf.hs <- AverageExpression(MOp.atac, assays = "chromvar", features = rownames(MOp.atac), slot = "data" )
ave.tf.hs <- ave.tf.hs$chromvar
hu.cl.motif <- hu.cl.motif[rownames(ave.tf.hs),]
ave.tf.hs$cl <- hu.cl.motif$cluster
ave.tf.hs <- na.omit(ave.tf.hs)

ave.tf.hs <- group_by(ave.tf.hs, cl) %>% summarize_if(is.numeric, mean)
ave.tf.hs <- data.frame(ave.tf.hs); rownames(ave.tf.hs) <- ave.tf.hs$cl; ave.tf.hs <- ave.tf.hs[,-1]
scaled <- t(scale(t(ave.tf.hs)))
scaled <- scale(scaled)
range(scaled)
scaled[scaled < 0] <- 0
scaled[scaled > 4] <- 4

tf.hs <- scaled

#Marmoset Average Scaled by Cluster
mar.cl.motif <- cl.motif[cl.motif$motif.name %in% rownames(marMOp.atac),]
rownames(mar.cl.motif) <- mar.cl.motif$motif.name
ave.tf.mar <- AverageExpression(marMOp.atac, assays = "chromvar", features = rownames(marMOp.atac), slot = "data" )
ave.tf.mar <- ave.tf.mar$chromvar
mar.cl.motif <- mar.cl.motif[rownames(ave.tf.mar),]
ave.tf.mar$cl <- mar.cl.motif$cluster
ave.tf.mar <- na.omit(ave.tf.mar)

ave.tf.mar <- group_by(ave.tf.mar, cl) %>% summarize_if(is.numeric, mean)
ave.tf.mar <- data.frame(ave.tf.mar); rownames(ave.tf.mar) <- ave.tf.mar$cl; ave.tf.mar <- ave.tf.mar[,-1]
scaled <- t(scale(t(ave.tf.mar)))
scaled <- scale(scaled)
range(scaled)
scaled[scaled < 0] <- 0
scaled[scaled > 4] <- 4

tf.mar <- scaled

colnames(tf.mar) <- paste("mar", colnames(tf.mar), sep = "_")
tf.hs.mar <- cbind(tf.hs, tf.mar[rownames(tf.hs),])
tf.hs.mar <- data.frame(tf.hs.mar)
tf.hs.mar$cluster <- gsub("cluster_","cl",rownames(tf.hs.mar))


#Generate plot for Lamp5 subclass
to.use <- c("cluster_60","cluster_34","cluster_12","cluster_33","cluster_98","cluster_11","cluster_48","cluster_47",
            "cluster_1","cluster_49","cluster_63","cluster_74","cluster_66","cluster_25","cluster_9","cluster_3",
            "cluster_89","cluster_21","cluster_38","cluster_10","cluster_37","cluster_14","cluster_27","cluster_23")
lamp5 <- tf.hs.mar[to.use,c("LAMP5","mar_Lamp5")]
rownames(lamp5) <- to.use
colnames(lamp5) <- c("human","marmoset")
ggHeat(lamp5, rescaling = "none", clustering = "none", x.lab.size = 8, y.lab.size = 11,
       heatscale = c(low = "skyblue", mid ="white", high = "#005a32")
) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  ggtitle("LAMP5")


#Generate plot for L5 IT subclass
L5IT <- tf.hs.mar[to.use,c("L5.IT","mar_L5.IT")]
rownames(L5IT) <- to.use
colnames(L5IT) <- c("human","marmoset")
ggHeat(L5IT, rescaling = "none", clustering = "none", x.lab.size = 8, y.lab.size = 11,
       heatscale = c(low = "skyblue", mid ="white", high = "#005a32")
) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  ggtitle("L5-IT")


