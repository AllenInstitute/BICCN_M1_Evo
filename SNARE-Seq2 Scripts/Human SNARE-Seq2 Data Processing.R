#Code Underlying SNARE-Seq2 M1 Data Analyses 
#Blue B. Lake
#b1lake@eng.ucsd.edu

# Human SNARE-Seq2 Data Processing ---------------------------------


####
##Generate Combined RNA count matrices, Apply QC filters
####

##Example code, raw library .rds files not provided
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

#process remaining libraries as above...


###DoubletDetection Filtering

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

#Process remaining samples as above



###Combine all Samples
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


###Apply Gene UMI Filters
library(Seurat)
library(dplyr)
library(Matrix)

MOp <- CreateSeuratObject(counts = MOp, project = "Human M1", min.cells = 3, min.features = 200)

#Gene Filter >200 and <7500 genes
MOp <- subset(x = MOp, subset = nFeature_RNA > 200 & nFeature_RNA < 7500)

#Gene/UMI Ratio Filter to remove low qualty barcodes
countMatrix <- GetAssayData(object = MOp, slot = "counts")
MOp.gmcf <- gene.vs.molecule.cell.filter(countMatrix,min.cell.size=200)
saveRDS(MOp.gmcf, file = paste0(dir3, "BICCN-h_hMOp_SNARE-R_UMI_dTn6_EmptyBC_MT_DD_Gene-UMI_Filter_07252019.rds"))









####
##Pagoda2 clustering
####

library(Seurat)
library(dplyr)
library(Matrix)
library(pagoda2)
require(parallel)

#Create Seurat object with subsampled oligodendrocytes (max 5000)
oli.sub.counts <- readRDS("~/NeMO_analysis_folder/SNARE/Analysis/Human/BICCN-h_hMOp_SNARE-R_UMI_dTn6_EmptyBC_MT_DD_Gene-UMI_Filter_Oli5000.rds")
MOp <- CreateSeuratObject(counts = oli.sub.counts, project = "Human M1", min.cells = 3, min.features = 200)

#Normalize counts
MOp <- NormalizeData(object = MOp, normalization.method = "LogNormalize", 
                     scale.factor = 10000)

#Detection of variable genes across the single cells
MOp <- FindVariableFeatures(object = MOp, selection.method = "vst", nfeatures = 3000)

#Scaling the data and removing unwanted sources of variation
all.genes <- rownames(x = MOp)
MOp <- ScaleData(object = MOp, features = all.genes, vars.to.regress = c("orig.ident"))

#Perform linear dimensional reduction
MOp <- RunPCA(object = MOp, features = VariableFeatures(object = MOp), npcs = 75)
MOp <- RunUMAP(object = MOp, dims = 1:75)

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

saveRDS(p2,file = "BICCN-h_hMOp_SPL-R_UMI_dTn6_EmptyBC_MT_DD_Gene_Filter_Pagoda_75PC_k100_2.RDS")

#Repeat makeKnnGraph and getKnnClusters for different k values

k50infomap <- p2$clusters$PCA$infomap
k100infomap <- p2$clusters$PCA$infomap
k200infomap <- p2$clusters$PCA$infomap
k500infomap <- p2$clusters$PCA$infomap


#Transfer Pagoda2 Clusters and PCA values (k = 100) to Seurat
p2 <- readRDS("~/NeMO_analysis_folder/SNARE/Analysis/Human/BICCN-h_hMOp_SPL-R_UMI_dTn6_EmptyBC_MT_DD_Gene_Filter_Pagoda_75PC_k100_2.RDS")
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

to.keep <- to.keep[!duplicated(to.keep)]
MOp <- subset(MOp, cells = to.keep)

#check for marker genes
p2 #k = 100
var.info <- p2$misc$varinfo; var.info <- var.info[order(var.info$lp, decreasing = F),];
sn.od.genes <- rownames(var.info)[1:6673]

MOp.markers <- FindAllMarkers(MOp, only.pos = TRUE, features = sn.od.genes,
                              logfc.threshold = 0.25, min.pct = 0.25)



####
##Compare local (Pagoda2) clusters with Smart-Seq v4 (SSv4) clusters
####

### Extended Data Fig. 4d
###Compare SMARTer Clusters with SNARE-Seq2 Pagoda2 Clusters - marker genes
load("~/NeMO_analysis_folder/SNARE/Analysis/Human/ABI_Smart-seq_MOp_Seuratv3.Robj")
abi.markers <- read.delim("~/NeMO_analysis_folder/SNARE/Analysis/Human/Lein_Smart-seq_Markers_betaprop1_0-4.txt", header = FALSE)
abi.markers <- abi.markers$V1

p2 #k = 100
var.info <- p2$misc$varinfo; var.info <- var.info[order(var.info$lp, decreasing = F),];
sn.od.genes <- rownames(var.info)[1:2000]

common.genes<-c(sn.od.genes, abi.markers)  
common.genes<-common.genes[common.genes %in% rownames(x = abi)]
common.genes<-common.genes[common.genes %in% rownames(x = MOp)]  #1808

ave.MOp<-AverageExpression(MOp, features = common.genes, slot = "scale.data")
ave.abi<-AverageExpression(abi, features = common.genes, slot = "scale.data")

library("corrplot")
ave.mop.cor<-cor(cbind(ave.MOp$RNA,ave.abi$RNA))
order.hc2 <- corrMatOrder(ave.mop.cor, order = "hclust", hclust.method = "ward.D")
ave.mop.cor<-ave.mop.cor[order.hc2,order.hc2]

x <- ave.mop.cor[rownames(ave.mop.cor) %in% colnames(ave.MOp$RNA),
                 colnames(ave.mop.cor) %in% colnames(ave.abi$RNA)]

#Refine ordering
order1 <- as.character(c(
  31,38,33,37,30,32,35,28,34,29,36,23,22,24,20,21,27,26,25,5,1,4,3,2,9,8,6,7,
  19,16,10,11,14,18,12,17,13,15,
  55,56,58,57,59,52,51,53,49,48,50,47,54,44,46,45,43,39,41,42,40,
  65,66,62,63,64,61,60,73,67,68,70,72,69,71
))
order2 <- colnames(x)[c(45:116,1:44,117:129)]

x=x[order1,order2]
write.table(x, file="Zhang_BICCN-H_20190523-20190611_huMOp_Corr_Matrix_P2-ABI.txt", sep = "\t", row.names=TRUE, col.names=TRUE)
col2 <- colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#FFFFFF","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))
corrplot(x, col = col2(500),  method="color",order = "original",
         hclust.method = "ward.D", cl.lim=c(-1,1))






####
##QC filter SNARE-seq2 AC data
####


library(SnapATAC)
dir <- "~/NeMO_analysis_folder/SNARE/Raw_Data/Human"

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
rna.table <- read.table("~/NeMO_analysis_folder/SNARE/Analysis/Human/BICCN-h_hMOp_Prelim_RNA_Consensus_Meta.txt", sep="\t",header=TRUE, row.names = 1)
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

promoter.df = read.table("~/NeMO_analysis_folder/SNARE/Analysis/Human/hg38.promoters.bed"); 
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






####
##Peak Calling - round 1
####

library(SnapATAC)
library(parallel)
library(GenomicRanges);
library(ggplot2)
library(dplyr)

meta <- meta[x.sp@barcode,]

###RNA-level Clusters
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






####
##RNA/AC - Integrative Analysis
####
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)

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
MOp.atac <- AddMetaData(MOp.atac, metadata = meta)
MOp.atac$tech <- "atac"

#pre-processing
DefaultAssay(MOp.atac) <- "ATAC"
VariableFeatures(MOp.atac) <- names(which(Matrix::rowSums(MOp.atac) > 100))
MOp.atac <- RunLSI(MOp.atac, n = 50, scale.max = NULL) #latent semantic indexing
MOp.atac <- RunUMAP(MOp.atac, reduction = "lsi", dims = 1:50)

#Update metadata for SNARE_Seq2 RNA Seurat object (MOp)
MOp$tech <- "rna"



###Generate gene activity matrix using cisero
library(cicero)

#Cicero analysis: generate predicted gene activity matrix from chrom data
Idents(object = MOp.atac) <- "consensus_cluster"

## Run cicero
hg38.chr.lengths <- read.table("~/NeMO_analysis_folder/SNARE/Analysis/Human/hg38.chr.lengths.txt", header = F, sep = "\t") 
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
input_cds <- suppressWarnings(new_cell_data_set(count.matrix,
                                                cell_metadata = pData,
                                                gene_metadata = fData))

cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = umap.coordinates, k = 50)
conns <- run_cicero(cicero_cds, hg38.chr.lengths) 

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
p2 <- readRDS("~/NeMO_analysis_folder/SNARE/Analysis/Human/BICCN-h_hMOp_SPL-R_UMI_dTn6_EmptyBC_MT_DD_Gene_Filter_Pagoda_75PC_k100_2.RDS")
var.info <- p2$misc$varinfo; var.info <- var.info[order(var.info$lp, decreasing = F),];
variable.genes <- rownames(var.info)[1:2000]

DefaultAssay(MOp.atac) <- "ACTIVITY"
abi.markers <- read.delim("~/NeMO_analysis_folder/SNARE/Analysis/Human/Lein_Smart-seq_Markers_betaprop1_0-4.txt", header = FALSE)
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


#Generate predictions for RNA consensus clusters
cl.predictions <- TransferData(anchorset = transfer.anchors, refdata = as.character(MOp$consensus_cluster), 
                               weight.reduction = MOp.atac[["lsi"]])
pred.meta <- cl.predictions[,c(1,130)]
colnames(pred.meta) <- c("cons-cl.predicted.id","cons-cl.prediction.score.max")
MOp.atac <- AddMetaData(MOp.atac, metadata = pred.meta)





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
coembed <- RenameIdents(coembed, '10' = '8', '18' = '8', '26' = '20', '29' = '27', '30' = '0', '31' = '11', '37' = '12',
                        '42' = '20', '45' = '12', '46' = '4', '5' = '20', '58' = '53', '59' = '0', '60' = '44', 
                        '72' = '47', '74' = '67', '76' = '20', '79' = '77', '82' = '6', '92' = '75', '93' = '0')
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


#Match RNA/AC integrated clusters with RNA-level clusters by jaccard similarity and frequency 
comp.1n <- compare.result$tb.df
comp.1n$cl <- as.character(comp.1n$cl)
comp.1n.id <- compare.result$cl.id.map
comp.1n.id$new <- as.character(comp.1n.id$new)
comp.1n <- left_join(comp.1n, comp.1n.id, 
                     by = c("cl" = "new"))
comp.1n$cl <- comp.1n$old
comp.1n %>% group_by(ref.cl) %>% top_n(2, jaccard) -> top2.jacc #select top two consensus clusters by jaccard similarity
top2.jacc %>% group_by(ref.cl) %>% top_n(1, Freq) -> top1 #select final top cluster by frequency
#Using these as a guide, consensus clusters were merged while accounting for subclass labels 


###Clean up coembed metadata
columns.keep <- c("library","nCount_RNA","nFeature_RNA","experiment_short", "patient",
                  "class","subclass","consensus_cluster","consensus_cluster_number",
                  "nCount_ATAC","nFeature_ATAC","nCount_ACTIVITY","nFeature_ACTIVITY",
                  "predicted.id","prediction.score.max","cons-cl.predicted.id",
                  "cons-cl.prediction.score.max","RNA_snn_res.4","res4_merged",
                  "atac_cons_clusters","tech")
coembed@meta.data <- coembed@meta.data[,columns.keep]

columns.renamed <- c("library","nCount_RNA","nFeature_RNA","experiment_short", "patient",
                     "class","subclass","RNA_cluster","RNA_cluster_ID",
                     "nCount_ATAC","nFeature_ATAC","nCount_ACTIVITY","nFeature_ACTIVITY",
                     "subclass.predicted.id","subclass.prediction.score.max","RNA-cl.predicted.id",
                     "RNA-cl.prediction.score.max","RNA-AC_snn_res.4","RNA-AC_snn_res.4.merged",
                     "AC_cluster","tech")

colnames(coembed@meta.data) <- columns.renamed

save(coembed, file = "Zhang_BICCN-H_20190523-20190611_MOp_RNA-AC_Coembed_Seurat_Final.rda")



####
##Peak Calling - round2 
####

library(SnapATAC)
library(parallel)
library(GenomicRanges);
library(ggplot2)
library(dplyr)


meta = read.table("~/NeMO_analysis_folder/SNARE/Raw_Data/Human/Zhang_BICCN-H_20190523-20190611_huMOp_Final_Sample_Metadata.txt", sep="\t",header=TRUE, row.names = 1)
meta <- meta[x.sp@barcode,]

###AC CLuster Peak calling
clust <- meta$AC_cluster_tree_order #AC Level clusters
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


###RNA-level Cluster Peak calling
clust <- meta$RNA_cluster
names(clust)<- rownames(meta)
clust <- factor(clust)
x.sp@cluster <- clust
clusters.sel = names(table(x.sp@cluster))[which(table(x.sp@cluster) > 100)];
peaks.ls = mclapply(seq(clusters.sel), function(i){
  print(clusters.sel[i]);
  peaks = runMACS(
    obj=x.sp[which(x.sp@cluster==clusters.sel[i]),], 
    output.prefix=paste0("RNA_clusters.", gsub(" ", "_", clusters.sel)[i]),
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





####
## Generate Final Combined Seurat object 
####

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

#final combined peak by cell barcode count matrix
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
meta = read.table("~/NeMO_analysis_folder/SNARE/Raw_Data/Human/Zhang_BICCN-H_20190523-20190611_huMOp_Final_Sample_Metadata.txt", sep="\t",header=TRUE, row.names = 1)
meta <- meta[colnames(MOp.atac), ]
MOp.atac <- AddMetaData(MOp.atac, metadata = meta)
MOp.atac$tech <- "atac"

#pre-processing
DefaultAssay(MOp.atac) <- "ATAC"
VariableFeatures(MOp.atac) <- names(which(Matrix::rowSums(MOp.atac) > 100))
MOp.atac <- RunLSI(MOp.atac, n = 50, scale.max = NULL) #latent semantic indexing
MOp.atac <- RunUMAP(MOp.atac, reduction = "lsi", dims = 1:50)


###Add in motif matrix
file = "~/NeMO_analysis_folder/SNARE/Analysis/Human/JASPAR2020_CORE_vertebrates_non-redundant_pfms_jaspar.txt"
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
Idents(object = MOp.atac) <- "AC_cluster_tree_order"
Idents(object = MOp.atac) <- factor(Idents(object = MOp.atac), levels = 1:43)

## Run cicero
hg38.chr.lengths <- read.table("~/NeMO_analysis_folder/SNARE/Raw_Data/Human/hg38.chr.lengths.txt", header = F, sep = "\t") ## A tab separated text file with chromosome lengths
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
input_cds <- suppressWarnings(new_cell_data_set(count.matrix,
                                                cell_metadata = pData,
                                                gene_metadata = fData))

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

save(MOp.atac, file = "Zhang_BICCN-H_20190523-20190611_huMOp_Seurat.rda")




