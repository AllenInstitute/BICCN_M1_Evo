#Code Underlying SNARE-Seq2 M1 Data Analyses 
#Blue B. Lake
#b1lake@eng.ucsd.edu

# Marmoset SNARE-Seq2 Data Processing --------------------------
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
#peak counts were performed according to human analyses with any 
#exceptions outlined in the manuscript methods


####
##Preparation of RNA Seurat object
####


marMOp.gmcf <- readRDS("~/NeMO_analysis_folder/SNARE/Analysis/Marmoset/BICCN-h_Marmoset_SNARE2-R_UMI_dTn6_EmptyBC_DD_Gene-UMI_Filter_12022019.rds")
marMOp <- CreateSeuratObject(counts = marMOp.gmcf, project = "Marmoset M1", min.cells = 3, min.features = 200)

#Normalize the data
marMOp <- NormalizeData(object = marMOp, normalization.method = "LogNormalize", 
                        scale.factor = 10000)

#Detection of variable genes across the single cells
marMOp <- FindVariableFeatures(object = marMOp, selection.method = "vst", nfeatures = 2000)

#Scaling the data and removing unwanted sources of variation
all.genes <- rownames(x = marMOp)
marMOp <- ScaleData(object = marMOp, features = all.genes, vars.to.regress = c("orig.ident"))

#Perform linear dimensional reduction
marMOp <- RunPCA(object = marMOp, features = VariableFeatures(object = marMOp), npcs = 75)
marMOp <- FindNeighbors(object = marMOp, dims = 1:75)
marMOp <- FindClusters(object = marMOp, resolution = 1)
marMOp <- RunUMAP(object = marMOp, dims = 1:75)
DimPlot(object = marMOp, reduction = "umap")

#Update metadata tables
meta <- read.delim("~/NeMO_analysis_folder/SNARE/Analysis/Marmoset/Marmoset_cell_assignment.txt", sep = "\t", header=TRUE, row.names = 1)
meta <- meta[rownames(marMOp@meta.data),]
marMOp <- AddMetaData(marMOp, metadata = meta)

#Identify possible multiplets from seurat clustering
library("scrattch.hicat")
final<-marMOp$consensus_cluster
names(final)<-rownames(marMOp@meta.data)
final <- factor(final)
final <- na.omit(final)
final <- factor(final)
final.df <- data.frame(cluster_id = unique(final),
                       cluster_label = paste0("cons_",unique(final)),
                       cluster_color = rainbow(length(unique(final))))
rownames(final.df) <- final.df$cluster_id
final.df.order<-final.df[order(final.df$cluster_id),]
prop<-Idents(object = marMOp)
compare.result <- compare_annotate(prop, final, final.df.order, reorder = TRUE)
compare.result$g + scale_x_discrete(name ="RNA Cluster", labels=compare.result$cl.id.map$old)

#Remove clusters that show mapping to multiple cell types
marMOp <- subset(marMOp, idents = c(7,9,20), invert = TRUE)



####
##Preparation of RNA/AC Seurat object
####

#load final AC peak by cell count matrix 
load("~/NeMO_analysis_folder/SNARE/Analysis/Marmoset/BICCN-h_20190730_20190903_Marmoset_MOp_peak_matrix_12062019.rda")
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

meta = read.table("~/NeMO_analysis_folder/SNARE/Raw_Data/Marmoset/Zhang_BICCN-H_20190730_20190903_marMOp_Final_Sample_Metadata.txt", sep="\t",header=TRUE, row.names = 1)
meta <- meta[colnames(marMOp.atac), ]
marMOp.atac <- AddMetaData(marMOp.atac, metadata = meta)

#pre-processing
DefaultAssay(marMOp.atac) <- "ATAC"
VariableFeatures(marMOp.atac) <- names(which(Matrix::rowSums(marMOp.atac) > 20))
marMOp.atac <- RunLSI(marMOp.atac, n = 50, scale.max = NULL) #latent semantic indexing
marMOp.atac <- RunUMAP(marMOp.atac, reduction = "lsi", dims = 1:40,
                       n.neighbors = 50L, min.dist = 0.5)


###Add in RNA data
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
file = "~/NeMO_analysis_folder/SNARE/Analysis/Human/JASPAR2020_CORE_vertebrates_non-redundant_pfms_jaspar.txt"
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
library(chromfunks)

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
calJac3.chr.lengths <- read.table("~/NeMO_analysis_folder/SNARE/Analysis/Marmoset/calJac3.chrom.sizes", header = F, sep = "\t") ## A tab separated text file with chromosome lengths
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
input_cds <-  suppressWarnings(new_cell_data_set(count.matrix,
                                                 cell_metadata = pData,
                                                 gene_metadata = fData))

cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = umap.coordinates, k = 50)
conns <- run_cicero(cicero_cds, calJac3.chr.lengths) 

## Assign peaks to modules
ccan.assigns <- generate_ccans(conns)

#generate annotation file
gtf <- rtracklayer::import('~/NeMO_analysis_folder/SNARE/Analysis/Marmoset/refdata-cellranger-Callithrix_jacchus-3.2-GCF_000004665.1_genes.gtf')
ch.conv <- read.delim("~/NeMO_analysis_folder/SNARE/Analysis/Marmoset/calJac3.2_ucscToRefSeq.txt.gz")
rownames(ch.conv) <- ch.conv$name
ch.conv <- ch.conv[seqlevels(gtf),]
seqlevels(gtf) <- as.character(ch.conv$X.chrom)
gene.coords <- gtf[gtf$type == 'gene']
seqlevelsStyle(gene.coords) <- 'UCSC'
gene.coords <- keepStandardChromosomes(gene.coords, pruning.mode = 'coarse')

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
  annotation = gene.coords,
  sep = c(':', '-')
)

peak_anno[peak_anno$distance > 5000,]$gene <- NA

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

save(marMOp.atac, file = "Zhang_BICCN-H_20190730_20190903_marMOp_Seurat.rda")

