#Code Underlying SNARE-Seq2 M1 Data Analyses 
#Blue B. Lake
#b1lake@eng.ucsd.edu

# Mouse ATAC TFBS Activities ------------------------------------
library(Signac)
library(Seurat)
library(dplyr)
library(Matrix)
library(methods)
library(swne)
library(cicero)
library(ggplot2)


####
##Prepare Seurat Object
####

counts <- readMM("~/NeMO_analysis_folder/SNARE/Analysis/Mouse_ATAC/MOp_intergration.subMGE.batchCorrect.pmat.mtx")
colnames <- read.delim("~/NeMO_analysis_folder/SNARE/Analysis/Mouse_ATAC/MOp_intergration.subMGE.batchCorrect.pmat.ygi", header = FALSE)
rownames <- read.delim("~/NeMO_analysis_folder/SNARE/Analysis/Mouse_ATAC/MOp_intergration.subMGE.batchCorrect.pmat.xgi", header = FALSE)
colnames(counts) <- colnames$V1; rownames(counts) <- rownames$V1
counts <- t(counts)
labels <- read.delim("~/NeMO_analysis_folder/SNARE/Analysis/Mouse_ATAC/MOp_intergration.subMGE.batchCorrect.labels.meta.txt")
rownames(labels) <- paste(labels$sample,labels$barcode, sep = ".")
labels <- labels[colnames(counts),]

msMOp.atac <- CreateSeuratObject(counts = counts, assay = "ATAC", project = "MOp_ATAC")
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
mm10.chr.lengths <- read.table("mm10.chrom.sizes", header = F, sep = "\t") 
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
input_cds <- suppressWarnings(new_cell_data_set(count.matrix,
                                                cell_metadata = pData,
                                                gene_metadata = fData))

cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = umap.coordinates, k = 50)
conns <- run_cicero(cicero_cds, mm10.chr.lengths) 

## Assign peaks to modules
ccan.assigns <- generate_ccans(conns, coaccess_cutoff = 0.08)


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

file = "~/NeMO_analysis_folder/SNARE/Analysis/Human/JASPAR2020_CORE_vertebrates_non-redundant_pfms_jaspar.txt"
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

save(msMOp.atac, file = "Mouse_MOp_scATAC_PVALB_seurat.rda")



####
##TFBS Activities ChC vs Basket
####

load("~/NeMO_analysis_folder/SNARE/Analysis/Mouse_ATAC/Mouse_MOp_scATAC_PVALB_seurat.rda")
DefaultAssay(msMOp.atac) <- "chromvar"
msMOp.atac <- RenameIdents(msMOp.atac, "Pv_Ntf3_Trim63" = "BC",
                           "Pv_Tac1" = "BC", "Pvalb_Vipr2" = "ChC")

tf.markers <- FindAllMarkers(
  object = msMOp.atac,
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'nCount_ATAC'
)
write.table(tf.markers, file="Zhang_Mouse_MOp_PVALB_TF_Activity_Marker_table.txt", sep = "\t", row.names=TRUE, col.names=TRUE)

#Plot Figure 6e
DotPlot(msMOp.atac, features = c("NFIB", "RORA")) + RotatedAxis()



