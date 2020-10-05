library(Seurat)

data <- readRDS("http://data.nemoarchive.org/biccn/lab/lein/2020_M1_study_analysis/Transcriptomics/cross_species_integration/sample.combined_exc_4_species_integration.RDS")
data <- RunUMAP(data, dims = 1:100)
DimPlot(data, reduction = "umap", group.by = "orig.ident")

Idents(data) <- data$subclass_label
data_et <- subset(data, idents = "L5 ET")
DimPlot(data_et, reduction = "umap", group.by = "orig.ident")

#integrate just L5 ET
Idents(data_et) <- data_et$orig.ident
human <- subset(data_et, idents = "human")
macaque <- subset(data_et, idents = "macaque")
marmoset <- subset(data_et, idents = "marmoset")
mouse <- subset(data_et, idents = "mouse")

human_et <- CreateSeuratObject(counts = human@assays$RNA@counts, meta.data = human@meta.data)
macaque_et <- CreateSeuratObject(counts = macaque@assays$RNA@counts, meta.data = macaque@meta.data)
marmoset_et <- CreateSeuratObject(counts = marmoset@assays$RNA@counts, meta.data = marmoset@meta.data)
mouse_et <- CreateSeuratObject(counts = mouse@assays$RNA@counts, meta.data = mouse@meta.data)

all.data <- merge(x = human_et, y = list(macaque_et, marmoset_et, mouse_et), add.cell.ids = c("human","macaque", "marmoset", "mouse"))
combined.list <- SplitObject(all.data, split.by = "orig.ident")

for (i in 1:length(combined.list)) {
  combined.list[[i]] <- SCTransform(combined.list[[i]], verbose = FALSE)
}

combined.features <- SelectIntegrationFeatures(object.list = combined.list, nfeatures = 3000)
combined.list <- PrepSCTIntegration(object.list = combined.list, anchor.features = combined.features, 
                                    verbose = T)

combined.anchors <- FindIntegrationAnchors(object.list = combined.list, normalization.method = "SCT", 
                                           anchor.features = combined.features, verbose = T)
combined.integrated <- IntegrateData(anchorset = combined.anchors, normalization.method = "SCT", 
                                     verbose = T)

combined.integrated <- RunPCA(combined.integrated, verbose = T)
combined.integrated <- RunUMAP(combined.integrated, dims = 1:30)
combined.integrated <- FindNeighbors(combined.integrated, dims = 1:30)
combined.integrated <- FindClusters(combined.integrated, resolution = 1)

#Visualize
DimPlot(combined.integrated, reduction = "umap", split.by = "orig.ident", group.by = "cluster_label") + NoLegend()
DimPlot(combined.integrated, reduction = "umap", split.by = "orig.ident", group.by = "seurat_clusters")
FeaturePlot(combined.integrated, reduction = "umap", split.by = "orig.ident", features = "nCount_RNA")
