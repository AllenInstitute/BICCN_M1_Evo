library(mgcv)
library(Seurat)
library(cowplot)
library(feather)
library(dplyr)
library(Matrix)
library(Hmisc)
library(tidyr)
library(gplots)
library(gdata)
library(stringr)
library(rlang)
library(scrattch.hicat)
library(scrattch.io)
library(matrixStats)
library(ggplot2)
library(igraph)
library(visNetwork)
library(riverplot)
library(scales)
library(dendextend)
library(reshape)
library(openxlsx)
library(rhdf5)
options(stringsAsFactors = FALSE)

#load integrated excitatory data
sample.combined <- readRDS("~/NeMO_analysis_folder/Transcriptomics/cross_species_integration/sample.combined_exc_4_species_integration.RDS")

#subset L5 ET and L5 IT cells
Idents(sample.combined) <- sample.combined$subclass_label
subclasses <- c("L5 ET", "L5 IT")
sample.combined <- subset(sample.combined, idents = subclasses)

#subset by species
Idents(sample.combined) <- sample.combined$orig.ident
human <- subset(sample.combined, idents = "human")
macaque <- subset(sample.combined, idents = "macaque")
marmoset <- subset(sample.combined, idents = "marmoset")
mouse <- subset(sample.combined, idents = "mouse")

#Find L5 ET specific genes for each species
Idents(human) <- human$subclass_label
Idents(macaque) <- macaque$subclass_label
Idents(marmoset) <- marmoset$subclass_label
Idents(mouse) <- mouse$subclass_label

#Find Log2FC L5 ET > L5 IT
tmp <- human@assays$SCT@counts
tmp@x <- log2(tmp@x + 1)
L5ET_human_avg <- tmp[ , which(colnames(tmp) %in% names(which(human$subclass_label == "L5 ET")))]
L5IT_human_avg <- tmp[ , which(colnames(tmp) %in% names(which(human$subclass_label == "L5 IT")))]
L5ET_human_avg <- rowMeans(L5ET_human_avg)
L5IT_human_avg <- rowMeans(L5IT_human_avg)
Log2FC_human <- L5ET_human_avg - L5IT_human_avg

tmp <- macaque@assays$SCT@counts
tmp@x <- log2(tmp@x + 1)
L5ET_macaque_avg <- tmp[ , which(colnames(tmp) %in% names(which(macaque$subclass_label == "L5 ET")))]
L5IT_macaque_avg <- tmp[ , which(colnames(tmp) %in% names(which(macaque$subclass_label == "L5 IT")))]
L5ET_macaque_avg <- rowMeans(L5ET_macaque_avg)
L5IT_macaque_avg <- rowMeans(L5IT_macaque_avg)
Log2FC_macaque <- L5ET_macaque_avg - L5IT_macaque_avg

tmp <- marmoset@assays$SCT@counts
tmp@x <- log2(tmp@x + 1)
L5ET_marmoset_avg <- tmp[ , which(colnames(tmp) %in% names(which(marmoset$subclass_label == "L5 ET")))]
L5IT_marmoset_avg <- tmp[ , which(colnames(tmp) %in% names(which(marmoset$subclass_label == "L5 IT")))]
L5ET_marmoset_avg <- rowMeans(L5ET_marmoset_avg)
L5IT_marmoset_avg <- rowMeans(L5IT_marmoset_avg)
Log2FC_marmoset <- L5ET_marmoset_avg - L5IT_marmoset_avg

tmp <- mouse@assays$SCT@counts
tmp@x <- log2(tmp@x + 1)
L5ET_mouse_avg <- tmp[ , which(colnames(tmp) %in% names(which(mouse$subclass_label == "L5 ET")))]
L5IT_mouse_avg <- tmp[ , which(colnames(tmp) %in% names(which(mouse$subclass_label == "L5 IT")))]
L5ET_mouse_avg <- rowMeans(L5ET_mouse_avg)
L5IT_mouse_avg <- rowMeans(L5IT_mouse_avg)
Log2FC_mouse <- L5ET_mouse_avg - L5IT_mouse_avg

#Find Genes that decrease with evolutionary distance from humans
All_species <- data.frame(human = Log2FC_human)
All_species$macaque <- Log2FC_macaque[match(rownames(All_species), names(Log2FC_macaque))]
All_species$marmoset <- Log2FC_marmoset[match(rownames(All_species), names(Log2FC_marmoset))]
All_species$mouse <- Log2FC_mouse[match(rownames(All_species), names(Log2FC_mouse))]

All_species$trend <- FALSE
for(i in 1:nrow(All_species)){
  if(All_species$human[i] > All_species$macaque[i]){
    if(All_species$macaque[i] > All_species$marmoset[i]){
      if(All_species$marmoset[i] > All_species$mouse[i]){
        All_species$trend[i] <- TRUE
        
      }
    }
  }
}
All_species <- All_species[which(All_species$trend == TRUE), ]
All_species <- All_species[which(All_species$human > 0.5), ]



sample.combined$heatmap_order <- paste(sample.combined$orig.ident, sample.combined$subclass_label, sep = "_")
Idents(sample.combined) <- sample.combined$heatmap_order
library(viridis)
heatmap_cells <- subset(sample.combined, downsample = 200)

heatmap_data <- data.frame(human = Log2FC_human)
heatmap_data$macaque = Log2FC_macaque
heatmap_data$marmoset = Log2FC_marmoset
heatmap_data$mouse = Log2FC_mouse
heatmap_data <- heatmap_data[which(rownames(heatmap_data) %in% rownames(All_species)), ]
heatmap_data <- as.matrix(heatmap_data)
heatmap.2(heatmap_data, scale = "none", col = viridis(100), Colv = NA, trace = "none", cexRow = .5) 

##########################################################################################################################
############################ spaghetti                ######################################################################
##########################################################################################################################
spaghetti_dat <- data.frame(genes = rep(rownames(All_species), 4)) 
spaghetti_dat$Log2FC <- 0 
spaghetti_dat$species <- NA 
spaghetti_dat$Log2FC[1:131] <- All_species$human
spaghetti_dat$species[1:131] <- "Human"
spaghetti_dat$Log2FC[132:262] <- All_species$macaque
spaghetti_dat$species[132:262] <- "Macaque"
spaghetti_dat$Log2FC[263:393] <- All_species$marmoset
spaghetti_dat$species[263:393] <- "Marmoset"
spaghetti_dat$Log2FC[394:524] <- All_species$mouse
spaghetti_dat$species[394:524] <- "Mouse"

#color genes by GO category
axon_genes <- c("ALCAM", "SLIT3", "SLIT2", "SEMA3D", "IRS2", "PTPRM", "ROBO2", "VSTM2L", "PRKCA", "EPHA6", "GFRA1", "NLGN1", "NYAP2", "TMEM108")
cell_adhesion_genes <- c("PCDH19", "ITGA7", "TENM1", "ALCAM", "PTPRM", "IL1RAP", "ROBO2", "FNDC3A", "CYP1B1", "VSTM2L", "NLGN1", "LMO7", "IL1RAPL1", "EPDR1", "SPON2", "DDR2", "ANGPT1", "PRKCA")
luteolysis_genes <- c("ROBO2", "SLIT3", "SLIT2")
all_gene_ontology <- data.frame(genes = unique(c(axon_genes, cell_adhesion_genes, luteolysis_genes)))
all_gene_ontology$axon <- FALSE
all_gene_ontology$adhesion <- FALSE
all_gene_ontology$luteolysis <- FALSE
all_gene_ontology$color <- NA
for(i in 1:nrow(all_gene_ontology)){
  if(all_gene_ontology$genes[i] %in% axon_genes){
    all_gene_ontology$axon[i] <- TRUE
  }  
  if(all_gene_ontology$genes[i] %in% cell_adhesion_genes){
    all_gene_ontology$adhesion[i] <- TRUE
  }  
  if(all_gene_ontology$genes[i] %in% luteolysis_genes){
    all_gene_ontology$luteolysis[i] <- TRUE
  }  
  if(all_gene_ontology$axon[i] == TRUE & all_gene_ontology$adhesion[i] == TRUE & all_gene_ontology$luteolysis[i] == TRUE){
    all_gene_ontology$color[i] <- "magenta"
  }
  if(all_gene_ontology$axon[i] == TRUE & all_gene_ontology$adhesion[i] == FALSE & all_gene_ontology$luteolysis[i] == TRUE){
    all_gene_ontology$color[i] <- "magenta4"
  }
  if(all_gene_ontology$axon[i] == FALSE & all_gene_ontology$adhesion[i] == TRUE & all_gene_ontology$luteolysis[i] == FALSE){
    all_gene_ontology$color[i] <- "green4"
  }
  if(all_gene_ontology$axon[i] == TRUE & all_gene_ontology$adhesion[i] == FALSE & all_gene_ontology$luteolysis[i] == FALSE){
    all_gene_ontology$color[i] <- "#2b569a"
  }
  if(all_gene_ontology$axon[i] == TRUE & all_gene_ontology$adhesion[i] == TRUE & all_gene_ontology$luteolysis[i] == FALSE){
    all_gene_ontology$color[i] <- "cyan"
  }
}

spaghetti_dat$color <- all_gene_ontology$color[match(spaghetti_dat$genes, all_gene_ontology$genes)]
spaghetti_dat$color[which(is.na(spaghetti_dat$color))] <- "grey74"

spaghetti_dat$label <- NA
spaghetti_dat$label[which(spaghetti_dat$genes %in% all_gene_ontology$genes)] <- spaghetti_dat$genes[which(spaghetti_dat$genes %in% all_gene_ontology$genes)]
spaghetti_dat$label[1:393] <- NA

size_line <- 1.5
library(ggrepel)
ggplot() +
  geom_line(data = spaghetti_dat[which(spaghetti_dat$color == "grey74"), ], aes(x = species, y = Log2FC, group = genes), color = "grey74", size = 1) +
  geom_line(data = spaghetti_dat[which(spaghetti_dat$color == "green4"), ], aes(x = species, y = Log2FC, group = genes), color = "green4", size = size_line) +
  geom_line(data = spaghetti_dat[which(spaghetti_dat$color == "magenta"), ], aes(x = species, y = Log2FC, group = genes), color = "magenta", size = size_line) +
  geom_line(data = spaghetti_dat[which(spaghetti_dat$color == "magenta4"), ], aes(x = species, y = Log2FC, group = genes), color = "magenta4", size = size_line) +
  geom_line(data = spaghetti_dat[which(spaghetti_dat$color == "#2b569a"), ], aes(x = species, y = Log2FC, group = genes), color = "#2b569a", size = size_line) +
  geom_line(data = spaghetti_dat[which(spaghetti_dat$color == "cyan"), ], aes(x = species, y = Log2FC, group = genes), color = "cyan", size = size_line) +
  geom_text_repel(data = spaghetti_dat, aes(x = species, y = Log2FC, label = label),
                  size = 2,
                  nudge_x = 40,
                  na.rm = TRUE,
                  angle = 45)+
  theme(axis.text.y  = element_text(size=20))+
  theme(axis.title.y = element_text(size=20))+
  theme(axis.text.x  = element_text(size=20))+
  theme(axis.title.x = element_text(size = 0)) 

#just axon guidance call outs
spaghetti_dat <- data.frame(genes = rep(rownames(All_species), 4)) 
spaghetti_dat$Log2FC <- 0 
spaghetti_dat$species <- NA 
spaghetti_dat$Log2FC[1:131] <- All_species$human
spaghetti_dat$species[1:131] <- "Human"
spaghetti_dat$Log2FC[132:262] <- All_species$macaque
spaghetti_dat$species[132:262] <- "Macaque"
spaghetti_dat$Log2FC[263:393] <- All_species$marmoset
spaghetti_dat$species[263:393] <- "Marmoset"
spaghetti_dat$Log2FC[394:524] <- All_species$mouse
spaghetti_dat$species[394:524] <- "Mouse"

#color genes by GO category
axon_genes <- c("ALCAM", "SLIT3", "SLIT2", "SEMA3D", "IRS2", "PTPRM", "ROBO2", "VSTM2L", "PRKCA", "EPHA6", "GFRA1", "NLGN1", "NYAP2", "TMEM108")

spaghetti_dat$color <- NA
spaghetti_dat$color[which(spaghetti_dat$genes %in% axon_genes)] <- "#2b569a"
spaghetti_dat$color[which(is.na(spaghetti_dat$color))] <- "grey74"

spaghetti_dat$label <- NA
spaghetti_dat$label[which(spaghetti_dat$genes %in% axon_genes)] <- spaghetti_dat$genes[which(spaghetti_dat$genes %in% axon_genes)]
spaghetti_dat$label[1:393] <- NA

spaghetti_dat$type <- "solid"
dashed <- c("ROBO2", "SLIT3", "SLIT2")
spaghetti_dat$type[which(spaghetti_dat$genes %in% dashed)] <- "dotted"
spaghetti_dat$plot_order <- paste(spaghetti_dat$color, spaghetti_dat$type, sep = "_")

size_line <- 1.5
library(ggrepel)
ggplot() +
  geom_line(data = spaghetti_dat[which(spaghetti_dat$color == "grey74"), ], aes(x = species, y = Log2FC, group = genes), color = "grey74", size = 1, linetype = "solid") +
  geom_line(data = spaghetti_dat[which(spaghetti_dat$plot_order == "#2b569a_solid"), ], aes(x = species, y = Log2FC, group = genes), color = "#2b569a", size = size_line, linetype = "solid") +
  geom_line(data = spaghetti_dat[which(spaghetti_dat$plot_order == "#2b569a_dotted"), ], aes(x = species, y = Log2FC, group = genes), color = "#2b569a", size = size_line, linetype = "dashed") +
  geom_text_repel(data = spaghetti_dat, aes(x = species, y = Log2FC, label = label),
                  size = 4,
                  nudge_x = 4,
                  na.rm = TRUE,
                  angle = 45)+
  theme(axis.text.y  = element_text(size=20))+
  theme(axis.title.y = element_text(size=20))+
  theme(axis.text.x  = element_text(size=20))+
  theme(axis.title.x = element_text(size = 0)) 


