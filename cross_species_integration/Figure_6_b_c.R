library(Seurat)
library(gridExtra)

sample.combined <- readRDS("http://data.nemoarchive.org/biccn/lab/lein/2020_M1_study_analysis/Transcriptomics/cross_species_integration/sample.combined_glia_integration.RDS")

#label Chcs and Basket cells in metadata
Idents(sample.combined) <- sample.combined$subclass_label
pvalb_data <- subset(sample.combined, idents = "Pvalb")

chandelier_clusters <- c("Inh L1-6 PVALB COL15A1", "Pvalb Vipr2_1", "Pvalb Vipr2_2", "Inh PVALB FAM19A4")
pvalb_data$comparison <- "Basket"
pvalb_data$comparison[which(pvalb_data$cluster_label %in% chandelier_clusters)] <- "Chandelier"

#Find DE genes for Chcs in each species
Idents(pvalb_data) <- pvalb_data$orig.ident
human <- subset(pvalb_data, idents = "human")
marmoset <- subset(pvalb_data, idents = "marmoset")
mouse <- subset(pvalb_data, idents = "mouse")

Idents(human) <- human$comparison
Idents(marmoset) <- marmoset$comparison
Idents(mouse) <- mouse$comparison

human_genes <- FindAllMarkers(human, assay = "SCT", slot = "data", test.use = "roc")
marmoset_genes <- FindAllMarkers(marmoset, assay = "SCT", slot = "data", test.use = "roc")
mouse_genes <- FindAllMarkers(mouse, assay = "SCT", slot = "data", test.use = "roc")


############################################################################################################################
########################## Venn diagram               #####################################################################
############################################################################################################################

human_genes <- human_genes[which(human_genes$avg_diff > 0), ]
marmoset_genes <- marmoset_genes[which(marmoset_genes$avg_diff > 0), ]
mouse_genes <- mouse_genes[which(mouse_genes$avg_diff > 0), ]

human_genes <- human_genes[which(human_genes$cluster == "Chandelier"), ]
marmoset_genes <- marmoset_genes[which(marmoset_genes$cluster == "Chandelier"), ]
mouse_genes <- mouse_genes[which(mouse_genes$cluster == "Chandelier"), ]

all.genes <- data.frame(genes = unique(c(human_genes$gene, marmoset_genes$gene, mouse_genes$gene)))
all.genes$Human <- as.character(match(all.genes$genes, human_genes$gene))
all.genes$Marmoset <- as.character(match(all.genes$genes, marmoset_genes$gene))
all.genes$Mouse <- as.character(match(all.genes$genes, mouse_genes$gene))
all.genes$Human[which(is.na(all.genes$Human))] <- FALSE 
all.genes$Human[which(all.genes$Human != FALSE)] <- TRUE
all.genes$Marmoset[which(is.na(all.genes$Marmoset))] <- FALSE 
all.genes$Marmoset[which(all.genes$Marmoset != FALSE)] <- TRUE
all.genes$Mouse[which(is.na(all.genes$Mouse))] <- FALSE 
all.genes$Mouse[which(all.genes$Mouse != FALSE)] <- TRUE
all.genes$Human <- as.logical(all.genes$Human)
all.genes$Marmoset <- as.logical(all.genes$Marmoset)
all.genes$Mouse <- as.logical(all.genes$Mouse)

library(eulerr)
plot(euler(
  all.genes[ ,2:4]),
  quantities = list(cex = 3),
  labels = NULL,
  main = paste0("Chandelier vs. Basket"),
  fills = c("royalblue1", "maroon4", "sienna2")
)

############################################################################################################################
########################## scatter plots               #####################################################################
############################################################################################################################

#use this code to change the default log expression to log2 expression
human@assays$SCT@data <- human@assays$SCT@counts
human@assays$SCT@data@x <- log2(human@assays$SCT@data@x + 1)
marmoset@assays$SCT@data <- marmoset@assays$SCT@counts
marmoset@assays$SCT@data@x <- log2(marmoset@assays$SCT@data@x + 1)
mouse@assays$SCT@data <- mouse@assays$SCT@counts
mouse@assays$SCT@data@x <- log2(mouse@assays$SCT@data@x + 1)

human_mat <- human@assays$SCT@data
colnames(human_mat) <- human$comparison
marmoset_mat <- marmoset@assays$SCT@data
colnames(marmoset_mat) <- marmoset$comparison
mouse_mat <- mouse@assays$SCT@data
colnames(mouse_mat) <- mouse$comparison

#find avg expression for Chc and Basket (log(umi + 1))
human_avg <- data.frame(Chandelier = rowMeans(human_mat[ , which(colnames(human_mat) == "Chandelier")]))
human_avg$Basket <- rowMeans(human_mat[ , which(colnames(human_mat) == "Basket")])

marmoset_avg <- data.frame(Chandelier = rowMeans(marmoset_mat[ , which(colnames(marmoset_mat) == "Chandelier")]))
marmoset_avg$Basket <- rowMeans(marmoset_mat[ , which(colnames(marmoset_mat) == "Basket")])

mouse_avg <- data.frame(Chandelier = rowMeans(mouse_mat[ , which(colnames(mouse_mat) == "Chandelier")]))
mouse_avg$Basket <- rowMeans(mouse_mat[ , which(colnames(mouse_mat) == "Basket")])

#plot scatter
conserved_genes <- all.genes$genes[which(all.genes$Human == T & all.genes$Marmoset ==T & all.genes$Mouse == T)]

human_avg$conserved <- F
human_avg$conserved[which(rownames(human_avg) %in% conserved_genes)] <- T

marmoset_avg$conserved <- F
marmoset_avg$conserved[which(rownames(marmoset_avg) %in% conserved_genes)] <- T

mouse_avg$conserved <- F
mouse_avg$conserved[which(rownames(mouse_avg) %in% conserved_genes)] <- T

TFs_of_interest <- c("RORA", "TRPS1", "NFIB")
human_avg$label <- NA
human_avg$label[which(rownames(human_avg) %in% TFs_of_interest)] <- rownames(human_avg)[which(rownames(human_avg) %in% TFs_of_interest)]
marmoset_avg$label <- NA
marmoset_avg$label[which(rownames(marmoset_avg) %in% TFs_of_interest)] <- rownames(marmoset_avg)[which(rownames(marmoset_avg) %in% TFs_of_interest)]
mouse_avg$label <- NA
mouse_avg$label[which(rownames(mouse_avg) %in% TFs_of_interest)] <- rownames(mouse_avg)[which(rownames(mouse_avg) %in% TFs_of_interest)]

library(ggrepel)
library(readxl)

tf_genelist <- "load TF gene list"
human_avg$tf <- F
human_avg$tf[which(rownames(human_avg) %in% tf_genelist$X__1)] <- T
marmoset_avg$tf <- F
marmoset_avg$tf[which(rownames(marmoset_avg) %in% tf_genelist$X__1)] <- T
mouse_avg$tf <- F
mouse_avg$tf[which(rownames(mouse_avg) %in% tf_genelist$X__1)] <- T

colnames(human_avg) <- paste("human", colnames(human_avg), sep = "_")
colnames(marmoset_avg) <- paste("marmoset", colnames(marmoset_avg), sep = "_")
colnames(mouse_avg) <- paste("mouse", colnames(mouse_avg), sep = "_")
to_plot <- cbind(human_avg, marmoset_avg, mouse_avg)

#plot  ChC vs. ChC
p1 <- ggplot() +
  geom_point(data = to_plot[which(to_plot$human_conserved == F & to_plot$human_tf == F), ], mapping = aes(x = human_Chandelier, y = marmoset_Chandelier), color = "grey80", size = 0.5) +
  geom_point(data = to_plot[which(to_plot$human_conserved == F & to_plot$human_tf == T), ], mapping = aes(x = human_Chandelier, y = marmoset_Chandelier), color = "cyan2", size = 1, shape = 15) +
  geom_point(data = to_plot[which(to_plot$human_conserved == T & to_plot$human_tf == F), ], mapping = aes(x = human_Chandelier, y = marmoset_Chandelier), color = "red", size = 1) +
  geom_point(data = to_plot[which(to_plot$human_conserved == T & is.na(to_plot$human_label) == F), ], mapping = aes(x = human_Chandelier, y = marmoset_Chandelier), color = "magenta", size = 3, shape = 15) +
  xlab("human_Chandelier (log2 expression)") +
  ylab("marmoset_Chandelier (log2 expression)") +
  xlim(c(0, 8)) +
  ylim(c(0, 8)) +
  theme(aspect.ratio = 1)
p1 <- p1 + geom_text_repel(data = to_plot[which(to_plot$human_conserved == T & is.na(to_plot$human_label) == F), ], aes(x = human_Chandelier, y = marmoset_Chandelier, label = human_label),
                           nudge_x = 5,
                           nudge_y = 0.5,
                           size = 3)

p2 <- ggplot() +
  geom_point(data = to_plot[which(to_plot$human_conserved == F & to_plot$human_tf == F), ], mapping = aes(x = human_Chandelier, y = mouse_Chandelier), color = "grey80", size = 0.5) +
  geom_point(data = to_plot[which(to_plot$human_conserved == F & to_plot$human_tf == T), ], mapping = aes(x = human_Chandelier, y = mouse_Chandelier), color = "cyan2", size = 1, shape = 15) +
  geom_point(data = to_plot[which(to_plot$human_conserved == T & to_plot$human_tf == F), ], mapping = aes(x = human_Chandelier, y = mouse_Chandelier), color = "red", size = 1) +
  geom_point(data = to_plot[which(to_plot$human_conserved == T & is.na(to_plot$human_label) == F), ], mapping = aes(x = human_Chandelier, y = mouse_Chandelier), color = "magenta", size = 3, shape = 15) +
  xlab("human_Chandelier (log2 expression)") +
  ylab("mouse_Chandelier (log2 expression)") +
  xlim(c(0, 8)) +
  ylim(c(0, 8)) +
  theme(aspect.ratio = 1)
p2 <- p2 + geom_text_repel(data = to_plot[which(to_plot$human_conserved == T & is.na(to_plot$human_label) == F), ], aes(x = human_Chandelier, y = mouse_Chandelier, label = human_label),
                           nudge_x = 5,
                           nudge_y = 0.5,
                           size = 3)

p3 <- ggplot() +
  geom_point(data = to_plot[which(to_plot$marmoset_conserved == F & to_plot$marmoset_tf == F), ], mapping = aes(x = marmoset_Chandelier, y = mouse_Chandelier), color = "grey80", size = 0.5) +
  geom_point(data = to_plot[which(to_plot$marmoset_conserved == F & to_plot$marmoset_tf == T), ], mapping = aes(x = marmoset_Chandelier, y = mouse_Chandelier), color = "cyan2", size = 1, shape = 15) +
  geom_point(data = to_plot[which(to_plot$marmoset_conserved == T & to_plot$marmoset_tf == F), ], mapping = aes(x = marmoset_Chandelier, y = mouse_Chandelier), color = "red", size = 1) +
  geom_point(data = to_plot[which(to_plot$marmoset_conserved == T & is.na(to_plot$marmoset_label) == F), ], mapping = aes(x = marmoset_Chandelier, y = mouse_Chandelier), color = "magenta", size = 3, shape = 15) +
  xlab("marmoset_Chandelier (log2 expression)") +
  ylab("mouse_Chandelier (log2 expression)") +
  xlim(c(0, 8)) +
  ylim(c(0, 8)) +
  theme(aspect.ratio = 1)
p3 <- p3 + geom_text_repel(data = to_plot[which(to_plot$marmoset_conserved == T & is.na(to_plot$marmoset_label) == F), ], aes(x = marmoset_Chandelier, y = mouse_Chandelier, label = marmoset_label),
                           nudge_x = 5,
                           nudge_y = 0.5,
                           size = 3)

grid.arrange(p1, p2, p3, nrow = 1)

###### plot ChC vs BC
p1 <- ggplot() +
  geom_point(data = to_plot[which(to_plot$human_conserved == F & to_plot$human_tf == F), ], mapping = aes(x = human_Chandelier, y = human_Basket), color = "grey80", size = 0.5) +
  geom_point(data = to_plot[which(to_plot$human_conserved == F & to_plot$human_tf == T), ], mapping = aes(x = human_Chandelier, y = human_Basket), color = "cyan2", size = 1, shape = 15) +
  geom_point(data = to_plot[which(to_plot$human_conserved == T & to_plot$human_tf == F), ], mapping = aes(x = human_Chandelier, y = human_Basket), color = "red", size = 1) +
  geom_point(data = to_plot[which(to_plot$human_conserved == T & is.na(to_plot$human_label) == F), ], mapping = aes(x = human_Chandelier, y = human_Basket), color = "magenta", size = 3, shape = 15) +
  xlab("human_Chandelier (log2 expression)") +
  ylab("human_Basket (log2 expression)") +
  xlim(c(0, 8)) +
  ylim(c(0, 8)) +
  theme(aspect.ratio = 1)
p1 <- p1 + geom_text_repel(data = to_plot[which(to_plot$human_conserved == T & is.na(to_plot$human_label) == F), ], aes(x = human_Chandelier, y = human_Basket, label = human_label),
                           nudge_x = 5,
                           nudge_y = 0.5,
                           size = 3)

p2 <- ggplot() +
  geom_point(data = to_plot[which(to_plot$marmoset_conserved == F & to_plot$marmoset_tf == F), ], mapping = aes(x = marmoset_Chandelier, y = marmoset_Basket), color = "grey80", size = 0.5) +
  geom_point(data = to_plot[which(to_plot$marmoset_conserved == F & to_plot$marmoset_tf == T), ], mapping = aes(x = marmoset_Chandelier, y = marmoset_Basket), color = "cyan2", size = 1, shape = 15) +
  geom_point(data = to_plot[which(to_plot$marmoset_conserved == T & to_plot$marmoset_tf == F), ], mapping = aes(x = marmoset_Chandelier, y = marmoset_Basket), color = "red", size = 1) +
  geom_point(data = to_plot[which(to_plot$marmoset_conserved == T & is.na(to_plot$marmoset_label) == F), ], mapping = aes(x = marmoset_Chandelier, y = marmoset_Basket), color = "magenta", size = 3, shape = 15) +
  xlab("marmoset_Chandelier (log2 expression)") +
  ylab("marmoset_Basket (log2 expression)") +
  xlim(c(0, 8)) +
  ylim(c(0, 8)) +
  theme(aspect.ratio = 1)
p2 <- p2 + geom_text_repel(data = to_plot[which(to_plot$marmoset_conserved == T & is.na(to_plot$marmoset_label) == F), ], aes(x = marmoset_Chandelier, y = marmoset_Basket, label = marmoset_label),
                           nudge_x = 5,
                           nudge_y = 0.5,
                           size = 3)

p3 <- ggplot() +
  geom_point(data = to_plot[which(to_plot$mouse_conserved == F & to_plot$mouse_tf == F), ], mapping = aes(x = mouse_Chandelier, y = mouse_Basket), color = "grey80", size = 0.5) +
  geom_point(data = to_plot[which(to_plot$mouse_conserved == F & to_plot$mouse_tf == T), ], mapping = aes(x = mouse_Chandelier, y = mouse_Basket), color = "cyan2", size = 1, shape = 15) +
  geom_point(data = to_plot[which(to_plot$mouse_conserved == T & to_plot$mouse_tf == F), ], mapping = aes(x = mouse_Chandelier, y = mouse_Basket), color = "red", size = 1) +
  geom_point(data = to_plot[which(to_plot$mouse_conserved == T & is.na(to_plot$mouse_label) == F), ], mapping = aes(x = mouse_Chandelier, y = mouse_Basket), color = "magenta", size = 3, shape = 15) +
  xlab("mouse_Chandelier (log2 expression)") +
  ylab("mouse_Basket (log2 expression)") +
  xlim(c(0, 8)) +
  ylim(c(0, 8)) +
  theme(aspect.ratio = 1)
p3 <- p3 + geom_text_repel(data = to_plot[which(to_plot$mouse_conserved == T & is.na(to_plot$mouse_label) == F), ], aes(x = mouse_Chandelier, y = mouse_Basket, label = mouse_label),
                           nudge_x = 5,
                           nudge_y = 0.5,
                           size = 3)

grid.arrange(p1, p2, p3, nrow = 1)

###############################################################
############### other point to address ########################
###############################################################

#how many and what % of TFs in each species have expression in ChC at least as high as NFIB?
NFIB_thresh_genes <- data.frame(human_n_genes = 0, human_pct_genes = 0, 
                                marmoset_n_genes = 0, marmoset_pct_genes = 0,
                                mouse_n_genes = 0, mouse_pct_genes = 0)

NFIB_thresh_genes$human_n_genes <-  nrow(to_plot[which(to_plot$human_tf & to_plot$human_Chandelier > to_plot["NFIB", "human_Chandelier"]), ])
NFIB_thresh_genes$human_pct_genes <- NFIB_thresh_genes$human_n_genes / length(which(to_plot$human_tf))
NFIB_thresh_genes$marmoset_n_genes <-  nrow(to_plot[which(to_plot$marmoset_tf & to_plot$marmoset_Chandelier > to_plot["NFIB", "marmoset_Chandelier"]), ])
NFIB_thresh_genes$marmoset_pct_genes <- NFIB_thresh_genes$marmoset_n_genes / length(which(to_plot$marmoset_tf))
NFIB_thresh_genes$mouse_n_genes <-  nrow(to_plot[which(to_plot$mouse_tf & to_plot$mouse_Chandelier > to_plot["NFIB", "mouse_Chandelier"]), ])
NFIB_thresh_genes$mouse_pct_genes <- NFIB_thresh_genes$mouse_n_genes / length(which(to_plot$mouse_tf))
NFIB_thresh_genes

#Color scatter to include human-specific TFs genes
human_specific <- as.character(all.genes$genes[which(all.genes$Human & all.genes$Marmoset == F & all.genes$Mouse == F)])

to_plot$human_specific <- F
to_plot$human_specific[which(rownames(to_plot) %in% human_specific & to_plot$human_tf)] <- T

#plot  ChC vs. ChC
p1 <- ggplot() +
  geom_point(data = to_plot[which(to_plot$human_conserved == F & to_plot$human_tf == F), ], mapping = aes(x = human_Chandelier, y = marmoset_Chandelier), color = "grey80", size = 0.5) +
  geom_point(data = to_plot[which(to_plot$human_conserved == F & to_plot$human_tf == T), ], mapping = aes(x = human_Chandelier, y = marmoset_Chandelier), color = "cyan2", size = 1, shape = 15) +
  geom_point(data = to_plot[which(to_plot$human_conserved == T & to_plot$human_tf == F), ], mapping = aes(x = human_Chandelier, y = marmoset_Chandelier), color = "red", size = 1) +
  geom_point(data = to_plot[which(to_plot$human_specific), ], mapping = aes(x = human_Chandelier, y = marmoset_Chandelier), color = "green4", size = 2, shape = 15) +
  geom_point(data = to_plot[which(to_plot$human_conserved == T & is.na(to_plot$human_label) == F), ], mapping = aes(x = human_Chandelier, y = marmoset_Chandelier), color = "magenta", size = 3, shape = 15) +
  xlab("human_Chandelier (log2 expression)") +
  ylab("marmoset_Chandelier (log2 expression)") +
  xlim(c(0, 8)) +
  ylim(c(0, 8)) 
p1 <- p1 + geom_text_repel(data = to_plot[which(to_plot$human_conserved == T & is.na(to_plot$human_label) == F), ], aes(x = human_Chandelier, y = marmoset_Chandelier, label = human_label),
                           nudge_x = 5,
                           nudge_y = 0.5,
                           size = 3)
p1 <- p1 +  geom_abline(intercept = 0)
p1 <- p1 + theme_bw() + theme(panel.border = element_blank(), aspect.ratio = 1, panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

p2 <- ggplot() +
  geom_point(data = to_plot[which(to_plot$human_conserved == F & to_plot$human_tf == F), ], mapping = aes(x = human_Chandelier, y = mouse_Chandelier), color = "grey80", size = 0.5) +
  geom_point(data = to_plot[which(to_plot$human_conserved == F & to_plot$human_tf == T), ], mapping = aes(x = human_Chandelier, y = mouse_Chandelier), color = "cyan2", size = 1, shape = 15) +
  geom_point(data = to_plot[which(to_plot$human_conserved == T & to_plot$human_tf == F), ], mapping = aes(x = human_Chandelier, y = mouse_Chandelier), color = "red", size = 1) +
  geom_point(data = to_plot[which(to_plot$human_specific), ], mapping = aes(x = human_Chandelier, y = mouse_Chandelier), color = "green4", size = 2, shape = 15) +
  geom_point(data = to_plot[which(to_plot$human_conserved == T & is.na(to_plot$human_label) == F), ], mapping = aes(x = human_Chandelier, y = mouse_Chandelier), color = "magenta", size = 3, shape = 15) +
  xlab("human_Chandelier (log2 expression)") +
  ylab("mouse_Chandelier (log2 expression)") +
  xlim(c(0, 8)) +
  ylim(c(0, 8)) 
p2 <- p2 + geom_text_repel(data = to_plot[which(to_plot$human_conserved == T & is.na(to_plot$human_label) == F), ], aes(x = human_Chandelier, y = mouse_Chandelier, label = human_label),
                           nudge_x = 5,
                           nudge_y = 0.5,
                           size = 3)
p2 <- p2 +  geom_abline(intercept = 0)
p2 <- p2 + theme_bw() + theme(panel.border = element_blank(), aspect.ratio = 1, panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

p3 <- ggplot() +
  geom_point(data = to_plot[which(to_plot$marmoset_conserved == F & to_plot$marmoset_tf == F), ], mapping = aes(x = marmoset_Chandelier, y = mouse_Chandelier), color = "grey80", size = 0.5) +
  geom_point(data = to_plot[which(to_plot$marmoset_conserved == F & to_plot$marmoset_tf == T), ], mapping = aes(x = marmoset_Chandelier, y = mouse_Chandelier), color = "cyan2", size = 1, shape = 15) +
  geom_point(data = to_plot[which(to_plot$marmoset_conserved == T & to_plot$marmoset_tf == F), ], mapping = aes(x = marmoset_Chandelier, y = mouse_Chandelier), color = "red", size = 1) +
  geom_point(data = to_plot[which(to_plot$human_specific), ], mapping = aes(x = marmoset_Chandelier, y = mouse_Chandelier), color = "green4", size = 2, shape = 15) +
  geom_point(data = to_plot[which(to_plot$marmoset_conserved == T & is.na(to_plot$marmoset_label) == F), ], mapping = aes(x = marmoset_Chandelier, y = mouse_Chandelier), color = "magenta", size = 3, shape = 15) +
  xlab("marmoset_Chandelier (log2 expression)") +
  ylab("mouse_Chandelier (log2 expression)") +
  xlim(c(0, 8)) +
  ylim(c(0, 8)) 
p3 <- p3 + geom_text_repel(data = to_plot[which(to_plot$marmoset_conserved == T & is.na(to_plot$marmoset_label) == F), ], aes(x = marmoset_Chandelier, y = mouse_Chandelier, label = marmoset_label),
                           nudge_x = 5,
                           nudge_y = 0.5,
                           size = 3)
p3 <- p3 +  geom_abline(intercept = 0)
p3 <- p3 + theme_bw() + theme(panel.border = element_blank(), aspect.ratio = 1, panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

grid.arrange(p1, p2, p3, nrow = 1)


#plot  ChC vs. BC
p1 <- ggplot() +
  geom_point(data = to_plot[which(to_plot$human_conserved == F & to_plot$human_tf == F), ], mapping = aes(x = human_Chandelier, y = human_Basket), color = "grey80", size = 0.5) +
  geom_point(data = to_plot[which(to_plot$human_conserved == F & to_plot$human_tf == T), ], mapping = aes(x = human_Chandelier, y = human_Basket), color = "cyan2", size = 1, shape = 15) +
  geom_point(data = to_plot[which(to_plot$human_conserved == T & to_plot$human_tf == F), ], mapping = aes(x = human_Chandelier, y = human_Basket), color = "red", size = 1) +
  geom_point(data = to_plot[which(to_plot$human_specific), ], mapping = aes(x = human_Chandelier, y = human_Basket), color = "green4", size = 2, shape = 15) +
  geom_point(data = to_plot[which(to_plot$human_conserved == T & is.na(to_plot$human_label) == F), ], mapping = aes(x = human_Chandelier, y = human_Basket), color = "magenta", size = 3, shape = 15) +
  xlab("human_Chandelier (log2 expression)") +
  ylab("human_Basket (log2 expression)") +
  xlim(c(0, 8)) +
  ylim(c(0, 8)) 
p1 <- p1 + geom_text_repel(data = to_plot[which(to_plot$human_conserved == T & is.na(to_plot$human_label) == F), ], aes(x = human_Chandelier, y = human_Basket, label = human_label),
                           nudge_x = 5,
                           nudge_y = 0.5,
                           size = 3)
p1 <- p1 +  geom_abline(intercept = 0)
p1 <- p1 + theme_bw() + theme(panel.border = element_blank(), aspect.ratio = 1, panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

p2 <- ggplot() +
  geom_point(data = to_plot[which(to_plot$marmoset_conserved == F & to_plot$marmoset_tf == F), ], mapping = aes(x = marmoset_Chandelier, y = marmoset_Basket), color = "grey80", size = 0.5) +
  geom_point(data = to_plot[which(to_plot$marmoset_conserved == F & to_plot$marmoset_tf == T), ], mapping = aes(x = marmoset_Chandelier, y = marmoset_Basket), color = "cyan2", size = 1, shape = 15) +
  geom_point(data = to_plot[which(to_plot$marmoset_conserved == T & to_plot$marmoset_tf == F), ], mapping = aes(x = marmoset_Chandelier, y = marmoset_Basket), color = "red", size = 1) +
  geom_point(data = to_plot[which(to_plot$human_specific), ], mapping = aes(x = marmoset_Chandelier, y = marmoset_Basket), color = "green4", size = 2, shape = 15) +
  geom_point(data = to_plot[which(to_plot$marmoset_conserved == T & is.na(to_plot$marmoset_label) == F), ], mapping = aes(x = marmoset_Chandelier, y = marmoset_Basket), color = "magenta", size = 3, shape = 15) +
  xlab("marmoset_Chandelier (log2 expression)") +
  ylab("marmoset_Basket (log2 expression)") +
  xlim(c(0, 8)) +
  ylim(c(0, 8)) 
p2 <- p2 + geom_text_repel(data = to_plot[which(to_plot$marmoset_conserved == T & is.na(to_plot$marmoset_label) == F), ], aes(x = marmoset_Chandelier, y = marmoset_Basket, label = marmoset_label),
                           nudge_x = 5,
                           nudge_y = 0.5,
                           size = 3)
p2 <- p2 +  geom_abline(intercept = 0)
p2 <- p2 + theme_bw() + theme(panel.border = element_blank(), aspect.ratio = 1, panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

p3 <- ggplot() +
  geom_point(data = to_plot[which(to_plot$mouse_conserved == F & to_plot$mouse_tf == F), ], mapping = aes(x = mouse_Chandelier, y = mouse_Basket), color = "grey80", size = 0.5) +
  geom_point(data = to_plot[which(to_plot$mouse_conserved == F & to_plot$mouse_tf == T), ], mapping = aes(x = mouse_Chandelier, y = mouse_Basket), color = "cyan2", size = 1, shape = 15) +
  geom_point(data = to_plot[which(to_plot$mouse_conserved == T & to_plot$mouse_tf == F), ], mapping = aes(x = mouse_Chandelier, y = mouse_Basket), color = "red", size = 1) +
  geom_point(data = to_plot[which(to_plot$human_specific), ], mapping = aes(x = mouse_Chandelier, y = mouse_Basket), color = "green4", size = 2, shape = 15) +
  geom_point(data = to_plot[which(to_plot$mouse_conserved == T & is.na(to_plot$mouse_label) == F), ], mapping = aes(x = mouse_Chandelier, y = mouse_Basket), color = "magenta", size = 3, shape = 15) +
  xlab("mouse_Chandelier (log2 expression)") +
  ylab("mouse_Basket (log2 expression)") +
  xlim(c(0, 8)) +
  ylim(c(0, 8)) 
p3 <- p3 + geom_text_repel(data = to_plot[which(to_plot$mouse_conserved == T & is.na(to_plot$mouse_label) == F), ], aes(x = mouse_Chandelier, y = mouse_Basket, label = mouse_label),
                           nudge_x = 5,
                           nudge_y = 0.5,
                           size = 3)
p3 <- p3 +  geom_abline(intercept = 0)
p3 <- p3 + theme_bw() + theme(panel.border = element_blank(), aspect.ratio = 1, panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


grid.arrange(p1, p2, p3, nrow = 1)
