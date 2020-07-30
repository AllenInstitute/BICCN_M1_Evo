#Code Underlying SNARE-Seq2 M1 Data Analyses 
#Blue B. Lake
#b1lake@eng.ucsd.edu

# Human Marmoset SNARE-Seq2 TF Motif Families -----------
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


####
##TF motif clusters/families
####

#Load Cluster Family table
load("~/NeMO_analysis_folder/SNARE/Analysis/Human/JASPAR_Motif_Family_Clusters.rda")

#Prepare human Seurat object
load("~/NeMO_analysis_folder/SNARE/Analysis/Human/Zhang_BICCN-H_20190523-20190611_huMOp_Seurat.rda")
DefaultAssay(MOp.atac) <- "chromvar"
Idents(object = MOp.atac) <- "subclass"
order <- c("LAMP5", "SNCG", "VIP", "SST CHODL", "SST", "PVALB", 
           "L2-3 IT", "L5 IT", "L5 ET", "L5-6 NP", "L6 IT", "L6 IT Car3", "L6 CT", "L6b",                                          
           "OPC", "Astro", "Oligo", "Micro-PVM", "VLMC", "Endo")  
Idents(object = MOp.atac) <- factor(Idents(object = MOp.atac), levels = order)


#Prepare marmoset Seurat object
load("~/NeMO_analysis_folder/SNARE/Analysis/Marmoset/Zhang_BICCN-H_20190730_20190903_marMOp_Seurat.rda")
DefaultAssay(marMOp.atac) <- "chromvar"
Idents(object = marMOp.atac) <- "subclass"
order <- c("Lamp5", "Sncg", "Meis2", "Vip", "Sst Chodl", "Sst", "Pvalb", 
           "L2-3 IT", "L5 IT", "L5 ET", "L5-6 NP", "L6 IT", "L6 IT Car3", "L6 CT", "L6b",                                          
           "OPC", "Astro", "Oligo", "Micro-PVM", "VLMC","Peri", "Endo")
Idents(object = marMOp.atac) <- factor(Idents(object = marMOp.atac), levels = order)



#Human Average by Subclass and TF Cluster
hu.cl.motif <- cl.motif[cl.motif$motif.name %in% rownames(MOp.atac),]
rownames(hu.cl.motif) <- hu.cl.motif$motif.name
ave.tf.hs <- AverageExpression(MOp.atac, assays = "chromvar", features = rownames(MOp.atac), slot = "data" )
ave.tf.hs <- ave.tf.hs$chromvar
hu.cl.motif <- hu.cl.motif[rownames(ave.tf.hs),]
ave.tf.hs$cl <- hu.cl.motif$cluster
ave.tf.hs <- na.omit(ave.tf.hs)

ave.tf.hs <- group_by(ave.tf.hs, cl) %>% summarize_if(is.numeric, mean)
ave.tf.hs <- data.frame(ave.tf.hs); rownames(ave.tf.hs) <- ave.tf.hs$cl; ave.tf.hs <- ave.tf.hs[,-1]
tf.hs <- ave.tf.hs

#Marmoset Average by Subclass and TF Cluster
mar.cl.motif <- cl.motif[cl.motif$motif.name %in% rownames(marMOp.atac),]
rownames(mar.cl.motif) <- mar.cl.motif$motif.name
ave.tf.mar <- AverageExpression(marMOp.atac, assays = "chromvar", features = rownames(marMOp.atac), slot = "data" )
ave.tf.mar <- ave.tf.mar$chromvar
mar.cl.motif <- mar.cl.motif[rownames(ave.tf.mar),]
ave.tf.mar$cl <- mar.cl.motif$cluster
ave.tf.mar <- na.omit(ave.tf.mar)

ave.tf.mar <- group_by(ave.tf.mar, cl) %>% summarize_if(is.numeric, mean)
ave.tf.mar <- data.frame(ave.tf.mar); rownames(ave.tf.mar) <- ave.tf.mar$cl; ave.tf.mar <- ave.tf.mar[,-1]
tf.mar <- ave.tf.mar

colnames(tf.mar) <- paste("mar", colnames(tf.mar), sep = "_")
tf.hs.mar <- cbind(tf.hs, tf.mar[rownames(tf.hs),])
tf.hs.mar <- data.frame(tf.hs.mar)
tf.hs.mar$cluster <- gsub("cluster_","cl",rownames(tf.hs.mar))
write.table(tf.hs.mar, file="hMOp_marMOp_TFBS_Family_Comparison_full_table.txt", sep = "\t", row.names=TRUE, col.names=TRUE)



#Lamp5
to.use <- c("cluster_3","cluster_9","cluster_14","cluster_21","cluster_23","cluster_25")
lamp5 <- tf.hs.mar[to.use,c("LAMP5","mar_Lamp5")]
rownames(lamp5) <- to.use
lamp5[is.na(lamp5)] <- 0
colnames(lamp5) <- c("human","marmoset")
range(lamp5)
write.table(lamp5, file="HuMOp_marMOp_TFBS_Family_Comparison_LAMP5_table.txt", sep = "\t", row.names=TRUE, col.names=TRUE)

#Plot Fig. 4h
ggHeat(lamp5, rescaling = "none", clustering = "none", x.lab.size = 8, y.lab.size = 11,
       heatscale = c(low = "skyblue", mid ="white", high = "#005a32")
) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  ggtitle("LAMP5")



#L5 IT
to.use <- c("cluster_1","cluster_8","cluster_12","cluster_34","cluster_46","cluster_48")
L5IT <- tf.hs.mar[to.use,c("L5.IT","mar_L5.IT")]
rownames(L5IT) <- to.use
L5IT[is.na(L5IT)] <- 0
colnames(L5IT) <- c("human","marmoset")
range(L5IT)
write.table(L5IT, file="HuMOp_marMOp_TFBS_Family_Comparison_L5-IT_table.txt", sep = "\t", row.names=TRUE, col.names=TRUE)

#Plot Fig. 4h
ggHeat(L5IT, rescaling = "none", clustering = "none", x.lab.size = 8, y.lab.size = 11,
       heatscale = c(low = "skyblue", mid ="white", high = "#005a32")
) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  ggtitle("L5-IT")






