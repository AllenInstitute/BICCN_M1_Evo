library(Seurat)
library(dplyr)
library(Matrix)
library(matrixStats)
library(gplots)
library(ggplot2)
library(feather)
##############################################################
############# colored correlation scatter ####################
##############################################################

#Custom color points
hum_marm <- read.csv("~/NeMO_analysis_folder/Transcriptomics/cross_species_integration/additional_intermediate_files/glial_files/human_marmoset_markers.csv")
hum_mouse <- read.csv("~/NeMO_analysis_folder/Transcriptomics/cross_species_integration/additional_intermediate_files/glial_files/human_mouse_markers.csv")
mouse_marm <- read.csv("~/NeMO_analysis_folder/Transcriptomics/cross_species_integration/additional_intermediate_files/glial_files/marmoset_mouse_markers.csv")


hum_marm <- hum_marm[grep("Astro_1", hum_marm$cluster), ]
hum_mouse <- hum_mouse[grep("Astro_1", hum_mouse$cluster), ]
mouse_marm <- mouse_marm[grep("Astro_1", mouse_marm$cluster), ]

FDR <- 0.01
avglogfc <- 2


hum_marm <- hum_marm[which(hum_marm$p_val_adj < FDR), ]
hum_marm <- hum_marm[which(hum_marm$avg_logFC > avglogfc), ]

hum_mouse <- hum_mouse[which(hum_mouse$p_val_adj < FDR), ]
hum_mouse <- hum_mouse[which(hum_mouse$avg_logFC > avglogfc), ]

mouse_marm <- mouse_marm[which(mouse_marm$p_val_adj < FDR), ]
mouse_marm <- mouse_marm[which(mouse_marm$avg_logFC > avglogfc), ]


avg_hum_marm <- read.csv("~/NeMO_analysis_folder/Transcriptomics/cross_species_integration/additional_intermediate_files/glial_files/glia_human_marmoset_avg_expr.csv")
avg_mouse <- read.csv("~/NeMO_analysis_folder/Transcriptomics/cross_species_integration/additional_intermediate_files/glial_files/glia_marmoset_mouse_avg_expr.csv")

avg_hum <- avg_hum_marm[ , grep("human", colnames(avg_hum_marm))]
avg_marm <- avg_hum_marm[ , grep("marmoset", colnames(avg_hum_marm))]
avg_mouse <- avg_mouse[ , grep("mouse", colnames(avg_mouse))]

colnames(avg_hum)
average_log2_expr <- data.frame(human = avg_hum[ , 6])
average_log2_expr$marmoset <- avg_marm[ , 6]
average_log2_expr$mouse <- avg_mouse[ , 6]
rownames(average_log2_expr) <- avg_hum_marm$X

average_log2_expr$hum_marm_col <- "black"
average_log2_expr$hum_marm_col[which(rownames(average_log2_expr) %in% hum_marm$gene[grep("human", hum_marm$cluster)])] <- "royalblue1"
average_log2_expr$hum_marm_col[which(rownames(average_log2_expr) %in% hum_marm$gene[grep("marmoset", hum_marm$cluster)])] <- "maroon3"

average_log2_expr$hum_mouse_col <- "black"
average_log2_expr$hum_mouse_col[which(rownames(average_log2_expr) %in% hum_mouse$gene[grep("human", hum_mouse$cluster)])] <- "royalblue1"
average_log2_expr$hum_mouse_col[which(rownames(average_log2_expr) %in% hum_mouse$gene[grep("mouse", hum_mouse$cluster)])] <- "sienna2"

average_log2_expr$mouse_marm_col <- "black"
average_log2_expr$mouse_marm_col[which(rownames(average_log2_expr) %in% mouse_marm$gene[grep("marmoset", mouse_marm$cluster)])] <- "maroon3"
average_log2_expr$mouse_marm_col[which(rownames(average_log2_expr) %in% mouse_marm$gene[grep("mouse", mouse_marm$cluster)])] <- "sienna2"

max(average_log2_expr[ , 1:3])
ggplot(data = average_log2_expr) +
  geom_point(mapping = aes(x = human, y = marmoset),  color = average_log2_expr$hum_marm_col)+
  xlim(0,8) + ylim(0,8)
ggplot(data = average_log2_expr) +
  geom_point(mapping = aes(x = human, y = mouse),  color = average_log2_expr$hum_mouse_col)+
  xlim(0,8) + ylim(0,8)
ggplot(data = average_log2_expr) +
  geom_point(mapping = aes(x = marmoset, y = mouse),  color = average_log2_expr$mouse_marm_col)+
  xlim(0,8) + ylim(0,8)

cor(average_log2_expr$human, average_log2_expr$marmoset, method = "spearman")
cor(average_log2_expr$human, average_log2_expr$mouse, method = "spearman")
cor(average_log2_expr$mouse, average_log2_expr$marmoset, method = "spearman")
