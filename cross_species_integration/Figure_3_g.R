library(dplyr)
library(matrixStats)
library(Matrix)
library(feather)

data <- read.csv("~/NeMO_analysis_folder/Transcriptomics/cross_species_integration/additional_intermediate_files/MediumHighIsoforms.tsv", sep = "\t")
head(data)
data$Relative_Human_Prop_diff <-  (data$Human_Isoform_Proportion - data$Mouse_Isoform_Proportion) / 
  (data$Human_Isoform_Proportion + data$Mouse_Isoform_Proportion) 

data$species_enriched_prop <- NA
data$species_enriched_prop[which(data$Relative_Human_Prop_diff > 0)] <- "Human"
data$species_enriched_prop[which(data$Relative_Human_Prop_diff < 0)] <- "Mouse"

#calculate data for box plot
anno <- read_feather("~/NeMO_analysis_folder/Transcriptomics/human_10x_SS_integration/Figure_ED1_source_data/anno.feather")
anno$subclass_label <- sub("/", "-", anno$subclass_label)
anno$subclass_label <- sub(" ", "_", anno$subclass_label)
data$Class <- anno$class_label[match(data$Subclass, anno$subclass_label)]
table(data$Class, data$Subclass)


# For each subclass and species:
#   R = (P1 - P2)/(P1 + P2) 
# S = number of isoforms where R < -0.8 (mouse) or R > 0.8 (human)
# I = number of isoforms for any value of R (i.e. all moderately expressed isoforms)
# Plot S / I

data$species_subclass <- paste(data$species_enriched_prop, data$Subclass, sep = "_")
to_plot <- data[which(data$Relative_Human_Prop_diff > 0.8 | data$Relative_Human_Prop_diff < -0.8), ]
s_table <- data.frame(table(to_plot$species_subclass))
i_table <- data.frame(table(data$species_subclass))
i_table <- i_table[-grep("NA", i_table$Var1), ]
colnames(i_table) <- paste(colnames(i_table), "all", sep = "_")
to_plot <- cbind(s_table, i_table)
to_plot$prop <- to_plot$Freq/to_plot$Freq_all
to_plot$Var1_all <- sub("Human_", "", to_plot$Var1_all)
to_plot$Var1_all <- sub("Mouse_", "", to_plot$Var1_all)
to_plot$species <- NA
to_plot$species[grep("Human", to_plot$Var1)] <- "Human"
to_plot$species[grep("Mouse", to_plot$Var1)] <- "Mouse"
to_plot$class <- anno$class_label[match(to_plot$Var1_all, anno$subclass_label)]
to_plot$species_class <- paste(to_plot$class, to_plot$species, sep = "_")
head(to_plot)

library(ggplot2)
ggplot(data = to_plot, mapping = aes(x = species_class, y= prop, color = species)) +
  geom_boxplot() + 
  geom_jitter(shape=16, position=position_jitter(0.1), color = "black")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  ylim(0, 0.65) +
  ylab("Proportion of number of enriched isoforms per subclass \n (9-fold enriched isoforms / all isoforms)")



