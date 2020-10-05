library(dplyr)
library(Seurat)

de_genes <- read.csv("Supplementary Table 7")
sample.combined <- readRDS("~/NeMO_analysis_folder/Transcriptomics/cross_species_integration/sample.combined_inh_integration.RDS")

subclasses <- c("Lamp5", "Sncg", "Vip", "Sst Chodl", "Sst", "Pvalb")

human <- de_genes[which(de_genes$conservation == "human"), ]
marmoset <- de_genes[which(de_genes$conservation == "marmoset"), ]
mouse <- de_genes[which(de_genes$conservation == "mouse"), ]

#make logFC plot
to_plot <- data.frame(species = NA, subclass = NA, prop_genes = 0, logFC = 0)
logFC <- 0

for(j in 1:8){
  logFC <- logFC + 0.25 
  for(i in 1:length(subclasses)){
    tmp <- data.frame(species = "human", subclass = subclasses[i])
    sub_loc <- which(human$cluster == subclasses[i])
    tmp$prop_genes <- length(which(human$human_avg.1[sub_loc] - human$marmoset_avg.1[sub_loc] > logFC & 
                                     human$human_avg.1[sub_loc] - human$mouse_avg.1[sub_loc] > logFC))       / length(sub_loc)
    tmp$logFC = logFC
    to_plot <- rbind(to_plot, tmp)
    
    tmp <- data.frame(species = "marmoset", subclass = subclasses[i])
    sub_loc <- which(marmoset$cluster == subclasses[i])
    tmp$prop_genes <- length(which(marmoset$marmoset_avg.1[sub_loc] - marmoset$human_avg.1[sub_loc] > logFC & 
                                     marmoset$marmoset_avg.1[sub_loc] - marmoset$mouse_avg.1[sub_loc] > logFC))       / length(sub_loc)
    tmp$logFC = logFC
    to_plot <- rbind(to_plot, tmp)
    
    tmp <- data.frame(species = "mouse", subclass = subclasses[i])
    sub_loc <- which(mouse$cluster == subclasses[i])
    tmp$prop_genes <- length(which(mouse$mouse_avg.1[sub_loc] - mouse$human_avg.1[sub_loc] > logFC & 
                                     mouse$mouse_avg.1[sub_loc] - mouse$marmoset_avg.1[sub_loc] > logFC))       / length(sub_loc)
    tmp$logFC = logFC
    to_plot <- rbind(to_plot, tmp)
  }
}
to_plot <- to_plot[-1, ]
to_plot$species <- paste(to_plot$species, ">", to_plot$logFC, sep = " ")


custom_color <- subclasses
names(custom_color) <- seurat_obj$subclass_color[match(custom_color, seurat_obj$subclass_label)]


ggplot(data = to_plot, mapping = aes(x = species, y = prop_genes, color = subclass)) +
  geom_point(size = 6, alpha = 0.5) +
  scale_color_manual(values = names(custom_color)) +
  ggtitle("Subclass logFC against other two species") +
  ylab("Proportion species-enriched DE genes") +
  xlab("LogFC enrichment")+
  ylim(c(0, 1))+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12, colour = c(rep("royalblue1", j), 
                                                                                rep("maroon4", j), 
                                                                                rep("sienna2", j))),
        axis.text.y = element_text(size = 12))

#make prop plot
to_plot <- data.frame(species = NA, subclass = NA, prop_genes = 0, prop.th = 0)
prop.th <- 0

for(j in 1:9){
  prop.th <- prop.th + 0.1
  for(i in 1:length(subclasses)){
    tmp <- data.frame(species = "human", subclass = subclasses[i])
    sub_loc <- which(human$cluster == subclasses[i])
    tmp$prop_genes <- length(which(human$human_prop.1[sub_loc] - human$marmoset_prop.1[sub_loc] > prop.th & 
                                     human$human_prop.1[sub_loc] - human$mouse_prop.1[sub_loc] > prop.th))       / length(sub_loc)
    tmp$prop.th = prop.th
    to_plot <- rbind(to_plot, tmp)
    
    tmp <- data.frame(species = "marmoset", subclass = subclasses[i])
    sub_loc <- which(marmoset$cluster == subclasses[i])
    tmp$prop_genes <- length(which(marmoset$marmoset_prop.1[sub_loc] - marmoset$human_prop.1[sub_loc] > prop.th & 
                                     marmoset$marmoset_prop.1[sub_loc] - marmoset$mouse_prop.1[sub_loc] > prop.th))       / length(sub_loc)
    tmp$prop.th = prop.th
    to_plot <- rbind(to_plot, tmp)
    
    tmp <- data.frame(species = "mouse", subclass = subclasses[i])
    sub_loc <- which(mouse$cluster == subclasses[i])
    tmp$prop_genes <- length(which(mouse$mouse_prop.1[sub_loc] - mouse$human_prop.1[sub_loc] > prop.th & 
                                     mouse$mouse_prop.1[sub_loc] - mouse$marmoset_prop.1[sub_loc] > prop.th))       / length(sub_loc)
    tmp$prop.th = prop.th
    to_plot <- rbind(to_plot, tmp)
  }
}
to_plot <- to_plot[-1, ]
to_plot$species <- paste(to_plot$species, ">", to_plot$prop.th, sep = " ")


custom_color <- subclasses
names(custom_color) <- seurat_obj$subclass_color[match(custom_color, seurat_obj$subclass_label)]


ggplot(data = to_plot, mapping = aes(x = species, y = prop_genes, color = subclass)) +
  geom_point(size = 3, alpha = 0.5) +
  scale_color_manual(values = names(custom_color)) +
  ggtitle("Subclass proportion difference against other two species") +
  ylab("Proportion species-specific DE genes") +
  xlab("Proportion enrichment") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7, colour = c(rep("royalblue1", j), 
                                                                               rep("maroon4", j), 
                                                                               rep("sienna2", j))))
