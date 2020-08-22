library(dendextend)
library(feather)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(scrattch.hicat)


calc_layer_prop <- function(roi, cluster) {
  cl.layer.cnt <- matrix(unlist(tapply(roi, cluster, 
                                       function(x) table(x))), 
                         ncol = nlevels(roi), byrow = TRUE, 
                         dimnames = list(levels(cluster), 
                                         levels(roi)))
  cl.layer.prop <- apply(apply(cl.layer.cnt, 2, 
                               function(x) x / sum(x)), 1, 
                         function(y) y / sum(y))
  return(cl.layer.prop)
}


order_rows <- function(m1) {
  max.row <- apply(m1, 1, function(x) {
    min(which(x > 0))
  })
  m2 <- m1[order(max.row), ]
  if (any(is.na(row.names(m2)))) {
    na.row <- is.na(row.names(m2))
    na.order <- c(which(! na.row), which(na.row))
    m2 <- m2[na.order, ]
  }
  return(m2)
}


fig.dir <- "misc_figures/output/Figure1/"
dir.create(fig.dir, recursive = TRUE, showWarnings = FALSE)


species <- c("human", "marmoset", "mouse")

anno.l <- readRDS("misc_figures/data/species.anno.rds")


#### Figure 1c-e. Plot dendrograms ####
for (species1 in species) {
  dend.fn <- paste0("misc_figures/data/", species1, "/", species1, "_dend.RData")
  dend <- readRDS(dend.fn)
  
  # Update dendrogram colors
  anno <- anno.l[[species1]]
  l.color <- data.frame(cluster = unique(anno$cluster_label))
  l.color$color <- anno$cluster_color[match(l.color$cluster, anno$cluster_label)]
  colorVectorNew <- makeColorsUnique(l.color$color)
  use.color <- setNames(colorVectorNew, l.color$cluster)
  
  dend <- dend %>% set("labels_col", use.color[labels(dend)])
  dend <- dend %>% set("leaves_col", use.color[labels(dend)])
  
  pdf(file = paste0(fig.dir, species1, "_dend.pdf"), width = 16, height = 5)
  par(mar = c(10, 2, 2, 2))
  plot(dend, main = species1)
  dev.off()
}



#### Figure 1c. Calc human layer distrib ####
cl.anno <- unique(anno.l[["human"]][, c("cluster_id", "cluster_label")])
cl.anno <- cl.anno[order(cl.anno$cluster_id), ]

anno <- as.data.frame(anno.l[["human_ss"]])
anno <- droplevels(subset(anno, cluster_label != "Astro L5-6 RORB CYBRD1"))
anno$layer_label <- as.factor(anno$layer_label)
anno$cluster_label <- factor(anno$cluster_label, levels = cl.anno$cluster_label)

cl.layer.prop <- calc_layer_prop(anno$layer_label, anno$cluster_label)


hm.colors <- colorRampPalette(c("white", brewer.pal(9, "YlOrRd")))(100)
pdf(file = paste0(fig.dir, "human_cl_layer_prop.pdf"), width = 12, height = 2, onefile = FALSE)
pheatmap(cl.layer.prop, cluster_rows = FALSE, cluster_cols = FALSE, 
         cellheight = 5, cellwidth = 5,
         color = hm.colors, show_colnames = FALSE)
dev.off()




#### Figure 1c-e. DNA-m/ATAC metadata ####
pdf(file = paste0(fig.dir, "epi_cluster_tracks.pdf"), width = 12, height = 6, onefile = TRUE)
for (species1 in species) {
  meta.fn <- paste0("misc_figures/data/", species1, "/Species_metadata_Cytosplore - ", species1, ".csv")
  meta <- read.csv(meta.fn)
  meta <- meta[order(meta$tree_order), ]
  
  for (epi.cl in c("ATAC_cluster_label", "DNAm_cluster_label")) {
    if (! all(is.na(meta[, epi.cl]))) {
      # conf1 <- table(meta$cluster_label, meta$DNAm_cluster_name)
      conf1 <- table(meta$cluster_label, meta[, epi.cl], useNA = "ifany")
      # conf1 <- sweep(conf1, 2, colSums(conf1), "/")
      conf2 <- as.data.frame.matrix(conf1)
      conf3 <- conf2[match(meta$cluster_label, row.names(conf2)), ]
      conf4 <- order_rows(t(conf3))
      if (species1 == "marmoset") {
        if (epi.cl == "ATAC_cluster_label") conf4 <- conf4[c(2:4,1,5:nrow(conf4)), ]
        if (epi.cl == "DNAm_cluster_label") conf4 <- conf4[c(2:11,1,12:nrow(conf4)), ]
      }
      pheatmap(conf4, cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 6, 
               color = c("white", "black"), border_color = NA,
               main = paste(species1, "-", epi.cl))
    }
  }
}
dev.off()



#### Figure 1f. Calc subclass proportions ####
subclasses <- c("Lamp5", "Sncg", "Vip", "Sst Chodl", "Sst", "Pvalb", 
                "L2/3 IT", "L5 IT", "L6 IT", "L6 IT Car3", "L5 ET", "L5/6 NP", 
                "L6 CT", "L6b", "Meis2", "OPC", "Astro", "Oligo", "Endo", 
                "VLMC", "SMC", "Peri", "Micro-PVM")


anno2.l <- lapply(anno.l[1:3], function(x) x[, c("species", "donor_id", 
                                                 "class_label", "subclass_label")])
anno.comb <- Reduce("rbind", anno2.l)

# Neuron class
g.neuron <- anno.comb %>% 
  filter(class_label != "Non-Neuronal") %>% 
  count(species, donor_id, class_label) %>% 
  group_by(species, donor_id) %>% 
  mutate(freq = n / sum(n)) %>% 
  ggplot(aes(x = factor(class_label), y = freq, color = species, fill = species)) +
  stat_summary(fun.data = mean_se, position = position_dodge(width = 0.75), 
               geom="errorbar", width = 0.5) +
  stat_summary(fun.y = mean, position = position_dodge(width = 0.75), geom="bar", width = 0.65) +
  geom_point(position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.15), 
             size = 0.5, color = "black", show.legend = FALSE) +
  xlab("") + ylab("Neuron proportion") + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(g.neuron, file = paste0(fig.dir, "neuron.prop.pdf"), width = 2.5, height = 3)

# Inhibitory
g.inh <- anno.comb %>% 
  mutate(subclass_label = factor(subclass_label, levels = subclasses)) %>% 
  filter(class_label == "GABAergic") %>% 
  count(species, donor_id, subclass_label) %>% 
  group_by(species, donor_id) %>% 
  mutate(freq = n / sum(n)) %>% 
  ggplot(aes(x = factor(subclass_label), y = freq, color = species, fill = species)) +
  stat_summary(fun.data = mean_se, position = position_dodge(width = 0.75), 
               geom="errorbar", width = 0.5) +
  stat_summary(fun.y = mean, position = position_dodge(width = 0.75), geom="bar", width = 0.65) +
  geom_point(position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.15), 
             size = 0.5, color = "black", show.legend = FALSE) +
  xlab("") + ylab("GABAergic neuron proportion") + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(g.inh, file = paste0(fig.dir, "inh.prop.pdf"), width = 4, height = 3)

# Excitatory
g.exc <- anno.comb %>% 
  mutate(subclass_label = factor(subclass_label, levels = subclasses)) %>% 
  filter(class_label == "Glutamatergic") %>% 
  count(species, donor_id, subclass_label) %>% 
  group_by(species, donor_id) %>% 
  mutate(freq = n / sum(n)) %>% 
  ggplot(aes(x = factor(subclass_label), y = freq, color = species, fill = species)) +
  stat_summary(fun.data = mean_se, position = position_dodge(width = 0.75), 
               geom="errorbar", width = 0.5) +
  stat_summary(fun.y = mean, position = position_dodge(width = 0.75), geom="bar", width = 0.65) +
  geom_point(position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.15), 
             size = 0.5, color = "black", show.legend = FALSE) +
  xlab("") + ylab("Glutamatergic neuron proportion") + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1))
# plot(g.exc)
ggsave(g.exc, file = paste0(fig.dir, "exc.prop.pdf"), width = 4.5, height = 3)

# Non-Neuronal
g.glia <- anno.comb %>% 
  mutate(subclass_label = factor(subclass_label, levels = subclasses)) %>%
  filter(class_label == "Non-Neuronal") %>% 
  count(species, donor_id, subclass_label) %>% 
  group_by(species, donor_id) %>% 
  mutate(freq = n / sum(n)) %>% 
  ggplot(aes(x = factor(subclass_label), y = freq, color = species, fill = species)) +
  stat_summary(fun.data = mean_se, position = position_dodge(width = 0.75), 
               geom="errorbar", width = 0.5) +
  stat_summary(fun.y = mean, position = position_dodge(width = 0.75), geom="bar", width = 0.65) +
  geom_point(position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.15), 
             size = 0.5, color = "black", show.legend = FALSE) +
  xlab("") + ylab("Non-neuronal proportion") + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(g.glia, file = paste0(fig.dir, "glia.prop.pdf"), width = 4.5, height = 3)



#### Compare subclass proportions - stats ####
anno.aov <- list()

anno1 <- anno.comb %>% 
  filter(class_label != "Non-Neuronal") %>% 
  count(species, donor_id, class_label) %>% 
  group_by(species, donor_id) %>% 
  mutate(freq = n / sum(n))
for (class1 in c("GABAergic", "Glutamatergic")) {
  anno.aov[[class1]] <- anno1 %>% filter(class_label == class1)
}

anno1 <- anno.comb %>% 
  filter(class_label == "GABAergic") %>% 
  count(species, donor_id, subclass_label) %>% 
  group_by(species, donor_id) %>% 
  mutate(freq = n / sum(n))
for (subclass1 in subclasses[1:6]) {
  anno.aov[[subclass1]] <- anno1 %>% filter(subclass_label == subclass1)
}

anno1 <- anno.comb %>% 
  filter(class_label == "Glutamatergic") %>% 
  count(species, donor_id, subclass_label) %>% 
  group_by(species, donor_id) %>% 
  mutate(freq = n / sum(n))
for (subclass1 in subclasses[7:14]) {
  anno.aov[[subclass1]] <- anno1 %>% filter(subclass_label == subclass1)
}

aov.th <- 0.05 / length(anno.aov)
anno.stats <- sapply(anno.aov, function(x) {
  aov1 <- aov(x$freq ~ x$species)
  pval1 <- summary(aov1)[[1]][["Pr(>F)"]][1]
  if (pval1 < aov.th) {
    tukey1 <- TukeyHSD(aov1)
    sign(tukey1[[1]][, "diff"]) * tukey1[[1]][, "p adj"]
  } else {
    rep(NA, 3)
  }
})

anno.stats <- na.omit(t(anno.stats))
anno.stats2 <- anno.stats * nrow(anno.stats)
anno.stats2[anno.stats2 > 1] <- 1
anno.stats2[anno.stats2 < -1] <- -1
write.csv(anno.stats2, file = paste0(fig.dir, "cell_class_prop_adjP.csv"))
