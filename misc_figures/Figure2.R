library(tidyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)


fig.dir <- "misc_figures/output/Figure2/"
dir.create(fig.dir, recursive = TRUE, showWarnings = FALSE)


#### Figure 2h. Bar plot stats ####
for (class1 in c("inh", "exc")) {
  leaves <- read.csv(paste0("misc_figures/data/", class1, ".nleaves.v2.txt"), sep = "\t")
  leaves$comp <- factor(sub(".", "_", leaves$comp, fixed = TRUE))
  
  aov1 <- aov(nleaves ~ celltype*comp, leaves)
  summary(aov1)
  
  subclasses <- levels(leaves$celltype)
  anno.aov <- list()
  
  for (subclass1 in subclasses) {
    anno.aov[[subclass1]] <- leaves %>% filter(celltype == subclass1)
  }
  
  aov.th <- 0.05 / length(anno.aov)
  anno.stats <- sapply(anno.aov, function(x) {
    aov1 <- aov(x$nleaves ~ x$comp)
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
  write.csv(anno.stats2, paste0(fig.dir, "Figure2_", class1, "_leaves_stats.csv"))
  
}



#### Figure 2i. MetaNeighbor gene sets ####
mn.df <- read.csv("misc_figures/data/Supplementary Table 9.csv")
mn.df$gene_set_type <- factor(mn.df$gene_set_type, 
                              levels = c("Other", "Signaling", 
                                         "Cell Adhesion", "Ion Channel"))

roc.tests <- c("human_meanROC", "marmoset_meanROC", "mouse_meanROC", 
               "human_marmoset_meanROC", "human_mouse_meanROC", 
               "marmoset_mouse_meanROC")

mn.l <- mn.df %>% 
  subset(cell_class == "GABAergic") %>%
  arrange(gene_set_type) %>% 
  gather(roc.tests, key = "comp", value = "ROC")

mn.l$comp <- factor(mn.l$comp, levels = roc.tests)
levels(mn.l$comp) <- c("Human", "Marmoset", "Mouse", 
                       "Human vs. Marmoset", "Human vs. Mouse", 
                       "Marmoset vs. Mouse")
mn.l$comp_type <- ifelse(mn.l$comp %in% c("Human", "Marmoset", "Mouse"), 
                         "Within-species", "Cross-species")

geneset.pal <- c("grey", "dark green", "orange", "purple")

g.scatter <- ggplot(mn.l, aes(x = within_species_meanROC, y = ROC)) +
  facet_grid(. ~ comp) +
  geom_hline(yintercept = 0.5, color = "grey", size = 0.25) +
  geom_vline(xintercept = 0.5, color = "grey", size = 0.25) +
  geom_abline(slope = 1, intercept = 0, color = "grey", size = 0.25) +
  geom_point(alpha = 0.5, aes(color = gene_set_type)) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  scale_colour_manual(values=geneset.pal) +
  coord_fixed(xlim = c(0.4, 1), ylim = c(0.4, 1)) +
  scale_x_continuous(breaks = c(0.5, 0.7, 0.9)) +
  scale_y_continuous(breaks = c(0.5, 0.7, 0.9)) +
  xlab("Within-species mean AUROC") +
  ylab("AUROC") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
plot(g.scatter)
ggsave(g.scatter, filename = paste0(fig.dir, "inh_species_AUROC.pdf"), 
       width = 10, height = 3)


mn.l %>% 
  group_by(comp) %>%
  group_modify(~ broom::tidy(lm(ROC ~ within_species_meanROC, data = .x)))

