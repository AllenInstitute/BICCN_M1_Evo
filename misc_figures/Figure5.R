library(ggplot2)


#### Figure 5d. Plot M1/MTG cell type layer distributions ####
layer.thk <- read.csv("misc_figures/data/human/mtg_m1_layer_thickness.csv")
layer.m1 <- subset(layer.thk, Area == "M1")
layer.mtg <- subset(layer.thk, Area == "MTG")

layer.map <- read.csv("misc_figures/data/human/layer_transfer.csv")
layer.map$M1_depth <- layer.m1$Mid.layer.relative.pial.depth[match(layer.map$M1_layer, layer.m1$Layer)]
layer.map$MTG_depth <- layer.mtg$Mid.layer.relative.pial.depth[match(layer.map$MTG_layer_transfer, layer.mtg$Layer)]

cl.df <- data.frame(m1.mean = tapply(layer.map$M1_depth, layer.map$cluster, mean), 
                    m1.sd = tapply(layer.map$M1_depth, layer.map$cluster, sd), 
                    mtg.mean = tapply(layer.map$MTG_depth, layer.map$cluster, mean), 
                    mtg.sd = tapply(layer.map$MTG_depth, layer.map$cluster, sd))

g.layer <- ggplot(cl.df, aes(x = m1.mean, y = mtg.mean)) +
  geom_vline(xintercept = c(layer.m1$Relative.pial.depth[1:4]), color = "grey80") +
  geom_hline(yintercept = c(layer.mtg$Relative.pial.depth[1:5]), color = "grey80") +
  geom_text(data = layer.m1[-1, ], aes(x = Mid.layer.relative.pial.depth, y = 0, 
                                       label = Layer), size = 5) +
  geom_text(data = layer.mtg[-1, ], aes(x = -0.05, y = Mid.layer.relative.pial.depth, 
                                        label = Layer), size = 5) +
  geom_errorbarh(aes(xmin = m1.mean - m1.sd, xmax = m1.mean + m1.sd)) +
  geom_errorbar(aes(ymin = mtg.mean - mtg.sd, ymax = mtg.mean + mtg.sd)) +
  geom_point(size = 2) +
  coord_equal() +
  scale_y_continuous(trans = "reverse") +
  xlab("M1 relative depth from pia") +
  ylab("MTG relative depth from pia") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot(g.layer)

ggsave(g.layer, filename = "output/Figure5/mtg_m1_exc_layer_scatter.pdf", width = 5, height = 5)
