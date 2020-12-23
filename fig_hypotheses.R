library(ggplot2)
library(qs)
library(data.table)

pa1 = qread("./output/cog-plot-0.qs")
pa2 = qread("./output/post-plot-0.qs")
pb1 = qread("./output/cog-plot-ImmEsc.qs")
pb2 = qread("./output/post-plot-ImmEsc.qs")
pc1 = qread("./output/cog-plot-S.qs")
pc2 = qread("./output/post-plot-S.qs")
pd1 = qread("./output/cog-plot-LatentPeriod.qs")
pd2 = qread("./output/post-plot-LatentPeriod.qs")

plot = cowplot::plot_grid(pa1, pa2, pb1, pb2, pc1, pc2, pd1, pd2, nrow = 4, labels = c("A", "", "B", "", "C", "", "D", ""), label_size = 10, align = "hv", axis = "bottom")
ggsave("./output/fig-hypotheses.pdf", plot, width = 34, height = 20, units = "cm", useDingbats = FALSE)
ggsave("./output/fig-hypotheses.png", plot, width = 34, height = 20, units = "cm")
plot
