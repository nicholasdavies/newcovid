p1 = qread("./output/ageplot-0.qs") + labs(title = "Increased infectiousness")
p2 = qread("./output/ageplot-ImmEsc.qs") + labs(title = "Immune escape")
p3 = qread("./output/ageplot-S.qs") + labs(title = "Increased susceptibility among children")
p4 = qread("./output/ageplot-LatentPeriod.qs") + labs(title = "Shorter generation time")


cowplot::plot_grid(p1, p2, p3, p4, nrow = 2, labels = LETTERS, label_size = 14)
ggsave("./output/fig_s2.pdf", width = 35, height = 20, units = "cm", useDingbats = FALSE)
ggsave("./output/fig_s2.png", width = 35, height = 20, units = "cm")
