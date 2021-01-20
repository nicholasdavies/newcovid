p1 = qread("./output/ageplot-infdur.qs") + labs(title = "Longer infectious period")
p2 = qread("./output/ageplot-ch_u.qs") + labs(title = "Increased susceptibility among children")


cowplot::plot_grid(p1, p2, nrow = 1, labels = LETTERS, label_size = 14)
ggsave("./output/fig_s2.pdf", width = 35, height = 14, units = "cm", useDingbats = FALSE)
ggsave("./output/fig_s2.png", width = 35, height = 14, units = "cm")
