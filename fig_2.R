library(data.table)
library(ggplot2)
library(cowplot)
library(lubridate)

# assumes g_rel and g_t from relativized_growth_rate_and_IPO... are loaded
# assumes p_rt_sgtf from rt_sgtf_scatter.R is loaded

p_muller = readRDS("./data/model 1_plot multinomial spline fit_muller plot fit.rds")
pal_variants = c(
    "B" = "#FF3A00",
    "B.1.98" = "#F78500",
    "B.40" = "#A6AA00",
    "B.1.1" = "#00C400",
    "B.1.1.257" = "#00D756",
    "B.1.1.1" = "#00DEE7",
    "B.1.1.315" = "#00CAFF",
    "B.1.177" = "#806CFF",
    "VOC 202012/01" = "#FF00FF",
    "minority variants" = "grey75",
    "B.1.1.301" = "#C3BC3F",
    "B.1.177.7" = "#453798"
)

pal_regions = c("#6388b4",  "#ffae34",  "#ef6f6a",  "#8cc2ca",  "#55ad89",  "#c3bc3f",  "#bb7693",  "#baa094",  "#a9b5ae",  "#767676")

p_2 = cowplot::plot_grid(
    cowplot::plot_grid(
        g_rel + theme_cowplot(font_size = 10) + theme(legend.position = "bottom") + 
            labs(title = NULL, colour = NULL) + scale_colour_manual(values = pal_variants) + guides(colour = guide_legend(nrow = 1, override.aes = list(size = 5))),
        g_t + theme_cowplot(font_size = 10) + theme(legend.position = "none") + labs(title = NULL, colour = NULL) + scale_colour_manual(values = pal_variants),
        ncol = 1, labels = LETTERS, label_size = 10),
    p_muller + scale_x_date(date_breaks = "1 month", date_labels = "%b", limits = ymd(c("2020-02-05", "2021-01-31"))) + 
        theme_cowplot(font_size = 10) + theme(legend.position = "bottom", strip.background = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5)) +
        labs(x = "Sample date", fill = NULL) + scale_fill_manual(values = pal_variants),
    p_rt_sgtf,
    nrow = 1, rel_widths = c(8, 10, 8), labels = c("", "C", "D"), label_size = 10
)

ggsave("./output/fig-2.pdf", p_2, width = 50, height = 17, units = "cm", useDingbats = FALSE)
ggsave("./output/fig-2.png", p_2, width = 50, height = 17, units = "cm")
