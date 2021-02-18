library(data.table)
library(ggplot2)
library(cowplot)
library(lubridate)
library(ggthemes)

pal_variants = c(
    "B" = "#FF3A00",
    "B.1.98" = "#F78500",
    "B.40" = "#A6AA00",
    "B.1" = "#00C400",
    "B.1.1" = "#00D756",
    "B.1.1.1" = "#00DEE7",
    "B.1.1.315" = "#00CAFF",
    "B.1.177" = "#806CFF",
    "VOC 202012/01" = "#FF00FF",
    "minority variants" = "grey75",
    "B.1.1.301" = "#C3BC3F",
    "B.1.177.7" = "#453798"
)

pal_regions = c("#6388b4", "#ffae34", "#ef6f6a", "#8cc2ca", "#55ad89", "#c3bc3f", "#bb7693", "#baa094", "#a9b5ae", "#767676")

# 1. Rebuild g_rel and g_t from relativized_growth_rate_and_IPO_analysis_VOC_within_B117.R
UK = fread("./data/UK_lineage_rt_estimates_VOC_within_B117_Dec_12_stop_7_lagging_zeros.csv")

UK[lineage!='B.1.1.7',relativized_rt:=(growth_rate-mean(growth_rate,na.rm=T))/sd(growth_rate,na.rm=T),by=sample_date]
UK[lineage == "VOC", lineage := "VOC 202012/01"]

IPO_length <- 31
IPOs <- UK[,list(initial_rt=mean(relativized_rt[1:IPO_length],na.rm=T),
                 sd_rt=sd(relativized_rt[1:IPO_length],na.rm=T),
                 start_date=min(sample_date)),by=lineage]
top_lineages <- c('B.1.177','B.1.177.7','VOC 202012/01','B.1.1.1','B.1.1.301','B.1.1')
g_rel <- ggplot(IPOs[start_date<as.Date('2020-11-13')],aes(start_date,initial_rt))+
    annotate("ribbon", ymin = qnorm(0.025), ymax = qnorm(0.975), x = as.Date(c("2020-02-05", "2020-11-08")), 
        fill='grey', col='black', alpha=0.2, size = 0.25, linetype = "33") +
    geom_point(cex = 0.6)+
    geom_point(data=IPOs[lineage %in% top_lineages],aes(col=lineage),cex=3)+
    geom_point(data=IPOs[lineage=='VOC 202012/01'],aes(col=lineage),cex=5)+
    theme_bw(base_size=15)+
    geom_hline(yintercept = 0)+
    scale_x_date('Date of initial phylogenetic observation',
        breaks = as.Date(c('2020-03-01','2020-05-01','2020-07-01','2020-09-01','2020-11-01')),
        labels = c('March','May','July','Sept','Nov'))+
    scale_y_continuous('Mean relativized growth rate', limits = c(-3,3))+
    geom_smooth(data=IPOs[!lineage %in% top_lineages & start_date<as.Date('2020-11-13')],col='black',fill='black',alpha=0.2)+
    guides(col=guide_legend(ncol=6))+
    theme(legend.position = 'bottom')

UK[,days_since_IPO:=1:.N,by=lineage]
g_t <- ggplot(UK,aes(days_since_IPO,relativized_rt))+
    annotate("ribbon", ymin = qnorm(0.025), ymax = qnorm(0.975), x = c(0, 296), fill='grey', col='black', alpha=0.2, size = 0.25, linetype = "33") +
    geom_line(lwd=1,alpha=0.02,aes(group=lineage))+
    geom_line(data=UK[lineage %in% top_lineages],aes(color=lineage),lwd=1)+
    theme_bw(base_size=15)+
    geom_hline(yintercept = 0, size = 0.25)+
    geom_smooth(size = 0.5, colour = "black", se = FALSE)+
    scale_y_continuous('Relativized growth rate',limits=c(-5,5))+
    scale_x_continuous('Age of lineage (days)')+
    theme(legend.position='none')


# 2. Rebuild Mueller plot
mfit_preds = read.csv(file = "./multinomial_logistic_fits/tables/model 1a_model predictions Fig2C.csv")
mfit_preds$sample_date = as.Date(mfit_preds$sample_date_num, origin = "1970-01-01")
mfit_preds$variant = factor(mfit_preds$variant, levels = names(pal_variants))
levels_nhs_name = c("South East", "London", "East of England",
                    "South West", "Midlands", "North East and Yorkshire",
                    "Scotland", "North West", "Wales")
mfit_preds$nhs_name = factor(mfit_preds$nhs_name, levels = levels_nhs_name)

p_muller = ggplot(mfit_preds, aes(x = sample_date, y = prob, group = variant)) +
    facet_wrap(~nhs_name) +
    geom_area(aes(lwd = I(1.2), colour = NULL, fill = variant), position = "stack") +
    annotate("rect", xmin = as.Date("2021-01-07"), xmax = as.Date("2021-03-01"), 
        ymin = 0, ymax = 1, alpha = 0.4, fill = "white")

# 3. Load p_rt_sgtf from rt_sgtf_scatter.R
source("./rt_sgtf_scatter.R")


p_2 = cowplot::plot_grid(
    cowplot::plot_grid(
        g_rel + theme_cowplot(font_size = 10) + theme(legend.position = "bottom") + 
            labs(title = NULL, colour = NULL) + scale_colour_manual(values = pal_variants) + 
            guides(colour = guide_legend(nrow = 2, override.aes = list(size = 3))),
        g_t + theme_cowplot(font_size = 10) + theme(legend.position = "none") + labs(title = NULL, colour = NULL) + 
            scale_colour_manual(values = pal_variants),
        ncol = 1, labels = LETTERS, label_size = 10, rel_heights = c(5, 4)),
    p_muller + scale_x_date(date_breaks = "2 months", date_labels = "%b", limits = ymd(c("2020-02-05", "2021-03-01"))) + 
        theme_cowplot(font_size = 10) + theme(legend.position = "bottom", strip.background = element_blank(), 
            axis.text.x = element_text(angle = 90, vjust = 0.5)) +
        labs(y = "Relative abundance", x = "Sample date", fill = NULL) + scale_fill_manual(values = pal_variants),
    p_rt_sgtf,
    nrow = 1, rel_widths = c(7, 11, 8), labels = c("", "C", "D"), label_size = 10
)

ggsave("./output/fig-2.pdf", p_2, width = 36, height = 14, units = "cm", useDingbats = FALSE)
ggsave("./output/fig-2.png", p_2, width = 36, height = 14, units = "cm")
