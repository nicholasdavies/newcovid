## Mobility, contacts, R and React. 

# Packages ----------------------------------------------------------------
library(data.table)
library(ggplot2)
library(patchwork)

# Source user written scripts ---------------------------------------------
source('get_react_data.R')

theme_set(cowplot::theme_cowplot(font_size = 8) + theme(strip.background = element_blank()))

cols = c("#8EC3F3", "#7EA59E")

# Mobility data -----------------------------------------------------------

gmob <- qs::qread('data/gm_for_analysis.qs')
tier_data <- read.csv('data/england_ltla_covid_tiers_2020_12_20.csv')

gmob_tier <- merge(gmob, tier_data, by.x = "lad_nm", by.y = "LAD20NM")

gmob_tier[, t4  := fifelse(Tier == 4, "Tier 4", "Not Tier 4")]
gmob_tier[, t4  := factor(t4, levels = c("Tier 4", "Not Tier 4"), label = c("Enters Tier 4", "Outside of Tier 4"))]

plot_a <- ggplot(gmob_tier[date > as.Date("2020-09-01")]) +
    geom_rect(aes(xmin = as.Date("2020-11-04"), xmax = as.Date("2020-12-05"), ymin=-Inf, ymax=Inf), col = "lightgrey", alpha = 0.01) +
    geom_vline(aes(xintercept = as.Date("2020-11-04")), linetype = 2) +
    geom_vline(aes(xintercept = as.Date("2020-12-05")), linetype = 2) +
    geom_smooth(aes(x = date, y = value, col = t4, fill = t4)) +
    facet_wrap(variable ~. , scales = "free_y", ncol = 1) +
    scale_x_date(breaks = "2 week", date_labels = "%d-%b", expand = expansion(0)) +
    scale_color_manual(values = cols, name = "") +
    scale_fill_manual(values = cols, name = "") +
    expand_limits(y = 0) +
    labs(title = "A", y = "% change in mobility", x = "") +
    geom_hline(aes(yintercept = 0), linetype = 2) +
    theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
    )

# Average contacts over time ----------------------------------------------


avg_cnts <- qs::qread("data/cmix_avg_cnts.qs")

avg_cnts[, area := factor(area,  levels = c("Tier 4", "Not Tier 4"), label = c("Entry to Tier 4", "Outside of Tier 4"))]

plot_b <-  
    ggplot(avg_cnts[!age %in% c("All", "Adult")], aes(x = start_date)) +
    geom_rect(aes(xmin = as.Date("2020-11-04"), xmax = as.Date("2020-12-05"), ymin=-Inf, ymax=Inf), col = "lightgrey", alpha = 0.01) +
    geom_vline(aes(xintercept = as.Date("2020-11-04")), linetype = 2) +
    geom_vline(aes(xintercept = as.Date("2020-12-05")), linetype = 2) +
    expand_limits(y = 0) +
    geom_ribbon(aes(ymin = lci, ymax = uci, fill =area), alpha = 0.3) +
    geom_line( aes(y = med, col = area)) +
    #geom_smooth( aes(y = med), method = "gam") +
    #geom_point( aes(y = avg)) +
    labs(title = "B", y = "Mean contacts", x = "") +
    facet_grid(setting ~ age , scales = "free_y") +
    scale_y_continuous(expand = expansion(0)) +
    scale_fill_manual(values = cols, name = "") +
    scale_color_manual(values = cols, name = "") +
    theme(text = element_text(size = 12)) +
    scale_x_date(breaks = "2 week", date_labels = "%d-%b") +
    expand_limits(x = as.Date("2020-10-01")) + 
    theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.spacing.y =  unit(1, "lines"), 
        
    )


# Input data ----------------------------------------------------------------

cmix_eng <- qs::qread("data/cmix_r0.qs")





# Compare to REACT --------------------------------------------------------
names(cmix_eng)
plot_c <- ggplot(cmix_eng) +
    geom_rect(data = react_Rt, aes(ymin = `0.05`, ymax = `0.95`, xmin = time_from, xmax = creation_date, fill = estimate), alpha = 0.9, col = "black") +
    geom_ribbon(aes(ymin = R_05, ymax = R_95, x = start_date), alpha = 0.2, color='DarkGrey') +
    geom_ribbon(aes(ymin = R_25, ymax = R_75, x = start_date), alpha = 0.3) +
    geom_hline(aes(yintercept = 1), color='black', linetype=2) +
    geom_point(aes(x = start_date, y = R_50), fill = "red", col = "black", pch =21) +
    scale_fill_manual(values = c("#8EC3F3", "#7EA59E"), labels = c("Single round", "Two rounds")) +
    scale_y_continuous(expand = expansion(0), limits = c(0,3)) +
    labs(fill = "REACT R estimate", y = "R", x = "", title = "C")  +
    scale_x_date(breaks = "month", date_labels = "%b") +
    theme(legend.position = c(0.8, 0.9),
          axis.text = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          strip.text = element_text(size = 12),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10)
    )


plot_ab <- (plot_a | plot_b) + plot_layout(guides = "collect")

plot_abc <- plot_ab / plot_c
plot_abc


#ggsave(filename = "outputs/Figure_1.png", plot_abc, width = 20, height = 12)
