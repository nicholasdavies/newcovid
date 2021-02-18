## Mobility, contacts, R and React. 

# Packages ----------------------------------------------------------------
library(data.table)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(lubridate)

# Source user written scripts ---------------------------------------------
source('get_react_data.R')

theme_set(cowplot::theme_cowplot(font_size = 10) + theme(strip.background = element_blank()))

cols_2 = c("#d575f0", "#6deda5")
cols_3 = c("#c83737", "#37abc8", "#000000")

# Mappings ----------------------------------------------------------------

map_mobility <- c(
    "grocery_and_pharmacy" = "Grocery and pharmacy",
    "parks" = "Parks",
    "residential" = "Residential",
    "retail_and_recreation" = "Retail and recreation",
    "transit_stations" = "Transit stations",
    "workplaces" = "Workplaces"  
)

# Mobility data -----------------------------------------------------------
gmob <- qs::qread('data/gm_for_analysis.qs')
tier_data <- read.csv('data/england_ltla_covid_tiers_2020_12_20.csv')
tier_data <-  as.data.table(tier_data)
tier_data[,  table(RGN19NM, Tier == 4)]
gmob_tier <- merge(gmob, tier_data, by.x = "lad_nm", by.y = "LAD20NM")


gmob_tier[, variable := map_mobility[variable]]
gmob_tier[, t4  := fifelse(Tier == 4, "Tier 4", "Not Tier 4")]
gmob_tier[, t4  := factor(t4, levels = c("Tier 4", "Not Tier 4"), label = c("Enters Tier 4", "Outside of Tier 4"))]

gmob_tier[RGN19NM %in% c("South East", "East of England", "London") | t4 == TRUE, area := "Tier 4 areas"]
gmob_tier[!RGN19NM %in% c("South East", "East of England", "London") | t4 == FALSE, area := "Rest of England"]

gmob_tier[, variable := factor(variable, levels = c(
                                "Residential",
                               "Workplaces" ,
                               "Grocery and pharmacy",
                               "Retail and recreation",
                               "Transit stations",
                               "Parks"))]

gmob_tier[, setattr(area, "levels", c("Tier 4 areas", "Rest of England"))]

p_gmob <- gmob_tier[date > as.Date("2020-09-01")] %>% 
  ggplot(aes(x = date)) +
  geom_smooth(aes(y = value, col = area, fill = area)) +
  annotate(geom = "rect", xmin = as.Date("2020-11-04"), xmax = as.Date("2020-12-02"), ymin=-Inf, ymax=Inf, col = "lightgrey", alpha = 0.2) +
  facet_wrap(variable ~. , scales = "free_y", ncol = 2) +
  scale_x_date(breaks = "2 week", date_labels = "%e %b", expand = expansion(0)) +
  scale_color_manual(values = cols_3, name = "", labels = c("East of England,\nLondon and South\nEast regions",
                                                            "Rest of England")) +
  scale_fill_manual(values = cols_3, name = "", labels = c("East of England,\nLondon and South\nEast regions",
                                                           "Rest of England")) +
  expand_limits(y = 0) +
  labs(y = "Relative change in mobility (%)", x = NULL) +
  geom_hline(aes(yintercept = 0), linetype = 2) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust=1),
    panel.spacing.y =  unit(0.25, "lines"),
    legend.key = element_rect(size = 6, fill = "white", colour = NA), legend.key.size = unit(1, "cm")
  )
  

p_gmob


# Average contacts over time ----------------------------------------------
avg_contacts <- qs::qread("data/cmix_avg_cnts.qs")

avg_contacts[, area := factor(area,  levels = c("Tier 4", "Not Tier 4"), label = c("Entry to Tier 4", "Outside of Tier 4"))]

avg_contacts[, setting2 := as.character(setting)]
avg_contacts[setting2 == "Work/Educ", setting2 := "Work / Education"]
avg_contacts[, table(setting2)]
avg_contacts[, setting2 := factor(setting2, levels = c("All", "Home", "Work / Education", "Other"))]
avg_contacts[, table(setting2)]

p_cmix <- avg_contacts[!age %in% c("All", "Adult") & !setting == "All"] %>%
  ggplot(aes(x = start_date)) +
  annotate(geom = "rect", xmin = as.Date("2020-11-04"), xmax = as.Date("2020-12-02"), ymin=-Inf, ymax=Inf, col = "lightgrey", alpha = 0.2) +
  expand_limits(y = 0) +
  geom_ribbon(aes(ymin = lci, ymax = uci, fill =area), alpha = 0.3) +
  geom_line( aes(y = med, col = area)) +
  labs(y = "Mean daily contacts", x = NULL) +
  facet_grid(setting2 ~ age , scales = "free_y", switch = "y") +
  scale_y_continuous(expand = expansion(0)) +
  scale_color_manual(values = cols_3, name = "", labels = c("Tier 4 areas",
                                                            "Rest of England")) +
  scale_fill_manual(values = cols_3, name = "", labels = c("Tier 4 areas",
                                                           "Rest of England")) +
  scale_x_date(breaks = "2 week", date_labels = "%e %b") +
  expand_limits(x = as.Date("2020-10-01")) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust=1),
    panel.spacing.y =  unit(0.4, "lines"),
    legend.position = c(0.68, 0.97),
    strip.placement = "outside")

p_cmix

# Input data ----------------------------------------------------------------

cmix_eng <- qs::qread("data/cmix_r0.qs")

# Compare to REACT --------------------------------------------------------
# p_r0 <- ggplot(cmix_eng) +
#     geom_rect(data = react_Rt, aes(ymin = `0.05`, ymax = `0.95`, xmin = time_from, xmax = creation_date, fill = estimate), alpha = 0.6) +
#     geom_ribbon(aes(ymin = R_05, ymax = R_95, x = start_date), alpha = 0.2) +
#     geom_ribbon(aes(ymin = R_25, ymax = R_75, x = start_date), alpha = 0.3) +
#     geom_hline(aes(yintercept = 1), color='black', linetype=2) +
#     geom_point(aes(x = start_date, y = R_50), fill = "red", col = "grey", pch =21) +
#     scale_fill_manual(values = c("#3b70bf", "#7dba5d"), labels = c("Single round", "Two rounds")) +
#     scale_y_continuous(expand = expansion(0), limits = c(0,3)) +
#     labs(fill = "REACT R estimate", y = "R", x = "")  +
#     scale_x_date(breaks = "month", date_labels = "%b") +
#     theme(
#         legend.position = c(.05, .85),
#         axis.text = element_text(size = 9),
#         axis.title.y = element_text(size = 10),
#         strip.text = element_text(size = 12),
#         legend.text = element_text(size = 10),
#         legend.title = element_text(size = 10)
#     )

p_r0 <- ggplot(cmix_eng) +
    geom_hline(aes(yintercept = 1), color='black', linetype=2) +
    geom_ribbon(aes(ymin = R_05, ymax = R_95, x = start_date), alpha = 0.2) +
    geom_line(aes(x = start_date, y = R_50), col = "black") +
    geom_point(data = react_Rt, aes(y = `0.5`, x = time_from + (creation_date - time_from)/2, colour = estimate)) +
    geom_linerange(data = react_Rt, aes(y = `0.5`, xmin = time_from, xmax = creation_date, colour = estimate)) +
    geom_linerange(data = react_Rt, aes(ymin = `0.05`, ymax = `0.95`, x = time_from + (creation_date - time_from)/2, colour = estimate)) +
    scale_colour_manual(values = c("#3b70bf", "#3dba1d"), labels = c("Single round", "Two rounds")) +
    scale_y_continuous(expand = expansion(0), limits = c(0,3)) +
    labs(colour = "REACT R estimate", y = "R", x = "")  +
    scale_x_date(breaks = "month", date_labels = "%b") +
    theme(
        legend.position = c(.05, .85),
        axis.text = element_text(size = 9),
        axis.title.y = element_text(size = 10),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10)
    )
p_r0

# Compare CoMix and REACT -------------------------------------------------


setkeyv(cmix_eng,"start_date")
setkeyv(react_Rt,"time_from")

cmix_react <- react_Rt[,nearest:=(time_from)][cmix_eng,roll = 'nearest'] #Closest _previous_ date


## Both rounds
# p_corr_plot <- ggplot(cmix_react, aes(R_50, `0.5`, col = estimate)) +
#   geom_jitter() +
#   scale_colour_manual(values = c("#3b70bf", "#7dba5d"), labels = c("Single round", "Two rounds")) +
#   scale_x_continuous(expand = expansion(0), limits = c(0,2)) +
#   scale_y_continuous(expand = expansion(0), limits = c(0,2)) +
#   labs(colour = "REACT R estimate", y = "R", x = "")  +
#   geom_abline(aes(slope = 1, intercept = 0)) +
#   ylab("REACT") +
#   xlab("CoMix")
# p_corr_plot

cmix_react2 = cmix_react[, .(`0.5` = mean(`0.5`), `0.05` = mean(`0.05`), `0.95` = mean(`0.95`),
    R_50 = mean(R_50), R_05 = mean(R_05), R_95 = mean(R_95)),
    by = .(estimate, nearest)]

p_corr_plot = ggplot(cmix_react2, aes(R_50, `0.5`, col = estimate)) +
    geom_abline(aes(slope = 1, intercept = 0), linetype = "33") + 
    geom_point() +
    geom_linerange(aes(x = R_50, y = `0.5`, ymin = `0.05`, ymax = `0.95`, colour = estimate)) +
    geom_linerange(aes(x = R_50, y = `0.5`, xmin = R_05, xmax = R_95, colour = estimate)) +
    scale_colour_manual(values = c("#3b70bf", "#3dba1d"), labels = c("Single round", "Two rounds")) +
    labs(colour = "REACT R estimate", y = "REACT R estimate", x = "CoMix R estimate")  +
    scale_x_log10(expand = expansion(0), limits = c(0.35, 2.3), breaks = c(0.5, 0.7, 1, 1.4, 2)) + 
    scale_y_log10(expand = expansion(0), limits = c(0.35, 2.3), breaks = c(0.5, 0.7, 1, 1.4, 2)) +
    theme(legend.position = "none")
p_corr_plot

## corr.tests
cmix_react[estimate == "two_rounds", cor.test(R_50, `0.5`)]
cmix_react[estimate == "per_round", cor.test(R_50, `0.5`)]

# Load UTLA raster plot
library(qs)
p_raster = qread("./output/fig_raster.qs")

# Load demographics plots
p_demo1 = qread("./output/plot_demographics_london_1.qs")
p_demo2 = qread("./output/plot_demographics_london_2.qs")
p_demo3 = qread("./output/plot_demographics_london_3.qs")

periods = stringr::str_trim(format(as.Date(as.character(p_demo3$data$time)), "%e %b %Y"));
p_demo3$data$time = factor(periods, levels = rev(unique(periods)))

# Combine plots -----------------------------------------------------------

plot_final = plot_grid(p_raster,
    plot_grid(
        plot_grid(p_demo1, p_demo2, p_demo3, nrow = 1, labels = c("B", "", ""), label_size = 10, rel_widths = c(10, 10, 6), align = "hv", axis = "bottom"),
        plot_grid(p_r0, p_corr_plot, nrow = 1, labels = c("C", "D"), label_size = 10, rel_widths = c(6.5, 3)),
        plot_grid(p_gmob, p_cmix, nrow = 1, labels = c("E", "F"), label_size = 10, rel_widths = c(2, 3)),
    nrow = 3, rel_heights = c(2, 2.1, 3)), 
    nrow = 1, labels = c("A", ""), label_size = 10, rel_widths = c(6, 8)
)

ggsave(filename = "./output/new_figure_1.png", 
       plot_final,
       width = 16,
       height = 9.6)

ggsave(filename = "./output/new_figure_1.pdf", 
       plot_final,
       width = 16,
       height = 9.6)




