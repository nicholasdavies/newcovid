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

theme_set(cowplot::theme_cowplot(font_size = 8) + theme(strip.background = element_blank()))

cols_2 = c("#d575f0", "#6deda5")
cols_3 = c("#d575f0", "#6deda5", "#000000")

# Mappings ----------------------------------------------------------------

map_mobility <- c(
    "grocery_and_pharmacy" = "Grocery and pharmacy",
    "parks" = "Parks",
    "residential" = "Residential",
    "retail_and_recreation" = "Retail and recreation",
    "transit_stations" = "Transit stations",
    "workplaces" = "Workplaces"  
)



# Spread in ZAF and England -----------------------------------------------

.debug <- "data"
.args <- if (interactive()) sprintf(c(
    "%s/nextstrain_groups_ngs-sa_COVID19-ZA-2020.12.17_metadata.tsv",
    "%s/newlin.csv",
    "comparison.png"
), .debug) else commandArgs(trailingOnly = TRUE)

rsa.dt <- fread(.args[1])[Country == "South Africa"][order(`Collection Data`)]
rsa.dt[, area := "ZAF"]
uk.dt <- fread(.args[2])[order(sample_date)]
names(uk.dt)

# Remove rest of UK and split England by region
uk.dt <- uk.dt[!nhs_name %in% c("Northern Ireland", "Scotland", "Wales")]
uk.dt[nhs_name%in% c("East of England", "London", "South East"), area := "EOE, Lon, SE"]
uk.dt[!nhs_name %in% c("East of England", "London", "South East"), area := "Rest of England"]



uk.dt[, table(country, nhs_name)]
plot.dt <- rbind(
    rsa.dt[,
           .(.N, iso3c = "ZAF"),
           by = .(date=`Collection Data`, newvariant = `Clade`=="20C", area)
    ],
    uk.dt[,
          .(.N, iso3c="GBR"),
          by=.(date=sample_date, newvariant = var2, area)
    ]
)

plot.dt[, total := sum(N), by=.(date, iso3c, area) ]


plot.dt[, area := factor(area, levels = c("EOE, Lon, SE", "Rest of England", "ZAF"))]

plot.dt[,
        c("lo95","hi95") := 
            as.data.table(t(mapply(
                function(x, n, p=x/n) binom.test(x, n, p, conf.level = .95)$conf.int,
                x = N, n = total
            )))
][,
  c("lo50","hi50") := 
      as.data.table(t(mapply(
          function(x, n, p=x/n) binom.test(x, n, p, conf.level = .50)$conf.int,
          x = N, n = total
      )))
]

p_var <- plot.dt[newvariant == TRUE & area != "ZAF" & date <= "2020-12-01"] %>%
  ggplot(aes(x = date)) +
  geom_line(aes(y=N/total, col = area)) +
  geom_ribbon(aes(ymin = lo95, ymax = hi95, fill = area), alpha = 0.1) +
  geom_ribbon(aes(ymin = lo50, ymax = hi50, fill = area), alpha = 0.2) +
  annotate(geom = "rect", xmin = as.Date("2020-11-04"), xmax = as.Date("2020-12-01"), 
           ymin = -Inf, ymax = Inf, col = "lightgrey", alpha = 0.2) +
  facet_grid() +
  scale_x_date(
    name = NULL,
    date_breaks = "months", date_minor_breaks = "weeks",
    date_labels = "%b"
  ) +
  scale_color_manual(values = cols_2, name = "", labels = c("East of England,\nLondon and South\nEast regions",
                                                            "Rest of England")) +
  scale_fill_manual(values = cols_2, name = "", labels = c("East of England,\nLondon and South\nEast regions",
                                                           "Rest of England")) +
  scale_y_continuous("Novel variant proportion", expand = expansion(0)) +
  coord_cartesian(xlim = c(as.Date("2020-10-01"), as.Date("2020-12-01"))) + 
  theme(legend.position = c(0.01, 0.98)) + 
  labs(title = "A")

p_var

# Rt corr plot ------------------------------------------------------------

week_start <- 5

english_pillars <- readRDS("data/english_pillars.rds") %>%
  filter(between(date_specimen,
                 as.Date("2020-09-01"), as.Date("2020-12-10"))) %>%
  mutate(week_specimen = floor_date(date_specimen, "week",
                                    week_start = week_start))

rt_estimates <-
  paste0("https://raw.githubusercontent.com/epiforecasts/covid-rt-estimates/",
         "master/subnational/united-kingdom-local/cases/summary/rt.csv")

suppressMessages(rt <- vroom::vroom(rt_estimates))

by_ltla <- english_pillars %>%
  select(-lower_age_limit, -positive, -total) %>%
  group_by(pillar, week_specimen, nhser_name, ltla_name, ltla_code) %>%
  summarise_if(is.numeric, sum) %>%
  ungroup() %>%
  pivot_longer(starts_with("sgene"), names_to = "sgene_result",
               values_to = "n") %>%
  mutate(sgene_result = sub("^sgene_", "", sgene_result))


rt_by_ltla <- rt %>%
  rename(ltla_name = region) %>%
  filter(type == "estimate")


rt_weekly <- rt_by_ltla %>%
  mutate(week_infection = floor_date(date, "week", week_start = week_start)) %>%
  group_by(ltla_name, week_infection) %>%
  summarise(mean = mean(mean), sd = mean(sd), n = n(), .groups = "drop") %>%
  filter(n == 7) %>%
  select(-n)

by_ltla_rt <- by_ltla %>%
  filter(pillar == "Pillar 2") %>%
  mutate(week_infection = week_specimen - 7) %>%
  select(-pillar, -negative) %>%
  pivot_wider(names_from = sgene_result, values_from = n) %>%
  mutate(prop_variant = negative / (positive + negative),
         cases = n_a + negative + positive) %>%
  inner_join(rt_weekly, by = c("week_infection", "ltla_name")) %>%
  select(week_infection, nhser_name, ltla_name, ltla_code, prop_variant, cases,
         rt_mean = mean, rt_sd = sd)

p_corr <- by_ltla_rt %>% 
  filter(week_infection == max(week_infection)) %>%
  #mutate(cases = cases*0.01) %>%
  ggplot(aes(x = prop_variant, y = rt_mean, fill = nhser_name, size = cases)) +
  geom_jitter(pch = 21) +
  scale_size_continuous(range = c(0.5, 3)) +
  scale_fill_brewer("", palette = "Set1") +
  xlab("Proportion with S gene dropped") +
  ylab("Mean reproduction number") +
  labs(title = "B", 
       size = "Cases since 9 October")


p_corr

# Mobility data -----------------------------------------------------------
gmob <- qs::qread('data/gm_for_analysis.qs')
tier_data <- read.csv('data/england_ltla_covid_tiers_2020_12_20.csv')
tier_data <-  as.data.table(tier_data)
tier_data[,  table(RGN19NM, Tier == 4)]
gmob_tier <- merge(gmob, tier_data, by.x = "lad_nm", by.y = "LAD20NM")


gmob_tier[, variable := map_mobility[variable]]
gmob_tier[, t4  := fifelse(Tier == 4, "Tier 4", "Not Tier 4")]
gmob_tier[, t4  := factor(t4, levels = c("Tier 4", "Not Tier 4"), label = c("Enters Tier 4", "Outside of Tier 4"))]

gmob_tier[RGN19NM %in% c("South East", "East of England", "London") | t4 == TRUE, area := "South East\nEast of England,\nand London\nTier 4 areas*"]
gmob_tier[!RGN19NM %in% c("South East", "East of England", "London") | t4 == FALSE, area := "Rest of England"]

gmob_tier[, variable := factor(variable, levels = c(
                                "Residential",
                               "Workplaces" ,
                               "Grocery and pharmacy",
                               "Retail and recreation",
                               "Transit stations",
                               "Parks"))]

gmob_tier[, setattr(area, "levels", c("South East\nEast of England,\nand London\nTier 4 areas*", "Rest of England"))]

p_gmob <- gmob_tier[date > as.Date("2020-09-01")] %>% 
  ggplot(aes(x = date)) +
  geom_smooth(aes(y = value, col = area, fill = area)) +
  annotate(geom = "rect", xmin = as.Date("2020-11-04"), xmax = as.Date("2020-12-02"), ymin=-Inf, ymax=Inf, col = "lightgrey", alpha = 0.2) +
  facet_wrap(variable ~. , scales = "free_y", ncol = 1) +
  scale_x_date(breaks = "2 week", date_labels = "%d-%b", expand = expansion(0)) +
  scale_color_manual(values = cols_3, name = "", labels = c("East of England,\nLondon and South\nEast regions",
                                                            "Rest of England")) +
  scale_fill_manual(values = cols_3, name = "", labels = c("East of England,\nLondon and South\nEast regions",
                                                           "Rest of England")) +
  expand_limits(y = 0) +
  labs(title = "C", y = "Relative change in mobility (%)", x = "") +
  geom_hline(aes(yintercept = 0), linetype = 2) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    legend.key = element_rect(size = 6, fill = "white", colour = NA), legend.key.size = unit(1, "cm")
  ) + 
  labs(title = "B")
  

p_gmob
# Average contacts over time ----------------------------------------------


avg_contacts <- qs::qread("data/cmix_avg_cnts.qs")

avg_contacts[, area := factor(area,  levels = c("Tier 4", "Not Tier 4"), label = c("Entry to Tier 4", "Outside of Tier 4"))]

p_cmix <- avg_contacts[!age %in% c("All", "Adult")] %>%
  ggplot(aes(x = start_date)) +
  annotate(geom = "rect", xmin = as.Date("2020-11-04"), xmax = as.Date("2020-12-02"), ymin=-Inf, ymax=Inf, col = "lightgrey", alpha = 0.2) +
  expand_limits(y = 0) +
  geom_ribbon(aes(ymin = lci, ymax = uci, fill =area), alpha = 0.3) +
  geom_line( aes(y = med, col = area)) +
  #geom_smooth( aes(y = med), method = "gam") +
  #geom_point( aes(y = avg)) +
  labs(title = "B", y = "Mean contacts", x = "") +
  facet_grid(setting ~ age , scales = "free_y") +
  scale_y_continuous(expand = expansion(0)) +
  scale_color_manual(values = cols_3, name = "", labels = c("East of England,\nLondon and South\nEast\nTier 4 areas*",
                                                            "Rest of England")) +
  scale_fill_manual(values = cols_3, name = "", labels = c("East of England,\nLondon and South\nEast\nTier 4 areas*",
                                                           "Rest of England")) +
  theme(text = element_text(size = 8)) +
  scale_x_date(breaks = "2 week", date_labels = "%d-%b") +
  expand_limits(x = as.Date("2020-10-01")) + 
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    panel.spacing.y =  unit(1, "lines")) + 
  labs(title = "C")


# Input data ----------------------------------------------------------------

cmix_eng <- qs::qread("data/cmix_r0.qs")

# Compare to REACT --------------------------------------------------------
names(cmix_eng)
names(cmix_eng)
p_r0 <- ggplot(cmix_eng) +
    geom_rect(data = react_Rt, aes(ymin = `0.05`, ymax = `0.95`, xmin = time_from, xmax = creation_date, fill = estimate), alpha = 0.6) +
    geom_ribbon(aes(ymin = R_05, ymax = R_95, x = start_date), alpha = 0.2) +
    geom_ribbon(aes(ymin = R_25, ymax = R_75, x = start_date), alpha = 0.3) +
    geom_hline(aes(yintercept = 1), color='black', linetype=2) +
    geom_point(aes(x = start_date, y = R_50), fill = "red", col = "grey", pch =21) +
    scale_fill_manual(values = c("#3b70bf", "#7dba5d"), labels = c("Single round", "Two rounds")) +
    scale_y_continuous(expand = expansion(0), limits = c(0,3)) +
    labs(fill = "REACT R estimate", y = "R", x = "", title = "D")  +
    scale_x_date(breaks = "month", date_labels = "%b") +
    theme(
        legend.position = c(.05, .85),
        axis.text = element_text(size = 9),
        axis.title.y = element_text(size = 10),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10)
    )


# Compare CoMix and REACT -------------------------------------------------


setkeyv(cmix_eng,"start_date")
setkeyv(react_Rt,"time_from")

cmix_react <- react_Rt[,nearest:=(time_from)][cmix_eng,roll = 'nearest'] #Closest _previous_ date


## Both rounds
p_corr_plot <- ggplot(cmix_react, aes(R_50, `0.5`, col = estimate)) +
  geom_jitter() +
  scale_colour_manual(values = c("#3b70bf", "#7dba5d"), labels = c("Single round", "Two rounds")) +
  scale_x_continuous(expand = expansion(0), limits = c(0,2)) +
  scale_y_continuous(expand = expansion(0), limits = c(0,2)) +
  labs(colour = "REACT R estimate", y = "R", x = "", title = "E")  +
  geom_abline(aes(slope = 1, intercept = 0)) +
  ylab("REACT") +
  xlab("CoMix")

## corr.tests
cmix_react[estimate == "two_rounds", cor.test(R_50, `0.5`)]
cmix_react[estimate == "per_round", cor.test(R_50, `0.5`)]


# Load UTLA raster plot
library(qs)
p_raster = qread("./output/fig_raster.qs") + labs(title = "A") + theme(plot.title = element_text(size = 9))


# Combine plots -----------------------------------------------------------

layout <- "
    AAAAAABBBBCCCC
    AAAAAABBBBCCCC
    AAAAAABBBBCCCC
    AAAAAABBBBCCCC
    AAAAAABBBBCCCC
    AAAAAADDDDDEEE
    AAAAAADDDDDEEE
    AAAAAADDDDDEEE
    "

#plot_final <- p_var + p_corr + p_gmob + p_cmix + p_r0 + p_corr_plot + plot_layout(design = layout)
plot_final <- p_raster + p_gmob + p_cmix + p_r0 + p_corr_plot + plot_layout(design = layout)

ggsave(filename = "./output/new_figure_1.png", 
       plot_final,
       width = 20,
       height = 12)

ggsave(filename = "./output/new_figure_1.pdf", 
       plot_final,
       width = 20,
       height = 12)




