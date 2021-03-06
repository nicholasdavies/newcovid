# Packages ----------------------------------------------------------------
library(here)
library(vroom)
library(ggplot2)
library(lemon)
library(dplyr)
library(cowplot)
library(lubridate)
library(stringr)
# devtools::install_github("epiforecasts/covidregionaldata")
library(covidregionaldata)
# Get data ----------------------------------------------------------------
utla_rt <- readRDS(url("https://raw.github.com/epiforecasts/covid19.sgene.utla.rt/main/data/utla_rt_with_covariates.rds"))
# link with cases
cases <- get_regional_data("UK", include_level_2_regions = TRUE) %>%
  filter(date >= min(utla_rt$week_infection) + 7) %>%
  mutate(week_infection = floor_date(date, "week", 1) - 7) %>%
  group_by(week_infection, utla_name = authority) %>%
  summarise(cases = sum(cases_new, na.rm = TRUE), .groups = "drop")
# join rt and cases
utla_rt_cases <- utla_rt %>%
  left_join(cases, by = c("week_infection", "utla_name")) %>%
  filter(week_infection >= "2020-10-01", !is.na(prop_sgtf))
# extract just Rt estimates made with a 5.5 day mean generation time
utla_rt  <- utla_rt_cases %>%
  rename_with(~ sub(paste0("_long_gt"), "", .x)) %>%
  filter(!is.na(rt_mean))
# filter for target dates
target_rt <- utla_rt %>%
  filter(as.character(week_infection) %in% c("2020-10-26", "2020-11-09",
                                             "2020-11-23", "2020-12-07",
                                             "2020-12-21"))

target_rt$week_infection = str_replace_all(format(target_rt$week_infection, "%b %e %Y"), "  ", " ")
target_rt$week_infection = factor(target_rt$week_infection, levels = unique(target_rt$week_infection))
target_rt$nhser_name[target_rt$nhser_name == "North East and Yorkshire"] = "North East\nand Yorkshire"

p_rt_sgtf <- ggplot(target_rt, aes(x = prop_sgtf, y = rt_mean,
                              fill = nhser_name, size = cases)) +
  geom_jitter(pch = 21) +
  facet_rep_wrap(. ~ week_infection, ncol = 2, repeat.tick.labels = "all") +
  scale_fill_brewer("", palette = "Set1") +
  xlab("Frequency of S gene target failure") +
  ylab("Mean reproduction number") +
  theme_cowplot(font_size = 10) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(size = paste("Cases")) +
  guides(fill = guide_legend(title = "NHS region", ncol = 1),
         size = guide_legend(ncol = 1)) +
  theme(legend.position = c(0.5, 0.15),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.key.height = unit(0.45, "cm"),
        strip.background = element_blank())
p_rt_sgtf
