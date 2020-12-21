library(data.table)
library(ggplot2)
library(zoo)
library(mgcv)
library(lubridate)
library(stringr)
library(cowplot)
library(qs)
library(ogwrangler)

theme_set(cowplot::theme_cowplot(font_size = 10) + theme(strip.background = element_blank()))

# Load Google Mobility data and interface with ogwrangler
gm = qread("./data/google_mobility_uk.qs");
CreateCache();
gm[, name := ifelse(sub_region_2 != "", paste0(sub_region_2, ", ", sub_region_1), sub_region_1)];
gm = gm[name != ""];
gm_match = data.table(code = ogcode("*", "gmcty"));
gm_match[, name := ogwhat(code, "name")];
gm_match[, rgn := ogwhat(code, "rgn")];
gm_match[, country := str_sub(rgn, 1, 1)];
gm_match[, rgn := NULL];
gm_match[, pop2019 := ogwhat(code, "pop2019")];
gm_match[, nhs := ogwhat(code, "nhser")];
gm_match[nhs %like% "^E", nhsname := ogwhat(nhs)];
gm = merge(gm, gm_match, by = "name");
gm = gm[, .SD, .SDcols = c(1, 16, 17, 18, 19, 9:15)];

# Format google mobility data
gm_melt = function(gm) {
    g = melt(gm, id.vars = 1:6)
    g[, variable := str_remove_all(variable, "_percent_change_from_baseline")]
    return (g)
}

#
# LOCKDOWN / CIRCUIT BREAKER ANALYSIS
#

# Northern Ireland
gm_ni = gm_melt(gm[country == "N"])
gm_ni = gm_ni[, .(value = weighted.mean(value, pop2019, na.rm = T)), keyby = .(variable, date, country)]
gm_ni[, rollval := rollmean(value, k = 7, fill = "extend"), by = variable]
gm_ni[, ld_start := ymd("2020-10-17")]; # start of first full day of lockdown
gm_ni[date >= ld_start & date <= "2020-10-30", ld_compdate := wday(date, week_start = 5) + ymd("2020-10-08")] # Compare to 9 to 15 October

# Plot NI google mobility data
pl_ni_mob = ggplot(gm_ni) +
    geom_line(aes(x = date, y = rollval, colour = variable)) +
    geom_vline(aes(xintercept = ymd("2020-10-16")), linetype = "22", size = 0.3) +
    facet_wrap(~variable) +
    theme(legend.position = "none") +
    labs(x = "Date", y = "Rolling mean - Northern Ireland")

# Compare post-lockdown values to pre-lockdown values
compare = gm_ni[!is.na(ld_compdate)]
compare = merge(compare,
    gm_ni[, .(ld_compdate = date, variable, baseline = value)],
    by = c("ld_compdate", "variable"))

compare = compare[order(date, variable)]
compare[, change := value - baseline]
compare[, .(mean = mean(change), se = sd(change)/sqrt(.N)), by = variable]
compare[, variable := str_replace_all(variable, "_", " ")]
compare[, variable := str_to_sentence(variable)]
compare[, variable := factor(variable)]

summaryN = compare[, .(mean = mean(change), se = sd(change)/sqrt(.N)), by = variable]

# IMPACT OF LOCKDOWN FOR NORTHERN IRELAND:
#                 variable        mean        se
# 1:  grocery_and_pharmacy  -1.3130269 0.9487168
# 2:                 parks   0.7009911 4.6486876
# 3:           residential   4.6109637 0.3439340
# 4: retail_and_recreation -14.9799870 1.4119360
# 5:      transit_stations -10.5919078 1.1598166
# 6:            workplaces -14.8809006 0.6219085

# Visual comparison
pl_ni = ggplot(compare) +
    geom_hline(aes(yintercept = 0), linetype = "23") + 
    geom_point(aes(x = date, y = change)) + 
    geom_smooth(aes(x = date, y = change), method = "lm", formula = y ~ 1) + 
    geom_label(data = summaryN, aes(x = ymd("2020-10-23") + 0.5, y = 15, 
        label = paste0(ifelse(mean > 0, "+", ""), round(mean, 2))), alpha = 0.75) +
    facet_wrap(~variable) +
    labs(x = "Date", y = "Change in Google Mobility index")





# Wales
gm_w = gm_melt(gm[country == "W"])
gm_w = gm_w[, .(value = weighted.mean(value, pop2019, na.rm = T)), keyby = .(variable, date, country)]
gm_w[, rollval := rollmean(value, k = 7, fill = "extend"), by = variable]
gm_w[, ld_start := ymd("2020-10-24")]; # start of first full day of lockdown
gm_w[date >= ld_start & date <= "2020-11-06", ld_compdate := wday(date, week_start = 5) + ymd("2020-10-15")] # Compare to 16 to 22 October

# Plot Wales google mobility data
pl_w_mob = ggplot(gm_w) + 
    geom_line(aes(x = date, y = rollval, colour = variable)) +
    geom_vline(aes(xintercept = ymd("2020-10-23")), linetype = "22", size = 0.3) +
    facet_wrap(~variable) +
    theme(legend.position = "none") +
    labs(x = "Date", y = "Rolling mean - Wales")

# Compare post-lockdown values to pre-lockdown values
compare = gm_w[!is.na(ld_compdate)]
compare = merge(compare,
    gm_w[, .(ld_compdate = date, variable, baseline = value)],
    by = c("ld_compdate", "variable"))

compare = compare[order(date, variable)]
compare[, change := value - baseline]
compare[, .(mean = mean(change), se = sd(change)/sqrt(.N)), by = variable]
compare[, variable := str_replace_all(variable, "_", " ")]
compare[, variable := str_to_sentence(variable)]
compare[, variable := factor(variable)]

summaryW = compare[, .(mean = mean(change), se = sd(change)/sqrt(.N)), by = variable]


# IMPACT OF LOCKDOWN FOR WALES:
#                 variable       mean       se
# 1:  grocery_and_pharmacy -16.131957 1.330679
# 2:                 parks -39.064812 6.721604
# 3:           residential   6.616149 0.452793
# 4: retail_and_recreation -37.521544 1.460017
# 5:      transit_stations -18.956221 1.154021
# 6:            workplaces -18.517148 1.565395

# Visual comparison
pl_w = ggplot(compare) +
    geom_hline(aes(yintercept = 0), linetype = "23") + 
    geom_point(aes(x = date, y = change)) + 
    geom_smooth(aes(x = date, y = change), method = "lm", formula = y ~ 1) + 
    geom_label(data = summaryW, aes(x = ymd("2020-10-30") + 0.5, y = -4, 
        label = paste0(ifelse(mean > 0, "+", ""), round(mean, 2))), alpha = 0.75) +
    facet_wrap(~variable) +
    labs(x = "Date", y = "Change in Google Mobility index")

cowplot::plot_grid(pl_ni_mob, pl_ni, pl_w_mob, pl_w, nrow = 2, labels = letters, label_size = 10)
ggsave("./figures/lockdown_analysis_2.png", width = 24, height = 18, units = "cm")
ggsave("./figures/lockdown_analysis_2.pdf", width = 24, height = 18, units = "cm", useDingbats = F)

cowplot::plot_grid(pl_ni, pl_w, nrow = 2, labels = letters, label_size = 10)

ggsave("./figures/lockdown_analysis.pdf", width = 18, height = 18, units = "cm", useDingbats = F)
ggsave("./figures/lockdown_analysis.png", width = 18, height = 18, units = "cm")




# Scotland - pub closures, coinciding with half term.
gm_s = gm_melt(gm[country == "S"])
gm_s = gm_s[, .(value = weighted.mean(value, pop2019, na.rm = T)), keyby = .(variable, date, country)]
gm_s[, rollval := rollmean(value, k = 7, fill = "extend"), by = variable]
gm_s[, ld_start := ymd("2020-10-10")]; # start of first full day of lockdown
gm_s[date >= ld_start & date <= "2020-10-23", ld_compdate := wday(date, week_start = 5) + ymd("2020-10-01")] # Compare to 2 to 8 October

# Plot Scotland google mobility data
ggplot(gm_s) + 
    geom_line(aes(x = date, y = rollmean(value, k = 7, fill = "extend"), colour = variable)) +
    geom_vline(aes(xintercept = ymd("2020-10-09"))) +
    geom_vline(aes(xintercept = ymd("2020-10-25"))) +
    facet_wrap(~variable) +
    theme(legend.position = "none")

# Compare post-lockdown values to pre-lockdown values
compare = gm_s[!is.na(ld_compdate)]
compare = merge(compare,
    gm_s[, .(ld_compdate = date, variable, baseline = value)],
    by = c("ld_compdate", "variable"))

compare = compare[order(date, variable)]
compare[, change := value - baseline]
compare[, .(mean = mean(change), se = sd(change)/sqrt(.N)), by = variable]
compare[, variable := str_replace_all(variable, "_", " ")]
compare[, variable := str_to_sentence(variable)]
compare[, variable := factor(variable)]

summaryS = compare[, .(mean = mean(change), se = sd(change)/sqrt(.N)), by = variable]

# IMPACT OF LOCKDOWN FOR SCOTLAND:
#                 variable       mean        se
# 1:  grocery_and_pharmacy  0.7392551 0.6692967
# 2:                 parks  6.3261325 7.3463237
# 3:           residential  1.9880707 0.4309069
# 4: retail_and_recreation -6.4052156 0.4742608
# 5:      transit_stations -5.7244764 1.1456064
# 6:            workplaces -8.4104283 1.1061221

# Visual comparison
pl_s = ggplot(compare) +
    geom_hline(aes(yintercept = 0), linetype = "23") + 
    geom_point(aes(x = date, y = change)) + 
    geom_smooth(aes(x = date, y = change), method = "lm", formula = y ~ 1) + 
    geom_label(data = summaryS, aes(x = ymd("2020-10-16") + 0.5, y = 28, 
        label = paste0(ifelse(mean > 0, "+", ""), round(mean, 2))), alpha = 0.75) +
    facet_wrap(~variable) +
    labs(x = "Date", y = "Change in Google Mobility index")




# NEW ADDITION 19 Nov 2020: England (Half term in most places, Sat 24 Oct to Sun 1 Nov)
gm_e_region = function(nhs_region)
{
    # England
    gm_e = gm_melt(gm[nhs == nhs_region])
    gm_e = gm_e[, .(value = weighted.mean(value, pop2019, na.rm = T)), keyby = .(variable, date, country)]
    gm_e[, rollval := rollmean(value, k = 7, fill = "extend"), by = variable]
    gm_e[, ld_start := ymd("2020-11-05")]; # start of first full day of lockdown
    gm_e[date >= ld_start & date <= "2020-12-01", ld_compdate := wday(date, week_start = 6) + ymd("2020-10-16")] # Compare to 17 to 23 October

    # Compare post-lockdown values to pre-lockdown values
    compare = gm_e[!is.na(ld_compdate)]
    compare = merge(compare,
        gm_e[, .(ld_compdate = date, variable, baseline = value)],
        by = c("ld_compdate", "variable"))
    
    compare = compare[order(date, variable)]
    compare[, change := value - baseline]
    compare[, variable := str_replace_all(variable, "_", " ")]
    compare[, variable := str_to_sentence(variable)]
    compare[, variable := factor(variable)]
    
    summaryE = compare[, .(mean = mean(change), se = sd(change)/sqrt(.N)), by = variable]
    
    cat(nhs_region, "\n")
    summaryE
}

reg = list()
reg[["East of England"]] =           gm_e_region("East of England")
reg[["London"]] =                    gm_e_region("London")
reg[["Midlands"]] =                  gm_e_region("Midlands")
reg[["North East and Yorkshire"]] =  gm_e_region("North East and Yorkshire")
reg[["North West"]] =                gm_e_region("North West")
reg[["South East"]] =                gm_e_region("South East")
reg[["South West"]] =                gm_e_region("South West")

regions = rbindlist(reg, idcol = "nhs")
regions[, .(mean = mean(mean), sd = sd(mean)), keyby = .(variable)]

# England
gm_e = gm_melt(gm[country == "E"])
gm_e = gm_e[, .(value = weighted.mean(value, pop2019, na.rm = T)), keyby = .(variable, date, country)]
gm_e[, rollval := rollmean(value, k = 7, fill = "extend"), by = variable]
gm_e[, ld_start := ymd("2020-11-05")]; # start of first full day of lockdown
gm_e[date >= ld_start & date <= "2020-12-01", ld_compdate := wday(date, week_start = 6) + ymd("2020-10-16")] # Compare to 17 to 23 October

# Plot England google mobility data
pl_e_mob = ggplot(gm_e) +
    geom_line(aes(x = date, y = rollval, colour = variable)) +
    geom_vline(aes(xintercept = ymd("2020-11-05")), linetype = "22", size = 0.3) +
    facet_wrap(~variable) +
    theme(legend.position = "none") +
    labs(x = "Date", y = "Rolling mean - England")

# Compare post-lockdown values to pre-lockdown values
compare = gm_e[!is.na(ld_compdate)]
compare = merge(compare,
    gm_e[, .(ld_compdate = date, variable, baseline = value)],
    by = c("ld_compdate", "variable"))

compare = compare[order(date, variable)]
compare[, change := value - baseline]
compare[, .(mean = mean(change), se = sd(change)/sqrt(.N)), by = variable]
compare[, variable := str_replace_all(variable, "_", " ")]
compare[, variable := str_to_sentence(variable)]
compare[, variable := factor(variable)]

summaryE = compare[, .(mean = mean(change), se = sd(change)/sqrt(.N)), by = variable]

# IMPACT OF LOCKDOWN FOR ENGLAND:
#                 variable       mean        se
# 1:  grocery_and_pharmacy  -4.782386 0.5578814
# 2:                 parks -19.335939 2.8611438
# 3:           residential   5.055462 0.2409448
# 4: retail_and_recreation -27.749227 0.8168795
# 5:      transit_stations -14.496543 1.1498847
# 6:            workplaces -10.295991 1.6405690

# Visual comparison
pl_e = ggplot(compare) +
    geom_hline(aes(yintercept = 0), linetype = "23") + 
    geom_point(aes(x = date, y = change)) + 
    geom_smooth(aes(x = date, y = change), method = "lm", formula = y ~ 1) + 
    geom_label(data = summaryE, aes(x = ymd("2020-11-14") + 0.5, y = -18, 
        label = paste0(ifelse(mean > 0, "+", ""), round(mean, 2))), alpha = 0.75) +
    facet_wrap(~variable) +
    labs(x = "Date", y = "Change in Google Mobility index")


# Plot showing overall effect of lockdown

ld = rbind(gm_e, gm_ni, gm_w)
ld = ld[date %between% list(ld_start - 21, ld_start + 28)]
ld[, date := as.Date(date)]
ld = ld[variable != "parks"]
ld[, variable := str_replace_all(variable, "_", " ")]
ld[, variable := str_to_sentence(variable)]
ld[, variable := factor(variable)]

ld[country == "E", country := "England"]
ld[country == "N", country := "Northern Ireland"]
ld[country == "W", country := "Wales"]

ld[, is_compare := date %in% ld_compdate, by = country]
ld[, is_lockdown := !is.na(ld_compdate)]

trend = ld[is_compare == TRUE | is_lockdown == TRUE, .(date, 
    comp_value = ifelse(is_compare, weighted.mean(value, is_compare), NA), 
    ld_value = ifelse(is_lockdown, weighted.mean(value, is_lockdown), NA),
    ld_start), by = .(country, variable)]

labels = trend[, 
    .(day_end = as.numeric(max(date[!is.na(ld_value)]) - ld_start), 
        value = unique(ld_value[!is.na(ld_value)]),
        change = unique(ld_value[!is.na(ld_value)]) - unique(comp_value[!is.na(comp_value)])), 
    by = .(country, variable)]

ggplot(ld) + 
    geom_line(aes(x = as.numeric(date - ld_start), y = value, colour = country), size = 0.2) +
    geom_line(data = trend, aes(x = as.numeric(date - ld_start), y = comp_value, colour = country), size = 1) +
    geom_line(data = trend, aes(x = as.numeric(date - ld_start), y = ld_value, colour = country), size = 1) +
    geom_vline(xintercept = 0, colour = "#666666", linetype = "22", size = 0.5) +
    geom_label(data = labels[!(variable == "Workplaces" & country %in% c("England", "Northern Ireland"))], aes(x = day_end, y = value, 
        label = paste0("(", ifelse(change > 0, "+", ""), round(change, 1), ")"), colour = country), 
        size = 3, hjust = 0, label.size = 0, label.padding = unit(0.1, "lines"), show.legend = FALSE) +
    geom_label(data = labels[variable == "Workplaces" & country == "England"], aes(x = day_end, y = value, 
        label = paste0("(", ifelse(change > 0, "+", ""), round(change, 1), ")"), colour = country), 
        size = 3, hjust = 0, vjust = 0, label.size = 0, label.padding = unit(0.1, "lines"), show.legend = FALSE) +
    geom_label(data = labels[variable == "Workplaces" & country == "Northern Ireland"], aes(x = day_end, y = value, 
        label = paste0("(", ifelse(change > 0, "+", ""), round(change, 1), ")"), colour = country), 
        size = 3, hjust = 0, vjust = 1.2, label.size = 0, label.padding = unit(0.1, "lines"), show.legend = FALSE) +
    facet_wrap(~variable, scales = "free") +
    theme(legend.position = c(0.76, 0.22)) +
    labs(x = "Days from lockdown", y = "Google Mobility index", colour = "Country") +
    scale_colour_manual(values = c("#bb0000", "#000000", "#0066ff")) +
    scale_x_continuous(breaks = c(-21, -14, -7, 0, 7, 14, 21, 28))

ggsave("./figures/lockdown_analysis_new.png", width = 24, height = 14, units = "cm")
ggsave("./figures/lockdown_analysis_new.pdf", width = 24, height = 14, units = "cm", useDingbats = F)

