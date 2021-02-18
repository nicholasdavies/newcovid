library(data.table)
library(ggplot2)
library(cowplot)
library(colourpal)
library(zoo)
library(lubridate)

ll = fread("~/Documents/uk_covid_data_sensitive/phe/20210118/Anonymised Combined Line List 20210118.csv")
sgtf = fread("~/Documents/uk_covid_data_sensitive/phe/20210118/SGTF_linelist_20210118.csv")
ll[, specimen_date := dmy(specimen_date)]
sgtf[, specimen_date := ymd(specimen_date)]
d = merge(ll, sgtf, by = c("FINALID", "specimen_date"), all = TRUE)

d[!is.na(age), agegroup := paste0(pmin((age %/% 10) * 10, 80), "-", pmin((age %/% 10) * 10, 80) + 9)]
d[agegroup == "80-89", agegroup := "80+"]

plotby = function(d, what)
{
    d2 = d[, .SD, .SDcols = c("specimen_date", "sgtf_under30CT", what)]
    setnames(d2, what, "what")
    d2 = d2[specimen_date >= "2020-11-01" & !is.na(what) & what != "",
        .(voc = sum(sgtf_under30CT == 1, na.rm = T), other = sum(sgtf_under30CT == 0, na.rm = T)), 
        keyby = .(date = specimen_date, what)];
    d2[, f_voc := rollmean(voc / (voc + other), 7, fill = NA), by = what]
    ggplot(d2) +
        geom_line(aes(x = date, y = f_voc, colour = what, group = what)) +
        scale_y_continuous(trans = scales::logit_trans(), limits = c(0.01, 0.99))
}

library(pals)

plotby2 = function(d, what, nhs_name = NULL, ylim = c(0.01, 0.85))
{
    d2 = copy(d)
    d2[, what := get(..what)]
    if (!is.null(nhs_name)) {
        d2 = d2[NHSER_name == nhs_name]
    }
    d2 = d2[!is.na(sgtf_under30CT) & !is.na(what) & what != "" & what != "Unknown" & specimen_date >= "2020-10-26", 
        .(sg = sum(sgtf_under30CT), oth = sum(1 - sgtf_under30CT)), 
        keyby = .(time = (as.numeric(specimen_date - ymd("2020-10-26")) %/% 14) * 14 + ymd("2020-10-26") , what)]
    
    breaks = c(0.01, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 0.99)
    breaks = breaks[breaks >= ylim[1] & breaks <= ylim[2]]
    
    d2[, time := factor(as.character(time), levels = rev(as.character(unique(time))))]
    
    ggplot(d2) +
        geom_line(aes(x = what, y = sg / (sg + oth), colour = time, group = time)) +
        scale_y_continuous(trans = scales::logit_trans(), limits = ylim, breaks = breaks)
}

plotby(d, "imd_decile")
plotby(d, "agegroup")
plotby(d, "NHSER_name")
plotby(d, "ethnicity_final")
plotby(d, "sex")
plotby(d, "UTLA_name") + theme(legend.position = "none")


# DEMOGRAPHICS PLOTS
pl1 = plotby2(d, "agegroup", "London", c(0.01, 0.98)) + labs(x = NULL, y = "Frequency of\nS gene target failure", colour = "Time") + 
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
pl2 = plotby2(d, "imd_decile", "London", c(0.01, 0.98)) + scale_x_continuous(breaks = 1:10) + 
    labs(x = "Index of multiple deprivation decile", y = NULL, colour = "Time") + theme(legend.position = "none")
pl3 = plotby2(d, "sex", "London", c(0.01, 0.98)) + labs(x = NULL, y = NULL, colour = "Period starting") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

qsave(pl1, "./output/plot_demographics_london_1.qs")
qsave(pl2, "./output/plot_demographics_london_2.qs")
qsave(pl3, "./output/plot_demographics_london_3.qs")


# MAKE OUTPUT
dd = d[!is.na(specimen_date) & !is.na(NHSER_name) & NHSER_name != "", 
    .(sgtf = sum(sgtf_under30CT == 1, na.rm = T), other = sum(sgtf_under30CT == 0, na.rm = T)), 
    keyby = .(date = specimen_date, pillar, nhs_name = NHSER_name)]

ggplot() +
    geom_line(data = dd[date > "2020-10-01" & pillar == "Pillar 2"], aes(x = date, y = sgtf / (other + sgtf))) + 
    facet_wrap(~nhs_name)

fwrite(dd[pillar == "Pillar 2" & date >= "2020-10-01", .(date, nhs_name, sgtf, other)], "./fitting_data/sgtf-2021-01-18.csv")


# RASTER PLOT
ll = fread("~/Documents/uk_covid_data_sensitive/phe/20210208/Anonymised Combined Line List 20210208.csv")
sgtf = fread("~/Documents/uk_covid_data_sensitive/phe/20210208/SGTF_linelist_20210208.csv")
ll[, specimen_date := dmy(specimen_date)]
sgtf[, specimen_date := ymd(specimen_date)]
d = merge(ll, sgtf, by = c("FINALID", "specimen_date"), all = TRUE)

# by UTLA (LTLA also avail), pillar 2 only
dd = d[!is.na(specimen_date) & !is.na(NHSER_name) & UTLA_name != "" & pillar == "Pillar 2", 
    .(sgtf = sum(sgtf_under30CT == 1, na.rm = T), other = sum(sgtf_under30CT == 0, na.rm = T)), 
    keyby = .(date = specimen_date, utla_name = UTLA_name, utla_code = UTLA_code)]

# Replace old Buckinghamshire county with new Buckinghamshire local authority
dd[utla_code == "E10000002", utla_code := "E06000060"]

# ONS Postcode Directory
onspd = fread("~/Downloads/ONSPD_AUG_2020_UK/Data/ONSPD_AUG_2020_UK.csv");
onspd = onspd[is.na(doterm), .(oslaua, oscty, lat, long)]
onspd = rbind(
    onspd[, .(lat = mean(lat, na.rm = T), long = mean(long, na.rm = T)), by = .(code = oslaua)],
    onspd[, .(lat = mean(lat, na.rm = T), long = mean(long, na.rm = T)), by = .(code = oscty)]
)
dd = merge(dd, onspd, by.x = "utla_code", by.y = "code", all = TRUE)
dd = dd[utla_code %like% "^E"]
dd = dd[!utla_code %like% "^E99"]

dd[, lat_rank := frank(lat, ties.method = "dense")]
dd = dd[order(lat_rank)]
dd[, utla_name := factor(utla_name, levels = unique(utla_name))]

library(ogwrangler)
CreateCache()
dd[, nhser := ogwhat(utla_code, "nhser")]
dd[, nhser_name := ogwhat(nhser)]

library(colourpal)

plotRaster = ggplot(dd[date >= "2020-10-15" & date <= "2021-02-01" & !is.na(sgtf) & sgtf + other > 0]) +
    geom_tile(aes(x = date, y = utla_name, fill = sgtf/(sgtf+other))) +
    theme_cowplot(font_size = 10) +
    theme(axis.text.y = element_text(size = unit(5, "pt")), panel.background = element_rect(fill = "#888888"), 
        legend.position = c(-0.24, 0.94), legend.background = element_rect(fill = "#ffffff"), legend.margin = margin(2,2,5,2)) +
    labs(x = NULL, y = NULL, fill = NULL) +
    scale_x_date(expand = expansion(0)) +
    ggpal("ocean.thermal")

plotRaster
ggsave("./output/utla_raster.pdf", plotRaster, width = 24, height = 30, units = "cm", useDingbats = FALSE)
ggsave("./output/utla_raster.png", plotRaster, width = 24, height = 30, units = "cm")
qsave(plotRaster, "./output/fig_raster.qs")

