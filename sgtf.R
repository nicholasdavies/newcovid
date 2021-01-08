ll = fread("~/Documents/uk_covid_data_sensitive/phe/20210106/Anonymised Combined Line List 20210106.csv")
sgtf = fread("~/Documents/uk_covid_data_sensitive/phe/20210106/SGTF_linelist_20210106.csv")
ll[, specimen_date := dmy(specimen_date)]
sgtf[, specimen_date := ymd(specimen_date)]
d = merge(ll, sgtf, by = c("FINALID", "specimen_date"), all = TRUE)

# IMD
ggplot(d[dmy(specimen_date) > "2020-10-01" & !is.na(imd_decile), 
    .(voc = sum(sgtf == 1, na.rm = T), other = sum(sgtf == 0, na.rm = T)), 
    keyby = .(date = dmy(specimen_date), imd_decile)]) + 
    geom_line(aes(x = date, y = voc / (voc + other), colour = imd_decile, group = imd_decile))

ggplot(d[dmy(specimen_date) > "2020-10-01" & !is.na(ethnicity_final), 
    .(voc = sum(sgtf == 1, na.rm = T), other = sum(sgtf == 0, na.rm = T)), 
    keyby = .(date = dmy(specimen_date), ethnicity_final)]) + 
    geom_line(aes(x = date, y = voc / (voc + other), colour = ethnicity_final, group = ethnicity_final))

ggplot(d[dmy(specimen_date) > "2020-10-01" & !is.na(sex), 
    .(voc = sum(sgtf == 1, na.rm = T), other = sum(sgtf == 0, na.rm = T)), 
    keyby = .(date = dmy(specimen_date), sex)]) + 
    geom_line(aes(x = date, y = voc / (voc + other), colour = sex, group = sex))

ggplot(d[dmy(specimen_date) > "2020-10-01" & !is.na(age), 
    .(voc = sum(sgtf == 1, na.rm = T), other = sum(sgtf == 0, na.rm = T)), 
    keyby = .(date = dmy(specimen_date), agegroup = cut(age, seq(0, 100, by = 10)))]) + 
    geom_line(aes(x = date, y = voc / (voc + other), colour = agegroup, group = agegroup)) +
    scale_y_continuous(trans = scales::logit_trans())

# by NHS region
dd = d[!is.na(specimen_date) & !is.na(NHSER_name) & NHSER_name != "", 
    .(sgtf = sum(sgtf == 1, na.rm = T), other = sum(sgtf == 0, na.rm = T)), 
    keyby = .(date = specimen_date, pillar, nhs_name = NHSER_name)]

ggplot(dd[date > "2020-10-01"]) +
    #geom_ribbon(data = raw[sample_date > "2020-10-01"], aes(x = sample_date, ymin = lower, ymax = upper), fill = "darkorchid", alpha = 0.4) +
    geom_line(aes(x = date, y = sgtf / (other + sgtf), colour = pillar)) + 
    facet_wrap(~nhs_name)

fwrite(dd[pillar == "Pillar 2" & date >= "2020-10-01", .(date, nhs_name, sgtf, other)], "./fitting_data/sgtf-2021-01-06.csv")

output = d[!is.na(specimen_date) & !is.na(NHSER_name) & NHSER_name != "" & pillar == "Pillar 2" & specimen_date >= "2020-09-01", 
    .(sgtf = sum(sgtf == 1, na.rm = T), other = sum(sgtf == 0, na.rm = T)), 
    keyby = .(date = specimen_date, nhser_name = NHSER_name, utla_name = UTLA_name, nhser_code = NHSER_code, utla_code = UTLA_code)]

saveRDS(output, "~/Desktop/sgtf.rds")

# Relative overrepresentation in under-20s
d[, age_group := ifelse(age < 20, "under_20", "20_plus")]
dage = d[pillar == "Pillar 2" & !is.na(sgtf) & !is.na(age_group) & !is.na(specimen_date) & NHSER_name != "", 
    .(under_20_sgtf = sum(age_group == "under_20" & sgtf == 1), over_20_sgtf = sum(age_group == "20_plus" & sgtf == 1),
      under_20_other = sum(age_group == "under_20" & sgtf == 0), over_20_other = sum(age_group == "20_plus" & sgtf == 0)), 
    keyby = .(specimen_date, NHSER_name)]

ggplot(dage) +
    geom_line(aes(x = specimen_date, y = log((under_20_sgtf/over_20_sgtf) / (under_20_other/over_20_other)))) +
    facet_wrap(~NHSER_name)

# Relative overrepresentation in over-70s
d[, age_group := ifelse(age < 70, "under_70", "70_plus")]
dage = d[pillar == "Pillar 2" & !is.na(sgtf) & !is.na(age_group) & !is.na(specimen_date) & NHSER_name != "", 
    .(under_70_sgtf = sum(age_group == "under_70" & sgtf == 1), over_70_sgtf = sum(age_group == "70_plus" & sgtf == 1),
      under_70_other = sum(age_group == "under_70" & sgtf == 0), over_70_other = sum(age_group == "70_plus" & sgtf == 0)), 
    keyby = .(specimen_date = floor_date(specimen_date, "week"), london = NHSER_name == "London")]

ggplot(dage) +
    geom_line(aes(x = specimen_date, y = log((over_70_sgtf/under_70_sgtf) / (over_70_other/under_70_other)), colour = london))


# by UTLA (LTLA also avail), pillar 2 only
dd = d[!is.na(specimen_date) & !is.na(NHSER_name) & UTLA_name != "" & pillar == "Pillar 2", 
    .(sgtf = sum(sgtf == 1, na.rm = T), other = sum(sgtf == 0, na.rm = T)), 
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

ggplot(dd[date >= "2020-10-01"]) +
    geom_line(aes(x = date, y = sgtf/(sgtf+other), colour = lat, group = lat)) +
    theme(legend.position = "none")

ggplot(dd[date >= "2020-10-01"]) +
    geom_smooth(aes(x = date, y = sgtf/(sgtf+other), colour = lat, group = lat), se = FALSE) +
    scale_y_continuous(trans = scales::logit_trans(), breaks = c(0.01, 0.1, 0.2, 0.5, 0.7, 0.9, 0.99), limits = c(0.01, 0.99)) +
    theme(legend.position = "none")

library(colourpal)

plotRaster = ggplot(dd[date >= "2020-10-01" & !is.na(sgtf) & sgtf + other > 0]) +
    geom_tile(aes(x = date, y = utla_name, fill = sgtf/(sgtf+other))) +
    theme_cowplot(font_size = 10) +
    theme(axis.text.y = element_text(size = unit(6, "pt")), panel.background = element_rect(fill = "#888888")) +
    labs(x = NULL, y = NULL, fill = "SGTF frequency") +
    scale_x_date(expand = expansion(0)) +
    ggpal("ocean.thermal")
ggsave("./output/utla_raster.pdf", plotRaster, width = 24, height = 30, units = "cm", useDingbats = FALSE)
ggsave("./output/utla_raster.png", plotRaster, width = 24, height = 30, units = "cm")


ggplot(dd[date >= "2020-10-01"]) +
    geom_smooth(aes(x = date, y = sgtf/(sgtf+other), colour = lat, group = lat), se = FALSE, method = "glm", method.args = list(family = "binomial")) +
    theme(legend.position = "none")

dd2 = dd[date >= "2020-11-02" & date < "2020-12-28"][order(date)]
dd2[, pop2019 := ogwhat(utla_code, "pop2019")]
dd2[, week := isoweek(date)]

dd3 = dd2[, .(sgtf = sum(sgtf), other = sum(other)), by = .(utla_code, utla_name, week, lat, long)]

plotXy = ggplot(dd3) +
    geom_line(aes(x = week, y = sgtf/(sgtf + other), group = utla_code, colour = lat)) +
    theme_cowplot(font_size = 10) +
    ggpal("ocean.thermal") +
    scale_y_continuous(trans = scales::logit_trans(), breaks = c(0.01, 0.1, 0.2, 0.5, 0.7, 0.9, 0.99), limits = c(0.01, 0.99)) +
    labs(x = "Week of 2020", y = "Relative frequency of S gene target failure", colour = "Latitude")

ggsave("./output/utla_xy.pdf", plotXy, width = 16, height = 12, units = "cm", useDingbats = FALSE)
ggsave("./output/utla_xy.png", plotXy, width = 16, height = 12, units = "cm")

raw
library(aod)
model = betabin(cbind(sgtf, other) ~ date + NHSER_name, ~ 1, data = dd[date >= "2020-10-01" & pillar == "Pillar 2"])

for (nhs in dd[, unique(nhs_name)]) {
    model = betabin(cbind(sgtf, other) ~ date, ~ 1, data = dd[date >= "2020-10-01" & pillar == "Pillar 2" & nhs_name == nhs])
    print(nhs)
    print(model@fixed.param[["date"]])
}

model
predicted2 = predict(model, newdata = raw[!nhs_name %in% c("Wales", "Northern Ireland", "Scotland")])
raw[!nhs_name %in% c("Wales", "Northern Ireland", "Scotland"), predicted2 := ..predicted2]

ggplot(d[!is.na(sgtf) & !is.na(imd_decile), .N, keyby = .(imd_decile, specimen_date, sgtf)]) +
    geom_line(aes(x = specimen_date, y = N, colour = as.factor(imd_decile), group = as.factor(imd_decile))) +
    facet_wrap(~as.factor(sgtf)) +
    scale_y_log10() +
    labs(colour = "IMD decile")

ggplot(d[!is.na(sgtf) & NHSER_name != "", .N, keyby = .(NHSER_name, specimen_date, sgtf)]) +
    geom_line(aes(x = specimen_date, y = N, colour = NHSER_name)) +
    facet_wrap(~as.factor(sgtf)) +
    scale_y_log10() +
    labs(colour = "NHS region")

ggplot(d[specimen_date > "2020-10-01" & NHSER_name != "" & pillar == "Pillar 2", .N, keyby = .(NHSER_name, specimen_date)]) +
    geom_line(aes(x = specimen_date, y = rollmean(N, 7, fill = NA), colour = NHSER_name)) +
    #scale_y_log10() +
    facet_wrap(~NHSER_name) +
    labs(colour = "NHS region") +
    scale_x_date(date_breaks = "1 week", date_label = "%b %d") +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))


d2 = copy(d[!is.na(sgtf)])
d2[, table(cat)]
d2[, age := (age %/% 5) * 5]
d2[, ethnicity_final := NULL]
d2

# delay to testing..?
delay = d[Onsetdate != "", .(delay = as.numeric(specimen_date - dmy(Onsetdate))), by = .(pillar, sgtf, asymptomatic_indicator)]
ggplot(delay[pillar == "Pillar 2"]) + geom_histogram(aes(x = delay, fill = asymptomatic_indicator))
