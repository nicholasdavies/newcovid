library(mgcv)
library(data.table)
library(readxl)
library(lubridate)
library(DescTools)

prefer = function(a, b)
{
    ifelse (!is.na(a), a, b)
}

# Load and assemble data
# Note: These files contain personally identifiable information, so they are not included with the repo.
d_death = read_excel("~/Documents/uk_covid_data_sensitive/phe/20210106/20210106 COVID19 Deaths.xlsx")
ll = fread("~/Documents/uk_covid_data_sensitive/phe/20210106/Anonymised Combined Line List 20210106.csv")
sgtf = fread("~/Documents/uk_covid_data_sensitive/phe/20210106/SGTF_linelist_20210106.csv")

setDT(d_death)
ll[, specimen_date := dmy(specimen_date)]
sgtf[, specimen_date := ymd(specimen_date)]

d = merge(ll, sgtf, by = c("FINALID", "specimen_date"), all = TRUE)
d = merge(d, d_death, by.x = "FINALID", by.y = "finalid", all = TRUE)
# sgtfvoc: from growth_stat.R
sgtfvoc = fread("./data/sgtfvoc.csv")
d = merge(d, sgtfvoc[, .(specimen_date.x = date, NHSER_name = nhs_name, sgtfv)], by = c("specimen_date.x", "NHSER_name"))

# For impact of hospital pressure: get avg time from pillar 2 test to hospital admission and to death
delay_admit = d[!is.na(dateadmission_NHSE) & !is.na(specimen_date.x) & !is.na(dod) & pillar == "Pillar 2", 
    as.numeric(as.Date(dateadmission_NHSE) - as.Date(specimen_date.x))]
delay_death = d[!is.na(dateadmission_NHSE) & !is.na(specimen_date.x) & !is.na(dod) & pillar == "Pillar 2", 
    as.numeric(as.Date(dod) - as.Date(specimen_date.x))]

mean(delay_admit[delay_admit >= 0])
median(delay_admit[delay_admit >= 0])
mean(delay_death[delay_death >= 0])
median(delay_death[delay_death >= 0])

d[, hosp_date := specimen_date.x + 9]

# pressure: from build_sitrep_data.R
d = merge(d, pressure, by.x = c("hosp_date", "NHSER_code"), by.y = c("date", "nhscd"))

# Extract data for model
max_death_date = d[, max(dod, na.rm = T)]
d2 = d[!is.na(sgtf), 
    .(p_voc = sgtf * sgtfv, age = age.x, sex, NHSER_name, imd = imd_decile, specimen_date = specimen_date.x,
        date_death = dod, date_admission = dateadmission_NHSE,
        mv_pressure, ni_pressure, os_pressure, ao_pressure, medstaff_abs, nursing_abs,
        status = ifelse(is.na(dod), 0, 1),
        time = ifelse(is.na(dod), 
            as.numeric(as.Date(max_death_date) - as.Date(specimen_date.x)), 
            as.numeric(as.Date(dod) - as.Date(specimen_date.x))))]

d2 = d2[!is.na(age) & !is.na(sex) & !is.na(NHSER_name) & !is.na(imd)]
d2[, specimen_t := as.numeric(specimen_date - ymd("2020-01-01"))]
d2[, sex := factor(sex)]
d2[, NHSER_name := factor(NHSER_name)]

# Notes
# Try 28 day cutoff for deaths
# Pressure on most likely admission date
# Also look at SGTF
# SGTF 1/0 and multiple fits
# Different false positives rate
# > Interaction between age and VOC/SGTF - evaluate using e.g. AIC
# NHS region as strata
# Cut off a week earlier? But if someone dies after cutoff, move their cutoff back & they're still alive
# Give an absolute risk as well - comparison without voc i.e. p_voc = 0 using predict
# Check whether NAs are missing by random (MCAR)
# Schoenfeld residuals; Schoenfeld test (P value for association between residuals and variable)

# Model 0
model0 = gam(time ~ p_voc + s(age) + NHSER_name,
         family = cox.ph(), data = d2, weights = status)

summary(model0) 
plot(model0, pages = 1, all.terms = TRUE)
plot(model0$linear.predictors, residuals(model))

# Model 1
model1 = gam(time ~ p_voc + s(age) + nursing_abs + medstaff_abs + mv_pressure + ni_pressure + NHSER_name, # + sex + imd + s(specimen_t) + NHSER_name + s(pressure, by = NHSER_name)
         family = cox.ph(), data = d2, weights = status)

model1 = gam(cbind(time, NHSER_name) ~ p_voc + s(age) + s(pressure) + NHSER_name, # + sex + imd + s(specimen_t) + NHSER_name + s(pressure, by = NHSER_name)
         family = cox.ph(), data = d2, weights = status)

summary(model1) 
plot(model1, pages = 1, all.terms = TRUE)
plot(model1$linear.predictors, residuals(model))

# Model 2
model2 = gam(time ~ p_voc + age + NHSER_name*pressure + s(specimen_t), # + sex + imd + s(specimen_t) + NHSER_name + s(pressure, by = NHSER_name)
         family = cox.ph(), data = d2, weights = status)

summary(model2) 
plot(model2, pages = 1, all.terms = TRUE)
plot(model$linear.predictors, residuals(model))

# Model 3
model3 = gam(time ~ p_voc + s(age) + s(pressure) + NHSER_name + s(imd), # + sex + imd + s(specimen_t) + NHSER_name + s(pressure, by = NHSER_name)
         family = cox.ph(), data = d2, weights = status)

summary(model3) 
plot(model3, pages = 4, all.terms = TRUE)
plot(model3$linear.predictors, residuals(model))


# Model
model = gam(time ~ sgtf + s(age) + NHSER_name + s(specimen_t), # + sex + NHSER_name + imd,
         family = cox.ph(), data = d2[specimen_date >= ymd("2020-12-01")], weights = status)

summary(model) 
plot(model, pages = 1, all.terms = TRUE)
plot(model$linear.predictors, residuals(model))

# Case-control?
d3 = d[!is.na(sgtf), 
    .(p_voc = sgtf * sgtfv, age = prefer(age.y, age.x), age_group = (prefer(age.y, age.x) %/% 5) * 5, sex, UTLA_name, 
        imd = imd_decile, imd_group = ((imd_decile - 1) %/% 2) * 2 + 1,
        specimen_date = specimen_date.x, specimen_week = week(specimen_date.x),
        died = !is.na(dod))]
voc_y = d3[p_voc > 0.9] # TODO which cutoff to use?
voc_n = d3[p_voc < 0.1]
voc_y[, pid := .I]
voc_n[, pid := -.I]

# Weird possibilities: because of matching to Pillar 2 only, if the variant causes different test-seeking behaviour
# could the analysis be biased?

# Need to do conditional logistic regression
# Give everyone same length of time of follow up

matched_cohort = function(grp1, grp2, match_by)
{
    m = merge(grp1, grp2, by = match_by, allow.cartesian = TRUE)
    m[, wt.x := 1/.N, by = pid.x]
    m[, wt.y := 1/.N, by = pid.y]
    d1 = m[, sum(died.x * wt.x)]
    d2 = m[, sum(died.y * wt.y)]
    s1 = m[, sum(!died.x * wt.x)]
    s2 = m[, sum(!died.y * wt.y)]
    tbl = matrix(c(s2, d2, s2, d1), nrow = 2, byrow = TRUE)
    
    cat("Group 1 ", d1, ", group 2 ", d2, " P = ", binom.test(d1, d1 + d2)$p.value, "\n");
    cat("Odds ratio:\n");
    print(OddsRatio(tbl, conf.level = 0.95));
}

matched_cohort(voc_y, voc_n, c("age", "sex", "UTLA_name", "specimen_date"))
matched_cohort(voc_y, voc_n, c("age_group", "sex", "UTLA_name", "specimen_date"))
matched_cohort(voc_y, voc_n, c("age_group", "imd", "sex", "UTLA_name", "specimen_date"))
matched_cohort(voc_y, voc_n, c("age", "imd", "sex", "UTLA_name", "specimen_date"))

matched_cohort(voc_y, voc_n, c("age", "sex", "UTLA_name", "specimen_week"))
matched_cohort(voc_y, voc_n, c("age_group", "sex", "UTLA_name", "specimen_week"))
matched_cohort(voc_y, voc_n, c("age_group", "imd", "sex", "UTLA_name", "specimen_week"))
matched_cohort(voc_y, voc_n, c("age_group", "imd_group", "sex", "UTLA_name", "specimen_week"))
matched_cohort(voc_y, voc_n, c("age", "imd", "sex", "UTLA_name", "specimen_week"))

voc_m1 = merge(voc_y, voc_n, by = c("age", "sex", "UTLA_name", "specimen_date"))
voc_m1s = voc_m1[, .SD[sample(.N, 1)], by = pid.x]
voc_m1s[, .(sum(died.x), sum(died.y))]

voc_m2 = merge(voc_y, voc_n, by = c("age_group", "sex", "UTLA_name", "specimen_date"))
voc_m2s = voc_m2[, .SD[sample(.N, 1)], by = pid.x]
voc_m2s[, .(sum(died.x), sum(died.y))]

voc_m3 = merge(voc_y, voc_n, by = c("age_group", "imd", "sex", "UTLA_name", "specimen_date"))
voc_m3s = voc_m3[, .SD[sample(.N, 1)], by = pid.x]
voc_m3s[, .(sum(died.x), sum(died.y))]

voc_m4 = merge(voc_y, voc_n, by = c("age", "imd", "sex", "UTLA_name", "specimen_date"))
voc_m4s = voc_m4[, .SD[sample(.N, 1)], by = pid.x]
voc_m4s[, .(sum(died.x), sum(died.y))]
