library(data.table)
library(readxl)
library(ggplot2)
library(lubridate)

source("./commit.R")

spim_data = read_excel("~/Documents/uk_covid_data_sensitive/spi-m_data/20201218_All_SPIM_trust_001.xlsx", "Extracted Data", col_type = "text")
setDT(spim_data)
spim_data[, DateVal := as.Date(DateVal)]
spim_data = cbind(
    spim_data[, .SD, .SDcols = c(1, 5, 6, 7, 8)],
    spim_data[, lapply(.SD, as.numeric), .SDcols = c(2, 3, 4, 9:ncol(spim_data))]
)

regions = c("East of England", "London", "Midlands", "North East and Yorkshire",
    "North West", "South West",  "South East", "Northern Ireland", "Scotland", "Wales")
rlevels = c("Region", "National")

# hospital_prev and icu_prev
hpEN = melt_annotate(spim_data[Geography %in% regions & !(Geography %in% c("Scotland", "Wales")) & ReportLevel %in% rlevels & !is.na(hospital_prev), hospital_prev,
  keyby = .(location = Geography, date = ymd(DateVal))])

hpW = melt_annotate(spim_data[Geography %in% regions & Geography == "Wales" & ReportLevel %in% rlevels & !is.na(hospital_prev), .(hospital_prev = hospital_prev + hospital_prev_recovering),
  keyby = .(location = Geography, date = ymd(DateVal))])

hpS = melt_annotate(spim_data[Geography == "Scotland" & ReportLevel %in% rlevels & !is.na(`hospital_prev_<28days`), .(hospital_prev = `hospital_prev_<28days`),
  keyby = .(location = Geography, date = ymd(DateVal))])

hp = rbind(hpEN, hpW, hpS)

ipENW = melt_annotate(spim_data[Geography %in% regions & Geography != "Scotland" & ReportLevel %in% rlevels & !is.na(icu_prev), icu_prev,
  keyby = .(location = Geography, date = ymd(DateVal))])

ipS = melt_annotate(spim_data[Geography == "Scotland" & ReportLevel %in% rlevels & !is.na(`icu_prev_<28days`), .(icu_prev = `icu_prev_<28days`),
  keyby = .(location = Geography, date = ymd(DateVal))])

ip = rbind(ipENW, ipS)

# hospital_inc
sumNA = function(x, y)
{
    ifelse(is.na(x) & is.na(y), NA, sum(c(x, y), na.rm = T))
}

hiENS = melt_annotate(spim_data[Geography %in% regions & Geography != "Wales" & ReportLevel %in% rlevels & !is.na(sumNA(hospital_inc, hospital_inc_new)),
  .(hospital_inc = sumNA(hospital_inc, hospital_inc_new)),
  keyby = .(location = Geography, date = ymd(DateVal))])

hiW = melt_annotate(spim_data[Geography %in% regions & ReportLevel %in% rlevels & !is.na(SPIM_hosp_inc_new),
  .(hospital_inc = SPIM_hosp_inc_new),
  keyby = .(location = Geography, date = ymd(DateVal))])

hi = rbindlist(list(hiENS, hiW[date >= "2020-02-01"]))

# death_inc_line
diE = melt_annotate(spim_data[Geography %in% regions & ReportLevel %in% rlevels & !is.na(PHE_type28_death_inc_line),
  .(death_inc_line = PHE_type28_death_inc_line),
  keyby = .(location = Geography, date = ymd(DateVal))])

# ### TEMP: just hospital deaths
# diE = melt_annotate(spim_data[Geography %in% regions & ReportLevel %in% rlevels & !is.na(CPNS_death_inc_line),
#   .(death_inc_line = CPNS_death_inc_line),
#   by = .(location = Geography, date = ymd(DateVal))])

# diE_0_4 = melt_annotate(spim_data[Geography %in% regions & ReportLevel %in% rlevels & !is.na(`PHE_type28_death_inc_line<1 age`),
#   .(death_inc_line_0_4 = `PHE_type28_death_inc_line<1 age` + `PHE_type28_death_inc_line1-4 age`),
#   by = .(location = Geography, date = ymd(DateVal))])
#
# diE_5_14 = melt_annotate(spim_data[Geography %in% regions & ReportLevel %in% rlevels & !is.na(PHE_type28_death_inc_line),
#   .(death_inc_line_5_14 = `PHE_type28_death_inc_line5-14 age`),
#   by = .(location = Geography, date = ymd(DateVal))])
#
# diE_15_24 = melt_annotate(spim_data[Geography %in% regions & ReportLevel %in% rlevels & !is.na(PHE_type28_death_inc_line),
#   .(death_inc_line_15_24 = `PHE_type28_death_inc_line15-24 age`),
#   by = .(location = Geography, date = ymd(DateVal))])
#
# diE_25_44 = melt_annotate(spim_data[Geography %in% regions & ReportLevel %in% rlevels & !is.na(PHE_type28_death_inc_line),
#   .(death_inc_line_25_44 = `PHE_type28_death_inc_line25-44 age`),
#   by = .(location = Geography, date = ymd(DateVal))])
#
# diE_45_54 = melt_annotate(spim_data[Geography %in% regions & ReportLevel %in% rlevels & !is.na(PHE_type28_death_inc_line),
#   .(death_inc_line_45_54 = `PHE_type28_death_inc_line45-54 age`),
#   by = .(location = Geography, date = ymd(DateVal))])
#
# diE_55_64 = melt_annotate(spim_data[Geography %in% regions & ReportLevel %in% rlevels & !is.na(PHE_type28_death_inc_line),
#   .(death_inc_line_55_64 = `PHE_type28_death_inc_line55-64 age`),
#   by = .(location = Geography, date = ymd(DateVal))])
#
# diE_65_74 = melt_annotate(spim_data[Geography %in% regions & ReportLevel %in% rlevels & !is.na(PHE_type28_death_inc_line),
#   .(death_inc_line_65_74 = `PHE_type28_death_inc_line65-74 age`),
#   by = .(location = Geography, date = ymd(DateVal))])
#
# diE_75_pl = melt_annotate(spim_data[Geography %in% regions & ReportLevel %in% rlevels & !is.na(PHE_type28_death_inc_line),
#   .(death_inc_line_75_pl = `PHE_type28_death_inc_line75-84 age` + `PHE_type28_death_inc_line>84 age`),
#   by = .(location = Geography, date = ymd(DateVal))])

diN = melt_annotate(spim_data[Geography == "Northern Ireland" & ReportLevel %in% rlevels & !is.na(SitRep_death_inc_line),
  .(death_inc_line = SitRep_death_inc_line),
  keyby = .(location = Geography, date = ymd(DateVal))])

diS = spim_data[Geography == "Scotland" & ReportLevel %in% rlevels & !is.na(SitRep_death_inc_line),
  .(death_inc_line = SitRep_death_inc_line),
  keyby = .(location = Geography, date = ymd(DateVal))]
zero_filler = data.table(location = "Scotland", date = min(diS$date) + 0:as.numeric(max(diS$date) - min(diS$date)))
diS = merge(zero_filler, diS, by = c("location", "date"), all = TRUE)
diS[is.na(death_inc_line), death_inc_line := 0]
diS = melt_annotate(diS)

diW = melt_annotate(spim_data[Geography == "Wales" & ReportLevel %in% rlevels & !is.na(PHW_Death_inc_line),
  .(death_inc_line = PHW_Death_inc_line),
  keyby = .(location = Geography, date = ymd(DateVal))])

di = rbindlist(list(diE, diN, diS, diW))

# Amalgamate
data = rbindlist(list(hp, ip, hi, di))

# Save
existing = fread("~/Dropbox/uk_covid_data/data-2020-12-15.csv")
existing[, date := as.Date(date)]
committed = commit(existing, data)
fwrite(committed, "~/Dropbox/uk_covid_data/data-2020-12-18.csv")

#######


w = spim_data[Geography %in% regions & ReportLevel %in% rlevels & !is.na(hospital_inc_acute1), hospital_inc_new,
  by = .(location = Geography, date = ymd(DateVal), ReportLevel)]

# hospital_inc
# English regions: hospital_inc
# Wales: SPIM_hosp_inc_new (?)
# NI: hospital_inc_new?
spim_data[Geography == "East of England", lapply(.SD, function(x) sum(!is.na(x))), .SDcols = 9:286]

# 28 day deaths
# Northern Ireland: SitRep_death_inc_dash SitRep_hospital_death_inc_line SitRep_death_inc_line ?

ggplot(w) + geom_line(aes(x = date, y = hospital_inc_new)) + facet_wrap(~location)

setwd("~/Documents/Git/covid19_automation/")
source("scripts/install_devel_packages.R")
reportfactory::install_deps()
reportfactory::rfh_load_scripts()

most_recent("~/Documents/uk_covid_data_sensitive/sitrep_n", "NHSitrep", ymd) # sheet 2 2-row header, sheet 3 normal
most_recent("~/Documents/uk_covid_data_sensitive/sitrep_n", "NIExport", ymd) # sheet 2 normal

most_recent("~/Documents/uk_covid_data_sensitive/sitrep_s", "C19Inpatients", ymd) # sheet 1 horizontal, sheet 2 horizontal
most_recent("~/Documents/uk_covid_data_sensitive/sitrep_s", "Positive_Hospitals", dmy) # csv normal
most_recent("~/Documents/uk_covid_data_sensitive/sitrep_s", "COVID19_deaths", ymd) # excel normal

most_recent("~/Documents/uk_covid_data_sensitive/sitrep_w", "Sitrep", ymd_hm) # csv normal with lookup


# Overall data
data = NULL;

# NORTHERN IRELAND
# source 1: sitrep, sheet 2
path = most_recent("~/Documents/uk_covid_data_sensitive/sitrep_n", "NHSitrep", ymd);
d1 = read_2_row_header(path, 2);
d1 = d1[, .(date = ymd(paste(2020, Date, Date.1)), death_inc_line = as.numeric(Dthsbydtfdthh), icu_prev = as.numeric(ICUc))];

# SKIPPING source 2: sitrep, sheet 3 (death by age & place of death)

# source 3: export, sheet 2
path = most_recent("~/Documents/uk_covid_data_sensitive/sitrep_n", "NIExport", ymd);
d2 = read_normal(path, 2);
d2 = d2[, .(date = ymd(Date), hospital_prev = as.numeric(Inpatients), hospital_inc = as.numeric(Admissions))];

# merge
d = merge(d1, d2, by = "date", all = T);
d[, location := "Northern Ireland"];
data = rbind(data, melt_annotate(d));


# SCOTLAND
# source 1: inpatients, sheet 1
path = most_recent("~/Documents/uk_covid_data_sensitive/sitrep_s", "C19Inpatients", ymd);
d1 = read_horiz_date(path, 1, cellranger::cell_rows(5:21), "location", "hospital_prev")[location == "Scotland"];

# source 2: inpatients, sheet 2
path = most_recent("~/Documents/uk_covid_data_sensitive/sitrep_s", "C19Inpatients", ymd);
d2 = read_horiz_date(path, 2, cellranger::cell_rows(5:21), "location", "icu_prev")[location == "Scotland"];

# source 3: tests in hospitals
path = most_recent("~/Documents/uk_covid_data_sensitive/sitrep_s", "Positive_Hospitals", dmy);
d3 = read_normal(path, NA);
d3 = d3[, .(date = ymd(EcossDateReceived), hospital_inc = N.Positive.In.Hospital)];

# source 4: COVID deaths
path = most_recent("~/Documents/uk_covid_data_sensitive/sitrep_s", "COVID19_deaths", ymd);
d4 = read_normal(path, 1);
d4 = d4[, .(death_inc_line = .N), keyby = .(date = ymd(NRS.Date.Death))];
d4 = blank_fill(d4, "date", 0);

# merge
d = merge(merge(merge(d1, d2, by = c("date", "location"), all = T), d3, by = "date", all = T), d4, by = "date", all = T);
data = rbind(data, melt_annotate(d));

# WALES
# export = read_wales_bulk_export("~/Documents/uk_covid_data_sensitive/sitrep_w");
# A A 8 5
#
# ggplot(export[Dataset == "A" & Section == "A" & Question == 8, sum(MeasureValue), by = .(Measure, date)]) +
#     geom_point(aes(x = date, y = V1, colour = Measure)) + facet_wrap(~Measure)
#
# export[date == "2020-06-01", sum(MeasureValue), by = .(date, Dataset, Section, Question, Measure)][V1 == 28]
#
# path = most_recent("~/Documents/uk_covid_data_sensitive/sitrep_w", "Sitrep", ymd_hm);
# d1 = read_normal(path, NA);
# d1[, date := ymd_hms(UpdateDateTime)]
#
# existing[location == "Wales" & date == "2020-06-01"]
# #          date location      indicator value                                  description
# # 1: 2020-06-01    Wales death_inc_line     9                All deaths (by date of death)
# # 2: 2020-06-01    Wales   hospital_inc    25 New and newly confirmed patients in hospital
# # 3: 2020-06-01    Wales  hospital_prev   712                          Total beds occupied
# # 4: 2020-06-01    Wales       icu_prev    28                            ICU beds occupied
#
# d1[, sum(MeasureValue), by = .(UpdateDateTime, Dataset, Section, Question, Measure)][V1 == 25]


# TODO

# ENGLAND
# source 1: sitrep list
# sl = import_sitrep_list_eng()
# s = rbindlist(sl, fill = T)

# TEMP source 1: sitrep series
d1 = data.table(import_sitrep_series_eng());
d1 = d1[, .(hospital_inc = sum(hospital_inc), hospital_prev = sum(hospital_prev), icu_prev = sum(icu_prev)), by = .(date = ymd(date), location = region)]

# source 2: deaths linelist
d2 = data.table(import_deaths_eng());
d2[, location := str_replace_all(nhser_name, "_", " ")];
d2[, location := str_to_title(location)];
d2[, location := str_replace_all(location, "Of", "of")];
d2[, location := str_replace_all(location, "And", "and")];
d2 = d2[location != "Missing" & death_type == "lab_confirmed" & death_type28 == 1, .(death_inc_line = .N), keyby = .(date = ymd(date_death), location)];
d2 = blank_fill(d2, "date", 0);

# merge
d = merge(d1, d2, by = c("date", "location"), all = T);
data = rbind(data, melt_annotate(d));

#existing = readRDS("~/Documents/covid-fit/data/forecast_data_2020-07-07.rds");
#names(existing) = c("date", "location", "indicator", "value", "type", "description")
#existing$type = NULL
existing = fread("~/Dropbox/uk_covid_data/test_data.csv")
existing$date = as.Date(existing$date)

goforit = commit(existing, existing)
goforit = commit(existing, data)
goforit = goforit[!is.na(date) & !is.na(location) & !is.na(value) & !(location %in% c("England", "Great Britain", "United Kingdom", "0"))]
ggplot(goforit) + geom_line(aes(x = date, y = value, colour = location)) + facet_wrap(~indicator, scales = "free")
fwrite(goforit, "~/Dropbox/uk_covid_data/test_data.csv")





# test...
file = most_recent("~/Documents/uk_covid_data_sensitive/sitrep_n", "NHSitrep", ymd);
data = read_excel(file, sheet = 2, skip = 2,
    col_names = c("month", "day", "pillar1_n", "pillar2_n", "pillar1_pos", "pillar2_pos", "pillar1_nonpos", "pillar2_nonpos",
        "icu_prev", "deaths_reported", "deaths_reported_today", "deaths", "pillar1_tests", "pillar2_tests"));
data = read_excel(file, sheet = 2);
setDT(data);

data = read_excel(file, sheet = 3, skip = 1,
    col_names = c("date", "lgd", "gender", "age", "deaths"));
setDT(data);
data


file = most_recent("~/Documents/uk_covid_data_sensitive/sitrep_n", "NIExport", ymd);
data = read_excel(file, sheet = 2, skip = 0,
    col_names = TRUE);
setDT(data);

?ymd
tail(test)
file

# temp
existing = rlang::duplicate(fd)
existing[, type := NULL]
names(existing) = c("date", "location", "indicator", "value", "description")
existing
fwrite(existing, "~/Dropbox/uk_covid_data/sitrep.csv")

# Visualize sitrep data
setwd("~/Documents/Git/covid19_automation/")
source("scripts/install_devel_packages.R")
reportfactory::install_deps()
reportfactory::rfh_load_scripts()

sl = import_sitrep_list_eng()
s = rbindlist(sl, fill = T)

types = s[, unique(value_type)]
plots = list()
date_limits = s[, c(min(date), max(date))]
for (i in seq_along(types)) {
    cat(".")
    plot = ggplot(s[value_type == types[i], .(value = sum(value, na.rm = T)), by = .(value_type, value_desc, date, region)]) +
        geom_point(aes(x = date, y = value, colour = str_sub(region, 1, 10)), size = 0.1) +
        geom_line(aes(x = date, y = value, colour = str_sub(region, 1, 10))) +
        labs(title = paste0(i, ". ", types[i]), subtitle = s[value_type == types[i], str_wrap(value_desc[1])], colour = "region") +
        scale_x_date(date_breaks = "1 month", date_labels = "%b", limits = as.Date(date_limits))
    ggsave(paste0("~/Dropbox/uk_covid_data/indicators/page", sprintf("%03i", i), ".pdf"), plot, width = 18, height = 12, units = "cm")
}
system('"/System/Library/Automator/Combine PDF Pages.action/Contents/Resources/join.py" -o ~/Dropbox/uk_covid_data/indicators/indicators.pdf ~/Dropbox/uk_covid_data/indicators/page*.pdf')

ss = import_sitrep_series_eng()
setDT(ss)

ggplot(ss[, .(v = sum(n_hdu_itu)), by = .(region, date)]) +
  geom_line(aes(x = date, y = v, colour = region))

ggplot(ss[, .(v = sum(icu_prev)), by = .(region, date)]) +
  geom_line(aes(x = date, y = v, colour = region))

# Visualize CO-CIN data
cocin = fread("~/Documents/covid-fit/data/cocin-7-july-2020/CCPUKSARI_DATA_2020-07-07_1244.csv")

cocin

# Deaths data
d = import_deaths_eng()
ggplot(d) + geom_histogram(aes(age), binwidth = 1)



library(shmanipulate)
shmanipulate({
    d = data.frame(age = seq(2.5, 77.5, by = 5));
    alpha = 0.2 * (concentration - 2) + 1;
    beta = 0.8 * (concentration - 2) + 1;
    d$multiplier = (1 - constant) * dbeta(d$age / 80, alpha, beta) + constant;
    ggplot(d) + geom_col(aes(age, multiplier))
}, concentration = c(1, 10), constant = c(0, 1))


# READ ENGLAND SITREPS

library(readxl)

#s0 = data.table(read_excel("~/Documents/uk_covid_data_sensitive/sitrep_e/Covid_sitrep_report_20200925_R_FINAL.xlsx", "R Data", skip = 1))
s1 = data.table(read_excel("~/Documents/uk_covid_data_sensitive/sitrep_e/Covid_sitrep_report_20200926_R_Final.xlsx", "R Data", skip = 1))
s2 = data.table(read_excel("~/Documents/uk_covid_data_sensitive/sitrep_e/Covid_sitrep_report_20200927_R_FINAL.xlsx", "R Data", skip = 1))
s3 = data.table(read_excel("~/Documents/uk_covid_data_sensitive/sitrep_e/Covid_sitrep_report_20200928_R_FINAL.xlsx", "R Data", skip = 1))

newdata = rbindlist(list(s1, s2, s3), use.names = T)
newdata = newdata[, .(hospital_inc = sum(SIT008_Total + SIT009_Total),
            hospital_prev = sum(SIT032_OccupiedCov + SIT033_OccupiedCov + SIT034_OccupiedCov + SIT058_OccupiedCov),
            icu_prev = sum(SIT054_OccupiedCov)),
  by = .(Region, period)]

newdata[, Region := str_to_title(Region)]
newdata[, Region := str_replace_all(Region, "Of", "of")]
newdata[, Region := str_replace_all(Region, "And", "and")]
names(newdata)[1:2] = c("location", "date")
newdata[, date := as.Date(date)]


existing = fread("~/Dropbox/uk_covid_data/test_data.csv")
existing = existing[!duplicated(existing, by = c("date", "location", "indicator"))]
existing$date = as.Date(existing$date)
data = rbind(existing, melt_annotate(newdata))

goforit = commit(existing, data)
fwrite(goforit, "~/Dropbox/uk_covid_data/test_data.csv")


goforit

newdata



# diagnosis time
d1 = data.table(import_sitrep_series_eng());

dd = d1[, .(
  adm = sum(n_patients_admitted),
  diag = sum(n_inpatients_diagnosed),
  diag_0_2  = sum(n_inpatients_diagnosed_0_2),
  diag_3_7  = sum(n_inpatients_diagnosed_3_7),
  diag_8_14 = sum(n_inpatients_diagnosed_8_14),
  diag_15_  = sum(n_inpatients_diagnosed_15_)
  ), by = .(dateweek = round_date(date, unit = "1 week"), region)]

ggplot(dd, aes(dateweek)) +
    geom_line(aes(y = diag, colour = "total")) +
    geom_line(aes(y = diag_0_2, colour = "0-2")) +
    geom_line(aes(y = diag_3_7, colour = "3-7")) +
    geom_line(aes(y = diag_8_14, colour = "8-14")) +
    geom_line(aes(y = diag_15_, colour = "15+")) +
    geom_line(aes(y = diag_0_2 + diag_3_7 + diag_8_14 + diag_15_, colour = "sum"))

ggplot(dd, aes(dateweek)) +
    geom_point(aes(y = (diag_8_14 + diag_15_) / (adm + diag))) +
    geom_smooth(aes(y = (diag_8_14 + diag_15_) / (adm + diag))) +
    facet_wrap(~region) +
    ylim(-0.1, 0.4)


# Avg death age
d2 = data.table(import_deaths_eng());
ggplot(d2[, age, keyby = date_death]) +
    geom_jitter(aes(x = date_death, y = age), size = 0.3)# +
#    geom_smooth(aes(x = date_death, y = age)) +
#    geom_hline(aes(yintercept = 76), colour = "blue") +
#    geom_hline(aes(yintercept = 83), colour = "red")



# Experimental: pull deaths from coronavirus.gov.uk web site
endpoint = 'https://api.coronavirus.data.gov.uk/v1/data?filters=areaType=nation&structure={"date":"date","deaths":"newDeaths28DaysByDeathDate","nation":"areaName"}'

response = httr::GET(url = endpoint, httr::timeout(10))

if (response$status_code >= 400) {
    err_msg = httr::http_status(response)
    stop(err_msg)
}

data = jsonlite::fromJSON(httr::content(response, "text"))
d = as.data.table(data$data)
d[, latest_date := max(date), by = nation]

govt_deaths = d
spim_deaths = spim_data[Geography %in% c("Northern Ireland", "Scotland", "Wales"), SitRep_death_inc_line, by = .(DateVal, Geography)]
phw_deaths = spim_data[Geography %in% c("Wales"), PHW_Death_inc_line, by = .(DateVal, Geography)]
govt_deaths
names(spim_deaths) = c("date", "nation", "deaths")
names(phw_deaths) = c("date", "nation", "deaths")
govt_deaths[, date := as.Date(date)]
spim_deaths[, date := as.Date(date)]
phw_deaths[, date := as.Date(date)]

ggplot() +
    geom_point(data = spim_deaths, aes(x = date, y = deaths+2, colour = "SitRep_death_inc_line"), size = 0.1) +
    geom_point(data = govt_deaths[nation != "England"], aes(x = date, y = deaths, colour = "Govt 28-day deaths data"), size = 0.1) +
    geom_point(data = phw_deaths, aes(x = date, y = deaths+1, colour = "PHW_Death_inc_line"), size = 0.1) +
    geom_line(data = spim_deaths, aes(x = date, y = deaths+2, colour = "SitRep_death_inc_line"), size = 0.1) +
    geom_line(data = govt_deaths[nation != "England"], aes(x = date, y = deaths, colour = "Govt 28-day deaths data"), size = 0.1) +
    geom_line(data = phw_deaths, aes(x = date, y = deaths+1, colour = "PHW_Death_inc_line"), size = 0.1) +
    facet_wrap(~nation, ncol = 1) +
    scale_x_date(date_breaks = "1 month", date_labels = "%b", limits = c(ymd("2020-03-01"), NA))

