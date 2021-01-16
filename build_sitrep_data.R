library(data.table)
library(readxl)
library(ggplot2)
library(lubridate)
library(ogwrangler)

source("./commit.R")

# Note: this file contains sensitive NHS data, so it is not included wth the repo.
spim_data = read_excel("~/Documents/uk_covid_data_sensitive/spi-m_data/20210115_All_SPIM_trust_002.xlsx", "Extracted Data", col_type = "text")
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
existing = fread("~/Dropbox/uk_covid_data/data-2021-01-08.csv")
existing[, date := as.Date(date)]
committed = commit(existing, data)
fwrite(committed, "~/Dropbox/uk_covid_data/data-2021-01-15.csv")
fwrite(committed, "./fitting_data/data-2021-01-15.csv")

#######

# Experimental: pull metrics from coronavirus.gov.uk web site
endpoint_hosp = 'https://api.coronavirus.data.gov.uk/v2/data?areaType=nhsRegion&metric=covidOccupiedMVBeds&metric=newAdmissions&metric=hospitalCases&format=csv'
response_hosp = httr::GET(url = endpoint_hosp, httr::timeout(10))
d_hosp = fread(httr::content(response_hosp, "text"))
d_hosp_all = rbind(
    melt_annotate(d_hosp[!is.na(newAdmissions), .(location = areaName, date = as.Date(date), hospital_inc = newAdmissions)]),
    melt_annotate(d_hosp[!is.na(hospitalCases), .(location = areaName, date = as.Date(date), hospital_prev = hospitalCases)]),
    melt_annotate(d_hosp[!is.na(covidOccupiedMVBeds), .(location = areaName, date = as.Date(date), icu_prev = covidOccupiedMVBeds)])
)

d_death = read_excel("~/Documents/uk_covid_data_sensitive/phe/20210101/20210101 COVID19 Deaths.xlsx")
setDT(d_death)
d_death = melt_annotate(d_death[!is.na(nhser_name) & death_type28 == 1 & !is.na(dod), .(death_inc_line = .N), keyby = .(location = nhser_name, date = as.Date(dod))])

# Save
existing = fread("~/Dropbox/uk_covid_data/data-2020-12-30.csv")
existing[, date := as.Date(date)]
committed = commit(existing, rbind(d_hosp_all, d_death))
fwrite(committed, "~/Dropbox/uk_covid_data/data-2021-01-01.csv")
fwrite(committed, "./fitting_data/data-2021-01-01.csv")

#######

# Available beds, staff absences, etc
dirs = list.files(path = "~/Documents/uk_covid_data_sensitive/folder_20127408/Archive/", include.dirs = TRUE)
alldat = NULL;
for (dir in dirs)
{
    cat(".")
    searchdir = paste0("~/Documents/uk_covid_data_sensitive/folder_20127408/Archive/", dir);
    file = list.files(path = searchdir, pattern = "^Covid.sitrep.report.*(xlsx|XLSX)$")
    dat = as.data.table(read_excel(paste0(searchdir, "/", file), "R Data", skip = 1))
    dat = dat[, .(date = as.Date(period), region = Region, sitecode = `Site/Org Code`, type = `Organisation Type`,
        mv_occ_cov = SIT032_OccupiedCov, mv_occ_sus = SIT032_OccupiedSuspected, mv_occ_non = SIT032_OccupiedNonCovNS, mv_unocc = SIT032_Unoccupied, # mechanical ventilation beds
        ni_occ_cov = SIT033_OccupiedCov, ni_occ_sus = SIT033_OccupiedSuspected, ni_occ_non = SIT033_OccupiedNonCovNS, ni_unocc = SIT033_Unoccupied, # noninvasive ventilation beds
        os_occ_cov = SIT034_OccupiedCov, os_occ_sus = SIT034_OccupiedSuspected, os_occ_non = SIT034_OccupiedNonCovNS, os_unocc = SIT034_Unoccupied, # oxygenation support beds
        ao_occ_cov = SIT058_OccupiedCov, ao_occ_sus = SIT058_OccupiedSuspected, ao_occ_non = SIT058_OccupiedNonCovNS, ao_unocc = SIT058_Unoccupied, # any other beds
        medstaff_absent = SIT048_TotalAbsent, nursing_absent = SIT049_TotalAbsent)]
    dat[, occ_all :=
        mv_occ_cov + mv_occ_sus + mv_occ_non +
        ni_occ_cov + ni_occ_sus + ni_occ_non +
        os_occ_cov + os_occ_sus + os_occ_non +
        ao_occ_cov + ao_occ_sus + ao_occ_non]
    dat[, occ_cov :=
        mv_occ_cov + mv_occ_sus +
        ni_occ_cov + ni_occ_sus +
        os_occ_cov + os_occ_sus +
        ao_occ_cov + ao_occ_sus]
    
    # Correct dates improperly entered as 2020 when they should be 2021
    dat[date <= "2020-06-01", date := date + 366]

    alldat = rbind(alldat, dat)
}

CreateCache()
etrust = fread("~/Dropbox/uk_covid_data/etrust/etrust.csv", header = FALSE)[, .(sitecode = V1, postcode = V10)]

alldat2 = merge(alldat[type == "Acute Trust"], etrust, by = "sitecode", all.x = TRUE)
alldat2[, ladcd := ogpost(postcode, geo = "lad")]
alldat2[, nhscd := ogpost(postcode, geo = "nhser")]

w = alldat2[, .(
    mv_pressure = 1 - sum(mv_unocc) / (sum(mv_occ_cov + mv_occ_sus + mv_occ_non + mv_unocc)),
    ni_pressure = 1 - sum(ni_unocc) / (sum(ni_occ_cov + ni_occ_sus + ni_occ_non + ni_unocc)),
    os_pressure = 1 - sum(os_unocc) / (sum(os_occ_cov + os_occ_sus + os_occ_non + os_unocc)),
    ao_pressure = 1 - sum(ao_unocc) / (sum(ao_occ_cov + ao_occ_sus + ao_occ_non + ao_unocc)),
    medstaff_abs = sum(medstaff_absent),
    nursing_abs = sum(nursing_absent)
    ), keyby = .(date, nhscd)]

# Plot data
w2 = melt(w, id.vars = 1:2)
ggplot(w2) + geom_line(aes(x = date, y = value, colour = nhscd)) + facet_wrap(~variable, scales = "free")

# Has outliers. Remove:
smooth_blips = function(x, date)
{
    x = ifelse(abs(x - zoo::rollmean(x, 7, "extend")) > 2 * sd(x) | date %between% ymd(c("2020-12-21", "2020-12-29")), NA, x)
    zoo::na.approx(x)
}

w2 = w2[, .(date, value = smooth_blips(value, date)), by = .(nhscd, variable)]
ggplot(w2) + geom_line(aes(x = date, y = value, colour = nhscd)) + facet_wrap(~variable, scales = "free")

# Standardise absences
w2[variable %in% c("medstaff_abs", "nursing_abs"), value := (value - mean(value)) / sd(value), by = .(nhscd, variable)]
w2[, value := zoo::rollmean(value, 7, fill = "extend"), by = .(nhscd, variable)]
ggplot(w2) + geom_line(aes(x = date, y = value, colour = nhscd)) + facet_wrap(~variable, scales = "free")

pressure = dcast(w2, nhscd + date ~ variable)
fwrite(pressure, "~/Documents/uk_covid_data_sensitive/pressure.csv")
