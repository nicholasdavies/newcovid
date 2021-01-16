library(data.table)
library(ggplot2)
library(lubridate)
library(here)
library(cowplot)
library(readxl)
library(sn)
library(qs)
library(stringr)
 
uk_covid_data_path = "./fitting_data/";
datapath = function(x) paste0(uk_covid_data_path, x)

# Sys.setenv(MP_DUPLICATE_LIB_OK=TRUE) # RCB needs to override this environment variable to run

#
# SETUP
#

# set up covidm
cm_path = "./covidm_for_fitting/";
cm_force_rebuild = F;
cm_build_verbose = T;
cm_version = 2;
source(paste0(cm_path, "/R/covidm.R"))
popUK = readRDS(datapath("popNHS.rds"));
matricesUK = readRDS(datapath("matricesNHS.rds"));

cm_populations = rbind(cm_populations[name != "United Kingdom"], popUK)
cm_matrices = c(cm_matrices, matricesUK)
source("./distribution_fit.R");
source("./spim_output.R");

#
# DATA
#

nhs_regions = popUK[, unique(name)]
pct = function(x) as.numeric(str_replace_all(x, "%", "")) / 100

# Process new data
new_data = fread(datapath("data-2021-01-15.csv"))

# removing problematic data with under reporting of hospital incidence
days_to_remove_scot <- 1  ## VISUALLY INSPECT DATA TO DETERMINE THESE VALUES
days_to_remove_wales <- 1 ## VISUALLY INSPECT DATA TO DETERMINE THESE VALUES
scottish_data <- new_data[new_data$location == "Scotland"]
scottish_hospital_inc <- scottish_data[scottish_data$indicator == "hospital_inc"]
max_date_scotland <- max(new_data[new_data$location=="Scotland" & 
                                    new_data$indicator=="hospital_inc"]$date)
dates_to_remove_scotland <- seq(scottish_hospital_inc$date[length(scottish_hospital_inc$date)
                                                           -(days_to_remove_scot-1)], 
                                by = "day", 
                                length.out = (max_date_scotland-
                                                scottish_hospital_inc$date[length(scottish_hospital_inc$date)
                                                                           -(days_to_remove_scot)]))
welsh_data <- new_data[new_data$location == "Wales"]
welsh_hospital_inc <- welsh_data[welsh_data$indicator == "hospital_inc"]
max_date_wales <- max(new_data[new_data$location=="Wales" & 
                                 new_data$indicator=="hospital_inc"]$date)
dates_to_remove_wales <- seq(welsh_hospital_inc$date[length(welsh_hospital_inc$date)
                                                     -(days_to_remove_wales-1)],
                            by = "day",
                             length.out = (max_date_wales-
                                             welsh_hospital_inc$date[length(welsh_hospital_inc$date)
                                                                     -(days_to_remove_wales)]))
new_data <- new_data[!(new_data$location == "Scotland" & 
                         new_data$indicator == "hospital_inc" & 
                         new_data$date %in% dates_to_remove_scotland)]
new_data <- new_data[!(new_data$location == "Wales" & 
                         new_data$indicator == "hospital_inc" & 
                         new_data$date %in% dates_to_remove_wales)]

new_data[, pid := match(location, nhs_regions) - 1]
ld = new_data[indicator == "death_inc_line", .(date, N = value, name = location, pid)];
sitreps = new_data[indicator != "death_inc_line", .(date = ymd(date), value_type = indicator, value, name = location, pid)]
sitreps = dcast(sitreps, date + name + pid ~ value_type, value.var = "value", fill = NA, fun.aggregate = function(x) x[1])
sitreps = sitreps[, .(date, n_in_itu = icu_prev, n_in_all_beds = hospital_prev, n_admitted_diagnosed = hospital_inc, name, pid)];


sero = fread(datapath("seroprev_nhs_regions_20201121.csv"))
sero = cbind(sero, sero[, skew_normal_solve(pct(Central.estimate), pct(Lower.bound), pct(Upper.bound), Data.source == "REACT-2")]);
sero[Data.source == "REACT-2", Central.estimate := mapply(function(x,o,a) paste0(round(100 * qsn(0.5, x, o, a), 2), "%"), xi, omega, alpha)]
sero[Data.source == "REACT-2", Lower.bound := mapply(function(x,o,a) paste0(round(100 * qsn(0.025, x, o, a), 2), "%"), xi, omega, alpha)]
sero[Data.source == "REACT-2", Upper.bound := mapply(function(x,o,a) paste0(round(100 * qsn(0.975, x, o, a), 2), "%"), xi, omega, alpha)]
sero[, Start.date := dmy(Start.date)]
sero[, End.date := dmy(End.date)]

virus = fread(datapath("virusprev_nhs_regions_20201219200923weightedREACTonly.csv"))
virus[, Start.date := dmy(Start.date)]
virus[, End.date := dmy(End.date)]
virus = cbind(virus, virus[, skew_normal_solve(pct(Central.estimate), pct(Lower.bound), pct(Upper.bound), Data.source == "REACT-2")]);

sero[, pid := match(NHS.region, nhs_regions) - 1]
virus[, pid := match(NHS.region, nhs_regions) - 1]


# Add England to deaths series
ld = rbind(ld,
    ld[!name %in% c("Northern Ireland", "Scotland", "Wales"), .(N = sum(N), name = "England", pid = 1), by = date]
)
ld[, date := as.Date(date)]

# Add England to sitrep series
sitreps = rbind(sitreps,
    sitreps[!name %in% c("Northern Ireland", "Scotland", "Wales"),
        .(n_in_itu = sum(n_in_itu, na.rm = T), n_in_all_beds = sum(n_in_all_beds, na.rm = T), n_admitted_diagnosed = sum(n_admitted_diagnosed, na.rm = T),
            name = "England", pid = 1),
        by = date]
)

# # Variant data, add England
# variant = fread(datapath("var2-2020-12-21.csv"))
# variant = rbind(variant, 
#     variant[!nhs_name %in% c("Northern Ireland", "Scotland", "Wales"),
#         .(all = sum(all, na.rm = T), var2 = sum(var2, na.rm = T), nhs_name = "England"),
#         by = sample_date], fill = TRUE
# )
# variant[, pid := match(nhs_name, nhs_regions) - 1]
# variant[, sample_date := as.Date(sample_date)]

# SGTF data, add England
sgtf = fread(datapath("sgtf-2021-01-15.csv"))
sgtf = rbind(sgtf, 
    sgtf[!nhs_name %in% c("Northern Ireland", "Scotland", "Wales"),
        .(sgtf = sum(sgtf, na.rm = T), other = sum(other, na.rm = T), nhs_name = "England"),
        by = date], fill = TRUE
)
sgtf[, pid := match(nhs_name, nhs_regions) - 1]
sgtf[, date := as.Date(date)]


qsave(list(ld, sitreps, virus, sero, sgtf), datapath("processed-data-2021-01-15.qs"))
