library(stringr)
library(zoo)
library(data.table)

# Google visits data
# We don't upload full Google Mobility Report files to the repo because they are 200+ MB, 
# but if you download the "global CSV" from https://www.google.com/covid19/mobility/, that's it.
googmo = fread("~/Dropbox/uk_covid_data/fitting/data/Global_Mobility_Report-2021-01-05.csv");
googmo = googmo[country_region_code == "GB" & sub_region_1 != "" & sub_region_2 == "" & metro_area == ""];

# Melt mobility data
googmo = melt(googmo, id.vars = 1:8, variable.name = "GPL_TYPE", value.name = "value")
googmo[, GPL_TYPE := str_replace(GPL_TYPE, "_percent_change_from_baseline", "")]

# Match GSS with Google placenames
process_name = function(x)
{
    x = str_replace_all(x, "Council", "");
    x = str_replace_all(x, "Principle Area", "");
    x = str_replace_all(x, "Principal Area", "");
    x = str_replace_all(x, "County of", "");
    x = str_replace_all(x, "City of", "");
    x = str_replace_all(x, "City", "");
    x = str_replace_all(x, "County Borough", "");
    x = str_replace_all(x, "Borough of", "");
    x = str_replace_all(x, "Greater", "");
    x = str_replace_all(x, "Islands", "");
    x = str_replace_all(x, ",", "");
    x = str_replace_all(x, "  ", " ");
    x = str_trim(x);
    x
}

# Create table of all GSS names relevant for us
wards = fread("~/Dropbox/uk_covid_data/fitting/data/Ward_to_Local_Authority_District_to_County_to_Region_to_Country_December_2019_Lookup_in_United_Kingdom.csv")

placematch = rbind(
    unique(data.table(gss_code = wards$LAD19CD, gss_name = wards$LAD19NM)),
    unique(data.table(gss_code = wards$CTY19CD, gss_name = wards$CTY19NM)),
    unique(data.table(gss_code = wards$RGN19CD, gss_name = wards$RGN19NM))
)[gss_name != ""];

placematch[, match_name := process_name(gss_name)];

# Create table of Google placenames
googmo[sub_region_1 %like% "Eileanan", sub_region_1 := "Na h-Eileanan Siar"] # Rename places with spelling differences
googmo[sub_region_1 %like% "Rhondda",  sub_region_1 := "Rhondda Cynon Taf"]  # in Google mobility data resp. ONS data
gplacematch = data.table(goog_name = googmo[, unique(sub_region_1)]);
gplacematch[, match_name := process_name(goog_name)];

# Merge
match_table = merge(placematch, gplacematch, by = "match_name", all.y = T);

# Check for duplicates
resolve_duplicates = function(match_table, goog_nm, gss_nm, gss_cd = "")
{
    if (is.na(gss_nm)) {
        match_table[!(goog_name == goog_nm & (gss_code != gss_cd))]
    } else {
        match_table[!(goog_name == goog_nm & (gss_name != gss_nm))]
    }
}
match_table = resolve_duplicates(match_table, "Greater London", "London") # London (Google) refers to Greater London (ONS), not City of London (ONS)
match_table = resolve_duplicates(match_table, "Greater Manchester", "Greater Manchester") # Greater Manchester (Google) refers to Greater Manchester (ONS), not Manchester (ONS)
match_table = resolve_duplicates(match_table, "West Midlands", NA, "E11000005") # West Midlands (Google) refers to the metropolitan county, not the metropolitan district

if (nrow(match_table[duplicated(goog_name)]) > 0) {
    stop("Duplicates in match_table: ", paste(match_table[, goog_name[duplicated(goog_name)]], collapse = ", "))
}

# Amalgamate by NHS England region / DA
lads = fread("~/Dropbox/uk_covid_data/fitting/data/Lower_Layer_Super_Output_Area_2011_to_Clinical_Commissioning_Group_to_Local_Authority_District_April_2020_Lookup_in_England.csv");
nhsers = fread("~/Dropbox/uk_covid_data/fitting/data/Clinical_Commissioning_Group_to_STP_and_NHS_England_Region_April_2020_Lookup_in_England.csv");

# Match LAD to CCG to NHS region
lad_to_ccg = unique(lads[, .(gss_code = LAD20CD, CCG20CD)]);
ccg_to_nhs = unique(nhsers[, .(CCG20CD, region_name = NHSER20NM)]);

match_table = merge(match_table, lad_to_ccg, by = "gss_code", all.x = T)
match_table = merge(match_table, ccg_to_nhs, by = "CCG20CD", all.x = T)

# Fill in DAs and some regions
match_table[gss_code %like% "^N", region_name := "Northern Ireland"]
match_table[gss_code %like% "^S", region_name := "Scotland"]
match_table[gss_code %like% "^W", region_name := "Wales"]
match_table[gss_name == "London", region_name := "London"]
match_table[gss_name == "Buckinghamshire", region_name := "South East"];

# Fill in missing counties E10 and metropolitan counties E11
ladcounty = unique(wards[, .(LAD19CD, LAD19NM, CTY19CD, CTY19NM)])
ladcounty = merge(ladcounty, unique(lads[, .(LAD19CD = LAD20CD, CCG20CD)]), all.x = T)
ladcounty = merge(ladcounty, unique(nhsers[, .(CCG20CD, region_name = NHSER20NM)]), by = "CCG20CD", all.x = T)
county_to_nhs = unique(ladcounty[CTY19CD != "", .(gss_code = CTY19CD, region_name)])
match_table = merge(match_table, unique(ladcounty[, .(gss_code = CTY19CD, region_name)]), by = "gss_code", all.x = T);

# Take region name x or y
if (nrow(match_table[!is.na(region_name.x) & !is.na(region_name.y)]) > 0) {
    stop("Match table has both region_name.x and region_name.y.");
}
region_names = unique(c(match_table[!is.na(region_name.x), region_name.x], match_table[!is.na(region_name.y), region_name.y]));
either = function(x, y) ifelse(is.na(x) & is.na(y), NULL, ifelse(is.na(x), y, x))
match_table[, region_name := either(region_name.x, region_name.y)];
match_table[, region_name.x := NULL];
match_table[, region_name.y := NULL];

# Remove entries that are duplicates because of two CCGs but which are in the same region
match_table[, CCG20CD := NULL];
match_table = unique(match_table);

# Count number of regions each google name matches to
match_table[, mappings := .N, by = goog_name];

# Get population info
match_table = merge(match_table, cm_structure_UK[, .(gss_code = Code, population = `All ages`)], by = "gss_code", all.x = T);
if (nrow(match_table[is.na(population)]) > 0) {
    stop("Missing population information.")
}

# Bring together and amalgamate by NHS region
googen2 = merge(googmo, match_table, by.x = "sub_region_1", by.y = "goog_name", allow.cartesian = T);

# Fill missing values in data sets, from median for date within region
googen2[, value := as.numeric(value)]
googen2[, value := na.aggregate(value, FUN = median), by = .(date, GPL_TYPE, region_name)]

# Amalgamate by NHS region
x = googen2[, .(value = weighted.mean(value, population/mappings, na.rm = T)), keyby = .(region_name, GPL_TYPE, date)]

# Add England and United Kingdom
x = rbind(x,
    googen2[, .(region_name = "United Kingdom", value = weighted.mean(value, population/mappings, na.rm = T)), keyby = .(GPL_TYPE, date)],
    googen2[!region_name %in% c("Northern Ireland", "Scotland", "Wales"), 
        .(region_name = "England", value = weighted.mean(value, population/mappings, na.rm = T)), keyby = .(GPL_TYPE, date)],
    use.names = TRUE
)

# Transform data and add date
x[, date := ymd(as.character(date))]
x[, value := 1 + value/100]
x[, t := as.numeric(date - ymd("2020-01-01"))]

# Fill in NAs within time series
x[, value := na.approx(value, na.rm = FALSE), by = .(region_name, GPL_TYPE)]

# Ensmoothen
x = x[, .(date, t, value = rollmean(value, 7, fill = NA)), by = .(region_name, GPL_TYPE)]
x = x[!is.na(value)]

# Plot to check
ggplot(x) +
    geom_line(aes(x = date, y = value, colour = region_name, group = region_name)) +
    geom_hline(aes(yintercept = 1), linetype = "22") +
    geom_vline(aes(xintercept = ymd("2020-12-25")), size = 0.2) +
    facet_wrap(~GPL_TYPE) +
    ylim(0, NA)

if (any(is.na(x$value))) {
    stop("NAs remain.")
}

# Reformat
y = cbind(
    x[GPL_TYPE == "grocery_and_pharmacy",           .(grocpharm = value), keyby = .(region_name, t)],
    parks = x[GPL_TYPE == "parks",                  .(parks = value), keyby = .(region_name, t)][, parks],
    residential = x[GPL_TYPE == "residential",      .(residential = value), keyby = .(region_name, t)][, residential],
    retrec = x[GPL_TYPE == "retail_and_recreation", .(retrec = value), keyby = .(region_name, t)][, retrec],
    transit = x[GPL_TYPE == "transit_stations",     .(transit = value), keyby = .(region_name, t)][, transit],
    workplace = x[GPL_TYPE == "workplaces",         .(workplaces = value), keyby = .(region_name, t)][, workplaces]
);

# Fill in from Jan 1
y = rbind(y, 
    data.table(region_name = rep(unique(y$region_name), each = 48), t = rep(0:47, length(unique(y$region_name))),
        grocpharm = 1, parks = 1, residential = 1, retrec = 1, transit = 1, workplace = 1)
)
y = y[order(region_name, t)]

# Carry forward last day
last_day = y[t == max(t)]
last_fortnight = rbindlist(rep(list(last_day), 14))
last_fortnight[, t := t + 0:13, by = region_name]

for (added_fortnights in 1:45) {
    y = rbind(y, last_fortnight);
    last_fortnight[, t := t + 14];
}

ggplot(y) + geom_line(aes(x = t, y = parks, colour = region_name, group = region_name))
ggplot(y) + geom_line(aes(x = t, y = workplace, colour = region_name, group = region_name))
region_names = y[, sort(unique(region_name))]
y[, pop := match(region_name, region_names) - 1]
y = y[order(pop, t)]

make_schedule = function(y, scenario)
{
    # Create schedule from Google data
    schedule = list()
    schedule[[1]]  = list(parameter = "contact", pops = 0,  mode = "assign", values = list(), times = numeric(0))
    schedule[[2]]  = list(parameter = "contact", pops = 1,  mode = "assign", values = list(), times = numeric(0))
    schedule[[3]]  = list(parameter = "contact", pops = 2,  mode = "assign", values = list(), times = numeric(0))
    schedule[[4]]  = list(parameter = "contact", pops = 3,  mode = "assign", values = list(), times = numeric(0))
    schedule[[5]]  = list(parameter = "contact", pops = 4,  mode = "assign", values = list(), times = numeric(0))
    schedule[[6]]  = list(parameter = "contact", pops = 5,  mode = "assign", values = list(), times = numeric(0))
    schedule[[7]]  = list(parameter = "contact", pops = 6,  mode = "assign", values = list(), times = numeric(0))
    schedule[[8]]  = list(parameter = "contact", pops = 7,  mode = "assign", values = list(), times = numeric(0))
    schedule[[9]]  = list(parameter = "contact", pops = 8,  mode = "assign", values = list(), times = numeric(0))
    schedule[[10]] = list(parameter = "contact", pops = 9,  mode = "assign", values = list(), times = numeric(0))
    schedule[[11]] = list(parameter = "contact", pops = 10, mode = "assign", values = list(), times = numeric(0))
    schedule[[12]] = list(parameter = "contact", pops = 11, mode = "assign", values = list(), times = numeric(0))
    schedule[[13]] = list(parameter = "fIs",     pops = 0,  mode = "multiply", values = list(), times = numeric(0))
    schedule[[14]] = list(parameter = "fIs",     pops = 1,  mode = "multiply", values = list(), times = numeric(0))
    schedule[[15]] = list(parameter = "fIs",     pops = 2,  mode = "multiply", values = list(), times = numeric(0))
    schedule[[16]] = list(parameter = "fIs",     pops = 3,  mode = "multiply", values = list(), times = numeric(0))
    schedule[[17]] = list(parameter = "fIs",     pops = 4,  mode = "multiply", values = list(), times = numeric(0))
    schedule[[18]] = list(parameter = "fIs",     pops = 5,  mode = "multiply", values = list(), times = numeric(0))
    schedule[[19]] = list(parameter = "fIs",     pops = 6,  mode = "multiply", values = list(), times = numeric(0))
    schedule[[20]] = list(parameter = "fIs",     pops = 7,  mode = "multiply", values = list(), times = numeric(0))
    schedule[[21]] = list(parameter = "fIs",     pops = 8,  mode = "multiply", values = list(), times = numeric(0))
    schedule[[22]] = list(parameter = "fIs",     pops = 9,  mode = "multiply", values = list(), times = numeric(0))
    schedule[[23]] = list(parameter = "fIs",     pops = 10, mode = "multiply", values = list(), times = numeric(0))
    schedule[[24]] = list(parameter = "fIs",     pops = 11, mode = "multiply", values = list(), times = numeric(0))

    if (scenario == 0) {
        W = rep(0, 713);
        S = rep(0, 713);
        L = rep(0, 713);
    }

    for (r in 1:nrow(y))
    {
        # schedule[[y[r, pop]+1]]$values = c(schedule[[y[r, pop]+1]]$values,
        #     list(y[r, c(0.5 + 0.5*transit, workplace + W[t], ifelse(t >= 81 & t <= 244, 0, 1) + S[t], 0.05*parks + 0.15*grocpharm + 0.4*retrec + 0.4*transit + L[t],
        #               c(0.5 + 0.5*transit, workplace       , ifelse(t >= 81 & t <= 244, 0, 1)       , 0.05*parks + 0.15*grocpharm + 0.4*retrec + 0.4*transit + L[t]/2) * 1)]));
        # schedule[[y[r, pop]+1]]$times = c(schedule[[y[r, pop]+1]]$times, y[r, t])
        # 
        # schedule[[y[r, pop]+11]]$values = c(schedule[[y[r, pop]+11]]$values,
        #     list(y[r, rep(0.5 + 0.5*transit, 16)]));
        # schedule[[y[r, pop]+11]]$times = c(schedule[[y[r, pop]+11]]$times, y[r, t])
        # below for schedule3
        schedule[[y[r, pop]+1]]$values = c(schedule[[y[r, pop]+1]]$values,
            list(y[r, c(residential, workplace, grocpharm, retrec, transit, 1, 1, 1)]));
        schedule[[y[r, pop]+1]]$times = c(schedule[[y[r, pop]+1]]$times, y[r, t])

        schedule[[y[r, pop]+13]]$values = c(schedule[[y[r, pop]+13]]$values,
            list(y[r, c(transit, rep(1, 15))]));
        schedule[[y[r, pop]+13]]$times = c(schedule[[y[r, pop]+13]]$times, y[r, t])
    }

    return (schedule)
}

schedule = make_schedule(y, 0);
saveRDS(schedule, "~/Documents/newcovid/fitting_data/schedule3-2021-01-05.rds")

