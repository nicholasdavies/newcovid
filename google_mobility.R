library(stringr)
library(zoo)
library(data.table)


# set up covidm
uk_covid_data_path = "./fitting_data/";
datapath = function(x) paste0(uk_covid_data_path, x)

cm_path = "./covidm_for_fitting/";
cm_force_rebuild = F;
cm_build_verbose = T;
cm_version = 2;
source(paste0(cm_path, "/R/covidm.R"))
popUK = readRDS(datapath("popNHS.rds"));
matricesUK = readRDS(datapath("matricesNHS.rds"));

# Google visits data
# We don't upload full Google Mobility Report files to the repo because they are 200+ MB, 
# but if you download the "global CSV" from https://www.google.com/covid19/mobility/, that's it.
googmo = fread("~/Dropbox/uk_covid_data/fitting/data/Global_Mobility_Report-2021-01-14.csv");
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
xs = x[, .(date, t, value = rollmean(value, 7, fill = NA)), by = .(region_name, GPL_TYPE)]
xs = xs[!is.na(value)]

# Plot to check
ggplot(xs) +
    geom_line(aes(x = date, y = value, colour = region_name, group = region_name)) +
    geom_hline(aes(yintercept = 1), linetype = "22") +
    geom_vline(aes(xintercept = ymd("2020-12-25")), size = 0.2) +
    facet_wrap(~GPL_TYPE) +
    ylim(0, NA)

if (any(is.na(xs$value))) {
    stop("NAs remain.")
}

# Reformat
y = cbind(
    xs[GPL_TYPE == "grocery_and_pharmacy",           .(grocpharm = value), keyby = .(region_name, t)],
    parks = xs[GPL_TYPE == "parks",                  .(parks = value), keyby = .(region_name, t)][, parks],
    residential = xs[GPL_TYPE == "residential",      .(residential = value), keyby = .(region_name, t)][, residential],
    retrec = xs[GPL_TYPE == "retail_and_recreation", .(retrec = value), keyby = .(region_name, t)][, retrec],
    transit = xs[GPL_TYPE == "transit_stations",     .(transit = value), keyby = .(region_name, t)][, transit],
    workplace = xs[GPL_TYPE == "workplaces",         .(workplaces = value), keyby = .(region_name, t)][, workplaces]
);

# Fill in from Jan 1 2020
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

# Add school terms: England
# school_close =  c("2020-2-16", "2020-4-05", "2020-5-24", "2020-7-22", "2020-10-25", "2020-12-20", "2021-02-14", "2021-04-01", "2021-05-30", "2021-07-25");
# school_reopen = c("2020-2-22", "2020-4-18", "2020-5-30", "2020-9-01", "2020-10-31", "2021-01-02", "2021-02-20", "2021-04-17", "2021-06-05", "2021-09-01");
school_close =  c("2020-2-16", "2020-3-22", "2020-10-25", "2020-12-20", "2021-02-14", "2021-04-01", "2021-05-30", "2021-07-25");
school_reopen = c("2020-2-22", "2020-9-01", "2020-10-31", "2021-01-02", "2021-02-20", "2021-04-17", "2021-06-05", "2021-09-01");
stop("See comment above.")
school_c = as.numeric(ymd(school_close) - ymd("2020-01-01"))
school_r = as.numeric(ymd(school_reopen) - ymd("2020-01-01"))

days = y[, unique(t)]

n_closures = rowSums(matrix(days, ncol = length(school_c), nrow = length(days), byrow = FALSE) >= 
        matrix(school_c, ncol = length(school_c), nrow = length(days), byrow = TRUE));
n_reopenings = rowSums(matrix(days, ncol = length(school_r), nrow = length(days), byrow = FALSE) > 
        matrix(school_r, ncol = length(school_r), nrow = length(days), byrow = TRUE));
y_school = data.table(t = days, school = 1 - (n_closures - n_reopenings))
y = merge(y, y_school, by = "t")
y[, date := ymd("2020-01-01") + t]

ggplot(y) + geom_line(aes(x = date, y = parks, colour = region_name, group = region_name))
ggplot(y[date %between% c("2020-12-01", "2021-01-15")]) + geom_line(aes(x = date, y = workplace, colour = region_name, group = region_name))
ggplot(y) + geom_line(aes(x = date, y = school, colour = region_name, group = region_name))
region_names = y[, sort(unique(region_name))]
y[, pop := match(region_name, region_names) - 1]
y = y[order(pop, t)]



make_schedule = function(y)
{
    # Create schedule from Google data and school closures
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

    for (r in 1:nrow(y))
    {
        schedule[[y[r, pop]+1]]$values = c(schedule[[y[r, pop]+1]]$values,
            list(y[r, c(residential, workplace, grocpharm, retrec, transit, school)]));
        schedule[[y[r, pop]+1]]$times = c(schedule[[y[r, pop]+1]]$times, y[r, t])
    }

    return (schedule)
}

schedule = make_schedule(y);
saveRDS(schedule, "~/Documents/newcovid/fitting_data/schedule3-2021-01-14.rds")




# Alternative schedules for projections

alt_y = function(ref1_start, ref1_end, ref2_start, ref2_end, dhold, d0, d1_start, d1_end, d2_start, school_open_date)
{
    # Reformat
    y = cbind(
        xs[GPL_TYPE == "grocery_and_pharmacy",           .(grocpharm = value), keyby = .(region_name, t)],
        parks = xs[GPL_TYPE == "parks",                  .(parks = value), keyby = .(region_name, t)][, parks],
        residential = xs[GPL_TYPE == "residential",      .(residential = value), keyby = .(region_name, t)][, residential],
        retrec = xs[GPL_TYPE == "retail_and_recreation", .(retrec = value), keyby = .(region_name, t)][, retrec],
        transit = xs[GPL_TYPE == "transit_stations",     .(transit = value), keyby = .(region_name, t)][, transit],
        workplace = xs[GPL_TYPE == "workplaces",         .(workplaces = value), keyby = .(region_name, t)][, workplaces]
    );

    # Fill in from Jan 1 2020
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
    
    y[, date := ymd("2020-01-01") + t]
    
    # Bring in contact changes
    measured1 = y[date >= ref1_start & date <= ref1_end, lapply(.SD, mean), 
        .SDcols = c("grocpharm", "parks", "residential", "retrec", "transit", "workplace"), by = region_name]
    measured2 = y[date >= ref2_start & date <= ref2_end, lapply(.SD, mean), 
        .SDcols = c("grocpharm", "parks", "residential", "retrec", "transit", "workplace"), by = region_name]
    
    m1 = lapply(split(measured1, by = "region_name"), function(x) as.matrix(x[, 2:7]))
    m2 = lapply(split(measured2, by = "region_name"), function(x) as.matrix(x[, 2:7]))
    
    yy = split(y, by = "region_name")
    ym = lapply(yy, function(x) as.matrix(x[, 3:8]))
    
    rhold = as.numeric(ymd(dhold) - ymd("2019-12-31"))
    r0 = as.numeric(ymd(d0) - ymd("2019-12-31"))
    r1_start = as.numeric(ymd(d1_start) - ymd("2019-12-31"))
    r1_end =   as.numeric(ymd(d1_end) - ymd("2019-12-31"))
    r2_start = as.numeric(ymd(d2_start) - ymd("2019-12-31"))
    
    for (reg in names(ym)) {
        # Leading up to r0
        ym[[reg]][rhold:r0, ] = 
            ym[[reg]][rep(rhold, r0 - rhold + 1), ]
        
        # Leading up to r1_start
        ym[[reg]][r0:r1_start, ] = 
            ym[[reg]][rep(r0, r1_start - r0 + 1), ] * seq(1, 0, length.out = r1_start - r0 + 1) +
            matrix(m1[[reg]], nrow = r1_start - r0 + 1, ncol = 6, byrow = T) * seq(0, 1, length.out = r1_start - r0 + 1)
        
        # r1_start to r1_end
        ym[[reg]][r1_start:r1_end, ] =
            matrix(m1[[reg]], nrow = r1_end - r1_start + 1, ncol = 6, byrow = T)
        
        # r1_end to r2_start
        ym[[reg]][r1_end:r2_start, ] = 
            matrix(m1[[reg]], nrow = r2_start - r1_end + 1, ncol = 6, byrow = T) * seq(1, 0, length.out = r2_start - r1_end + 1) +
            matrix(m2[[reg]], nrow = r2_start - r1_end + 1, ncol = 6, byrow = T) * seq(0, 1, length.out = r2_start - r1_end + 1)
        
        # r2_start to end
        rf = nrow(ym[[reg]]);
        ym[[reg]][r2_start:rf, ] =
            matrix(m2[[reg]], nrow = rf - r2_start + 1, ncol = 6, byrow = T)
        
        for (j in 1:6) {
            yy[[reg]][[j + 2]] = ym[[reg]][, j]
        }
    }
    
    y = rbindlist(yy)
    
    # Add school terms: England
    school_close =  c("2020-2-16", "2020-3-22", "2020-10-25", "2020-12-20", "2021-02-14", "2021-04-01", "2021-05-30", "2021-07-25");
    school_reopen = c("2020-2-22", "2020-9-01", "2020-10-31", "2021-01-02", "2021-02-20", "2021-04-17", "2021-06-05", "2021-09-01");
    
    school_c = as.numeric(ymd(school_close) - ymd("2020-01-01"))
    school_r = as.numeric(ymd(school_reopen) - ymd("2020-01-01"))

    days = y[, unique(t)]

    n_closures = rowSums(matrix(days, ncol = length(school_c), nrow = length(days), byrow = FALSE) >= 
            matrix(school_c, ncol = length(school_c), nrow = length(days), byrow = TRUE));
    n_reopenings = rowSums(matrix(days, ncol = length(school_r), nrow = length(days), byrow = FALSE) > 
            matrix(school_r, ncol = length(school_r), nrow = length(days), byrow = TRUE));
    y_school = data.table(t = days, school = 1 - (n_closures - n_reopenings))
    y = merge(y, y_school, by = "t")
    
    # Adjust school opening date
    y[date >= "2020-12-25" & date < school_open_date, school := 0]

    region_names = y[, sort(unique(region_name))]
    y[, pop := match(region_name, region_names) - 1]
    y = y[order(pop, t)]
    
    return (y)
}

# Create alternative series for projections
# 1) return to mid-October (12-18 October)
# 2) lockdown from Christmas (9-15 Nov, Lockdown 2)
# 3) same but with school closures

y1 = alt_y(ref1_start = "2020-10-05", ref1_end = "2020-10-11", ref2_start = "2020-10-05", ref2_end = "2020-10-11",
    dhold = "2020-12-16", d0 = "2021-01-03", d1_start = "2021-01-07", d1_end = "2021-02-28", d2_start = "2021-04-01", school_open_date = "2021-01-02")

y2 = alt_y(ref1_start = "2020-11-09", ref1_end = "2020-11-15", ref2_start = "2020-10-05", ref2_end = "2020-10-11",
    dhold = "2020-12-16", d0 = "2021-01-03", d1_start = "2021-01-07", d1_end = "2021-02-28", d2_start = "2021-04-01", school_open_date = "2021-01-02")

y3 = alt_y(ref1_start = "2020-11-09", ref1_end = "2020-11-15", ref2_start = "2020-10-05", ref2_end = "2020-10-11",
    dhold = "2020-12-16", d0 = "2021-01-03", d1_start = "2021-01-07", d1_end = "2021-02-28", d2_start = "2021-04-01", school_open_date = "2021-02-19")

y4 = alt_y(ref1_start = "2020-04-06", ref1_end = "2020-04-12", ref2_start = "2020-10-05", ref2_end = "2020-10-11",
    dhold = "2020-12-16", d0 = "2021-01-03", d1_start = "2021-01-07", d1_end = "2021-02-28", d2_start = "2021-04-01", school_open_date = "2021-01-02")

y5 = alt_y(ref1_start = "2020-04-06", ref1_end = "2020-04-12", ref2_start = "2020-10-05", ref2_end = "2020-10-11",
    dhold = "2020-12-16", d0 = "2021-01-03", d1_start = "2021-01-07", d1_end = "2021-02-28", d2_start = "2021-04-01", school_open_date = "2021-02-19")

plot_y = function(y)
{
    melted = melt(y, id.vars = c("t", "region_name", "date", "pop"))
    ggplot(melted[date <= "2021-06-30" & !variable %in% c("parks", "residential") & 
            !region_name %in% c("Northern Ireland", "Scotland", "Wales", "United Kingdom", "England")]) +
        geom_line(aes(x = date, y = value, colour = region_name)) +
        facet_wrap(~variable, scales = "free") +
        theme(legend.position = c(0.7, 0.2)) +
        scale_x_date(date_labels = "%b %y")
}

theme_set(theme_cowplot(font_size = 10) + theme(strip.background = element_blank()))


py1 = plot_y(y1) + labs(x = "Date", y = "Mobility index", colour = "NHS region", title = "Moderate stringency (October 2020)")
py2 = plot_y(y2) + labs(x = "Date", y = "Mobility index", colour = "NHS region", title = "High stringency (November 2020) with schools open")
py3 = plot_y(y3) + labs(x = "Date", y = "Mobility index", colour = "NHS region", title = "High stringency (November 2020) with schools closed")
py4 = plot_y(y5) + labs(x = "Date", y = "Mobility index", colour = "NHS region", title = "Very high stringency (March 2020)")
cowplot::plot_grid(py1, py2, py3, py4, nrow = 2, labels = LETTERS, label_size = 10)
ggsave("./output/scenarios.pdf", width = 35, height = 25, units = "cm", useDingbats = TRUE)
ggsave("./output/scenarios.png", width = 35, height = 25, units = "cm")

plot_y(y1) + labs(x = "Date", y = "Mobility index", title = "Scenario 1")
ggsave("./output/scenario1.pdf", width = 15, height = 10, units = "cm", useDingbats = TRUE)
plot_y(y2) + labs(x = "Date", y = "Mobility index", title = "Scenario 2")
ggsave("./output/scenario2.pdf", width = 15, height = 10, units = "cm", useDingbats = TRUE)
plot_y(y3) + labs(x = "Date", y = "Mobility index", title = "Scenario 3")
ggsave("./output/scenario3.pdf", width = 15, height = 10, units = "cm", useDingbats = TRUE)
plot_y(y4) + labs(x = "Date", y = "Mobility index", title = "Scenario 4")
ggsave("./output/scenario4.pdf", width = 15, height = 10, units = "cm", useDingbats = TRUE)
plot_y(y5) + labs(x = "Date", y = "Mobility index", title = "Scenario 5")
ggsave("./output/scenario5.pdf", width = 15, height = 10, units = "cm", useDingbats = TRUE)


sc1 = make_schedule(y1);
sc2 = make_schedule(y2);
sc3 = make_schedule(y3);
sc4 = make_schedule(y4);
sc5 = make_schedule(y5);
saveRDS(sc1, "~/Documents/newcovid/fitting_data/scenario1.rds")
saveRDS(sc2, "~/Documents/newcovid/fitting_data/scenario2.rds")
saveRDS(sc3, "~/Documents/newcovid/fitting_data/scenario3.rds")
saveRDS(sc4, "~/Documents/newcovid/fitting_data/scenario4.rds")
saveRDS(sc5, "~/Documents/newcovid/fitting_data/scenario5.rds")

ggplot(y) + geom_line(aes(x = date, y = parks, colour = region_name, group = region_name))
ggplot(y[date %between% c("2020-12-01", "2021-01-15")]) + geom_line(aes(x = date, y = workplace, colour = region_name, group = region_name))
ggplot(y) + geom_line(aes(x = date, y = school, colour = region_name, group = region_name))




