library(ggplot2)
library(data.table)
library(zoo)
library(ogwrangler)
library(lubridate)
CreateCache()

# ONS Postcode Directory
onspd = fread("~/Downloads/ONSPD_AUG_2020_UK/Data/ONSPD_AUG_2020_UK.csv");
onspd = onspd[is.na(doterm), .(pcd, nhser, oslaua, lat, long)]
setkey(onspd, lat, long)

# Find ONSPD entry closest to a given latitude and longitude.
latlon = function(latitude, longitude, tol = 1)
{
    onspd2 = onspd[lat %between% c(latitude - tol, latitude + tol) & long %between% c(longitude - tol, longitude + tol)];
    onspd2[, dist := (lat - latitude)^2 + (long - longitude)^2];
    return (onspd2[which.min(dist)])
}

variance_beta = function(a, b) (a * b) / ((a + b)^2 * (a + b + 1))
mean_beta = function(a, b) a / (a + b)
variance_betas = function(w, a, b) sum(w^2 * mapply(variance_beta, a, b))
mean_betas = function(w, a, b) sum(w * mapply(mean_beta, a, b))

# mode of a beta dist
beta_mode = function(a, b)
{
    ifelse(a > 1 & b > 1, (a - 1) / (a + b - 2),
        ifelse(a < b, 0, 1))
}

# For identifying HDI of a beta distribution, when mode is not 0 or 1
p_range_beta = function(d, a, b)
{
    min = uniroot(function(x) dbeta(x, a, b) - d, lower = 0, upper = mode)$root;
    max = uniroot(function(x) dbeta(x, a, b) - d, lower = mode, upper = 1)$root;
    
    # return probability between min and max
    pbeta(max, a, b) - pbeta(min, a, b)
}

# Identify HDI of a beta distribution, when mode is not 0 or 1
beta_hdi_find = function(m, a, b)
{
    mode = beta_mode(a, b);
    if (mode == 0 || mode == 1)
        stop("beta distribution must have 0 < mode < 1.")
    peak = dbeta(mode, a, b);

    d = log(uniroot(function(x) p_range_beta(x, a, b) - m, lower = 0, upper = peak)$root);

    min = uniroot(function(x) dbeta(x, a, b, log = TRUE) - d, lower = 0, upper = mode)$root;
    max = uniroot(function(x) dbeta(x, a, b, log = TRUE) - d, lower = mode, upper = 1)$root;

    return (c(min = min, max = max))
}

# HDI of a beta dist
beta_hdi = function(a, b, m = 0.95)
{
    if (a == 1 && b == 1) {
        return (c(0.5 - m/2, 0.5 + m/2))
    }
    
    mode = beta_mode(a, b)
    
    if (mode == 0) {
        min = 0;
        max = qbeta(m, a, b);
    } else if (mode == 1) {
        min = qbeta(1 - m, a, b);
        max = 1;
    } else {
        return (beta_hdi_find(m, a, b))
    }
    
    return (c(min = min, max = max))
}

# Approximate weighted mean of several beta distributions
approx_mean_betas = function(w, a, b)
{
    mean = mean_betas(w, a, b);
    variance = variance_betas(w, a, b);
    alpha = mean * ((mean * (1 - mean)) / variance - 1)
    beta = (1 - mean) * ((mean * (1 - mean)) / variance - 1)
    return (c(shape1 = alpha, shape2 = beta))
}

# Get new lineage data
nl = function()
{
    #newlin = fread("./data/cog_metadata_microreact_public-2020-12-18.csv")
    newlin = fread("./data/cog_metadata_microreact_public-2020-12-22.csv")#
    newlin = newlin[country == "UK"]
    newlin[, site := .GRP, by = .(longitude, latitude)]
    newlin[, B117 := lineage == "B.1.1.7"]
    newlin[, N_B117 := sum(B117), by = site]
    newlin[, var2 := B117 & n501y == "Y"]
    newlin[, N_var2 := sum(var2), by = site]
    newlin[, spec := paste0(d614g, n439k, p323l, a222v, y453f, n501y, del_21765_6)]
    
    # Assign localities
    sites = newlin[, unique(site)]
    for (s in sites) {
        cat(".");
        latitude = newlin[site == s, latitude[1]];
        longitude = newlin[site == s, longitude[1]];
        loc = latlon(latitude, longitude);
        if (nrow(loc) > 0) {
            newlin[site == s, pcd := loc$pcd];
            newlin[site == s, lad := loc$oslaua];
            newlin[site == s, nhs := loc$nhser];
        }
    }
    
    # Omit Gibraltar and unknown sites
    newlin = newlin[!is.na(nhs)]
    
    # NHS regions
    newlin[nhs %like% "E", nhs_name := ogwhat(nhs)]
    newlin[nhs %like% "N", nhs_name := "Northern Ireland"]
    newlin[nhs %like% "S", nhs_name := "Scotland"]
    newlin[nhs %like% "W", nhs_name := "Wales"]
    
    # Reassign site id to remove missing sites
    newlin[, site := .GRP, by = .(longitude, latitude)]
    
    return (newlin[])
}

newlin = nl()


# View sites
ggplot(unique(newlin[, .(latitude, longitude, nhs_name)])) + 
    geom_point(aes(x = longitude, y = latitude, colour = nhs_name))

# Build site-frequency-corrected variant frequency table
sitefreq = newlin[, .N, by = .(nhs_name, site)]
sitefreq[, freq := N / sum(N), by = .(nhs_name)]
date_min = newlin[, min(sample_date)]
date_max = newlin[, max(sample_date)]
ndate = as.numeric(date_max - date_min) + 1
nsite = sitefreq[, uniqueN(site)]
prior_a = 0.5
prior_b = 0.5

varfreq = data.table(site = rep(1:nsite, each = ndate))
varfreq[, date := rep(date_min + (0:(ndate - 1)), nsite)]
varfreq = merge(varfreq, 
    newlin[, .(var1 = sum(!var2), var2 = sum(var2)), keyby = .(date = sample_date, site)],
    by = c("date", "site"), all = TRUE)
varfreq[is.na(var1), var1 := 0]
varfreq[is.na(var2), var2 := 0]
varfreq = merge(varfreq, sitefreq[, .(site, nhs_name, sitefreq = freq)], by = "site")

# Find cutoff point
ggplot(varfreq[, sum(var1 > 0 | var2 > 0) / .N, by = .(date, nhs_name)]) + 
    geom_line(aes(date, V1, colour = nhs_name)) +
    geom_vline(aes(xintercept = ymd("2020-11-25"))) +
    facet_wrap(~nhs_name)

data_site = newlin[, .(all = .N, var2 = sum(var2)), keyby = .(site, sample_date, nhs_name)]
data_site[, site2 := as.numeric(match(site, unique(site))), by = nhs_name]
data_site[, site2 := site2 / max(site2), by = nhs_name]
ggplot(data_site[sample_date > "2020-10-01"]) +
    geom_col(aes(x = sample_date, y = all, fill = site2), colour = "black", size = 0.2, position = "fill") +
    facet_wrap(~nhs_name) +
    theme(legend.position = "none") +
    scale_fill_gradientn(colours = c("red", "yellow", "green", "blue", "violet")) +
    geom_vline(aes(xintercept = ymd("2020-12-01"))) +
    labs(x = "Sample date", y = "Proportion of COG-UK samples from each site")


ggplot(varfreq[, sum(var1 + var2), by = .(date, nhs_name)]) + 
    geom_line(aes(date, V1, colour = nhs_name)) +
    geom_vline(aes(xintercept = ymd("2020-11-25"))) +
    facet_wrap(~nhs_name, scales = "free")

v = varfreq[, .(var1 = sum(var1 * sitefreq) + 1, var2 = sum(var2 * sitefreq) + 0.5), by = .(nhs_name, date)]

for (i in 1:nrow(v)) {
    a = v[i, var2] # inverted because we are interested
    b = v[i, var1] # in the frequency of variant 2
    v[i, q_lo := qbeta(0.025, a, b)]
    v[i, q_hi := qbeta(0.975, a, b)]
    v[i, mode := beta_mode(v[i, var2], v[i, var1])];
}

ggplot(v) +
    geom_ribbon(aes(x = date, ymin = q_lo, ymax = q_hi), alpha = 0.2) +
    geom_line(aes(x = date, y = mode), alpha = 0.5) +
    facet_wrap(~nhs_name) + scale_y_log10()


data = newlin[, .(all = .N, var2 = sum(var2)), keyby = .(sample_date, nhs_name)]
ggplot(data[sample_date <= "2020-12-31"]) +
    geom_line(aes(x = sample_date, y = var2 / all, colour = nhs_name)) +
    facet_wrap(~nhs_name) +
    theme(legend.position = "none")

fwrite(data[sample_date <= "2020-12-01"], "./fitting_data/var2-2020-12-16.csv")

data_site = newlin[, .(all = .N, var2 = sum(var2)), keyby = .(site, sample_date)]
data_site
ggplot(data_site) +
    geom_line(aes(x = sample_date, y = var2 / all, colour = site)) +
    facet_wrap(~site) +
    theme(legend.position = "none") +
    geom_vline(aes(xintercept = ymd("2020-12-01")))
