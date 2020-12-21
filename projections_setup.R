library(data.table)
library(ggplot2)
library(lubridate)
library(here)
library(cowplot)
library(readxl)
library(sn)
library(qs)
library(stringr)
library(mgcv)
library(binom)

N_THREADS = 36
which_pops = c(3,1,9)

uk_covid_data_path = "./fitting_data/";
datapath = function(x) paste0(uk_covid_data_path, x)

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
nhs_regions = popUK[, unique(name)]

cm_populations = rbind(cm_populations[name != "United Kingdom"], popUK)
cm_matrices = c(cm_matrices, matricesUK)
source("./distribution_fit.R");
source("./spim_output.R");
source("./check_fit.R")


#
# Load components
#

# Health burden processes
source("./processes.R")

# CPP functions for covidm
source("./cpp_funcs.R")

# LOAD FITS
source("./latest_fit.R")

# Determine population sizes and ensure no undefined vaccine parameters
popsize = NULL
for (i in seq_along(parametersI)) {
    if (!is.null(parametersI[[i]])) {
        popsize = rbind(popsize,
            data.table(Geography = parametersI[[i]]$pop[[1]]$name, population_size = sum(parametersI[[i]]$pop[[1]]$size))
        )
    }
}
popsize = rbind(data.table(Geography = "England", population_size = sum(popsize$population_size)), popsize)
setkey(popsize, Geography)



# 2. Bringing in of interventions

# Extract first column from each of tables_list with scenario names
england_only = function(tables_list, names)
{
    t = tables_list[[1]][, 1]
    for (i in seq_along(tables_list)) {
        t = cbind(t, tables_list[[i]][, 2]);
        names(t)[i + 1] = names[i];    
    }
    
    return (t)
}

arrange_projection = function(proj, cumulative_deaths = FALSE, from_date = NULL, england = FALSE)
{
    ft = 0;
    if (!is.null(from_date)) {
        ft = as.numeric(ymd(from_date) - ymd("2020-01-01"))
    }
    
    if (england == FALSE) {
        w = proj[t >= ft, .(deaths = sum(death_o + death2_o), admissions = sum(hosp_undetected_o + hosp_undetected2_o),
            beds = sum(hosp_p + hosp2_p - hosp_undetected_p - hosp_undetected2_p), icu = sum(icu_p + icu2_p), 
            Rt = obs0[1], tier = obs0[2], cb = obs0[3]), keyby = .(run, t, population)]
    } else {
        w = proj[t >= ft, .(population = "England", deaths = sum(death_o + death2_o), admissions = sum(hosp_undetected_o + hosp_undetected2_o),
            beds = sum(hosp_p + hosp2_p - hosp_undetected_p - hosp_undetected2_p), icu = sum(icu_p + icu2_p), 
            Rt = mean(obs0[seq(1, .N, by = 16)]), tier = mean(obs0[seq(2, .N, by = 16)]), cb = mean(obs0[seq(3, .N, by = 16)])), keyby = .(run, t)]
    }
    
    if (cumulative_deaths) {
        w[, cum_deaths := cumsum(deaths), by = .(run, population)]
    }
    w = melt(w, id.vars = 1:3)
    w = w[, as.list(quantile(value, c(0.025, 0.5, 0.975))), by = .(population, variable, t)]
    w[variable == "deaths", variable := "Deaths"]
    w[variable == "cum_deaths", variable := "Cumulative deaths"]
    w[variable == "admissions", variable := "Admissions"]
    w[variable == "beds", variable := "Hospital beds"]
    w[variable == "icu", variable := "ICU beds"]
}

plot_projection = function(proj_list, proj_names, from_date, pal = "Accent")
{
    p = NULL
    for (i in seq_along(proj_list)) {
        pp = arrange_projection(proj_list[[i]]);
        pp[, name := proj_names[[i]]];
        p = rbind(p, pp);
    }
    p[, name := factor(name, unique(name))]

    plot = ggplot(p[!variable %in% c("tier", "cb") & ymd("2020-01-01") + t >= from_date]);
    
    cbs = p[variable == "cb", .(t = ymd("2020-01-01") + t[`50%` > 0.5]), by = .(population)][, 
        .(tmin = min(t), tmax = max(t)), by = population]
    if (nrow(cbs) > 0) {
        plot = plot +
            geom_rect(data = cbs, aes(xmin = tmin, xmax = tmax, ymin = -Inf, ymax = Inf), fill = "black", alpha = 0.1)
    }
    
    rline = data.table(variable = factor("Rt", levels(p$variable)), y = 1)
    
    plot +
        geom_ribbon(aes(x = ymd("2020-01-01") + t, ymin = `2.5%`, ymax = `97.5%`, fill = name), alpha = 0.5) +
        geom_hline(data = rline, aes(yintercept = y), size = 0.3) +
        geom_line(aes(x = ymd("2020-01-01") + t, y = `50%`, colour = name)) +
        facet_grid(variable ~ population, switch = "y", scales = "free") +
        cowplot::theme_cowplot(font_size = 11) + 
        theme(strip.background = element_blank(), strip.placement = "outside", 
            legend.position = "bottom",
            panel.background = element_rect(fill = "#f4f4f4"),
            panel.grid.major = element_line(colour = "#ffffff", size = 0.4)) +
        labs(x = NULL, y = NULL, colour = NULL, fill = NULL) +
        scale_x_date(date_breaks = "1 month", date_labels = "%b") +
        scale_color_brewer(aesthetics = c("colour", "fill"), palette = pal) +
        ylim(0, NA)
}

plot_cum_deaths = function(proj_list, proj_names, from_date, ystart = 26, ydiff = -1.2, titl = "Cumulative deaths")
{
    stop("Check for both virus strains.")
    p = NULL
    for (i in seq_along(proj_list)) {
        w = arrange_projection(proj_list[[i]], cumulative_deaths = TRUE, from_date = from_date, england = TRUE)
        p = rbind(p,
            cbind(w[variable == "Cumulative deaths" | variable == "cb"], name = proj_names[i])
        )
    }
    
    p[, name := factor(name, levels = unique(name))]
    
    cbs = p[variable == "cb", .(t = ymd("2020-01-01") + t[`50%` > 0.5]), by = .(population, name)][, 
        .(tmin = min(t), tmax = max(t)), by = .(population, name)]
    
    plot = ggplot(p[variable != "cb"]) +
        geom_ribbon(aes(x = ymd("2020-01-01") + t, ymin = `2.5%` / 1000, ymax = `97.5%` / 1000, fill = name), alpha = 0.5) +
        geom_line(aes(x = ymd("2020-01-01") + t, y = `50%` / 1000, colour = name), size = 0.2) +
        facet_wrap(~population, scales = "free") +
        theme(strip.background = element_blank(), strip.text = element_blank(), legend.position = c(0.02, 0.8),
            panel.background = element_rect(fill = "#f4f4f4"),
            panel.grid.major = element_line(colour = "#ffffff", size = 0.2)) +
        labs(x = NULL, y = "Cumulative deaths\n(thousands)", colour = NULL, fill = NULL, title = titl) +
        scale_x_date(date_breaks = "1 month", date_labels = "%b")
    
    if (nrow(cbs) > 0) {
        yy = ystart
        for (nm in cbs[, unique(name)]) {
            plot = plot +
                geom_linerange(data = cbs[name == nm], aes(xmin = tmin, xmax = tmax, colour = name), y = yy, show.legend = FALSE, linetype = "32", size = 0.2) +
                geom_linerange(data = cbs[name == nm], aes(x = tmin, colour = name), ymin = yy - (ydiff * 0.2), ymax = yy + (ydiff * 0.2), show.legend = FALSE, size = 0.2) +
                geom_linerange(data = cbs[name == nm], aes(x = tmax, colour = name), ymin = yy - (ydiff * 0.2), ymax = yy + (ydiff * 0.2), show.legend = FALSE, size = 0.2)
            yy = yy + ydiff
        }
    }

    return (plot)
}

plot_icu = function(proj_list, proj_names, from_date)
{
    stop("Check for both virus strains.")

    p = NULL
    for (i in seq_along(proj_list)) {
        w = arrange_projection(proj_list[[i]], cumulative_deaths = TRUE, from_date = from_date, england = TRUE)
        p = rbind(p,
            cbind(w[variable == "ICU beds"], name = proj_names[i])
        )
    }
    
    p[, name := factor(name, levels = unique(name))]
    
    plot = ggplot(p) +
        geom_ribbon(aes(x = ymd("2020-01-01") + t, ymin = `2.5%`, ymax = `97.5%`, fill = name), alpha = 0.5) +
        geom_line(aes(x = ymd("2020-01-01") + t, y = `50%`, colour = name), size = 0.2) +
        facet_wrap(~population, scales = "free") +
        theme(strip.background = element_blank(), strip.text = element_blank(), legend.position = "none",
            panel.background = element_rect(fill = "#f4f4f4"),
            panel.grid.major = element_line(colour = "#ffffff", size = 0.2)) +
        labs(x = NULL, y = "ICU beds occupied", colour = NULL, fill = NULL) +
        scale_x_date(date_breaks = "1 month", date_labels = "%b")

    return (plot)
}


summarize_projection = function(proj, from_date, popsize, scenario_output = NULL, to_date = NULL, wh = "after")
{
    stop("Check for both virus strains.")

    if (!is.null(to_date)) {
        proj = rlang::duplicate(proj)
        proj = proj[ymd("2020-01-01") + t <= to_date]
    }
    
    geo = c("England", proj[, unique(population)])
    
    # Build summaries
    totals = proj[, .(cases = sum(cases), infections = sum(infections_i), 
        deaths = sum(death_o), admissions = sum(hosp_undetected_o)), 
        by = .(population, run, when = ifelse(ymd("2020-01-01") + t < from_date, "before", "after"))]
    totals = rbind(totals[, .(population = "England", cases = sum(cases), infections = sum(infections),
        deaths = sum(deaths), admissions = sum(admissions)), by = .(run, when)], totals, fill = TRUE);
    
    prev_peaks = proj[ymd("2020-01-01") + t < from_date,
        .(all = sum(hosp_p - hosp_undetected_p), icu = sum(icu_p)),
        by = .(population, run, t)][,
        .(all = max(all), icu = max(icu)),
        by = .(population, run)][,
        .(peak_all = mean(all), peak_icu = mean(icu)), keyby = population]
    prev_peaks = rbind(prev_peaks[, .(population = "England", peak_all = sum(peak_all), peak_icu = sum(peak_icu))], prev_peaks, fill = TRUE);
    setkey(prev_peaks, population)
    
    if (wh == "before") {
        post_peaks = proj[ymd("2020-01-01") + t < from_date,
            .(all = sum(hosp_p - hosp_undetected_p), icu = sum(icu_p)),
            by = .(population, run, t)][,
            .(all = max(all), icu = max(icu)),
            keyby = .(population, run)]
    } else {
        post_peaks = proj[ymd("2020-01-01") + t >= from_date,
            .(all = sum(hosp_p - hosp_undetected_p), icu = sum(icu_p)),
            by = .(population, run, t)][,
            .(all = max(all), icu = max(icu)),
            keyby = .(population, run)]
    }
    post_peaks = rbind(post_peaks[, .(population = "England", all = sum(all), icu = sum(icu)), by = run], post_peaks, fill = TRUE);
    setkey(post_peaks, population)
    
    peaks = post_peaks[,
        .(all = all, rel_all = all / prev_peaks[population, peak_all],
          icu = icu, rel_icu = icu / prev_peaks[population, peak_icu]),
        keyby = .(population, run)]

    hi = proj[ymd("2020-01-01") + t >= from_date,
        .(all = sum(hosp_p - hosp_undetected_p), icu = sum(icu_p)),
        keyby = .(population, run, t)][,
        .(all_hi = sum(all > prev_peaks[population, peak_all] * 0.5) / 7,
          icu_hi = sum(icu > prev_peaks[population, peak_icu] * 0.5) / 7), 
        keyby = .(population, run)]
    
    tier = proj[ymd("2020-01-01") + t >= from_date,
        .(tier = obs0[2], cb = obs0[3]), keyby = .(population, run, t)][,
        .(tier2 = sum(tier == 2 & cb == 0) / 7, tier3 = sum(tier == 3 & cb == 0) / 7, cb = sum(cb == 1) / 7),
        keyby = .(population, run)]
    
    hi2 = merge(hi, popsize, by.x = "population", by.y = "Geography")
    hi2 = hi2[, .(population = "England", all_hi = weighted.mean(all_hi, population_size), 
        icu_hi = weighted.mean(icu_hi, population_size)), by = run]
    hi = rbind(hi2, hi)
    
    tier2 = merge(tier, popsize, by.x = "population", by.y = "Geography")
    tier2 = tier2[, .(population = "England", tier2 = weighted.mean(tier2, population_size),
        tier3 = weighted.mean(tier3, population_size), cb = weighted.mean(cb, population_size)), by = run]
    tier = rbind(tier2, tier)

    # Factor
    totals[, population := factor(population, geo)]
    peaks[, population := factor(population, geo)]
    hi[, population := factor(population, geo)]
    tier[, population := factor(population, geo)]
    
    if (is.null(scenario_output)) { # Nice table output
        tab = rbind(
            totals[when == wh, .(indicator = "Deaths", value = niceq(deaths)), keyby = population],
            totals[when == wh, .(indicator = "Admissions", value = niceq(admissions)), keyby = population],
            peaks[, .(indicator = "Peak ICU requirement", value = niceq(icu)), keyby = population],
            peaks[, .(indicator = "Peak ICU (rel. W1)", value = nicepc(rel_icu)), keyby = population],
            hi[, .(indicator = "Weeks of high ICU occupancy", value = niceq(icu_hi)), keyby = population],
            tier[, .(indicator = "Weeks in Tier 2", value = niceq(tier2)), keyby = population],
            tier[, .(indicator = "Weeks in Tier 3", value = niceq(tier3)), keyby = population],
            tier[, .(indicator = "Weeks in lockdown", value = niceq(cb)), keyby = population]
        )
        return (dcast(tab, indicator ~ population))
    } else { # Data output
        tab = rbind(
            totals[when == "after", .(indicator = "Cases",      value = quantile(cases,      c(0.025, 0.5, 0.975)), q = c("lo", "mid", "hi")), keyby = population],
            totals[when == "after", .(indicator = "Infections", value = quantile(infections, c(0.025, 0.5, 0.975)), q = c("lo", "mid", "hi")), keyby = population],
            totals[when == "after", .(indicator = "Deaths",     value = quantile(deaths,     c(0.025, 0.5, 0.975)), q = c("lo", "mid", "hi")), keyby = population],
            totals[when == "after", .(indicator = "Admissions", value = quantile(admissions, c(0.025, 0.5, 0.975)), q = c("lo", "mid", "hi")), keyby = population],
            
            prev_peaks[, .(indicator = "Peak hospital beds", value = peak_all, q = "ref"), keyby = population],
            prev_peaks[, .(indicator = "Peak ICU beds", value = peak_icu, q = "ref"), keyby = population],
            peaks[, .(indicator = "Peak hospital beds", value = quantile(all, c(0.025, 0.5, 0.975)), q = c("lo", "mid", "hi")), keyby = population],
            peaks[, .(indicator = "Peak ICU beds", value = quantile(icu, c(0.025, 0.5, 0.975)), q = c("lo", "mid", "hi")), keyby = population],
            
            hi[, .(indicator = "Weeks of high ICU occupancy", value = quantile(icu_hi, c(0.025, 0.5, 0.975)), q = c("lo", "mid", "hi")), keyby = population],
            hi[, .(indicator = "Weeks of high hospital occupancy", value = quantile(all_hi, c(0.025, 0.5, 0.975)), q = c("lo", "mid", "hi")), keyby = population],
            
            tier[, .(indicator = "Weeks in Tier 2", value = quantile(tier2, c(0.025, 0.5, 0.975)), q = c("lo", "mid", "hi")), keyby = population],
            tier[, .(indicator = "Weeks in Tier 3", value = quantile(tier3, c(0.025, 0.5, 0.975)), q = c("lo", "mid", "hi")), keyby = population],
            tier[, .(indicator = "Weeks in lockdown",  value = quantile(cb, c(0.025, 0.5, 0.975)), q = c("lo", "mid", "hi")), keyby = population]
        )
        
        return (cbind(dcast(tab, population + indicator ~ q), scenario = scenario_output))
    }
}




cpp_vec = function(x) 
{
    paste("{", paste(x, collapse = ", "), "}")
}

cpp_obsI_tiers = function(tier_start_date, tier2_threshold, tier3_threshold, pop_size)
{
    tier_start_time = as.numeric(ymd(tier_start_date) - ymd("2020-01-01"));
    
    glue::glue(
        'if (t == 0) { dyn.scratch["tier"] = 0; dyn.scratch["revt"] = 0; dyn.scratch["cb"] = 0; }', # current tier and review time, circuit breaker yes/no
        'if (dyn.scratch["cb"] == 1) {',
        '        P.changes.ch[2].mode = Change::Bypass;',
        '        P.changes.ch[3].mode = Change::Bypass;',
        '}',
        'else if (t >= ${tier_start_time}) {',
        '    double inf_per100k_perweek = dyn("infections_p", t, {}, {}) / ${pop_size / 100000};',
        '    if (inf_per100k_perweek >= ${tier3_threshold} && dyn.scratch["tier"] < 3) {',
        '        P.changes.ch[2].mode = Change::Bypass;',
        '        P.changes.ch[3].mode = Change::Assign;',
        '        dyn.scratch["tier"] = 3;',
        '        dyn.scratch["revt"] = t + 28;',
        '    }',
        '    else if (inf_per100k_perweek >= ${tier2_threshold} && (dyn.scratch["tier"] < 2',
        '             || (t >= dyn.scratch["revt"] && dyn.scratch["tier"] >= 2))) {',
        '        P.changes.ch[2].mode = Change::Assign;',
        '        P.changes.ch[3].mode = Change::Bypass;',
        '        dyn.scratch["tier"] = 2;',
        '        dyn.scratch["revt"] = t + 28;',
        '    }',
        '    else if (inf_per100k_perweek < ${tier2_threshold} && (dyn.scratch["tier"] < 1 || t >= dyn.scratch["revt"])) {',
        '        P.changes.ch[2].mode = Change::Bypass;',
        '        P.changes.ch[3].mode = Change::Bypass;',
        '        dyn.scratch["tier"] = 1;',
        '        dyn.scratch["revt"] = t + 28;',
        '     }',
        '}',
        'dyn.Obs(t, 0, 1, 0) = dyn.scratch["tier"];',
        .sep = "\n", .open = "${", .close = "}")
}

cpp_obsI_cb = function(cb_dates, cb_durations, behaviour)
{
    cb_durations = rep_len(cb_durations, length(cb_dates))
    ret = 'if (t == 0) { dyn.scratch["cb"] = 0; }'; # cb yes/no;
    
    for (k in seq_along(cb_dates)) {
        cb_date = cb_dates[k];
        cb_duration = cb_durations[k];
        
        cb_t0 = as.integer(ymd(cb_date) - ymd("2020-01-01"));
        cb_t1 = cb_t0 + cb_duration;
        
        ret = c(ret,
            glue::glue(
                'if (t == ${cb_t0}) {',
                switch(behaviour,
                    default      = 'dyn.scratch["cb"] = 1;',
                    stay_in_tier = 'dyn.scratch["revt"] = t + 42; dyn.scratch["cb"] = 1;',
                    go_to_tier_3 = 'dyn.scratch["revt"] = t + 42; dyn.scratch["cb"] = 1;'
                ),
                '}',
                'if (t == ${cb_t1}) {',
                '    dyn.scratch["cb"] = 0;',
                if (behaviour == "go_to_tier_3") 'P.changes.ch[2].mode = Change::Bypass; P.changes.ch[3].mode = Change::Assign; dyn.scratch["tier"] = 3;' else '',
                '}',
                .sep = "\n", .open = "${", .close = "}")
        )
    }
    
    ret = c(ret,
        'dyn.Obs(t, 0, 2, 0) = dyn.scratch["cb"];'
    )
    
    return (ret)
}

cpp_chgI_close_schools = function(cb_dates, cb_durations)
{
    cb_durations = rep_len(cb_durations, length(cb_dates))
    ret = c(
        'P.changes.ch[5].values = { vector<double>(8, x[15]) };',
        'P.changes.ch[5].times = { x[17] };'
    )
    
    for (k in seq_along(cb_dates)) {
        cb_date = cb_dates[k];
        cb_duration = cb_durations[k];
        
        cb_t0 = as.integer(ymd(cb_date) - ymd("2020-01-01"));
        cb_t1 = cb_t0 + cb_duration;
        
        ret = c(ret,
            '{',
            glue::glue('double school_close = {cb_t0}, school_open = {cb_t1};'),
            
            'P.changes.ch[5].values.push_back(vector<double>(8, 1.0));',
            'P.changes.ch[5].values.push_back(vector<double>(8, x[15]));',
            'P.changes.ch[5].times.push_back(school_close);',
            'P.changes.ch[5].times.push_back(school_open + 1);',
    
            # fitting of google mobility indices
            'for (unsigned int k : vector<unsigned int> { 0, 2, 3 }) {',
            '    for (unsigned int i = 0; i < P.changes.ch[k].times.size(); ++i) {',
            '        double t = P.changes.ch[k].times[i];',
            '        if (t >= school_close && t <= school_open) {',
            '            P.changes.ch[k].values[i][2] = 0;',
            '            P.changes.ch[k].values[i][6] = 0;',
            '        }',
            '    }',
            '}',
            '}'
        )
    }
    
    return (ret)
}

project = function(popset, tiers = FALSE, tier2 = NULL, tier3 = NULL, 
    cb_date = NA, cb_duration = 14, lockdown = NULL, se = rep(0, 12), close_schools = FALSE, cb_behaviour = "default", 
    waning_duration = -1, seasonality = 0, n_run = 100, expire = NULL, parameters_only = FALSE)
{
    dynamicsO = list()
    paramsO = list()
    for (p in popset) {
        # calculate population size
        pop_size = sum(parametersI[[p]]$pop[[1]]$size);
        
        # Set parameters for this run . . .
        paramsI = rlang::duplicate(parametersI[[p]])
        paramsI$time1 = "2021-03-31";
        
        # get lockdown days in schedule1
        cb_days = rep(0, length(paramsI$schedule[[1]]$times))
        if (!is.na(cb_date[1])) {
            for (k in seq_along(cb_date)) {
                for (i in seq_along(paramsI$schedule[[1]]$values)) {
                    today = ymd(paramsI$date0) + paramsI$schedule[[1]]$times[i];
                    if (today %between% c(ymd(cb_date[k]), ymd(cb_date[k]) + rep_len(cb_duration, length(cb_date))[k] - 1)) {
                        cb_days[i] = 1;
                    }
                }
            }
        }
        
        # load user defined functions
        cm_source_backend(
            user_defined = list(
                model_v2 = list(
                    cpp_changes = c(
                        cpp_chgI_proj(cb_days, se),
                        if (close_schools) cpp_chgI_close_schools(cb_date, cb_duration) else NULL
                    ),
                    cpp_loglikelihood = "",
                    cpp_observer = c(
                        cpp_obsI(P.death), 
                        if (!is.na(cb_date[1])) cpp_obsI_cb(cb_date, cb_duration, cb_behaviour) else NULL,
                        if (tiers) cpp_obsI_tiers("2020-10-14", 700, 2100, pop_size) else NULL,
                        if (waning_duration > 0) paste("if (t == 274) { P.pop[0].wn = vector<double>(16,", 1 / waning_duration, "); }") else NULL,
                        if (seasonality > 0) paste("if (t == 274) { P.pop[0].season_A[0] = ", seasonality, "; }") else NULL
                    )
                )
            )
        )
        
        # Expire mobility changes after date "expire"
        if (!is.null(expire)) {
            expire_t = as.numeric(ymd(expire) - ymd(paramsI$date0))
            index = which(paramsI$schedule[[1]]$times == expire_t)
            for (i in index:length(paramsI$schedule[[1]]$values)) {
                paramsI$schedule[[1]]$values[[i]] = paramsI$schedule[[1]]$values[[index]]
            }
        }
        
        # adjust schedule1 for lockdown 
        if (!is.na(cb_date[1])) {
            for (k in seq_along(cb_date)) {
                for (i in seq_along(paramsI$schedule[[1]]$values)) {
                    today = ymd(paramsI$date0) + paramsI$schedule[[1]]$times[i];
                    if (today %between% c(ymd(cb_date[k]), ymd(cb_date[k]) + rep_len(cb_duration, length(cb_date))[k] - 1)) {
                        paramsI$schedule[[1]]$values[[i]][1] = max(0.0, paramsI$schedule[[1]]$values[[i]][1] + lockdown$res / 100); # res
                        paramsI$schedule[[1]]$values[[i]][2] = max(0.0, paramsI$schedule[[1]]$values[[i]][2] + lockdown$wor / 100); # work
                        paramsI$schedule[[1]]$values[[i]][3] = max(0.0, paramsI$schedule[[1]]$values[[i]][3] + lockdown$gro / 100); # groc
                        paramsI$schedule[[1]]$values[[i]][4] = max(0.0, paramsI$schedule[[1]]$values[[i]][4] + lockdown$ret / 100); # retail
                        paramsI$schedule[[1]]$values[[i]][5] = max(0.0, paramsI$schedule[[1]]$values[[i]][5] + lockdown$tra / 100); # transit
                    }
                }
            }
        }
        
        # contacts for tier 2
        for (i in seq_along(paramsI$schedule[[3]]$values)) {
            paramsI$schedule[[3]]$values[[i]][1] = max(0.0, paramsI$schedule[[1]]$values[[i]][1] + tier2$res / 100); # res
            paramsI$schedule[[3]]$values[[i]][2] = max(0.0, paramsI$schedule[[1]]$values[[i]][2] + tier2$wor / 100); # work
            paramsI$schedule[[3]]$values[[i]][3] = max(0.0, paramsI$schedule[[1]]$values[[i]][3] + tier2$gro / 100); # groc
            paramsI$schedule[[3]]$values[[i]][4] = max(0.0, paramsI$schedule[[1]]$values[[i]][4] + tier2$ret / 100); # retail
            paramsI$schedule[[3]]$values[[i]][5] = max(0.0, paramsI$schedule[[1]]$values[[i]][5] + tier2$tra / 100); # transit
        }
    
        # contacts for tier 3
        for (i in seq_along(paramsI$schedule[[4]]$values)) {
            paramsI$schedule[[4]]$values[[i]][1] = max(0.0, paramsI$schedule[[1]]$values[[i]][1] + tier3$res / 100); # res
            paramsI$schedule[[4]]$values[[i]][2] = max(0.0, paramsI$schedule[[1]]$values[[i]][2] + tier3$wor / 100); # work
            paramsI$schedule[[4]]$values[[i]][3] = max(0.0, paramsI$schedule[[1]]$values[[i]][3] + tier3$gro / 100); # groc
            paramsI$schedule[[4]]$values[[i]][4] = max(0.0, paramsI$schedule[[1]]$values[[i]][4] + tier3$ret / 100); # retail
            paramsI$schedule[[4]]$values[[i]][5] = max(0.0, paramsI$schedule[[1]]$values[[i]][5] + tier3$tra / 100); # transit
        }
        
        paramsI$processes[[length(paramsI$processes) + 1]] = list(
            source = "newII2",
            type = "multinomial",
            names = "infections",
            report = "ip",
            prob = matrix(1, nrow = 1, ncol = 16),
            delays = matrix(c(rep(0, 7 * 4 - 1), 1), nrow = 1, byrow = T)
        )
        
        # Sampling fits
        if (!parameters_only) {
            postI = rlang::duplicate(posteriorsI[[p]])
            postI[, stdnorm_ld := rnorm(.N)]
            postI[, stdnorm_t2 := rnorm(.N)]
            postI[, stdnorm_t3 := rnorm(.N)]
            test = cm_backend_sample_fit_test(cm_translate_parameters(paramsI), postI, n_run, seed = 0);
    
            test = rbindlist(test)
            test[, population := p]
            dynamicsO[[p]] = test
        }
        
        paramsO[[p]] = paramsI
        print(p)
    }
    
    if (parameters_only) {
        return (paramsO)
    }

    testX = rbindlist(dynamicsO)
    testX[, population := nhs_regions[population]]

    return (testX)
}

make_data = function(ld, sitreps)
{
    rbind(
        ld[, .(ValueType = "type28_death_inc_line", Geography = name,
            dmin = NA, d = as.Date(date), dmax = NA, ymin = NA, y = N, ymax = NA)],
        sitreps[, .(ValueType = "icu_prev", Geography = name,
            dmin = NA, d = as.Date(date), dmax = NA, ymin = NA, y = n_in_itu, ymax = NA)],
        sitreps[, .(ValueType = "hospital_prev", Geography = name,
            dmin = NA, d = as.Date(date), dmax = NA, ymin = NA, y = n_in_all_beds, ymax = NA)],
        sitreps[, .(ValueType = "hospital_inc", Geography = name,
            dmin = NA, d = as.Date(date), dmax = NA, ymin = NA, y = n_admitted_diagnosed, ymax = NA)]
    )
}

niceq = function(x)
{
    q = quantile(x, c(0.05, 0.5, 0.95));
    f = function(y) prettyNum(signif(y, 3), big.mark = ",")
    paste0(f(q[2]), " (", f(q[1]), " - ", f(q[3]), ")")
}

nicepc = function(x)
{
    q = quantile(x, c(0.05, 0.5, 0.95));
    f = function(y) round(y * 100, 0)
    paste0(f(q[2]), "% (", f(q[1]), " - ", f(q[3]), "%)")
}

combine_R = function(R_ldN4o, R_ldN4c) {
    R_ldN4 = cbind(R_ldN4o[, 1:3], R_ldN4c[, 3], R_ldN4o[, 4], R_ldN4c[, 4])
    names(R_ldN4)[c(3,5)] = paste0(names(R_ldN4)[c(3,5)], ", schools open")
    names(R_ldN4)[c(4,6)] = paste0(names(R_ldN4)[c(4,6)], ", schools closed")
    R_ldN4
}

arrs_ld = function(proj) {
    stop("Update.")

    # group == 1 is Rt, group == 2 is tier, group == 3 is lockdown
    trace = proj[, .(ld_0 = obs0[3], Rt_0 = obs0[1]), by = .(run, population, t)]
    trace[, ld_1 := shift(ld_0, -1), by = .(run, population)]
    trace[, ld_2 := shift(ld_0, -2), by = .(run, population)]
    trace[, ld_3 := shift(ld_0, -3), by = .(run, population)]
    trace[, Rt_1 := shift(Rt_0, -1), by = .(run, population)]
    trace[, Rt_2 := shift(Rt_0, -2), by = .(run, population)]
    trace[, Rt_3 := shift(Rt_0, -3), by = .(run, population)]
    
    lock = rbind(
        trace[ld_0 == 0 & ld_1 == 0 & ld_2 == 1 & ld_3 == 1, .(Rt_pre = Rt_0, Rt_in = Rt_3, effect = Rt_3 / Rt_0), by = .(run, population)]
    )

    rbind(
        lock[, .(population = "England", `Rt pre-lockdown` = niceq(Rt_pre), `Rt in lockdown` = niceq(Rt_in), `Reduction` = nicepc(1 - effect))],
        lock[, .(`Rt pre-lockdown` = niceq(Rt_pre), `Rt in lockdown` = niceq(Rt_in), `Reduction` = nicepc(1 - effect)), keyby = population]
    )
}

arrs_tier = function(proj) {
    stop("Update.")

    # group == 1 is Rt, group == 2 is tier, group == 3 is lockdown
    trace = proj[, .(tier_0 = obs0[2], Rt_0 = obs0[1]), by = .(run, population, t)]
    trace = trace[tier_0 > 0]
    trace[, tier_1 := shift(tier_0, -1), by = .(run, population)]
    trace[, tier_2 := shift(tier_0, -2), by = .(run, population)]
    trace[, tier_3 := shift(tier_0, -3), by = .(run, population)]
    trace[, Rt_1 := shift(Rt_0, -1), by = .(run, population)]
    trace[, Rt_2 := shift(Rt_0, -2), by = .(run, population)]
    trace[, Rt_3 := shift(Rt_0, -3), by = .(run, population)]
    
    tier_1_to_2 = rbind(
        trace[tier_0 == 1 & tier_1 == 1 & tier_2 == 2 & tier_3 == 2, .(effect = Rt_3 / Rt_0), by = .(run, population)],
        trace[tier_0 == 2 & tier_1 == 2 & tier_2 == 1 & tier_3 == 1, .(effect = Rt_0 / Rt_3), by = .(run, population)]
    )

    tier_2_to_3 = rbind(
        trace[tier_0 == 2 & tier_1 == 2 & tier_2 == 3 & tier_3 == 3, .(effect = Rt_3 / Rt_0), by = .(run, population)],
        trace[tier_0 == 3 & tier_1 == 3 & tier_2 == 2 & tier_3 == 2, .(effect = Rt_0 / Rt_3), by = .(run, population)]
    )
    
    tier_1_to_3 = cbind(
        tier_1_to_2[, .(t1_2 = sample(effect, 1000, replace = TRUE)), by = population],
        tier_2_to_3[, .(t2_3 = sample(effect, 1000, replace = TRUE)), by = population]
    )
    tier_1_to_3[, effect := t1_2 * t2_3]
    
    table_t2 = rbind(
        tier_1_to_2[, .(population = "England", `Tier 1 to Tier 2` = nicepc(1 - effect))],
        tier_1_to_2[, .(`Tier 1 to Tier 2` = nicepc(1 - effect)), keyby = population]
    )

    table_t3 = rbind(
        tier_1_to_3[, .(population = "England", `Tier 1 to Tier 3` = nicepc(1 - effect))],
        tier_1_to_3[, .(`Tier 1 to Tier 3` = nicepc(1 - effect)), keyby = population]
    )
    
    merge(table_t2, table_t3, by = "population", all = TRUE, sort = FALSE)
}

plot_indicator = function(series, indicator_name)
{
    ggplot(series[indicator == indicator_name]) +
        geom_pointrange(aes(x = population, ymin = lo, y = mid, ymax = hi, colour = scenario), position = position_dodge(width = 0.8), shape = 20, size = 0.0) +
        geom_linerange(aes(x = population, ymin = lo, y = mid, ymax = hi, colour = scenario), position = position_dodge(width = 0.8), size = 0.2) +
        labs(x = NULL, y = indicator_name) + ylim(0, NA)
}

plot_indicators_england = function(series, indicator_names, y_axis_title, unit = 1, legpos = "none", pal = "Dark2")
{
    ggplot(series[population == "England" & indicator %in% indicator_names]) +
        geom_ribbon(aes(x = scenario, ymin = lo / unit, ymax = hi / unit, fill = indicator, group = indicator), alpha = 0.5) +
        geom_line(aes(x = scenario, y = mid / unit, colour = indicator, group = indicator), size = 0.2) +
        labs(x = NULL, y = y_axis_title, fill = NULL, colour = NULL) + ylim(0, NA) +
        theme(panel.background = element_rect(fill = "#f4f4f4"),
            panel.grid.major = element_line(colour = "#ffffff", size = 0.2),
            legend.position = legpos) +
        scale_color_brewer(aesthetics = c("colour", "fill"), palette = pal)
}

plot_indicators_england2 = function(series, indicator_names, y_axis_title, unit = 1, legpos = "none", pal = "Accent")
{
    ggplot(series[population == "England" & indicator %in% indicator_names]) +
        geom_linerange(aes(x = scenario, ymin = lo / unit, ymax = hi / unit, colour = indicator, group = indicator),
            position = position_dodge(width = 0.5), size = 0.4) +
        geom_point(aes(x = scenario, y = mid / unit, colour = indicator, group = indicator),
            position = position_dodge(width = 0.5), size = 0.2, shape = 3) +
        labs(x = NULL, y = y_axis_title, colour = NULL) + ylim(0, NA) +
        theme(panel.background = element_rect(fill = "#f4f4f4"),
            panel.grid.major = element_line(colour = "#ffffff", size = 0.4),
            legend.position = legpos) +
        scale_color_brewer(aesthetics = c("colour", "fill"), palette = pal)
}

plot_indicators_england3 = function(series, indicator_names, y_axis_title, unit = 1, legpos = "none", pal = "Accent")
{
    ggplot(series[population == "England" & indicator %in% indicator_names]) +
        geom_col(aes(x = scenario, y = mid / unit, fill = indicator, group = indicator), 
            position = position_dodge(width = 0.6), size = 0.25, colour = "black", width = 0.4) +
        geom_linerange(aes(x = scenario, ymin = lo / unit, ymax = hi / unit, group = indicator),
            position = position_dodge(width = 0.6), size = 0.25) +
        labs(x = NULL, y = y_axis_title, fill = NULL) + ylim(0, NA) +
        theme(panel.background = element_rect(fill = "#f4f4f4"),
            panel.grid.major = element_line(colour = "#ffffff", size = 0.4),
            legend.position = legpos) +
        scale_color_brewer(aesthetics = c("colour", "fill"), palette = pal)
}

plot_indicators_england_stack2 = function(series, indicator_names, y_axis_title, stack_order, legpos = "none", pal = "Set2")
{
    ser2 = copy(series[population == "England" & indicator %in% indicator_names])
    ser2[, indicator := factor(indicator, stack_order)]
    ser3 = copy(ser2)
    ser3 = ser3[order(scenario, rev(indicator))]
    ser3[, lo := lo - mid]
    ser3[, hi := hi - mid]
    ser3[, mid := cumsum(mid), by = scenario]
    ser3[, lo := mid + lo]
    ser3[, hi := mid + hi]
    
    ggplot() +
        geom_col(data = ser2, aes(x = scenario, y = mid, fill = indicator), position = "stack", colour = "black", size = 0.25, width = 0.5) +
        geom_linerange(data = ser3, aes(x = scenario, ymin = lo, ymax = hi), size = 0.25) +
        labs(x = NULL, y = y_axis_title, fill = NULL, colour = NULL) + ylim(0, NA) +
        theme(panel.background = element_rect(fill = "#f4f4f4"),
            panel.grid.major = element_line(colour = "#ffffff", size = 0.2),
            legend.position = legpos) +
        scale_color_brewer(aesthetics = c("colour", "fill"), palette = pal)
}

plot_indicators_england_stack = function(series, indicator_names, y_axis_title, stack_order, legpos = "none", pal = "Set2")
{
    ser2 = rlang::duplicate(series[population == "England" & indicator %in% indicator_names])
    ser2[, indicator := factor(indicator, stack_order)]
    ser3 = rlang::duplicate(ser2)
    ser3 = ser3[order(scenario, rev(indicator))]
    ser3[, lo := lo - mid]
    ser3[, hi := hi - mid]
    ser3[, mid := cumsum(mid), by = scenario]
    ser3[, lo := mid + lo]
    ser3[, hi := mid + hi]

    ggplot() +
        geom_area(data = ser2, aes(x = scenario, y = mid, fill = indicator, group = indicator)) +
        geom_ribbon(data = ser3, aes(x = scenario, ymin = lo, ymax = hi, group = indicator), alpha = 0.25) +
        labs(x = NULL, y = y_axis_title, fill = NULL, colour = NULL) + ylim(0, NA) +
        theme(panel.background = element_rect(fill = "#f4f4f4"),
            panel.grid.major = element_line(colour = "#ffffff", size = 0.2),
            legend.position = legpos) +
        scale_color_brewer(aesthetics = c("colour", "fill"), palette = pal)
}

plot_indicator_regions = function(series, indicator_name, label = indicator_name, suffix = "\n(Thousands)", unit = 1000, legpos = "none", bigtitle = NULL)
{
    ser2 = rlang::duplicate(series[population != "England" & indicator == indicator_name])
    ser2[population == "East of England", population := "EE"]
    ser2[population == "London", population := "Ldn"]
    ser2[population == "Midlands", population := "Mlds"]
    ser2[population == "North East and Yorkshire", population := "NE&Y"]
    ser2[population == "North West", population := "NW"]
    ser2[population == "South East", population := "SE"]
    ser2[population == "South West", population := "SW"]
    ggplot(ser2) +
        geom_point(aes(x = population, y = mid / unit, colour = scenario), position = position_dodge(width = 0.8), shape = 20, size = 0.5) +
        geom_linerange(aes(x = population, ymin = lo / unit, ymax = hi / unit, colour = scenario), position = position_dodge(width = 0.8), size = 0.3) +
        labs(x = NULL, y = paste0(label, suffix), colour = NULL, title = bigtitle) + 
        ylim(0, NA) +
        theme(legend.position = legpos)
}


england_pops = c(1, 3, 4, 5, 6, 9, 10)
