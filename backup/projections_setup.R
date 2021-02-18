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
which_pops = c(1, 3, 4, 5, 6, 9, 10)

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

# Extract England column from each of tables_list with scenario names
england_only = function(tables_list, names)
{
    t = tables_list[[1]][, 1]
    for (i in seq_along(tables_list)) {
        t = cbind(t, tables_list[[i]][, England]);
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

plot_projection_2parts = function(proj_list, proj_names, from_date, colours, regions_a)
{
    p = NULL
    for (i in seq_along(proj_list)) {
        pp = arrange_projection(proj_list[[i]]);
        pp[, name := proj_names[[i]]];
        p = rbind(p, pp);
    }
    p[, name := factor(name, unique(name))]

    plotA = ggplot(p[!variable %in% c("tier", "cb") & ymd("2020-01-01") + t >= from_date & population %in% regions_a]);
    
    cbs = p[variable == "cb" & population %in% regions_a, .(t = ymd("2020-01-01") + t[`50%` > 0.5]), by = .(population)][, 
        .(tmin = min(t), tmax = max(t)), by = population]
    if (nrow(cbs) > 0) {
        plotA = plotA +
            geom_rect(data = cbs, aes(xmin = tmin, xmax = tmax, ymin = -Inf, ymax = Inf), fill = "black", alpha = 0.1)
    }
    
    rline = data.table(variable = factor("Rt", levels(p$variable)), y = 1)
    
    plotA = plotA +
        geom_ribbon(aes(x = ymd("2020-01-01") + t, ymin = `2.5%`, ymax = `97.5%`, fill = name), alpha = 0.4) +
        geom_hline(data = rline, aes(yintercept = y), size = 0.3) +
        geom_line(aes(x = ymd("2020-01-01") + t, y = `50%`, colour = name)) +
        facet_grid(variable ~ population, switch = "y", scales = "free") +
        cowplot::theme_cowplot(font_size = 11) + 
        theme(strip.background = element_blank(), strip.placement = "outside", 
            legend.position = "bottom",
            panel.background = element_rect(fill = "#f4f4f4"),
            panel.grid.major = element_line(colour = "#ffffff", size = 0.4),
            axis.text.x = element_text(angle = 90, vjust = 0.5)) +
        labs(x = NULL, y = NULL, colour = NULL, fill = NULL) +
        scale_x_date(date_breaks = "1 month", date_labels = "%b") +
        scale_color_manual(aesthetics = c("colour", "fill"), values = colours) +
        ylim(0, NA)
    
    # B

    plotB = ggplot(p[!variable %in% c("tier", "cb") & ymd("2020-01-01") + t >= from_date & !population %in% regions_a]);
    
    cbs = p[variable == "cb" & !population %in% regions_a, .(t = ymd("2020-01-01") + t[`50%` > 0.5]), by = .(population)][, 
        .(tmin = min(t), tmax = max(t)), by = population]
    if (nrow(cbs) > 0) {
        plotB = plotB +
            geom_rect(data = cbs, aes(xmin = tmin, xmax = tmax, ymin = -Inf, ymax = Inf), fill = "black", alpha = 0.1)
    }
    
    rline = data.table(variable = factor("Rt", levels(p$variable)), y = 1)
    
    plotB = plotB +
        geom_ribbon(aes(x = ymd("2020-01-01") + t, ymin = `2.5%`, ymax = `97.5%`, fill = name), alpha = 0.4) +
        geom_hline(data = rline, aes(yintercept = y), size = 0.3) +
        geom_line(aes(x = ymd("2020-01-01") + t, y = `50%`, colour = name)) +
        facet_grid(variable ~ population, switch = "y", scales = "free") +
        cowplot::theme_cowplot(font_size = 11) + 
        theme(strip.background = element_blank(), strip.placement = "outside", 
            legend.position = "bottom",
            panel.background = element_rect(fill = "#f4f4f4"),
            panel.grid.major = element_line(colour = "#ffffff", size = 0.4),
            axis.text.x = element_text(angle = 90, vjust = 0.5)) +
        labs(x = NULL, y = NULL, colour = NULL, fill = NULL) +
        scale_x_date(date_breaks = "1 month", date_labels = "%b") +
        scale_color_manual(aesthetics = c("colour", "fill"), values = colours) +
        ylim(0, NA)
    
    cowplot::plot_grid(plotA, plotB, nrow = 1, rel_widths = c(0.48, 0.52), labels = LETTERS, label_size = 10)
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


summarize_projection = function(proj, from_date, popsize, to_date = NULL, wh = "after")
{
    proj = copy(proj)
    if (!is.null(to_date)) {
        proj = proj[ymd("2020-01-01") + t <= to_date]
    }
    
    england = proj[, lapply(.SD, sum), .SDcols = 5:ncol(proj), by = .(run, t, group)]
    england[, population := "England"]
    proj = rbind(proj, england, fill = TRUE)
    
    geo = proj[, unique(population)]
    
    # Build summaries
    totals = proj[, .(
        deaths = sum(death_o + death2_o), admissions = sum(hosp_undetected_o + hosp_undetected2_o)), 
        by = .(population, run, when = ifelse(ymd("2020-01-01") + t < from_date, "before", "after"))]

    prev_peaks = proj[ymd("2020-01-01") + t < from_date,
        .(all = sum(hosp_p - hosp_undetected_p + hosp2_p - hosp_undetected2_p), icu = sum(icu_p + icu2_p), dea = sum(death_o + death2_o)),
        by = .(population, run, t)][,
        .(all = max(all), icu = max(icu), dea = max(dea)),
        by = .(population, run)][,
        .(peak_all = mean(all), peak_icu = mean(icu), peak_dea = mean(dea)), keyby = population]
    setkey(prev_peaks, population)
    
    if (wh == "before") {
        post_peaks = proj[ymd("2020-01-01") + t < from_date,
            .(all = sum(hosp_p - hosp_undetected_p + hosp2_p - hosp_undetected2_p), icu = sum(icu_p + icu2_p), dea = sum(death_o + death2_o)),
            by = .(population, run, t)][,
            .(all = max(all), icu = max(icu), dea = max(dea)),
            keyby = .(population, run)]
    } else {
        post_peaks = proj[ymd("2020-01-01") + t >= from_date,
            .(all = sum(hosp_p - hosp_undetected_p + hosp2_p - hosp_undetected2_p), icu = sum(icu_p + icu2_p), dea = sum(death_o + death2_o)),
            by = .(population, run, t)][,
            .(all = max(all), icu = max(icu), dea = max(dea)),
            keyby = .(population, run)]
    }
    setkey(post_peaks, population)
    
    peaks = post_peaks[,
        .(all = all, rel_all = all / prev_peaks[population, peak_all],
          dea = dea, rel_dea = dea / prev_peaks[population, peak_dea],
          icu = icu, rel_icu = icu / prev_peaks[population, peak_icu]),
        keyby = .(population, run)]
    
    
    tier = proj[ymd("2020-01-01") + t >= from_date,
        .(tier = obs0[2], cb = obs0[3]), keyby = .(population, run, t)][,
        .(tier2 = sum(tier == 2 & cb == 0) / 7, tier3 = sum(tier == 3 & cb == 0) / 7, tier4 = sum(cb == 1) / 7),
        keyby = .(population, run)]

    tier2 = merge(tier[population != "England"], popsize, by.x = "population", by.y = "Geography")
    tier2 = tier2[, .(population = "England", tier2 = weighted.mean(tier2, population_size),
        tier3 = weighted.mean(tier3, population_size), tier4 = weighted.mean(tier4, population_size)), by = run]
    tier = rbind(tier2, tier[population != "England"])

    # Factor
    totals[, population := factor(population, geo)]
    peaks[, population := factor(population, geo)]
    tier[, population := factor(population, geo)]
    
    tab = rbind(
        totals[when == wh, .(indicator = "Total deaths", value = niceq(deaths)), keyby = population],
        totals[when == wh, .(indicator = "Total admissions", value = niceq(admissions)), keyby = population],
        peaks[, .(indicator = "Peak deaths", value = niceq(dea)), keyby = population],
        peaks[, .(indicator = "Peak ICU requirement", value = niceq(icu)), keyby = population],
        peaks[, .(indicator = "Peak ICU (rel. to 1st wave)", value = nicepc(rel_icu)), keyby = population],
        tier[, .(indicator = "Weeks in Tier 2", value = niceq(tier2)), keyby = population],
        tier[, .(indicator = "Weeks in Tier 3", value = niceq(tier3)), keyby = population],
        tier[, .(indicator = "Weeks in Tier 4", value = niceq(tier4)), keyby = population]
    )
    return (dcast(tab, indicator ~ population))
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

cpp_obsI_partclose_schools = function(params, d_school_date0 = NULL, d_school_date1 = NULL, d_school_contact = NULL)
{
    if (is.null(d_school_date0)) {
        return (NULL);
    }
    
    old_matrix = c(t(params$pop[[1]]$matrices[[3]]));
    adj_matrix = c(t(params$pop[[1]]$matrices[[3]] * d_school_contact));
    
    t0 = as.integer(ymd(d_school_date0) - ymd("2020-01-01"));
    t1 = as.integer(ymd(d_school_date1) - ymd("2020-01-01"));
    
    glue::glue(
        'if (t == ${t0}) {',
        '    P.pop[0].matrices[2].x = ${ cpp_vec(adj_matrix) };',
        '    P.pop[0].needs_recalc = true;',
        '    P.pop[0].Recalculate();',
        '}',
        'if (t == ${t1}) {',
        '    P.pop[0].matrices[2].x = ${ cpp_vec(old_matrix) };',
        '    P.pop[0].needs_recalc = true;',
        '    P.pop[0].Recalculate();',
        '}',
        .sep = "\n", .open = "${", .close = "}")
}

cpp_chgI_close_schools = function(school_breaks, school_factors_r = NULL, school_factors_c = NULL)
{
    ret = c(
        'P.changes.ch[5].values = { vector<double>(8, x[15]) };',
        'P.changes.ch[5].times = { x[17] };'
    )
    
    if (is.null(school_factors_r)) {
        school_factors_r = rep(0, length(school_breaks) %/% 2);
    }
    if (is.null(school_factors_c)) {
        school_factors_c = rep(0, length(school_breaks) %/% 2);
    }
    
    for (k in seq(1, length(school_breaks), by = 2)) {
        cb_t0 = as.integer(ymd(school_breaks[k]) - ymd("2020-01-01"));
        cb_t1 = as.integer(ymd(school_breaks[k + 1]) - ymd("2020-01-01"));
        school_factor_r = school_factors_r[(k + 1) %/% 2];
        school_factor_c = school_factors_c[(k + 1) %/% 2];

        ret = c(ret,
            '{',
            glue::glue('double school_close = {cb_t0}, school_open = {cb_t1}, school_factor_r = {school_factor_r}, school_factor_c = {school_factor_c};'),
            
            'P.changes.ch[5].values.push_back(vector<double>(8, (1 - school_factor_r) * 1.0 + school_factor_r * x[15]));',
            'P.changes.ch[5].values.push_back(vector<double>(8, x[15]));',
            'P.changes.ch[5].times.push_back(school_close);',
            'P.changes.ch[5].times.push_back(school_open);',
    
            # fitting of google mobility indices
            'for (unsigned int k : vector<unsigned int> { 0, 2, 3 }) {',
            '    for (unsigned int i = 0; i < P.changes.ch[k].times.size(); ++i) {',
            '        double t = P.changes.ch[k].times[i];',
            '        if (t >= school_close && t < school_open) {',
            '            P.changes.ch[k].values[i][2] *= school_factor_c;',
            '            P.changes.ch[k].values[i][6] *= school_factor_c;',
            '        }',
            '    }',
            '}',
            '}'
        )
    }
    
    return (ret)
}

project = function(popset, tiers = FALSE, tier2 = NULL, tier3 = NULL, 
    cb_date = NA, cb_duration = 14, lockdown = NULL, se = rep(0, 12), school_breaks = NULL, cb_behaviour = "default", 
    waning_duration = -1, seasonality = 0, n_run = 100, expire = NULL, vacc = NULL, ei_v = rep(1, 16), ed_vi = rep(1, 16), parameters_only = FALSE, 
    school_factors_r = NULL, school_factors_c = NULL, d_school_date0 = NULL, d_school_date1 = NULL, d_school_contact = NULL)
{
    dynamicsO = list()
    paramsO = list()
    for (p in popset) {
        # calculate population size
        pop_size = sum(parametersI[[p]]$pop[[1]]$size);
        
        # Set parameters for this run . . .
        paramsI = rlang::duplicate(parametersI[[p]])
        paramsI$time1 = "2021-06-30";
        
        paramsI$pop[[1]]$ei_v = ei_v;
        paramsI$pop[[1]]$ei2_v = ei_v;
        paramsI$pop[[1]]$ed_vi = ed_vi;
        paramsI$pop[[1]]$ed_vi2 = ed_vi;
        
        # get lockdown days in schedule1
        # Old version with multiple circuit breakers:
        # cb_days = rep(0, length(paramsI$schedule[[1]]$times))
        # if (!is.na(cb_date[1])) {
        #     for (k in seq_along(cb_date)) {
        #         for (i in seq_along(paramsI$schedule[[1]]$values)) {
        #             today = ymd(paramsI$date0) + paramsI$schedule[[1]]$times[i];
        #             if (today %between% c(ymd(cb_date[k]), ymd(cb_date[k]) + rep_len(cb_duration, length(cb_date))[k] - 1)) {
        #                 cb_days[i] = 1;
        #             }
        #         }
        #     }
        # }

        # New version with diff circuit breakers in different regions:
        cb_days = rep(0, length(paramsI$schedule[[1]]$times))
        if (!is.na(cb_date[p])) {
            for (i in seq_along(paramsI$schedule[[1]]$values)) {
                today = ymd(paramsI$date0) + paramsI$schedule[[1]]$times[i];
                if (today %between% c(ymd(cb_date[p]), ymd(cb_date[p]) + cb_duration[p] - 1)) {
                    cb_days[i] = 1;
                }
            }
        }
        
        # load user defined functions
        cm_source_backend(
            user_defined = list(
                model_v2 = list(
                    cpp_changes = c(
                        cpp_chgI_proj(cb_days, se, ncol(posteriorsI[[p]]) - 4),
                        cpp_chgI_close_schools(school_breaks, school_factors_r, school_factors_c)
                    ),
                    cpp_loglikelihood = "",
                    cpp_observer = c(
                        cpp_obsI(P.death), 
                        if (!is.na(cb_date[p])) cpp_obsI_cb(cb_date[p], cb_duration[p], cb_behaviour) else NULL,
                        if (tiers) cpp_obsI_tiers("2020-12-15", 700, 1400, pop_size) else NULL,
                        if (waning_duration > 0) paste("if (t == 274) { P.pop[0].wn = vector<double>(16,", 1 / waning_duration, "); }") else NULL,
                        if (seasonality > 0) paste("if (t == 274) { P.pop[0].season_A[0] = ", seasonality, "; }") else NULL,
                        if (!is.null(vacc)) cpp_obsI_vax(parametersI[[p]], vacc[[p]]) else NULL,
                        if (!is.null(d_school_contact)) cpp_obsI_partclose_schools(paramsI, d_school_date0, d_school_date1, d_school_contact) else NULL
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
        if (!is.na(cb_date[p])) {
            for (i in seq_along(paramsI$schedule[[1]]$values)) {
                today = ymd(paramsI$date0) + paramsI$schedule[[1]]$times[i];
                if (today %between% c(ymd(cb_date[p]), ymd(cb_date[p]) + cb_duration[p] - 1)) {
                    paramsI$schedule[[1]]$values[[i]][1] = max(0.0, paramsI$schedule[[1]]$values[[i]][1] + lockdown$res / 100); # res
                    paramsI$schedule[[1]]$values[[i]][2] = max(0.0, paramsI$schedule[[1]]$values[[i]][2] + lockdown$wor / 100); # work
                    paramsI$schedule[[1]]$values[[i]][3] = max(0.0, paramsI$schedule[[1]]$values[[i]][3] + lockdown$gro / 100); # groc
                    paramsI$schedule[[1]]$values[[i]][4] = max(0.0, paramsI$schedule[[1]]$values[[i]][4] + lockdown$ret / 100); # retail
                    paramsI$schedule[[1]]$values[[i]][5] = max(0.0, paramsI$schedule[[1]]$values[[i]][5] + lockdown$tra / 100); # transit
                }
            }
        }

        # Old version with multiple CBs
        # if (!is.na(cb_date[1])) {
        #     for (k in seq_along(cb_date)) {
        #         for (i in seq_along(paramsI$schedule[[1]]$values)) {
        #             today = ymd(paramsI$date0) + paramsI$schedule[[1]]$times[i];
        #             if (today %between% c(ymd(cb_date[k]), ymd(cb_date[k]) + rep_len(cb_duration, length(cb_date))[k] - 1)) {
        #                 paramsI$schedule[[1]]$values[[i]][1] = max(0.0, paramsI$schedule[[1]]$values[[i]][1] + lockdown$res / 100); # res
        #                 paramsI$schedule[[1]]$values[[i]][2] = max(0.0, paramsI$schedule[[1]]$values[[i]][2] + lockdown$wor / 100); # work
        #                 paramsI$schedule[[1]]$values[[i]][3] = max(0.0, paramsI$schedule[[1]]$values[[i]][3] + lockdown$gro / 100); # groc
        #                 paramsI$schedule[[1]]$values[[i]][4] = max(0.0, paramsI$schedule[[1]]$values[[i]][4] + lockdown$ret / 100); # retail
        #                 paramsI$schedule[[1]]$values[[i]][5] = max(0.0, paramsI$schedule[[1]]$values[[i]][5] + lockdown$tra / 100); # transit
        #             }
        #         }
        #     }
        # }
        
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
            postI = copy(posteriorsI[[p]])
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

make_vaccine_schedule = function(ymd_start, weekly_v, targeting, popset)
{
    # Get total population
    popsize2 = NULL
    for (i in popset) {
        if (!is.null(parametersI[[i]])) {
            popsize2 = rbind(popsize2,
                data.table(p = i, 
                    age = parametersI[[i]]$pop[[1]]$group_names,
                    population_size = parametersI[[i]]$pop[[1]]$size)
            )
        }
    }
    total_popsize = popsize2[, .(population_size = sum(population_size)), by = age]$population_size;

    # Make vaccination schedule
    npops = length(popset)
    
    schedule = data.table(date = rep(ymd(ymd_start) + 0:(length(weekly_v) - 1), each = npops), 
        pop = rep(popset, length(weekly_v)),
        daily_v = rep(weekly_v / 7, each = npops),
        v = rep(list(rep(0, 16)), length(weekly_v) * length(npops)));
    setkey(schedule, date);
    
    # Allocate to each region by population and targeting
    target_prev = rep(0, 16);
    target_i = 1;

    vacc = rep(0, 16);

    for (d in schedule[, unique(date)])
    {
        alloc_pattern = targeting[[target_i]] - target_prev;
        alloc_target = total_popsize * alloc_pattern;
        alloc_today = schedule[date == d, daily_v[1]] * alloc_target / sum(alloc_target);
        
        s = cbind(popsize2[, 1:3], alloc_today = rep(alloc_today, npops));
        s[, alloc_today := alloc_today * population_size / sum(population_size), by = age];
        for (p in popset) {
            schedule[date == d & pop == p, v := s[p == ..p, alloc_today]]
        }
        
        vacc = vacc + alloc_today;
        if (all(vacc >= alloc_target)) {
            target_prev = targeting[[target_i]];
            target_i = target_i + 1;
            if (target_i > length(targeting)) {
                break;
            }
            if (!any(targeting[[target_i]] > target_prev)) {
                stop("Malformed targeting")
            }
        }
    }
    
    # Keep only rows in schedule that change vaccination allotment
    ret = list()

    for (i in popset) {
        # Keep only rows in schedule that change vaccination allotment
        sched = schedule[pop == i];
        for (r in nrow(sched):2)
        {
            if (all(sched[r, v][[1]]   == sched[r - 1, v][[1]])) {
                sched = sched[-r]
            }
        }
        sched = rbind(sched,
            data.table(date = ymd(ymd_start) + length(weekly_v), pop = i, 
                daily_v = 0, v = list(rep(0, 16))))
        
        ret[[i]] = list(
            vt = sched[, date],
            v = sched[, v]
        );
    }
    
    return (ret)
}


england_pops = c(1, 3, 4, 5, 6, 9, 10)
