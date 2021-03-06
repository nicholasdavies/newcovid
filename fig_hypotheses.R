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
library(extraDistr)
library(Hmisc)

theme_set(cowplot::theme_cowplot(font_size = 10) + theme(strip.background = element_blank()))

which_pops = c(1, 3, 9)
england_pops = c(1, 3, 4, 5, 6, 9, 10)

uk_covid_data_path = "./fitting_data/";
datapath = function(x) paste0(uk_covid_data_path, x)
pct = function(x) as.numeric(str_replace_all(x, "%", "")) / 100
data_file = "processed-data-2021-01-18.qs"
date_fitting = "2020-12-24"
date_predict = "2021-01-07"

#
# SETUP
#

# set up covidm
cm_path = "./covidm_for_fitting/";
cm_force_rebuild = F;
cm_build_verbose = F;
cm_version = 2;
source(paste0(cm_path, "/R/covidm.R"))
popUK = readRDS(datapath("popNHS.rds"));
matricesUK = readRDS(datapath("matricesNHS.rds"));

cm_populations = rbind(cm_populations[name != "United Kingdom"], popUK)
cm_matrices = c(cm_matrices, matricesUK)
source("./distribution_fit.R");
source("./spim_output.R");
source("./check_fit.R");
source("./project-helper.R");
source("./cpp_funcs.R")
source("./processes.R")

nhs_regions = popUK[, unique(name)]

all_data = qread(datapath(data_file))
ld = all_data[[1]]
sitreps = all_data[[2]]
virus = all_data[[3]][!Data.source %like% "7a|7b|6a|6b" & omega > 1e-9]
sero = all_data[[4]]
sgtf = all_data[[5]]


load_hyp = function(filename, check = FALSE, withvoc = TRUE)
{
    # Load fit
    fit = qread(filename);
    
    # Run fit
    posteriorsI = fit[[1]]
    parametersI = fit[[2]]
    post1 = posteriorsI[[which_pops[1]]]
    priorsI = names(post1)[5:(ncol(post1))];
    names(priorsI) = priorsI;
    opt_seas = FALSE;
    opt_conc = TRUE;
    opt_v2 = "v2_when" %in% priorsI;
    opt_relu = "v2_relu" %in% priorsI;
    opt_latdur = "v2_latdur" %in% priorsI;
    opt_serial = "v2_serial" %in% priorsI;
    opt_infdur = "v2_infdur" %in% priorsI;
    opt_immesc = "v2_immesc" %in% priorsI;
    opt_ch_u = "v2_ch_u" %in% priorsI;
    
    if (!withvoc) {
        opt_v2 = FALSE;
    }
    
    # Generate SPI-M output
    # Sample dynamics from fit
    dynamicsI = list()
    ll_mean_theta = list()
    prev_ll = list()
    new_ll = list()
    for (p in which_pops)  {
        cat(paste0("Sampling fit for population ", p, "...\n"))
        
        ldI = rlang::duplicate(ld);
        ldI = ldI[pid == p - 1];
        sitrepsI = rlang::duplicate(sitreps);
        sitrepsI = sitrepsI[pid == p - 1];
        seroI = rlang::duplicate(sero);
        seroI = seroI[pid == p - 1 & Data.source != "NHSBT"];   # sero: all but NHSBT
        virusI = rlang::duplicate(virus);
        virusI = virusI[pid == p - 1 & Data.source %like% "REACT"]; # virus: REACT only
        sgtfI = copy(sgtf);
        sgtfI = sgtfI[pid == p - 1];
        
        # Source backend
        cm_source_backend(
            user_defined = list(
                model_v2 = list(
                    cpp_changes = cpp_chgI_voc(priorsI, seasonality = opt_seas, 
                        v2 = opt_v2, v2_relu = opt_relu, v2_latdur = opt_latdur, v2_serial = opt_serial, v2_infdur = opt_infdur, v2_immesc = opt_immesc, v2_ch_u = opt_ch_u),
                    cpp_loglikelihood = cpp_likI_voc(parametersI[[p]], ldI, sitrepsI, seroI, virusI, sgtfI, p, date_predict, priorsI, death_cutoff = 7, use_sgtf = opt_v2),
                    cpp_observer = cpp_obsI_voc(concentration = opt_conc, v2 = opt_v2, P.death, P.critical, priorsI)
                )
            )
        )

        # Sampling fits
        paramsI2 = rlang::duplicate(parametersI[[p]])
        paramsI2$time1 = as.character(ymd(parametersI[[p]]$time1) + 56);
        test = cm_backend_sample_fit_test(cm_translate_parameters(paramsI2), posteriorsI[[p]], 100, seed = 0);
        rows = cm_backend_sample_fit_rows(cm_translate_parameters(paramsI2), posteriorsI[[p]], 100, seed = 0);
        
        test = rbindlist(test)
        test[, population := p]
        
        # Add dispersion parameters
        disp = posteriorsI[[p]][rows, .SD, .SDcols = patterns("^disp|v2_conc|v2_disp|v2_sgtf0")]
        disp[, run := .I]
        test = merge(test, disp, by = "run")

        dynamicsI[[p]] = test
        
        # Posterior mean sampling
        post_mean = posteriorsI[[p]][, lapply(.SD, mean)]
        theta = unname(unlist(c(post_mean[, 5:ncol(post_mean)])))
        ll_mean_theta[[p]] = cm_backend_loglik(cm_translate_parameters(paramsI2), theta, seed = 0);
        
        # Recover previous and predictive LL
        prev_ll[[p]] = posteriorsI[[p]][rows, ll]
        new_ll[[p]] = rep(0, 100);
        for (i in 1:100) {
            theta = unname(unlist(c(posteriorsI[[p]][rows[i], 5:ncol(posteriorsI[[p]])])))
            new_ll[[p]][i] = cm_backend_loglik(cm_translate_parameters(paramsI2), theta, seed = 0);
        }
    }
    
    # Concatenate dynamics for SPI-M output
    test = rbindlist(dynamicsI, fill = TRUE)
    test[, population := nhs_regions[population]]

    # Visually inspect fit
    if (check) {
        plot = check_fit(test, parametersI, ld, sitreps, virus, sero, nhs_regions[which_pops], death_cutoff = 0, "2020-12-31")
        plot = plot + geom_vline(aes(xintercept = ymd("2020-12-24")), size = 0.25, linetype = "42")
        return (plot)
    }

    # Posteriors
    post = rbindlist(posteriorsI, idcol = "population", fill = TRUE)
    post[, pop := nhs_regions[population]]
    post = post[population %in% which_pops]
    post_melted = melt(post, id.vars = c(1:5, ncol(post)))

    # Calculate DIC
    post[, D := -2 * ll]
    dic1 = post[, 0.5 * var(D) + mean(D), by = population][, mean(V1)];

    # Assess fit to SGTF data
    sgtf[, qlo := qbeta(0.025, sgtf + 1, other + 1)]
    sgtf[, qhi := qbeta(0.975, sgtf + 1, other + 1)]
    if ("v2_conc" %in% names(test)) {
        vmodel = test[, .(I1 = sum(test_o), I2 = sum(test2_o), sgtf0 = v2_sgtf0[1], conc = v2_conc[1]), by = .(t, population, run)]
    } else {
        vmodel = test[, .(I1 = sum(test_o), I2 = sum(test2_o), sgtf0 = v2_sgtf0[1], conc = 1 / (v2_disp[1] * v2_disp[1])), by = .(t, population, run)]
    }
    vmodel[, p2 := I2 / (I1 + I2)]
    vmodel[is.nan(p2), p2 := 0]
    vmodel[, sgtf := (1 - p2) * sgtf0 + p2];
    vmodel[, alpha := sgtf * (conc - 2) + 1];
    vmodel[, beta := (1 - sgtf) * (conc - 2) + 1];
    
    times = vmodel[population == nhs_regions[which_pops[1]] & run == 1, t];
    runs = vmodel[, max(run)];
    reps = 100;
    
    quants = NULL
    for (p in which_pops)
    {
        cat("Quantiles for", nhs_regions[p], "...\n");
        vmat = matrix(0, nrow = length(times), ncol = runs * reps);
        c = 1;
        for (REP in 1:reps) {
            for (RUN in 1:runs) {
                vmat[, c] = vmodel[population == nhs_regions[p] & run == RUN, qbeta(runif(1), alpha, beta)];
                c = c + 1;
            }
        }
        quants0 = as.data.table(t(apply(vmat, 1, quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))))
        quants0[, population := nhs_regions[p]];
        quants0[, t := times];
        quants = rbind(quants, quants0);
    }
    
    # Predictive value (strain frequencies)
    pv = merge(vmodel[, .(date = ymd("2020-01-01") + t, population, run, alpha, beta)], sgtf[, .(date, population = nhs_name, sgtf, other)], 
        by = c("population", "date"))
    pv = pv[, sum(dbbinom(sgtf, sgtf + other, alpha, beta, log = TRUE)), by = .(population, run)]
    
    # Assess growth rate / Rt
    gr = growth_rate(test)

    return (list(post = post_melted, dic1 = dic1, vmodel = quants, 
        ll_mean_theta = ll_mean_theta, pv = pv, prev_ll = prev_ll, new_ll = new_ll,
        dr = gr$dr, dRt = gr$dRt))
}

south_east = function(filename, keep_voc = TRUE, plot_until = "2020-12-24")
{
    # Load fit
    fit = qread(filename);
    
    # Run fit
    posteriorsI = fit[[1]]
    parametersI = fit[[2]]
    post1 = posteriorsI[[which_pops[1]]]
    priorsI = names(post1)[5:(ncol(post1))];
    names(priorsI) = priorsI;
    opt_seas = FALSE;
    opt_conc = TRUE;
    opt_v2 = "v2_when" %in% priorsI;
    opt_relu = "v2_relu" %in% priorsI;
    opt_latdur = "v2_latdur" %in% priorsI;
    opt_serial = "v2_serial" %in% priorsI;
    opt_infdur = "v2_infdur" %in% priorsI;
    opt_immesc = "v2_immesc" %in% priorsI;
    opt_ch_u = "v2_ch_u" %in% priorsI;
    
    # Generate SPI-M output
    # Sample dynamics from fit
    dynamicsI = list()

    new_ll = list()
    for (p in 9)  {
        cat(paste0("Sampling fit for population ", p, "...\n"))
        
        ldI = rlang::duplicate(ld);
        ldI = ldI[pid == p - 1];
        sitrepsI = rlang::duplicate(sitreps);
        sitrepsI = sitrepsI[pid == p - 1];
        seroI = rlang::duplicate(sero);
        seroI = seroI[pid == p - 1 & Data.source != "NHSBT"];   # sero: all but NHSBT
        virusI = rlang::duplicate(virus);
        virusI = virusI[pid == p - 1 & Data.source %like% "REACT"]; # virus: REACT only
        sgtfI = copy(sgtf);
        sgtfI = sgtfI[pid == p - 1];
        
        # Source backend
        cm_source_backend(
            user_defined = list(
                model_v2 = list(
                    cpp_changes = cpp_chgI_voc(priorsI, seasonality = opt_seas, 
                        v2 = opt_v2, v2_relu = opt_relu, v2_latdur = opt_latdur, v2_serial = opt_serial, v2_infdur = opt_infdur, v2_immesc = opt_immesc, v2_ch_u = opt_ch_u),
                    cpp_loglikelihood = cpp_likI_voc(parametersI[[p]], ldI, sitrepsI, seroI, virusI, sgtfI, p, date_predict, priorsI, death_cutoff = 7, use_sgtf = opt_v2),
                    cpp_observer = cpp_obsI_voc(concentration = opt_conc, v2 = opt_v2, P.death, P.critical, priorsI)
                )
            )
        )
        
        # Sampling fits
        paramsI2 = rlang::duplicate(parametersI[[p]])
        paramsI2$time1 = as.character(ymd(parametersI[[p]]$time1) + 56);
        
        if (keep_voc == FALSE) {
            posteriorsI[[p]][, v2_when := 99999];
        }
        
        test = cm_backend_sample_fit_test(cm_translate_parameters(paramsI2), posteriorsI[[p]], 100, seed = 0);
        rows = cm_backend_sample_fit_rows(cm_translate_parameters(paramsI2), posteriorsI[[p]], 100, seed = 0);
        
        test = rbindlist(test)
        test[, population := p]
        
        # Add dispersion parameters
        disp = posteriorsI[[p]][rows, .SD, .SDcols = patterns("^disp|v2_conc|v2_disp|v2_sgtf0")]
        disp[, run := .I]
        test = merge(test, disp, by = "run")

        dynamicsI[[p]] = test
    }
    
    # Concatenate dynamics for SPI-M output
    test = rbindlist(dynamicsI, fill = TRUE)
    test[, population := nhs_regions[population]]

    # Visually inspect fit
    p = check_fit(test, parametersI, ld[date <= plot_until], sitreps[date <= plot_until], virus, sero, nhs_regions[9], death_cutoff = 0, plot_until)

    return (list(dynamics = test, plot = p + geom_vline(aes(xintercept = ymd("2020-12-24")), size = 0.25, linetype = "42")))
}

predictive = function(hyp)
{
    diff = 0;
    for (p in which_pops) {
        diff = diff + hyp$new_ll[[p]] - hyp$prev_ll[[p]]
    }
    cat("Mean", mean(diff), "sd", sd(diff), "deviance", -2 * mean(diff), "\n")
    return (-2 * mean(diff))
}

# Plot posteriors
plot_posterior = function(filename, basename, pops)
{
    post = qread(filename)[[1]]
    post = rbindlist(post, idcol = "population", fill = TRUE)
    post[, pop := nhs_regions[population]]
    post = post[population %in% pops]
    post = melt(post, id.vars = c(1:5, ncol(post)))
    post = post[!is.na(value)]
    
    plot_posterior = ggplot(post) +
        geom_density(aes(value, colour = pop)) + 
        facet_wrap(~variable, scales = "free") +
        labs(x = "Parameter value", y = "Density", colour = "NHS England\nregion")
    
    ggsave(paste0("./output/posterior_", basename, ".pdf"), plot_posterior, width = 30, height = 20, units = "cm", useDingbats = FALSE)
    ggsave(paste0("./output/posterior_", basename, ".png"), plot_posterior, width = 30, height = 20, units = "cm")
}

# Estimate posterior distribution for parameter
overall = function(fit, varname, log_transform, pops)
{
    post = copy(fit$post)
    post = merge(post, popsize, by.x = "pop", by.y = "Geography")
    post = post[variable == varname]
    
    if (log_transform) {
        post[, value := log(value)];
    }
    
    calc = function(populations)
    {
        samples = post[population %in% populations, .(
                mu = wtd.mean(value, weights = population_size, normwt = TRUE),
                sd = sqrt(wtd.var(value, weights = population_size, normwt = TRUE))
            ), by = .(trial, chain)];
        samples = samples[!is.nan(sd)];
        samples[, y := rnorm(.N, mu, sd)];
        
        if (log_transform) {
            samples[, y := exp(y)];
        }
        o_mean = samples[, mean(y)];
        o_025 = samples[, quantile(y, 0.025)];
        o_975 = samples[, quantile(y, 0.975)];
        o_above = samples[, mean(y > ifelse(log_transform, 1, 0))];

        cat(o_mean, " (", o_025, " - ", o_975, ") [P(y increased) = ", o_above, "]\n");
        return (c(o_mean, o_025, o_975))
    }
    
    calc(pops)
}

hyp_plots = function(h, pops, vars, lt, varnames)
{
    ps = ggplot(h$vmodel[t > 274 & population %in% nhs_regions[pops]]) +
        geom_ribbon(aes(x = t + ymd("2020-01-01"), ymin = `2.5%`, ymax = `97.5%`), fill = "darkorchid", alpha = 0.4) +
        geom_line(aes(x = t + ymd("2020-01-01"), y = `50%`), colour = "darkorchid") +
        geom_ribbon(data = sgtf2, aes(x = date, ymin = lo, ymax = hi), alpha = 0.4) +
        geom_line(data = sgtf2, aes(x = date, y = mean)) +
        geom_vline(aes(xintercept = ymd("2020-12-24")), size = 0.25, linetype = "22") +
        # scale_y_continuous(trans = scales::logit_trans(), breaks = c(0.01, 0.1, 0.2, 0.3, 0.5, 0.7, 0.8, 0.9, 0.99), limits = c(0.01, 0.99)) +
        facet_wrap(~population)
    
    dv = NULL
    for (v in seq_along(vars)) {
        ci = overall(h, vars[v], lt[v], pops)
        dv = rbind(dv,
            data.table(var = varnames[v], mean = ci[1], lo = ci[2], hi = ci[3])
        )
    }
    
    # For calculating delta-r
    dr_fit = list(post = copy(h$dr))
    dr_fit$post[, trial := run]
    dr_fit$post[, chain := 0]
    dr_fit$post[, variable := "dr"]
    dr_fit$post[, value := dr * 7]
    dr_fit$post[, pop := population]
    dr_fit$post[, population := as.numeric(match(population, nhs_regions))]
    
    cidr = exp(overall(dr_fit, "dr", FALSE, c(1, 3, 9)))
    cih = exp(overall(h, "v2_hosp_rlo", FALSE, pops))
    cii = exp(overall(h, "v2_icu_rlo",  FALSE, pops))
    cid = exp(overall(h, "v2_cfr_rlo",  FALSE, pops))

    dv = rbind(dv,
        data.table(var = "Growth rate", mean = cidr[1], lo = cidr[2], hi = cidr[3]),
        data.table(var = "Severe", mean = cih[1], lo = cih[2], hi = cih[3]),
        data.table(var = "Critical",  mean = cii[1], lo = cii[2], hi = cii[3]),
        data.table(var = "Fatal",  mean = cid[1], lo = cid[2], hi = cid[3])
    )
    
    dv[, var := factor(var, rev(unique(var)))]
    
    pp = ggplot(dv) +
        geom_pointrange(aes(x = mean, xmin = lo, xmax = hi, y = var, colour = var), shape = 20, size = 0.5)
    
    return (list(ps, pp))
}

style = function(p, titl, n, yl1, region_titles, max_x = 3.25)
{
    palette = c("#6388b4", "#ffae34", "#ef6f6a", "#8cc2ca", "#55ad89", "#c3bc3f", "#bb7693", "#baa094", "#a9b5ae", "#767676")

    p1 = p[[1]] + labs(x = NULL, y = ifelse(yl1, "Frequency of S gene target failure", ""), title = titl) + ylim(0, 1) +
            theme(plot.title = element_text(size = 9, margin = margin(0, 0, 0, 0)))
    if (region_titles == FALSE) {
        p1 = p1 + theme(strip.text = element_blank())
    }
    
    p[[2]]$layers[[1]]$aes_params$size = 0.4
    p[[2]]$layers[[1]]$aes_params$shape = 20
    
    p2 = p[[2]] + scale_colour_manual(values = c(palette[4:1], palette[4 + n])) + 
            theme(legend.position = "none", axis.text.y = element_text(size = 8)) + 
            geom_vline(aes(xintercept = 1.0), linetype = "22", size = 0.25) + labs(x = NULL, y = NULL) + xlim(0, max_x)
    
    list(p1, p2)
}

growth_rate = function(test)
{
    # Extract daily growth rate
    decnov_r = test[t %between% c(306, 353), .(incidence_1 = sum(test_o), incidence_2 = sum(test2_o)), by = .(population, run, t)]
    decnov_r = decnov_r[, .(t = tail(t, -1), r1 = diff(log(incidence_1)), r2 = diff(log(incidence_2))), by = .(population, run)]
    dr_df = decnov_r[, .(r1 = mean(r1), r2 = mean(r2), dr = mean(r2 - r1)), by = .(population, run)]
    
    # Extract difference in Rt
    decnov_Rt = test[t %between% c(306, 353) & group %in% c(11, 12), .(population, run, t, obs0)]
    decnov_Rt[, strain := rep_len(c("Rt1", "Rt2"), .N)]
    decnov_Rt = dcast(decnov_Rt, population + run + t ~ strain, value.var = "obs0")
    dRt_df = decnov_Rt[, .(Rt1 = exp(mean(log(Rt1))), Rt2 = exp(mean(log(Rt2)))), by = .(population, run)]
    dRt_df[, dRt := Rt2 / Rt1]
    
    return (list(dr = dr_df, dRt = dRt_df))
}

grs = function(h, nm)
{
    # For calculating delta-r
    dr_fit = list(post = rbind(
        melt(h$dr, id.vars = c("population", "run")),
        melt(h$dRt, id.vars = c("population", "run"))
    ))
    dr_fit$post[, trial := run]
    dr_fit$post[, chain := 0]
    dr_fit$post[, pop := population]
    dr_fit$post[, population := as.numeric(match(population, nhs_regions))]
    
    ci_r1 = overall(dr_fit, "r1", FALSE, c(1, 3, 9))
    ci_r2 = overall(dr_fit, "r2", FALSE, c(1, 3, 9))
    ci_dr = overall(dr_fit, "dr", FALSE, c(1, 3, 9))
    ci_Rt1 = overall(dr_fit, "Rt1", TRUE, c(1, 3, 9))
    ci_Rt2 = overall(dr_fit, "Rt2", TRUE, c(1, 3, 9))
    ci_dRt = overall(dr_fit, "dRt", TRUE, c(1, 3, 9))
    
    rpw = function(r) {
        100 * (exp(r * 7) - 1)
    }
    
    g = cbind(c("r1", "r2", "dr", "rpw1", "rpw2", "drpw", "Rt1", "Rt2", "dRt"), 
        data.table(rbind(ci_r1, ci_r2, ci_dr, rpw(ci_r1), rpw(ci_r2), rpw(ci_dr), ci_Rt1, ci_Rt2, ci_dRt)))
    names(g) = c("variable", "mean", "q025", "q975")
    g[, name := nm]
    return (g[])
}


sgtf2 = sgtf[date >= "2020-10-01" & nhs_name %in% nhs_regions[which_pops], .(date, population = nhs_name, sgtf, other)]
sgtf2[, c("mean", "lo", "hi") := binom.confint(sgtf, sgtf + other, methods = "bayes", prior.shape1 = 0.01, prior.shape2 = 0.01)[, 4:6]]

# Determine population sizes
nhs_regions = c("East of England", "England", "London", "Midlands", "North East and Yorkshire", "North West", "Northern Ireland", "Scotland", "South East", "South West", "United Kingdom", "Wales")
parametersI = qread("./fits/final/relu.qs")[[2]]
popsize = NULL
for (i in seq_along(parametersI)) {
    if (!is.null(parametersI[[i]])) {
        popsize = rbind(popsize,
            data.table(Geography = parametersI[[i]]$pop[[1]]$name, population_size = sum(parametersI[[i]]$pop[[1]]$size))
        )
    }
}

plot_posterior("./fits/final/relu.qs", "relu", c(1, 3, 4, 5, 6, 9, 10))
plot_posterior("./fits/final/immesc.qs", "immesc", c(1, 3, 9))
plot_posterior("./fits/final/ch_u.qs", "ch_u", c(1, 3, 9))
plot_posterior("./fits/final/infdur.qs", "infdur", c(1, 3, 9))
plot_posterior("./fits/final/serial.qs", "serial", c(1, 3, 9))
plot_posterior("./fits/final/combined.qs", "combined", c(1, 3, 9))

which_pops = england_pops
hu7 = load_hyp("./fits/final/relu.qs")
which_pops = c(1, 3, 9)

hu = load_hyp("./fits/final/relu.qs")    # Note - this just loads for which_pops = c(1, 3, 9) despite the fit file having all NHSE regions
hm = load_hyp("./fits/final/immesc.qs")
hc = load_hyp("./fits/final/ch_u.qs")
hi = load_hyp("./fits/final/infdur.qs")
hs = load_hyp("./fits/final/serial.qs")
ha = load_hyp("./fits/final/combined.qs")

# Assess fit of each hypothesis
fit_dt = data.table(
    hypothesis = c("Transmissibility", "Duration of infectiousness", "Cross protection", 
        "Susceptibility in children", "Generation time", "Combined"),
    dic = c(hu$dic1, hi$dic1, hm$dic1, hc$dic1, hs$dic1, ha$dic1),
    pred = c(predictive(hu), predictive(hi), predictive(hm), predictive(hc), predictive(hs), predictive(ha))
)

fit_dt[, Ddic := dic - min(dic)]
fit_dt[, Dpred := pred - min(pred)]

fwrite(fit_dt, "./output/model_performance.csv")

# Growth rates / R
gr = rbind(
    grs(hu7, "Increased transmissibility (all England)"),
    grs(hu, "Increased transmissibility"),
    grs(hi, "Increased duration of infectiousness"),
    grs(hm, "Immune escape"),
    grs(hc, "Increased susceptibility in children"),
    grs(hs, "Shorter generation time"),
    grs(ha, "Combined model")
)
gr[, name := factor(name, unique(name))]
gr[, variable := factor(variable, unique(variable))]
gr[variable %like% "R", digits := 2]
gr[variable %like% "r", digits := 3]
gr[variable %like% "rpw", digits := 0]
gr[, summary := paste0(round(mean, digits), " (", round(q025, digits), "-", round(q975, digits), ")")]
gr2 = dcast(gr, name ~ variable, value.var = "summary")
fwrite(gr2, "./output/model_growth.csv")

# Best fitting model check params
u_u = overall(hu, "v2_relu", TRUE, c(1, 3, 9))
u_h = overall(hu, "v2_hosp_rlo", FALSE, c(1, 3, 9))
u_i = overall(hu, "v2_icu_rlo", FALSE, c(1, 3, 9))
u_c = overall(hu, "v2_cfr_rlo", FALSE, c(1, 3, 9))

exp(u_h)
exp(u_i)
exp(u_c)

u_u7 = overall(hu7, "v2_relu", TRUE, england_pops)
u_h7 = overall(hu7, "v2_hosp_rlo", FALSE, england_pops)
u_i7 = overall(hu7, "v2_icu_rlo", FALSE, england_pops)
u_c7 = overall(hu7, "v2_cfr_rlo", FALSE, england_pops)

exp(u_h7)
exp(u_i7)
exp(u_c7)

pu = hyp_plots(hu, c(1, 3, 9), "v2_relu", TRUE, "Transmissibility")
pi = hyp_plots(hi, c(1, 3, 9), "v2_infdur", TRUE, "Duration of\ninfectiousness")
pm = hyp_plots(hm, c(1, 3, 9), "v2_immesc", TRUE, "Cross\nprotection")
pc = hyp_plots(hc, c(1, 3, 9), "v2_ch_u", TRUE, "Susceptibility\nin children")
ps = hyp_plots(hs, c(1, 3, 9), "v2_serial", TRUE, "Generation\ntime")
    
su = style(pu, "Increased transmissibility", 1, FALSE, TRUE)
si = style(pi, "Increased duration of infectiousness", 2, FALSE, FALSE)
sm = style(pm, "Immune escape", 3, TRUE, FALSE)
sc = style(pc, "Increased susceptibility in children", 4, FALSE, FALSE)
ss = style(ps, "Shorter generation time", 5, FALSE, FALSE)

se_novoc =  south_east("./fits/novoc_ELSE4.qs", TRUE, "2020-12-31")
se_trans0 = south_east("./fits/final/relu.qs", FALSE, "2020-12-31")
se_trans1 = south_east("./fits/final/relu.qs", TRUE, "2020-12-31")

extenders = data.table(ValueType = c("Deaths", "Hospital\nadmissions", "Hospital beds\noccupied",
    "ICU beds\noccupied", "Infection\nincidence", "PCR\nprevalence (%)", "Seroprevalence\n(%)"), 
    x = rep(as.Date("2020-06-01"), 7), 
    y = c(200, 800, 6000, 500, 40000, 4, 12.5))

se1 = se_trans1$plot + theme(strip.text.x = element_blank()) + geom_blank(data = extenders, aes(x, y))
se0 = se_trans0$plot + theme(axis.title.y = element_blank(), strip.text = element_blank()) + geom_blank(data = extenders, aes(x, y))
sen = se_novoc$plot + theme(axis.title.y = element_blank(), strip.text = element_blank()) + geom_blank(data = extenders, aes(x, y))

plot_hyp = cowplot::plot_grid(
    cowplot::plot_grid(su[[1]], si[[1]], sm[[1]], sc[[1]], ss[[1]], ncol = 1, rel_heights = c(1.15, 1, 1, 1, 1)),
    cowplot::plot_grid(su[[2]] + labs(title = ""), si[[2]], sm[[2]], sc[[2]], ss[[2]], ncol = 1, align = "v", rel_heights = c(1.15, 1, 1, 1, 1)),
    se1, se0, sen, 
    labels = LETTERS, nrow = 1, label_size = 10, rel_widths = c(6, 3.6, 3.6, 3.1, 3.1))

ggsave("./output/plot_hyp.pdf", plot_hyp, width = 35, height = 17, units = "cm", useDingbats = FALSE)
ggsave("./output/plot_hyp.png", plot_hyp, width = 35, height = 17, units = "cm")





# Determine population sizes by age group
nhs_regions = c("East of England", "England", "London", "Midlands", "North East and Yorkshire", "North West", "Northern Ireland", "Scotland", "South East", "South West", "United Kingdom", "Wales")
parametersI = qread("./fits/final/relu.qs")[[2]]
popsize2 = NULL
for (i in seq_along(parametersI)) {
    if (!is.null(parametersI[[i]])) {
        popsize2 = rbind(popsize2,
            data.table(Geography = parametersI[[i]]$pop[[1]]$name, population_size = parametersI[[i]]$pop[[1]]$size, age = seq(0, 75, by = 5))
        )
    }
}


# For age plots
testC = south_east("./fits/final/ch_u.qs")
testU = south_east("./fits/final/infdur.qs")

a_u = overall(ha, "v2_relu", TRUE, c(1, 3, 9))
