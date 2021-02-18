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

# Lines with updated data: ###
# Need to apply these also to other fits.R versions.

N_THREADS = 32
REP_START = 19
REP_END = 19
BURN_IN = 2500
BURN_IN_FINAL = 2500
ITER = 500
REF_FIT = "./fits/pp18.qs"

reference_posterior = qread(REF_FIT)[[1]]

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

cm_populations = rbind(cm_populations[name != "United Kingdom"], popUK)
cm_matrices = c(cm_matrices, matricesUK)
source("./distribution_fit.R");
source("./spim_output.R");
source("./check_fit.R")


#
# DATA
#

nhs_regions = popUK[, unique(name)]
pct = function(x) as.numeric(str_replace_all(x, "%", "")) / 100

all_data = qread(datapath("processed-data-2021-01-01.qs")) ###
ld = all_data[[1]]
sitreps = all_data[[2]]
virus = all_data[[3]][!Data.source %like% "7a|7b|6a|6b"]
sero = all_data[[4]]
sgtf = all_data[[5]]

#
# FITTING
#

# NUMBER OF REGIONS TO FIT
N_REG = 12;

# Build parameters for NHS regions ###
params = cm_parameters_SEI3R(nhs_regions[1:N_REG], deterministic = T, 
                             date_start = "2020-01-01", 
                             date_end = as.character(max(ld$date, sitreps$date, virus$End.date, sero$End.date, sgtf$date) + 1), ###
    dE  = cm_delay_gamma(2.5, 2.5, t_max = 15, t_step = 0.25)$p,
    dIp = cm_delay_gamma(2.5, 4.0, t_max = 15, t_step = 0.25)$p,
    dIs = cm_delay_gamma(2.5, 4.0, t_max = 15, t_step = 0.25)$p,
    dIa = cm_delay_gamma(5.0, 4.0, t_max = 15, t_step = 0.25)$p)
params = cm_split_matrices_ex_in(params, 15)

# school terms
school_close =  c("2020-2-16", "2020-4-05", "2020-5-24", "2020-7-22", "2020-10-25", "2020-12-20", "2021-02-14", "2021-04-01", "2021-05-30", "2021-07-25");
school_reopen = c("2020-2-22", "2020-4-18", "2020-5-30", "2020-9-01", "2020-10-31", "2021-01-02", "2021-02-20", "2021-04-17", "2021-06-05", "2021-09-01");

# Load age-varying symptomatic rate
covid_scenario = qread(datapath("2-linelist_both_fit_fIa0.5-rbzvih.qs"));
covu = unname(rep(colMeans(covid_scenario[,  5:12]), each = 2));
covy = unname(rep(colMeans(covid_scenario[, 13:20]), each = 2));

for (i in seq_along(params$pop)) {
    params$pop[[i]]$u = covu / mean(covu);
    params$pop[[i]]$u2 = covu / mean(covu);
    params$pop[[i]]$y = covy;
    params$pop[[i]]$y2 = covy;
}

# Health burden processes
source("./processes.R")
params$processes = burden_processes

# changes
schedule_all = readRDS(datapath("schedule3-2020-12-31.rds")); ###
schedule = list();
for (i in seq_along(schedule_all)) {
    if (schedule_all[[i]]$pops < N_REG) {
        schedule[[length(schedule) + 1]] = schedule_all[[i]]
    }
}

# Remove NAs
for (i in seq_along(schedule)) {
    for (j in seq_along(schedule[[i]]$values)) {
        if (any(is.na(schedule[[i]]$values[[j]]))) {
            schedule[[i]]$values[[j]] = ifelse(is.na(schedule[[i]]$values[[j]]), prev, schedule[[i]]$values[[j]])
        }
        prev = schedule[[i]]$values[[j]];
    }
}
params$schedule = schedule


#
# Individual fits
#

source("./cpp_funcs.R")

# Fitting
priorsI = list(
    tS = "U 0 60",
    u = "N 0.07 0.01 T 0.04 0.2",
    death_mean = "N 15 2 T 5 30",    # <<< co-cin
    death_shape = "N 1.9 0.2 T 0.1 3", # <<< co-cin
    admission = "N 8 1 T 4 20", # <<< co-cin
    cfr_rel = "N 1 0.1 T 0.1 4", # <<< co-cin
    icu_rlo = "N 0 0.1 T -2 2",
    hosp_rlo = "N 0 0.1 T -2 2", 
    icu_admission = "N 12.5 1 T 8 14", # <<< co-cin
    contact_final = "N 1 0.1 T 0 1",
    contact_s0 = "E 0.1 0.1",
    contact_s1 = "E 0.1 0.1",
    concentration1 = "N 2 .3 T 2 10", # was .5
    concentration2 = "N 2 .2 T 2 10", # was .4
    concentration3 = "N 2 .1 T 2 10", # was .2
    sep_boost = "N 1 0.05",
    cfr_rel2 = "N 0.45 0.1 T 0 1", # <<<
    sep_when = "U 224 264",
    v2_when = "U 144 365",
    v2_relu = "L 0.0 0.4 T 0.25 4",
    v2_hosp_rlo = "N 0 0.1 T -2 2", # hosp x[20]
    v2_cfr_rel = "N 1 0.1 T 0.1 4", # cfr_rel x[21]
    v2_conc = "E 0.1 0.1 T 2 1000", # x[22] conc ###
    v2_sgtf0 = "B 1.5 15",  # x[23] sgtf0 ###
    v2_infdur = "L 0.0 0.4 T 0.25 4", # x[24]
    test_delay = "N 2.5 5 T 0 15", # x[25]
    v2_immesc = "B 2 1", # x[26]
    v2_ch_u = "N 0 0.25 T 0 1"
)

priors_multi = list(
    tS = "fix",
    u = "fix",
    death_mean = "fix",
    death_shape = "fix",
    admission = "fix",
    cfr_rel = "fix",
    icu_rlo = "fix",
    hosp_rlo = "fix",
    icu_admission = "fix",
    contact_final = "fix",
    contact_s0 = "fix",
    contact_s1 = "fix",
    concentration1 = "fix",
    concentration2 = "fix",
    concentration3 = "fix",
    sep_boost = "fix",
    cfr_rel2 = "fix",
    sep_when = "fix",
    v2_when = "multi",
    v2_relu = "single",
    v2_hosp_rlo = "single",
    v2_cfr_rel = "single",
    v2_conc = "single",
    v2_sgtf0 = "fix",
    v2_infdur = "single",
    test_delay = "single",
    v2_immesc = "single",
    v2_ch_u = "single"
)

posteriorsI = list()
dynamicsI = list()
parametersI = list()

# Remove problematic virus entries
virus = virus[omega > 1e-9]

existing_file = paste0("./fits/pp", REP_START - 1, ".qs");
if (file.exists(existing_file)) {
    saved = qread(existing_file)
    posteriorsI = saved[[1]]
    parametersI = saved[[2]]
    rm(saved)
    
    # Ensure newly added parameters are in parametersI
    for (p in seq_along(parametersI))
    {
        if (!is.null(parametersI[[p]])) {
            if (is.null(parametersI[[p]]$pop[[1]]$ed_vi)) {
                parametersI[[p]]$pop[[1]]$ed_vi = rep(0, 16);
            }
            if (is.null(parametersI[[p]]$pop[[1]]$ed_vi2)) {
                parametersI[[p]]$pop[[1]]$ed_vi2 = rep(0, 16);
            }
        }
    }
}

for (replic in REP_START:REP_END)
{
    init_previous = TRUE
    init_previous_amount = 1
    
    # Loop through regions
    pi = 1;
    details = list();
    for (p in which_pops) {
        paramsI = rlang::duplicate(params);
        paramsI$pop = list(rlang::duplicate(params$pop[[p]]));
        paramsI$travel = matrix(1, nrow = 1, ncol = 1);
        paramsI$schedule = list();
        j = 1;
        for (i in seq_along(params$schedule)) {
            if (p - 1 == params$schedule[[i]]$pops) {
                paramsI$schedule[[j]] = rlang::duplicate(params$schedule[[i]]);
                paramsI$schedule[[j]]$pops = 0;
                j = j + 1;
            }
        }
    
        # contact placeholder for tier 2
        paramsI$schedule[[3]] = rlang::duplicate(paramsI$schedule[[1]]);
        for (i in seq_along(paramsI$schedule[[3]]$values)) {
            paramsI$schedule[[3]]$values[[i]][1] = paramsI$schedule[[1]]$values[[i]][1] +  0.2497655 / 100;
            paramsI$schedule[[3]]$values[[i]][2] = paramsI$schedule[[1]]$values[[i]][2] + -0.2307939 / 100;
            paramsI$schedule[[3]]$values[[i]][3] = paramsI$schedule[[1]]$values[[i]][3] + -1.5907698 / 100;
            paramsI$schedule[[3]]$values[[i]][4] = paramsI$schedule[[1]]$values[[i]][4] + -3.4866544 / 100;
            paramsI$schedule[[3]]$values[[i]][5] = paramsI$schedule[[1]]$values[[i]][5] + -3.4524518 / 100;
        }
        paramsI$schedule[[3]]$mode = "bypass";
    
        # contact placeholder for tier 3
        paramsI$schedule[[4]] = rlang::duplicate(paramsI$schedule[[1]]);
        for (i in seq_along(paramsI$schedule[[4]]$values)) {
            paramsI$schedule[[4]]$values[[i]][1] = paramsI$schedule[[1]]$values[[i]][1] +  2.080457 / 100;
            paramsI$schedule[[4]]$values[[i]][2] = paramsI$schedule[[1]]$values[[i]][2] + -8.045226 / 100;
            paramsI$schedule[[4]]$values[[i]][3] = paramsI$schedule[[1]]$values[[i]][3] + -2.476266 / 100;
            paramsI$schedule[[4]]$values[[i]][4] = paramsI$schedule[[1]]$values[[i]][4] + -10.144043 / 100;
            paramsI$schedule[[4]]$values[[i]][5] = paramsI$schedule[[1]]$values[[i]][5] + -7.681244 / 100;
        }
        paramsI$schedule[[4]]$mode = "bypass";
    
        # contact multiplier for gradual contact change
        paramsI$schedule[[5]] = list(
            parameter = "contact",
            pops = 0,
            mode = "multiply",
            values = rep(list(rep(1, 8)), 366),
            times = 0:365
        )
    
        # contact multiplier for september boost
        paramsI$schedule[[6]] = list(
            parameter = "contact",
            pops = 0,
            mode = "multiply",
            values = list(rep(1, 8)),
            times = c(244)
        )
    
        ldI = copy(ld);
        ldI = ldI[pid == p - 1];
        sitrepsI = copy(sitreps);
        sitrepsI = sitrepsI[pid == p - 1];
        seroI = copy(sero);
        seroI = seroI[pid == p - 1 & Data.source != "NHSBT"];   # sero: all but NHSBT
        virusI = copy(virus);
        virusI = virusI[pid == p - 1 & Data.source %like% "REACT"]; # virus: REACT only
        sgtfI = copy(sgtf); ###
        sgtfI = sgtfI[pid == p - 1]; ###
    
        # specify user defined functions
        model_v2I = list(
            cpp_changes = cpp_chgI_everything(TRUE),
            cpp_loglikelihood = cpp_likI(paramsI, ldI, sitrepsI, seroI, virusI, sgtfI, p),
            cpp_observer = cpp_obsI(P.death)
        )

        details$parameters[[pi]] = rlang::duplicate(paramsI);
        details$parameters_translated[[pi]] = rlang::duplicate(cm_translate_parameters(paramsI));
        details$model_v2[[pi]] = rlang::duplicate(model_v2I);
        pi = pi + 1;
    
        # postI = cm_backend_mcmc_test(cm_translate_parameters(paramsI), priorsI2,
        #     seed = 0, burn_in = ifelse(replic == REP_END, BURN_IN_FINAL, BURN_IN), 
        #     iterations = ITER, n_threads = N_THREADS, classic_gamma = T);
        # setDT(postI)
        # posteriorsI[[p]] = postI
        # 
        # parametersI[[p]] = rlang::duplicate(paramsI)
        # qsave(rlang::duplicate(list(posteriorsI, parametersI)), paste0("./fits/pp", replic, "-progress.qs"))
        # 
        # print(p)
    }
    
    # load user defined functions
    details$pops = nhs_regions[which_pops];
    cm_source_backend(
        user_defined = list(
            model_v2 = list(
                cpp_changes = wrap_region(details$model_v2, "cpp_changes", details$pops),
                cpp_loglikelihood = wrap_region(details$model_v2, "cpp_loglikelihood", details$pops),
                cpp_observer = wrap_region(details$model_v2, "cpp_observer", details$pops)
            )
        )
    )
    
    priors_new = list();
    for (i in seq_along(priors_multi))
    {
        nm = names(priors_multi)[i];
        pr = priors_multi[[i]];
        
        if (pr == "fix") { # Fixed parameter
            posterior_mean = rep(0, length(which_pops));
            for (j in seq_along(which_pops))
            {
                posterior_mean[j] = mean(reference_posterior[[which_pops[j]]][[nm]])
            }
            priors_new[[nm]] = posterior_mean;
        } else { # Fitted parameter
            if (nm %in% names(reference_posterior[[which_pops[1]]])) {
                if (pr == "single") {
                    values = NULL;
                    for (j in which_pops) {
                        values = c(values, reference_posterior[[j]][[nm]]);
                    }
                    range = quantile(values, c(0.025, 0.975));
                    priors_new[[nm]] = paste0(priorsI[[nm]], " I ", range[1], " ", range[2], " ", init_previous_amount);
                } else if (pr == "multi") {
                    for (j in seq_along(which_pops)) {
                        values = reference_posterior[[which_pops[j]]][[nm]];
                        range = quantile(values, c(0.025, 0.975));
                        priors_new[[nm]][j] = paste0(priorsI[[nm]], " I ", range[1], " ", range[2], " ", init_previous_amount);
                    }
                } else {
                    stop("Cannot interpret ", pr);
                }
            } else {
                if (pr == "single") {
                    priors_new[[nm]] = priorsI[[nm]];
                } else if (pr == "multi") {
                    for (j in seq_along(which_pops)) {
                        priors_new[[nm]][j] = priorsI[[nm]];
                    }
                } else {
                    stop("Cannot interpret ", pr);
                }
            }
        }
    }
    
    details$priors = priors_new

    postI = cm_backend_mcmc_multi(details$parameters_translated, details$priors, 
        seed = 0, burn_in = ifelse(replic == REP_END, BURN_IN_FINAL, BURN_IN), 
        iterations = ITER, n_threads = N_THREADS, classic_gamma = T)

    qsave(rlang::duplicate(list(postI, details)), paste0("./fits/multiEverything", replic, ".qs"))
    
    # Generate SPI-M output
    # Sample dynamics from fit
    # saved = qread("./fits/multiB19.qs")
    # postI = saved[[1]]
    # details = saved[[2]]
    # rm(saved)
    ##dynamicsI = list()
    ##dynamics0 = list()

    # # Source backend
    # cm_source_backend(
    #     user_defined = list(
    #         model_v2 = list(
    #             cpp_changes = wrap_region(details$model_v2, "cpp_changes", details$pops),
    #             cpp_loglikelihood = wrap_region(details$model_v2, "cpp_loglikelihood", details$pops),
    #             cpp_observer = wrap_region(details$model_v2, "cpp_observer", details$pops)
    #         )
    #     )
    # )
    
    # Sampling fits
    #paramsI2 = rlang::duplicate(parametersI[[p]])
    #paramsI2$time1 = as.character(ymd(parametersI[[p]]$time1) + 56);
    test = cm_backend_sample_fit_multi(details$parameters_translated, details$priors, postI, 20, seed = 0)
    test = rbindlist(test, idcol = "simulation")
    test[, population := as.character(details$pops[(simulation - 1) %% length(details$pops) + 1])]
    
    # Again, without new variant
    posteriors0 = copy(posteriorsI[[p]]);
    posteriors0[, v2_when := 9999];
    test0 = cm_backend_sample_fit_test(cm_translate_parameters(paramsI2), posteriors0, 100, seed = 0);
    
    #dynamics0[[p]] = test0

    # Concatenate dynamics for SPI-M output
    #test = rbindlist(dynamicsI, fill = TRUE)
    #test[, population := nhs_regions[population]]
    #test0 = rbindlist(dynamics0, fill = TRUE)
    #test0[, population := nhs_regions[population]]
    
    # Fit to SGTF data ###
    sgtf[, qlo := qbeta(0.025, sgtf + 1, other + 1)]
    sgtf[, qhi := qbeta(0.975, sgtf + 1, other + 1)]
    #vmodel = test[, .(I1 = sum(Ip + Is + Ia), I2 = sum(Ip2 + Is2 + Ia2)), by = .(t, population, run)]
    vmodel = test[, .(I1 = sum(test_o), I2 = sum(test2_o)), by = .(t, population, run)]
    vmodel[is.nan(p2), p2 := 0]
    sgtf0 = data.table(population = details$pops, f0 = details$priors$v2_sgtf0)
    vmodel = merge(vmodel, sgtf0, by = "population")
    vmodel = vmodel[, as.list(quantile(p2, c(0.025, 0.5, 0.975))), by = .(t, nhs_name = population)]
    plotS = ggplot(sgtf[(pid + 1) %in% which_pops]) +
        geom_ribbon(aes(x = date, ymin = qlo, ymax = qhi), fill = "black", alpha = 0.1) +
        geom_ribbon(data = vmodel[t + ymd("2020-01-01") >= "2020-10-01"], 
            aes(x = ymd("2020-01-01") + t, ymin = `2.5%`, ymax = `97.5%`), fill = "darkorchid", alpha = 0.5) +
        geom_line(aes(x = date, y = sgtf / (sgtf + other)), size = 0.25) +
        #scale_y_continuous(trans = scales::logit_trans(), breaks = c(0.01, 0.1, 0.2, 0.3, 0.5, 0.7, 0.8, 0.9, 0.99), limits = c(0.01, 0.99)) +
        facet_wrap(~nhs_name) +
        labs(x = NULL, y = "Relative frequency of\nS gene target failure") +
        scale_x_date(date_breaks = "1 month", date_labels = "%b")
    ggsave(paste0("./output/sgtf_check_", replic, ".pdf"), plotS, width = 20, height = 6, units = "cm", useDingbats = FALSE)
    #ggsave(paste0("./output/sgtf_check_", replic, ".png"), plotS, width = 20, height = 6, units = "cm")
    
    # # Fit to COG data ###
    # variant[, qlo := qbeta(0.025, var2 + 1, all - var2 + 1)]
    # variant[, qhi := qbeta(0.975, var2 + 1, all - var2 + 1)]
    # vmodel = test[, .(p2 = sum(Ip2 + Is2 + Ia2) / sum(Ip + Is + Ia + Ip2 + Is2 + Ia2)), by = .(t, population, run)]
    # vmodel[is.nan(p2), p2 := 0]
    # vmodel = vmodel[, as.list(quantile(p2, c(0.025, 0.5, 0.975))), by = .(t, nhs_name = population)]
    # plot1 = ggplot(variant[(pid + 1) %in% which_pops]) +
    #     geom_ribbon(aes(x = sample_date, ymin = qlo, ymax = qhi), fill = "black", alpha = 0.1) +
    #     geom_ribbon(data = vmodel[t + ymd("2020-01-01") >= "2020-10-01"], 
    #         aes(x = ymd("2020-01-01") + t, ymin = `2.5%`, ymax = `97.5%`), fill = "darkorchid", alpha = 0.5) +
    #     geom_line(aes(x = sample_date, y = var2 / all), size = 0.25) +
    #     facet_wrap(~nhs_name) +
    #     labs(x = NULL, y = "Relative frequency of\nVOC 202012/01") +
    #     scale_x_date(date_breaks = "1 month", date_labels = "%b")
    # ggsave(paste0("./output/variant_check_", replic, ".pdf"), plot1, width = 20, height = 6, units = "cm", useDingbats = FALSE)
    # #ggsave(paste0("./output/variant_check_", replic, ".png"), plot1, width = 20, height = 6, units = "cm")
    
    # Posteriors of interest
    post = rbindlist(posteriorsI[which_pops], idcol = "population")
    post[, D := -2 * ll]
    post[, 0.5 * var(D) + mean(D), by = population][, mean(V1)]
    
    post[, population := nhs_regions[which_pops[population]]]
    post = melt(post, id.vars = 1, measure.vars = c("v2_relu", "v2_hosp_rlo", "v2_cfr_rel"))
    post[variable == "v2_relu", variable := "Relative transmission rate"]
    post[variable == "v2_hosp_rlo", value := exp(value)]
    post[variable == "v2_hosp_rlo", variable := "Associated OR of hospitalisation"]
    post[variable == "v2_cfr_rel", variable := "Associated RR of death"]
    
    prior = data.table(x = rep(seq(0.7, 1.3, 0.03), 3), 
        variable = rep(c("Relative transmission rate", "Associated OR of hospitalisation", "Associated RR of death"), each = 21))
    prior[variable == "Relative transmission rate", y := dlnorm(x, 0, 0.2)]
    prior[variable == "Associated OR of hospitalisation", y := dlnorm(x, 0, 0.1)]
    prior[variable == "Associated RR of death", y := dnorm(x, 1, 0.1)]
    plot2 = ggplot(post) +
        geom_line(data = prior, aes(x, y), colour = "#888888") +
        geom_density(aes(value, colour = population), adjust = 4) +
        #geom_histogram(aes(value, colour = population), bins = 20) +
        facet_wrap(~variable, scales = "free") +
        theme(legend.position = c(0.06, 0.9)) +
        labs(x = NULL, y = NULL, colour = NULL) +
        expand_limits(x = 1)
    ggsave(paste0("./output/variant_stats_", replic, ".pdf"), plot2, width = 20, height = 6, units = "cm", useDingbats = FALSE)
    #ggsave(paste0("./output/variant_stats_", replic, ".png"), plot2, width = 20, height = 6, units = "cm")

    qsave(plotS, "./output/sgtf-plot-0.qs")
    # qsave(plot1, "./output/cog-plot-0.qs") ###
    qsave(plot2, "./output/post-plot-0.qs")
    
    post[variable == "Relative transmission rate", quantile(value, c(0.025, 0.5, 0.975))]

    # Visually inspect fit
    plot_a = check_fit(test0, ld, sitreps, virus, sero, nhs_regions[which_pops], "2021-01-02")
    plot_b = check_fit(test, ld, sitreps, virus, sero, nhs_regions[which_pops], "2021-01-02")
    plot3 = cowplot::plot_grid(plot_a, plot_b, nrow = 1, labels = LETTERS)
    ggsave(paste0("./output/fit_", replic, ".pdf"), plot3, width = 30, height = 25, units = "cm", useDingbats = FALSE)
    #ggsave(paste0("./output/fit_", replic, ".png"), plot3, width = 30, height = 25, units = "cm")

    plot_a = check_fit(test0, ld, sitreps, virus, sero, nhs_regions[which_pops], "2021-01-02", "2020-09-01")
    plot_b = check_fit(test, ld, sitreps, virus, sero, nhs_regions[which_pops], "2021-01-02", "2020-09-01")
    plot3L = cowplot::plot_grid(plot_a, plot_b, nrow = 1, labels = LETTERS)
    ggsave(paste0("./output/fitL_", replic, ".pdf"), plot3L, width = 30, height = 25, units = "cm", useDingbats = FALSE)
    #ggsave(paste0("./output/fitL_", replic, ".png"), plot3L, width = 30, height = 25, units = "cm")
    
    # england_pops = c(1, 3, 4, 5, 6, 9, 10)
    # plot = compare_fit(test, test0, ld, sitreps, virus, sero, nhs_regions[which_pops], nhs_regions[england_pops], "2020-12-18")
    # ggsave(paste0("./output/fitboth_", replic, ".pdf"), plot, width = 30, height = 25, units = "cm", useDingbats = FALSE)
    # ggsave(paste0("./output/fitboth_", replic, ".png"), plot, width = 30, height = 25, units = "cm")
}
