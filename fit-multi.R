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

N_THREADS = 24
REP_START = 4
REP_END = 4
BURN_IN = 2500
BURN_IN_FINAL = 2500
ITER = 500

which_pops = c(1, 3, 4, 5, 6, 9, 10)

REF_FIT = "./fits/baseline16.qs"
data_file = "processed-data-2021-01-06.qs"
mobility_file = "schedule3-2021-01-05.rds"

reference_posterior = qread(REF_FIT)[[1]]
reference_sgtf0 = c(0.02571052, 0.03440736, 0.01339611, 0.01631563, 0.01018319, 0.01554367, 0.03621880) # from growth_stat.R
for (i in seq_along(which_pops))
{
    reference_posterior[[which_pops[i]]][, v2_sgtf0 := reference_sgtf0[i]]
}

# Command line
FIT_TYPE = commandArgs(trailingOnly = TRUE)[[1]];

opt_relu = FALSE;
opt_latdur = FALSE;
opt_infdur = FALSE;
opt_immesc = FALSE;
opt_ch_u = FALSE;

if (FIT_TYPE == "relu") {
    extra_priors_1 = list(v2_relu = "L 0.0 0.4 T 0.25 4");
    extra_priors_2 = list(v2_relu = "single");
    opt_relu = TRUE;
} else if (FIT_TYPE == "latdur") {
    extra_priors_1 = list(v2_latdur = "L 0.0 0.4 T 0.01 8");
    extra_priors_2 = list(v2_latdur = "single");
    opt_latdur = TRUE;
} else if (FIT_TYPE == "infdur") {
    extra_priors_1 = list(v2_infdur = "L 0.0 0.4 T 0.01 8");
    extra_priors_2 = list(v2_infdur = "single");
    opt_infdur = TRUE;
} else if (FIT_TYPE == "immesc") {
    extra_priors_1 = list(v2_immesc = "B 2 1");
    extra_priors_2 = list(v2_immesc = "single");
    opt_immesc = TRUE;
} else if (FIT_TYPE == "ch_u") {
    extra_priors_1 = list(v2_ch_u = "N 0 0.25 T 0 2");
    extra_priors_2 = list(v2_ch_u = "single");
    opt_ch_u = TRUE;
} else if (FIT_TYPE == "everything") {
    extra_priors_1 = list(
        v2_relu = "L 0.0 0.4 T 0.25 4", 
        v2_latdur = "L 0.0 0.4 T 0.01 8", 
        v2_infdur = "L 0.0 0.4 T 0.01 8", 
        v2_immesc = "B 2 1", 
        v2_ch_u = "N 0 0.25 T 0 2");
    extra_priors_2 = list(v2_relu = "single", v2_latdur = "single", v2_infdur = "single", v2_immesc = "single", v2_ch_u = "single");
    opt_relu = TRUE;
    opt_latdur = TRUE;
    opt_infdur = TRUE;
    opt_immesc = TRUE;
    opt_ch_u = TRUE;
} else if (FIT_TYPE == "infec_profile") {
    extra_priors_1 = list(
        v2_relu = "L 0.0 0.4 T 0.125 8", 
        v2_latdur = "L 0.0 0.4 T 0.01 8", 
        v2_infdur = "L 0.0 0.4 T 0.01 8");
    extra_priors_2 = list(v2_relu = "single", v2_latdur = "single", v2_infdur = "single");
    opt_relu = TRUE;
    opt_latdur = TRUE;
    opt_infdur = TRUE;
} else {
    stop("Need to specify fit type at command line.");
}



# Adjust number of threads for fitting
N_THREADS = 22 + length(extra_priors_1) * 2;


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

all_data = qread(datapath(data_file))
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
                             date_end = as.character(max(ld$date, sitreps$date, virus$End.date, sero$End.date, sgtf$date) + 1),
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
schedule_all = readRDS(datapath(mobility_file));
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
priorsI = c(list(
    tS = "U 0 60",
    u = "N 0.1 0.01 T 0.04 0.2",
    death_mean = "N 15 2 T 5 30",
    hosp_admission = "N 8 1 T 4 20",
    icu_admission = "N 12.5 1 T 8 14",
    cfr_rlo = "N 0 0.1 T -2 2",
    cfr_rlo2 = "N 0 0.1 T -2 2",
    cfr_rlo3 = "N 0 0.1 T -2 2",
    hosp_rlo = "N 0 0.1 T -2 2", 
    icu_rlo = "N 0 0.1 T -2 2",
    icu_rlo2 = "N 0 0.1 T -2 2",
    contact_final = "N 1 0.1 T 0 1",
    contact_s0 = "E 0.1 0.1",
    contact_s1 = "E 0.1 0.1",
    concentration1 = "N 2 .3 T 2 10",
    concentration2 = "N 2 .2 T 2 10",
    concentration3 = "N 2 .1 T 2 10",
    sep_boost = "N 1 0.05",
    sep_when = "U 214 274",
    disp_deaths = "E 10 10",
    disp_hosp_inc = "E 10 10",
    disp_hosp_prev = "E 10 10",
    disp_icu_prev = "E 10 10",

    v2_when = "U 144 365",
    v2_sgtf0 = "B 1.5 15",
    v2_conc = "E 0.1 0.1 T 2 1000",
    v2_hosp_rlo = "N 0 0.1 T -4 4",
    v2_icu_rlo = "N 0 0.1 T -4 4",
    v2_cfr_rlo = "N 0 0.1 T -4 4"
), extra_priors_1)

priors_multi = c(list(
    tS = "fix",
    u = "fix",
    death_mean = "fix",
    hosp_admission = "fix",
    icu_admission = "fix",
    cfr_rlo = "fix",
    cfr_rlo2 = "fix",
    cfr_rlo3 = "fix",
    hosp_rlo = "fix", 
    icu_rlo = "fix",
    icu_rlo2 = "fix",
    contact_final = "fix",
    contact_s0 = "fix",
    contact_s1 = "fix",
    concentration1 = "fix",
    concentration2 = "fix",
    concentration3 = "fix",
    sep_boost = "fix",
    sep_when = "fix",
    disp_deaths = "fix",
    disp_hosp_inc = "fix",
    disp_hosp_prev = "fix",
    disp_icu_prev = "fix",
    
    v2_when = "multi",
    v2_sgtf0 = "fix",
    v2_conc = "single",
    v2_hosp_rlo = "single",
    v2_icu_rlo = "single",
    v2_cfr_rlo = "single"
), extra_priors_2)

posteriorsI = list()
dynamicsI = list()

# Remove problematic virus entries
virus = virus[omega > 1e-9]

existing_file = paste0("./fits/multi-", FIT_TYPE, REP_START - 1, ".qs");
prev_posterior = NULL;
if (file.exists(existing_file)) {
    saved = qread(existing_file)
    prev_posterior = saved[[1]]
    details = saved[[2]]
    rm(saved)
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
        sgtfI = copy(sgtf);
        sgtfI = sgtfI[pid == p - 1];
    
        # specify user defined functions
        model_v2I = list(
            cpp_changes = cpp_chgI_voc(priorsI, v2 = TRUE, v2_relu = opt_relu, v2_latdur = opt_latdur, v2_infdur = opt_infdur, v2_immesc = opt_immesc, v2_ch_u = opt_ch_u),
            cpp_loglikelihood = cpp_likI_voc(paramsI, ldI, sitrepsI, seroI, virusI, sgtfI, p, "2100-01-01", priorsI, death_cutoff = 7, use_sgtf = TRUE),
            cpp_observer = cpp_obsI_voc(v2 = TRUE, P.death, P.critical, priorsI)
        )

        details$parameters[[pi]] = rlang::duplicate(paramsI);
        details$parameters_translated[[pi]] = rlang::duplicate(cm_translate_parameters(paramsI));
        details$model_v2[[pi]] = rlang::duplicate(model_v2I);
        pi = pi + 1;
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
            if (!is.null(prev_posterior) && nm %in% names(prev_posterior)) {
                values = prev_posterior[[nm]];
                range = quantile(values, c(0.025, 0.975));
                priors_new[[nm]] = paste0(priorsI[[nm]], " I ", range[1], " ", range[2], " ", init_previous_amount);
                cat("Restarting from ", priors_new[[nm]], " for ", nm, "\n");
            } else if (!is.null(prev_posterior) && all(paste0(nm, seq_along(which_pops) - 1) %in% names(prev_posterior))) {
                for (j in seq_along(which_pops)) {
                    values = prev_posterior[[paste0(nm, j - 1)]];
                    range = quantile(values, c(0.025, 0.975));
                    priors_new[[nm]][j] = paste0(priorsI[[nm]], " I ", range[1], " ", range[2], " ", init_previous_amount);
                    cat("Restarting from ", priors_new[[nm]][j], " for ", nm, "\n");
                }
            } else if (nm %in% names(reference_posterior[[which_pops[1]]])) {
                if (pr == "single") {
                    values = NULL;
                    for (j in which_pops) {
                        values = c(values, reference_posterior[[j]][[nm]]);
                    }
                    range = quantile(values, c(0.025, 0.975));
                    priors_new[[nm]] = paste0(priorsI[[nm]], " I ", range[1], " ", range[2], " ", init_previous_amount);
                    cat("Starting from ", priors_new[[nm]], " for ", nm, "\n");
                } else if (pr == "multi") {
                    for (j in seq_along(which_pops)) {
                        values = reference_posterior[[which_pops[j]]][[nm]];
                        range = quantile(values, c(0.025, 0.975));
                        priors_new[[nm]][j] = paste0(priorsI[[nm]], " I ", range[1], " ", range[2], " ", init_previous_amount);
                        cat("Starting from ", priors_new[[nm]][j], " for ", nm, "\n");
                    }
                } else {
                    stop("Cannot interpret ", pr);
                }
            } else {
                if (pr == "single") {
                    priors_new[[nm]] = priorsI[[nm]];
                    cat("Initialising from ", priors_new[[nm]], " for ", nm, "\n");
                } else if (pr == "multi") {
                    for (j in seq_along(which_pops)) {
                        priors_new[[nm]][j] = priorsI[[nm]];
                        cat("Initialising from ", priors_new[[nm]][j], " for ", nm, "\n");
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

    qsave(rlang::duplicate(list(postI, details)), paste0("./fits/multi-", FIT_TYPE, replic, ".qs"))
    
    # Sampling fits
    test = cm_backend_sample_fit_multi(details$parameters_translated, details$priors, postI, 20, seed = 0)
    test = rbindlist(test, idcol = "simulation")
    test[, population := as.character(details$pops[(simulation - 1) %% length(details$pops) + 1])]
    rows = cm_backend_sample_fit_rows(details$parameters_translated[[1]], postI, 20, seed = 0);
    
    # Add dispersion parameters
    disp = NULL
    for (p in which_pops) {
        disp_p = as.data.table(as.list(reference_posterior[[p]][, colMeans(.SD), .SDcols = patterns("^disp")]))
        disp_p[, population := nhs_regions[p]]
        disp = rbind(disp, disp_p);
    }
    test = merge(test, disp, by = "population")
    
    # Fit to SGTF data ###
    sgtf[, qlo := qbeta(0.025, sgtf + 1, other + 1)]
    sgtf[, qhi := qbeta(0.975, sgtf + 1, other + 1)]
    #vmodel = test[, .(I1 = sum(Ip + Is + Ia), I2 = sum(Ip2 + Is2 + Ia2)), by = .(t, population, run)]
    vmodel = test[, .(I1 = sum(test_o), I2 = sum(test2_o)), by = .(t, population, run)]
    vmodel[, p2 := I2 / (I1 + I2)]
    vmodel[is.nan(p2), p2 := 0]
    sgtf0 = data.table(population = details$pops, f0 = reference_sgtf0)
    vmodel = merge(vmodel, sgtf0, by = "population")
    vmodel[, sgtf := (1 - p2) * f0 + p2];
    conc = data.table(v2_conc = as.data.table(postI)[rows, v2_conc])
    conc[, run := 1:.N]
    vmodel = merge(vmodel, conc, by = "run")
    vmodel[, alpha := sgtf * (v2_conc - 2) + 1]
    vmodel[, beta := (1 - sgtf) * (v2_conc - 2) + 1]
    vmodel[, q025 := qbeta(0.025, alpha, beta)]
    vmodel[, q500 := qbeta(0.500, alpha, beta)]
    vmodel[, q975 := qbeta(0.975, alpha, beta)]
    vmodel = vmodel[, lapply(.SD, mean), .SDcols = c("q025", "q500", "q975"), by = .(nhs_name = population, t)]
    plotS = ggplot(sgtf[(pid + 1) %in% which_pops]) +
        geom_ribbon(aes(x = date, ymin = qlo, ymax = qhi), fill = "black", alpha = 0.1) +
        geom_ribbon(data = vmodel[t + ymd("2020-01-01") >= "2020-10-01"], 
            aes(x = ymd("2020-01-01") + t, ymin = q025, ymax = q975), fill = "darkorchid", alpha = 0.5) +
        geom_line(data = vmodel[t + ymd("2020-01-01") >= "2020-10-01"], 
            aes(x = ymd("2020-01-01") + t, y = q500), colour = "darkorchid") +
        geom_line(aes(x = date, y = sgtf / (sgtf + other)), size = 0.25) +
        #scale_y_continuous(trans = scales::logit_trans(), breaks = c(0.01, 0.1, 0.2, 0.3, 0.5, 0.7, 0.8, 0.9, 0.99), limits = c(0.01, 0.99)) +
        facet_wrap(~nhs_name) +
        labs(x = NULL, y = "Relative frequency of\nS gene target failure") +
        scale_x_date(date_breaks = "1 month", date_labels = "%b")
    ggsave(paste0("./output/sgtf_check_", FIT_TYPE, replic, ".pdf"), plotS, width = 20, height = 6, units = "cm", useDingbats = FALSE)
    #ggsave(paste0("./output/sgtf_check_", replic, ".png"), plotS, width = 20, height = 6, units = "cm")
    
    # # Posteriors of interest
    # post[, D := -2 * ll]
    # post[, 0.5 * var(D) + mean(D), by = population][, mean(V1)]
    # 
    # post[, population := nhs_regions[which_pops[population]]]
    # post = melt(post, id.vars = 1, measure.vars = c("v2_relu", "v2_hosp_rlo", "v2_cfr_rel"))
    # post[variable == "v2_relu", variable := "Relative transmission rate"]
    # post[variable == "v2_hosp_rlo", value := exp(value)]
    # post[variable == "v2_hosp_rlo", variable := "Associated OR of hospitalisation"]
    # post[variable == "v2_cfr_rel", variable := "Associated RR of death"]
    # 
    # prior = data.table(x = rep(seq(0.7, 1.3, 0.03), 3), 
    #     variable = rep(c("Relative transmission rate", "Associated OR of hospitalisation", "Associated RR of death"), each = 21))
    # prior[variable == "Relative transmission rate", y := dlnorm(x, 0, 0.2)]
    # prior[variable == "Associated OR of hospitalisation", y := dlnorm(x, 0, 0.1)]
    # prior[variable == "Associated RR of death", y := dnorm(x, 1, 0.1)]
    # plot2 = ggplot(post) +
    #     geom_line(data = prior, aes(x, y), colour = "#888888") +
    #     geom_density(aes(value, colour = population), adjust = 4) +
    #     #geom_histogram(aes(value, colour = population), bins = 20) +
    #     facet_wrap(~variable, scales = "free") +
    #     theme(legend.position = c(0.06, 0.9)) +
    #     labs(x = NULL, y = NULL, colour = NULL) +
    #     expand_limits(x = 1)
    # ggsave(paste0("./output/variant_stats_", replic, ".pdf"), plot2, width = 20, height = 6, units = "cm", useDingbats = FALSE)
    # 
    # qsave(plotS, "./output/sgtf-plot-0.qs")
    # qsave(plot2, "./output/post-plot-0.qs")
    # 
    # post[variable == "Relative transmission rate", quantile(value, c(0.025, 0.5, 0.975))]

    # Visually inspect fit
    #plot_a = check_fit(test0, ld, sitreps, virus, sero, nhs_regions[which_pops], "2021-01-02")
    plot_b = check_fit(test, details$parameters, ld, sitreps, virus, sero, nhs_regions[which_pops], death_cutoff = 7, "2021-01-02")
    #plot3 = cowplot::plot_grid(plot_a, plot_b, nrow = 1, labels = LETTERS)
    ggsave(paste0("./output/multi_fit_", FIT_TYPE, replic, ".pdf"), plot_b, width = 30, height = 25, units = "cm", useDingbats = FALSE)

    #plot_a = check_fit(test0, ld, sitreps, virus, sero, nhs_regions[which_pops], "2021-01-02", "2020-09-01")
    plot_b = check_fit(test, details$parameters, ld, sitreps, virus, sero, nhs_regions[which_pops], death_cutoff = 7, "2021-01-02", "2020-09-01")
    #plot3L = cowplot::plot_grid(plot_a, plot_b, nrow = 1, labels = LETTERS)
    ggsave(paste0("./output/multi_fitL_", FIT_TYPE, replic, ".pdf"), plot_b, width = 30, height = 25, units = "cm", useDingbats = FALSE)
}


stop("Add back in measured sgtf0")

