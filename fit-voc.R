# fit-voc: For VOC fits.

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

theme_set(cowplot::theme_cowplot(font_size = 10) + theme(strip.background = element_blank()))

N_THREADS = 46
REP_START = 1
REP_END = 2
BURN_IN = 3000
ITER = 5000
BURN_IN_FINAL = 3000
ITER_FINAL = 5000

which_pops = c(1, 3, 4, 5, 6, 9, 10)
set_id = ""

data_file = "processed-data-2021-01-18.qs"
mobility_file = "schedule3-2021-01-07.rds"
date_fitting = "2020-12-24"
# data_file = "processed-data-2021-01-15.qs"
# mobility_file = "schedule3-2021-01-14.rds"
# date_fitting = "2021-01-15"


# Command line
c_args = commandArgs(trailingOnly = TRUE);
FIT_TYPE = c_args[[1]];
POP_SET = c_args[[2]];
if (length(c_args) > 2) {
    REP_START = as.numeric(c_args[[3]]);
    REP_END = as.numeric(c_args[[4]]);
}

opt_conc = TRUE;
opt_seas = FALSE;

opt_v2 = TRUE;
opt_relu = FALSE;
opt_latdur = FALSE;
opt_serial = FALSE;
opt_infdur = FALSE;
opt_immesc = FALSE;
opt_ch_u = FALSE;
extra_priors = list();

if (FIT_TYPE == "relu") {
    extra_priors = list(v2_relu = "L 0.0 0.4 T 0.25 4");
    opt_relu = TRUE;
} else if (FIT_TYPE == "latdur") {
    extra_priors = list(v2_latdur = "L 0.0 0.4 T 0.01 8");
    opt_latdur = TRUE;
} else if (FIT_TYPE == "serial") {
    extra_priors = list(v2_serial = "L 0.0 0.4 T 0.01 1");
    opt_serial = TRUE;
} else if (FIT_TYPE == "infdur") {
    extra_priors = list(v2_infdur = "L 0.0 0.4 T 0.01 8");
    opt_infdur = TRUE;
} else if (FIT_TYPE == "immesc") {
    extra_priors = list(v2_immesc = "B 3 1");
    opt_immesc = TRUE;
} else if (FIT_TYPE == "ch_u") {
    extra_priors = list(v2_ch_u = "L 0.0 0.4 T 0.04 24");
    opt_ch_u = TRUE;
} else if (FIT_TYPE == "combined") {
    extra_priors = list(
        v2_relu = "L 0.0 0.4 T 0.25 4", 
        v2_serial = "L 0.0 0.4 T 0.01 8", 
        v2_immesc = "B 3 1", 
        v2_ch_u = "L 0.0 0.4 T 0.04 24");
    opt_relu = TRUE;
    opt_serial = TRUE;
    opt_immesc = TRUE;
    opt_ch_u = TRUE;
} else if (FIT_TYPE == "everything") {
    extra_priors = list(
        v2_relu = "L 0.0 0.4 T 0.25 4", 
        v2_latdur = "L 0.0 0.4 T 0.01 8", 
        v2_infdur = "L 0.0 0.4 T 0.01 8", 
        v2_immesc = "B 3 1", 
        v2_ch_u = "L 0.0 0.4 T 0.04 24");
    opt_relu = TRUE;
    opt_latdur = TRUE;
    opt_infdur = TRUE;
    opt_immesc = TRUE;
    opt_ch_u = TRUE;
} else if (FIT_TYPE == "infec_profile" || FIT_TYPE == "ip_test") {
    extra_priors = list(
        v2_relu = "L 0.0 0.4 T 0.125 8",
        v2_serial = "L 0.0 0.4 T 0.01 1");
    opt_relu = TRUE;
    opt_serial = TRUE;
} else if (FIT_TYPE == "novoc") {
    opt_v2 = FALSE;
} else {
    stop("Need to specify fit type at command line.");
}

if (POP_SET == "else") {
    which_pops = c(1, 3, 9)
    pop_letter = "ELSE"
} else if (POP_SET == "mnsw") {
    which_pops = c(4, 5, 6, 10)
    pop_letter = "MNSW"
} else if (POP_SET == "all") {
    which_pops = c(1, 3, 4, 5, 6, 9, 10)
    pop_letter = "ALL"
} else {
    which_pops = as.numeric(POP_SET)
    pop_letter = paste0(POP_SET, "p")
    if (is.na(which_pops)) {
        stop("POP_SET must be else or all.");
    }
}

set_id = paste0(FIT_TYPE, "_", pop_letter);


# Adjust number of threads for fitting
N_THREADS = 54 + length(extra_priors) * 2;

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
                             date_end = date_fitting,
    dE  = cm_delay_gamma(2.5, 2.5, t_max = 15, t_step = 0.25)$p,
    dIp = cm_delay_gamma(2.5, 4.0, t_max = 15, t_step = 0.25)$p,
    dIs = cm_delay_gamma(2.5, 4.0, t_max = 15, t_step = 0.25)$p,
    dIa = cm_delay_gamma(5.0, 4.0, t_max = 15, t_step = 0.25)$p)
params = cm_split_matrices_ex_in(params, 15)

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
priorsI = list(
    tS = "U 0 60",
    u = "N 0.09 0.02 T 0.04 0.2",
    death_mean = "N 15 2 T 5 30",    # <<< co-cin
    hosp_admission = "N 8 1 T 4 20", # <<< co-cin
    icu_admission = "N 12.5 1 T 8 14", # <<< co-cin
    cfr_rlo = "N 0 0.1 T -2 2",
    cfr_rlo2 = "N 0 0.1 T -2 2",
    cfr_rlo3 = "N 0 0.1 T -2 2",
    hosp_rlo = "N 0 0.1 T -2 2", 
    icu_rlo = "N 0 0.1 T -2 2",
    icu_rlo2 = "N 0 0.1 T -2 2",
    contact_final = "N 1 0.1 T 0 1",
    contact_s0 = "E 0.1 0.1",
    contact_s1 = "E 0.1 0.1",
    disp_deaths = "E 10 10",
    disp_hosp_inc = "E 10 10",
    disp_hosp_prev = "E 10 10",
    disp_icu_prev = "E 10 10",
    concentration1 = "N 2 .3 T 2 10",
    concentration2 = "N 2 .2 T 2 10",
    concentration3 = "N 2 .1 T 2 10"
);

if (opt_v2) {
    priorsI = c(priorsI, list(
        v2_when = "U 144 365",
        v2_sgtf0 = "B 1.5 15",
        v2_disp = "E 10 10 T 0 0.7",
        v2_hosp_rlo = "N 0 0.1 T -4 4",
        v2_icu_rlo = "N 0 0.1 T -4 4",
        v2_cfr_rlo = "N 0 0.1 T -4 4"
    ));
}

priorsI = c(priorsI, extra_priors);


posteriorsI = list()
dynamicsI = list()
parametersI = list()

# Remove problematic virus entries
virus = virus[omega > 1e-9]

init_previous = TRUE
init_previous_amount = 1

existing_file = paste0("./fits/", set_id, REP_START - 1, ".qs");

if (file.exists(existing_file)) {
    saved = qread(existing_file)
    posteriorsI = saved[[1]]
    parametersI = saved[[2]]
    rm(saved)
}

for (i in seq_along(posteriorsI)) {
    if (!is.null(posteriorsI[[i]]) && "v2_conc" %in% names(posteriorsI[[i]])) {
        posteriorsI[[i]][, v2_disp := 1 / sqrt(v2_conc)];
    }
}

for (replic in REP_START:REP_END)
{
    # RCB checking execution time to test multithreading
    time1 <- Sys.time()
    
    # Loop through regions
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
        paramsI$schedule[[2]] = rlang::duplicate(paramsI$schedule[[1]]);
        for (i in seq_along(paramsI$schedule[[2]]$values)) {
            paramsI$schedule[[2]]$values[[i]][1] = paramsI$schedule[[1]]$values[[i]][1] +  0.2497655 / 100;
            paramsI$schedule[[2]]$values[[i]][2] = paramsI$schedule[[1]]$values[[i]][2] + -0.2307939 / 100;
            paramsI$schedule[[2]]$values[[i]][3] = paramsI$schedule[[1]]$values[[i]][3] + -1.5907698 / 100;
            paramsI$schedule[[2]]$values[[i]][4] = paramsI$schedule[[1]]$values[[i]][4] + -3.4866544 / 100;
            paramsI$schedule[[2]]$values[[i]][5] = paramsI$schedule[[1]]$values[[i]][5] + -3.4524518 / 100;
        }
        paramsI$schedule[[2]]$mode = "bypass";
    
        # contact placeholder for tier 3
        paramsI$schedule[[3]] = rlang::duplicate(paramsI$schedule[[1]]);
        for (i in seq_along(paramsI$schedule[[3]]$values)) {
            paramsI$schedule[[3]]$values[[i]][1] = paramsI$schedule[[1]]$values[[i]][1] +  2.080457 / 100;
            paramsI$schedule[[3]]$values[[i]][2] = paramsI$schedule[[1]]$values[[i]][2] + -8.045226 / 100;
            paramsI$schedule[[3]]$values[[i]][3] = paramsI$schedule[[1]]$values[[i]][3] + -2.476266 / 100;
            paramsI$schedule[[3]]$values[[i]][4] = paramsI$schedule[[1]]$values[[i]][4] + -10.144043 / 100;
            paramsI$schedule[[3]]$values[[i]][5] = paramsI$schedule[[1]]$values[[i]][5] + -7.681244 / 100;
        }
        paramsI$schedule[[3]]$mode = "bypass";
    
        # contact multiplier for gradual contact change
        paramsI$schedule[[4]] = list(
            parameter = "contact",
            pops = 0,
            mode = "multiply",
            values = rep(list(rep(1, 8)), 366),
            times = 0:365
        )
    
        # contact multiplier for september boost
        paramsI$schedule[[5]] = list(
            parameter = "contact",
            pops = 0,
            mode = "multiply",
            values = list(rep(1, 8)),
            times = c(244)
        )
        
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
    
        # load user defined functions
        cm_source_backend(
            user_defined = list(
                model_v2 = list(
                    cpp_changes = cpp_chgI_voc(priorsI, seasonality = opt_seas, 
                        v2 = opt_v2, v2_relu = opt_relu, v2_latdur = opt_latdur, v2_serial = opt_serial, v2_infdur = opt_infdur, v2_immesc = opt_immesc, v2_ch_u = opt_ch_u),
                    cpp_loglikelihood = cpp_likI_voc(paramsI, ldI, sitrepsI, seroI, virusI, sgtfI, p, date_fitting, priorsI, death_cutoff = 0, use_sgtf = opt_v2),
                    cpp_observer = cpp_obsI_voc(concentration = opt_conc, v2 = opt_v2, P.death, P.critical, priorsI)
                )
            )
        )

        priorsI2 = rlang::duplicate(priorsI)
        if (init_previous) {
            for (k in seq_along(priorsI2)) {
                pname = names(priorsI2)[k];
                if (length(posteriorsI) >= p && pname %in% names(posteriorsI[[p]])) {
                    init_values = quantile(posteriorsI[[p]][trial > mean(trial), get(pname)], c(0.025, 0.975));
                    cat(paste0("Using 95% CI ", init_values[1], " - ", init_values[2], " for initial values of parameter ", pname, 
                        " with probability ", init_previous_amount, "\n"));
                    priorsI2[[pname]] = paste0(priorsI2[[pname]], " I ", init_values[1], " ", init_values[2], " ", init_previous_amount);
                    cat(paste0(priorsI2[[pname]], "\n"));
                } else {
                    cat(paste0("Could not find init values for parameter ", pname, "\n"));
                    cat(paste0(priorsI2[[pname]], "\n"));
                }
            }
        }
    
        postI = cm_backend_mcmc_test(cm_translate_parameters(paramsI), priorsI2,
            seed = 0, 
            burn_in = ifelse(replic == REP_END, BURN_IN_FINAL, BURN_IN), 
            iterations = ifelse(replic == REP_END, ITER_FINAL, ITER), 
            n_threads = N_THREADS, classic_gamma = F, do_migration = T);
        setDT(postI)
        posteriorsI[[p]] = postI
    
        parametersI[[p]] = rlang::duplicate(paramsI)
        qsave(rlang::duplicate(list(posteriorsI, parametersI)), paste0("./fits/", set_id, replic, "-progress.qs"))
    
        print(p)
    }
    
    # RCB timing check again
    time2 <- Sys.time()
    print(time2-time1)
    # 45 mins for England

    qsave(rlang::duplicate(list(posteriorsI, parametersI)), paste0("./fits/", set_id, replic, ".qs"))
    
    # Generate SPI-M output
    # Sample dynamics from fit
    dynamicsI = list()
    for (p in which_pops)  {
        cat(paste0("Sampling fit for population ", p, "...\n"))
        
        # Source backend
        cm_source_backend(
            user_defined = list(
                model_v2 = list(
                    cpp_changes = cpp_chgI_voc(priorsI, seasonality = opt_seas, 
                        v2 = opt_v2, v2_relu = opt_relu, v2_latdur = opt_latdur, v2_serial = opt_serial, v2_infdur = opt_infdur, v2_immesc = opt_immesc, v2_ch_u = opt_ch_u),
                    cpp_loglikelihood = "",
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
    }
    
    # Concatenate dynamics for SPI-M output
    test = rbindlist(dynamicsI, fill = TRUE)
    test[, population := nhs_regions[population]]

    # Visually inspect fit
    plot = check_fit(test, parametersI, ld, sitreps, virus, sero, nhs_regions[which_pops], death_cutoff = 0, "2020-12-31")
    plot = plot + geom_vline(aes(xintercept = ymd("2020-12-24")), size = 0.25, linetype = "42")
    ggsave(paste0("./output/fit_", set_id, replic, ".pdf"), plot, width = 40, height = 25, units = "cm", useDingbats = FALSE)
    
    # Posteriors
    post = rbindlist(posteriorsI, idcol = "population")
    post[, pop := nhs_regions[population]]
    melted = melt(post, id.vars = c(1:5, ncol(post)))
    plot = ggplot(melted) + geom_density(aes(x = value, colour = pop)) + facet_wrap(~variable, scales = "free")
    ggsave(paste0("./output/post_", set_id, replic, ".pdf"), plot, width = 20, height = 15, units = "cm", useDingbats = FALSE)

    # Fit to SGTF data
    if (opt_v2) {
        sgtf[, qlo := qbeta(0.025, sgtf + 1, other + 1)]
        sgtf[, qhi := qbeta(0.975, sgtf + 1, other + 1)]
        vmodel = test[, .(I1 = sum(test_o), I2 = sum(test2_o), sgtf0 = v2_sgtf0[1], conc = 1/(v2_disp[1]*v2_disp[1])), by = .(t, population, run)]
        vmodel[, p2 := I2 / (I1 + I2)]
        vmodel[is.nan(p2), p2 := 0]
        vmodel[, sgtf := (1 - p2) * sgtf0 + p2];
        vmodel[, alpha := sgtf * (conc - 2) + 1]
        vmodel[, beta := (1 - sgtf) * (conc - 2) + 1]
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
        ggsave(paste0("./output/sgtf_check_", set_id, replic, ".pdf"), plotS, width = 20, height = 6, units = "cm", useDingbats = FALSE)
    }
}
