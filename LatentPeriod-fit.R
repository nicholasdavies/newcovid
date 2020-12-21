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
REP_START = 10
REP_END = 10
BURN_IN = 1500
BURN_IN_FINAL = 2500
ITER = 500
which_pops = c(3, 1, 9)

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

all_data = qread(datapath("processed-data-2020-12-19.qs"))
ld = all_data[[1]]
sitreps = all_data[[2]]
virus = all_data[[3]][!Data.source %like% "7a|7b|6a|6b"]
sero = all_data[[4]]
variant = all_data[[5]][sample_date >= "2020-10-01"]

#
# FITTING
#

# NUMBER OF REGIONS TO FIT
N_REG = 12;

# Build parameters for NHS regions
params = cm_parameters_SEI3R(nhs_regions[1:N_REG], deterministic = T, 
                             date_start = "2020-01-01", 
                             date_end = as.character(max(ld$date, sitreps$date, virus$End.date, sero$End.date) + 1),
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
schedule_all = readRDS(datapath("schedule3-2020-12-18.rds"));
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
    v2_latent = "N 1 0.2 T 0.1 2",
    v2_hosp_rlo = "N 0 0.1 T -2 2", # hosp x[20]
    v2_cfr_rel = "N 1 0.1 T 0.1 4" # cfr_rel x[21]
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
}

for (replic in REP_START:REP_END)
{
    init_previous = TRUE
    init_previous_amount = 1
    
    # RCB checking execution time to test multithreading
    time1 <- Sys.time()
    
    # Loop through regions
    for (p in which_pops) { ### 1:12) {
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
    
        ldI = rlang::duplicate(ld);
        ldI = ldI[pid == p - 1];
        sitrepsI = rlang::duplicate(sitreps);
        sitrepsI = sitrepsI[pid == p - 1];
        seroI = rlang::duplicate(sero);
        seroI = seroI[pid == p - 1 & Data.source != "NHSBT"];   # sero: all but NHSBT
        virusI = rlang::duplicate(virus);
        virusI = virusI[pid == p - 1 & Data.source %like% "REACT"]; # virus: REACT only
        variantI = copy(variant);
        variantI = variantI[pid == p - 1];
    
        # load user defined functions
        cm_source_backend(
            user_defined = list(
                model_v2 = list(
                    cpp_changes = cpp_chgI_LatentPeriod(),
                    cpp_loglikelihood = cpp_likI(paramsI, ldI, sitrepsI, seroI, virusI, variantI, p),
                    cpp_observer = cpp_obsI(P.death)
                )
            )
        )
        
        priorsI2 = rlang::duplicate(priorsI)
        if (init_previous) {
            for (k in seq_along(priorsI2)) {
                pname = names(priorsI2)[k];
                if (length(posteriorsI) >= p && pname %in% names(posteriorsI[[p]])) {
                    init_values = quantile(posteriorsI[[p]][[pname]], c(0.25, 0.75));
                    cat(paste0("Using IQR ", init_values[1], " - ", init_values[2], " for initial values of parameter ", pname, 
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
            seed = 0, burn_in = ifelse(replic == REP_END, BURN_IN_FINAL, BURN_IN), 
            iterations = ITER, n_threads = N_THREADS, classic_gamma = T);
        setDT(postI)
        posteriorsI[[p]] = postI
    
        # # Sampling fits
        # paramsI2 = rlang::duplicate(paramsI)
        # paramsI2$time1 = as.character(ymd(paramsI$time1) + 56);
        # test = cm_backend_sample_fit_test(cm_translate_parameters(paramsI2), postI, 500, seed = 0);
        # 
        # test = rbindlist(test)
        # test = test[, population := p]
        # dynamicsI[[p]] = test
    
        parametersI[[p]] = rlang::duplicate(paramsI)
        qsave(rlang::duplicate(list(posteriorsI, parametersI)), paste0("./fits/LatentPeriod-pp", replic, "-progress.qs"))
    
        print(p)
    }
    
    # RCB timing check again
    time2 <- Sys.time()
    print(time2-time1)
    # 45 mins for England

    qsave(rlang::duplicate(list(posteriorsI, parametersI)), paste0("./fits/LatentPeriod-pp", replic, ".qs"))
    
    # Generate SPI-M output
    # Sample dynamics from fit
    dynamicsI = list()
    dynamics0 = list()
    for (p in which_pops)  { ### 1:12
        cat(paste0("Sampling fit for population ", p, "...\n"))
        
        # Source backend
        cm_source_backend(
            user_defined = list(
                model_v2 = list(
                    cpp_changes = cpp_chgI_LatentPeriod(),
                    cpp_loglikelihood = "",
                    cpp_observer = cpp_obsI(P.death)
                )
            )
        )
        
        # Sampling fits
        paramsI2 = rlang::duplicate(parametersI[[p]])
        paramsI2$time1 = as.character(ymd(parametersI[[p]]$time1) + 56);
        test = cm_backend_sample_fit_test(cm_translate_parameters(paramsI2), posteriorsI[[p]], 100, seed = 0);
        
        test = rbindlist(test)
        test[, population := p]
        dynamicsI[[p]] = test
        
        # Again, without new variant
        posteriors0 = copy(posteriorsI[[p]]);
        posteriors0[, v2_when := 9999];
        test0 = cm_backend_sample_fit_test(cm_translate_parameters(paramsI2), posteriors0, 100, seed = 0);
        
        test0 = rbindlist(test0)
        test0[, population := p]
        dynamics0[[p]] = test0
    }
    
    # Concatenate dynamics for SPI-M output
    test = rbindlist(dynamicsI, fill = TRUE)
    test[, population := nhs_regions[population]]
    test0 = rbindlist(dynamics0, fill = TRUE)
    test0[, population := nhs_regions[population]]
    
    # Fit to COG data
    variant[, qlo := qbeta(0.025, var2 + 1, all - var2 + 1)]
    variant[, qhi := qbeta(0.975, var2 + 1, all - var2 + 1)]
    vmodel = test[, .(p2 = sum(Ip2 + Is2 + Ia2) / sum(Ip + Is + Ia + Ip2 + Is2 + Ia2)), by = .(t, population, run)]
    vmodel[is.nan(p2), p2 := 0]
    vmodel = vmodel[, as.list(quantile(p2, c(0.025, 0.5, 0.975))), by = .(t, nhs_name = population)]
    plot = ggplot(variant[(pid + 1) %in% which_pops]) +
        geom_ribbon(aes(x = sample_date, ymin = qlo, ymax = qhi), fill = "black", alpha = 0.1) +
        geom_ribbon(data = vmodel[t + ymd("2020-01-01") >= "2020-10-01"], 
            aes(x = ymd("2020-01-01") + t, ymin = `2.5%`, ymax = `97.5%`), fill = "darkorchid", alpha = 0.5) +
        geom_line(aes(x = sample_date, y = var2 / all), size = 0.25) +
        facet_wrap(~nhs_name) +
        labs(x = NULL, y = "Frequency of\nSARS-CoV-2 variant") +
        scale_x_date(date_breaks = "1 month", date_labels = "%b")
    ggsave(paste0("./output/LatentPeriod-variant_check_", replic, ".pdf"), plot, width = 20, height = 6, units = "cm", useDingbats = FALSE)
    
    # Posteriors of interest
    post = rbindlist(posteriorsI[which_pops], idcol = "population")
    post[, population := nhs_regions[which_pops[population]]]
    post = melt(post, id.vars = 1, measure.vars = c("v2_immesc", "v2_hosp_rlo", "v2_cfr_rel"))
    post[variable == "v2_latent", variable := "Relative length of\nlatent period for B.1.1.7 N501Y"]
    post[variable == "v2_hosp_rlo", value := exp(value)]
    post[variable == "v2_hosp_rlo", variable := "Associated OR of hospitalisation"]
    post[variable == "v2_cfr_rel", variable := "Associated RR of death"]
    
    prior = data.table(x = rep(seq(0.7, 1.3, 0.03), 3), 
        variable = rep(c("Relative length of\nlatent period for B.1.1.7 N501Y", "Associated OR of hospitalisation", "Associated RR of death"), each = 21))
    prior[variable == "Relative length of\nlatent period for B.1.1.7 N501Y", y := dnorm(x, 1, 0.2)]
    prior[variable == "Associated OR of hospitalisation", y := dlnorm(x, 0, 0.1)]
    prior[variable == "Associated RR of death", y := dnorm(x, 1, 0.1)]
    plot = ggplot(post) +
        geom_line(data = prior, aes(x, y), colour = "#888888") +
        geom_density(aes(value, colour = population), adjust = 4) +
        #geom_histogram(aes(value, colour = population), bins = 20) +
        facet_wrap(~variable, scales = "free") +
        theme(legend.position = c(0.01, 0.9)) +
        labs(x = NULL, y = NULL, colour = NULL) +
        expand_limits(x = 1)
    ggsave(paste0("./output/LatentPeriod-variant_stats_", replic, ".pdf"), plot, width = 20, height = 6, units = "cm", useDingbats = FALSE)

    # Visually inspect fit
    plot_a = check_fit(test0, ld, sitreps, virus, sero, nhs_regions[which_pops])
    plot_b = check_fit(test, ld, sitreps, virus, sero, nhs_regions[which_pops])
    plot = cowplot::plot_grid(plot_a, plot_b, nrow = 1, labels = letters)
    ggsave(paste0("./output/LatentPeriod-fit_", replic, ".pdf"), plot, width = 30, height = 25, units = "cm", useDingbats = FALSE)

    # # Generate SPI-M output
    # ## UPDATE (creation year, creation month, creation day, forecast start date, forecast end date)
    # mtp_output = SPIM_output_full(test[!population %in% c("England", "United Kingdom")], 2020, 12, 15, "2020-12-15", "2021-01-28")
    # 
    # # Inspect output
    # plot = ggplot(mtp_output[AgeBand == "All"], aes(x = make_date(`Year of Value`, `Month of Value`, `Day of Value`))) +
    #     geom_ribbon(aes(ymin = `Quantile 0.05`, ymax = `Quantile 0.95`, fill = `Geography`)) +
    #     facet_grid(ValueType ~ Geography, scales = "free") +
    #     theme(legend.position = "none")
    # 
    # ggsave(paste0("./output/mtp_check_spim_", replic, ".pdf"), plot, width = 25, height = 20, units = "cm", useDingbats = FALSE)
    # 
    # 
    # # Save output
    # fwrite(mtp_output, paste0("./output/SPIM_mtp_", replic, ".csv"))
}
