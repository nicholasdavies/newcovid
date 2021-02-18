# project-voc: For VOC projections.

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

which_pops = c(1, 3, 4, 5, 6, 9, 10)

uk_covid_data_path = "./fitting_data/";
datapath = function(x) paste0(uk_covid_data_path, x)

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

project = function(FIT_TYPE, FIT_FILE, SCH_FILE, VAX_TYPE, forced_seasonality = 0, include_voc = TRUE)
{
    data_file = "processed-data-2021-01-08.qs"
    mobility_file = SCH_FILE
    project_date = "2021-06-30"
    
    output_id = paste0(FIT_TYPE, "_", str_remove_all(SCH_FILE, "\\.rds"), "_", VAX_TYPE);
    
    # Model options
    opt_conc = TRUE;
    opt_seas = FALSE;
    
    opt_relu = FALSE;
    opt_latdur = FALSE;
    opt_serial = FALSE;
    opt_infdur = FALSE;
    opt_immesc = FALSE;
    opt_ch_u = FALSE;
    
    if (FIT_TYPE == "relu") {
        extra_priors = list(v2_relu = "L 0.0 0.4 T 0.25 4");
        opt_relu = TRUE;
    } else if (FIT_TYPE == "latdur") {
        extra_priors = list(v2_latdur = "L 0.0 0.4 T 0.01 8");
        opt_latdur = TRUE;
    } else if (FIT_TYPE == "serial") {
        extra_priors = list(v2_serial = "L 0.0 0.4 T 0.01 8");
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
    } else if (FIT_TYPE == "infec_profile") {
        extra_priors = list(
            v2_relu = "L 0.0 0.4 T 0.125 8", 
            v2_latdur = "L 0.0 0.4 T 0.01 8", 
            v2_infdur = "L 0.0 0.4 T 0.01 8");
        opt_relu = TRUE;
        opt_latdur = TRUE;
        opt_infdur = TRUE;
    } else {
        stop("Need to specify fit type at command line.");
    }
    set_id = FIT_TYPE;
    
    if (VAX_TYPE == "0") {
        print("No vaccination.");
    } else if (VAX_TYPE == "v200k") {
        vax_rate = 200000;
        vax_ei_v = rep(0.6, 16);
        vax_ed_vi = rep(0.875, 16);
        vax_targeting = targeting_old_to_young;
    } else if (VAX_TYPE == "v2M") {
        vax_rate = 2000000;
        vax_ei_v = rep(0.6, 16);
        vax_ed_vi = rep(0.875, 16);
        vax_targeting = targeting_old_to_young;
    } else {
        stop("Unsupported VAX_TYPE.");
    }
    
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
    # Individual fits
    #
    
    
    # Fitting
    priorsI = c(list(
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
        concentration3 = "N 2 .1 T 2 10",
        
        v2_when = "U 144 365",
        v2_sgtf0 = "B 1.5 15",
        v2_conc = "E 0.1 0.1 T 2 1000",
        v2_hosp_rlo = "N 0 0.1 T -4 4",
        v2_icu_rlo = "N 0 0.1 T -4 4",
        v2_cfr_rlo = "N 0 0.1 T -4 4"
    ), extra_priors)
    
    
    posteriorsI = list()
    dynamicsI = list()
    parametersI = list()
    
    # Remove problematic virus entries
    virus = virus[omega > 1e-9]
    
    init_previous = TRUE
    init_previous_amount = 1
    
    saved = qread(FIT_FILE)
    posteriorsI = saved[[1]]
    parametersI = saved[[2]]
    rm(saved)
    
    N_REG = 12
    
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
    
    # Loop through regions to reassign schedule
    for (p in which_pops) {
        paramsI = rlang::duplicate(parametersI[[p]]);
        paramsI$schedule = list();
        j = 1;
        for (i in seq_along(schedule)) {
            if (p - 1 == schedule[[i]]$pops) {
                paramsI$schedule[[j]] = rlang::duplicate(schedule[[i]]);
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
        
        parametersI[[p]] = rlang::duplicate(paramsI)
        print(p)
    }
    
    # Assign vaccination
    if (VAX_TYPE == "0") {
        VS = NULL;
    } else {
        cat("Building vaccine schedule...\n");
        VS = make_vaccine_schedule(parametersI, "2021-01-01", rep(vax_rate, 365), vax_targeting, which_pops)
        
        for (p in which_pops)
        {
            parametersI[[p]]$pop[[1]]$ei_v = vax_ei_v;
            parametersI[[p]]$pop[[1]]$ed_vi = vax_ed_vi;
        }
    }
    
    # Generate SPI-M output
    # Sample dynamics from fit
    dynamicsI = list()
    for (p in which_pops)  {
        cat(paste0("Sampling fit for population ", p, "...\n"))
        
        # Make parameters
        paramsI2 = rlang::duplicate(parametersI[[p]])
        paramsI2$time1 = project_date;
        
        # Source backend
        cm_source_backend(
            user_defined = list(
                model_v2 = list(
                    cpp_changes = cpp_chgI_voc(priorsI, seasonality = opt_seas, 
                        v2 = include_voc, v2_relu = opt_relu, v2_latdur = opt_latdur, v2_serial = opt_serial, v2_infdur = opt_infdur, v2_immesc = opt_immesc, v2_ch_u = opt_ch_u),
                    cpp_loglikelihood = "",
                    cpp_observer = c(
                        cpp_obsI_voc(concentration = opt_conc, v2 = include_voc, P.death, P.critical, priorsI),
                        if (!is.null(VS)) cpp_obsI_vax(paramsI2, VS[[p]]) else NULL,
                        if (forced_seasonality > 0) cpp_obsI_seasonality(forced_seasonality, 366) else NULL
                    )
                )
            )
        )
        
        # Sampling fits
        test = cm_backend_sample_fit_test(cm_translate_parameters(paramsI2), posteriorsI[[p]], 100, seed = 0);
        rows = cm_backend_sample_fit_rows(cm_translate_parameters(paramsI2), posteriorsI[[p]], 100, seed = 0);
        
        test = rbindlist(test)
        test[, population := p]
        
        # Add dispersion parameters
        disp = posteriorsI[[p]][rows, .SD, .SDcols = patterns("^disp|v2_conc|v2_sgtf0")]
        disp[, run := .I]
        test = merge(test, disp, by = "run")
    
        dynamicsI[[p]] = test
    }
    
    # Concatenate dynamics for SPI-M output
    test = rbindlist(dynamicsI, fill = TRUE)
    test[, population := nhs_regions[population]]
    
    return (test)
}


# Determine population sizes
popsize = NULL
parametersI = qread("./fits/relu_ALL18.qs")[[2]]
for (i in seq_along(parametersI)) {
    if (!is.null(parametersI[[i]])) {
        popsize = rbind(popsize,
            data.table(Geography = parametersI[[i]]$pop[[1]]$name, population_size = sum(parametersI[[i]]$pop[[1]]$size))
        )
    }
}

# Load data
data_file = "processed-data-2021-01-08.qs"

all_data = qread(datapath(data_file))
ld = all_data[[1]]
sitreps = all_data[[2]]
virus = all_data[[3]][!Data.source %like% "7a|7b|6a|6b"]
sero = all_data[[4]]
sgtf = all_data[[5]]

max_beds_w1_england = sitreps[name == "England" & date <= "2020-09-01", max(n_in_all_beds, na.rm = T)]
max_deaths_w1_england = ld[name == "England" & date <= "2020-09-01", max(N, na.rm = T)]


for (VAX in c("0", "v200k", "v2M")) {

    p1 = project("relu", "./fits/relu_ALL18.qs", "scenario1.rds", VAX)
    p2 = project("relu", "./fits/relu_ALL18.qs", "scenario2.rds", VAX)
    p3 = project("relu", "./fits/relu_ALL18.qs", "scenario3.rds", VAX)
    p5 = project("relu", "./fits/relu_ALL18.qs", "scenario5.rds", VAX)
    
    plot_england = plot_projection(list(p1, p2, p3, p5), 
        c("Moderate (October 2020)", "High (November 2020) with schools open", "High with schools closed", "Very high (March 2020)"), 
        "2020-12-01", england = TRUE, reduced = TRUE, pal = "Dark2", 
        hosp_line = max_beds_w1_england / 1000, deaths_line = max_deaths_w1_england, colour_label = "Stringency of NPIs")
    
    plot_regions = plot_projection(list(p1, p2, p3, p5), 
        c("Moderate (October 2020)", "High (November 2020) with schools open", "High with schools closed", "Very high (March 2020)"), 
        "2020-12-01", england = FALSE, reduced = FALSE, pal = "Dark2")
    
    tb_1 = summarize_projection(p1, "2020-12-15", popsize)
    tb_2 = summarize_projection(p2, "2020-12-15", popsize)
    tb_3 = summarize_projection(p3, "2020-12-15", popsize)
    tb_5 = summarize_projection(p5, "2020-12-15", popsize)
    
    tbEngland = england_only(list(tb_1, tb_2, tb_3, tb_5), 
        c("Moderate (October 2020)", "High (November 2020) with schools open", "High with schools closed", "Very high (March 2020)"))
    
    plot_england$plot_env$proj_list = NULL
    plot_regions$plot_env$proj_list = NULL
    
    qsave(plot_england, paste0("./output/plot_england_", VAX, ".qs"))
    qsave(plot_regions, paste0("./output/plot_regions_", VAX, ".qs"))
    fwrite(tb_1, paste0("./output/table_1_", VAX, ".csv"));
    fwrite(tb_2, paste0("./output/table_2_", VAX, ".csv"));
    fwrite(tb_3, paste0("./output/table_3_", VAX, ".csv"));
    fwrite(tb_5, paste0("./output/table_5_", VAX, ".csv"));
    fwrite(tbEngland, paste0("./output/table_england_", VAX, ".csv"));
}


# Assemble plots
pe_0 = qread("./output/plot_england_0.qs")
pe_v200k = qread("./output/plot_england_v200k.qs")
pe_v2M = qread("./output/plot_england_v2M.qs")

data_blank = data.table(variable = 
        factor(c("Hospital beds\noccupied (thousands)", "Deaths"), 
            levels = c("Rt", "Hospital beds\noccupied (thousands)", "Deaths")), x = ymd("2021-01-01"), y = c(145, 4000))

plot_e = cowplot::plot_grid(
    pe_0 + geom_blank(data = data_blank, aes(x, y)) + 
        theme(legend.position = "none", strip.text.x = element_blank(), plot.title = element_text(size = 10)) + 
        labs(title = "No vaccination"),
    pe_v200k + geom_blank(data = data_blank, aes(x, y)) + 
        theme(legend.spacing.x = unit(0.5, "cm"), strip.text.x = element_blank(), plot.title = element_text(size = 10)) + 
        labs(title = "200,000 vaccinations per week"),
    pe_v2M + geom_blank(data = data_blank, aes(x, y)) + 
        theme(legend.position = "none", strip.text.x = element_blank(), plot.title = element_text(size = 10)) + 
        labs(title = "2 million vaccinations per week"),
    nrow = 1, labels = LETTERS, label_size = 10, align = "hv", axis = "bottom")

ggsave("./output/projections_england.pdf", plot_e, width = 30, height = 14, units = "cm", useDingbats = FALSE)
ggsave("./output/projections_england.png", plot_e, width = 30, height = 14, units = "cm")

pr_0 = qread("./output/plot_regions_0.qs")
pr_v200k = qread("./output/plot_regions_v200k.qs")
pr_v2M = qread("./output/plot_regions_v2M.qs")

ggsave("./output/projections_regions_0.pdf", pr_0, width = 38, height = 22, units = "cm", useDingbats = FALSE)
ggsave("./output/projections_regions_0.png", pr_0, width = 38, height = 22, units = "cm")
ggsave("./output/projections_regions_v200k.pdf", pr_v200k, width = 38, height = 22, units = "cm", useDingbats = FALSE)
ggsave("./output/projections_regions_v200k.png", pr_v200k, width = 38, height = 22, units = "cm")
ggsave("./output/projections_regions_v2M.pdf", pr_v2M, width = 38, height = 22, units = "cm", useDingbats = FALSE)
ggsave("./output/projections_regions_v2M.png", pr_v2M, width = 38, height = 22, units = "cm")




###############
# Seasonality #
###############


for (VAX in c("0", "v200k", "v2M")) {

    p1 = project("relu", "./fits/relu_ALL18.qs", "scenario1.rds", VAX, 0.1)
    p2 = project("relu", "./fits/relu_ALL18.qs", "scenario2.rds", VAX, 0.1)
    p3 = project("relu", "./fits/relu_ALL18.qs", "scenario3.rds", VAX, 0.1)
    p5 = project("relu", "./fits/relu_ALL18.qs", "scenario5.rds", VAX, 0.1)
    
    plot_england = plot_projection(list(p1, p2, p3, p5), 
        c("Moderate (October 2020)", "High (November 2020) with schools open", "High with schools closed", "Very high (March 2020)"), 
        "2020-12-01", england = TRUE, reduced = TRUE, pal = "Dark2", 
        hosp_line = max_beds_w1_england / 1000, deaths_line = max_deaths_w1_england, colour_label = "Stringency of NPIs")
    
    plot_regions = plot_projection(list(p1, p2, p3, p5), 
        c("Moderate (October 2020)", "High (November 2020) with schools open", "High with schools closed", "Very high (March 2020)"), 
        "2020-12-01", england = TRUE, reduced = TRUE, pal = "Dark2")
    
    tb_1 = summarize_projection(p1, "2020-12-15", popsize)
    tb_2 = summarize_projection(p2, "2020-12-15", popsize)
    tb_3 = summarize_projection(p3, "2020-12-15", popsize)
    tb_5 = summarize_projection(p5, "2020-12-15", popsize)
    
    tbEngland = england_only(list(tb_1, tb_2, tb_3, tb_5), 
        c("Moderate (October 2020)", "High (November 2020) with schools open", "High with schools closed", "Very high (March 2020)"))
    
    plot_england$plot_env$proj_list = NULL
    plot_regions$plot_env$proj_list = NULL
    
    qsave(plot_england, paste0("./output/seas_plot_england_", VAX, ".qs"))
    qsave(plot_regions, paste0("./output/seas_plot_regions_", VAX, ".qs"))
    fwrite(tb_1, paste0("./output/seas_table_1_", VAX, ".csv"));
    fwrite(tb_2, paste0("./output/seas_table_2_", VAX, ".csv"));
    fwrite(tb_3, paste0("./output/seas_table_3_", VAX, ".csv"));
    fwrite(tb_5, paste0("./output/seas_table_5_", VAX, ".csv"));
    fwrite(tbEngland, paste0("./output/seas_table_england_", VAX, ".csv"));
}


# Assemble plots
pe_0 = qread("./output/seas_plot_england_0.qs")
pe_v200k = qread("./output/seas_plot_england_v200k.qs")
pe_v2M = qread("./output/seas_plot_england_v2M.qs")

data_blank = data.table(variable = 
        factor(c("Hospital beds\noccupied (thousands)", "Deaths"), 
            levels = c("Rt", "Hospital beds\noccupied (thousands)", "Deaths")), x = ymd("2021-01-01"), y = c(145, 4000))

plot_e = cowplot::plot_grid(
    pe_0 + geom_blank(data = data_blank, aes(x, y)) + 
        theme(legend.position = "none", strip.text.x = element_blank(), plot.title = element_text(size = 10)) + 
        labs(title = "No vaccination"),
    pe_v200k + geom_blank(data = data_blank, aes(x, y)) + 
        theme(legend.spacing.x = unit(0.5, "cm"), strip.text.x = element_blank(), plot.title = element_text(size = 10)) + 
        labs(title = "200,000 vaccinations per week"),
    pe_v2M + geom_blank(data = data_blank, aes(x, y)) + 
        theme(legend.position = "none", strip.text.x = element_blank(), plot.title = element_text(size = 10)) + 
        labs(title = "2 million vaccinations per week"),
    nrow = 1, labels = LETTERS, label_size = 10, align = "hv", axis = "bottom")

ggsave("./output/seas_projections_england.pdf", plot_e, width = 30, height = 14, units = "cm", useDingbats = FALSE)
ggsave("./output/seas_projections_england.png", plot_e, width = 30, height = 14, units = "cm")

pr_0 = qread("./output/seas_plot_regions_0.qs")
pr_v200k = qread("./output/seas_plot_regions_v200k.qs")
pr_v2M = qread("./output/seas_plot_regions_v2M.qs")

ggsave("./output/seas_projections_regions_0.pdf", pr_0, width = 30, height = 12, units = "cm", useDingbats = FALSE)
ggsave("./output/seas_projections_regions_0.png", pr_0, width = 30, height = 12, units = "cm")
ggsave("./output/seas_projections_regions_v200k.pdf", pr_v200k, width = 30, height = 12, units = "cm", useDingbats = FALSE)
ggsave("./output/seas_projections_regions_v200k.png", pr_v200k, width = 30, height = 12, units = "cm")
ggsave("./output/seas_projections_regions_v2M.pdf", pr_v2M, width = 30, height = 12, units = "cm", useDingbats = FALSE)
ggsave("./output/seas_projections_regions_v2M.png", pr_v2M, width = 30, height = 12, units = "cm")




#####################
# Infectious period #
#####################

which_pops = c(1, 3, 9)

plot_projection(list(p1), 
        c("Moderate (October 2020)"), 
        "2020-12-01", england = FALSE, reduced = FALSE, pal = "Dark2")

for (VAX in c("0", "v200k", "v2M")) {

    p1 = project("infdur", "./fits/infdur_ELSE11.qs", "scenario1.rds", VAX)
    p2 = project("infdur", "./fits/infdur_ELSE11.qs", "scenario2.rds", VAX)
    p3 = project("infdur", "./fits/infdur_ELSE11.qs", "scenario3.rds", VAX)
    p5 = project("infdur", "./fits/infdur_ELSE11.qs", "scenario5.rds", VAX)
    
    plot_regions = plot_projection(list(p1, p2, p3, p5), 
        c("Moderate (October 2020)", "High (November 2020) with schools open", "High with schools closed", "Very high (March 2020)"), 
        "2020-12-01", england = FALSE, reduced = FALSE, pal = "Dark2")
    
    tb_1 = summarize_projection(p1, "2020-12-15", popsize)
    tb_2 = summarize_projection(p2, "2020-12-15", popsize)
    tb_3 = summarize_projection(p3, "2020-12-15", popsize)
    tb_5 = summarize_projection(p5, "2020-12-15", popsize)
    
    plot_regions$plot_env$proj_list = NULL
    
    qsave(plot_regions, paste0("./output/infdur_plot_regions_", VAX, ".qs"))
    fwrite(tb_1, paste0("./output/infdur_table_1_", VAX, ".csv"));
    fwrite(tb_2, paste0("./output/infdur_table_2_", VAX, ".csv"));
    fwrite(tb_3, paste0("./output/infdur_table_3_", VAX, ".csv"));
    fwrite(tb_5, paste0("./output/infdur_table_5_", VAX, ".csv"));
}


# Assemble plots
pr_0 = qread("./output/infdur_plot_regions_0.qs")
pr_v200k = qread("./output/infdur_plot_regions_v200k.qs")
pr_v2M = qread("./output/infdur_plot_regions_v2M.qs")

ggsave("./output/infdur_projections_regions_0.pdf", pr_0, width = 38, height = 22, units = "cm", useDingbats = FALSE)
ggsave("./output/infdur_projections_regions_0.png", pr_0, width = 38, height = 22, units = "cm")
ggsave("./output/infdur_projections_regions_v200k.pdf", pr_v200k, width = 38, height = 22, units = "cm", useDingbats = FALSE)
ggsave("./output/infdur_projections_regions_v200k.png", pr_v200k, width = 38, height = 22, units = "cm")
ggsave("./output/infdur_projections_regions_v2M.pdf", pr_v2M, width = 38, height = 22, units = "cm", useDingbats = FALSE)
ggsave("./output/infdur_projections_regions_v2M.png", pr_v2M, width = 38, height = 22, units = "cm")




###############
# Without VOC #
###############

which_pops = england_pops

for (VAX in c("0", "v200k", "v2M")) {

    p1 = project("relu", "./fits/relu_ALL18.qs", "scenario1.rds", VAX, 0.1, include_voc = FALSE)
    p2 = project("relu", "./fits/relu_ALL18.qs", "scenario2.rds", VAX, 0.1, include_voc = FALSE)
    p3 = project("relu", "./fits/relu_ALL18.qs", "scenario3.rds", VAX, 0.1, include_voc = FALSE)
    p5 = project("relu", "./fits/relu_ALL18.qs", "scenario5.rds", VAX, 0.1, include_voc = FALSE)
    
    plot_england = plot_projection(list(p1, p2, p3, p5), 
        c("Moderate (October 2020)", "High (November 2020) with schools open", "High with schools closed", "Very high (March 2020)"), 
        "2020-12-01", england = TRUE, reduced = TRUE, pal = "Dark2", 
        hosp_line = max_beds_w1_england / 1000, deaths_line = max_deaths_w1_england, colour_label = "Stringency of NPIs")
    
    plot_regions = plot_projection(list(p1, p2, p3, p5), 
        c("Moderate (October 2020)", "High (November 2020) with schools open", "High with schools closed", "Very high (March 2020)"), 
        "2020-12-01", england = TRUE, reduced = TRUE, pal = "Dark2")
    
    tb_1 = summarize_projection(p1, "2020-12-15", popsize)
    tb_2 = summarize_projection(p2, "2020-12-15", popsize)
    tb_3 = summarize_projection(p3, "2020-12-15", popsize)
    tb_5 = summarize_projection(p5, "2020-12-15", popsize)
    
    tbEngland = england_only(list(tb_1, tb_2, tb_3, tb_5), 
        c("Moderate (October 2020)", "High (November 2020) with schools open", "High with schools closed", "Very high (March 2020)"))
    
    plot_england$plot_env$proj_list = NULL
    plot_regions$plot_env$proj_list = NULL
    
    qsave(plot_england, paste0("./output/0voc_plot_england_", VAX, ".qs"))
    qsave(plot_regions, paste0("./output/0voc_plot_regions_", VAX, ".qs"))
    fwrite(tb_1, paste0("./output/0voc_table_1_", VAX, ".csv"));
    fwrite(tb_2, paste0("./output/0voc_table_2_", VAX, ".csv"));
    fwrite(tb_3, paste0("./output/0voc_table_3_", VAX, ".csv"));
    fwrite(tb_5, paste0("./output/0voc_table_5_", VAX, ".csv"));
    fwrite(tbEngland, paste0("./output/0voc_table_england_", VAX, ".csv"));
}


# Assemble plots
pe_0 = qread("./output/0voc_plot_england_0.qs")
pe_v200k = qread("./output/0voc_plot_england_v200k.qs")
pe_v2M = qread("./output/0voc_plot_england_v2M.qs")

data_blank = data.table(variable = 
        factor(c("Hospital beds\noccupied (thousands)", "Deaths"), 
            levels = c("Rt", "Hospital beds\noccupied (thousands)", "Deaths")), x = ymd("2021-01-01"), y = c(145, 4000))

plot_e = cowplot::plot_grid(
    pe_0 + geom_blank(data = data_blank, aes(x, y)) + 
        theme(legend.position = "none", strip.text.x = element_blank(), plot.title = element_text(size = 10)) + 
        labs(title = "No vaccination"),
    pe_v200k + geom_blank(data = data_blank, aes(x, y)) + 
        theme(legend.spacing.x = unit(0.5, "cm"), strip.text.x = element_blank(), plot.title = element_text(size = 10)) + 
        labs(title = "200,000 vaccinations per week"),
    pe_v2M + geom_blank(data = data_blank, aes(x, y)) + 
        theme(legend.position = "none", strip.text.x = element_blank(), plot.title = element_text(size = 10)) + 
        labs(title = "2 million vaccinations per week"),
    nrow = 1, labels = LETTERS, label_size = 10, align = "hv", axis = "bottom")

ggsave("./output/0voc_projections_england.pdf", plot_e, width = 30, height = 14, units = "cm", useDingbats = FALSE)
ggsave("./output/0voc_projections_england.png", plot_e, width = 30, height = 14, units = "cm")

pr_0 = qread("./output/0voc_plot_regions_0.qs")
pr_v200k = qread("./output/0voc_plot_regions_v200k.qs")
pr_v2M = qread("./output/0voc_plot_regions_v2M.qs")

ggsave("./output/0voc_projections_regions_0.pdf", pr_0, width = 30, height = 12, units = "cm", useDingbats = FALSE)
ggsave("./output/0voc_projections_regions_0.png", pr_0, width = 30, height = 12, units = "cm")
ggsave("./output/0voc_projections_regions_v200k.pdf", pr_v200k, width = 30, height = 12, units = "cm", useDingbats = FALSE)
ggsave("./output/0voc_projections_regions_v200k.png", pr_v200k, width = 30, height = 12, units = "cm")
ggsave("./output/0voc_projections_regions_v2M.pdf", pr_v2M, width = 30, height = 12, units = "cm", useDingbats = FALSE)
ggsave("./output/0voc_projections_regions_v2M.png", pr_v2M, width = 30, height = 12, units = "cm")
