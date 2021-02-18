
plot_fits = function(fit_filename, FIT_TYPE, which_pops, FLAG = "")
{
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
    
    fit = qread(fit_filename)
    posteriorsI = fit[[1]]
    parametersI = fit[[2]]

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
    ggsave(paste0("./output/check_fit_", FIT_TYPE, FLAG, ".pdf"), plot, width = 20, height = 25, units = "cm", useDingbats = FALSE)
}

plot_fits("./fits/relu_ALL18.qs", "relu", england_pops)
plot_fits("./fits/relu_ALL18.qs", "novoc", england_pops, "_removed")
plot_fits("./fits/novoc_ELSE4.qs", "novoc", c(1, 3, 9), "_never_present")


