library(data.table)

SPIM_output_1var = function(summ, varname, outname, cyear, cmonth, cday, ymd_from, ymd_to, size)
{
    quants = seq(0.05, 0.95, by = 0.05)

    # new code...
    rq = data.table(run = 1:summ[, max(run)], q = runif(summ[, max(run)]))
    summ = merge(summ, rq, by = "run")
    if (is.character(size)) {
        summ[, rv := qnbinom(q, size = 1.0 / get(size)^2, mu = get(varname))]
    } else if (size == 0) { # poisson
        summ[, rv := qpois(q, lambda = get(varname))]
    } else if (size == -1) { # no variation
        summ[, rv := get(varname)]
    } else if (size > 0) { # neg binom
        summ[, rv := qnbinom(q, size = size, mu = get(varname))]
    } else {
        stop("Size must be -1, 0, or a positive number")
    }
    qsumm = summ[, as.list(quantile(rv, quants)), by = .(t, population, age_group)]

    #qsumm = summ[, as.list(qnbinom(quants, size = 20, mu = quantile(get(varname), quants, na.rm = T))), by = .(t, population, age_group)]
    qsumm[, date := t + ymd("2020-01-01")]
    qsumm = rbind(qsumm,
        qsumm[, lapply(.SD, sum), .SDcols = `5%`:`95%`, by = .(t, population, date, age_group = rep("All", nrow(qsumm)))],
        use.names = TRUE, fill = TRUE)

    data_out = qsumm[date %between% c(ymd_from, ymd_to),
        .(Group = "LSHTM", Model = "Transmission", Scenario = "MTP", ModelType = "Multiple", Version = 2,
        `Creation Day` = cday, `Creation Month` = cmonth, `Creation Year` = cyear,
        `Day of Value` = day(date), `Month of Value` = month(date), `Year of Value` = year(date),
        AgeBand = age_group, Geography = population, ValueType = outname, Value = `50%`,
        `5%`, `10%`, `15%`, `20%`, `25%`, `30%`, `35%`, `40%`, `45%`, `50%`, `55%`, `60%`, `65%`, `70%`, `75%`, `80%`, `85%`, `90%`, `95%`)]

    names(data_out)[names(data_out) %like% "\\%$"] = paste("Quantile", quants);
    data_out
}

SPIM_output_full = function(test0, cyear, cmonth, cday, ymd_from, ymd_to)
{
    cat("Restricting time span...\n")
    t0 = as.numeric(ymd(ymd_from) - ymd("2020-01-01"))
    t1 = as.numeric(ymd(ymd_to) - ymd("2020-01-01"))
    test = test0[t %between% c(t0, t1)]

    cat("Adding age groups...\n");
    test[group %between% c(1, 1), age_group := "0-4"]
    test[group %between% c(2, 3), age_group := "5-14"]
    test[group %between% c(4, 5), age_group := "15-24"]
    test[group %between% c(6, 9), age_group := "25-44"]
    test[group %between% c(10, 13), age_group := "45-64"]
    test[group %between% c(14, 15), age_group := "65-74"]
    test[group %between% c(16, 16), age_group := "75+"]
    test[, age_group := factor(age_group, levels = unique(age_group))]

    cat("Summarizing variables...\n");
    summ = test[, .(death_o = sum(death_o + death2_o), disp = mean(disp_deaths)), by = .(run, t, population, age_group)]
    summ3 = test[, .(icu_p = sum(icu_p + icu2_p), disp = mean(disp_icu_prev)), by = .(run, t, population, age_group)]
    summ4 = test[, .(bed_p = sum(pmax(0, hosp_p + hosp2_p - hosp_undetected_p - hosp_undetected2_p)), disp = mean(disp_hosp_prev)), by = .(run, t, population, age_group)]
    summ34 = test[, .(admissions = sum(hosp_undetected_o + hosp_undetected2_o), disp = disp_hosp_inc), by = .(run, t, population, age_group)]
    summ_inf_i = test[, .(infections_i = sum(pcr_positive_i)), by = .(run, t, population, age_group)]
    summ_inf_p = test[, .(infections_p = sum(pcr_positive_p)), by = .(run, t, population, age_group)];
    summ_sero_p = test[, .(sero_p = sum(lfia_positive_p)), by = .(run, t, population, age_group)];

    cat("Running quantiles...\n");
    w = rbind(
        SPIM_output_1var(summ, "death_o", "type28_death_inc_line", cyear, cmonth, cday, ymd_from, ymd_to, "disp"),
        SPIM_output_1var(summ3, "icu_p", "icu_prev", cyear, cmonth, cday, ymd_from, ymd_to, "disp"),
        SPIM_output_1var(summ4, "bed_p", "hospital_prev", cyear, cmonth, cday, ymd_from, ymd_to, "disp"),
        SPIM_output_1var(summ34, "admissions", "hospital_inc", cyear, cmonth, cday, ymd_from, ymd_to, "disp"),
        SPIM_output_1var(summ_inf_p, "infections_p", "prevalence_mtp", cyear, cmonth, cday, ymd_from, ymd_to, -1),
        SPIM_output_1var(summ_inf_i, "infections_i", "infections_inc", cyear, cmonth, cday, ymd_from, ymd_to, -1),
        SPIM_output_1var(summ_sero_p, "sero_p", "sero_prev", cyear, cmonth, cday, ymd_from, ymd_to, -1)
    )
    return (w)
}

make_data = function(ld, sitreps, virus, sero)
{
    rbind(
        ld[, .(ValueType = "type28_death_inc_line", Geography = name,
            dmin = as.Date(NA), d = as.Date(date), dmax = as.Date(NA), ymin = NA, y = N, ymax = NA)],
        sitreps[, .(ValueType = "icu_prev", Geography = name,
            dmin = as.Date(NA), d = as.Date(date), dmax = as.Date(NA), ymin = NA, y = n_in_itu, ymax = NA)],
        sitreps[, .(ValueType = "hospital_prev", Geography = name,
            dmin = as.Date(NA), d = as.Date(date), dmax = as.Date(NA), ymin = NA, y = n_in_all_beds, ymax = NA)],
        sitreps[, .(ValueType = "hospital_inc", Geography = name,
            dmin = as.Date(NA), d = as.Date(date), dmax = as.Date(NA), ymin = NA, y = n_admitted_diagnosed, ymax = NA)],
        virus[Data.source %like% "REACT", .(ValueType = "prevalence_mtp", Geography = NHS.region,
            dmin = as.Date(Start.date), d = as.Date(Start.date) + (as.Date(End.date) - as.Date(Start.date)) / 2, dmax = as.Date(End.date), 
            ymin = pct(Lower.bound), y = pct(Central.estimate), ymax = pct(Upper.bound))],
        sero[!Data.source %like% "NHS BT", .(ValueType = "sero_prev", Geography = NHS.region,
            dmin = as.Date(Start.date), d = as.Date(Start.date) + (as.Date(End.date) - as.Date(Start.date)) / 2, dmax = as.Date(End.date), 
            ymin = pct(Lower.bound), y = pct(Central.estimate), ymax = pct(Upper.bound))]
    )
}
