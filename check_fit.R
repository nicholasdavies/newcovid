#
# PLOT FITS.
#

gen_fit = function(test, ld, sitreps, virus, sero, populations)
{
    test = copy(test)[population %in% populations]
    sero = copy(sero)[NHS.region %in% populations]
    virus = copy(virus)[NHS.region %in% populations]
    ld = copy(ld)[name %in% populations]
    sitreps = copy(sitreps)[name %in% populations]
    
    # Calculate total population
    popsize = NULL
    for (i in seq_along(parametersI)) {
        if (!is.null(parametersI[[i]])) {
            popsize = rbind(popsize,
                data.table(Geography = parametersI[[i]]$pop[[1]]$name, population_size = sum(parametersI[[i]]$pop[[1]]$size))
            )
        }
    }
    
    # Create formatted output
    output = SPIM_output_full(test, 2020, 11, 1, "2020-01-01", as.character(ymd("2020-01-01") + max(test$t)))
    output[, d := make_date(`Year of Value`, `Month of Value`, `Day of Value`)]
    output = merge(output, popsize, by = "Geography")
    
    adj_output = function(output, val_type, div, pop = 0) {
        output[ValueType == val_type, Value := Value / (div + population_size * pop)]
        output[ValueType == val_type, `Quantile 0.05` := `Quantile 0.05` / (div + population_size * pop)]
        output[ValueType == val_type, `Quantile 0.1`  := `Quantile 0.1`  / (div + population_size * pop)]
        output[ValueType == val_type, `Quantile 0.15` := `Quantile 0.15` / (div + population_size * pop)]
        output[ValueType == val_type, `Quantile 0.2`  := `Quantile 0.2`  / (div + population_size * pop)]
        output[ValueType == val_type, `Quantile 0.25` := `Quantile 0.25` / (div + population_size * pop)]
        output[ValueType == val_type, `Quantile 0.3`  := `Quantile 0.3`  / (div + population_size * pop)]
        output[ValueType == val_type, `Quantile 0.35` := `Quantile 0.35` / (div + population_size * pop)]
        output[ValueType == val_type, `Quantile 0.4`  := `Quantile 0.4`  / (div + population_size * pop)]
        output[ValueType == val_type, `Quantile 0.45` := `Quantile 0.45` / (div + population_size * pop)]
        output[ValueType == val_type, `Quantile 0.5`  := `Quantile 0.5`  / (div + population_size * pop)]
        output[ValueType == val_type, `Quantile 0.55` := `Quantile 0.55` / (div + population_size * pop)]
        output[ValueType == val_type, `Quantile 0.6`  := `Quantile 0.6`  / (div + population_size * pop)]
        output[ValueType == val_type, `Quantile 0.65` := `Quantile 0.65` / (div + population_size * pop)]
        output[ValueType == val_type, `Quantile 0.7`  := `Quantile 0.7`  / (div + population_size * pop)]
        output[ValueType == val_type, `Quantile 0.75` := `Quantile 0.75` / (div + population_size * pop)]
        output[ValueType == val_type, `Quantile 0.8`  := `Quantile 0.8`  / (div + population_size * pop)]
        output[ValueType == val_type, `Quantile 0.85` := `Quantile 0.85` / (div + population_size * pop)]
        output[ValueType == val_type, `Quantile 0.9`  := `Quantile 0.9`  / (div + population_size * pop)]
        output[ValueType == val_type, `Quantile 0.95` := `Quantile 0.95` / (div + population_size * pop)]
    }
    
    adj_output(output, "hospital_inc", 1)
    adj_output(output, "hospital_prev", 1)
    adj_output(output, "icu_prev", 1)
    adj_output(output, "prevalence_mtp", 0, 0.01)
    adj_output(output, "sero_prev", 0, 0.01)
    adj_output(output, "type28_death_inc_line", 1)
    
    # Make data to output
    data = make_data(ld, sitreps, virus, sero)
    data = merge(data, popsize, by = "Geography")
    
    adj_data = function(data, val_type, div, pop = 0) {
        data[ValueType == val_type, ymin := ymin / (div + population_size * pop)]
        data[ValueType == val_type, y    := y    / (div + population_size * pop)]
        data[ValueType == val_type, ymax := ymax / (div + population_size * pop)]
    }
    
    adj_data(data, "hospital_inc", 1)
    adj_data(data, "hospital_prev", 1)
    adj_data(data, "icu_prev", 1)
    adj_data(data, "prevalence_mtp", 0.01)
    adj_data(data, "sero_prev", 0.01)
    adj_data(data, "type28_death_inc_line", 1)
    
    output[ValueType == "hospital_inc", ValueType := "Hospital admissions"]
    output[ValueType == "hospital_prev", ValueType := "Hospital beds occupied"]
    output[ValueType == "icu_prev", ValueType := "ICU beds occupied"]
    output[ValueType == "infections_inc", ValueType := "Infection incidence"]
    output[ValueType == "prevalence_mtp", ValueType := "PCR prevalence (%)"]
    output[ValueType == "sero_prev", ValueType := "Seroprevalence (%)"]
    output[ValueType == "type28_death_inc_line", ValueType := "Deaths"]

    data[ValueType == "hospital_inc", ValueType := "Hospital admissions"]
    data[ValueType == "hospital_prev", ValueType := "Hospital beds occupied"]
    data[ValueType == "icu_prev", ValueType := "ICU beds occupied"]
    data[ValueType == "infections_inc", ValueType := "Infection incidence"]
    data[ValueType == "prevalence_mtp", ValueType := "PCR prevalence (%)"]
    data[ValueType == "sero_prev", ValueType := "Seroprevalence (%)"]
    data[ValueType == "type28_death_inc_line", ValueType := "Deaths"]
    
    return (list(data, output))
}

check_fit = function(test, ld, sitreps, virus, sero, populations, max_date, min_date = NULL)
{
    fit = gen_fit(test, ld, sitreps, virus, sero, populations)
    data = fit[[1]]
    output = fit[[2]]
    output = output[d <= max_date]
    
    if (!is.null(min_date)) {
        output = output[d >= min_date]
        data = data[d >= min_date]
    }
    
    # Make plot
    theme_set(cowplot::theme_cowplot(font_size = 10) + theme(strip.background = element_blank()))
    
    linetypes = c("Deaths", "Hospital admissions", "Hospital beds occupied", "ICU beds occupied")
    
    plot = ggplot(output[d > "2020-03-01" & AgeBand == "All"]) +
        geom_ribbon(aes(x = d, ymin = `Quantile 0.05`, ymax = `Quantile 0.95`, fill = ValueType), alpha = 0.5) +
        geom_line(aes(x = d, y = Value, colour = ValueType)) +
        geom_line(data = data[ValueType %in% linetypes], aes(x = d, y = y), size = 0.2) +
        geom_point(data = data[!ValueType %in% linetypes], aes(x = d, y = y), size = 0.01, shape = 20) +
        geom_linerange(data = data, aes(x = d, ymin = ymin, ymax = ymax), size = 0.2) +
        geom_linerange(data = data, aes(xmin = dmin, xmax = dmax, y = y), size = 0.2) +
        facet_grid(ValueType ~ Geography, scales = "free", switch = "y") +
        theme(legend.position = "none", strip.placement = "outside", strip.background = element_blank()) +
        scale_x_date(date_breaks = "2 months", date_labels = "%b") +
        labs(x = NULL, y = NULL)
    
    return (plot)
}



compare_fit = function(test, test0, ld, sitreps, virus, sero, populations, populations0, max_date)
{
    fit = gen_fit(test, ld, sitreps, virus, sero, populations)
    output = fit[[2]]

    fit = gen_fit(test0, ld, sitreps, virus, sero, populations0)
    data = fit[[1]]
    output0 = fit[[2]]
    
    output[, kind := "With VOC"]
    output0[, kind := "Without VOC"]
    output = rbind(output, output0)
    output = output[d <= max_date]
    
    # Make plot
    theme_set(cowplot::theme_cowplot(font_size = 10) + theme(strip.background = element_blank()))
    
    linetypes = c("Deaths", "Hospital admissions", "Hospital beds occupied", "ICU beds occupied")
    
    plot = ggplot(output[d > "2020-03-01" & AgeBand == "All"]) +
        geom_ribbon(aes(x = d, ymin = `Quantile 0.05`, ymax = `Quantile 0.95`, fill = ValueType, group = kind), alpha = 0.5) +
        geom_line(aes(x = d, y = Value, colour = ValueType, linetype = kind)) +
        geom_line(data = data[ValueType %in% linetypes], aes(x = d, y = y), size = 0.2) +
        geom_point(data = data[!ValueType %in% linetypes], aes(x = d, y = y), size = 0.01, shape = 20) +
        geom_linerange(data = data, aes(x = d, ymin = ymin, ymax = ymax), size = 0.2) +
        geom_linerange(data = data, aes(xmin = dmin, xmax = dmax, y = y), size = 0.2) +
        facet_grid(ValueType ~ Geography, scales = "free", switch = "y") +
        theme(legend.position = "none", strip.placement = "outside", strip.background = element_blank()) +
        scale_x_date(date_breaks = "2 months", date_labels = "%b") +
        labs(x = NULL, y = NULL)
    
    return (plot)
}
