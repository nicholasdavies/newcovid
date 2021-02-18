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

arrange_projection = function(proj, from_date = NULL, england = FALSE)
{
    ft = 0;
    if (!is.null(from_date)) {
        ft = as.numeric(ymd(from_date) - ymd("2020-01-01"))
    }
    
    if (england == FALSE) {
        w = proj[t >= ft, .(deaths = sum(death_o + death2_o), admissions = sum(hosp_undetected_o + hosp_undetected2_o),
            beds = sum(hosp_p + hosp2_p - hosp_undetected_p - hosp_undetected2_p), icu = sum(icu_p + icu2_p), 
            disp_deaths = mean(disp_deaths), disp_hosp_inc = mean(disp_hosp_inc), disp_hosp_prev = mean(disp_hosp_prev), disp_icu_prev = mean(disp_icu_prev),
            Rt = obs0[1], tier = obs0[2], cb = obs0[3]), keyby = .(run, t, population)]
    } else {
        w = proj[t >= ft, .(population = "England", deaths = sum(death_o + death2_o), admissions = sum(hosp_undetected_o + hosp_undetected2_o),
            beds = sum(hosp_p + hosp2_p - hosp_undetected_p - hosp_undetected2_p), icu = sum(icu_p + icu2_p), 
            disp_deaths = mean(disp_deaths), disp_hosp_inc = mean(disp_hosp_inc), disp_hosp_prev = mean(disp_hosp_prev), disp_icu_prev = mean(disp_icu_prev),
            Rt = mean(obs0[seq(1, .N, by = 16)]), tier = mean(obs0[seq(2, .N, by = 16)]), cb = mean(obs0[seq(3, .N, by = 16)])), keyby = .(run, t)]
    }
    
    quants = function(x, qu, disp)
    {
        size = 1/(disp^2);
        qnbinom(qu, size = size, mu = x)
    }
    
    w[, deaths_lo := quants(deaths, 0.025, disp_deaths)]
    w[, deaths_md := quants(deaths, 0.500, disp_deaths)]
    w[, deaths_hi := quants(deaths, 0.975, disp_deaths)]
    w[, admissions_lo := quants(admissions, 0.025, disp_hosp_inc)]
    w[, admissions_md := quants(admissions, 0.500, disp_hosp_inc)]
    w[, admissions_hi := quants(admissions, 0.975, disp_hosp_inc)]
    w[, bedsK_lo := quants(beds, 0.025, disp_hosp_prev) / 1000]
    w[, bedsK_md := quants(beds, 0.500, disp_hosp_prev) / 1000]
    w[, bedsK_hi := quants(beds, 0.975, disp_hosp_prev) / 1000]
    w[, icu_lo := quants(icu, 0.025, disp_icu_prev)]
    w[, icu_md := quants(icu, 0.500, disp_icu_prev)]
    w[, icu_hi := quants(icu, 0.975, disp_icu_prev)]
    
    w = cbind(
        w[, lapply(.SD, mean), .SDcols = patterns("_lo$|_md$|_hi$"), keyby = .(t, population)],
        w[, .(Rt_lo = quantile(Rt, 0.025)), keyby = .(t, population)][, .(Rt_lo)],
        w[, .(Rt_md = quantile(Rt, 0.500)), keyby = .(t, population)][, .(Rt_md)],
        w[, .(Rt_hi = quantile(Rt, 0.975)), keyby = .(t, population)][, .(Rt_hi)]
    )

    w = melt(w, id.vars = 1:2)
    w[, quant := str_remove_all(variable, ".*_")];
    w[, variable := str_remove_all(variable, "_lo$|_md$|_hi$")];
    w = dcast(w, t + population + variable ~ quant, value.var = "value")

    w[variable == "deaths", variable := "Deaths"]
    w[variable == "admissions", variable := "Hospital\nadmissions"]
    w[variable == "bedsK", variable := "Hospital beds\noccupied (thousands)"]
    w[variable == "icu", variable := "ICU beds\noccupied"]
    
    setnames(w, c("t", "population", "variable", "97.5%", "2.5%", "50%"))
    return (w)
}

plot_projection = function(proj_list, proj_names, from_date, england = FALSE, reduced = FALSE, pal = "Accent", hosp_line = NA, deaths_line = NA, colour_label = NULL)
{
    p = NULL
    for (i in seq_along(proj_list)) {
        pp = arrange_projection(proj_list[[i]], england = england);
        pp[, name := proj_names[[i]]];
        p = rbind(p, pp);
    }
    
    if (reduced) {
        p = p[variable %in% c("cb", "Rt", "Hospital beds\noccupied (thousands)", "Deaths")];
    }
    
    p[, name := factor(name, unique(name))]
    p[, variable := factor(variable, unique(variable))]

    plot = ggplot(p[!variable %in% c("tier", "cb") & ymd("2020-01-01") + t >= from_date]);
    
    cbs = p[variable == "cb", .(t = ymd("2020-01-01") + t[`50%` > 0.5]), by = .(population)][, 
        .(tmin = min(t), tmax = max(t)), by = population]
    if (nrow(cbs) > 0) {
        plot = plot +
            geom_rect(data = cbs, aes(xmin = tmin, xmax = tmax, ymin = -Inf, ymax = Inf), fill = "black", alpha = 0.1)
    }
    
    rline = data.table(variable = factor(c("Rt", "Hospital beds\noccupied (thousands)", "Deaths"), levels(p$variable)), y = c(1, hosp_line, deaths_line))
    
    plot +
        geom_ribbon(aes(x = ymd("2020-01-01") + t, ymin = `2.5%`, ymax = `97.5%`, fill = name), alpha = 0.5) +
        geom_hline(data = rline, aes(yintercept = y), size = 0.3, linetype = "33") +
        geom_line(aes(x = ymd("2020-01-01") + t, y = `50%`, colour = name)) +
        facet_grid(variable ~ population, switch = "y", scales = "free") +
        cowplot::theme_cowplot(font_size = 11) + 
        theme(strip.background = element_blank(), strip.placement = "outside", 
            legend.position = "bottom",
            panel.background = element_rect(fill = "#f4f4f4"),
            panel.grid.major = element_line(colour = "#ffffff", size = 0.4)) +
        labs(x = NULL, y = NULL, colour = colour_label, fill = colour_label) +
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

# plot_cum_deaths = function(proj_list, proj_names, from_date, ystart = 26, ydiff = -1.2, titl = "Cumulative deaths")
# {
#     stop("Check for both virus strains.")
#     p = NULL
#     for (i in seq_along(proj_list)) {
#         w = arrange_projection(proj_list[[i]], cumulative_deaths = TRUE, from_date = from_date, england = TRUE)
#         p = rbind(p,
#             cbind(w[variable == "Cumulative deaths" | variable == "cb"], name = proj_names[i])
#         )
#     }
#     
#     p[, name := factor(name, levels = unique(name))]
#     
#     cbs = p[variable == "cb", .(t = ymd("2020-01-01") + t[`50%` > 0.5]), by = .(population, name)][, 
#         .(tmin = min(t), tmax = max(t)), by = .(population, name)]
#     
#     plot = ggplot(p[variable != "cb"]) +
#         geom_ribbon(aes(x = ymd("2020-01-01") + t, ymin = `2.5%` / 1000, ymax = `97.5%` / 1000, fill = name), alpha = 0.5) +
#         geom_line(aes(x = ymd("2020-01-01") + t, y = `50%` / 1000, colour = name), size = 0.2) +
#         facet_wrap(~population, scales = "free") +
#         theme(strip.background = element_blank(), strip.text = element_blank(), legend.position = c(0.02, 0.8),
#             panel.background = element_rect(fill = "#f4f4f4"),
#             panel.grid.major = element_line(colour = "#ffffff", size = 0.2)) +
#         labs(x = NULL, y = "Cumulative deaths\n(thousands)", colour = NULL, fill = NULL, title = titl) +
#         scale_x_date(date_breaks = "1 month", date_labels = "%b")
#     
#     if (nrow(cbs) > 0) {
#         yy = ystart
#         for (nm in cbs[, unique(name)]) {
#             plot = plot +
#                 geom_linerange(data = cbs[name == nm], aes(xmin = tmin, xmax = tmax, colour = name), y = yy, show.legend = FALSE, linetype = "32", size = 0.2) +
#                 geom_linerange(data = cbs[name == nm], aes(x = tmin, colour = name), ymin = yy - (ydiff * 0.2), ymax = yy + (ydiff * 0.2), show.legend = FALSE, size = 0.2) +
#                 geom_linerange(data = cbs[name == nm], aes(x = tmax, colour = name), ymin = yy - (ydiff * 0.2), ymax = yy + (ydiff * 0.2), show.legend = FALSE, size = 0.2)
#             yy = yy + ydiff
#         }
#     }
# 
#     return (plot)
# }

# plot_icu = function(proj_list, proj_names, from_date)
# {
#     stop("Check for both virus strains.")
# 
#     p = NULL
#     for (i in seq_along(proj_list)) {
#         w = arrange_projection(proj_list[[i]], cumulative_deaths = TRUE, from_date = from_date, england = TRUE)
#         p = rbind(p,
#             cbind(w[variable == "ICU beds"], name = proj_names[i])
#         )
#     }
#     
#     p[, name := factor(name, levels = unique(name))]
#     
#     plot = ggplot(p) +
#         geom_ribbon(aes(x = ymd("2020-01-01") + t, ymin = `2.5%`, ymax = `97.5%`, fill = name), alpha = 0.5) +
#         geom_line(aes(x = ymd("2020-01-01") + t, y = `50%`, colour = name), size = 0.2) +
#         facet_wrap(~population, scales = "free") +
#         theme(strip.background = element_blank(), strip.text = element_blank(), legend.position = "none",
#             panel.background = element_rect(fill = "#f4f4f4"),
#             panel.grid.major = element_line(colour = "#ffffff", size = 0.2)) +
#         labs(x = NULL, y = "ICU beds occupied", colour = NULL, fill = NULL) +
#         scale_x_date(date_breaks = "1 month", date_labels = "%b")
# 
#     return (plot)
# }


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

    # Factor
    totals[, population := factor(population, geo)]
    peaks[, population := factor(population, geo)]

    tab = rbind(
        totals[when == wh, .(indicator = "Total deaths", value = niceq(deaths)), keyby = population],
        totals[when == wh, .(indicator = "Total admissions", value = niceq(admissions)), keyby = population],
        peaks[, .(indicator = "Peak deaths", value = niceq(dea)), keyby = population],
        peaks[, .(indicator = "Peak ICU requirement", value = niceq(icu)), keyby = population],
        peaks[, .(indicator = "Peak ICU (rel. to 1st wave)", value = nicepc(rel_icu)), keyby = population]
    )
    return (dcast(tab, indicator ~ population))
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

# combine_R = function(R_ldN4o, R_ldN4c) {
#     R_ldN4 = cbind(R_ldN4o[, 1:3], R_ldN4c[, 3], R_ldN4o[, 4], R_ldN4c[, 4])
#     names(R_ldN4)[c(3,5)] = paste0(names(R_ldN4)[c(3,5)], ", schools open")
#     names(R_ldN4)[c(4,6)] = paste0(names(R_ldN4)[c(4,6)], ", schools closed")
#     R_ldN4
# }

# plot_indicator = function(series, indicator_name)
# {
#     ggplot(series[indicator == indicator_name]) +
#         geom_pointrange(aes(x = population, ymin = lo, y = mid, ymax = hi, colour = scenario), position = position_dodge(width = 0.8), shape = 20, size = 0.0) +
#         geom_linerange(aes(x = population, ymin = lo, y = mid, ymax = hi, colour = scenario), position = position_dodge(width = 0.8), size = 0.2) +
#         labs(x = NULL, y = indicator_name) + ylim(0, NA)
# }
# 
# plot_indicators_england = function(series, indicator_names, y_axis_title, unit = 1, legpos = "none", pal = "Dark2")
# {
#     ggplot(series[population == "England" & indicator %in% indicator_names]) +
#         geom_ribbon(aes(x = scenario, ymin = lo / unit, ymax = hi / unit, fill = indicator, group = indicator), alpha = 0.5) +
#         geom_line(aes(x = scenario, y = mid / unit, colour = indicator, group = indicator), size = 0.2) +
#         labs(x = NULL, y = y_axis_title, fill = NULL, colour = NULL) + ylim(0, NA) +
#         theme(panel.background = element_rect(fill = "#f4f4f4"),
#             panel.grid.major = element_line(colour = "#ffffff", size = 0.2),
#             legend.position = legpos) +
#         scale_color_brewer(aesthetics = c("colour", "fill"), palette = pal)
# }
# 
# plot_indicators_england2 = function(series, indicator_names, y_axis_title, unit = 1, legpos = "none", pal = "Accent")
# {
#     ggplot(series[population == "England" & indicator %in% indicator_names]) +
#         geom_linerange(aes(x = scenario, ymin = lo / unit, ymax = hi / unit, colour = indicator, group = indicator),
#             position = position_dodge(width = 0.5), size = 0.4) +
#         geom_point(aes(x = scenario, y = mid / unit, colour = indicator, group = indicator),
#             position = position_dodge(width = 0.5), size = 0.2, shape = 3) +
#         labs(x = NULL, y = y_axis_title, colour = NULL) + ylim(0, NA) +
#         theme(panel.background = element_rect(fill = "#f4f4f4"),
#             panel.grid.major = element_line(colour = "#ffffff", size = 0.4),
#             legend.position = legpos) +
#         scale_color_brewer(aesthetics = c("colour", "fill"), palette = pal)
# }
# 
# plot_indicators_england3 = function(series, indicator_names, y_axis_title, unit = 1, legpos = "none", pal = "Accent")
# {
#     ggplot(series[population == "England" & indicator %in% indicator_names]) +
#         geom_col(aes(x = scenario, y = mid / unit, fill = indicator, group = indicator), 
#             position = position_dodge(width = 0.6), size = 0.25, colour = "black", width = 0.4) +
#         geom_linerange(aes(x = scenario, ymin = lo / unit, ymax = hi / unit, group = indicator),
#             position = position_dodge(width = 0.6), size = 0.25) +
#         labs(x = NULL, y = y_axis_title, fill = NULL) + ylim(0, NA) +
#         theme(panel.background = element_rect(fill = "#f4f4f4"),
#             panel.grid.major = element_line(colour = "#ffffff", size = 0.4),
#             legend.position = legpos) +
#         scale_color_brewer(aesthetics = c("colour", "fill"), palette = pal)
# }
# 
# plot_indicators_england_stack2 = function(series, indicator_names, y_axis_title, stack_order, legpos = "none", pal = "Set2")
# {
#     ser2 = copy(series[population == "England" & indicator %in% indicator_names])
#     ser2[, indicator := factor(indicator, stack_order)]
#     ser3 = copy(ser2)
#     ser3 = ser3[order(scenario, rev(indicator))]
#     ser3[, lo := lo - mid]
#     ser3[, hi := hi - mid]
#     ser3[, mid := cumsum(mid), by = scenario]
#     ser3[, lo := mid + lo]
#     ser3[, hi := mid + hi]
#     
#     ggplot() +
#         geom_col(data = ser2, aes(x = scenario, y = mid, fill = indicator), position = "stack", colour = "black", size = 0.25, width = 0.5) +
#         geom_linerange(data = ser3, aes(x = scenario, ymin = lo, ymax = hi), size = 0.25) +
#         labs(x = NULL, y = y_axis_title, fill = NULL, colour = NULL) + ylim(0, NA) +
#         theme(panel.background = element_rect(fill = "#f4f4f4"),
#             panel.grid.major = element_line(colour = "#ffffff", size = 0.2),
#             legend.position = legpos) +
#         scale_color_brewer(aesthetics = c("colour", "fill"), palette = pal)
# }
# 
# plot_indicators_england_stack = function(series, indicator_names, y_axis_title, stack_order, legpos = "none", pal = "Set2")
# {
#     ser2 = rlang::duplicate(series[population == "England" & indicator %in% indicator_names])
#     ser2[, indicator := factor(indicator, stack_order)]
#     ser3 = rlang::duplicate(ser2)
#     ser3 = ser3[order(scenario, rev(indicator))]
#     ser3[, lo := lo - mid]
#     ser3[, hi := hi - mid]
#     ser3[, mid := cumsum(mid), by = scenario]
#     ser3[, lo := mid + lo]
#     ser3[, hi := mid + hi]
# 
#     ggplot() +
#         geom_area(data = ser2, aes(x = scenario, y = mid, fill = indicator, group = indicator)) +
#         geom_ribbon(data = ser3, aes(x = scenario, ymin = lo, ymax = hi, group = indicator), alpha = 0.25) +
#         labs(x = NULL, y = y_axis_title, fill = NULL, colour = NULL) + ylim(0, NA) +
#         theme(panel.background = element_rect(fill = "#f4f4f4"),
#             panel.grid.major = element_line(colour = "#ffffff", size = 0.2),
#             legend.position = legpos) +
#         scale_color_brewer(aesthetics = c("colour", "fill"), palette = pal)
# }
# 
# plot_indicator_regions = function(series, indicator_name, label = indicator_name, suffix = "\n(Thousands)", unit = 1000, legpos = "none", bigtitle = NULL)
# {
#     ser2 = rlang::duplicate(series[population != "England" & indicator == indicator_name])
#     ser2[population == "East of England", population := "EE"]
#     ser2[population == "London", population := "Ldn"]
#     ser2[population == "Midlands", population := "Mlds"]
#     ser2[population == "North East and Yorkshire", population := "NE&Y"]
#     ser2[population == "North West", population := "NW"]
#     ser2[population == "South East", population := "SE"]
#     ser2[population == "South West", population := "SW"]
#     ggplot(ser2) +
#         geom_point(aes(x = population, y = mid / unit, colour = scenario), position = position_dodge(width = 0.8), shape = 20, size = 0.5) +
#         geom_linerange(aes(x = population, ymin = lo / unit, ymax = hi / unit, colour = scenario), position = position_dodge(width = 0.8), size = 0.3) +
#         labs(x = NULL, y = paste0(label, suffix), colour = NULL, title = bigtitle) + 
#         ylim(0, NA) +
#         theme(legend.position = legpos)
# }

make_vaccine_schedule = function(parametersI, ymd_start, weekly_v, targeting, popset)
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

# Vaccine targeting
targeting_old_to_young = 
    list(c(0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,1,1)*0.85,
         c(0,0,0,0, 0,0,0,0, 0,0,0,0, 1,1,1,1)*0.85,
         c(0,0,0,0, 0,0,0,0, 0,0,1,1, 1,1,1,1)*0.85,
         c(0,0,0,0, 0,0,0,0, 1,1,1,1, 1,1,1,1)*0.85,
         c(0,0,0,0, 0,0,1,1, 1,1,1,1, 1,1,1,1)*0.85,
         c(0,0,0,0, 1,1,1,1, 1,1,1,1, 1,1,1,1)*0.85,
         c(0,0,1,1, 1,1,1,1, 1,1,1,1, 1,1,1,1)*0.85,
         c(1,1,1,1, 1,1,1,1, 1,1,1,1, 1,1,1,1)*0.85)

targeting_young_to_old = 
    list(c(0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,1,1)*0.85,
         c(0,0,0,0, 0,0,0,0, 0,0,0,0, 1,1,1,1)*0.85,
         c(1,1,0,0, 0,0,0,0, 0,0,0,0, 1,1,1,1)*0.85,
         c(1,1,1,1, 0,0,0,0, 0,0,0,0, 1,1,1,1)*0.85,
         c(1,1,1,1, 1,1,0,0, 0,0,0,0, 1,1,1,1)*0.85,
         c(1,1,1,1, 1,1,1,1, 0,0,0,0, 1,1,1,1)*0.85,
         c(1,1,1,1, 1,1,1,1, 1,1,0,0, 1,1,1,1)*0.85,
         c(1,1,1,1, 1,1,1,1, 1,1,1,1, 1,1,1,1)*0.85)


