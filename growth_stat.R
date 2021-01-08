library(data.table)
library(RCppMCMC)
library(ggplot2)
library(lubridate)
library(cowplot)
library(stringr)

var2 = fread("./fitting_data/var2-2020-12-21.csv")
sgtf = fread("./fitting_data/sgtf-2021-01-06.csv")

# Convert to same format and order
var2 = var2[, .(date = as.Date(sample_date), nhs_name, voc = var2, other = all - var2)][order(nhs_name, date)];
sgtf = sgtf[, .(date = as.Date(date), nhs_name, voc = sgtf, other)][order(nhs_name, date)];

# Helper functions
nameval = function(names, values)
{
    x = values;
    names(x) = names;
    return (x)
}

logistic = function(x, a, b)
{
    xx = a * (x - b);
    ifelse(xx < -200, 0, ifelse(xx > 200, 1, exp(xx) / (1 + exp(xx))))
}
    
cpp_vec = function(x) paste("{", paste(x, collapse = ", "), "}")
cpp_bbinom = 
'auto bbinom = [](double k, double n, double mode, double conc)
{
    auto lgamma = [](double x) { return gsl_sf_lngamma(x); };
    double a = mode * (conc - 2) + 1;
    double b = (1 - mode) * (conc - 2) + 1;
    
    return (lgamma(n + 1) + lgamma(k + a) + lgamma(n - k + b) + lgamma(a + b))
        - (lgamma(k + 1) + lgamma(n - k + 1) + lgamma(n + a + b) + lgamma(a) + lgamma(b));
};';
cpp_logistic = 
'auto logistic = [](double x, double a, double b)
{
    double xx = a * (x - b);
    if (xx < -200) return 0.0;
    if (xx > 200)  return 1.0;
    return exp(xx) / (1 + exp(xx));
};';
    
build_model = function(data, falsepos)
{
    regions = data[, unique(nhs_name)];
    n_regions = length(regions);
    
    # Set up priors
    prior_intercept_name = paste0("intercept", 1:n_regions);
    prior_intercept_dist = rep("N 0 1000", n_regions);
    prior_intercept_code = paste0("double intercept[] = {", paste(prior_intercept_name, collapse = ", "), "};");
    prior_intercept = nameval(prior_intercept_name, prior_intercept_dist);
    
    prior_conc_name = paste0("conc", 1:n_regions);
    prior_conc_dist = rep("N 0 500 T 2 2000", n_regions);
    prior_conc_code = paste0("double conc[] = {", paste(prior_conc_name, collapse = ", "), "};");
    prior_conc = nameval(prior_conc_name, prior_conc_dist);
    
    if (is.character(falsepos)) {
        prior_falsepos_name = paste0("falsepos", 1:n_regions);
        prior_falsepos_code = paste0("double falsepos[] = {", paste(prior_falsepos_name, collapse = ", "), "};");
        prior_falsepos_dist = rep(falsepos, n_regions);
        prior_falsepos = nameval(prior_falsepos_name, prior_falsepos_dist);
    } else {
        prior_falsepos_code = paste0("double falsepos[] = ", cpp_vec(rep_len(falsepos, n_regions)), ";");
        prior_falsepos = NULL;
    }
    
    # Make priors
    priors = c(
        slope = "N 0 1 T 0 100",
        prior_intercept,
        prior_conc,
        prior_falsepos
    );
    
    # Make code
    code = glue::glue(
        cpp_bbinom,
        cpp_logistic,
        prior_intercept_code,
        prior_conc_code,
        prior_falsepos_code,
        'std::vector<double> t = ${cpp_vec(as.numeric(data$date - ymd("2020-01-01")))};',
        'std::vector<double> s = ${cpp_vec(data$voc)};',
        'std::vector<double> f = ${cpp_vec(data$other)};',
        'std::vector<unsigned int> r = ${cpp_vec(match(data$nhs_name, regions) - 1)};',
        'for (unsigned int i = 0; i < t.size(); ++i) {',
        '    double frequency = logistic(t[i], slope, intercept[r[i]]);',
        '    double predicted = frequency + (1 - frequency) * falsepos[r[i]];',
        '    ll += bbinom(s[i], s[i] + f[i], predicted, conc[r[i]]);',
        '}',
    .sep = "\n", .open = "${", .close = "}")
    
    
    make_model("bbinom", priors, code)
}

# Shared slope
sgtf_model = build_model(sgtf, "B 1.5 15")
results = RCppMCMC(sgtf_model, 30000, 1000, threads = 6, verbose = TRUE)
setDT(results)
results

results[, 0.5 * var(-2 * ll) + mean(-2 * ll)]

# Independent slopes
nhs_regions = sgtf[, unique(nhs_name)]
results_1 = NULL;
ll = NULL;
for (i in seq_along(nhs_regions))
{
    sgtf_model_1 = build_model(sgtf[nhs_name == nhs_regions[i]], "B 1.5 15")
    results_1_0 = RCppMCMC(sgtf_model_1, 10000, 1000, threads = 6, verbose = TRUE)
    setnames(results_1_0, "slope", paste0("slope", i))
    setnames(results_1_0, "intercept1", paste0("intercept", i))
    setnames(results_1_0, "conc1", paste0("conc", i))
    setnames(results_1_0, "falsepos1", paste0("falsepos", i))
    ll = cbind(ll, results_1_0$ll)
    if (is.null(results_1)) {
        results_1 = results_1_0
    } else {
        results_1 = cbind(results_1, results_1_0[, 5:ncol(results_1_0)])
    }
}

# Get false positives by region
setDT(results_1)
results_1[, colMeans(.SD), .SDcols = patterns("^falsepos")]
results_1[, lapply(.SD, sd), .SDcols = patterns("^falsepos")]


ll_1 = rowSums(ll)
ll_1
setDT(results_1)
0.5 * var(-2 * ll_1) + mean(-2 * ll_1)
mean(ll*7)
results[, mean(ll)]


# Plotting of results
plot_results = function(data_list, data_names, data_colours, results, date_min = "1900-01-01", date_max = "2100-01-01")
{
    # Data
    data_all = NULL;
    nhs_names = data_list[[1]][, unique(nhs_name)];
    
    for (i in seq_along(data_list))
    {
        data = copy(data_list[[i]][date >= date_min & date <= date_max & !nhs_name %in% c("Wales", "Scotland", "Northern Ireland")]);
        data[, qlo := qbeta(0.025, voc + 1, other + 1)];
        data[, qhi := qbeta(0.975, voc + 1, other + 1)];
        data[, name := data_names[i]];
        data_all = rbind(data_all, data);
    }
    
    plot = ggplot(data_all) +
        geom_ribbon(aes(x = date, ymin = qlo, ymax = qhi, fill = name), alpha = 0.4) +
        geom_line(aes(x = date, y = voc / (voc + other), colour = name)) +
        scale_colour_manual(values = nameval(data_names, data_colours), aesthetics = c("fill", "colour")) +
        scale_x_date(date_breaks = "1 month", date_labels = "%b") +
        facet_wrap(~nhs_name, nrow = 1) +
        labs(x = NULL, y = "Variant-\ncompatible frequency", colour = NULL, fill = NULL) +
        theme_cowplot(font_size = 10) +
        theme(strip.background = element_blank(), legend.position = c(0.01, 0.9))
    
    # If present, results
    if (!is.null(results))
    {
        res_data_all = NULL;
        
        res = melt(results, id.vars = numeric(), measure.vars = 5:ncol(results));
        res[, nhs_name := nhs_names[as.numeric(str_remove_all(variable, "[a-z]*"))]];
        res[, var := str_remove_all(variable, "[0-9]*")];
        
        extract = function(res, nhs, varname)
        {
            ex = res[nhs_name == nhs & var == varname, value];
            if (length(ex) > 0)
                return (ex);
            return (res[is.na(nhs_name) & var == varname, value]);
        }
        
        set.seed(12345);
        for (i in seq_along(nhs_names))
        {
            # Extract posterior for this NHS region
            res_data = data.table(nhs_name = nhs_names[i], date = data_all[, min(date)] + 0:data_all[, as.numeric(max(date) - min(date))]);
            slope = extract(res, nhs_names[i], "slope");
            intercept = extract(res, nhs_names[i], "intercept");
            conc = extract(res, nhs_names[i], "conc");
            falsepos = extract(res, nhs_names[i], "falsepos");

            nsamp = 2000;
            rows = sample(length(slope), nsamp);
            slope = slope[rows];
            intercept = intercept[rows];
            conc = conc[rows];
            falsepos = falsepos[rows];
            
            pred  = matrix(0, nrow = nrow(res_data), ncol = nsamp);
            predr = matrix(0, nrow = nrow(res_data), ncol = nsamp);
            sgtfv = matrix(0, nrow = nrow(res_data), ncol = nsamp);
            
            for (j in 1:nsamp)
            {
                pred[, j] = logistic(as.numeric(res_data$date - ymd("2020-01-01")), slope[j], intercept[j]);
                sgtfv[, j] = pred[, j] / (pred[, j] + (1 - pred[, j]) * falsepos[j]);
                pred[, j] = pred[, j] + (1 - pred[, j]) * falsepos[j];
                predr[, j] = rbeta(length(pred[, j]), pred[, j] * (conc[j] - 2) + 1, (1 - pred[, j]) * (conc[j] - 2) + 1);
            }
            
            res_data[, mlo := apply(pred,  1, function(x) quantile(x, 0.025))];
            res_data[, mmd := apply(pred,  1, function(x) quantile(x, 0.500))];
            res_data[, mhi := apply(pred,  1, function(x) quantile(x, 0.975))];
            res_data[, rlo := apply(predr, 1, function(x) quantile(x, 0.025))];
            res_data[, rmd := apply(predr, 1, function(x) quantile(x, 0.500))];
            res_data[, rhi := apply(predr, 1, function(x) quantile(x, 0.975))];
            res_data[, sgtfv := apply(sgtfv, 1, mean)];
            
            res_data_all = rbind(res_data_all, res_data);
        }
        
        plot = plot +
            geom_ribbon(data = res_data_all, aes(x = date, ymin = rlo, ymax = rhi), alpha = 0.4) +
            geom_ribbon(data = res_data_all, aes(x = date, ymin = mlo, ymax = mhi), alpha = 0.4)
    }
    
    plot
}

plot0 = plot_results(list(sgtf, var2), c("SGTF", "Sequence"), c("darkorchid", "darkorange2"), NULL, date_min = "2020-10-01")
plot1 = plot_results(list(sgtf, var2), c("SGTF", "Sequence"), c("darkorchid", "darkorange2"), results, date_min = "2020-10-01")
plot2 = plot_results(list(sgtf, var2), c("SGTF", "Sequence"), c("darkorchid", "darkorange2"), results_1, date_min = "2020-10-01")

# sgtfv = plot_results(list(sgtf, var2), c("SGTF", "Sequence"), c("darkorchid", "darkorange2"), results, date_min = "2020-08-01") # with func returning res_data_all instead of plot
# fwrite(sgtfv[, .(nhs_name, date, sgtfv)], "./data/sgtfvoc.csv")

plot012 = plot_grid(plot0, plot1, plot2, ncol = 1, labels = LETTERS, label_size = 10)
ggsave("./output/growth_stat.pdf", plot012, width = 34, height = 12, units = "cm", useDingbats = FALSE)
ggsave("./output/growth_stat.png", plot012, width = 34, height = 12, units = "cm")

ggsave("./output/growth_data.pdf", plot0, width = 20, height = 15, units = "cm", useDingbats = FALSE)
ggsave("./output/growth_data.png", plot0, width = 20, height = 15, units = "cm")
ggsave("./output/growth_data_model.pdf", plot1, width = 20, height = 15, units = "cm", useDingbats = FALSE)
ggsave("./output/growth_data_model.png", plot1, width = 20, height = 15, units = "cm")

# Table of results
niceq = function(x)
{
    q = quantile(x, c(0.05, 0.5, 0.95));
    f = function(y) prettyNum(signif(y, 3), big.mark = ",")
    paste0(f(q[2]), " (", f(q[1]), " - ", f(q[3]), ")")
}

nicedate = function(x, origin = "2020-01-01")
{
    q = quantile(x, c(0.05, 0.5, 0.95));
    f = function(y) format(as.Date(y + ymd(origin)), "%d %b")
    paste0(f(q[2]), " (", f(q[1]), " - ", f(q[3]), ")")
}

table_results = function(results)
{
    res = melt(results, id.vars = numeric(), measure.vars = 5:ncol(results));
    res[, nhs_name := nhs_regions[as.numeric(str_remove_all(variable, "[a-z]*"))]];
    res[, var := str_remove_all(variable, "[0-9]*")];
    
    extract = function(res, nhs, varname)
    {
        ex = res[nhs_name == nhs & var == varname, value];
        if (length(ex) > 0)
            return (ex);
        return (res[is.na(nhs_name) & var == varname, value]);
    }
    
    tb = NULL;
    
    for (i in seq_along(nhs_regions))
    {
        # Extract posterior for this NHS region
        slope = extract(res, nhs_regions[i], "slope");
        intercept = extract(res, nhs_regions[i], "intercept");
        conc = extract(res, nhs_regions[i], "conc");
        falsepos = extract(res, nhs_regions[i], "falsepos");
        
        # Put in table
        tb = rbind(tb,
            data.table(`NHS region` = nhs_regions[i], 
                `Relative growth rate` = niceq(slope), 
                `Intercept (f_VOC = 50%)` = nicedate(intercept),
                `SGTF false positive rate` = niceq(falsepos),
                `Data precision` = niceq(conc)
            )
        )
    }
    
    tb
}

tb = table_results(results)
tb_1 = table_results(results_1)
fwrite(tb, "./output/growth_stat_tb_single.csv")
fwrite(tb_1, "./output/growth_stat_tb_multiple.csv")
tb_1
# TODO Plot counts as well
plot
