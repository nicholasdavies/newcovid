# 4-fitting.R

# covidm options
cm_path = "~/Dropbox/nCoV/covidm/"; ### CHANGE THIS to reflect the path to covidm.
cm_force_rebuild = F;
cm_build_verbose = T;
cm_version = 1;
source(paste0(cm_path, "/R/covidm.R"))

# define data to fit to.
# note, these are just to illustrate the fitting process -- no empirical data here.
incidence = fread(
"date new_cases
2020-03-16 1
2020-03-18 2
2020-03-20 2
2020-03-22 4
2020-03-24 7
2020-03-26 10
2020-03-28 14
2020-03-30 21");
incidence[, date := ymd(date)];

# define likelihood function
likelihood = function(parameters, dynamics, data, x)
{
    inc = data;
    inc[, t := as.numeric(date - ymd(parameters$date0))];
    eval = merge(dynamics[compartment == "cases_reported", .(model_case = sum(value)), by = t], inc, by = "t");
    ll = sum(dpois(eval$new_cases, lambda = pmax(0.1, eval$model_case), log = T));
    
    return (ll)
}

# define parameters func, which interprets a proposal for the posterior distribution as a parameter set usable by the underlying model.
pf = function(parameters, x)
{
    x = as.list(x);

    n_groups = length(parameters$pop[[1]]$size);
    parameters$pop[[1]]$u = rep(x$u, n_groups);
    parameters$pop[[1]]$seed_times = floor(x$t_intro) + 0:6; # seed for 7 days

    return (parameters);
}

# build parameters for England
params = cm_parameters_SEI3R(cm_uk_locations("E", 1), deterministic = T);

# build priors
priors = list(
    u = "N 0.1 0.025 T 0 0.2",
    t_intro = "U 0 10"
);

# run fit.
# note, the MCMC fitter will report an acceptance rate of -1 for the first few hundred iterations.
# don't worry about this, it's only because it waits until there are 1000 samples overall to 
# start calculating the acceptance rate.
fit = cm_fit(
    base_parameters = params,
    priors = priors,
    parameters_func = pf,
    likelihood_func = likelihood,
    data = incidence,
    mcmc_burn_in = 500, mcmc_samples = 2000, mcmc_init_opt = F, opt_maxeval = 25
);

# show posterior
cm_plot_posterior(fit);
cm_plot_pairwise(fit);

# use posterior to generate sample dynamics from the model
dyn = cm_sample_fit(fit, 25)

# summarize these runs
summ = dyn[compartment == "cases_reported", .(cases_reported = sum(value)), by = .(t, run)]
summ = summ[, cm_mean_hdi(cases_reported), by = t]

# show model fit
ggplot(summ[t <= 30]) +
    geom_ribbon(aes(x = t, ymin = lower, ymax = upper), fill = "#eeeeff") +
    geom_line(aes(x = t, y = mean), colour = "#8888ff") +
    geom_point(data = incidence, aes(x = t, y = new_cases))

