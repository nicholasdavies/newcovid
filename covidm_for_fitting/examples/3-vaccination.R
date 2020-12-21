# 2-interventions.R

# covidm options
cm_path = "~/Dropbox/nCoV/covidm/"; ### CHANGE THIS to reflect the path to covidm.
cm_force_rebuild = F;
cm_build_verbose = T;
cm_version = 2;
source(paste0(cm_path, "/R/covidm.R"))

# build parameters to start from
params = cm_parameters_SEI3R("Peru", deterministic = T, date_start = "2020-03-01", date_end = "2021-03-01", u = 0.03);
n_age_groups = length(params$pop[[1]]$size);
cm_calc_R0(params, 1)

# There are 3 parameters for controlling vaccination.
# v -- number of vaccines administered per day, for each age group.
# ev -- vaccine effectiveness, assuming all-or-nothing protection. (So effectively, v[i] * ev[i] people get vaccine protection each day in age group i.)
# wv -- rate at which vaccine protection wanes (per day), for each age group.

params$pop[[1]]$v = rep(0, n_age_groups);        # no vaccines administered
params$pop[[1]]$ev = rep(0.8, n_age_groups);     # 80% efficacy
params$pop[[1]]$wv = rep(1 / 365, n_age_groups); # vaccine protection lasts for 1 year on average
params$schedule = list();                        # no scheduled changes to parameters

# So let's contrast three scenarios:
# (1) no vaccination

run1 = cm_simulate(params)

# (2) 10,000 vaccines per day in each age group

params$pop[[1]]$v = rep(10000, n_age_groups);
run2 = cm_simulate(params)

# (3) 100,000 vaccines per day in each age group, just in under-50s, and only for the first 20 days of the simulation

params$pop[[1]]$v = rep(0, n_age_groups);
params$schedule = list(
    list(
        parameter = "v",
        pops = 0,
        mode = "assign",
        values = list(rep(c(100000, 0), times = c(10, n_age_groups - 10)), rep(0, n_age_groups)),
        times = c(0, 20)
    )
);
run3 = cm_simulate(params)

# show results

ggplot(run1$dynamics[compartment == "cases"]) +
    geom_line(aes(ymd("2020-03-01") + t, value, colour = group, group = group)) + 
    facet_wrap(~population)

ggplot(run2$dynamics[compartment == "cases"]) +
    geom_line(aes(ymd("2020-03-01") + t, value, colour = group, group = group)) + 
    facet_wrap(~population)

ggplot(run3$dynamics[compartment == "cases"]) +
    geom_line(aes(ymd("2020-03-01") + t, value, colour = group, group = group)) + 
    facet_wrap(~population)

# Total number of cases in each scenario
run1$dynamics[compartment == "cases", sum(value)]
run2$dynamics[compartment == "cases", sum(value)]
run3$dynamics[compartment == "cases", sum(value)]
run1$dynamics[compartment == "cases", sum(value), by = group]
run2$dynamics[compartment == "cases", sum(value), by = group]
run3$dynamics[compartment == "cases", sum(value), by = group]
