# 2-interventions.R

# covidm options
cm_path = "~/Dropbox/nCoV/covidm/"; ### CHANGE THIS to reflect the path to covidm.
cm_force_rebuild = F;
cm_build_verbose = T;
cm_version = 1;
source(paste0(cm_path, "/R/covidm.R"))

# build parameters for all of England, regional level (9 regions)
params = cm_parameters_SEI3R(cm_uk_locations("E", 2), deterministic = T, date_start = "2020-03-01", date_end = "2021-03-01");

# introduce some simple interventions. 
# These are school terms in England; school terms in Wales, Scotland, and Northern Ireland differ slightly.
# (updated on 2nd April to start in February; note this only looks ahead to 1st September 2021)
school_close  = c("2020-2-16", "2020-4-05", "2020-5-24", "2020-7-22", "2020-10-25", "2020-12-20", "2021-02-14", "2021-04-01", "2021-05-30", "2021-07-25");
school_reopen = c("2020-2-22", "2020-4-18", "2020-5-30", "2020-9-01", "2020-10-31", "2021-01-02", "2021-02-20", "2021-04-17", "2021-06-05", "2021-09-01");

iv = cm_iv_build(params) # this sets up a data structure for doing interventions
cm_iv_school_breaks(iv, school_close, school_reopen) # populates iv with school terms -- uses cm_iv_set under the hood.
cm_iv_contact(iv, "2020-05-01", "2020-06-18", c(0.1, 0.1, 1, 0.1)) # changes contact matrix components. For May 1st to June 18th inclusive, 
# home, work, and other contacts 90% less. By default the matrix components are home, work, school, and other.
cm_iv_set(iv, "2020-03-01", "2020-04-01", travel = 0.1 + 0.9 * diag(9), fun = "assign")
params = cm_iv_apply(params, iv) # sets the "schedule" parameter to follow interventions in iv.

# Set seeds to control start of outbreak
params$pop[[1]]$seed_times = rep(0:6, each = 5) # 5 new infections each day for first 7 days
params$pop[[1]]$dist_seed_ages = cm_age_coefficients(20, 50, 5 * (0:length(params$pop[[1]]$size))) # infections start in individuals aged 20-50

# One useful thing to do is to set parameters to achieve a desired R0.
# You can calculate the starting R0 for a population (i.e. before any scheduled interventions take effect)
# using cm_calc_R0.
current_R0 = cm_calc_R0(params, 1); # calculate R0 in population 1 of params

# This is a little high, we want to try R0 = 2.4. To do this we can reduce u (susceptibility) by the appropriate amount.
target_R0 = 2.4
params$pop[[1]]$u = params$pop[[1]]$u * target_R0 / current_R0

# check to see it's worked...
cm_calc_R0(params, 1)

# Finally - model is deterministic by default, but can also be stochastic.
params$deterministic = F

# run the model
run = cm_simulate(params, 1, 1515) # second parameter is number of runs w same parameters but different seed. third parameter is starting seed. 0 is default

# show results
ggplot(run$dynamics[compartment == "cases"]) +
    geom_line(aes(ymd("2020-03-01") + t, value, colour = group, group = group)) + 
    facet_wrap(~population)
