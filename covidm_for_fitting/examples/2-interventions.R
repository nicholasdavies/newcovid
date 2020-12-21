# 2-interventions.R

# covidm options
cm_path = "~/Dropbox/nCoV/covidm/"; ### CHANGE THIS to reflect the path to covidm.
cm_force_rebuild = F;
cm_build_verbose = T;
cm_version = 2;
source(paste0(cm_path, "/R/covidm.R"))

# build parameters for all of England, regional level (9 regions)
params = cm_parameters_SEI3R(cm_uk_locations("E", 2), deterministic = T, date_start = "2020-03-01", date_end = "2021-03-01");

# Introduce some simple interventions. 

# These are school holidays in England; school holidays in Wales, Scotland, and Northern Ireland differ slightly.
school_holidays  = c(
  "2020-02-16", "2020-02-23",
  "2020-04-05", "2020-04-19",
  "2020-05-24", "2020-05-31",
  "2020-07-22", "2020-09-02",
  "2020-10-25", "2020-11-01",
  "2020-12-20", "2021-01-03",
  "2021-02-14", "2021-02-21",
  "2021-04-01", "2021-04-18",
  "2021-05-30", "2021-06-06",
  "2021-07-25", "2021-09-02");

# To schedule changes to parameters in covidm2, you can use the schedule parameter.
# It's a list of lists.
params$schedule = list()

# Each inner list has 5 components, as illustrated here in setting school holidays:
params$schedule[[1]] = list(
    parameter = "contact",   # Which parameter to change
    pops = numeric(),        # Which populations to apply the parameter change to (0 to N-1); or leave empty for all populations
    mode = "assign",         # How parameter changes interact with the existing parameter. Here, assign means just overwrite whatever is already there with the new values.
    values = rep(list(c(1, 1, 0, 1), c(1, 1, 1, 1)), 10), # Turn school contacts off and then on 10 times. By default the matrix components are home, work, school, and other.
    times = school_holidays  # Times at which to change parameter
);

# Let's also set home, work, and other contacts to 90% less from May 1st to June 18th.
# We use mode = "multiply" here to sit well with the other changes to the contact parameter.
params$schedule[[2]] = list(
    parameter = "contact",
    pops = numeric(),
    mode = "multiply",
    values = list(c(0.1, 0.1, 1, 0.1), c(1, 1, 1, 1)),
    times = c("2020-05-01", "2020-06-19")
);

# Also, set the travel matrix from March to April.
params$schedule[[3]] = list(
    parameter = "travel",
    pops = numeric(),
    mode = "assign",
    values = list(0.1 + 0.9 * diag(9), numeric(0)),
    times = c("2020-03-01", "2020-04-02")
);

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
