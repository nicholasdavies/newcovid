# 1-getting-started.R
# for covidm version 2

# covidm options
cm_path = "~/Dropbox/nCoV/covidm/"; ### CHANGE THIS to reflect the path to covidm.
cm_force_rebuild = F;
cm_build_verbose = T;
cm_version = 2;
source(paste0(cm_path, "/R/covidm.R"))

# build parameters for all of UK, down to the regional level (level 2).
params = cm_parameters_SEI3R(cm_uk_locations("UK", 2), deterministic = T);

# # alternatively: build parameters for another country.
# params = cm_parameters_SEI3R("Italy");

# run the model
run = cm_simulate(params, 1)

# show results
ggplot(run$dynamics[compartment == "cases"]) +
    geom_line(aes(t, value, colour = group, group = group)) +
    facet_wrap(~population)

