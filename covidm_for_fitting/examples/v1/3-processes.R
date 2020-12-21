# 1-getting-started.R

# covidm options
cm_path = "~/Dropbox/nCoV/covidm/"; ### CHANGE THIS to reflect the path to covidm.
cm_force_rebuild = F;
cm_build_verbose = T;
cm_version = 1;
source(paste0(cm_path, "/R/covidm.R"))

# Build parameters for a country 
params = cm_parameters_SEI3R("Spain");

# Set an observation process
# TODO explain this better and make relevant
params$processes = list(
    list(source = "S", type = "multinomial", names = c("test1", "test2"), report = c("ipo", "ipo"),
        prob = matrix(c(rep(0.1, 16), rep(0.9, 16)), nrow = 2, ncol = 16, byrow = T),
        delays = matrix(c(cm_delay_gamma(7, 7, 60, 0.25)$p, cm_delay_gamma(14, 14, 60, 0.25)$p), nrow = 2, byrow = T)
    )
)

# what the above does.
# this says to set up an observation process. Individuals enter the process when they *leave* compartment [source].
# So e.g. if you want to track symptomatics (e.g. to track burdens) you should set source = "Is".
# When individuals enter this process, they get split (that's the type = "multinomial" bit -- only thing implemented so far) into
# two subprocesses, test1 and test2. prob is a [number of subprocesses] x [number of age groups] matrix with age-specific probabilities of
# entering either the "test1" or "test2" subprocess. These are unnormalised weights actually, so you can't just leave one blank at the moment.
# report = "ipo" means report the incidence, prevalence, and "outcidence" of individuals in the subprocesses. Can just specify e.g. "p" or "ip" or whatever.
# delays is how long individuals stay in the subprocesses and is not split by age. So far I just have cm_delay_gamma (for a gamma distribution,
# with parameters mean, shape, max_t and t_step, which needs to be the same as the simulation's time step).
# You can chain these together by setting the source of another process to one of the names of the subprocesses.

# run the model
run = cm_simulate(params)

# show results
ggplot(run$dynamics[compartment %like% "test"]) +
    geom_line(aes(t, value, colour = group, group = group)) + 
    facet_wrap(~compartment)
