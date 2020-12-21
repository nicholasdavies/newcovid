# 5-observer.R

# covidm options
cm_path = "~/Dropbox/nCoV/covidm/"; ### CHANGE THIS to reflect the path to covidm.
cm_force_rebuild = F;
cm_build_verbose = T;
cm_version = 1;
source(paste0(cm_path, "/R/covidm.R"))

# build parameters for Italy
params = cm_parameters_SEI3R("Italy", deterministic = T);

# Add an observer. The observer function is called each day.
params$pop[[1]]$observer = function(time, dynamics)
{
    # The dynamics data.table contains the same information as the dynamics returned
    # in cm_simulate, but it is formatted differently. There are no "run" or "compartments"
    # columns; instead, compartments are each their own separate column. So dynamics contains
    # columns t, population, group, S, E, Ia, and so on. Also, since the dynamics have not
    # been "annotated" at this point, populations and groups are integers starting at 1, rather
    # than factors with the specified names of populations and groups. This may change in
    # future releases.
    
    # We can use dynamics to, for example, calculate the number of cases reported today.
    cases_today = dynamics[t == time, sum(cases_reported)];

    # The observer needs to return either nothing (NULL) or a list of actions to do.
    actions = list()
    
    # Actions can do four things. First, they can change parameters in the same way as schedule.
    if (cases_today > 25000) {
        actions$changes = list(contact = c(0.1, 0.1, 0.1, 0.1));
    }
    
    # They can print a message to R console output. Note, this needs to be a character vector.
    if (time == 75) {
        actions$print = paste("Cases on day 75:", cases_today);
    }
    
    # They can halt the simulation early with halt = T.
    if (time == 90) {
        actions$print = "Halting early.";
        actions$halt = T;
    }
    
    # Finally, they can save lines to an internal buffer. This is called csv because the intention
    # is for you to format these as comma-separated values, but you can put any character vector here.
    # A newline is added to the end if there isn't one there.
    actions$csv = paste(time, dynamics[t == time & group == 1, cases_reported], sep = ",");
    
    return (actions)
}

# run the model
run = cm_simulate(params, 1)

# get csv output
under5 = fread(run$csv[[1]])
plot(under5)

# show full results
ggplot(run$dynamics[compartment == "cases", .(total = sum(value)), by = .(population, t)]) +
    geom_line(aes(t, total)) + 
    facet_wrap(~population)
