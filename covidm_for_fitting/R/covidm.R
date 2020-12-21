# covidm.R
# main file to source for users of the covidm API.

# TODO
#  - recover better from user halting of execution of cm_fit
#  - recover better from crashes during cm_fit
# TODO
#  - make sure seq_len and seq_along are used instead of 1:length() for 0-length things

packageStartupMessage("Loading covidm.")

# Required libraries
suppressPackageStartupMessages({
    library(data.table)   # for data.table, an enhanced (and faster) data.frame
    library(ggplot2)      # for plotting
    library(Rcpp)         # for running the C++ model backend
    library(RcppGSL)      # for linking with the GNU Scientific Library for the C++ backend
    library(qs)           # for qsave and qread, faster equivalents of saveRDS and readRDS
    library(lubridate)    # for manipulating dates and times. NB requires stringr
    library(HDInterval)   # for summarizing results
    library(cowplot)      # for plotting grids
})

# Settings for covidm
cm_option_ = function(name, default_value) {
    if (!exists(name)) {
        packageStartupMessage("Using default value for ", name, ": ", default_value)
    }
    return (get0(name, ifnotfound = default_value))
}

cm_path_           = cm_option_("cm_path", "~/Dropbox/nCoV/covidm/")
cm_version_        = cm_option_("cm_version", 2)
cm_build_          = cm_option_("cm_build", T)
cm_force_rebuild_  = cm_option_("cm_force_rebuild", F)
cm_build_dir_      = cm_option_("cm_build_dir", paste0(cm_path_, "/build/"))
cm_build_verbose_  = cm_option_("cm_build_verbose", T)
cm_force_shared_   = cm_option_("cm_force_shared", F)

# Attach code
source(paste0(cm_path_, "/R/shared/cmS_misc.R"))
source(paste0(cm_path_, "/R/shared/cmS_plot.R"))

if (cm_version_ == 1) {
    source(paste0(cm_path_, "/R/v1/cm1_backend.R"))
    source(paste0(cm_path_, "/R/v1/cm1_run.R"))
    source(paste0(cm_path_, "/R/v1/cm1_params.R"))
    source(paste0(cm_path_, "/R/v1/cm1_interventions.R"))
    source(paste0(cm_path_, "/R/v1/cm1_fit.R"))
} else if (cm_version_ == 2) {
    source(paste0(cm_path_, "/R/v2/cm2_backend.R"))
    source(paste0(cm_path_, "/R/v2/cm2_run.R"))
    source(paste0(cm_path_, "/R/v2/cm2_params.R"))
    source(paste0(cm_path_, "/R/v2/cm2_interventions.R"))
    source(paste0(cm_path_, "/R/v2/cm2_fit.R"))
} else {
    stop(paste("Requested covidm version", cm_version_, "but there are only versions 1 and 2."));
}

if (cm_build_) {
    packageStartupMessage("Building backend code...");
    cm_source_backend();
}

# Attach data
packageStartupMessage("Loading data...");
cm_matrices     = readRDS(paste0(cm_path_, "/data/all_matrices.rds"));
cm_populations  = readRDS(paste0(cm_path_, "/data/wpp2019_pop2020.rds"));
cm_structure_UK = readRDS(paste0(cm_path_, "/data/structure_UK.rds"));
cm_highrisk     = readRDS(paste0(cm_path_, "/data/prevalence_morbidities.rds"))

