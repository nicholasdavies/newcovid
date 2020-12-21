# v1/cm1_backend.R
# for building the backend

# Build C++ code
cm_source_backend = function(user_defined = NULL)
{
    # Build C++ code
    if (!cm_force_shared_) {
        sourceCpp(paste0(cm_path_, "/model_v1/corona.cpp"), 
                  rebuild  = cm_force_rebuild_,
                  cacheDir = cm_build_dir_,
                  verbose  = cm_build_verbose_)
    } else {
        # WARNING: will use old code even if newer available
        for (src in list.files(cm_build_dir_, pattern = "*\\.R", recursive = T, full.names = T)) source(src)
    }
}
