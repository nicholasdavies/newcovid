load_fit = function(filename)
{
    saved = qread(filename)
    posteriorsI = saved[[1]]
    parametersI = saved[[2]]
    rm(saved)
    
    # Ensure newly added parameters are in parametersI
    for (p in seq_along(parametersI))
    {
        if (!is.null(parametersI[[p]])) {
            if (is.null(parametersI[[p]]$pop[[1]]$ed_vi)) {
                parametersI[[p]]$pop[[1]]$ed_vi = rep(0, 16);
            }
            if (is.null(parametersI[[p]]$pop[[1]]$ed_vi2)) {
                parametersI[[p]]$pop[[1]]$ed_vi2 = rep(0, 16);
            }
        }
    }
    
    # Remove leftover columns from posteriorsI from earlier fitting
    for (p in seq_along(posteriorsI))
    {
        if (!is.null(posteriorsI[[p]])) {
            if (!is.null(posteriorsI[[p]]$v2_ab)) {
                posteriorsI[[p]][, v2_ab := NULL];
            }
        }
    }
    
    posteriorsI <<- posteriorsI
    parametersI <<- parametersI
}

load_fit("./fits/pp10.qs")

