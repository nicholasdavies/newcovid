library(data.table)
library(qs)
library(ggplot2)
library(cowplot)

nhs_regions = c("East of England", "England", "London", "Midlands", "North East and Yorkshire",
    "North West", "Northern Ireland", "Scotland", "South East", "South West", "United Kingdom", "Wales")

# detect outlier chains
outliers = function(post_list, pops)
{
    r = NULL;
    for (p in pops) {
        vars = names(post_list[[p]])[-1:-4]
        post = post_list[[p]]
        for (v in vars) {
            mu = post[, mean(get(v))]
            sd = post[, sd(get(v))]
            
            ind = post[, mean(get(v)), by = chain][, .(chain, div = abs(V1 - mean(V1)), indicator = abs(V1 - mu) > 1 * sd)]
            if (ind[, any(indicator == TRUE)]) {
                r = rbind(r, 
                    data.table(pop = p, variable = v, chains = ind[indicator == TRUE, chain])
                )
            }
        }
    }
    return (r)
}

# remove outliers and burn-in
strip_burn = function(post_list, pops, out, burn_in)
{
    for (p in pops)
    {
        post_list[[p]][, trial := trial - burn_in]
        if (is.null(out)) {
            outchains = numeric()
        } else {
            outchains = out[pop == p, unique(chains)]
        }
        post_list[[p]] = post_list[[p]][trial >= 0 & !chain %in% outchains]
    }
    post_list
}

# gelman-rubin diagnostic, Rhat
gelman_rubin = function(post, varname)
{
    # number of chains
    M = post[, uniqueN(chain)]
    
    # number of iterations
    N = post[, max(trial + 1), by = chain][, min(V1)]
    
    # posterior mean for each chain
    thm = post[, mean(get(varname)), by = chain][, V1]
    
    # posterior variance for each chain
    s2m = post[, var(get(varname)), by = chain][, V1]
    
    # grand posterior mean
    th = mean(thm)
    
    # how the individual means vary around the grand mean
    B_over_n = var(thm)
    
    # averaged variance of the chains
    W = mean(s2m)
    
    Vh = (N - 1) * W / N + (M + 1) * B_over_n / M;
    
    sqrt(Vh/W)
}

# Calculate all Rhats for a model fit
gr_all = function(post_list, pops)
{
    r = NULL
    for (p in pops) {
        vars = names(post_list[[p]])[-1:-4]
        for (v in vars) {
            g = gelman_rubin(post_list[[p]], v);
            r = rbind(r, data.table(pop = p, variable = v, Rhat = g))
        }
    }
    r
}

# i.e. trace plot for posteriors
caterpillar_plot = function(post_list, pops)
{
    plots = list()
    for (p in pops) {
        pm = melt(post_list[[p]], id.vars = c("trial", "chain"))
        plots[[length(plots) + 1]] = ggplot(pm) +
            geom_line(aes(x = trial, y = value, group = chain)) +
            facet_wrap(~variable, scales = "free_y", ncol = 2) +
            theme_cowplot(font_size = 8) + labs(y = NULL, title = nhs_regions[p])
    }
    plot_grid(plotlist = plots, nrow = 1)
}

# Generate output
fits = list(
    # list(file = "./fits/relu_4p22.qs", id = "relu", name = "Increased transmissibility", pops = c(1, 3, 4, 5, 6, 9, 10)),
    # list(file = "./fits/infdur_ELSE13.qs", id = "infdur", name = "Duration of infectiousness", pops = c(1, 3, 9)),
    # list(file = "./fits/immesc_ELSE13.qs", id = "immesc", name = "Immune escape", pops = c(1, 3, 9)),
    # list(file = "./fits/ch_u_ELSE12.qs", id = "ch_u", name = "Increased susceptibility in children", pops = c(1, 3, 9)),
    # list(file = "./fits/serial_3p12.qs", id = "serial", name = "Shorter generation time", pops = c(1, 3, 9))
    list(file = "./fits/combined_3p11.qs", id = "combined", name = "Combined model", pops = c(1, 3, 9))
)

#Rhat = NULL
for (f in fits) {
    cat("Loading", f$id, "...\n");
    w = qread(f$file)
    post = w[[1]]
    cat("Detecting outlier chains...\n");
    out = outliers(post, f$pops)
    cat("Stripping", nrow(out[, unique(chains), by = pop]), "chains...\n");
    post = strip_burn(post, f$pops, out, if (f$id == "relu") 0 else 1000)
    cat("Calculating convergence diagnostics...\n");
    Rhat = rbind(Rhat,
        cbind(model = f$name, gr_all(post, f$pops))
    )
    cat("Building trace plot...\n");
    cp = caterpillar_plot(post, f$pops)
    cat("Saving trace plot...\n");
    ggsave(paste0("./output/trace3_", f$id, ".png"), cp, width = 20 * length(f$pops) / 3, height = 30, units = "cm")
    cat("Saving fit...\n");
    w[[1]] = post
    qsave(w, paste0("./fits/final/", f$id, ".qs"));
}

Rhat[Rhat > 1.1]
fwrite(dcast(Rhat, variable ~ model + pop), "./output/Rhat.csv")
