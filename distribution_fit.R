skewness = function(x)
{
    mean((x - mean(x))^3) / (mean((x - mean(x))^2)^1.5)
}

skew_normal_objective = function(x, target)
{
    central = target[1]
    lower = target[2]
    upper = target[3]
    
    xi = x[1]
    omega = x[2]
    alpha = x[3]

    delta = alpha / sqrt(1 + alpha^2)
    
    mean = xi + omega * delta * sqrt(2 / pi)
    lq = sn::psn(lower, xi, omega, alpha)
    uq = sn::psn(upper, xi, omega, alpha)
    
    return ((central - mean) ^ 2 + (lq - 0.025) ^ 2 + (uq - 0.975) ^ 2)
}

skew_normal_solve = function(central, lower, upper, adj_react2)
{
    d = NULL;
    set.seed(12345)
    for (i in seq_along(central))
    {
        if (adj_react2[i]) {
            # Get point estimate of q, crude seroprevalence (non-adjusted)
            p = central[i]
            sensitivity = 0.844
            specificity = 0.986
            q = p * (sensitivity) + (p - 1) * (specificity - 1)
            
            # Model p with uncertainty over sensitivity and specificity
            sensitivity = rbeta(100000, shape1 = 31.65, shape2 = 6.3)
            specificity = rbeta(100000, 4*98.6, 4*1.4)[sensitivity > q]
            sensitivity = sensitivity[sensitivity > q]
            p = (q + specificity - 1) / (sensitivity + specificity - 1);
            
            # Get skew-normal distribution that matches these parameters
            par = cp2dp(c(mean(p), sd(p), min(0.98, skewness(p))), "SN");
            d = rbind(d,
                data.table(xi = par[1], omega = par[2], alpha = par[3]));
        } else {
            solution = optim(par = c(central[i], (upper[i] - lower[i]) / 4, 0), fn = skew_normal_objective, 
                target = c(central[i], lower[i], upper[i]), 
                method = "Nelder-Mead", control = list(maxit = 10000, abstol = 0, reltol = 0));
            d = rbind(d, 
                data.table(xi = solution$par[1], omega = solution$par[2], alpha = solution$par[3])
            );
        }
    }
    
    return (d)
}
