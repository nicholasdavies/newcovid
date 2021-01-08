library(deSolve)
library(ggplot2)
library(data.table)
library(stringr)
library(shmanipulate)

# Simple 2-strain model with 2 groups
strain = function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
        dS1 = - S1 * (b1 * I1 + b12 * I2) - S1 * (b1 * J1 + b12 * J2) * a
        dI1 =   S1 * (b1 * I1 + b12 * I2)                                 - g * I1 
        dJ1 =                             + S1 * (b1 * J1 + b12 * J2) * a - g * J1 
        dR1 =                                                             + g * (I1 + J1)
        dS2 = - S2 * (b2 * I2 + b12 * I1) - S2 * (b2 * J2 + b12 * J1) * a
        dI2 =   S2 * (b2 * I2 + b12 * I1)                                 - g * I2 
        dJ2 =                             + S2 * (b2 * J2 + b12 * J1) * a - g * J2 
        dR2 =                                                             + g * (I2 + J2)
        list(c(dS1, dI1, dJ1, dR1, dS2, dI2, dJ2, dR2))
    })
}

parameters = c(b1 = 0.4, b2 = 0.3, b12 = 0.2, a = 1.25, g = 0.35)
state      = c(S1 = 0.999, I1 = 0.001, J1 = 0, R1 = 0, S2 = 0.999, I2 = 0.00999, J2 = 0.00001, R2 = 0)
times      = seq(0, 250, by = 0.1)

out = ode(y = state, times = times, func = strain, parms = parameters)
plot(out)

out = as.data.table(out)

ggplot(out) +
    geom_line(aes(x = time, y = (J1 + J2) / (I1 + I2 + J1 + J2))) +
    scale_y_continuous(trans = scales::logit_trans(), breaks = c(0.01, 0.1, 0.2, 0.5, 0.7, 0.9, 0.99))



# 2-strain model with N groups
n_group = 4
strainN = function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
        
        S = state[seq(from = 1, by = 4, length.out = n_group)]
        I = state[seq(from = 2, by = 4, length.out = n_group)]
        J = state[seq(from = 3, by = 4, length.out = n_group)]
        R = state[seq(from = 4, by = 4, length.out = n_group)]
        
        dS = rep(0, 2)
        dI = rep(0, 2)
        dJ = rep(0, 2)
        dR = rep(0, 2)
        
        a = parameters[1]; # transmission advantage of strain J relative to strain I
        g = parameters[2]; # recovery rate
        b = diag(parameters[4:(4 + n_group - 1)], nrow = n_group) + (1 - diag(n_group)) * parameters[3]; # contact matrix
        
        foiI = b %*% I      # I strain FOI
        foiJ = b %*% J * a  # J strain FOI
        
        for (i in 1:n_group)
        {
            dS[i] = - S[i] * foiI[i] - S[i] * foiJ[i]
            dI[i] =   S[i] * foiI[i]                  - g * I[i] 
            dJ[i] =                  + S[i] * foiJ[i] - g * J[i]
            dR[i] =                                   + g * (I[i] + J[i])
        }
        
        list(c(rbind(dS, dI, dJ, dR)))
    })
}

make_params = function(advantage, duration, b_off, b_highest, b_lowest)
{
    c(advantage, 1 / duration, b_off, seq(b_highest, b_lowest, length.out = n_group))
}

make_initial = function(I0, J0, whereJ)
{
    whereJ = unique(pmin(whereJ, n_group))
    
    S = rep(1/n_group, n_group);
    I = S * I0;
    J = rep(0, n_group);
    J[whereJ] = S[whereJ] * J0;
    S = S - I - J;
    R = rep(0, n_group);
    
    state = c(rbind(S, I, J, R))
    names(state) = paste0(rep(c("S", "I", "J", "R"), n_group), rep(1:n_group, each = 4))
    
    return (state)
}

shmanipulate({
    n_group <<- n_groups
    
    parameters = make_params(advantage, duration, b_off, b_highest, b_lowest)
    state      = make_initial(I0, J0, start_group)
    times      = seq(0, tmax, by = 0.1)

    out = ode(y = state, times = times, func = strainN, parms = parameters)
    out = as.data.table(out)

    melted = melt(out, id.vars = "time")
    melted[, state := str_remove_all(variable, "[0-9]*")]
    melted[, group := str_remove_all(variable, "[A-Z]*")]
    freq = merge(melted[state == "I", .(I = sum(value)), by = time], melted[state == "J", .(J = sum(value)), by = time])

    plot_states = ggplot(melted) +
        geom_line(aes(x = time, y = value, colour = group)) +
        facet_wrap(~state, scales = "free")
    
    plot_freq = ggplot(freq) +
        geom_line(aes(x = time, y = J / (I + J))) +
        scale_y_continuous(trans = scales::logit_trans(), breaks = c(0.01, 0.1, 0.2, 0.3, 0.5, 0.7, 0.8, 0.9, 0.99), limits = c(0.01, 0.99))
    
    cowplot::plot_grid(plot_states, plot_freq, nrow = 1)
}, n_groups = c(4, 1, 8, 1), advantage = c(1.25, 1, 2), duration = c(6, 1, 10), b_off = c(0.1, 0, 1), b_highest = c(0.8, 0, 1), b_lowest = c(0.4, 0, 1),
    start_group = c(4, 1, 8, 1), I0 = c(0.01, 0, 0.01), J0 = c(0.0001, 0, 0.001), tmax = 250, options = list(ncol = 3))

