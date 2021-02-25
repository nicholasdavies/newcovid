# Calculate Tg for transmission model
#
#

rgamma2 = function(n, mean, shape)
{
    rate = shape / mean
    rgamma(n, shape = shape, rate = rate)
}

ns = 100000
symp = data.table(start = rgamma2(ns, 2.5, 2.5))
symp[, end := start + rgamma2(ns, 2.5, 4.0) + rgamma2(ns, 2.5, 4.0)]
symp[, n := rpois(ns, end - start)]
symp[, Tg := start + (end - start) / 2]
symp[, weighted.mean(Tg, n)] # 5.3

rows = symp[, sample(.N, ns, replace = TRUE, prob = n)]
sd(symp[rows, runif(.N, start, end)]) # 2.5

asymp = data.table(start = rgamma2(ns, 2.5, 2.5))
asymp[, end := start + rgamma2(ns, 5.0, 4.0)]
asymp[, n := rpois(ns, end - start)]
asymp[, Tg := start + (end - start) / 2]
asymp[, weighted.mean(Tg, n)] # 5.6

rows = asymp[, sample(.N, ns, replace = TRUE, prob = n)]
sd(asymp[rows, runif(.N, start, end)]) # 2.9

pc = 1/3
ps = 1 - pc

(pc * 5.3 + 0.5 * ps * 5.6) / (pc + 0.5 * ps)
(pc * 2.5 + 0.5 * ps * 2.9) / (pc + 0.5 * ps)

    # dE  = cm_delay_gamma(2.5, 2.5, t_max = 15, t_step = 0.25)$p,
    # dIp = cm_delay_gamma(2.5, 4.0, t_max = 15, t_step = 0.25)$p,
    # dIs = cm_delay_gamma(2.5, 4.0, t_max = 15, t_step = 0.25)$p,
    # dIa = cm_delay_gamma(5.0, 4.0, t_max = 15, t_step = 0.25)$p)


Tg = 5.5 # 5.5
sigma = 2.7 # 2.1
k = (sigma/Tg)^2

toRt = function(r)
{
    (1 + k * r * Tg) ^ (1/k)
}

tor = function(Rt)
{
    (Rt^k - 1) / (k * Tg)
}
