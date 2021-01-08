# Probability of ICU given hospitalisation (derived from CO-CIN)
picu_cocin_func = function(age)
{
    model = qread(datapath("icu_model.qs"))
    pred = data.table(patient_age = age)
    p = predict(model, pred)
    exp(p) / (1 + exp(p))
}
picu_cocin = picu_cocin_func(0:85)

# Infection fatality rate (derived from Levin et al., preprint)
ifr_levin = 100 * exp(-7.56 + 0.121 * 0:85) / (100 + exp(-7.56 + 0.121 * 0:85)) / 100;

# Infection hospitalisation rate (derived from Salje et al., Science)
ihr_salje = exp(-7.37 + 0.068 * 0:85) / (1 + exp(-7.37 + 0.068 * 0:85));

# Amalgamate probabilities
probabilities = data.table(age = 0:85, ihr = ihr_salje, ifr = ifr_levin, picu = picu_cocin)
probabilities[, age_group := pmin(15, age %/% 5)]
probabilities = probabilities[, lapply(.SD, mean), by = age_group, .SDcols = 2:4]

# Create model burden processes
P.hosp     = probabilities[, ihr];
P.critical = probabilities[, ihr * picu];
P.severe   = probabilities[, ihr * (1 - picu)];
P.death    = probabilities[, ifr];

delay_normal = function(mu, sd, t_max, t_step)
{
    t_points = seq(0, t_max, by = t_step);
    heights = pnorm(t_points + t_step/2, mu, sd) -
        pnorm(pmax(0, t_points - t_step/2), mu, sd);
    return (data.table(t = t_points, p = heights / sum(heights)))
}


burden_processes = list(
    # Strain 1 burden processes
    list(source = "newI", type = "multinomial", names = c("to_hosp", "null"), report = c("", ""),
        prob = matrix(c(P.hosp, 1 - P.hosp), nrow = 2, ncol = 16, byrow = T),
        delays = matrix(c(cm_delay_gamma(6.0 + 2.5, 0.71, 60, 0.25)$p, cm_delay_skip(60, 0.25)$p), nrow = 2, byrow = T)),

    list(source = "newI", type = "multinomial", names = c("to_icu", "null"), report = c("", ""),
        prob = matrix(c(P.critical, 1 - P.critical), nrow = 2, ncol = 16, byrow = T),
        delays = matrix(c(cm_delay_gamma(9.6 + 2.5, 1.91, 60, 0.25)$p, cm_delay_skip(60, 0.25)$p), nrow = 2, byrow = T)),

    list(source = "to_hosp", type = "multinomial", names = "hosp", report = "ip",
        prob = matrix(1, nrow = 1, ncol = 16, byrow = T),
        delays = matrix(cm_delay_lnorm(11.08, 1.202, 60, 0.25)$p, nrow = 1, byrow = T)),

    list(source = "to_icu", type = "multinomial", names = "icu", report = "ip",
        prob = matrix(1, nrow = 1, ncol = 16, byrow = T),
        delays = matrix(cm_delay_lnorm(13.33, 1.25, 60, 0.25)$p, nrow = 1, byrow = T)),

    list(source = "newI", type = "multinomial", names = c("death", "null"), report = c("o", ""),
        prob = matrix(c(P.death, 1 - P.death), nrow = 2, ncol = 16, byrow = T),
        delays = matrix(c(cm_delay_lnorm(15, 0.9, 60, 0.25)$p, cm_delay_skip(60, 0.25)$p), nrow = 2, byrow = T)),

    # Sero and PCR positivity
    # note -- delay should actually be gamma (symptom onset) plus normal -- to be fixed.
    # from Borremans et al
    list(source = "newII2", type = "multinomial", names = "to_lfia_positive", report = c(""),
        prob = matrix(1, nrow = 1, ncol = 16, byrow = T),
        delays = matrix(delay_normal(11.9 + 2.5, 5.3, 60, 0.25)$p, nrow = 1, byrow = T)),

    list(source = "to_lfia_positive", type = "multinomial", names = "lfia_positive", report = c("ip"),
        prob = matrix(1, nrow = 1, ncol = 16, byrow = T),
        delays = matrix(cm_delay_gamma(224, 1, 730, 0.25)$p, nrow = 1, byrow = T)),
    
    list(source = "to_hosp", type = "multinomial", names = "hosp_undetected", report = c("po"),
        prob = matrix(1, nrow = 1, ncol = 16, byrow = T),
        delays = matrix(cm_delay_gamma(6, 1, 60, 0.25)$p, nrow = 1, byrow = T)),
    
    list(source = "newEE2", type = "multinomial", names = "to_pcr_positive", report = c(""),
        prob = matrix(1, nrow = 1, ncol = 16, byrow = T),
        delays = matrix(c(cm_delay_gamma(2.76, 4.79, 60, 0.25)$p), nrow = 1, byrow = T)),

    list(source = "to_pcr_positive", type = "multinomial", names = "pcr_positive", report = c("ip"),
        prob = matrix(1, nrow = 1, ncol = 16, byrow = T),
        delays = matrix(cm_delay_gamma(8.47, 1.96, 730, 0.25)$p, nrow = 1, byrow = T)),
    
    # Strain 2 burden processes
    list(source = "newI2", type = "multinomial", names = c("to_hosp2", "null"), report = c("", ""),
        prob = matrix(c(P.hosp, 1 - P.hosp), nrow = 2, ncol = 16, byrow = T),
        delays = matrix(c(cm_delay_gamma(6.0 + 2.5, 0.71, 60, 0.25)$p, cm_delay_skip(60, 0.25)$p), nrow = 2, byrow = T)),

    list(source = "newI2", type = "multinomial", names = c("to_icu2", "null"), report = c("", ""),
        prob = matrix(c(P.critical, 1 - P.critical), nrow = 2, ncol = 16, byrow = T),
        delays = matrix(c(cm_delay_gamma(9.6 + 2.5, 1.91, 60, 0.25)$p, cm_delay_skip(60, 0.25)$p), nrow = 2, byrow = T)),

    list(source = "to_hosp2", type = "multinomial", names = "hosp2", report = "ip",
        prob = matrix(1, nrow = 1, ncol = 16, byrow = T),
        delays = matrix(cm_delay_lnorm(11.08, 1.202, 60, 0.25)$p, nrow = 1, byrow = T)),

    list(source = "to_icu2", type = "multinomial", names = "icu2", report = "ip",
        prob = matrix(1, nrow = 1, ncol = 16, byrow = T),
        delays = matrix(cm_delay_lnorm(13.33, 1.25, 60, 0.25)$p, nrow = 1, byrow = T)),

    list(source = "newI2", type = "multinomial", names = c("death2", "null"), report = c("o", ""),
        prob = matrix(c(P.death, 1 - P.death), nrow = 2, ncol = 16, byrow = T),
        delays = matrix(c(cm_delay_lnorm(15, 0.9, 60, 0.25)$p, cm_delay_skip(60, 0.25)$p), nrow = 2, byrow = T)),
    
    list(source = "to_hosp2", type = "multinomial", names = "hosp_undetected2", report = c("po"),
        prob = matrix(1, nrow = 1, ncol = 16, byrow = T),
        delays = matrix(cm_delay_gamma(6, 1, 60, 0.25)$p, nrow = 1, byrow = T)),
    
    # Pillar 2
    list(source = "newI", type = "multinomial", names = "test", report = c("o"),
        prob = matrix(1, nrow = 1, ncol = 16, byrow = T),
        delays = matrix(cm_delay_gamma(2.5 + 3.9, 2, 60, 0.25)$p, nrow = 1, byrow = T)),

    list(source = "newI2", type = "multinomial", names = "test2", report = c("o"),
        prob = matrix(1, nrow = 1, ncol = 16, byrow = T),
        delays = matrix(cm_delay_gamma(2.5 + 3.9, 2, 60, 0.25)$p, nrow = 1, byrow = T))
)

