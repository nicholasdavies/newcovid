# For fitting

# Give names to elements of x from priors
named_params = function(priors)
{
    paste(mapply(function(nm, i) paste0("double x_", nm, " = x[", i - 1, "]; (void) x_", nm, ";"), names(priors), seq_along(priors)), collapse = "\n");
}

# initialiser list with values from x
cpp_vec = function(x) 
{
    paste("{", paste(x, collapse = ", "), "}")
}

asc = function(x, y0, y1, s0, s1)
{
    xx = s0 + x * (s1 - s0);
    h0 = exp(s0) / (1 + exp(s0));
    h1 = exp(s1) / (1 + exp(s1));
    h = (exp(xx) / (1 + exp(xx)) - h0) / (h1 - h0);
    return (y0 + (y1 - y0) * h);
}


# create cpp changes
cpp_chgI_voc = function(priors, v2, v2_relu, v2_latdur, v2_infdur, v2_immesc, v2_ch_u)
{
    glue::glue(
        named_params(priors),
        'vector<double> work_curve = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.008, 0.021, 0.033, 0.046, 0.058, 0.071, 0.083, 0.096, 0.108, 0.121, 0.133, 0.146, 0.158, 0.171, 0.183, 0.196, 0.208, 0.221, 0.233, 0.246, 0.258, 0.271, 0.283, 0.296, 0.308, 0.321, 0.334, 0.346, 0.359, 0.371, 0.384, 0.397, 0.41, 0.422, 0.435, 0.448, 0.461, 0.474, 0.487, 0.5, 0.513, 0.526, 0.539, 0.552, 0.566, 0.579, 0.592, 0.606, 0.619, 0.633, 0.646, 0.66, 0.674, 0.687, 0.701, 0.715, 0.729, 0.743, 0.757, 0.771, 0.785, 0.799, 0.813, 0.828, 0.842, 0.856, 0.87, 0.885, 0.899, 0.914, 0.928, 0.942, 0.957, 0.971, 0.986, 1, 1.014, 1.029, 1.043, 1.058, 1.072, 1.087, 1.101, 1.115, 1.13, 1.144, 1.159, 1.173, 1.188, 1.202, 1.216, 1.231, 1.245, 1.26, 1.274, 1.289, 1.303, 1.317, 1.332, 1.346, 1.361 };',
        'vector<double> other_curve = { 0.064, 0.066, 0.067, 0.068, 0.069, 0.071, 0.072, 0.073, 0.075, 0.076, 0.077, 0.078, 0.08, 0.081, 0.082, 0.084, 0.085, 0.086, 0.087, 0.089, 0.09, 0.091, 0.092, 0.094, 0.095, 0.096, 0.098, 0.099, 0.1, 0.101, 0.103, 0.104, 0.105, 0.106, 0.108, 0.109, 0.11, 0.112, 0.113, 0.114, 0.116, 0.118, 0.119, 0.121, 0.123, 0.125, 0.128, 0.13, 0.132, 0.135, 0.137, 0.14, 0.143, 0.146, 0.15, 0.154, 0.159, 0.164, 0.169, 0.175, 0.182, 0.19, 0.198, 0.207, 0.217, 0.228, 0.24, 0.252, 0.266, 0.28, 0.295, 0.31, 0.327, 0.344, 0.361, 0.379, 0.398, 0.418, 0.438, 0.459, 0.48, 0.502, 0.525, 0.549, 0.572, 0.597, 0.621, 0.647, 0.672, 0.698, 0.725, 0.751, 0.778, 0.805, 0.833, 0.86, 0.888, 0.916, 0.944, 0.972, 1, 1.028, 1.056, 1.084, 1.112, 1.14, 1.168, 1.196, 1.224, 1.252, 1.28, 1.308, 1.337, 1.365, 1.393, 1.421, 1.449, 1.477, 1.505, 1.533, 1.561, 1.589, 1.617, 1.645, 1.673, 1.701 };',
        'auto interp = [&](double y, vector<double>& curve) {',
        '    if (y < 0) return curve[0];',
        '    unsigned int i = (unsigned int)(y * 100);',
        '    if (i >= curve.size()) return curve.back();',
        '    double f = y * 100 - i;',
        '    return (1 - f) * curve[i] + f * curve[i + 1];',
        '};',
        
        'auto odds = [&](double v, double lo) {',
        '    double a = v / (1 - v);',
        '    return a * exp(lo) / (a * exp(lo) + 1);',
        '};',

        'for (unsigned int g = 0; g < P.processes[0].prob.size(); ++g) {',
        '    // To hospital',
        '    P.processes[0].prob[g][0]  = odds(P.processes[0].prob[g][0], x_hosp_rlo);',
        '    P.processes[0].prob[g][1]  = 1 - P.processes[0].prob[g][0];',
        '    P.processes[10].prob[g][0] = odds(P.processes[10].prob[g][0], x_hosp_rlo${ if (v2) "+ x_v2_hosp_rlo" else "" });',
        '    P.processes[10].prob[g][1] = 1 - P.processes[10].prob[g][0];',
        '    // To ICU',
        '    P.processes[1].prob[g][0]  = odds(P.processes[1].prob[g][0], x_icu_rlo);',
        '    P.processes[1].prob[g][1]  = 1 - P.processes[1].prob[g][0];',
        '    P.processes[11].prob[g][0] = odds(P.processes[11].prob[g][0], x_icu_rlo${ if (v2) "+ x_v2_icu_rlo" else "" });',
        '    P.processes[11].prob[g][1] = 1 - P.processes[11].prob[g][0];',
        '    // To death',
        '    P.processes[4].prob[g][0]  = odds(P.processes[4].prob[g][0], x_cfr_rlo);',
        '    P.processes[4].prob[g][1]  = 1 - P.processes[4].prob[g][0];',
        '    P.processes[14].prob[g][0] = odds(P.processes[14].prob[g][0], x_cfr_rlo${ if (v2) "+ x_v2_cfr_rlo" else "" });',
        '    P.processes[14].prob[g][1] = 1 - P.processes[14].prob[g][0];',
        '    // Relative susceptibility',
        '    P.pop[0].u[g]  = P.pop[0].u[g]  * x_u;',
        if (v2_relu) {
        '    P.pop[0].u2[g] = P.pop[0].u2[g] * x_u * x_v2_relu;'
        } else {
        '    P.pop[0].u2[g] = P.pop[0].u2[g] * x_u;'
        },
        if (v2_immesc) {
        '    P.pop[0].pi2_r[g] = x_v2_immesc;'
        } else '',
        '}',

        if (v2_ch_u && v2_relu) {
        'P.pop[0].u2[0] = P.pop[0].u2[0] + x_v2_ch_u * x_u * x_v2_relu;
        P.pop[0].u2[1] = P.pop[0].u2[1] + x_v2_ch_u * x_u * x_v2_relu;
        P.pop[0].u2[2] = P.pop[0].u2[2] + x_v2_ch_u * x_u * x_v2_relu;
        P.pop[0].u2[3] = P.pop[0].u2[3] + x_v2_ch_u * x_u * x_v2_relu;'
        } else if (v2_ch_u && !v2_relu) {
        'P.pop[0].u2[0] = P.pop[0].u2[0] + x_v2_ch_u * x_u;
        P.pop[0].u2[1] = P.pop[0].u2[1] + x_v2_ch_u * x_u;
        P.pop[0].u2[2] = P.pop[0].u2[2] + x_v2_ch_u * x_u;
        P.pop[0].u2[3] = P.pop[0].u2[3] + x_v2_ch_u * x_u;'
        } else '',

        if (v2_latdur) {
        'P.pop[0].dE2 = delay_gamma(x_v2_latdur * 2.5, 2.5, 15, 0.25);'
        } else '',

        if (v2_infdur) {
        'P.pop[0].dIp2 = delay_gamma(x_v2_infdur * 2.5, 4.0, 15, 0.25);
        P.pop[0].dIs2 = delay_gamma(x_v2_infdur * 2.5, 4.0, 15, 0.25);
        P.pop[0].dIa2 = delay_gamma(x_v2_infdur * 5.0, 4.0, 15, 0.25);'
        } else '',
        
        '// Delays to admission to hospital & ICU',
        'P.processes[0]. delays[0] = delay_gamma(x_hosp_admission, 0.71, 60, 0.25);',
        'P.processes[10].delays[0] = delay_gamma(x_hosp_admission, 0.71, 60, 0.25);',
        'P.processes[1]. delays[0] = delay_gamma(x_icu_admission, 1.91, 60, 0.25);',
        'P.processes[11].delays[0] = delay_gamma(x_icu_admission, 1.91, 60, 0.25);',
        '// Death',
        'P.processes[4] .delays[0] = delay_gamma(x_death_mean, 2.2, 60, 0.25);',
        'P.processes[14].delays[0] = delay_gamma(x_death_mean, 2.2, 60, 0.25);',
        '// Seeding of original and variant strain',
        'P.pop[0].seed_times = seq((int)x_tS, (int)x_tS + 27);',
        if (v2) {
        'P.pop[0].seed_times2 = vector<double>(10, x_v2_when);'
        } else {
        'P.pop[0].seed_times2 = vector<double>(1, 99999);'
        },

        'auto asc = [&](double x, double y0, double y1, double s0, double s1) {',
        '    double xx = s0 + x * (s1 - s0);',
        '    double h0 = exp(s0) / (1 + exp(s0));',
        '    double h1 = exp(s1) / (1 + exp(s1));',
        '    double h = (exp(xx) / (1 + exp(xx)) - h0) / (h1 - h0);',
        '    return y0 + (y1 - y0) * h;',
        '};',

        # Contact adjustment
        'for (unsigned int i = 0; i < P.changes.ch[4].times.size(); ++i) {',
        '    double tx = double(i) / (P.changes.ch[4].times.size() - 1.0);',
        '    P.changes.ch[4].values[i] = vector<double>(8, asc(tx, 1.0, x_contact_final, -x_contact_s0, x_contact_s1));',
        '}',

        # Sep boost
        'P.changes.ch[5].values[0] = vector<double>(8, x_sep_boost);',
        'P.changes.ch[5].times[0] = x_sep_when;',

        # fitting of google mobility indices
        'for (unsigned int k : vector<unsigned int> { 0, 2, 3 }) {',
        '    for (unsigned int i = 0; i < P.changes.ch[k].times.size(); ++i) {',
        '        //double resi = P.changes.ch[k].values[i][0];',
        '        double wplc = P.changes.ch[k].values[i][1];',
        '        double groc = P.changes.ch[k].values[i][2];',
        '        double rtrc = P.changes.ch[k].values[i][3];',
        '        double trns = P.changes.ch[k].values[i][4];',
        '        double othx = rtrc * 0.345 + trns * 0.445 + groc * 0.210;',
        '        double t = P.changes.ch[k].times[i];',

        # from CoMix analysis
        '        double home = asc(min(1.0, t / 365.0), 1.0, 1.545019 / 3.875622, -79 * 0.6, 286 * 0.6);',
        '        double work = interp(wplc, work_curve);',
        '        double scho = (t >= 81 && t <= 244) ? 0 : 1;',
        '        double othe = interp(othx, other_curve);',
        '        P.changes.ch[k].values[i] = { home, work, scho, othe, home, work, scho, othe };',
        '    }',
        '}',
        
        # old... for fIs changes - overwriting changes in fIs here.
        'for (unsigned int i = 0; i < P.changes.ch[1].times.size(); ++i) {',
        '    P.changes.ch[1].values[i] = vector<double>(16, 1.0);',
        '}',
        .sep = "\n    ", .open = "${", .close = "}"
    )
}

# create c++ likelihood components
cpp_likI_voc = function(params, ld, sitreps, sero, virus, sgtfd, popid, max_date, priors, death_cutoff, use_sgtf)
{
    refdate = ld[, max(date)] + 1;
    ld = copy(ld[date <= max_date]);
    sitreps = copy(sitreps[date <= max_date]);
    sero = copy(sero[End.date <= max_date]);
    virus = copy(virus[End.date <= max_date]);
    sgtfd = copy(sgtfd[date <= max_date]);

    glue::glue(
        named_params(priors),
        'std::vector<double> death_v =  ${ cpp_vec(ld[!is.na(N), N]) };',
        'std::vector<double> death_t =  ${ cpp_vec(ld[!is.na(N), as.numeric(ymd(date) - ymd(params$date0))]) };',
        'std::vector<double> death_tr = ${ cpp_vec(ld[!is.na(N), as.numeric(refdate - ymd(date))]) };',
        'std::vector<double> hosp_i_v = ${ cpp_vec(sitreps[!is.na(n_admitted_diagnosed), n_admitted_diagnosed]) };',
        'std::vector<double> hosp_i_t = ${ cpp_vec(sitreps[!is.na(n_admitted_diagnosed), as.numeric(ymd(date) - ymd(params$date0))]) };',
        'std::vector<double> hosp_p_v = ${ cpp_vec(sitreps[!is.na(n_in_all_beds), n_in_all_beds]) };',
        'std::vector<double> hosp_p_t = ${ cpp_vec(sitreps[!is.na(n_in_all_beds), as.numeric(ymd(date) - ymd(params$date0))]) };',
        'std::vector<double> icu_p_v =  ${ cpp_vec(sitreps[!is.na(n_in_itu), n_in_itu]) };',
        'std::vector<double> icu_p_t =  ${ cpp_vec(sitreps[!is.na(n_in_itu), as.numeric(ymd(date) - ymd(params$date0))]) };',
        'std::vector<double> sero_xi =   ${ cpp_vec(sero[, xi]) };',
        'std::vector<double> sero_om =   ${ cpp_vec(sero[, omega]) };',
        'std::vector<double> sero_al =   ${ cpp_vec(sero[, alpha]) };',
        'std::vector<double> sero_t0 =   ${ cpp_vec(sero[, as.numeric(ymd(Start.date) - ymd(params$date0))]) };',
        'std::vector<double> sero_t1 =   ${ cpp_vec(sero[, as.numeric(ymd(End.date) - ymd(params$date0))]) };',
        'std::vector<double> virus_xi =  ${ cpp_vec(virus[, xi]) };',
        'std::vector<double> virus_om =  ${ cpp_vec(virus[, omega]) };',
        'std::vector<double> virus_al =  ${ cpp_vec(virus[, alpha]) };',
        'std::vector<double> virus_t0 =  ${ cpp_vec(virus[, as.numeric(ymd(Start.date) - ymd(params$date0))]) };',
        'std::vector<double> virus_t1 =  ${ cpp_vec(virus[, as.numeric(ymd(End.date) - ymd(params$date0))]) };',
        if (use_sgtf) {
        'std::vector<double> sgtf_s = ${ cpp_vec(sgtfd[!is.na(sgtf), sgtf]) };
        std::vector<double> sgtf_o = ${ cpp_vec(sgtfd[!is.na(other), other]) };
        std::vector<double> sgtf_t = ${ cpp_vec(sgtfd[!is.na(sgtf), as.numeric(ymd(date) - ymd(params$date0))]) };'
        } else '',
        
        'auto size_param = [](double disp) { return 1.0 / (disp * disp); };',

        'for (int i = 0; i < int(death_t.size()) - ${death_cutoff}; ++i)',
        '    ll += nbinom(death_v[i], max(0.1, dyn("death_o", death_t[i], {}, {}) + dyn("death2_o", death_t[i], {}, {})), size_param(x_disp_deaths));',
        'for (unsigned int i = 0; i < hosp_i_t.size(); ++i)',
        '    ll += nbinom(hosp_i_v[i], max(0.1, dyn("hosp_undetected_o", hosp_i_t[i], {}, {}) + dyn("hosp_undetected2_o", hosp_i_t[i], {}, {})), size_param(x_disp_hosp_inc));',
        'for (unsigned int i = 0; i < hosp_p_t.size(); ++i)',
        '    ll += nbinom(hosp_p_v[i], max(0.1, dyn("hosp_p", hosp_p_t[i], {}, {}) + dyn("hosp2_p", hosp_p_t[i], {}, {}) - dyn("hosp_undetected_p", hosp_p_t[i], {}, {}) - dyn("hosp_undetected2_p", hosp_p_t[i], {}, {})), size_param(x_disp_hosp_prev));',
        'for (unsigned int i = 0; i < icu_p_t.size(); ++i)',
        '    ll += nbinom(icu_p_v[i], max(0.1, dyn("icu_p", icu_p_t[i], {}, {}) + dyn("icu2_p", icu_p_t[i], {}, {})), size_param(x_disp_icu_prev));',
        'for (unsigned int i = 0; i < virus_t0.size(); ++i) {',
        '    double virus_prev = 0;',
        '    for (double tt = virus_t0[i]; tt <= virus_t1[i]; ++tt)',
        '        virus_prev += dyn("pcr_positive_p", tt, {}, {}) / (virus_t1[i] - virus_t0[i] + 1);',
        '    ll += 10 * skewnorm(virus_prev / ${ sum(params$pop[[1]]$size) }, virus_xi[i], virus_om[i], virus_al[i]);',
        '}',
        'for (unsigned int i = 0; i < sero_t0.size(); ++i) {',
        '    double sero_prev = 0;',
        '    for (double tt = sero_t0[i]; tt <= sero_t1[i]; ++tt)',
        '        sero_prev += dyn("lfia_positive_p", tt, {}, {}) / (sero_t1[i] - sero_t0[i] + 1);',
        '    ll += skewnorm(sero_prev / ${ sum(params$pop[[1]]$size) }, sero_xi[i], sero_om[i], sero_al[i]);',
        '}',
        if (use_sgtf) {
        'for (unsigned int i = 0; i < sgtf_t.size(); ++i) {
            double s1 = dyn("test_o",  sgtf_t[i], {}, {});
            double s2 = dyn("test2_o", sgtf_t[i], {}, {});
            double p2 = (s1 + s2) > 0 ? s2 / (s1 + s2) : 0;
            double model_sgtf = p2 + (1 - p2) * x_v2_sgtf0;
            ll += 10 * bbinom(sgtf_s[i], sgtf_s[i] + sgtf_o[i], model_sgtf, x_v2_conc);
        }'
        } else '',

        .sep = "\n    ", .open = "${", .close = "}"
    )
}

# create c++ observer function
cpp_obsI_voc = function(v2, P.death, P.critical, priors)
{
    glue::glue(
        named_params(priors),
        'auto asc = [&](double x, double y0, double y1, double s0, double s1) {',
        '    double xx = s0 + x * (s1 - s0);',
        '    double h0 = exp(s0) / (1 + exp(s0));',
        '    double h1 = exp(s1) / (1 + exp(s1));',
        '    double h = (exp(xx) / (1 + exp(xx)) - h0) / (h1 - h0);',
        '    return y0 + (y1 - y0) * h;',
        '};',
        
        'auto odds = [&](double v, double lo) {',
        '    double a = v / (1 - v);',
        '    return a * exp(lo) / (a * exp(lo) + 1);',
        '};',
        
        'auto clamp = [&](double v) {',
        '    return max(0.0, min(1.0, v));',
        '};',

        'dyn.Obs(t, 0, 0, 0) = estimate_Rt(P, dyn, t, 0, 50);',
        'dyn.Obs(t, 0, 3, 0) = estimate_R0(P, dyn, t, 0, 50);',
        
        # increase in young person mobility
        'if (t == 182) {',
        '    double mode = 0.2;',
        '    double conc = x_concentration1;',
        '    double constant = 0.2;',
        '    for (unsigned int a = 0; a < P.pop[0].u.size(); ++a) {',
        '        P.pop[0].u[a]  *= (1 - constant) * dbeta((a + 0.5) / P.pop[0].u.size(), mode * (conc - 2) + 1, (1 - mode) * (conc - 2) + 1) + constant;',
        '        P.pop[0].u2[a] *= (1 - constant) * dbeta((a + 0.5) / P.pop[0].u.size(), mode * (conc - 2) + 1, (1 - mode) * (conc - 2) + 1) + constant;',
        '    }',
        '}',
        'if (t == 213) {',
        '    double mode = 0.2;',
        '    double conc_prev = x_concentration1;',
        '    double conc = x_concentration2;',
        '    double constant = 0.2;',
        '    for (unsigned int a = 0; a < P.pop[0].u.size(); ++a) {',
        '        P.pop[0].u[a]  /= (1 - constant) * dbeta((a + 0.5) / P.pop[0].u.size(), mode * (conc_prev - 2) + 1, (1 - mode) * (conc_prev - 2) + 1) + constant;',
        '        P.pop[0].u[a]  *= (1 - constant) * dbeta((a + 0.5) / P.pop[0].u.size(), mode * (conc - 2) + 1, (1 - mode) * (conc - 2) + 1) + constant;',
        '        P.pop[0].u2[a] /= (1 - constant) * dbeta((a + 0.5) / P.pop[0].u.size(), mode * (conc_prev - 2) + 1, (1 - mode) * (conc_prev - 2) + 1) + constant;',
        '        P.pop[0].u2[a] *= (1 - constant) * dbeta((a + 0.5) / P.pop[0].u.size(), mode * (conc - 2) + 1, (1 - mode) * (conc - 2) + 1) + constant;',
        '    }',
        '}',
        'if (t == 244) {',
        '    double mode = 0.2;',
        '    double conc_prev = x_concentration2;',
        '    double conc = x_concentration3;',
        '    double constant = 0.2;',
        '    for (unsigned int a = 0; a < P.pop[0].u.size(); ++a) {',
        '        P.pop[0].u[a]  /= (1 - constant) * dbeta((a + 0.5) / P.pop[0].u.size(), mode * (conc_prev - 2) + 1, (1 - mode) * (conc_prev - 2) + 1) + constant;',
        '        P.pop[0].u[a]  *= (1 - constant) * dbeta((a + 0.5) / P.pop[0].u.size(), mode * (conc - 2) + 1, (1 - mode) * (conc - 2) + 1) + constant;',
        '        P.pop[0].u2[a] /= (1 - constant) * dbeta((a + 0.5) / P.pop[0].u.size(), mode * (conc_prev - 2) + 1, (1 - mode) * (conc_prev - 2) + 1) + constant;',
        '        P.pop[0].u2[a] *= (1 - constant) * dbeta((a + 0.5) / P.pop[0].u.size(), mode * (conc - 2) + 1, (1 - mode) * (conc - 2) + 1) + constant;',
        '    }',
        '}',

        # changing CFR and ICU admission
        'if ((int)t % 7 == 0) {',
        '    double ifr_table[] = ${ cpp_vec(P.death) };',
        '    double icr_table[] = ${ cpp_vec(P.critical) };',
        '    double adj_f1 = asc(clamp(t / 190.0),          1.0, 0.0, -4.0, 1.0);',
        '    double adj_f3 = asc(clamp((t - 240.0) / 90.0), 0.0, 1.0, -5.0, 5.0);',
        '    double adj_f2 = 1.0 - adj_f1 - adj_f3;',
        '    double adj_c = asc(clamp(t/366.0), 0.0, 1.0, -6, 6);',
        '    for (unsigned int g = 0; g < P.processes[0].prob.size(); ++g) {',
        '        double ifr = ifr_table[g];',
        '        double icr = icr_table[g];',
        '        P.processes[4].prob[g][0]  = odds(ifr, adj_f1 * x_cfr_rlo + adj_f2 * x_cfr_rlo2 + adj_f3 * x_cfr_rlo3);',
        '        P.processes[4].prob[g][1]  = 1 - P.processes[4].prob[g][0];',
        '        P.processes[14].prob[g][0] = odds(ifr, adj_f1 * x_cfr_rlo + adj_f2 * x_cfr_rlo2 + adj_f3 * x_cfr_rlo3${ if (v2) "+ x_v2_cfr_rlo" else "" });',
        '        P.processes[14].prob[g][1] = 1 - P.processes[14].prob[g][0];',
        '        P.processes[1].prob[g][0]  = odds(icr, (1 - adj_c) * x_icu_rlo + adj_c * x_icu_rlo2);',
        '        P.processes[1].prob[g][1]  = 1 - P.processes[1].prob[g][0];',
        '        P.processes[11].prob[g][0] = odds(icr, (1 - adj_c) * x_icu_rlo + adj_c * x_icu_rlo2${ if (v2) "+ x_v2_icu_rlo" else "" });',
        '        P.processes[11].prob[g][1] = 1 - P.processes[11].prob[g][0];',
        '    }',
        '}',
        # changing detection rate in patients admitted to hospital
        'double detection = asc(min(t / 365.0, 1.0), 14, 1, -5.86, 33.4);',
        'P.processes[7].delays[0]  = delay_gamma(detection, 0.59, 60, 0.25);',
        'P.processes[15].delays[0] = delay_gamma(detection, 0.59, 60, 0.25);',
        .sep = "\n    ", .open = "${", .close = "}"
    )
}






# create cpp changes
cpp_chgI_baseline = function(priors)
{
    c(
        named_params(priors),
        'vector<double> work_curve = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.008, 0.021, 0.033, 0.046, 0.058, 0.071, 0.083, 0.096, 0.108, 0.121, 0.133, 0.146, 0.158, 0.171, 0.183, 0.196, 0.208, 0.221, 0.233, 0.246, 0.258, 0.271, 0.283, 0.296, 0.308, 0.321, 0.334, 0.346, 0.359, 0.371, 0.384, 0.397, 0.41, 0.422, 0.435, 0.448, 0.461, 0.474, 0.487, 0.5, 0.513, 0.526, 0.539, 0.552, 0.566, 0.579, 0.592, 0.606, 0.619, 0.633, 0.646, 0.66, 0.674, 0.687, 0.701, 0.715, 0.729, 0.743, 0.757, 0.771, 0.785, 0.799, 0.813, 0.828, 0.842, 0.856, 0.87, 0.885, 0.899, 0.914, 0.928, 0.942, 0.957, 0.971, 0.986, 1, 1.014, 1.029, 1.043, 1.058, 1.072, 1.087, 1.101, 1.115, 1.13, 1.144, 1.159, 1.173, 1.188, 1.202, 1.216, 1.231, 1.245, 1.26, 1.274, 1.289, 1.303, 1.317, 1.332, 1.346, 1.361 };',
        'vector<double> other_curve = { 0.064, 0.066, 0.067, 0.068, 0.069, 0.071, 0.072, 0.073, 0.075, 0.076, 0.077, 0.078, 0.08, 0.081, 0.082, 0.084, 0.085, 0.086, 0.087, 0.089, 0.09, 0.091, 0.092, 0.094, 0.095, 0.096, 0.098, 0.099, 0.1, 0.101, 0.103, 0.104, 0.105, 0.106, 0.108, 0.109, 0.11, 0.112, 0.113, 0.114, 0.116, 0.118, 0.119, 0.121, 0.123, 0.125, 0.128, 0.13, 0.132, 0.135, 0.137, 0.14, 0.143, 0.146, 0.15, 0.154, 0.159, 0.164, 0.169, 0.175, 0.182, 0.19, 0.198, 0.207, 0.217, 0.228, 0.24, 0.252, 0.266, 0.28, 0.295, 0.31, 0.327, 0.344, 0.361, 0.379, 0.398, 0.418, 0.438, 0.459, 0.48, 0.502, 0.525, 0.549, 0.572, 0.597, 0.621, 0.647, 0.672, 0.698, 0.725, 0.751, 0.778, 0.805, 0.833, 0.86, 0.888, 0.916, 0.944, 0.972, 1, 1.028, 1.056, 1.084, 1.112, 1.14, 1.168, 1.196, 1.224, 1.252, 1.28, 1.308, 1.337, 1.365, 1.393, 1.421, 1.449, 1.477, 1.505, 1.533, 1.561, 1.589, 1.617, 1.645, 1.673, 1.701 };',
        'auto interp = [&](double y, vector<double>& curve) {',
        '    if (y < 0) return curve[0];',
        '    unsigned int i = (unsigned int)(y * 100);',
        '    if (i >= curve.size()) return curve.back();',
        '    double f = y * 100 - i;',
        '    return (1 - f) * curve[i] + f * curve[i + 1];',
        '};',
        
        'auto odds = [&](double v, double lo) {',
        '    double a = v / (1 - v);',
        '    return a * exp(lo) / (a * exp(lo) + 1);',
        '};',

        'for (unsigned int g = 0; g < P.processes[0].prob.size(); ++g) {',
        '    // To hospital',
        '    P.processes[0].prob[g][0]  = odds(P.processes[0].prob[g][0], x_hosp_rlo);',
        '    P.processes[0].prob[g][1]  = 1 - P.processes[0].prob[g][0];',
        '    P.processes[10].prob[g][0] = odds(P.processes[10].prob[g][0], x_hosp_rlo /* + v2 */);',
        '    P.processes[10].prob[g][1] = 1 - P.processes[10].prob[g][0];',
        '    // To ICU',
        '    P.processes[1].prob[g][0]  = odds(P.processes[1].prob[g][0], x_icu_rlo);',
        '    P.processes[1].prob[g][1]  = 1 - P.processes[1].prob[g][0];',
        '    P.processes[11].prob[g][0] = odds(P.processes[11].prob[g][0], x_icu_rlo /* + v2 */);',
        '    P.processes[11].prob[g][1] = 1 - P.processes[11].prob[g][0];',
        '    // To death',
        '    P.processes[4].prob[g][0]  = odds(P.processes[4].prob[g][0], x_cfr_rlo);',
        '    P.processes[4].prob[g][1]  = 1 - P.processes[4].prob[g][0];',
        '    P.processes[14].prob[g][0] = odds(P.processes[14].prob[g][0], x_cfr_rlo /* + v2 */);',
        '    P.processes[14].prob[g][1] = 1 - P.processes[14].prob[g][0];',
        '    // Relative susceptibility',
        '    P.pop[0].u[g]  = P.pop[0].u[g]  * x_u;',
        '    P.pop[0].u2[g] = P.pop[0].u2[g] * x_u /* * v2 */;',
        '}',
        
        '// Delays to admission to hospital & ICU',
        'P.processes[0]. delays[0] = delay_gamma(x_hosp_admission, 0.71, 60, 0.25);',
        'P.processes[10].delays[0] = delay_gamma(x_hosp_admission, 0.71, 60, 0.25);',
        'P.processes[1]. delays[0] = delay_gamma(x_icu_admission, 1.91, 60, 0.25);',
        'P.processes[11].delays[0] = delay_gamma(x_icu_admission, 1.91, 60, 0.25);',
        '// Death',
        'P.processes[4] .delays[0] = delay_gamma(x_death_mean, 2.2, 60, 0.25);',
        'P.processes[14].delays[0] = delay_gamma(x_death_mean, 2.2, 60, 0.25);',
        '// Seeding of original and variant strain',
        'P.pop[0].seed_times = seq((int)x_tS, (int)x_tS + 27);',
        'P.pop[0].seed_times2 = vector<double>(1, 99999);',

        'auto asc = [&](double x, double y0, double y1, double s0, double s1) {',
        '    double xx = s0 + x * (s1 - s0);',
        '    double h0 = exp(s0) / (1 + exp(s0));',
        '    double h1 = exp(s1) / (1 + exp(s1));',
        '    double h = (exp(xx) / (1 + exp(xx)) - h0) / (h1 - h0);',
        '    return y0 + (y1 - y0) * h;',
        '};',

        # Contact adjustment
        'for (unsigned int i = 0; i < P.changes.ch[4].times.size(); ++i) {',
        '    double tx = double(i) / (P.changes.ch[4].times.size() - 1.0);',
        '    P.changes.ch[4].values[i] = vector<double>(8, asc(tx, 1.0, x_contact_final, -x_contact_s0, x_contact_s1));',
        '}',

        # Sep boost
        'P.changes.ch[5].values[0] = vector<double>(8, x_sep_boost);',
        'P.changes.ch[5].times[0] = x_sep_when;',

        # fitting of google mobility indices
        'for (unsigned int k : vector<unsigned int> { 0, 2, 3 }) {',
        '    for (unsigned int i = 0; i < P.changes.ch[k].times.size(); ++i) {',
        '        //double resi = P.changes.ch[k].values[i][0];',
        '        double wplc = P.changes.ch[k].values[i][1];',
        '        double groc = P.changes.ch[k].values[i][2];',
        '        double rtrc = P.changes.ch[k].values[i][3];',
        '        double trns = P.changes.ch[k].values[i][4];',
        '        double othx = rtrc * 0.345 + trns * 0.445 + groc * 0.210;',
        '        double t = P.changes.ch[k].times[i];',

        # from CoMix analysis
        '        double home = asc(min(1.0, t / 365.0), 1.0, 1.545019 / 3.875622, -79 * 0.6, 286 * 0.6);',
        '        double work = interp(wplc, work_curve);',
        '        double scho = (t >= 81 && t <= 244) ? 0 : 1;',
        '        double othe = interp(othx, other_curve);',
        '        P.changes.ch[k].values[i] = { home, work, scho, othe, home, work, scho, othe };',
        '    }',
        '}',
        
        # old... for fIs changes - overwriting changes in fIs here.
        'for (unsigned int i = 0; i < P.changes.ch[1].times.size(); ++i) {',
        '    P.changes.ch[1].values[i] = vector<double>(16, 1.0);',
        '}'
    )
}


# create c++ likelihood components
cpp_likI_baseline = function(params, ld, sitreps, sero, virus, sgtfd, popid, max_date, priors, death_confirm = FALSE)
{
    refdate = ld[, max(date)] + 1;
    ld = copy(ld[date <= max_date]);
    sitreps = copy(sitreps[date <= max_date]);
    sero = copy(sero[End.date <= max_date]);
    virus = copy(virus[End.date <= max_date]);
    sgtfd = copy(sgtfd[date <= max_date]);

    glue::glue(
        named_params(priors),
        'std::vector<double> death_v =  ${ cpp_vec(ld[!is.na(N), N]) };',
        'std::vector<double> death_t =  ${ cpp_vec(ld[!is.na(N), as.numeric(ymd(date) - ymd(params$date0))]) };',
        'std::vector<double> death_tr = ${ cpp_vec(ld[!is.na(N), as.numeric(refdate - ymd(date))]) };',
        'std::vector<double> hosp_i_v = ${ cpp_vec(sitreps[!is.na(n_admitted_diagnosed), n_admitted_diagnosed]) };',
        'std::vector<double> hosp_i_t = ${ cpp_vec(sitreps[!is.na(n_admitted_diagnosed), as.numeric(ymd(date) - ymd(params$date0))]) };',
        'std::vector<double> hosp_p_v = ${ cpp_vec(sitreps[!is.na(n_in_all_beds), n_in_all_beds]) };',
        'std::vector<double> hosp_p_t = ${ cpp_vec(sitreps[!is.na(n_in_all_beds), as.numeric(ymd(date) - ymd(params$date0))]) };',
        'std::vector<double> icu_p_v =  ${ cpp_vec(sitreps[!is.na(n_in_itu), n_in_itu]) };',
        'std::vector<double> icu_p_t =  ${ cpp_vec(sitreps[!is.na(n_in_itu), as.numeric(ymd(date) - ymd(params$date0))]) };',
        'std::vector<double> sero_xi =   ${ cpp_vec(sero[, xi]) };',
        'std::vector<double> sero_om =   ${ cpp_vec(sero[, omega]) };',
        'std::vector<double> sero_al =   ${ cpp_vec(sero[, alpha]) };',
        'std::vector<double> sero_t0 =   ${ cpp_vec(sero[, as.numeric(ymd(Start.date) - ymd(params$date0))]) };',
        'std::vector<double> sero_t1 =   ${ cpp_vec(sero[, as.numeric(ymd(End.date) - ymd(params$date0))]) };',
        'std::vector<double> virus_xi =  ${ cpp_vec(virus[, xi]) };',
        'std::vector<double> virus_om =  ${ cpp_vec(virus[, omega]) };',
        'std::vector<double> virus_al =  ${ cpp_vec(virus[, alpha]) };',
        'std::vector<double> virus_t0 =  ${ cpp_vec(virus[, as.numeric(ymd(Start.date) - ymd(params$date0))]) };',
        'std::vector<double> virus_t1 =  ${ cpp_vec(virus[, as.numeric(ymd(End.date) - ymd(params$date0))]) };',
        
        'auto size_param = [](double disp) { return 1.0 / (disp * disp); };',

        'for (unsigned int i = 0; i < death_t.size(); ++i)',
        if (death_confirm) {
        '    ll += nbinom_gammaconf(death_v[i], max(0.1, dyn("death_o", death_t[i], {}, {}) + dyn("death2_o", death_t[i], {}, {})), size_param(x_disp_deaths), death_tr[i], 4.5, 2.7);'
        } else {
        '    ll += nbinom(death_v[i], max(0.1, dyn("death_o", death_t[i], {}, {}) + dyn("death2_o", death_t[i], {}, {})), size_param(x_disp_deaths));'
        },
        'for (unsigned int i = 0; i < hosp_i_t.size(); ++i)',
        '    ll += nbinom(hosp_i_v[i], max(0.1, dyn("hosp_undetected_o", hosp_i_t[i], {}, {}) + dyn("hosp_undetected2_o", hosp_i_t[i], {}, {})), size_param(x_disp_hosp_inc));',
        'for (unsigned int i = 0; i < hosp_p_t.size(); ++i)',
        '    ll += nbinom(hosp_p_v[i], max(0.1, dyn("hosp_p", hosp_p_t[i], {}, {}) + dyn("hosp2_p", hosp_p_t[i], {}, {}) - dyn("hosp_undetected_p", hosp_p_t[i], {}, {}) - dyn("hosp_undetected2_p", hosp_p_t[i], {}, {})), size_param(x_disp_hosp_prev));',
        'for (unsigned int i = 0; i < icu_p_t.size(); ++i)',
        '    ll += nbinom(icu_p_v[i], max(0.1, dyn("icu_p", icu_p_t[i], {}, {}) + dyn("icu2_p", icu_p_t[i], {}, {})), size_param(x_disp_icu_prev));',
        'for (unsigned int i = 0; i < virus_t0.size(); ++i) {',
        '    double virus_prev = 0;',
        '    for (double tt = virus_t0[i]; tt <= virus_t1[i]; ++tt)',
        '        virus_prev += dyn("pcr_positive_p", tt, {}, {}) / (virus_t1[i] - virus_t0[i] + 1);',
        '    ll += 10 * skewnorm(virus_prev / ${ sum(params$pop[[1]]$size) }, virus_xi[i], virus_om[i], virus_al[i]);',
        '}',
        'for (unsigned int i = 0; i < sero_t0.size(); ++i) {',
        '    double sero_prev = 0;',
        '    for (double tt = sero_t0[i]; tt <= sero_t1[i]; ++tt)',
        '        sero_prev += dyn("lfia_positive_p", tt, {}, {}) / (sero_t1[i] - sero_t0[i] + 1);',
        '    ll += skewnorm(sero_prev / ${ sum(params$pop[[1]]$size) }, sero_xi[i], sero_om[i], sero_al[i]);',
        '}',

        .sep = "\n    ", .open = "${", .close = "}"
    )
}


# Observer function
cpp_obsI_baseline = function(P.death, P.critical, priors)
{
    c(
        named_params(priors),
        'auto asc = [&](double x, double y0, double y1, double s0, double s1) {',
        '    double xx = s0 + x * (s1 - s0);',
        '    double h0 = exp(s0) / (1 + exp(s0));',
        '    double h1 = exp(s1) / (1 + exp(s1));',
        '    double h = (exp(xx) / (1 + exp(xx)) - h0) / (h1 - h0);',
        '    return y0 + (y1 - y0) * h;',
        '};',
        
        'auto odds = [&](double v, double lo) {',
        '    double a = v / (1 - v);',
        '    return a * exp(lo) / (a * exp(lo) + 1);',
        '};',
        
        'auto clamp = [&](double v) {',
        '    return max(0.0, min(1.0, v));',
        '};',

        'dyn.Obs(t, 0, 0, 0) = estimate_Rt(P, dyn, t, 0, 50);',
        'dyn.Obs(t, 0, 3, 0) = estimate_R0(P, dyn, t, 0, 50);',
        
        # increase in young person mobility
        'if (t == 182) {',
        '    double mode = 0.2;',
        '    double conc = x_concentration1;',
        '    double constant = 0.2;',
        '    for (unsigned int a = 0; a < P.pop[0].u.size(); ++a) {',
        '        P.pop[0].u[a]  *= (1 - constant) * dbeta((a + 0.5) / P.pop[0].u.size(), mode * (conc - 2) + 1, (1 - mode) * (conc - 2) + 1) + constant;',
        '        P.pop[0].u2[a] *= (1 - constant) * dbeta((a + 0.5) / P.pop[0].u.size(), mode * (conc - 2) + 1, (1 - mode) * (conc - 2) + 1) + constant;',
        '    }',
        '}',
        'if (t == 213) {',
        '    double mode = 0.2;',
        '    double conc_prev = x_concentration1;',
        '    double conc = x_concentration2;',
        '    double constant = 0.2;',
        '    for (unsigned int a = 0; a < P.pop[0].u.size(); ++a) {',
        '        P.pop[0].u[a]  /= (1 - constant) * dbeta((a + 0.5) / P.pop[0].u.size(), mode * (conc_prev - 2) + 1, (1 - mode) * (conc_prev - 2) + 1) + constant;',
        '        P.pop[0].u[a]  *= (1 - constant) * dbeta((a + 0.5) / P.pop[0].u.size(), mode * (conc - 2) + 1, (1 - mode) * (conc - 2) + 1) + constant;',
        '        P.pop[0].u2[a] /= (1 - constant) * dbeta((a + 0.5) / P.pop[0].u.size(), mode * (conc_prev - 2) + 1, (1 - mode) * (conc_prev - 2) + 1) + constant;',
        '        P.pop[0].u2[a] *= (1 - constant) * dbeta((a + 0.5) / P.pop[0].u.size(), mode * (conc - 2) + 1, (1 - mode) * (conc - 2) + 1) + constant;',
        '    }',
        '}',
        'if (t == 244) {',
        '    double mode = 0.2;',
        '    double conc_prev = x_concentration2;',
        '    double conc = x_concentration3;',
        '    double constant = 0.2;',
        '    for (unsigned int a = 0; a < P.pop[0].u.size(); ++a) {',
        '        P.pop[0].u[a]  /= (1 - constant) * dbeta((a + 0.5) / P.pop[0].u.size(), mode * (conc_prev - 2) + 1, (1 - mode) * (conc_prev - 2) + 1) + constant;',
        '        P.pop[0].u[a]  *= (1 - constant) * dbeta((a + 0.5) / P.pop[0].u.size(), mode * (conc - 2) + 1, (1 - mode) * (conc - 2) + 1) + constant;',
        '        P.pop[0].u2[a] /= (1 - constant) * dbeta((a + 0.5) / P.pop[0].u.size(), mode * (conc_prev - 2) + 1, (1 - mode) * (conc_prev - 2) + 1) + constant;',
        '        P.pop[0].u2[a] *= (1 - constant) * dbeta((a + 0.5) / P.pop[0].u.size(), mode * (conc - 2) + 1, (1 - mode) * (conc - 2) + 1) + constant;',
        '    }',
        '}',

        # changing CFR and ICU admission
        'if ((int)t % 7 == 0) {',
        '    double ifr_table[] = {', paste(P.death, collapse = ", "), '};',
        '    double icr_table[] = {', paste(P.critical, collapse = ", "), '};',
        '    double adj_f1 = asc(clamp(t / 190.0),          1.0, 0.0, -4.0, 1.0);',
        '    double adj_f3 = asc(clamp((t - 240.0) / 90.0), 0.0, 1.0, -5.0, 5.0);',
        '    double adj_f2 = 1.0 - adj_f1 - adj_f3;',
        '    double adj_c = asc(clamp(t/366.0), 0.0, 1.0, -6, 6);',
        '    for (unsigned int g = 0; g < P.processes[0].prob.size(); ++g) {',
        '        double ifr = ifr_table[g];',
        '        double icr = icr_table[g];',
        '        P.processes[4].prob[g][0]  = odds(ifr, adj_f1 * x_cfr_rlo + adj_f2 * x_cfr_rlo2 + adj_f3 * x_cfr_rlo3);',
        '        P.processes[4].prob[g][1]  = 1 - P.processes[4].prob[g][0];',
        '        P.processes[14].prob[g][0] = odds(ifr, adj_f1 * x_cfr_rlo + adj_f2 * x_cfr_rlo2 + adj_f3 * x_cfr_rlo3 /* + v2 */);',
        '        P.processes[14].prob[g][1] = 1 - P.processes[14].prob[g][0];',
        '        P.processes[1].prob[g][0]  = odds(icr, (1 - adj_c) * x_icu_rlo + adj_c * x_icu_rlo2);',
        '        P.processes[1].prob[g][1]  = 1 - P.processes[1].prob[g][0];',
        '        P.processes[11].prob[g][0] = odds(icr, (1 - adj_c) * x_icu_rlo + adj_c * x_icu_rlo2 /* + v2 */);',
        '        P.processes[11].prob[g][1] = 1 - P.processes[11].prob[g][0];',
        '    }',
        '}',
        # changing detection rate in patients admitted to hospital
        'double detection = asc(min(t / 365.0, 1.0), 14, 1, -5.86, 33.4);',
        'P.processes[7].delays[0]  = delay_gamma(detection, 0.59, 60, 0.25);',
        'P.processes[15].delays[0] = delay_gamma(detection, 0.59, 60, 0.25);'
    )
}





















# create cpp changes
cpp_chgI = function(variant2 = TRUE)
{
    c(
        'vector<double> work_curve = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.008, 0.021, 0.033, 0.046, 0.058, 0.071, 0.083, 0.096, 0.108, 0.121, 0.133, 0.146, 0.158, 0.171, 0.183, 0.196, 0.208, 0.221, 0.233, 0.246, 0.258, 0.271, 0.283, 0.296, 0.308, 0.321, 0.334, 0.346, 0.359, 0.371, 0.384, 0.397, 0.41, 0.422, 0.435, 0.448, 0.461, 0.474, 0.487, 0.5, 0.513, 0.526, 0.539, 0.552, 0.566, 0.579, 0.592, 0.606, 0.619, 0.633, 0.646, 0.66, 0.674, 0.687, 0.701, 0.715, 0.729, 0.743, 0.757, 0.771, 0.785, 0.799, 0.813, 0.828, 0.842, 0.856, 0.87, 0.885, 0.899, 0.914, 0.928, 0.942, 0.957, 0.971, 0.986, 1, 1.014, 1.029, 1.043, 1.058, 1.072, 1.087, 1.101, 1.115, 1.13, 1.144, 1.159, 1.173, 1.188, 1.202, 1.216, 1.231, 1.245, 1.26, 1.274, 1.289, 1.303, 1.317, 1.332, 1.346, 1.361 };',
        'vector<double> other_curve = { 0.064, 0.066, 0.067, 0.068, 0.069, 0.071, 0.072, 0.073, 0.075, 0.076, 0.077, 0.078, 0.08, 0.081, 0.082, 0.084, 0.085, 0.086, 0.087, 0.089, 0.09, 0.091, 0.092, 0.094, 0.095, 0.096, 0.098, 0.099, 0.1, 0.101, 0.103, 0.104, 0.105, 0.106, 0.108, 0.109, 0.11, 0.112, 0.113, 0.114, 0.116, 0.118, 0.119, 0.121, 0.123, 0.125, 0.128, 0.13, 0.132, 0.135, 0.137, 0.14, 0.143, 0.146, 0.15, 0.154, 0.159, 0.164, 0.169, 0.175, 0.182, 0.19, 0.198, 0.207, 0.217, 0.228, 0.24, 0.252, 0.266, 0.28, 0.295, 0.31, 0.327, 0.344, 0.361, 0.379, 0.398, 0.418, 0.438, 0.459, 0.48, 0.502, 0.525, 0.549, 0.572, 0.597, 0.621, 0.647, 0.672, 0.698, 0.725, 0.751, 0.778, 0.805, 0.833, 0.86, 0.888, 0.916, 0.944, 0.972, 1, 1.028, 1.056, 1.084, 1.112, 1.14, 1.168, 1.196, 1.224, 1.252, 1.28, 1.308, 1.337, 1.365, 1.393, 1.421, 1.449, 1.477, 1.505, 1.533, 1.561, 1.589, 1.617, 1.645, 1.673, 1.701 };',
        'auto interp = [&](double x, vector<double>& curve) {',
        '    if (x < 0) return curve[0];',
        '    unsigned int i = (unsigned int)(x * 100);',
        '    if (i >= curve.size()) return curve.back();',
        '    double f = x * 100 - i;',
        '    return (1 - f) * curve[i] + f * curve[i + 1];',
        '};',
        
        'auto odds = [&](double x, double lo) {',
        '    double a = x / (1 - x);',
        '    return a * exp(lo) / (a * exp(lo) + 1);',
        '};',

        'for (unsigned int g = 0; g < P.processes[0].prob.size(); ++g) {',
        '    // To hospital',
        '    P.processes[0].prob[g][0]  = odds(P.processes[0].prob[g][0], x[7]);',
        '    P.processes[0].prob[g][1]  = 1 - P.processes[0].prob[g][0];',
        '    P.processes[10].prob[g][0] = odds(P.processes[10].prob[g][0], x[7] + x[20]);',
        '    P.processes[10].prob[g][1] = 1 - P.processes[10].prob[g][0];',
        '    // To ICU',
        '    P.processes[1].prob[g][0]  = odds(P.processes[1].prob[g][0], x[7] + x[6]);',
        '    P.processes[1].prob[g][1]  = 1 - P.processes[1].prob[g][0];',
        '    P.processes[11].prob[g][0] = odds(P.processes[11].prob[g][0], x[7] + x[6] + x[20]);',
        '    P.processes[11].prob[g][1] = 1 - P.processes[11].prob[g][0];',
        '    // To death',
        '    P.processes[4].prob[g][0]  = min(1.0, P.processes[4].prob[g][0] * x[5]);',
        '    P.processes[4].prob[g][1]  = 1 - P.processes[4].prob[g][0];',
        '    P.processes[14].prob[g][0] = min(1.0, P.processes[14].prob[g][0] * x[5] * x[21]);',
        '    P.processes[14].prob[g][1] = 1 - P.processes[14].prob[g][0];',
        '    // Relative susceptibility',
        '    P.pop[0].u[g]  = P.pop[0].u[g]  * x[1];',
        '    P.pop[0].u2[g] = P.pop[0].u2[g] * x[1] * x[19];',
        '}',
        
        '// Delays to admission to hospital & ICU',
        'P.processes[0]. delays[0] = delay_gamma(x[4], 0.71, 60, 0.25);',
        'P.processes[10].delays[0] = delay_gamma(x[4], 0.71, 60, 0.25);',
        'P.processes[1]. delays[0] = delay_gamma(x[8], 1.91, 60, 0.25);',
        'P.processes[11].delays[0] = delay_gamma(x[8], 1.91, 60, 0.25);',
        '// Death',
        'P.processes[4] .delays[0] = delay_gamma(x[2], x[3], 60, 0.25);',
        'P.processes[14].delays[0] = delay_gamma(x[2], x[3], 60, 0.25);',
        '// Seeding of original and variant strain',
        'P.pop[0].seed_times = seq((int)x[0], (int)x[0] + 27);',
        if (variant2) {
            'P.pop[0].seed_times2 = vector<double>(10, x[18]);'
        } else {
            'P.pop[0].seed_times2 = vector<double>(10, 99999);'
        },

        'auto asc = [&](double x, double y0, double y1, double s0, double s1) {',
        '    double xx = s0 + x * (s1 - s0);',
        '    double h0 = exp(s0) / (1 + exp(s0));',
        '    double h1 = exp(s1) / (1 + exp(s1));',
        '    double h = (exp(xx) / (1 + exp(xx)) - h0) / (h1 - h0);',
        '    return y0 + (y1 - y0) * h;',
        '};',

        # Contact adjustment
        'for (unsigned int i = 0; i < P.changes.ch[4].times.size(); ++i) {',
        '    double tx = double(i) / (P.changes.ch[4].times.size() - 1.0);',
        '    P.changes.ch[4].values[i] = vector<double>(8, asc(tx, 1.0, x[9], -x[10], x[11]));',
        '}',

        # Sep boost
        'P.changes.ch[5].values[0] = vector<double>(8, x[15]);',
        'P.changes.ch[5].times[0] = x[17];',

        # fitting of google mobility indices
        'for (unsigned int k : vector<unsigned int> { 0, 2, 3 }) {',
        '    for (unsigned int i = 0; i < P.changes.ch[k].times.size(); ++i) {',
        '        //double resi = P.changes.ch[k].values[i][0];',
        '        double wplc = P.changes.ch[k].values[i][1];',
        '        double groc = P.changes.ch[k].values[i][2];',
        '        double rtrc = P.changes.ch[k].values[i][3];',
        '        double trns = P.changes.ch[k].values[i][4];',
        '        double othx = rtrc * 0.345 + trns * 0.445 + groc * 0.210;',
        '        double t = P.changes.ch[k].times[i];',

        # from CoMix analysis
        '        double home = asc(min(1.0, t / 365.0), 1.0, 1.545019 / 3.875622, -79 * 0.6, 286 * 0.6);',
        '        double work = interp(wplc, work_curve);',
        '        double scho = (t >= 81 && t <= 244) ? 0 : 1;',
        '        double othe = interp(othx, other_curve);',
        '        P.changes.ch[k].values[i] = { home, work, scho, othe, home, work, scho, othe };',
        '    }',
        '}',
        
        # old... for fIs changes - overwriting changes in fIs here.
        'for (unsigned int i = 0; i < P.changes.ch[1].times.size(); ++i) {',
        '    P.changes.ch[1].values[i] = vector<double>(16, 1.0);',
        '}'
    )
}


# create cpp changes: multiple region fit
cpp_chgI_multi = function(variant2 = TRUE)
{
    c(
        'vector<double> work_curve = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.008, 0.021, 0.033, 0.046, 0.058, 0.071, 0.083, 0.096, 0.108, 0.121, 0.133, 0.146, 0.158, 0.171, 0.183, 0.196, 0.208, 0.221, 0.233, 0.246, 0.258, 0.271, 0.283, 0.296, 0.308, 0.321, 0.334, 0.346, 0.359, 0.371, 0.384, 0.397, 0.41, 0.422, 0.435, 0.448, 0.461, 0.474, 0.487, 0.5, 0.513, 0.526, 0.539, 0.552, 0.566, 0.579, 0.592, 0.606, 0.619, 0.633, 0.646, 0.66, 0.674, 0.687, 0.701, 0.715, 0.729, 0.743, 0.757, 0.771, 0.785, 0.799, 0.813, 0.828, 0.842, 0.856, 0.87, 0.885, 0.899, 0.914, 0.928, 0.942, 0.957, 0.971, 0.986, 1, 1.014, 1.029, 1.043, 1.058, 1.072, 1.087, 1.101, 1.115, 1.13, 1.144, 1.159, 1.173, 1.188, 1.202, 1.216, 1.231, 1.245, 1.26, 1.274, 1.289, 1.303, 1.317, 1.332, 1.346, 1.361 };',
        'vector<double> other_curve = { 0.064, 0.066, 0.067, 0.068, 0.069, 0.071, 0.072, 0.073, 0.075, 0.076, 0.077, 0.078, 0.08, 0.081, 0.082, 0.084, 0.085, 0.086, 0.087, 0.089, 0.09, 0.091, 0.092, 0.094, 0.095, 0.096, 0.098, 0.099, 0.1, 0.101, 0.103, 0.104, 0.105, 0.106, 0.108, 0.109, 0.11, 0.112, 0.113, 0.114, 0.116, 0.118, 0.119, 0.121, 0.123, 0.125, 0.128, 0.13, 0.132, 0.135, 0.137, 0.14, 0.143, 0.146, 0.15, 0.154, 0.159, 0.164, 0.169, 0.175, 0.182, 0.19, 0.198, 0.207, 0.217, 0.228, 0.24, 0.252, 0.266, 0.28, 0.295, 0.31, 0.327, 0.344, 0.361, 0.379, 0.398, 0.418, 0.438, 0.459, 0.48, 0.502, 0.525, 0.549, 0.572, 0.597, 0.621, 0.647, 0.672, 0.698, 0.725, 0.751, 0.778, 0.805, 0.833, 0.86, 0.888, 0.916, 0.944, 0.972, 1, 1.028, 1.056, 1.084, 1.112, 1.14, 1.168, 1.196, 1.224, 1.252, 1.28, 1.308, 1.337, 1.365, 1.393, 1.421, 1.449, 1.477, 1.505, 1.533, 1.561, 1.589, 1.617, 1.645, 1.673, 1.701 };',
        'auto interp = [&](double x, vector<double>& curve) {',
        '    if (x < 0) return curve[0];',
        '    unsigned int i = (unsigned int)(x * 100);',
        '    if (i >= curve.size()) return curve.back();',
        '    double f = x * 100 - i;',
        '    return (1 - f) * curve[i] + f * curve[i + 1];',
        '};',
        
        'auto odds = [&](double x, double lo) {',
        '    double a = x / (1 - x);',
        '    return a * exp(lo) / (a * exp(lo) + 1);',
        '};',

        'for (unsigned int g = 0; g < P.processes[0].prob.size(); ++g) {',
        '    // To hospital',
        '    P.processes[0].prob[g][0]  = odds(P.processes[0].prob[g][0], x[7]);',
        '    P.processes[0].prob[g][1]  = 1 - P.processes[0].prob[g][0];',
        '    P.processes[10].prob[g][0] = odds(P.processes[10].prob[g][0], x[7] + x[20]);',
        '    P.processes[10].prob[g][1] = 1 - P.processes[10].prob[g][0];',
        '    // To ICU',
        '    P.processes[1].prob[g][0]  = odds(P.processes[1].prob[g][0], x[7] + x[6]);',
        '    P.processes[1].prob[g][1]  = 1 - P.processes[1].prob[g][0];',
        '    P.processes[11].prob[g][0] = odds(P.processes[11].prob[g][0], x[7] + x[6] + x[20]);',
        '    P.processes[11].prob[g][1] = 1 - P.processes[11].prob[g][0];',
        '    // To death',
        '    P.processes[4].prob[g][0]  = min(1.0, P.processes[4].prob[g][0] * x[5]);',
        '    P.processes[4].prob[g][1]  = 1 - P.processes[4].prob[g][0];',
        '    P.processes[14].prob[g][0] = min(1.0, P.processes[14].prob[g][0] * x[5] * x[21]);',
        '    P.processes[14].prob[g][1] = 1 - P.processes[14].prob[g][0];',
        '    // Relative susceptibility',
        '    P.pop[0].u[g]  = P.pop[0].u[g]  * x[1];',
        '    P.pop[0].u2[g] = P.pop[0].u2[g] * x[1] * x[19];',
        '}',
        
        # different duration of infectiousness
        'P.pop[0].dIp2 = delay_gamma(x[24] * 2.5, 4.0, 15, 0.25);',
        'P.pop[0].dIs2 = delay_gamma(x[24] * 2.5, 4.0, 15, 0.25);',
        'P.pop[0].dIa2 = delay_gamma(x[24] * 5.0, 4.0, 15, 0.25);',
        # delay to getting tested
        'P.processes[16].delays[0] = delay_gamma(x[25], 2, 60, 0.25);',
        'P.processes[17].delays[0] = delay_gamma(x[25], 2, 60, 0.25);',

        '// Delays to admission to hospital & ICU',
        'P.processes[0]. delays[0] = delay_gamma(x[4], 0.71, 60, 0.25);',
        'P.processes[10].delays[0] = delay_gamma(x[4], 0.71, 60, 0.25);',
        'P.processes[1]. delays[0] = delay_gamma(x[8], 1.91, 60, 0.25);',
        'P.processes[11].delays[0] = delay_gamma(x[8], 1.91, 60, 0.25);',
        '// Death',
        'P.processes[4] .delays[0] = delay_gamma(x[2], x[3], 60, 0.25);',
        'P.processes[14].delays[0] = delay_gamma(x[2], x[3], 60, 0.25);',
        '// Seeding of original and variant strain',
        'P.pop[0].seed_times = seq((int)x[0], (int)x[0] + 27);',
        if (variant2) {
            'P.pop[0].seed_times2 = vector<double>(10, x[18]);'
        } else {
            'P.pop[0].seed_times2 = vector<double>(10, 99999);'
        },

        'auto asc = [&](double x, double y0, double y1, double s0, double s1) {',
        '    double xx = s0 + x * (s1 - s0);',
        '    double h0 = exp(s0) / (1 + exp(s0));',
        '    double h1 = exp(s1) / (1 + exp(s1));',
        '    double h = (exp(xx) / (1 + exp(xx)) - h0) / (h1 - h0);',
        '    return y0 + (y1 - y0) * h;',
        '};',

        # Contact adjustment
        'for (unsigned int i = 0; i < P.changes.ch[4].times.size(); ++i) {',
        '    double tx = double(i) / (P.changes.ch[4].times.size() - 1.0);',
        '    P.changes.ch[4].values[i] = vector<double>(8, asc(tx, 1.0, x[9], -x[10], x[11]));',
        '}',

        # Sep boost
        'P.changes.ch[5].values[0] = vector<double>(8, x[15]);',
        'P.changes.ch[5].times[0] = x[17];',

        # fitting of google mobility indices
        'for (unsigned int k : vector<unsigned int> { 0, 2, 3 }) {',
        '    for (unsigned int i = 0; i < P.changes.ch[k].times.size(); ++i) {',
        '        //double resi = P.changes.ch[k].values[i][0];',
        '        double wplc = P.changes.ch[k].values[i][1];',
        '        double groc = P.changes.ch[k].values[i][2];',
        '        double rtrc = P.changes.ch[k].values[i][3];',
        '        double trns = P.changes.ch[k].values[i][4];',
        '        double othx = rtrc * 0.345 + trns * 0.445 + groc * 0.210;',
        '        double t = P.changes.ch[k].times[i];',

        # from CoMix analysis
        '        double home = asc(min(1.0, t / 365.0), 1.0, 1.545019 / 3.875622, -79 * 0.6, 286 * 0.6);',
        '        double work = interp(wplc, work_curve);',
        '        double scho = (t >= 81 && t <= 244) ? 0 : 1;',
        '        double othe = interp(othx, other_curve);',
        '        P.changes.ch[k].values[i] = { home, work, scho, othe, home, work, scho, othe };',
        '    }',
        '}',
        
        # old... for fIs changes - overwriting changes in fIs here.
        'for (unsigned int i = 0; i < P.changes.ch[1].times.size(); ++i) {',
        '    P.changes.ch[1].values[i] = vector<double>(16, 1.0);',
        '}'
    )
}


# create cpp changes: multiple region fit with everything
cpp_chgI_everything = function(variant2 = TRUE)
{
    c(
        'vector<double> work_curve = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.008, 0.021, 0.033, 0.046, 0.058, 0.071, 0.083, 0.096, 0.108, 0.121, 0.133, 0.146, 0.158, 0.171, 0.183, 0.196, 0.208, 0.221, 0.233, 0.246, 0.258, 0.271, 0.283, 0.296, 0.308, 0.321, 0.334, 0.346, 0.359, 0.371, 0.384, 0.397, 0.41, 0.422, 0.435, 0.448, 0.461, 0.474, 0.487, 0.5, 0.513, 0.526, 0.539, 0.552, 0.566, 0.579, 0.592, 0.606, 0.619, 0.633, 0.646, 0.66, 0.674, 0.687, 0.701, 0.715, 0.729, 0.743, 0.757, 0.771, 0.785, 0.799, 0.813, 0.828, 0.842, 0.856, 0.87, 0.885, 0.899, 0.914, 0.928, 0.942, 0.957, 0.971, 0.986, 1, 1.014, 1.029, 1.043, 1.058, 1.072, 1.087, 1.101, 1.115, 1.13, 1.144, 1.159, 1.173, 1.188, 1.202, 1.216, 1.231, 1.245, 1.26, 1.274, 1.289, 1.303, 1.317, 1.332, 1.346, 1.361 };',
        'vector<double> other_curve = { 0.064, 0.066, 0.067, 0.068, 0.069, 0.071, 0.072, 0.073, 0.075, 0.076, 0.077, 0.078, 0.08, 0.081, 0.082, 0.084, 0.085, 0.086, 0.087, 0.089, 0.09, 0.091, 0.092, 0.094, 0.095, 0.096, 0.098, 0.099, 0.1, 0.101, 0.103, 0.104, 0.105, 0.106, 0.108, 0.109, 0.11, 0.112, 0.113, 0.114, 0.116, 0.118, 0.119, 0.121, 0.123, 0.125, 0.128, 0.13, 0.132, 0.135, 0.137, 0.14, 0.143, 0.146, 0.15, 0.154, 0.159, 0.164, 0.169, 0.175, 0.182, 0.19, 0.198, 0.207, 0.217, 0.228, 0.24, 0.252, 0.266, 0.28, 0.295, 0.31, 0.327, 0.344, 0.361, 0.379, 0.398, 0.418, 0.438, 0.459, 0.48, 0.502, 0.525, 0.549, 0.572, 0.597, 0.621, 0.647, 0.672, 0.698, 0.725, 0.751, 0.778, 0.805, 0.833, 0.86, 0.888, 0.916, 0.944, 0.972, 1, 1.028, 1.056, 1.084, 1.112, 1.14, 1.168, 1.196, 1.224, 1.252, 1.28, 1.308, 1.337, 1.365, 1.393, 1.421, 1.449, 1.477, 1.505, 1.533, 1.561, 1.589, 1.617, 1.645, 1.673, 1.701 };',
        'auto interp = [&](double x, vector<double>& curve) {',
        '    if (x < 0) return curve[0];',
        '    unsigned int i = (unsigned int)(x * 100);',
        '    if (i >= curve.size()) return curve.back();',
        '    double f = x * 100 - i;',
        '    return (1 - f) * curve[i] + f * curve[i + 1];',
        '};',
        
        'auto odds = [&](double x, double lo) {',
        '    double a = x / (1 - x);',
        '    return a * exp(lo) / (a * exp(lo) + 1);',
        '};',

        'for (unsigned int g = 0; g < P.processes[0].prob.size(); ++g) {',
        '    // To hospital',
        '    P.processes[0].prob[g][0]  = odds(P.processes[0].prob[g][0], x[7]);',
        '    P.processes[0].prob[g][1]  = 1 - P.processes[0].prob[g][0];',
        '    P.processes[10].prob[g][0] = odds(P.processes[10].prob[g][0], x[7] + x[20]);',
        '    P.processes[10].prob[g][1] = 1 - P.processes[10].prob[g][0];',
        '    // To ICU',
        '    P.processes[1].prob[g][0]  = odds(P.processes[1].prob[g][0], x[7] + x[6]);',
        '    P.processes[1].prob[g][1]  = 1 - P.processes[1].prob[g][0];',
        '    P.processes[11].prob[g][0] = odds(P.processes[11].prob[g][0], x[7] + x[6] + x[20]);',
        '    P.processes[11].prob[g][1] = 1 - P.processes[11].prob[g][0];',
        '    // To death',
        '    P.processes[4].prob[g][0]  = min(1.0, P.processes[4].prob[g][0] * x[5]);',
        '    P.processes[4].prob[g][1]  = 1 - P.processes[4].prob[g][0];',
        '    P.processes[14].prob[g][0] = min(1.0, P.processes[14].prob[g][0] * x[5] * x[21]);',
        '    P.processes[14].prob[g][1] = 1 - P.processes[14].prob[g][0];',
        '    // Relative susceptibility',
        '    P.pop[0].u[g]  = P.pop[0].u[g]  * x[1];',
        '    P.pop[0].u2[g] = P.pop[0].u2[g] * x[1] * x[19];',
        '    // Immune escape',
        '    P.pop[0].pi2_r[g] = x[26];',
        '}',
        
        '// Increased susceptibility in children',
        'P.pop[0].u2[0] = P.pop[0].u2[0] + x[27] * x[1];',
        'P.pop[0].u2[1] = P.pop[0].u2[1] + x[27] * x[1];',
        'P.pop[0].u2[2] = P.pop[0].u2[2] + x[27] * x[1];',
        'P.pop[0].u2[3] = P.pop[0].u2[3] + x[27] * x[1];',
        
        # different duration of infectiousness
        'P.pop[0].dIp2 = delay_gamma(x[24] * 2.5, 4.0, 15, 0.25);',
        'P.pop[0].dIs2 = delay_gamma(x[24] * 2.5, 4.0, 15, 0.25);',
        'P.pop[0].dIa2 = delay_gamma(x[24] * 5.0, 4.0, 15, 0.25);',
        # delay to getting tested
        'P.processes[16].delays[0] = delay_gamma(x[25], 2, 60, 0.25);',
        'P.processes[17].delays[0] = delay_gamma(x[25], 2, 60, 0.25);',

        '// Delays to admission to hospital & ICU',
        'P.processes[0]. delays[0] = delay_gamma(x[4], 0.71, 60, 0.25);',
        'P.processes[10].delays[0] = delay_gamma(x[4], 0.71, 60, 0.25);',
        'P.processes[1]. delays[0] = delay_gamma(x[8], 1.91, 60, 0.25);',
        'P.processes[11].delays[0] = delay_gamma(x[8], 1.91, 60, 0.25);',
        '// Death',
        'P.processes[4] .delays[0] = delay_gamma(x[2], x[3], 60, 0.25);',
        'P.processes[14].delays[0] = delay_gamma(x[2], x[3], 60, 0.25);',
        '// Seeding of original and variant strain',
        'P.pop[0].seed_times = seq((int)x[0], (int)x[0] + 27);',
        if (variant2) {
            'P.pop[0].seed_times2 = vector<double>(10, x[18]);'
        } else {
            'P.pop[0].seed_times2 = vector<double>(10, 99999);'
        },

        'auto asc = [&](double x, double y0, double y1, double s0, double s1) {',
        '    double xx = s0 + x * (s1 - s0);',
        '    double h0 = exp(s0) / (1 + exp(s0));',
        '    double h1 = exp(s1) / (1 + exp(s1));',
        '    double h = (exp(xx) / (1 + exp(xx)) - h0) / (h1 - h0);',
        '    return y0 + (y1 - y0) * h;',
        '};',

        # Contact adjustment
        'for (unsigned int i = 0; i < P.changes.ch[4].times.size(); ++i) {',
        '    double tx = double(i) / (P.changes.ch[4].times.size() - 1.0);',
        '    P.changes.ch[4].values[i] = vector<double>(8, asc(tx, 1.0, x[9], -x[10], x[11]));',
        '}',

        # Sep boost
        'P.changes.ch[5].values[0] = vector<double>(8, x[15]);',
        'P.changes.ch[5].times[0] = x[17];',

        # fitting of google mobility indices
        'for (unsigned int k : vector<unsigned int> { 0, 2, 3 }) {',
        '    for (unsigned int i = 0; i < P.changes.ch[k].times.size(); ++i) {',
        '        //double resi = P.changes.ch[k].values[i][0];',
        '        double wplc = P.changes.ch[k].values[i][1];',
        '        double groc = P.changes.ch[k].values[i][2];',
        '        double rtrc = P.changes.ch[k].values[i][3];',
        '        double trns = P.changes.ch[k].values[i][4];',
        '        double othx = rtrc * 0.345 + trns * 0.445 + groc * 0.210;',
        '        double t = P.changes.ch[k].times[i];',

        # from CoMix analysis
        '        double home = asc(min(1.0, t / 365.0), 1.0, 1.545019 / 3.875622, -79 * 0.6, 286 * 0.6);',
        '        double work = interp(wplc, work_curve);',
        '        double scho = (t >= 81 && t <= 244) ? 0 : 1;',
        '        double othe = interp(othx, other_curve);',
        '        P.changes.ch[k].values[i] = { home, work, scho, othe, home, work, scho, othe };',
        '    }',
        '}',
        
        # old... for fIs changes - overwriting changes in fIs here.
        'for (unsigned int i = 0; i < P.changes.ch[1].times.size(); ++i) {',
        '    P.changes.ch[1].values[i] = vector<double>(16, 1.0);',
        '}'
    )
}

# create c++ likelihood components
cpp_likI = function(params, ld, sitreps, sero, virus, sgtfd, popid)
{
    refdate = ld[, max(date)] + 1;

    glue::glue(
        'std::vector<double> death_v =  ${ cpp_vec(ld[!is.na(N), N]) };',
        'std::vector<double> death_t =  ${ cpp_vec(ld[!is.na(N), as.numeric(ymd(date) - ymd(params$date0))]) };',
        'std::vector<double> death_tr = ${ cpp_vec(ld[!is.na(N), as.numeric(refdate - ymd(date))]) };',
        'std::vector<double> hosp_i_v = ${ cpp_vec(sitreps[!is.na(n_admitted_diagnosed), n_admitted_diagnosed]) };',
        'std::vector<double> hosp_i_t = ${ cpp_vec(sitreps[!is.na(n_admitted_diagnosed), as.numeric(ymd(date) - ymd(params$date0))]) };',
        'std::vector<double> hosp_p_v = ${ cpp_vec(sitreps[!is.na(n_in_all_beds), n_in_all_beds]) };',
        'std::vector<double> hosp_p_t = ${ cpp_vec(sitreps[!is.na(n_in_all_beds), as.numeric(ymd(date) - ymd(params$date0))]) };',
        'std::vector<double> icu_p_v =  ${ cpp_vec(sitreps[!is.na(n_in_itu), n_in_itu]) };',
        'std::vector<double> icu_p_t =  ${ cpp_vec(sitreps[!is.na(n_in_itu), as.numeric(ymd(date) - ymd(params$date0))]) };',
        'std::vector<double> sero_xi =   ${ cpp_vec(sero[, xi]) };',
        'std::vector<double> sero_om =   ${ cpp_vec(sero[, omega]) };',
        'std::vector<double> sero_al =   ${ cpp_vec(sero[, alpha]) };',
        'std::vector<double> sero_t0 =   ${ cpp_vec(sero[, as.numeric(ymd(Start.date) - ymd(params$date0))]) };',
        'std::vector<double> sero_t1 =   ${ cpp_vec(sero[, as.numeric(ymd(End.date) - ymd(params$date0))]) };',
        'std::vector<double> virus_xi =  ${ cpp_vec(virus[, xi]) };',
        'std::vector<double> virus_om =  ${ cpp_vec(virus[, omega]) };',
        'std::vector<double> virus_al =  ${ cpp_vec(virus[, alpha]) };',
        'std::vector<double> virus_t0 =  ${ cpp_vec(virus[, as.numeric(ymd(Start.date) - ymd(params$date0))]) };',
        'std::vector<double> virus_t1 =  ${ cpp_vec(virus[, as.numeric(ymd(End.date) - ymd(params$date0))]) };',
        'std::vector<double> sgtf_s = ${ cpp_vec(sgtfd[!is.na(sgtf), sgtf]) };',
        'std::vector<double> sgtf_o = ${ cpp_vec(sgtfd[!is.na(other), other]) };',
        'std::vector<double> sgtf_t = ${ cpp_vec(sgtfd[!is.na(sgtf), as.numeric(ymd(date) - ymd(params$date0))]) };',
        
        'for (unsigned int i = 0; i < death_t.size(); ++i)',
        '    ll += nbinom_gammaconf(death_v[i], max(0.1, dyn("death_o", death_t[i], {}, {}) + dyn("death2_o", death_t[i], {}, {})), 10, death_tr[i], 4.5, 2.7);',
        'for (unsigned int i = 0; i < hosp_i_t.size(); ++i)',
        '    ll += nbinom(hosp_i_v[i], max(0.1, dyn("hosp_undetected_o", hosp_i_t[i], {}, {}) + dyn("hosp_undetected2_o", hosp_i_t[i], {}, {})), 20);',
        'for (unsigned int i = 0; i < hosp_p_t.size(); ++i)',
        '    ll += nbinom(hosp_p_v[i], max(0.1, dyn("hosp_p", hosp_p_t[i], {}, {}) + dyn("hosp2_p", hosp_p_t[i], {}, {}) - dyn("hosp_undetected_p", hosp_p_t[i], {}, {}) - dyn("hosp_undetected2_p", hosp_p_t[i], {}, {})), 20);',
        'for (unsigned int i = 0; i < icu_p_t.size(); ++i)',
        '    ll += nbinom(icu_p_v[i], max(0.1, dyn("icu_p", icu_p_t[i], {}, {}) + dyn("icu2_p", icu_p_t[i], {}, {})), 20);',
        'for (unsigned int i = 0; i < virus_t0.size(); ++i) {',
        '    double virus_prev = 0;',
        '    for (double tt = virus_t0[i]; tt <= virus_t1[i]; ++tt)',
        '        virus_prev += dyn("pcr_positive_p", tt, {}, {}) / (virus_t1[i] - virus_t0[i] + 1);',
        '    ll += 10 * skewnorm(virus_prev / ${ sum(params$pop[[1]]$size) }, virus_xi[i], virus_om[i], virus_al[i]);',
        '}',
        'for (unsigned int i = 0; i < sero_t0.size(); ++i) {',
        '    double sero_prev = 0;',
        '    for (double tt = sero_t0[i]; tt <= sero_t1[i]; ++tt)',
        '        sero_prev += dyn("lfia_positive_p", tt, {}, {}) / (sero_t1[i] - sero_t0[i] + 1);',
        '    ll += skewnorm(sero_prev / ${ sum(params$pop[[1]]$size) }, sero_xi[i], sero_om[i], sero_al[i]);',
        '}',
        'for (unsigned int i = 0; i < sgtf_t.size(); ++i) {',
        '    double s1 = dyn("test_o",  sgtf_t[i], {}, {});',
        '    double s2 = dyn("test2_o", sgtf_t[i], {}, {});',
        '    double p2 = (s1 + s2) > 0 ? s2 / (s1 + s2) : 0;',
        '    double model_sgtf = p2 + (1 - p2) * x[23];',
        '    ll += bbinom(sgtf_s[i], sgtf_s[i] + sgtf_o[i], model_sgtf, x[22]);',
        '}',
        
        .sep = "\n    ", .open = "${", .close = "}"
    )
}


# Observer function (basic)
cpp_obsI = function(P.death)
{
    c(
        'auto asc = [&](double x, double y0, double y1, double s0, double s1) {',
        '    double xx = s0 + x * (s1 - s0);',
        '    double h0 = exp(s0) / (1 + exp(s0));',
        '    double h1 = exp(s1) / (1 + exp(s1));',
        '    double h = (exp(xx) / (1 + exp(xx)) - h0) / (h1 - h0);',
        '    return y0 + (y1 - y0) * h;',
        '};',

        'dyn.Obs(t, 0, 0, 0) = estimate_Rt(P, dyn, t, 0, 50);',
        'dyn.Obs(t, 0, 3, 0) = estimate_R0(P, dyn, t, 0, 50);',

        # increase in young person mobility
        'if (t == 182) {',
        '    double mode = 0.2;',
        '    double conc = x[12];',
        '    double constant = 0.2;',
        '    for (unsigned int a = 0; a < P.pop[0].u.size(); ++a) {',
        '        P.pop[0].u[a]  *= (1 - constant) * dbeta((a + 0.5) / P.pop[0].u.size(), mode * (conc - 2) + 1, (1 - mode) * (conc - 2) + 1) + constant;',
        '        P.pop[0].u2[a] *= (1 - constant) * dbeta((a + 0.5) / P.pop[0].u.size(), mode * (conc - 2) + 1, (1 - mode) * (conc - 2) + 1) + constant;',
        '    }',
        '}',
        'if (t == 213) {',
        '    double mode = 0.2;',
        '    double conc_prev = x[12];',
        '    double conc = x[13];',
        '    double constant = 0.2;',
        '    for (unsigned int a = 0; a < P.pop[0].u.size(); ++a) {',
        '        P.pop[0].u[a]  /= (1 - constant) * dbeta((a + 0.5) / P.pop[0].u.size(), mode * (conc_prev - 2) + 1, (1 - mode) * (conc_prev - 2) + 1) + constant;',
        '        P.pop[0].u[a]  *= (1 - constant) * dbeta((a + 0.5) / P.pop[0].u.size(), mode * (conc - 2) + 1, (1 - mode) * (conc - 2) + 1) + constant;',
        '        P.pop[0].u2[a] /= (1 - constant) * dbeta((a + 0.5) / P.pop[0].u.size(), mode * (conc_prev - 2) + 1, (1 - mode) * (conc_prev - 2) + 1) + constant;',
        '        P.pop[0].u2[a] *= (1 - constant) * dbeta((a + 0.5) / P.pop[0].u.size(), mode * (conc - 2) + 1, (1 - mode) * (conc - 2) + 1) + constant;',
        '    }',
        '}',
        'if (t == 244) {',
        '    double mode = 0.2;',
        '    double conc_prev = x[13];',
        '    double conc = x[14];',
        '    double constant = 0.2;',
        '    for (unsigned int a = 0; a < P.pop[0].u.size(); ++a) {',
        '        P.pop[0].u[a]  /= (1 - constant) * dbeta((a + 0.5) / P.pop[0].u.size(), mode * (conc_prev - 2) + 1, (1 - mode) * (conc_prev - 2) + 1) + constant;',
        '        P.pop[0].u[a]  *= (1 - constant) * dbeta((a + 0.5) / P.pop[0].u.size(), mode * (conc - 2) + 1, (1 - mode) * (conc - 2) + 1) + constant;',
        '        P.pop[0].u2[a] /= (1 - constant) * dbeta((a + 0.5) / P.pop[0].u.size(), mode * (conc_prev - 2) + 1, (1 - mode) * (conc_prev - 2) + 1) + constant;',
        '        P.pop[0].u2[a] *= (1 - constant) * dbeta((a + 0.5) / P.pop[0].u.size(), mode * (conc - 2) + 1, (1 - mode) * (conc - 2) + 1) + constant;',
        '    }',
        '}',
        # changing CFR
        'if ((int)t % 7 == 0) {',
        '    double ifr_table[] = {', paste(P.death, collapse = ", "), '};',
        '    for (unsigned int g = 0; g < P.processes[0].prob.size(); ++g) {',
        '        double ifr = ifr_table[g];',
        '        P.processes[4].prob[g][0]  = asc(min(365.0, t/365.0), ifr * x[5], ifr * x[5] * x[16], -2.9, 7.8);',
        '        P.processes[4].prob[g][1]  = 1 - P.processes[4].prob[g][0];',
        '        P.processes[14].prob[g][0] = asc(min(365.0, t/365.0), ifr * x[5] * x[21], ifr * x[5] * x[16] * x[21], -2.9, 7.8);',
        '        P.processes[14].prob[g][1] = 1 - P.processes[14].prob[g][0];',
        '    }',
        '}',
        # changing detection rate in patients admitted to hospital
        'double detection = asc(min(t / 365.0, 1.0), 14, 1, -5.86, 33.4);',
        'P.processes[7].delays[0]  = delay_gamma(detection, 0.59, 60, 0.25);',
        'P.processes[15].delays[0] = delay_gamma(detection, 0.59, 60, 0.25);'
    )
}

# create cpp changes - for projections (i.e. adding SE of tiers effects)
cpp_chgI_proj = function(cb_days, se, nx)
{
    glue::glue(
        'vector<double> work_curve = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.008, 0.021, 0.033, 0.046, 0.058, 0.071, 0.083, 0.096, 0.108, 0.121, 0.133, 0.146, 0.158, 0.171, 0.183, 0.196, 0.208, 0.221, 0.233, 0.246, 0.258, 0.271, 0.283, 0.296, 0.308, 0.321, 0.334, 0.346, 0.359, 0.371, 0.384, 0.397, 0.41, 0.422, 0.435, 0.448, 0.461, 0.474, 0.487, 0.5, 0.513, 0.526, 0.539, 0.552, 0.566, 0.579, 0.592, 0.606, 0.619, 0.633, 0.646, 0.66, 0.674, 0.687, 0.701, 0.715, 0.729, 0.743, 0.757, 0.771, 0.785, 0.799, 0.813, 0.828, 0.842, 0.856, 0.87, 0.885, 0.899, 0.914, 0.928, 0.942, 0.957, 0.971, 0.986, 1, 1.014, 1.029, 1.043, 1.058, 1.072, 1.087, 1.101, 1.115, 1.13, 1.144, 1.159, 1.173, 1.188, 1.202, 1.216, 1.231, 1.245, 1.26, 1.274, 1.289, 1.303, 1.317, 1.332, 1.346, 1.361 };',
        'vector<double> other_curve = { 0.064, 0.066, 0.067, 0.068, 0.069, 0.071, 0.072, 0.073, 0.075, 0.076, 0.077, 0.078, 0.08, 0.081, 0.082, 0.084, 0.085, 0.086, 0.087, 0.089, 0.09, 0.091, 0.092, 0.094, 0.095, 0.096, 0.098, 0.099, 0.1, 0.101, 0.103, 0.104, 0.105, 0.106, 0.108, 0.109, 0.11, 0.112, 0.113, 0.114, 0.116, 0.118, 0.119, 0.121, 0.123, 0.125, 0.128, 0.13, 0.132, 0.135, 0.137, 0.14, 0.143, 0.146, 0.15, 0.154, 0.159, 0.164, 0.169, 0.175, 0.182, 0.19, 0.198, 0.207, 0.217, 0.228, 0.24, 0.252, 0.266, 0.28, 0.295, 0.31, 0.327, 0.344, 0.361, 0.379, 0.398, 0.418, 0.438, 0.459, 0.48, 0.502, 0.525, 0.549, 0.572, 0.597, 0.621, 0.647, 0.672, 0.698, 0.725, 0.751, 0.778, 0.805, 0.833, 0.86, 0.888, 0.916, 0.944, 0.972, 1, 1.028, 1.056, 1.084, 1.112, 1.14, 1.168, 1.196, 1.224, 1.252, 1.28, 1.308, 1.337, 1.365, 1.393, 1.421, 1.449, 1.477, 1.505, 1.533, 1.561, 1.589, 1.617, 1.645, 1.673, 1.701 };',
        'vector<double> cb_days = ${ cpp_vec(cb_days) };',
        'vector<double> se = ${ cpp_vec(se) };',
        'auto interp = [&](double x, vector<double>& curve) {',
        '    if (x < 0) return curve[0];',
        '    unsigned int i = (unsigned int)(x * 100);',
        '    if (i >= curve.size()) return curve.back();',
        '    double f = x * 100 - i;',
        '    return (1 - f) * curve[i] + f * curve[i + 1];',
        '};',
        
        'auto odds = [&](double x, double lo) {',
        '    double a = x / (1 - x);',
        '    return a * exp(lo) / (a * exp(lo) + 1);',
        '};',

        'for (unsigned int g = 0; g < P.processes[0].prob.size(); ++g) {',
        '    // To hospital',
        '    P.processes[0].prob[g][0]  = odds(P.processes[0].prob[g][0], x[7]);',
        '    P.processes[0].prob[g][1]  = 1 - P.processes[0].prob[g][0];',
        '    P.processes[10].prob[g][0] = odds(P.processes[10].prob[g][0], x[7] + x[20]);',
        '    P.processes[10].prob[g][1] = 1 - P.processes[10].prob[g][0];',
        '    // To ICU',
        '    P.processes[1].prob[g][0]  = odds(P.processes[1].prob[g][0], x[7] + x[6]);',
        '    P.processes[1].prob[g][1]  = 1 - P.processes[1].prob[g][0];',
        '    P.processes[11].prob[g][0] = odds(P.processes[11].prob[g][0], x[7] + x[6] + x[20]);',
        '    P.processes[11].prob[g][1] = 1 - P.processes[11].prob[g][0];',
        '    // To death',
        '    P.processes[4].prob[g][0]  = min(1.0, P.processes[4].prob[g][0] * x[5]);',
        '    P.processes[4].prob[g][1]  = 1 - P.processes[4].prob[g][0];',
        '    P.processes[14].prob[g][0] = min(1.0, P.processes[14].prob[g][0] * x[5] * x[21]);',
        '    P.processes[14].prob[g][1] = 1 - P.processes[14].prob[g][0];',
        '    // Relative susceptibility',
        '    P.pop[0].u[g]  = P.pop[0].u[g]  * x[1];',
        '    P.pop[0].u2[g] = P.pop[0].u2[g] * x[1] * x[19];',
        '}',

        '// Delays to admission to hospital & ICU',
        'P.processes[0]. delays[0] = delay_gamma(x[4], 0.71, 60, 0.25);',
        'P.processes[10].delays[0] = delay_gamma(x[4], 0.71, 60, 0.25);',
        'P.processes[1]. delays[0] = delay_gamma(x[8], 1.91, 60, 0.25);',
        'P.processes[11].delays[0] = delay_gamma(x[8], 1.91, 60, 0.25);',
        '// Death',
        'P.processes[4] .delays[0] = delay_gamma(x[2], x[3], 60, 0.25);',
        'P.processes[14].delays[0] = delay_gamma(x[2], x[3], 60, 0.25);',
        '// Seeding of original and variant strain',
        'P.pop[0].seed_times = seq((int)x[0], (int)x[0] + 27);',
        'P.pop[0].seed_times2 = vector<double>(10, x[18]);',

        'auto asc = [&](double x, double y0, double y1, double s0, double s1) {',
        '    double xx = s0 + x * (s1 - s0);',
        '    double h0 = exp(s0) / (1 + exp(s0));',
        '    double h1 = exp(s1) / (1 + exp(s1));',
        '    double h = (exp(xx) / (1 + exp(xx)) - h0) / (h1 - h0);',
        '    return y0 + (y1 - y0) * h;',
        '};',

        # Contact adjustment
        'for (unsigned int i = 0; i < P.changes.ch[4].times.size(); ++i) {',
        '    double tx = double(i) / (P.changes.ch[4].times.size() - 1.0);',
        '    P.changes.ch[4].values[i] = vector<double>(8, asc(tx, 1.0, x[9], -x[10], x[11]));',
        '}',

        # Sep boost
        'P.changes.ch[5].values[0] = vector<double>(8, x[15]);',
        'P.changes.ch[5].times[0] = x[17];',

        # Processing of google mobility indices into contact rates. 0 = regular/lockdown, 2 = tier 2, 3 = tier 3
        'unsigned int stdnorm_index = 0;',
        'for (unsigned int k : vector<unsigned int> { 0, 2, 3 }) {',
        '    for (unsigned int i = 0; i < P.changes.ch[k].times.size(); ++i) {',
        '        //double resi = P.changes.ch[k].values[i][0];',
        '        double wplc = P.changes.ch[k].values[i][1];',
        '        double groc = P.changes.ch[k].values[i][2];',
        '        double rtrc = P.changes.ch[k].values[i][3];',
        '        double trns = P.changes.ch[k].values[i][4];',
        '        if ((k == 0 && cb_days[i] == 1) || k > 0)',
        '        {',
        '            wplc += x[stdnorm_index + ${nx}] * se[stdnorm_index * 4 + 0] / 100;',
        '            groc += x[stdnorm_index + ${nx}] * se[stdnorm_index * 4 + 1] / 100;',
        '            rtrc += x[stdnorm_index + ${nx}] * se[stdnorm_index * 4 + 2] / 100;',
        '            trns += x[stdnorm_index + ${nx}] * se[stdnorm_index * 4 + 3] / 100;',
        '        }',
        '        double othx = rtrc * 0.345 + trns * 0.445 + groc * 0.210;',
        '        double t = P.changes.ch[k].times[i];',

        # from CoMix analysis
        '        double home = asc(min(1.0, t / 365.0), 1.0, 1.545019 / 3.875622, -79 * 0.6, 286 * 0.6);',
        '        double work = interp(wplc, work_curve);',
        '        double scho = (t >= 81 && t <= 244) ? 0 : 1;',
        '        double othe = interp(othx, other_curve);',
        '        P.changes.ch[k].values[i] = { home, work, scho, othe, home, work, scho, othe };',
        '    }',
        '    ++stdnorm_index;',
        '}',
        
        # old... for fIs changes - overwriting changes in fIs here.
        'for (unsigned int i = 0; i < P.changes.ch[1].times.size(); ++i) {',
        '    P.changes.ch[1].values[i] = vector<double>(16, 1.0);',
        '}',
    .sep = "\n", .open = "${", .close = "}")
}




# Immune escape version
# create cpp changes
cpp_chgI_ImmEsc = function()
{
    c(
        'vector<double> work_curve = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.008, 0.021, 0.033, 0.046, 0.058, 0.071, 0.083, 0.096, 0.108, 0.121, 0.133, 0.146, 0.158, 0.171, 0.183, 0.196, 0.208, 0.221, 0.233, 0.246, 0.258, 0.271, 0.283, 0.296, 0.308, 0.321, 0.334, 0.346, 0.359, 0.371, 0.384, 0.397, 0.41, 0.422, 0.435, 0.448, 0.461, 0.474, 0.487, 0.5, 0.513, 0.526, 0.539, 0.552, 0.566, 0.579, 0.592, 0.606, 0.619, 0.633, 0.646, 0.66, 0.674, 0.687, 0.701, 0.715, 0.729, 0.743, 0.757, 0.771, 0.785, 0.799, 0.813, 0.828, 0.842, 0.856, 0.87, 0.885, 0.899, 0.914, 0.928, 0.942, 0.957, 0.971, 0.986, 1, 1.014, 1.029, 1.043, 1.058, 1.072, 1.087, 1.101, 1.115, 1.13, 1.144, 1.159, 1.173, 1.188, 1.202, 1.216, 1.231, 1.245, 1.26, 1.274, 1.289, 1.303, 1.317, 1.332, 1.346, 1.361 };',
        'vector<double> other_curve = { 0.064, 0.066, 0.067, 0.068, 0.069, 0.071, 0.072, 0.073, 0.075, 0.076, 0.077, 0.078, 0.08, 0.081, 0.082, 0.084, 0.085, 0.086, 0.087, 0.089, 0.09, 0.091, 0.092, 0.094, 0.095, 0.096, 0.098, 0.099, 0.1, 0.101, 0.103, 0.104, 0.105, 0.106, 0.108, 0.109, 0.11, 0.112, 0.113, 0.114, 0.116, 0.118, 0.119, 0.121, 0.123, 0.125, 0.128, 0.13, 0.132, 0.135, 0.137, 0.14, 0.143, 0.146, 0.15, 0.154, 0.159, 0.164, 0.169, 0.175, 0.182, 0.19, 0.198, 0.207, 0.217, 0.228, 0.24, 0.252, 0.266, 0.28, 0.295, 0.31, 0.327, 0.344, 0.361, 0.379, 0.398, 0.418, 0.438, 0.459, 0.48, 0.502, 0.525, 0.549, 0.572, 0.597, 0.621, 0.647, 0.672, 0.698, 0.725, 0.751, 0.778, 0.805, 0.833, 0.86, 0.888, 0.916, 0.944, 0.972, 1, 1.028, 1.056, 1.084, 1.112, 1.14, 1.168, 1.196, 1.224, 1.252, 1.28, 1.308, 1.337, 1.365, 1.393, 1.421, 1.449, 1.477, 1.505, 1.533, 1.561, 1.589, 1.617, 1.645, 1.673, 1.701 };',
        'auto interp = [&](double x, vector<double>& curve) {',
        '    if (x < 0) return curve[0];',
        '    unsigned int i = (unsigned int)(x * 100);',
        '    if (i >= curve.size()) return curve.back();',
        '    double f = x * 100 - i;',
        '    return (1 - f) * curve[i] + f * curve[i + 1];',
        '};',
        
        'auto odds = [&](double x, double lo) {',
        '    double a = x / (1 - x);',
        '    return a * exp(lo) / (a * exp(lo) + 1);',
        '};',

        'for (unsigned int g = 0; g < P.processes[0].prob.size(); ++g) {',
        '    // To hospital',
        '    P.processes[0].prob[g][0]  = odds(P.processes[0].prob[g][0], x[7]);',
        '    P.processes[0].prob[g][1]  = 1 - P.processes[0].prob[g][0];',
        '    P.processes[10].prob[g][0] = odds(P.processes[10].prob[g][0], x[7] + x[20]);',
        '    P.processes[10].prob[g][1] = 1 - P.processes[10].prob[g][0];',
        '    // To ICU',
        '    P.processes[1].prob[g][0]  = odds(P.processes[1].prob[g][0], x[7] + x[6]);',
        '    P.processes[1].prob[g][1]  = 1 - P.processes[1].prob[g][0];',
        '    P.processes[11].prob[g][0] = odds(P.processes[11].prob[g][0], x[7] + x[6] + x[20]);',
        '    P.processes[11].prob[g][1] = 1 - P.processes[11].prob[g][0];',
        '    // To death',
        '    P.processes[4].prob[g][0]  = min(1.0, P.processes[4].prob[g][0] * x[5]);',
        '    P.processes[4].prob[g][1]  = 1 - P.processes[4].prob[g][0];',
        '    P.processes[14].prob[g][0] = min(1.0, P.processes[14].prob[g][0] * x[5] * x[21]);',
        '    P.processes[14].prob[g][1] = 1 - P.processes[14].prob[g][0];',
        '    // Relative susceptibility',
        '    P.pop[0].u[g]  = P.pop[0].u[g]  * x[1];',
        '    P.pop[0].u2[g] = P.pop[0].u2[g] * x[1];',
        '    // Immune escape',
        '    P.pop[0].pi2_r[g] = x[19];',
        '}',

        '// Delays to admission to hospital & ICU',
        'P.processes[0]. delays[0] = delay_gamma(x[4], 0.71, 60, 0.25);',
        'P.processes[10].delays[0] = delay_gamma(x[4], 0.71, 60, 0.25);',
        'P.processes[1]. delays[0] = delay_gamma(x[8], 1.91, 60, 0.25);',
        'P.processes[11].delays[0] = delay_gamma(x[8], 1.91, 60, 0.25);',
        '// Death',
        'P.processes[4] .delays[0] = delay_gamma(x[2], x[3], 60, 0.25);',
        'P.processes[14].delays[0] = delay_gamma(x[2], x[3], 60, 0.25);',
        '// Seeding of original and variant strain',
        'P.pop[0].seed_times = seq((int)x[0], (int)x[0] + 27);',
        'P.pop[0].seed_times2 = vector<double>(10, x[18]);',

        'auto asc = [&](double x, double y0, double y1, double s0, double s1) {',
        '    double xx = s0 + x * (s1 - s0);',
        '    double h0 = exp(s0) / (1 + exp(s0));',
        '    double h1 = exp(s1) / (1 + exp(s1));',
        '    double h = (exp(xx) / (1 + exp(xx)) - h0) / (h1 - h0);',
        '    return y0 + (y1 - y0) * h;',
        '};',

        # Contact adjustment
        'for (unsigned int i = 0; i < P.changes.ch[4].times.size(); ++i) {',
        '    double tx = double(i) / (P.changes.ch[4].times.size() - 1.0);',
        '    P.changes.ch[4].values[i] = vector<double>(8, asc(tx, 1.0, x[9], -x[10], x[11]));',
        '}',

        # Sep boost
        'P.changes.ch[5].values[0] = vector<double>(8, x[15]);',
        'P.changes.ch[5].times[0] = x[17];',

        # fitting of google mobility indices
        'for (unsigned int k : vector<unsigned int> { 0, 2, 3 }) {',
        '    for (unsigned int i = 0; i < P.changes.ch[k].times.size(); ++i) {',
        '        //double resi = P.changes.ch[k].values[i][0];',
        '        double wplc = P.changes.ch[k].values[i][1];',
        '        double groc = P.changes.ch[k].values[i][2];',
        '        double rtrc = P.changes.ch[k].values[i][3];',
        '        double trns = P.changes.ch[k].values[i][4];',
        '        double othx = rtrc * 0.345 + trns * 0.445 + groc * 0.210;',
        '        double t = P.changes.ch[k].times[i];',

        # from CoMix analysis
        '        double home = asc(min(1.0, t / 365.0), 1.0, 1.545019 / 3.875622, -79 * 0.6, 286 * 0.6);',
        '        double work = interp(wplc, work_curve);',
        '        double scho = (t >= 81 && t <= 244) ? 0 : 1;',
        '        double othe = interp(othx, other_curve);',
        '        P.changes.ch[k].values[i] = { home, work, scho, othe, home, work, scho, othe };',
        '    }',
        '}',
        
        # old... for fIs changes - overwriting changes in fIs here.
        'for (unsigned int i = 0; i < P.changes.ch[1].times.size(); ++i) {',
        '    P.changes.ch[1].values[i] = vector<double>(16, 1.0);',
        '}'
    )
}



# Latent period version
# create cpp changes
cpp_chgI_LatentPeriod = function()
{
    c(
        'vector<double> work_curve = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.008, 0.021, 0.033, 0.046, 0.058, 0.071, 0.083, 0.096, 0.108, 0.121, 0.133, 0.146, 0.158, 0.171, 0.183, 0.196, 0.208, 0.221, 0.233, 0.246, 0.258, 0.271, 0.283, 0.296, 0.308, 0.321, 0.334, 0.346, 0.359, 0.371, 0.384, 0.397, 0.41, 0.422, 0.435, 0.448, 0.461, 0.474, 0.487, 0.5, 0.513, 0.526, 0.539, 0.552, 0.566, 0.579, 0.592, 0.606, 0.619, 0.633, 0.646, 0.66, 0.674, 0.687, 0.701, 0.715, 0.729, 0.743, 0.757, 0.771, 0.785, 0.799, 0.813, 0.828, 0.842, 0.856, 0.87, 0.885, 0.899, 0.914, 0.928, 0.942, 0.957, 0.971, 0.986, 1, 1.014, 1.029, 1.043, 1.058, 1.072, 1.087, 1.101, 1.115, 1.13, 1.144, 1.159, 1.173, 1.188, 1.202, 1.216, 1.231, 1.245, 1.26, 1.274, 1.289, 1.303, 1.317, 1.332, 1.346, 1.361 };',
        'vector<double> other_curve = { 0.064, 0.066, 0.067, 0.068, 0.069, 0.071, 0.072, 0.073, 0.075, 0.076, 0.077, 0.078, 0.08, 0.081, 0.082, 0.084, 0.085, 0.086, 0.087, 0.089, 0.09, 0.091, 0.092, 0.094, 0.095, 0.096, 0.098, 0.099, 0.1, 0.101, 0.103, 0.104, 0.105, 0.106, 0.108, 0.109, 0.11, 0.112, 0.113, 0.114, 0.116, 0.118, 0.119, 0.121, 0.123, 0.125, 0.128, 0.13, 0.132, 0.135, 0.137, 0.14, 0.143, 0.146, 0.15, 0.154, 0.159, 0.164, 0.169, 0.175, 0.182, 0.19, 0.198, 0.207, 0.217, 0.228, 0.24, 0.252, 0.266, 0.28, 0.295, 0.31, 0.327, 0.344, 0.361, 0.379, 0.398, 0.418, 0.438, 0.459, 0.48, 0.502, 0.525, 0.549, 0.572, 0.597, 0.621, 0.647, 0.672, 0.698, 0.725, 0.751, 0.778, 0.805, 0.833, 0.86, 0.888, 0.916, 0.944, 0.972, 1, 1.028, 1.056, 1.084, 1.112, 1.14, 1.168, 1.196, 1.224, 1.252, 1.28, 1.308, 1.337, 1.365, 1.393, 1.421, 1.449, 1.477, 1.505, 1.533, 1.561, 1.589, 1.617, 1.645, 1.673, 1.701 };',
        'auto interp = [&](double x, vector<double>& curve) {',
        '    if (x < 0) return curve[0];',
        '    unsigned int i = (unsigned int)(x * 100);',
        '    if (i >= curve.size()) return curve.back();',
        '    double f = x * 100 - i;',
        '    return (1 - f) * curve[i] + f * curve[i + 1];',
        '};',
        
        'auto odds = [&](double x, double lo) {',
        '    double a = x / (1 - x);',
        '    return a * exp(lo) / (a * exp(lo) + 1);',
        '};',

        'for (unsigned int g = 0; g < P.processes[0].prob.size(); ++g) {',
        '    // To hospital',
        '    P.processes[0].prob[g][0]  = odds(P.processes[0].prob[g][0], x[7]);',
        '    P.processes[0].prob[g][1]  = 1 - P.processes[0].prob[g][0];',
        '    P.processes[10].prob[g][0] = odds(P.processes[10].prob[g][0], x[7] + x[20]);',
        '    P.processes[10].prob[g][1] = 1 - P.processes[10].prob[g][0];',
        '    // To ICU',
        '    P.processes[1].prob[g][0]  = odds(P.processes[1].prob[g][0], x[7] + x[6]);',
        '    P.processes[1].prob[g][1]  = 1 - P.processes[1].prob[g][0];',
        '    P.processes[11].prob[g][0] = odds(P.processes[11].prob[g][0], x[7] + x[6] + x[20]);',
        '    P.processes[11].prob[g][1] = 1 - P.processes[11].prob[g][0];',
        '    // To death',
        '    P.processes[4].prob[g][0]  = min(1.0, P.processes[4].prob[g][0] * x[5]);',
        '    P.processes[4].prob[g][1]  = 1 - P.processes[4].prob[g][0];',
        '    P.processes[14].prob[g][0] = min(1.0, P.processes[14].prob[g][0] * x[5] * x[21]);',
        '    P.processes[14].prob[g][1] = 1 - P.processes[14].prob[g][0];',
        '    // Relative susceptibility',
        '    P.pop[0].u[g]  = P.pop[0].u[g]  * x[1];',
        '    P.pop[0].u2[g] = P.pop[0].u2[g] * x[1];',
        '}',
        
        '// Latent period',
        'P.pop[0].dE2 = delay_gamma(x[19] * 2.5, 2.5, 15, 0.25);',

        '// Delays to admission to hospital & ICU',
        'P.processes[0]. delays[0] = delay_gamma(x[4], 0.71, 60, 0.25);',
        'P.processes[10].delays[0] = delay_gamma(x[4], 0.71, 60, 0.25);',
        'P.processes[1]. delays[0] = delay_gamma(x[8], 1.91, 60, 0.25);',
        'P.processes[11].delays[0] = delay_gamma(x[8], 1.91, 60, 0.25);',
        '// Death',
        'P.processes[4] .delays[0] = delay_gamma(x[2], x[3], 60, 0.25);',
        'P.processes[14].delays[0] = delay_gamma(x[2], x[3], 60, 0.25);',
        '// Seeding of original and variant strain',
        'P.pop[0].seed_times = seq((int)x[0], (int)x[0] + 27);',
        'P.pop[0].seed_times2 = vector<double>(10, x[18]);',

        'auto asc = [&](double x, double y0, double y1, double s0, double s1) {',
        '    double xx = s0 + x * (s1 - s0);',
        '    double h0 = exp(s0) / (1 + exp(s0));',
        '    double h1 = exp(s1) / (1 + exp(s1));',
        '    double h = (exp(xx) / (1 + exp(xx)) - h0) / (h1 - h0);',
        '    return y0 + (y1 - y0) * h;',
        '};',

        # Contact adjustment
        'for (unsigned int i = 0; i < P.changes.ch[4].times.size(); ++i) {',
        '    double tx = double(i) / (P.changes.ch[4].times.size() - 1.0);',
        '    P.changes.ch[4].values[i] = vector<double>(8, asc(tx, 1.0, x[9], -x[10], x[11]));',
        '}',

        # Sep boost
        'P.changes.ch[5].values[0] = vector<double>(8, x[15]);',
        'P.changes.ch[5].times[0] = x[17];',

        # fitting of google mobility indices
        'for (unsigned int k : vector<unsigned int> { 0, 2, 3 }) {',
        '    for (unsigned int i = 0; i < P.changes.ch[k].times.size(); ++i) {',
        '        //double resi = P.changes.ch[k].values[i][0];',
        '        double wplc = P.changes.ch[k].values[i][1];',
        '        double groc = P.changes.ch[k].values[i][2];',
        '        double rtrc = P.changes.ch[k].values[i][3];',
        '        double trns = P.changes.ch[k].values[i][4];',
        '        double othx = rtrc * 0.345 + trns * 0.445 + groc * 0.210;',
        '        double t = P.changes.ch[k].times[i];',

        # from CoMix analysis
        '        double home = asc(min(1.0, t / 365.0), 1.0, 1.545019 / 3.875622, -79 * 0.6, 286 * 0.6);',
        '        double work = interp(wplc, work_curve);',
        '        double scho = (t >= 81 && t <= 244) ? 0 : 1;',
        '        double othe = interp(othx, other_curve);',
        '        P.changes.ch[k].values[i] = { home, work, scho, othe, home, work, scho, othe };',
        '    }',
        '}',
        
        # old... for fIs changes - overwriting changes in fIs here.
        'for (unsigned int i = 0; i < P.changes.ch[1].times.size(); ++i) {',
        '    P.changes.ch[1].values[i] = vector<double>(16, 1.0);',
        '}'
    )
}


# Children susceptibility version
# create cpp changes
cpp_chgI_Kids = function()
{
    c(
        'vector<double> work_curve = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.008, 0.021, 0.033, 0.046, 0.058, 0.071, 0.083, 0.096, 0.108, 0.121, 0.133, 0.146, 0.158, 0.171, 0.183, 0.196, 0.208, 0.221, 0.233, 0.246, 0.258, 0.271, 0.283, 0.296, 0.308, 0.321, 0.334, 0.346, 0.359, 0.371, 0.384, 0.397, 0.41, 0.422, 0.435, 0.448, 0.461, 0.474, 0.487, 0.5, 0.513, 0.526, 0.539, 0.552, 0.566, 0.579, 0.592, 0.606, 0.619, 0.633, 0.646, 0.66, 0.674, 0.687, 0.701, 0.715, 0.729, 0.743, 0.757, 0.771, 0.785, 0.799, 0.813, 0.828, 0.842, 0.856, 0.87, 0.885, 0.899, 0.914, 0.928, 0.942, 0.957, 0.971, 0.986, 1, 1.014, 1.029, 1.043, 1.058, 1.072, 1.087, 1.101, 1.115, 1.13, 1.144, 1.159, 1.173, 1.188, 1.202, 1.216, 1.231, 1.245, 1.26, 1.274, 1.289, 1.303, 1.317, 1.332, 1.346, 1.361 };',
        'vector<double> other_curve = { 0.064, 0.066, 0.067, 0.068, 0.069, 0.071, 0.072, 0.073, 0.075, 0.076, 0.077, 0.078, 0.08, 0.081, 0.082, 0.084, 0.085, 0.086, 0.087, 0.089, 0.09, 0.091, 0.092, 0.094, 0.095, 0.096, 0.098, 0.099, 0.1, 0.101, 0.103, 0.104, 0.105, 0.106, 0.108, 0.109, 0.11, 0.112, 0.113, 0.114, 0.116, 0.118, 0.119, 0.121, 0.123, 0.125, 0.128, 0.13, 0.132, 0.135, 0.137, 0.14, 0.143, 0.146, 0.15, 0.154, 0.159, 0.164, 0.169, 0.175, 0.182, 0.19, 0.198, 0.207, 0.217, 0.228, 0.24, 0.252, 0.266, 0.28, 0.295, 0.31, 0.327, 0.344, 0.361, 0.379, 0.398, 0.418, 0.438, 0.459, 0.48, 0.502, 0.525, 0.549, 0.572, 0.597, 0.621, 0.647, 0.672, 0.698, 0.725, 0.751, 0.778, 0.805, 0.833, 0.86, 0.888, 0.916, 0.944, 0.972, 1, 1.028, 1.056, 1.084, 1.112, 1.14, 1.168, 1.196, 1.224, 1.252, 1.28, 1.308, 1.337, 1.365, 1.393, 1.421, 1.449, 1.477, 1.505, 1.533, 1.561, 1.589, 1.617, 1.645, 1.673, 1.701 };',
        'auto interp = [&](double x, vector<double>& curve) {',
        '    if (x < 0) return curve[0];',
        '    unsigned int i = (unsigned int)(x * 100);',
        '    if (i >= curve.size()) return curve.back();',
        '    double f = x * 100 - i;',
        '    return (1 - f) * curve[i] + f * curve[i + 1];',
        '};',
        
        'auto odds = [&](double x, double lo) {',
        '    double a = x / (1 - x);',
        '    return a * exp(lo) / (a * exp(lo) + 1);',
        '};',

        'for (unsigned int g = 0; g < P.processes[0].prob.size(); ++g) {',
        '    // To hospital',
        '    P.processes[0].prob[g][0]  = odds(P.processes[0].prob[g][0], x[7]);',
        '    P.processes[0].prob[g][1]  = 1 - P.processes[0].prob[g][0];',
        '    P.processes[10].prob[g][0] = odds(P.processes[10].prob[g][0], x[7] + x[20]);',
        '    P.processes[10].prob[g][1] = 1 - P.processes[10].prob[g][0];',
        '    // To ICU',
        '    P.processes[1].prob[g][0]  = odds(P.processes[1].prob[g][0], x[7] + x[6]);',
        '    P.processes[1].prob[g][1]  = 1 - P.processes[1].prob[g][0];',
        '    P.processes[11].prob[g][0] = odds(P.processes[11].prob[g][0], x[7] + x[6] + x[20]);',
        '    P.processes[11].prob[g][1] = 1 - P.processes[11].prob[g][0];',
        '    // To death',
        '    P.processes[4].prob[g][0]  = min(1.0, P.processes[4].prob[g][0] * x[5]);',
        '    P.processes[4].prob[g][1]  = 1 - P.processes[4].prob[g][0];',
        '    P.processes[14].prob[g][0] = min(1.0, P.processes[14].prob[g][0] * x[5] * x[21]);',
        '    P.processes[14].prob[g][1] = 1 - P.processes[14].prob[g][0];',
        '    // Relative susceptibility',
        '    P.pop[0].u[g]  = P.pop[0].u[g]  * x[1];',
        '    P.pop[0].u2[g] = P.pop[0].u2[g] * x[1];',
        '}',
        
        '// Increased susceptibility in children',
        'double umax = 1.15 * x[1];',
        'P.pop[0].u2[0] = P.pop[0].u2[0] + x[19] * (umax - P.pop[0].u2[0]);',
        'P.pop[0].u2[1] = P.pop[0].u2[1] + x[19] * (umax - P.pop[0].u2[1]);',
        'P.pop[0].u2[2] = P.pop[0].u2[2] + x[19] * (umax - P.pop[0].u2[2]);',
        'P.pop[0].u2[3] = P.pop[0].u2[3] + x[19] * (umax - P.pop[0].u2[3]);',
        
        '// Delays to admission to hospital & ICU',
        'P.processes[0]. delays[0] = delay_gamma(x[4], 0.71, 60, 0.25);',
        'P.processes[10].delays[0] = delay_gamma(x[4], 0.71, 60, 0.25);',
        'P.processes[1]. delays[0] = delay_gamma(x[8], 1.91, 60, 0.25);',
        'P.processes[11].delays[0] = delay_gamma(x[8], 1.91, 60, 0.25);',
        '// Death',
        'P.processes[4] .delays[0] = delay_gamma(x[2], x[3], 60, 0.25);',
        'P.processes[14].delays[0] = delay_gamma(x[2], x[3], 60, 0.25);',
        '// Seeding of original and variant strain',
        'P.pop[0].seed_times = seq((int)x[0], (int)x[0] + 27);',
        'P.pop[0].seed_times2 = vector<double>(10, x[18]);',

        'auto asc = [&](double x, double y0, double y1, double s0, double s1) {',
        '    double xx = s0 + x * (s1 - s0);',
        '    double h0 = exp(s0) / (1 + exp(s0));',
        '    double h1 = exp(s1) / (1 + exp(s1));',
        '    double h = (exp(xx) / (1 + exp(xx)) - h0) / (h1 - h0);',
        '    return y0 + (y1 - y0) * h;',
        '};',

        # Contact adjustment
        'for (unsigned int i = 0; i < P.changes.ch[4].times.size(); ++i) {',
        '    double tx = double(i) / (P.changes.ch[4].times.size() - 1.0);',
        '    P.changes.ch[4].values[i] = vector<double>(8, asc(tx, 1.0, x[9], -x[10], x[11]));',
        '}',

        # Sep boost
        'P.changes.ch[5].values[0] = vector<double>(8, x[15]);',
        'P.changes.ch[5].times[0] = x[17];',

        # fitting of google mobility indices
        'for (unsigned int k : vector<unsigned int> { 0, 2, 3 }) {',
        '    for (unsigned int i = 0; i < P.changes.ch[k].times.size(); ++i) {',
        '        //double resi = P.changes.ch[k].values[i][0];',
        '        double wplc = P.changes.ch[k].values[i][1];',
        '        double groc = P.changes.ch[k].values[i][2];',
        '        double rtrc = P.changes.ch[k].values[i][3];',
        '        double trns = P.changes.ch[k].values[i][4];',
        '        double othx = rtrc * 0.345 + trns * 0.445 + groc * 0.210;',
        '        double t = P.changes.ch[k].times[i];',

        # from CoMix analysis
        '        double home = asc(min(1.0, t / 365.0), 1.0, 1.545019 / 3.875622, -79 * 0.6, 286 * 0.6);',
        '        double work = interp(wplc, work_curve);',
        '        double scho = (t >= 81 && t <= 244) ? 0 : 1;',
        '        double othe = interp(othx, other_curve);',
        '        P.changes.ch[k].values[i] = { home, work, scho, othe, home, work, scho, othe };',
        '    }',
        '}',
        
        # old... for fIs changes - overwriting changes in fIs here.
        'for (unsigned int i = 0; i < P.changes.ch[1].times.size(); ++i) {',
        '    P.changes.ch[1].values[i] = vector<double>(16, 1.0);',
        '}'
    )
}


# Observer for vaccination programmes
cpp_obsI_vax = function(params, vacc)
{
    glue::glue(
        'if (t == 0) { dyn.scratch["vaxphase"] = 0; }',
        '{',
        'vector<double> vt  = ${ cpp_vec(as.numeric(ymd(vacc$vt) - ymd(params$date0))) };',
        'vector<double> v   = ${ cpp_vec(unlist(vacc$v)) };',
        'unsigned int phase = dyn.scratch["vaxphase"];',
        'if (vt.size() > phase && t >= vt[phase]) {',
        '    int age_groups = P.pop[0].size.size();',
        '    P.pop[0].v   = vector<double>(v.begin()   + age_groups * phase, v.begin()   + age_groups * (phase + 1));',
        '    dyn.scratch["vaxphase"] = phase + 1;',
        '}',
        'dyn.Obs(t, 0, 5, 0) = dyn.scratch["vaxphase"];',
        '}',
        .sep = "\n", .open = "${", .close = "}")
}


# Wrap Cpp funcs for multi region fitting
wrap_region = function(model_v2, name, region_names)
{
    ret = 'if (false) {}';
    for (i in seq_along(region_names))
    {
        ret = c(ret,
            glue::glue('else if (P.pop[0].name == "{region_names[i]}")'),
            '{',
            model_v2[[i]][[name]],
            '}'
        )
    }
    return (ret);
}
