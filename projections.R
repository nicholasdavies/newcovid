# Run projections for control scenarios
source("./projections_setup.R")
source("./define_tiers.R")

# Get variant 2 quantities for other regions which were not fitted because of insufficient data.
variant_pops = c(1, 3, 9)
variant_posterior = rbindlist(posteriorsI[variant_pops], idcol = "population")

variant_posterior[, mean(v2_when)]

for (p in setdiff(england_pops, variant_pops))
{
    which_rows = sample(nrow(variant_posterior), nrow(posteriorsI[[p]]), replace = TRUE)
    posteriorsI[[p]][, v2_when     := rnorm(.N, 290, 5)]
    posteriorsI[[p]]$v2_relu     = variant_posterior[which_rows, v2_relu]
    posteriorsI[[p]]$v2_hosp_rlo = variant_posterior[which_rows, v2_hosp_rlo]
    posteriorsI[[p]]$v2_cfr_rel  = variant_posterior[which_rows, v2_cfr_rel]
}

# Vaccine targeting
targeting_old_to_young = 
    list(c(0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,1,1)*0.85,
         c(0,0,0,0, 0,0,0,0, 0,0,0,0, 1,1,1,1)*0.85,
         c(0,0,0,0, 0,0,0,0, 0,0,1,1, 1,1,1,1)*0.85,
         c(0,0,0,0, 0,0,0,0, 1,1,1,1, 1,1,1,1)*0.85,
         c(0,0,0,0, 0,0,1,1, 1,1,1,1, 1,1,1,1)*0.85,
         c(0,0,0,0, 1,1,1,1, 1,1,1,1, 1,1,1,1)*0.85,
         c(0,0,1,1, 1,1,1,1, 1,1,1,1, 1,1,1,1)*0.85,
         c(1,1,1,1, 1,1,1,1, 1,1,1,1, 1,1,1,1)*0.85)

targeting_young_to_old = 
    list(c(0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,1,1)*0.85,
         c(0,0,0,0, 0,0,0,0, 0,0,0,0, 1,1,1,1)*0.85,
         c(1,1,0,0, 0,0,0,0, 0,0,0,0, 1,1,1,1)*0.85,
         c(1,1,1,1, 0,0,0,0, 0,0,0,0, 1,1,1,1)*0.85,
         c(1,1,1,1, 1,1,0,0, 0,0,0,0, 1,1,1,1)*0.85,
         c(1,1,1,1, 1,1,1,1, 0,0,0,0, 1,1,1,1)*0.85,
         c(1,1,1,1, 1,1,1,1, 1,1,0,0, 1,1,1,1)*0.85,
         c(1,1,1,1, 1,1,1,1, 1,1,1,1, 1,1,1,1)*0.85)

lockdown_dates = rep("2020-12-26", 12)
lockdown_dates[variant_pops] = "2020-12-20"
lockdown_durations = rep(36, 12)
lockdown_durations[variant_pops] = 42

school_breaks_o = c("2020-12-19", "2021-01-04")
school_breaks_c = c("2020-12-19", "2021-02-01")

proj_00    = project(england_pops, tiers = TRUE, tier2 = tier2, tier3 = tier3, cb_date = NA,             cb_duration = NA,                 lockdown = NA, school_breaks = school_breaks_o)
proj_Eo    = project(england_pops, tiers = TRUE, tier2 = tier2, tier3 = tier3, cb_date = lockdown_dates, cb_duration = lockdown_durations, lockdown = lockdownER, se = seER, school_breaks = school_breaks_o, cb_behaviour = "default")
proj_Ec    = project(england_pops, tiers = TRUE, tier2 = tier2, tier3 = tier3, cb_date = lockdown_dates, cb_duration = lockdown_durations, lockdown = lockdownER, se = seER, school_breaks = school_breaks_c, cb_behaviour = "default")
proj_Ec_Vl = project(england_pops, tiers = TRUE, tier2 = tier2, tier3 = tier3, cb_date = lockdown_dates, cb_duration = lockdown_durations, lockdown = lockdownER, se = seER, school_breaks = school_breaks_c, cb_behaviour = "default",
    vacc = make_vaccine_schedule("2021-01-01", rep( 200000, 365), targeting_old_to_young, england_pops), ei_v = rep(0.6, 16), ed_vi = rep(0.875, 16))
proj_Ec_Vh = project(england_pops, tiers = TRUE, tier2 = tier2, tier3 = tier3, cb_date = lockdown_dates, cb_duration = lockdown_durations, lockdown = lockdownER, se = seER, school_breaks = school_breaks_c, cb_behaviour = "default",
    vacc = make_vaccine_schedule("2021-01-01", rep(2000000, 365), targeting_old_to_young, england_pops), ei_v = rep(0.6, 16), ed_vi = rep(0.875, 16))

tb00 = summarize_projection(proj_00, "2020-12-15", popsize)
fwrite(tb00, "./output/table_00.csv");
tbEo = summarize_projection(proj_Eo, "2020-12-15", popsize)
fwrite(tbEo, "./output/table_Eo.csv");
tbEc = summarize_projection(proj_Ec, "2020-12-15", popsize)
fwrite(tbEc, "./output/table_Ec.csv");
tbEcVl = summarize_projection(proj_Ec_Vl, "2020-12-15", popsize)
fwrite(tbEcVl, "./output/table_EcVl.csv");
tbEcVh = summarize_projection(proj_Ec_Vh, "2020-12-15", popsize)
fwrite(tbEcVh, "./output/table_EcVh.csv");

tbEngland = england_only(list(tb00, tbEo, tbEc, tbEcVl, tbEcVh), 
    c("Tiers 1-3 only", "Tier 4, schools open", "Tier 4, schools closed", "Plus 200k immunised per week", "Plus 2M immunised per week"))
fwrite(tbEngland, "./output/table_england.csv");

# plot_projection(list(proj_00), list("Nothing"), "2020-10-01", "Dark2") + geom_vline(aes(xintercept = ymd("2020-12-19")))
# plot_projection(list(proj_Eo), list("England lockdown / schools open"), "2020-10-01", "Dark2")

plot_projection_2parts(list(proj_00, proj_Eo, proj_Ec), list("Tiers 1-3 only", "Tier 4, schools open", "Tier 4, schools closed"), 
    "2020-10-01", c("#b66353", "#5679a3", "#31a1b3"), c("East of England", "London", "South East"))
replic = 10
ggsave(paste0("./output/proj_", replic, ".pdf"), width = 34, height = 16, units = "cm", useDingbats = FALSE)
ggsave(paste0("./output/proj_", replic, ".png"), width = 34, height = 16, units = "cm")

plot_projection_2parts(list(proj_Ec, proj_Ec_Vl, proj_Ec_Vh), list("Tier 4, schools closed", "Plus 200k immunised per week", "Plus 2M immunised per week"), 
    "2020-10-01", c("#31a1b3", "#f47942", "#59a14f"), c("East of England", "London", "South East"))
replic = 10
ggsave(paste0("./output/proj_vacc_", replic, ".pdf"), width = 34, height = 16, units = "cm", useDingbats = FALSE)
ggsave(paste0("./output/proj_vacc_", replic, ".png"), width = 34, height = 16, units = "cm")

# TO HERE - replace Table 1, Fig 4, Fig 5


# New school closure analyses
primary = c(0.0, 1.0, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
sndary0 = c(0.0, 0.0, 0.8, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
sndary1 = c(0.0, 0.0, 0.0, 0.4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
univers = c(0.0, 0.0, 0.0, 0.4, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
teacher = c(1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    
dsc_Sa = 1.0 * primary + 0.5 * sndary0 + 1.0 * sndary1 + 0.0 * univers + 1.0 * teacher;
dsc_Sb = 0.5 * primary + 0.5 * sndary0 + 1.0 * sndary1 + 0.0 * univers + 1.0 * teacher;

school_breaks_2 = c("2020-12-19", "2021-01-03", "2021-01-04", "2021-04-01")


proj_SO = project(england_pops, tiers = TRUE, tier2 = tier2, tier3 = tier3, cb_date = lockdown_dates, cb_duration = lockdown_durations, lockdown = lockdownER, se = seER, cb_behaviour = "default",
    school_breaks = school_breaks_2, school_factors_r = c(0.0, 1.0), school_factors_c = c(0.0, 1.0),
    d_school_date0 = "2021-01-04", d_school_date1 = "2021-04-01", d_school_contact = rep(1, 16))
proj_Sa = project(england_pops, tiers = TRUE, tier2 = tier2, tier3 = tier3, cb_date = lockdown_dates, cb_duration = lockdown_durations, lockdown = lockdownER, se = seER, cb_behaviour = "default",
    school_breaks = school_breaks_2, school_factors_r = c(0.0, 0.79), school_factors_c = c(0.0, 1.0),
    d_school_date0 = "2021-01-04", d_school_date1 = "2021-04-01", d_school_contact = dsc_Sa)
proj_Sb = project(england_pops, tiers = TRUE, tier2 = tier2, tier3 = tier3, cb_date = lockdown_dates, cb_duration = lockdown_durations, lockdown = lockdownER, se = seER, cb_behaviour = "default",
    school_breaks = school_breaks_2, school_factors_r = c(0.0, 0.54), school_factors_c = c(0.0, 1.0),
    d_school_date0 = "2021-01-04", d_school_date1 = "2021-04-01", d_school_contact = dsc_Sb)
proj_SC = project(england_pops, tiers = TRUE, tier2 = tier2, tier3 = tier3, cb_date = lockdown_dates, cb_duration = lockdown_durations, lockdown = lockdownER, se = seER, cb_behaviour = "default",
    school_breaks = school_breaks_2, school_factors_r = c(0.0, 0.0), school_factors_c = c(0.0, 1.0),
    d_school_date0 = "2021-01-04", d_school_date1 = "2021-04-01", d_school_contact = rep(0, 16))


plll = plot_projection_2parts(list(proj_SO, proj_Sa, proj_Sb, proj_SC), list("All open", "M1: P open, SK open, SN alt, U closed", "M2: P alt, SK open, SN alt, U closed", "All closed"), 
    "2020-10-01", c("#5679a3", "#f47942", "#59a14f", "#31a1b3"), c("East of England", "London", "South East"))
replic = 10
ggsave(paste0("./output/28Dec_proj_", replic, ".pdf"), plll, width = 34, height = 16, units = "cm", useDingbats = FALSE)
ggsave(paste0("./output/28Dec_proj_", replic, ".png"), plll, width = 34, height = 16, units = "cm")

tbSO = summarize_projection(proj_SO, "2020-12-15", popsize)
tbSa = summarize_projection(proj_Sa, "2020-12-15", popsize)
tbSb = summarize_projection(proj_Sb, "2020-12-15", popsize)
tbSC = summarize_projection(proj_SC, "2020-12-15", popsize)
fwrite(tbSO, "./output/28Dec_table_SO.csv");
fwrite(tbSa, "./output/28Dec_table_Sa.csv");
fwrite(tbSb, "./output/28Dec_table_Sb.csv");
fwrite(tbSC, "./output/28Dec_table_SC.csv");
tbS = england_only(list(tbSO, tbSa, tbSb, tbSC), c("All open", "P open, SK open, SN 1/1, U closed", "P 1/1, SK open, SN 1/1, U closed", "All closed"))
fwrite(tbS, "./output/28Dec_table_S.csv");


proj_SOv = project(england_pops, tiers = TRUE, tier2 = tier2, tier3 = tier3, cb_date = lockdown_dates, cb_duration = lockdown_durations, lockdown = lockdownER, se = seER, cb_behaviour = "default",
    school_breaks = school_breaks_2, school_factors_r = c(0.0, 1.0), school_factors_c = c(0.0, 1.0),
    d_school_date0 = "2021-01-04", d_school_date1 = "2021-04-01", d_school_contact = rep(1, 16),
    vacc = make_vaccine_schedule("2021-01-01", rep(1500000, 365), targeting_old_to_young, england_pops), ei_v = rep(0.375, 16), ed_vi = rep(0.6, 16))
proj_Sav = project(england_pops, tiers = TRUE, tier2 = tier2, tier3 = tier3, cb_date = lockdown_dates, cb_duration = lockdown_durations, lockdown = lockdownER, se = seER, cb_behaviour = "default",
    school_breaks = school_breaks_2, school_factors_r = c(0.0, 0.79), school_factors_c = c(0.0, 1.0),
    d_school_date0 = "2021-01-04", d_school_date1 = "2021-04-01", d_school_contact = dsc_Sa,
    vacc = make_vaccine_schedule("2021-01-01", rep(1500000, 365), targeting_old_to_young, england_pops), ei_v = rep(0.375, 16), ed_vi = rep(0.6, 16))
proj_Sbv = project(england_pops, tiers = TRUE, tier2 = tier2, tier3 = tier3, cb_date = lockdown_dates, cb_duration = lockdown_durations, lockdown = lockdownER, se = seER, cb_behaviour = "default",
    school_breaks = school_breaks_2, school_factors_r = c(0.0, 0.54), school_factors_c = c(0.0, 1.0),
    d_school_date0 = "2021-01-04", d_school_date1 = "2021-04-01", d_school_contact = dsc_Sb,
    vacc = make_vaccine_schedule("2021-01-01", rep(1500000, 365), targeting_old_to_young, england_pops), ei_v = rep(0.375, 16), ed_vi = rep(0.6, 16))
proj_SCv = project(england_pops, tiers = TRUE, tier2 = tier2, tier3 = tier3, cb_date = lockdown_dates, cb_duration = lockdown_durations, lockdown = lockdownER, se = seER, cb_behaviour = "default",
    school_breaks = school_breaks_2, school_factors_r = c(0.0, 0.0), school_factors_c = c(0.0, 1.0),
    d_school_date0 = "2021-01-04", d_school_date1 = "2021-04-01", d_school_contact = rep(0, 16),
    vacc = make_vaccine_schedule("2021-01-01", rep(1500000, 365), targeting_old_to_young, england_pops), ei_v = rep(0.375, 16), ed_vi = rep(0.6, 16))

pllv = plot_projection_2parts(list(proj_SOv, proj_Sav, proj_Sbv, proj_SCv), list("All open", "M1: P open, SK open, SN alt, U closed", "M2: P alt, SK open, SN alt, U closed", "All closed"), 
    "2020-10-01", c("#5679a3", "#f47942", "#59a14f", "#31a1b3"), c("East of England", "London", "South East"))
replic = 10
ggsave(paste0("./output/28Dec_proj_vacc_", replic, ".pdf"), pllv, width = 34, height = 16, units = "cm", useDingbats = FALSE)
ggsave(paste0("./output/28Dec_proj_vacc_", replic, ".png"), pllv, width = 34, height = 16, units = "cm")

tbSOv = summarize_projection(proj_SOv, "2020-12-15", popsize)
tbSav = summarize_projection(proj_Sav, "2020-12-15", popsize)
tbSbv = summarize_projection(proj_Sbv, "2020-12-15", popsize)
tbSCv = summarize_projection(proj_SCv, "2020-12-15", popsize)
fwrite(tbSOv, "./output/28Dec_table_SOv.csv");
fwrite(tbSav, "./output/28Dec_table_Sav.csv");
fwrite(tbSbv, "./output/28Dec_table_Sbv.csv");
fwrite(tbSCv, "./output/28Dec_table_SCv.csv");
tbSv = england_only(list(tbSOv, tbSav, tbSbv, tbSCv), c("All open", "P open, SK open, SN 1/1, U closed", "P 1/1, SK open, SN 1/1, U closed", "All closed"))
fwrite(tbSv, "./output/28Dec_table_Sv.csv");

### South Eastern England only
southeast = c("London", "East of England", "South East")
tbSO = summarize_projection(proj_SO[population %in% southeast], "2020-12-15", popsize)
tbSa = summarize_projection(proj_Sa[population %in% southeast], "2020-12-15", popsize)
tbSb = summarize_projection(proj_Sb[population %in% southeast], "2020-12-15", popsize)
tbSC = summarize_projection(proj_SC[population %in% southeast], "2020-12-15", popsize)
tbS = england_only(list(tbSO, tbSa, tbSb, tbSC), c("All open", "P open, SK open, SN 1/1, U closed", "P 1/1, SK open, SN 1/1, U closed", "All closed"))
fwrite(tbS, "./output/28Dec_table_S_SE.csv");


tbSOv = summarize_projection(proj_SOv[population %in% southeast], "2020-12-15", popsize)
tbSav = summarize_projection(proj_Sav[population %in% southeast], "2020-12-15", popsize)
tbSbv = summarize_projection(proj_Sbv[population %in% southeast], "2020-12-15", popsize)
tbSCv = summarize_projection(proj_SCv[population %in% southeast], "2020-12-15", popsize)
tbSv = england_only(list(tbSOv, tbSav, tbSbv, tbSCv), c("All open", "P open, SK open, SN 1/1, U closed", "P 1/1, SK open, SN 1/1, U closed", "All closed"))
fwrite(tbSv, "./output/28Dec_table_Sv_SE.csv");



# Note for vaccination without new strain, or scenarios without new strain, can set all v2_when in posteriorsI to 9999

trace = function(proj)
{
    dat = proj[,
              .(dth1 = sum(death_o),
                dth2 = sum(death2_o), 
                adm1 = sum(hosp_undetected_o),
                adm2 = sum(hosp_undetected2_o),
                bed1 = sum(hosp_p - hosp_undetected_p),
                bed2 = sum(hosp2_p - hosp_undetected2_p),
                icu1 = sum(icu_p),
                icu2 = sum(icu2_p)), by = .(population, run, t)]
    dat = dat[, lapply(.SD, mean), .SDcols = 4:11, by = .(population, t)]
    
    dat = melt(dat, id.vars = 1:2)
    dat[, var := str_sub(variable, 1, 3)]
    ggplot(dat) + geom_line(aes(x = t, y = value, colour = variable)) + facet_grid(population~var)
}

trace(proj_00)

proj_Ec_Vo = project(england_pops, tiers = TRUE, tier2 = tier2, tier3 = tier3, cb_date = lockdown_dates, cb_duration = lockdown_durations, lockdown = lockdownER, se = seER, school_breaks = school_breaks_c, cb_behaviour = "default",
    vacc = make_vaccine_schedule("2021-01-01", rep(2000000, 365), targeting_old_to_young, england_pops), ei_v = rep(0.6, 16), ed_vi = rep(0.875, 16))
proj_Ec_Vy = project(england_pops, tiers = TRUE, tier2 = tier2, tier3 = tier3, cb_date = lockdown_dates, cb_duration = lockdown_durations, lockdown = lockdownER, se = seER, school_breaks = school_breaks_c, cb_behaviour = "default",
    vacc = make_vaccine_schedule("2021-01-01", rep(2000000, 365), targeting_young_to_old, england_pops), ei_v = rep(0.6, 16), ed_vi = rep(0.875, 16))


# Test
proj_basic = project(3, tiers = FALSE, tier2 = tier2, tier3 = tier3, cb_date = NA,           cb_duration = NA, lockdown = NA,                  close_schools = FALSE, cb_behaviour = "default")

proj_ldE4c = project(which_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = "2020-12-20", cb_duration = 14, lockdown = lockdownE, se = seE, close_schools = TRUE,  cb_behaviour = "default")


# TO DO.
# account for school closures reducing R
# 





# OLD CODE


# RUN PROJECTIONS

# ldN4o = lockdown, like Northern Ireland, 4 weeks, schools open
proj_basic = project(england_pops, tiers = FALSE, tier2 = tier2, tier3 = tier3, cb_date = NA,           cb_duration = NA, lockdown = NA,                  close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13")
proj_tiers = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = NA,           cb_duration = NA, lockdown = NA,        se = seN, close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13")
proj_ldN4o = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = "2020-11-05", cb_duration = 27, lockdown = lockdownN, se = seN, close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13")
proj_ldW4o = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = "2020-11-05", cb_duration = 27, lockdown = lockdownW, se = seW, close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13")
proj_ldN4c = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = "2020-11-05", cb_duration = 27, lockdown = lockdownN, se = seN, close_schools = TRUE,  cb_behaviour = "default", expire = "2020-10-13")
proj_ldW4c = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = "2020-11-05", cb_duration = 27, lockdown = lockdownW, se = seW, close_schools = TRUE,  cb_behaviour = "default", expire = "2020-10-13")

# total infected
proj_basic[t == 455, nicepc(1 - sum(S)/sum(S+E+Ip+Is+Ia+R)), by = .(population)]

plot_projection(list(proj_basic, proj_tiers), list("No tiers", "Tiers"), "2020-10-01", "Dark2")
ggsave("./figures/proj_tiers.png", width = 24, height = 12, units = "cm")
ggsave("./figures/proj_tiers.pdf", width = 24, height = 12, units = "cm", useDingbats = FALSE)
plot_projection(list(proj_tiers, proj_ldN4o, proj_ldN4c), list("Tiers only", "NI lockdown, schools open", "NI lockdown, schools closed"), "2020-10-01", "Set2")
ggsave("./figures/proj_ld_N.png", width = 24, height = 12, units = "cm")
ggsave("./figures/proj_ld_N.pdf", width = 24, height = 12, units = "cm", useDingbats = FALSE)
plot_projection(list(proj_tiers, proj_ldW4o, proj_ldW4c), list("Tiers only", "Wales lockdown, schools open", "Wales lockdown, schools closed"), "2020-10-01", "Set2")
ggsave("./figures/proj_ld_W.png", width = 24, height = 12, units = "cm")
ggsave("./figures/proj_ld_W.pdf", width = 24, height = 12, units = "cm", useDingbats = FALSE)

# Different duration of lockdown
proj_ldW1o = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = "2020-11-05", cb_duration =  6, lockdown = lockdownW, se = seW, close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13")
proj_ldW2o = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = "2020-11-05", cb_duration = 13, lockdown = lockdownW, se = seW, close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13")
proj_ldW3o = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = "2020-11-05", cb_duration = 20, lockdown = lockdownW, se = seW, close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13")
proj_ldW4o = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = "2020-11-05", cb_duration = 27, lockdown = lockdownW, se = seW, close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13")
proj_ldW5o = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = "2020-11-05", cb_duration = 34, lockdown = lockdownW, se = seW, close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13")
proj_ldW6o = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = "2020-11-05", cb_duration = 41, lockdown = lockdownW, se = seW, close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13")

# Different timing of lockdown
proj_ldW4o_10_08 = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = "2020-10-08", cb_duration = 27, lockdown = lockdownW, se = seW, close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13")
proj_ldW4o_10_15 = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = "2020-10-15", cb_duration = 27, lockdown = lockdownW, se = seW, close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13")
proj_ldW4o_10_22 = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = "2020-10-22", cb_duration = 27, lockdown = lockdownW, se = seW, close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13")
proj_ldW4o_10_29 = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = "2020-10-29", cb_duration = 27, lockdown = lockdownW, se = seW, close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13")
proj_ldW4o_11_05 = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = "2020-11-05", cb_duration = 27, lockdown = lockdownW, se = seW, close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13")
proj_ldW4o_11_12 = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = "2020-11-12", cb_duration = 27, lockdown = lockdownW, se = seW, close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13")
proj_ldW4o_11_19 = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = "2020-11-19", cb_duration = 27, lockdown = lockdownW, se = seW, close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13")

# Sensitivity analyses with waning and seasonality
proj_tiers_w  = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = NA,           cb_duration = NA, lockdown = NA,        se = seW, close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13",
    waning_duration = 40*7, seasonality = 0)
proj_ldW4o_w  = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = "2020-11-05", cb_duration = 27, lockdown = lockdownW, se = seW, close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13",
    waning_duration = 40*7, seasonality = 0)
proj_tiers_s  = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = NA,           cb_duration = NA, lockdown = NA,        se = seW, close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13",
    waning_duration = -1, seasonality = 0.1)
proj_ldW4o_s  = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = "2020-11-05", cb_duration = 27, lockdown = lockdownW, se = seW, close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13",
    waning_duration = -1, seasonality = 0.1)
proj_tiers_ws = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = NA,           cb_duration = NA, lockdown = NA,        se = seW, close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13",
    waning_duration = 40*7, seasonality = 0.1)
proj_ldW4o_ws = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = "2020-11-05", cb_duration = 27, lockdown = lockdownW, se = seW, close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13",
    waning_duration = 40*7, seasonality = 0.1)

plot_projection(list(proj_tiers_w, proj_ldW4o_w), list("Tiers + waning immunity", "Lockdown + waning immunity"), "2020-10-01", "Set1")
ggsave("./figures/proj_sensitivity_w.png", width = 24, height = 12, units = "cm")
ggsave("./figures/proj_sensitivity_w.pdf", width = 24, height = 12, units = "cm", useDingbats = FALSE)
plot_projection(list(proj_tiers_s, proj_ldW4o_s), list("Tiers + seasonality", "Lockdown + seasonality"), "2020-10-01", "Set1")
ggsave("./figures/proj_sensitivity_s.png", width = 24, height = 12, units = "cm")
ggsave("./figures/proj_sensitivity_s.pdf", width = 24, height = 12, units = "cm", useDingbats = FALSE)
plot_projection(list(proj_tiers_ws, proj_ldW4o_ws), list("Tiers + waning immunity + seasonality", "Lockdown + waning immunity + seasonality"), "2020-10-01", "Set1")
ggsave("./figures/proj_sensitivity_ws.png", width = 24, height = 12, units = "cm")
ggsave("./figures/proj_sensitivity_ws.pdf", width = 24, height = 12, units = "cm", useDingbats = FALSE)


# Summary tables
tbBa_before = summarize_projection(proj_basic, "2020-10-01", popsize, wh = "before")

tbBa = summarize_projection(proj_basic, "2020-10-01", popsize)
fwrite(tbBa, "./figures/proj_basic.csv");
tbTi = summarize_projection(proj_tiers, "2020-10-01", popsize)
fwrite(tbTi, "./figures/proj_tiers.csv");
tbNo = summarize_projection(proj_ldN4o, "2020-10-01", popsize)
fwrite(tbNo, "./figures/proj_ldN4o.csv");
tbNc = summarize_projection(proj_ldN4c, "2020-10-01", popsize)
fwrite(tbNc, "./figures/proj_ldN4c.csv");
tbWo = summarize_projection(proj_ldW4o, "2020-10-01", popsize)
fwrite(tbWo, "./figures/proj_ldW4o.csv");
tbWc = summarize_projection(proj_ldW4c, "2020-10-01", popsize)
fwrite(tbWc, "./figures/proj_ldW4c.csv");

tbEo = summarize_projection(proj_ldE4o, "2020-10-01", popsize)
fwrite(tbEo, "./figures/proj_ldE4o.csv");
tbEc = summarize_projection(proj_ldE4c, "2020-10-01", popsize)
fwrite(tbEc, "./figures/proj_ldE4c.csv");


tbEngland = england_only(list(tbBa, tbTi, tbNo, tbNc, tbWo, tbWc), 
    c("Baseline", "Tiers only", "NI-style lockdown, schools open", "NI-style lockdown, schools closed", "Wales-style lockdown, schools open", "Wales-style lockdown, schools closed"))
fwrite(tbEngland, "./figures/table_england.csv");

tbs_tiers_w = summarize_projection(proj_tiers_w, "2020-10-01", popsize)
tbs_ldown_w = summarize_projection(proj_ldW4o_w, "2020-10-01", popsize)
tbs_tiers_s = summarize_projection(proj_tiers_s, "2020-10-01", popsize)
tbs_ldown_s = summarize_projection(proj_ldW4o_s, "2020-10-01", popsize)
tbs_tiers_ws = summarize_projection(proj_tiers_ws, "2020-10-01", popsize)
tbs_ldown_ws = summarize_projection(proj_ldW4o_ws, "2020-10-01", popsize)
tbSensitivity = england_only(list(tbTi, tbWo, tbs_tiers_s, tbs_ldown_s, tbs_tiers_w, tbs_ldown_w, tbs_tiers_ws, tbs_ldown_ws), 
    c("Tiers only", "Lockdown", "Tiers only + seasonality", "Lockdown + seasonality", "Tiers only + waning", "Lockdown + waning", "Tiers only + seasonality + waning", "Lockdown + seasonality + waning"))
fwrite(tbSensitivity, "./figures/table_sensitivity.csv");

R_tiers = arrs_tier(proj_tiers)
R_ldN4o = arrs_ld(proj_ldN4o)
R_ldN4c = arrs_ld(proj_ldN4c)
R_ldW4o = arrs_ld(proj_ldW4o)
R_ldW4c = arrs_ld(proj_ldW4c)

fwrite(R_tiers, "./figures/R0_table_tiers.csv")
fwrite(R_ldN4o, "./figures/R0_table_N4o.csv")
fwrite(R_ldN4c, "./figures/R0_table_N4c.csv")
fwrite(R_ldW4o, "./figures/R0_table_W4o.csv")
fwrite(R_ldW4c, "./figures/R0_table_W4c.csv")

series_weeks = rbind(
    summarize_projection(proj_tiers, "2020-10-01", popsize, "No lockdown"),
    summarize_projection(proj_ldW1o, "2020-10-01", popsize, "1 week"),
    summarize_projection(proj_ldW2o, "2020-10-01", popsize, "2 weeks"),
    summarize_projection(proj_ldW3o, "2020-10-01", popsize, "3 weeks"),
    summarize_projection(proj_ldW4o, "2020-10-01", popsize, "4 weeks"),
    summarize_projection(proj_ldW5o, "2020-10-01", popsize, "5 weeks"),
    summarize_projection(proj_ldW6o, "2020-10-01", popsize, "6 weeks")
)
series_weeks[, scenario := factor(scenario, unique(scenario))]

series_timing = rbind(
    summarize_projection(proj_ldW4o_10_08, "2020-10-01", popsize, "8 Oct"),
    summarize_projection(proj_ldW4o_10_15, "2020-10-01", popsize, "15 Oct"),
    summarize_projection(proj_ldW4o_10_22, "2020-10-01", popsize, "22 Oct"),
    summarize_projection(proj_ldW4o_10_29, "2020-10-01", popsize, "29 Oct"),
    summarize_projection(proj_ldW4o_11_05, "2020-10-01", popsize, "5 Nov"),
    summarize_projection(proj_ldW4o_11_12, "2020-10-01", popsize, "12 Nov"),
    summarize_projection(proj_ldW4o_11_19, "2020-10-01", popsize, "19 Nov")
)
series_timing[, scenario := factor(scenario, unique(scenario))]

series_what = rbind(
    summarize_projection(proj_basic, "2020-10-01", popsize, "Nothing"),
    summarize_projection(proj_tiers, "2020-10-01", popsize, "Tiers"),
    summarize_projection(proj_ldN4o, "2020-10-01", popsize, "Ld N/o"),
    summarize_projection(proj_ldN4c, "2020-10-01", popsize, "Ld N/c"),
    summarize_projection(proj_ldW4o, "2020-10-01", popsize, "Ld W/o"),
    summarize_projection(proj_ldW4c, "2020-10-01", popsize, "Ld W/c")
)
series_what[, scenario := factor(scenario, unique(scenario))]


# Regional plots: Type of intervention
theme_set(theme_cowplot(font_size = 6))
pr1  = plot_indicator_regions(series_what, "Deaths", legpos = "top", bigtitle = "Type of intervention")
pr2  = plot_indicator_regions(series_what, "Admissions")
pr3  = plot_indicator_regions(series_what, "Cases")
pr4  = plot_indicator_regions(series_what, "Infections")
pr5  = plot_indicator_regions(series_what, "Peak ICU beds")
pr6  = plot_indicator_regions(series_what, "Weeks of high ICU occupancy", "High ICU weeks", suffix = "", unit = 1)
pr7  = plot_indicator_regions(series_what, "Peak hospital beds", "Peak hosp beds")
pr8  = plot_indicator_regions(series_what, "Weeks of high hospital occupancy", "High hosp weeks", suffix = "", unit = 1)
pr9  = plot_indicator_regions(series_what, "Weeks in Tier 2", "Weeks in\nTier 2", suffix = "", unit = 1)
pr10 = plot_indicator_regions(series_what, "Weeks in Tier 3", "Weeks in\nTier 3", suffix = "", unit = 1)
pr11 = plot_indicator_regions(series_what, "Weeks in lockdown", "Weeks in\nlockdown", suffix = "", unit = 1)
ppr1 = plot_grid(pr1, pr2, pr3, pr4, pr5, pr6, pr7, pr8, pr9, pr10, pr11, ncol = 1, align = "v", rel_heights = c(1.5, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1))
ggsave("./figures/regional_what.pdf", ppr1, width = 7, height = 21, units = "cm")

# Regional plots: Duration of lockdown
pr1  = plot_indicator_regions(series_weeks, "Deaths", legpos = "top", bigtitle = "Duration of lockdown")
pr2  = plot_indicator_regions(series_weeks, "Admissions")
pr3  = plot_indicator_regions(series_weeks, "Cases")
pr4  = plot_indicator_regions(series_weeks, "Infections")
pr5  = plot_indicator_regions(series_weeks, "Peak ICU beds")
pr6  = plot_indicator_regions(series_weeks, "Weeks of high ICU occupancy", "High ICU weeks", suffix = "", unit = 1)
pr7  = plot_indicator_regions(series_weeks, "Peak hospital beds", "Peak hosp beds")
pr8  = plot_indicator_regions(series_weeks, "Weeks of high hospital occupancy", "High hosp weeks", suffix = "", unit = 1)
pr9  = plot_indicator_regions(series_weeks, "Weeks in Tier 2", "Weeks in\nTier 2", suffix = "", unit = 1)
pr10 = plot_indicator_regions(series_weeks, "Weeks in Tier 3", "Weeks in\nTier 3", suffix = "", unit = 1)
pr11 = plot_indicator_regions(series_weeks, "Weeks in lockdown", "Weeks in\nlockdown", suffix = "", unit = 1)
ppr2 = plot_grid(pr1, pr2, pr3, pr4, pr5, pr6, pr7, pr8, pr9, pr10, pr11, ncol = 1, align = "v", rel_heights = c(1.5, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1))
ggsave("./figures/regional_weeks.pdf", ppr2, width = 7, height = 21, units = "cm")

# Regional plots: Timing of lockdown
pr1  = plot_indicator_regions(series_timing, "Deaths", legpos = "top", bigtitle = "Timing of lockdown")
pr2  = plot_indicator_regions(series_timing, "Admissions")
pr3  = plot_indicator_regions(series_timing, "Cases")
pr4  = plot_indicator_regions(series_timing, "Infections")
pr5  = plot_indicator_regions(series_timing, "Peak ICU beds")
pr6  = plot_indicator_regions(series_timing, "Weeks of high ICU occupancy", "High ICU weeks", suffix = "", unit = 1)
pr7  = plot_indicator_regions(series_timing, "Peak hospital beds", "Peak hosp beds")
pr8  = plot_indicator_regions(series_timing, "Weeks of high hospital occupancy", "High hosp weeks", suffix = "", unit = 1)
pr9  = plot_indicator_regions(series_timing, "Weeks in Tier 2", "Weeks in\nTier 2", suffix = "", unit = 1)
pr10 = plot_indicator_regions(series_timing, "Weeks in Tier 3", "Weeks in\nTier 3", suffix = "", unit = 1)
pr11 = plot_indicator_regions(series_timing, "Weeks in lockdown", "Weeks in\nlockdown", suffix = "", unit = 1)
ppr3 = plot_grid(pr1, pr2, pr3, pr4, pr5, pr6, pr7, pr8, pr9, pr10, pr11, ncol = 1, align = "v", rel_heights = c(1.5, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1))
ggsave("./figures/regional_timing.pdf", ppr3, width = 7, height = 21, units = "cm")

ppr = plot_grid(ppr1, ppr2, ppr3, nrow = 1, label_size = 8, labels = letters)
ggsave("./figures/regional.pdf", ppr, width = 21, height = 21, units = "cm", useDingbats = FALSE)
ggsave("./figures/regional.png", ppr, width = 21, height = 21, units = "cm")



# Comparative plots: Type of intervention
p0 = plot_cum_deaths(list(proj_basic, proj_tiers, proj_ldN4o, proj_ldN4c, proj_ldW4o, proj_ldW4c), 
    c("Baseline", "Tiers", "Ld N/o", "Ld N/c", "Ld W/o", "Ld W/c"), 
    "2020-10-01", ystart = 52, ydiff = -3.7, titl = "Type of intervention")
p7 = plot_indicators_england2(series_what, c("Weeks of high ICU occupancy",
    "Weeks of high hospital occupancy"), "Hospital pressure\n(weeks)", legpos = "bottom")
p8 = plot_indicators_england_stack2(series_what, c(
    "Weeks in Tier 2",
    "Weeks in Tier 3",
    "Weeks in lockdown"), "Time under\nrestrictions (weeks)", stack_order = c("Weeks in Tier 2", "Weeks in Tier 3", "Weeks in lockdown"), legpos = "bottom")

theme_set(theme_cowplot(font_size = 6))
pp1 = plot_grid(p0, p7, p8, ncol = 1, rel_heights = c(2.8, 2, 2), align = "v", axis = "bottom")
ggsave("./figures/big_whatA.pdf", pp1, width = 7, height = 12, units = "cm")


# Comparative plots: Duration of lockdown
p0 = plot_cum_deaths(list(proj_tiers, proj_ldW1o, proj_ldW2o, proj_ldW3o, proj_ldW4o, proj_ldW5o, proj_ldW6o), 
    c("No lockdown", "1 week", "2 weeks", "3 weeks", "4 weeks", "5 weeks", "6 weeks"), 
    "2020-10-01", ystart = 48.3, ydiff = -3.2, titl = "Duration of lockdown")
p7 = plot_indicators_england2(series_weeks, c("Weeks of high ICU occupancy",
    "Weeks of high hospital occupancy"), "Hospital pressure\n(weeks)", legpos = "bottom")
p8 = plot_indicators_england_stack2(series_weeks, c(
    "Weeks in Tier 2",
    "Weeks in Tier 3",
    "Weeks in lockdown"), "Time under\nrestrictions (weeks)", stack_order = c("Weeks in Tier 2", "Weeks in Tier 3", "Weeks in lockdown"), legpos = "bottom")

theme_set(theme_cowplot(font_size = 6))
pp2 = plot_grid(p0, p7, p8, ncol = 1, rel_heights = c(2.8, 2, 2), align = "v", axis = "bottom")
ggsave("./figures/big_durationA.pdf", pp2, width = 7, height = 12, units = "cm")


# Comparative plots: Timing of lockdown
p0 = plot_cum_deaths(list(proj_ldW4o_10_08, proj_ldW4o_10_15, proj_ldW4o_10_22, proj_ldW4o_10_29, proj_ldW4o_11_05, proj_ldW4o_11_12, proj_ldW4o_11_19), 
    c(  "                  8 Oct", 
        "                      15 Oct", 
        "                           22 Oct", 
        "29 Oct", 
        "5 Nov", 
        "12 Nov", 
        "19 Nov"), 
    "2020-10-01", ystart = 45.8, ydiff = -2.85, titl = "Timing of lockdown")
p7 = plot_indicators_england2(series_timing, c("Weeks of high ICU occupancy",
    "Weeks of high hospital occupancy"), "Hospital pressure\n(weeks)", legpos = "bottom")
p8 = plot_indicators_england_stack2(series_timing, c(
    "Weeks in Tier 2",
    "Weeks in Tier 3",
    "Weeks in lockdown"), "Time under\nrestrictions (weeks)", stack_order = c("Weeks in Tier 2", "Weeks in Tier 3", "Weeks in lockdown"), legpos = "bottom")

theme_set(theme_cowplot(font_size = 6))
pp3 = plot_grid(p0, p7, p8, ncol = 1, rel_heights = c(2.8, 2, 2), align = "v", axis = "bottom")
ggsave("./figures/big_timingA.pdf", pp3, width = 7, height = 12, units = "cm")


pp = plot_grid(pp1, pp2, pp3, nrow = 1, label_size = 8, labels = letters)
ggsave("./figures/bigA.pdf", pp, width = 21, height = 12, units = "cm", useDingbats = FALSE)
ggsave("./figures/bigA.png", pp, width = 21, height = 12, units = "cm")










# New 29 Nov 2020 - England lockdown
# RUN PROJECTIONS


# ldN4o = lockdown, like Northern Ireland, 4 weeks, schools open
#proj_basic = project(england_pops, tiers = FALSE, tier2 = tier2, tier3 = tier3, cb_date = NA,           cb_duration = NA, lockdown = NA,                  close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13")
#proj_tiers = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = NA,           cb_duration = NA, lockdown = NA,        se = seN, close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13")
#proj_ldN4o = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = "2020-11-05", cb_duration = 27, lockdown = lockdownN, se = seN, close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13")
#proj_ldW4o = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = "2020-11-05", cb_duration = 27, lockdown = lockdownW, se = seW, close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13")
#proj_ldN4c = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = "2020-11-05", cb_duration = 27, lockdown = lockdownN, se = seN, close_schools = TRUE,  cb_behaviour = "default", expire = "2020-10-13")
#proj_ldW4c = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = "2020-11-05", cb_duration = 27, lockdown = lockdownW, se = seW, close_schools = TRUE,  cb_behaviour = "default", expire = "2020-10-13")
proj_ldE4o = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = "2020-11-05", cb_duration = 27, lockdown = lockdownE, se = seE, close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13")
proj_ldE4c = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = "2020-11-05", cb_duration = 27, lockdown = lockdownE, se = seE, close_schools = TRUE,  cb_behaviour = "default", expire = "2020-10-13")

plot_projection(list(proj_tiers, proj_ldE4o, proj_ldE4c), list("Tiers only", "England lockdown, schools open", "England lockdown, schools closed"), "2020-10-01", "Set2")
ggsave("./figures/proj_ld_E.png", width = 24, height = 12, units = "cm")
ggsave("./figures/proj_ld_E.pdf", width = 24, height = 12, units = "cm", useDingbats = FALSE)

# Different duration of lockdown
proj_ldE1o = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = "2020-11-05", cb_duration =  6, lockdown = lockdownE, se = seE, close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13")
proj_ldE2o = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = "2020-11-05", cb_duration = 13, lockdown = lockdownE, se = seE, close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13")
proj_ldE3o = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = "2020-11-05", cb_duration = 20, lockdown = lockdownE, se = seE, close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13")
proj_ldE4o = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = "2020-11-05", cb_duration = 27, lockdown = lockdownE, se = seE, close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13")
proj_ldE5o = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = "2020-11-05", cb_duration = 34, lockdown = lockdownE, se = seE, close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13")
proj_ldE6o = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = "2020-11-05", cb_duration = 41, lockdown = lockdownE, se = seE, close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13")

# Different timing of lockdown
proj_ldE4o_10_08 = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = "2020-10-08", cb_duration = 27, lockdown = lockdownE, se = seE, close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13")
proj_ldE4o_10_15 = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = "2020-10-15", cb_duration = 27, lockdown = lockdownE, se = seE, close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13")
proj_ldE4o_10_22 = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = "2020-10-22", cb_duration = 27, lockdown = lockdownE, se = seE, close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13")
proj_ldE4o_10_29 = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = "2020-10-29", cb_duration = 27, lockdown = lockdownE, se = seE, close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13")
proj_ldE4o_11_05 = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = "2020-11-05", cb_duration = 27, lockdown = lockdownE, se = seE, close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13")
proj_ldE4o_11_12 = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = "2020-11-12", cb_duration = 27, lockdown = lockdownE, se = seE, close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13")
proj_ldE4o_11_19 = project(england_pops, tiers =  TRUE, tier2 = tier2, tier3 = tier3, cb_date = "2020-11-19", cb_duration = 27, lockdown = lockdownE, se = seE, close_schools = FALSE, cb_behaviour = "default", expire = "2020-10-13")


R_ldE4o = arrs_ld(proj_ldE4o)
R_ldE4c = arrs_ld(proj_ldE4c)

fwrite(R_ldE4o, "./figures/R0_table_E4o.csv")
fwrite(R_ldE4c, "./figures/R0_table_E4c.csv")

series_weeks = rbind(
    summarize_projection(proj_tiers, "2020-10-01", popsize, "No lockdown"),
    summarize_projection(proj_ldE1o, "2020-10-01", popsize, "1 week"),
    summarize_projection(proj_ldE2o, "2020-10-01", popsize, "2 weeks"),
    summarize_projection(proj_ldE3o, "2020-10-01", popsize, "3 weeks"),
    summarize_projection(proj_ldE4o, "2020-10-01", popsize, "4 weeks"),
    summarize_projection(proj_ldE5o, "2020-10-01", popsize, "5 weeks"),
    summarize_projection(proj_ldE6o, "2020-10-01", popsize, "6 weeks")
)
series_weeks[, scenario := factor(scenario, unique(scenario))]

series_timing = rbind(
    summarize_projection(proj_ldE4o_10_08, "2020-10-01", popsize, "8 Oct"),
    summarize_projection(proj_ldE4o_10_15, "2020-10-01", popsize, "15 Oct"),
    summarize_projection(proj_ldE4o_10_22, "2020-10-01", popsize, "22 Oct"),
    summarize_projection(proj_ldE4o_10_29, "2020-10-01", popsize, "29 Oct"),
    summarize_projection(proj_ldE4o_11_05, "2020-10-01", popsize, "5 Nov"),
    summarize_projection(proj_ldE4o_11_12, "2020-10-01", popsize, "12 Nov"),
    summarize_projection(proj_ldE4o_11_19, "2020-10-01", popsize, "19 Nov")
)
series_timing[, scenario := factor(scenario, unique(scenario))]

series_what = rbind(
    summarize_projection(proj_basic, "2020-10-01", popsize, "Nothing"),
    summarize_projection(proj_tiers, "2020-10-01", popsize, "Tiers"),
    summarize_projection(proj_ldN4o, "2020-10-01", popsize, "Ld N/o"),
    summarize_projection(proj_ldN4c, "2020-10-01", popsize, "Ld N/c"),
    summarize_projection(proj_ldW4o, "2020-10-01", popsize, "Ld W/o"),
    summarize_projection(proj_ldW4c, "2020-10-01", popsize, "Ld W/c"),
    summarize_projection(proj_ldE4o, "2020-10-01", popsize, "Ld E/o"),
    summarize_projection(proj_ldE4c, "2020-10-01", popsize, "Ld E/c")
)
series_what[, scenario := factor(scenario, unique(scenario))]


# Comparative plots: Type of intervention
p0 = plot_cum_deaths(list(proj_basic, proj_tiers, proj_ldN4o, proj_ldN4c, proj_ldW4o, proj_ldW4c, proj_ldE4o, proj_ldE4c), 
    c("Nothing", "Tiers", "Ld N/o", "Ld N/c", "Ld W/o", "Ld W/c", "Ld E/o", "Ld E/c"), 
    "2020-10-01", ystart = 56, ydiff = -3.8, titl = "Type of intervention")
p1 = plot_icu(list(proj_basic, proj_tiers, proj_ldN4o, proj_ldN4c, proj_ldW4o, proj_ldW4c, proj_ldE4o, proj_ldE4c), 
    c("Nothing", "Tiers", "Ld N/o", "Ld N/c", "Ld W/o", "Ld W/c", "Ld E/o", "Ld E/c"), 
    "2020-10-01")
p7 = plot_indicators_england(series_what, c("Weeks of high ICU occupancy",
    "Weeks of high hospital occupancy"), "Hospital pressure\n(weeks)", legpos = "bottom")
p8 = plot_indicators_england_stack(series_what, c(
    "Weeks in Tier 2",
    "Weeks in Tier 3",
    "Weeks in lockdown"), "Time under\nrestrictions (weeks)", stack_order = c("Weeks in Tier 2", "Weeks in Tier 3", "Weeks in lockdown"), legpos = "bottom")

theme_set(theme_cowplot(font_size = 6))
pp1 = plot_grid(p0, p7, p8, ncol = 1, rel_heights = c(2.8, 2, 2), align = "v", axis = "bottom")
ggsave("./figures/ENG_big_what.pdf", pp1, width = 7, height = 12, units = "cm")


# Comparative plots: Duration of lockdown
p0 = plot_cum_deaths(list(proj_tiers, proj_ldE1o, proj_ldE2o, proj_ldE3o, proj_ldE4o, proj_ldE5o, proj_ldE6o), 
    c("No lockdown", "1 week", "2 weeks", "3 weeks", "4 weeks", "5 weeks", "6 weeks"), 
    "2020-10-01", ystart = 48.5, ydiff = -3.2, titl = "Duration of lockdown")
p1 = plot_icu(list(proj_tiers, proj_ldE1o, proj_ldE2o, proj_ldE3o, proj_ldE4o, proj_ldE5o, proj_ldE6o), 
    c("No lockdown", "1 week", "2 weeks", "3 weeks", "4 weeks", "5 weeks", "6 weeks"), 
    "2020-10-01")
p7 = plot_indicators_england(series_weeks, c("Weeks of high ICU occupancy",
    "Weeks of high hospital occupancy"), "Hospital pressure\n(weeks)", legpos = "bottom")
p8 = plot_indicators_england_stack(series_weeks, c(
    "Weeks in Tier 2",
    "Weeks in Tier 3",
    "Weeks in lockdown"), "Time under\nrestrictions (weeks)", stack_order = c("Weeks in Tier 2", "Weeks in Tier 3", "Weeks in lockdown"), legpos = "bottom")

theme_set(theme_cowplot(font_size = 6))
pp2 = plot_grid(p0, p7, p8, ncol = 1, rel_heights = c(2.8, 2, 2), align = "v", axis = "bottom")
ggsave("./figures/ENG_big_duration.pdf", pp2, width = 7, height = 12, units = "cm")


# Comparative plots: Timing of lockdown
p0 = plot_cum_deaths(list(proj_ldE4o_10_08, proj_ldE4o_10_15, proj_ldE4o_10_22, proj_ldE4o_10_29, proj_ldE4o_11_05, proj_ldE4o_11_12, proj_ldE4o_11_19), 
    c(  "                  8 Oct", 
        "                      15 Oct", 
        "                          22 Oct", 
        "29 Oct", 
        "5 Nov", 
        "12 Nov", 
        "19 Nov"), 
    "2020-10-01", ystart = 48.5, ydiff = -2.9, titl = "Timing of lockdown")
p1 = plot_icu(list(proj_ldE4o_10_08, proj_ldE4o_10_15, proj_ldE4o_10_22, proj_ldE4o_10_29, proj_ldE4o_11_05, proj_ldE4o_11_12, proj_ldE4o_11_19), 
    c(  "                  8 Oct", 
        "                      15 Oct", 
        "                          22 Oct", 
        "29 Oct", 
        "5 Nov", 
        "12 Nov", 
        "19 Nov"), 
    "2020-10-01")
p7 = plot_indicators_england(series_timing, c("Weeks of high ICU occupancy",
    "Weeks of high hospital occupancy"), "Hospital pressure\n(weeks)", legpos = "bottom")
p8 = plot_indicators_england_stack(series_timing, c(
    "Weeks in Tier 2",
    "Weeks in Tier 3",
    "Weeks in lockdown"), "Time under\nrestrictions (weeks)", stack_order = c("Weeks in Tier 2", "Weeks in Tier 3", "Weeks in lockdown"), legpos = "bottom")

theme_set(theme_cowplot(font_size = 6))
pp3 = plot_grid(p0, p7, p8, ncol = 1, rel_heights = c(2.8, 2, 2), align = "v", axis = "bottom")
ggsave("./figures/ENG_big_timing.pdf", pp3, width = 7, height = 12, units = "cm")


pp = plot_grid(pp1, pp2, pp3, nrow = 1, label_size = 8, labels = letters)
ggsave("./figures/ENG_big.pdf", pp, width = 21, height = 12, units = "cm", useDingbats = FALSE)
ggsave("./figures/ENG_big.png", pp, width = 21, height = 12, units = "cm")
