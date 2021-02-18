# Requires testC and testU from fig_hypotheses.R

age_distC = dcast(testC$dynamics[population %in% c("East of England", "South East", "London") & t %between% c(335, 366),
    .(s1 = sum(Ip + Is + Ia), s2 = sum(Ip2 + Is2 + Ia2)), by = .(population, group, run)][,
    .(age = group * 5 - 2.5, s1 = s1 / sum(s1), s2 = s2 / sum(s2)), by = .(population, run)][,
        .(variant = rep(c("Preexisting", "VOC 202012/01"), each = 3), 
          value = c(quantile(s1, c(0.025, 0.5, 0.975)), quantile(s2, c(0.025, 0.5, 0.975))),
          q = rep(c("lo", "mid", "hi"), 2)), by = .(age)], age + variant ~ q)

age_distU = dcast(testU$dynamics[population %in% c("East of England", "South East", "London") & t %between% c(335, 366),
    .(s1 = sum(Ip + Is + Ia), s2 = sum(Ip2 + Is2 + Ia2)), by = .(population, group, run)][,
    .(age = group * 5 - 2.5, s1 = s1 / sum(s1), s2 = s2 / sum(s2)), by = .(population, run)][,
        .(variant = rep(c("Preexisting", "VOC 202012/01"), each = 3), 
          value = c(quantile(s1, c(0.025, 0.5, 0.975)), quantile(s2, c(0.025, 0.5, 0.975))),
          q = rep(c("lo", "mid", "hi"), 2)), by = .(age)], age + variant ~ q)

age_distC[, age := (age %/% 10) * 10 + 5]
age_distC = age_distC[, lapply(.SD, sum), by = .(age, variant)]

age_distU[, age := (age %/% 10) * 10 + 5]
age_distU = age_distU[, lapply(.SD, sum), by = .(age, variant)]

plot_ageC = ggplot(age_distC) +
    geom_col(aes(age, mid, fill = variant, group = variant), position = "dodge") +
    geom_linerange(aes(age, ymin = lo, ymax = hi, group = variant), position = position_dodge(width = 10)) +
    labs(x = "Age", y = "Proportion of infections", fill = "Variant") +
    theme_cowplot() + theme(legend.position = c(0.05, 0.9)) + ylim(0, 0.35)

plot_ageU = ggplot(age_distU) +
    geom_col(aes(age, mid, fill = variant, group = variant), position = "dodge") +
    geom_linerange(aes(age, ymin = lo, ymax = hi, group = variant), position = position_dodge(width = 10)) +
    labs(x = "Age", y = "Proportion of infections", fill = "Variant") +
    theme_cowplot() + theme(legend.position = c(0.05, 0.9)) + ylim(0, 0.35)

dS = fread("~/Documents/newcovid/data/sgtf_agedist.csv") # from PHE Pillar 2 testing data
dS[, variant := c("S gene present", "S gene target failure")[sgtf + 1]]

plot_ageE = ggplot(dS) +
    geom_col(aes(ag0, mean, fill = variant, group = variant), position = "dodge") +
    geom_linerange(aes(ag0, ymin = lo, ymax = hi, group = variant), position = position_dodge(width = 10)) +
    labs(x = "Age", y = "Proportion of cases", fill = "Variant") +
    theme_cowplot() + theme(legend.position = c(0.05, 0.9)) + ylim(0, 0.35)

p1 = plot_ageU + labs(title = "Longer infectious period")
p2 = plot_ageC + labs(title = "Increased susceptibility among children")
p3 = plot_ageE + labs(title = "Pillar 2 testing data")

cowplot::plot_grid(p1, p2, p3, nrow = 2, labels = LETTERS, label_size = 14)
ggsave("./output/age_dist_model.pdf", width = 35, height = 28, units = "cm", useDingbats = FALSE)
ggsave("./output/age_dist_model.png", width = 35, height = 28, units = "cm")
