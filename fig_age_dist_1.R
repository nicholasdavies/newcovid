
age_dist = dcast(testC[population %in% c("East of England", "South East", "London") & t %between% c(305, 340),
    .(s1 = sum(Ip + Is + Ia), s2 = sum(Ip2 + Is2 + Ia2)), by = .(population, group, run)][,
    .(age = group * 5 - 2.5, s1 = s1 / sum(s1), s2 = s2 / sum(s2)), by = .(population, run)][,
        .(variant = rep(c("Preexisting", "VOC 202012/01"), each = 3), 
          value = c(quantile(s1, c(0.025, 0.5, 0.975)), quantile(s2, c(0.025, 0.5, 0.975))),
          q = rep(c("lo", "mid", "hi"), 2)), by = .(age)], age + variant ~ q)

age_dist[, age := (age %/% 10) * 10 + 5]
age_dist = age_dist[, lapply(.SD, sum), by = .(age, variant)]

plot_age = ggplot(age_dist) +
    geom_col(aes(age, mid, fill = variant, group = variant), position = "dodge") +
    geom_linerange(aes(age, ymin = lo, ymax = hi, group = variant), position = position_dodge(width = 10)) +
    labs(x = "Age", y = "Proportion of infections", fill = "Variant") +
    theme_cowplot() + theme(legend.position = c(0.05, 0.9)) + ylim(0, 0.35)

plot_age

qsave(plot_age, "./output/ageplot-ch_u.qs")
