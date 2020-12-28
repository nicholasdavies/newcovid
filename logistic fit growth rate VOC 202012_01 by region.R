library(lme4)
library(emmeans)
library(effects)
library(ggplot2)
library(ggthemes)
library(MASS)
library(nlme)
library(mgcv)

# 1. LOAD & PREPROCESS DATA ####
data = read.csv(".\\data\\cog_metadata_microreact_public-2020-12-22-annotated.csv")
data$sample_date = as.Date(data$sample_date)
data$sample_date_num = as.numeric(data$sample_date)
data$nhs_name = as.factor(data$nhs_name)
# var2 == TRUE for samples that are VOC 202012/01
datasubs = data[data$sample_date >= "2020-09-01"&data$nhs_name != "Northern Ireland",]
# we leave out Northern Ireland since VOC 202012/01 has not been reported there yet
datasubs$nhs_name = droplevels(datasubs$nhs_name)

# aggregate binomial data in 2 separate columns
datasubs_aggr = aggregate(datasubs[,c("var2")], by=list(sample_date = datasubs$sample_date, nhs_name = datasubs$nhs_name), FUN="sum")
colnames(datasubs_aggr)[3] = "VOC_202012_01"
datasubs_aggr2 = aggregate(datasubs[,c("var2")], by=list(sample_date = datasubs$sample_date, nhs_name = datasubs$nhs_name), FUN=function(x) sum(!x))
datasubs_aggr$OTHER_STRAINS = datasubs_aggr2$x
datasubs_aggr$sample_date = as.Date(datasubs_aggr$sample_date)
datasubs_aggr$sample_date_num = as.numeric(datasubs_aggr$sample_date)
datasubs_aggr$nhs_name = as.factor(datasubs_aggr$nhs_name)
datasubs_aggr$propVOC = datasubs_aggr$VOC_202012_01/(datasubs_aggr$VOC_202012_01+datasubs_aggr$OTHER_STRAINS)
datasubs_aggr$obs = factor(1:nrow(datasubs_aggr))
head(datasubs_aggr)


# 2. LOGISTIC FIT ####

# logistic fit to model growth in proportion of strains that are VOC 202012/01 
# with or without overdispersion
fit = glm(cbind(VOC_202012_01, OTHER_STRAINS) ~  nhs_name*sample_date_num, family=binomial(logit), data=datasubs_aggr) 
summary(fit)
fit.overdisp = glm(cbind(VOC_202012_01, OTHER_STRAINS) ~  nhs_name*sample_date_num, family=quasibinomial(logit), data=datasubs_aggr) 
summary(fit.overdisp)$dispersion # dispersion coefficient = 1.056009
# likelihood-ratio test show that a quasibinomial GLM with overdispersion does not fit significantly better 
# than a regular binomial GLM without overdispersion
pchisq(summary(fit.overdisp)$dispersion * fit$df.residual, fit$df.residual, lower = F) # p=0.12

# # PS1 in a mixed model we would take into account overdispersion using an observation level random effect,
# # but since there is no overdispersion not much point of going down this route
# fit = glmer(cbind(VOC_202012_01, OTHER_STRAINS) ~  (1|obs) + nhs_name*scale(sample_date_num), family=binomial(logit), data=datasubs_aggr) 
# summary(fit) 

# # PS2 a model with overdispersion and temporally autocorrelated residuals could be
# fit_tempautocor = glmmPQL(cbind(VOC_202012_01, OTHER_STRAINS) ~ nhs_name*sample_date_num,
#                           random = ~ 1 | obs, # observation-level random effect to deal with possible overdispersion
#                           family = binomial(logit),
#                           correlation = corCAR1(form = ~ sample_date_num), # lag-1 temporal autocorrelation 
#                           data=datasubs_aggr)
# summary(fit_tempautocor) # temp autocorrelation PHI=0.2
# # but this gives similar results that the simple binomial GLM above

# PS3 an explicitly spatial model could be done using a binomial generalized additive model:
# (with a gamm temporal autocorrelation and overdispersion could also potentially be added)
# fit_gam = gam(var2 ~ ti(latitude,bs="cr")+ti(longitude,bs="cr")+ti(longitude,latitude,bs="cr") +
#                                sample_date_num, family=binomial(logit), data=datasubs)
# summary(fit_gam)

# PS4 effect plot would be produced like this, but below I make a prettier one using emmeans and ggplot2
# plot(allEffects(fit, x.var="sample_date_num", residuals=TRUE), residuals.pch=16, 
#     smooth.residuals=FALSE, residuals.color=alpha("steelblue",0.3))



# 3. PLOT MODEL FIT ####

extrapolate = 60 # nr of days to extrapolate fit into the future
df = as.data.frame(emmeans(fit, ~ sample_date_num*nhs_name, at=list(sample_date_num=seq(min(datasubs$sample_date_num),
                                                                          max(datasubs$sample_date_num)+extrapolate))), type="response")
df$sample_date = as.Date(df$sample_date_num, origin="1970-01-01")
df$nhs_name = factor(df$nhs_name, 
                     levels=unique(df$nhs_name)[order(df$prob[df$sample_date_num==max(df$sample_date_num)],decreasing = TRUE)])
qplot(data=df, x=sample_date, y=prob, geom="blank") +
    facet_wrap(~nhs_name) +
    geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL, fill=nhs_name), alpha=I(0.3))+
    geom_line(aes(y=prob, colour=nhs_name), alpha=I(0.8))+
    labs(tag = "@TWenseleers\ndata COG-UK") +
    theme(plot.tag.position = "bottomright",
          plot.tag = element_text(vjust = 1, hjust = 1, size=8)) +
    ylab("Representation of VOC strain 202012/01 (%)") +
    theme_hc() + xlab("") + ggtitle("GROWTH OF VOC STRAIN 202012/01 BY NHS REGION IN UK") +
    scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01")),
                       labels=c("M","A","M","J","J","A","S","O","N","D","J","F")) +
    scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                        labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
    coord_cartesian(xlim=c(as.Date("2020-09-01"),as.Date("2021-02-01")), ylim=c(0.0001,0.999)) +
    scale_color_discrete("", h=c(0, 280), c=150) +
    scale_fill_discrete("", h=c(0, 280), c=150) +
    geom_point(data=datasubs_aggr, aes(x=sample_date, y=propVOC, colour=nhs_name)) +
    theme(legend.position = "none")
ggsave(file="growth VOC strain 202012_01 UK.png", width=8, height=6)


# 4. CALCULATE GROWTH RATES & SELECTIVE ADVANTAGE OF VOC STRAIN IN DIFFERENT REGIONS ####

# growth rate of strain = slope of logistic regression = selective advantage s / generation time (6.5 days)

growth_rates = as.data.frame(emtrends(fit, ~ nhs_name, var="sample_date_num"))[,-c(3,4)] # sample_date_num.trend = growth rate of VOC 202012/01 strain
colnames(growth_rates)[2] = "growth_rate"
growth_rates
#                   nhs_name growth_rate  asymp.LCL  asymp.UCL
# 1          East of England  0.08623527 0.07967674 0.09279381
# 2                   London  0.10452845 0.09926676 0.10979013
# 3                 Midlands  0.11568813 0.09976412 0.13161214
# 4 North East and Yorkshire  0.13374421 0.11404233 0.15344609
# 5               North West  0.13046064 0.10271897 0.15820231
# 6                 Scotland  0.09387689 0.06327649 0.12447728
# 7               South East  0.08774132 0.08225427 0.09322838
# 8               South West  0.11158065 0.08342929 0.13973202
# 9                    Wales  0.10317128 0.07354453 0.13279802

# corresponding selective advantage of VOC 202012/01 
s = growth_rates
generation_time = 6.5
s[,c(2:4)] = s[,c(2:4)]*generation_time
colnames(s)[2] = "s"
s
#                   nhs_name         s asymp.LCL asymp.UCL
# 1          East of England 0.5605293 0.5178988 0.6031597
# 2                   London 0.6794349 0.6452339 0.7136359
# 3                 Midlands 0.7519728 0.6484668 0.8554789
# 4 North East and Yorkshire 0.8693374 0.7412751 0.9973996
# 5               North West 0.8479942 0.6676733 1.0283150
# 6                 Scotland 0.6101998 0.4112972 0.8091023
# 7               South East 0.5703186 0.5346527 0.6059845
# 8               South West 0.7252742 0.5422904 0.9082581
# 9                    Wales 0.6706133 0.4780394 0.8631872
range(s$s) # selective advantage in different NHS regions = 56%-87%

growth_rate_avg = as.data.frame(emtrends(fit, ~ 1, var="sample_date_num"))[,-c(3,4)][2:4]
colnames(growth_rate_avg)[1] = "growth_rate"
growth_rate_avg
#   growth_rate asymp.LCL asymp.UCL
# 1   0.1074474 0.1003166 0.1145783

s_avg = growth_rate_avg*generation_time
colnames(s_avg)[1] = "s"
s_avg # avg selective advantage s = 69.84% [65.20%,74.48%] 95% CLs
#            s asymp.LCL asymp.UCL
# 1 0.6984083 0.6520577 0.7447588


# 5. TUKEY POSTHOC TESTS TO TEST FOR DIFFERENCES IN GROWTH RATE ACROSS DIFFERENT NHS REGIONS ####

# as.data.frame(emtrends(fit, pairwise ~ nhs_name, var="sample_date_num", adjust="tukey")$contrasts)[,-4]
# contrast     estimate          SE     z.ratio      p.value
# 1                    East of England - London -0.018293174 0.004290035 -4.26410843 0.0006770443
# 2                  East of England - Midlands -0.029452859 0.008786765 -3.35195685 0.0227825437
# 3  East of England - North East and Yorkshire -0.047508938 0.010594499 -4.48430236 0.0002515467
# 4                East of England - North West -0.044225370 0.014544346 -3.04072594 0.0597994588
# 5                  East of England - Scotland -0.007641615 0.015967304 -0.47857887 0.9999295789
# 6                East of England - South East -0.001506052 0.004362912 -0.34519408 0.9999943633
# 7                East of England - South West -0.025345380 0.014747850 -1.71858133 0.7352347839
# 8                     East of England - Wales -0.016936005 0.015481919 -1.09392151 0.9752668580
# 9                           London - Midlands -0.011159685 0.008556683 -1.30420687 0.9304572491
# 10          London - North East and Yorkshire -0.029215764 0.010404470 -2.80800122 0.1128696255
# 11                        London - North West -0.025932196 0.014406512 -1.80003299 0.6823068473
# 12                          London - Scotland  0.010651559 0.015841856  0.67236814 0.9991051614
# 13                        London - South East  0.016787122 0.003878735  4.32798884 0.0005107594
# 14                        London - South West -0.007052206 0.014611936 -0.48263328 0.9999248975
# 15                             London - Wales  0.001357169 0.015352505  0.08840049 0.9999999999
# 16        Midlands - North East and Yorkshire -0.018056079 0.012925008 -1.39698780 0.8993871047
# 17                      Midlands - North West -0.014772511 0.016320246 -0.90516474 0.9927542168
# 18                        Midlands - Scotland  0.021811244 0.017600206  1.23926075 0.9478571984
# 19                      Midlands - South East  0.027946807 0.008593453  3.25210449 0.0314849221
# 20                      Midlands - South West  0.004107479 0.016501864  0.24890998 0.9999995685
# 21                           Midlands - Wales  0.012516854 0.017161068  0.72937499 0.9983889205
# 22      North East and Yorkshire - North West  0.003283568 0.017360489  0.18914028 0.9999999510
# 23        North East and Yorkshire - Scotland  0.039867323 0.018568883  2.14699634 0.4409801160
# 24      North East and Yorkshire - South East  0.046002886 0.010434730  4.40863195 0.0003555908
# 25      North East and Yorkshire - South West  0.022163558 0.017531334  1.26422537 0.9415765239
# 26           North East and Yorkshire - Wales  0.030572933 0.018153193  1.68416281 0.7565123865
# 27                      North West - Scotland  0.036583755 0.021073633  1.73599664 0.7242035780
# 28                    North West - South East  0.042719318 0.014428381  2.96078386 0.0750225897
# 29                    North West - South West  0.018879990 0.020165373  0.93625790 0.9909160639
# 30                         North West - Wales  0.027289365 0.020708284  1.31779947 0.9263715463
# 31                      Scotland - South East  0.006135563 0.015861746  0.38681510 0.9999863272
# 32                      Scotland - South West -0.017703766 0.021214596 -0.83450873 0.9958401014
# 33                           Scotland - Wales -0.009294390 0.021731310 -0.42769581 0.9999702659
# 34                    South East - South West -0.023839328 0.014633498 -1.62909292 0.7889543577
# 35                         South East - Wales -0.015429953 0.015373029 -1.00370286 0.9856630129
# 36                         South West - Wales  0.008409375 0.020851717  0.40329415 0.9999811094

# slight evidence for slower spread in the East of England, for the rest mostly no sign differences
