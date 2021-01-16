# ANALYSIS OF SECONDARY ATTACK RATES AND SUSCEPTIBILITY TO VOC 202012/01 IFO AGE BASED ON
# PUBLIC HEALTH ENGLAND DATA REPORTED IN REPORT
# "Investigation of novel SARS-CoV-2 variant Variant of Concern 202012/01, Technical briefing 3", Tables 6 & 7 in
# https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/950823/Variant_of_Concern_VOC_202012_01_Technical_Briefing_3_-_England.pdf

# T. Wenseleers, last updated: 15 January 2021

library(export)
library(lme4)
library(MASS)
library(car)
library(afex)
library(effects)
library(emmeans)
library(ggplot2)
library(ggthemes)
# install from https://github.com/tomwenseleers/export
# library(devtools)
# devtools::install_github("tomwenseleers/export")
library(export) 

# READ IN DATA ####

data_SA = read.csv(".\\multinomial_logistic_fits\\data\\PHE_Variant_of_Concern_VOC_202012_01_Technical_Briefing_3_England_SECATTACKRATES_BYAGE.csv")
head(data_SA)
str(data_SA)
data_SA$data_type = factor(data_SA$data_type) # SGTF (S-gene dropout) or sequencing data
data_SA$age_group = factor(data_SA$age_group) # age of person getting infected
data_SA$variant = factor(data_SA$variant)
data_SA$variant = relevel(data_SA$variant, ref="non_VOC") # variant carried by index patient (VOC or non-VOC)
data_SA$variant_num = as.numeric(data_SA$variant=="VOC") # numeric version of VOC factor
data_SA$obs = factor(1:nrow(data_SA)) # for observation-level random effect to take into account overdispersion in mixed effect models
aggregate(total ~ data_type, data=data_SA, FUN=sum) # total contacts
#    data_type  total
# 1 sequencing  17701
# 2       SGTF 456086
aggregate(count_became_cases  ~ data_type, data=data_SA, FUN=sum) # contacts that became cases
#    data_type count_became_cases
# 1 sequencing               2455
# 2       SGTF              64325
# secondary attack rates:
# 2455/17701 # 13.87%
# 64325/456086 # 14.10%
data_SA = data_SA[complete.cases(data_SA),]
data_SA$prop = data_SA$count_became_cases / data_SA$total


# 1. ANALYSIS OF SECONDARY ATTACK RATES AND SUSCEPTIBILITY TO VOC 202012/01 IFO AGE ####

# 1.1 FIT BINOMIAL GLM TO GIVEN AGGREGATED DATA ####
# we take into account data type (SGTF or sequencing data), variant carried by index case (VOC or not)
# and age of people getting infected plus possible 1st order interaction effects between these

set_sum_contrasts()
fitdSA_bGLM = glm(cbind(count_became_cases, total-count_became_cases) ~ 
                    (data_type + age_group + variant)^2, family=binomial(logit), data=data_SA)
fitdSA_bGLM_od = glm(cbind(count_became_cases, total-count_became_cases) ~ 
                  (data_type + age_group + variant)^2, family=quasibinomial(logit), data=data_SA)
summary(fitdSA_bGLM_od)$dispersion 
# 0.63, <1, so definitely no overdispersion here, so we continue with fitdSA_bGLM


# # PS one could also contemplate fitting binomial GLMMs like those below, but models
# # fitdSA2, fitdSA3 & fitdSA4 give almost identical results as the fixed effect modl above
# fitdSA1 = glmer(cbind(count_became_cases, total-count_became_cases) ~ (1|obs) +
#                     data_type + age_group + variant, family=binomial(logit), data=data_SA)
# fitdSA2 = glmer(cbind(count_became_cases, total-count_became_cases) ~ (1|obs) +
#                     data_type + age_group * variant, family=binomial(logit), data=data_SA)
# fitdSA3 = glmer(cbind(count_became_cases, total-count_became_cases) ~ (1|obs) +
#                     data_type + age_group * variant + data_type : variant, family=binomial(logit), data=data_SA)
# fitdSA4 = glmer(cbind(count_became_cases, total-count_became_cases) ~ (1|data_type/obs) +
#                     age_group * variant, family=binomial(logit), data=data_SA)

# we continue with the fixed effect binomial GLM above, as we are specifically interested in testing
# for effects of data_type, age_group, variant plus possible 1st order interactions

# saveRDS(fitdSA_bGLM, file = ".\\multinomial_logistic_fits\\fits\\bGLMfit_secondaryattack_no overdispersion_best model.rds")
# saveRDS(fitdSA_bGLM_od, file = ".\\multinomial_logistic_fits\\fits\\bGLMfit_secondaryattack_with overdispersion.rds")
# or to directly load previously fitted models
fitdSA_bGLM = readRDS(file = ".\\multinomial_logistic_fits\\fits\\bGLMfit_secondaryattack_no overdispersion_best model.rds")
fitdSA_bGLM_od = readRDS(file = ".\\multinomial_logistic_fits\\fits\\bGLMfit_secondaryattack_with overdispersion.rds")


fit = fitdSA_bGLM 

summary(fit)
# Coefficients:
#                    Estimate Std. Error z value Pr(>|z|)    
# (Intercept)           -1.7287078  0.0173402 -99.694  < 2e-16 ***
# data_type1            -0.0001277  0.0173293  -0.007 0.994121    
# age_group1            -0.7935586  0.0364098 -21.795  < 2e-16 ***
# age_group2            -0.3969570  0.0298829 -13.284  < 2e-16 ***
# age_group3             0.0237290  0.0310826   0.763 0.445215    
# age_group4             0.1075901  0.0300246   3.583 0.000339 ***
# age_group5             0.1456478  0.0288678   5.045 4.53e-07 ***
# age_group6             0.1582106  0.0309902   5.105 3.30e-07 ***
# age_group7             0.2021944  0.0421366   4.799 1.60e-06 ***
# age_group8             0.3495595  0.0623936   5.602 2.11e-08 ***
# variant1              -0.1721365  0.0123307 -13.960  < 2e-16 ***
# data_type1:age_group1 -0.0041364  0.0363757  -0.114 0.909465    
# data_type1:age_group2 -0.0029125  0.0298714  -0.098 0.922327    
# data_type1:age_group3 -0.0051792  0.0310816  -0.167 0.867659    
# data_type1:age_group4 -0.0244671  0.0300100  -0.815 0.414902    
# data_type1:age_group5  0.0118847  0.0288502   0.412 0.680380    
# data_type1:age_group6 -0.0058587  0.0309842  -0.189 0.850024    
# data_type1:age_group7 -0.0248670  0.0421565  -0.590 0.555276    
# data_type1:age_group8  0.0857293  0.0624058   1.374 0.169523    
# data_type1:variant1   -0.0033250  0.0111998  -0.297 0.766555    
# age_group1:variant1   -0.0251657  0.0137667  -1.828 0.067549 .  
# age_group2:variant1    0.0258927  0.0112398   2.304 0.021242 *  
# age_group3:variant1    0.0096242  0.0118996   0.809 0.418639    
# age_group4:variant1   -0.0294065  0.0115195  -2.553 0.010687 *  
# age_group5:variant1   -0.0380698  0.0110391  -3.449 0.000563 ***
# age_group6:variant1   -0.0289793  0.0118573  -2.444 0.014526 *  
# age_group7:variant1    0.0052756  0.0161309   0.327 0.743629    
# age_group8:variant1   -0.0078859  0.0246847  -0.319 0.749375    
# ---
#     Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)


Anova(fit, type="III")
# Analysis of Deviance Table (Type III tests)
# 
# Response: cbind(count_became_cases, total - count_became_cases)
#                     LR Chisq Df Pr(>Chisq)    
# data_type               0.00  1     0.9941    
# age_group             955.35  8  < 2.2e-16 ***
# variant               194.28  1  < 2.2e-16 ***
# data_type:age_group     2.90  8     0.9407    
# data_type:variant       0.09  1     0.7665    
# age_group:variant      36.74  8  1.283e-05 ***
table2csv(Anova(fit, type="III"), file=".\\multinomial_logistic_fits\\tables\\bGLMfit_secondaryattack_Anova type III tests.csv")

# effect plots
plot(allEffects(fit), ylab="Prop. of contacts with index case that became cases")



# 1.2 PLOT MODEL FIT ####

fitdSA_preds = as.data.frame(emmeans(fit, ~ variant, by=c("age_group","data_type"), 
                                    type="response"))
fitdSA_preds
fitdSA_preds_avg = as.data.frame(emmeans(fit, ~ variant, by=c("data_type"), 
                                         type="response"))
fitdSA_preds_avg = data.frame(variant=fitdSA_preds_avg[,1], age_group="avg",
                                       fitdSA_preds_avg[,2:7])

# fitdSA_preds$total = data_SA$total[match(interaction(fitdSA_preds$variant, fitdSA_preds$age_group, fitdSA_preds$data_type),
#                                         interaction(data_SA$variant, data_SA$age_group, data_SA$data_type))]
fitdSA_preds = rbind(fitdSA_preds, fitdSA_preds_avg)
fitdSA_preds$variant = as.character(fitdSA_preds$variant)
fitdSA_preds$variant = gsub("non_VOC","other",fitdSA_preds$variant)
fitdSA_preds$variant = factor(fitdSA_preds$variant, levels=c("VOC","other"))
table2csv(fitdSA_preds, file=".\\multinomial_logistic_fits\\tables\\bGLM_secondary_attack_rates_model predictions.csv")

plot_fitdSA = qplot(data=fitdSA_preds, x=age_group, y=prob, group=variant, geom="blank") +
    facet_wrap(~data_type, ncol=1) +
    # geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL, fill=variant), alpha=I(0.3)) +
    # geom_col(aes(y=prob, colour=NULL, 
    #             fill=variant), position = position_dodge(width = 0.9), alpha=I(1)) +
    geom_point(aes(y=prob, colour=variant, # size=total, 
                 fill=variant), position = position_dodge(width = 0.2), alpha=I(20)) + # position = position_dodge(width = 0.9)
    geom_linerange(aes(ymin=asymp.LCL, ymax=asymp.UCL, colour=variant), position = position_dodge(width = 0.2)) + # position = position_dodge(width = 0.9)
    geom_line(data=fitdSA_preds[fitdSA_preds$age_group!="avg",], aes(y=prob, colour=variant), 
              position = position_dodge(width = 0.2), alpha=I(1)) +
    ylab("Contacts with index patient that became cases (%)") +
    xlab("Age of person being infected") +
    theme_hc() + 
    # ggtitle("SECONDARY ATTACK RATES") +
    scale_y_continuous( trans="logit",   
                        breaks=c(0.01, 0.05, 0.1, 0.15, 0.20, 0.25, 0.30, 0.35),
                        labels = c("1","5","10","15","20","25", "30", "35")) +
    coord_cartesian(ylim=c(0.05,0.30), xlim=c(0.5,10.5), expand=c(0,0)) +
    scale_color_manual("variant carried\nby index patient", values=c("red","blue","red")) +
    scale_fill_manual("variant carried\nby index patient", values=c("red","blue")) +
    # geom_point(data=data_SA, 
    #           aes(y=prop, colour=variant, size=total), alpha=I(0.5)) +
    scale_size_continuous("all contacts", trans="log10" 
                          # range=c(0.01, 4), limits=c(100,1600), breaks=c(100,200,400,800,1600)
    ) +
    guides(size=FALSE) + 
    # guides(colour=FALSE) +
    theme(legend.position = c(0.85, 0.1)) + 
    theme(
        legend.title = element_text(color = "black", size = 8+1),
        legend.text = element_text(color = "black", size = 7+1),
        legend.key.size = unit(0.4, "cm")
    ) # +
    # labs(tag = "@TWenseleers\ndata PHE SGSS, NHS T&T, COG-UK") +
    # theme(plot.tag.position = "bottomright",
    #      plot.tag = element_text(vjust = 1, hjust = 1, size=8))
# theme(legend.position = "right") # theme(legend.position = "none")
plot_fitdSA
ggsave(file=".\\multinomial_logistic_fits\\plots\\bGLM_plot secondary_attack_rates_VOC nonVOC_byage.png", width=7, height=5)
saveRDS(plot_fitdSA, file = ".\\multinomial_logistic_fits\\plots\\bGLM_plot secondary_attack_rates_VOC nonVOC_byage.rds")
graph2ppt(file=".\\multinomial_logistic_fits\\plots\\bGLM_plot secondary_attack_rates_VOC nonVOC_byage.pptx", width=6, height=6)


# 1.3 SIDAK POSTHOC TESTS TO TEST FOR DIFFERENCES IN SECONDARY INFECTIONS BEING BY VOC VS NON-VOC IFO AGE ####

# Sidak posthoc tests to compare odds ratio to be infected by VOC ifo age against H0 of OR=1 :

# all p values for odds ratios to be sign different from 1 are highly significant, except for 80+
# which implies that VOC is indeed sign more infectious
# p values are Sidak corrected for multiple testing
# odds ratios:
fitSA_contrasts = data.frame(as.data.frame(emmeans(fit, revpairwise ~ variant|age_group, by="data_type",
                               adjust="sidak", type="response")$contrasts),
                             as.data.frame(confint(emmeans(fit, revpairwise ~ variant|age_group, by="data_type",
                                                           adjust="sidak", type="response")$contrasts))[,c("asymp.LCL","asymp.UCL")]) 
fitSA_contrasts_avg = data.frame(as.data.frame(emmeans(fit, revpairwise ~ variant, by="data_type",
                                                            adjust="sidak", type="response")$contrasts),
                                 as.data.frame(confint(emmeans(fit, revpairwise ~ variant, by="data_type",
                                                           adjust="sidak", type="response")$contrasts))[,c("asymp.LCL","asymp.UCL")]) 
fitSA_contrasts_avg = data.frame(contrast=fitSA_contrasts_avg[,1], age_group="avg",
                                 fitSA_contrasts_avg[,2:9])

fitSA_contrasts = rbind(fitSA_contrasts, fitSA_contrasts_avg)

fitSA_contrasts_avg_bothdatatypes = data.frame(as.data.frame(emmeans(fit, revpairwise ~ variant, 
                                                       adjust="sidak", type="response")$contrasts),
                                 as.data.frame(confint(emmeans(fit, revpairwise ~ variant, 
                                                               adjust="sidak", type="response")$contrasts))[,c("asymp.LCL","asymp.UCL")]) 
fitSA_contrasts_avg_bothdatatypes
#        contrast odds.ratio         SE  df  z.ratio     p.value asymp.LCL asymp.UCL
# 1 VOC / non_VOC   1.410964 0.03479632 Inf 13.96001 2.73413e-44  1.344386  1.480838

# difference in log(odds):
fitSA_contrasts
#         contrast age_group  data_type odds.ratio         SE  df   z.ratio       p.value asymp.LCL asymp.UCL
# 1  VOC / non_VOC       0-9 sequencing   1.493697 0.07628370 Inf  7.856892  3.937820e-15 1.3514234  1.650950
# 2  VOC / non_VOC     10-19 sequencing   1.348695 0.06448192 Inf  6.256723  3.931507e-10 1.2280537  1.481188
# 3  VOC / non_VOC     20-29 sequencing   1.393300 0.06763318 Inf  6.832774  8.328803e-12 1.2668515  1.532369
# 4  VOC / non_VOC     30-39 sequencing   1.506420 0.07239708 Inf  8.525688  1.519041e-17 1.3710025  1.655214
# 5  VOC / non_VOC     40-49 sequencing   1.532749 0.07279668 Inf  8.991894  2.430070e-19 1.3965094  1.682279
# 6  VOC / non_VOC     50-59 sequencing   1.505134 0.07298370 Inf  8.432319  3.388830e-17 1.3686754  1.655197
# 7  VOC / non_VOC     60-69 sequencing   1.405470 0.07657451 Inf  6.247281  4.176592e-10 1.2631224  1.563860
# 8  VOC / non_VOC     70-79 sequencing   1.442958 0.09944993 Inf  5.320517  1.034725e-07 1.2606310  1.651654
# 9  VOC / non_VOC       80+ sequencing   1.189453 0.11497549 Inf  1.794841  7.267902e-02 0.9841655  1.437562
# 10 VOC / non_VOC       0-9       SGTF   1.473962 0.04029133 Inf 14.192389  1.021203e-45 1.3970711  1.555086
# 11 VOC / non_VOC     10-19       SGTF   1.330876 0.02734437 Inf 13.911979  5.357605e-44 1.2783469  1.385564
# 12 VOC / non_VOC     20-29       SGTF   1.374891 0.03077337 Inf 14.224321  6.472913e-46 1.3158803  1.436548
# 13 VOC / non_VOC     30-39       SGTF   1.486517 0.03171638 Inf 18.580586  4.614410e-77 1.4256362  1.549998
# 14 VOC / non_VOC     40-49       SGTF   1.512498 0.03022817 Inf 20.703046  3.251550e-95 1.4543973  1.572920
# 15 VOC / non_VOC     50-59       SGTF   1.485248 0.03307166 Inf 17.765566  1.306169e-70 1.4218226  1.551502
# 16 VOC / non_VOC     60-69       SGTF   1.386901 0.04622484 Inf  9.813253  9.873368e-23 1.2991977  1.480525
# 17 VOC / non_VOC     70-79       SGTF   1.423893 0.07677047 Inf  6.554554  5.580841e-11 1.2811030  1.582598
# 18 VOC / non_VOC       80+       SGTF   1.173738 0.10154342 Inf  1.851674  6.407261e-02 0.9906753  1.390628
# 19 VOC / non_VOC       avg sequencing   1.420378 0.06412250 Inf  7.773301  7.646654e-15 1.3001000  1.551784
# 20 VOC / non_VOC       avg       SGTF   1.401612 0.01889904 Inf 25.039168 2.290923e-138 1.3650557  1.439147

table2csv(fitSA_contrasts, file=".\\multinomial_logistic_fits\\tables\\bGLMfit_secondaryattack_Sidak posthoc tests odds VOC vs nonVOC by age.csv")



# plot of odds ratio to be infected by VOC ifo age ####
plot_fitSA_oddsratios = qplot(data=fitSA_contrasts, x=age_group, y=odds.ratio, group=data_type, 
                                     geom="blank") +
    facet_wrap(~data_type, ncol=1) +
    # geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL, fill=variant), alpha=I(0.3)) +
    # geom_col(aes(y=prob, colour=NULL, 
    #             fill=variant), position = position_dodge(width = 0.9), alpha=I(1)) +
    geom_point(aes(y=odds.ratio), fill=I("red"), colour=I("red"), alpha=I(1)) + # position = position_dodge(width = 0.9)
    geom_linerange(aes(ymin=asymp.LCL, ymax=asymp.UCL), fill=I("red"), colour=I("red")) + # position = position_dodge(width = 0.9)
    geom_line(data=fitSA_contrasts[fitSA_contrasts$age_group!="avg",],
              aes(x=age_group, y=odds.ratio, fill=I("red"), colour=I("red")), alpha=I(1)) +
    geom_hline(yintercept=1, color=alpha("black", 0.5)) +
    ylab("Odds ratio to be infected by VOC") +
    xlab("Age of person being infected") +
    theme_hc() + 
    # ggtitle("ODDS RATIO TO BE INFECTED BY VOC 202012/01 BY AGE") +
    scale_y_continuous( # trans="log" # ,   
        # breaks=c(0.01, 0.05, 0.1, 0.15, 0.20, 0.25),
        # labels = c("1","5","10","15","20","25")
    ) +
    theme(legend.position = "right") + # theme(legend.position = "none")
    coord_cartesian(ylim=c(1, 1.7), 
                     xlim=c(0.5, 10.5), expand=c(0,0)) # +
    # labs(tag = "@TWenseleers\ndata PHE SGSS, NHS T&T, COG-UK") +
    # theme(plot.tag.position = "bottomright",
    #        plot.tag = element_text(vjust = 1, hjust = 1, size=8))
plot_fitSA_oddsratios
ggsave(file=".\\multinomial_logistic_fits\\plots\\bGLM_plot secondary_attack_rates_odds ratios VOC_byage.png", width=7, height=5)
saveRDS(plot_fitdSA, file = ".\\multinomial_logistic_fits\\plots\\bGLM_plot secondary_attack_rates_odds ratios VOC_byage.rds")
graph2ppt(file=".\\multinomial_logistic_fits\\plots\\bGLM_plot secondary_attack_rates_odds ratios VOC_byage.pptx", width=6, height=6)


# plot of relative risk to be infected by VOC ifo age ####
# cf e.g. https://www.flutterbys.com.au/stats/tut/tut10.5a.html

# function to convert ODDS RATIO to RELATIVE RISK 
# with PO=mean prop that became infected with VOC for each 
# particular age category, or on average
OR_to_RR = function(OR, P0) { RR = OR / (1 - P0 + (P0 * OR))
                       return(RR) }
newdat = data.frame(fitSA_contrasts[-tail(1:nrow(fitSA_contrasts),2),
                                    c("age_group","data_type")],
                    variant="VOC")
P0 = predict(fit, newdata = newdat, type="response")
P0_link = predict(fit, newdata = newdat, type="link")
fitSA_contrasts$P0 = c(P0, plogis(aggregate(P0_link, by=list(data_type=newdat$data_type), mean)[,2])) 
fitSA_contrasts$RR = OR_to_RR(fitSA_contrasts$odds.ratio, fitSA_contrasts$P0)
fitSA_contrasts$RR.LCL = OR_to_RR(fitSA_contrasts$asymp.LCL, fitSA_contrasts$P0)
fitSA_contrasts$RR.UCL = OR_to_RR(fitSA_contrasts$asymp.UCL, fitSA_contrasts$P0)
fitSA_contrasts
#         contrast age_group  data_type odds.ratio         SE  df   z.ratio       p.value asymp.LCL asymp.UCL         P0       RR    RR.LCL
# 1  VOC / non_VOC       0-9 sequencing   1.493697 0.07628370 Inf  7.856892  3.937820e-15 1.3514234  1.650950 0.08900027 1.430828 1.3104371
# 2  VOC / non_VOC     10-19 sequencing   1.348695 0.06448192 Inf  6.256723  3.931507e-10 1.2280537  1.481188 0.12141094 1.293917 1.1949672
# 3  VOC / non_VOC     20-29 sequencing   1.393300 0.06763318 Inf  6.832774  8.328803e-12 1.2668515  1.532369 0.17588954 1.303151 1.2100558
# 4  VOC / non_VOC     30-39 sequencing   1.506420 0.07239708 Inf  8.525688  1.519041e-17 1.3710025  1.655214 0.19141462 1.373298 1.2800961
# 5  VOC / non_VOC     40-49 sequencing   1.532749 0.07279668 Inf  8.991894  2.430070e-19 1.3965094  1.682279 0.20460271 1.382098 1.2917165
# 6  VOC / non_VOC     50-59 sequencing   1.505134 0.07298370 Inf  8.432319  3.388830e-17 1.3686754  1.655197 0.20229001 1.365593 1.2736848
# 7  VOC / non_VOC     60-69 sequencing   1.405470 0.07657451 Inf  6.247281  4.176592e-10 1.2631224  1.563860 0.20079675 1.299656 1.1997354
# 8  VOC / non_VOC     70-79 sequencing   1.442958 0.09944993 Inf  5.320517  1.034725e-07 1.2606310  1.651654 0.24783375 1.300220 1.1841434
# 9  VOC / non_VOC       80+ sequencing   1.189453 0.11497549 Inf  1.794841  7.267902e-02 0.9841655  1.437562 0.18714069 1.148726 0.9870905
# 10 VOC / non_VOC       0-9       SGTF   1.473962 0.04029133 Inf 14.192389  1.021203e-45 1.3970711  1.555086 0.08915266 1.414205 1.3493058
# 11 VOC / non_VOC     10-19       SGTF   1.330876 0.02734437 Inf 13.911979  5.357605e-44 1.2783469  1.385564 0.12135019 1.279502 1.2365784
# 12 VOC / non_VOC     20-29       SGTF   1.374891 0.03077337 Inf 14.224321  6.472913e-46 1.3158803  1.436548 0.17646483 1.289579 1.2464035
# 13 VOC / non_VOC     30-39       SGTF   1.486517 0.03171638 Inf 18.580586  4.614410e-77 1.4256362  1.549998 0.19808524 1.355851 1.3147838
# 14 VOC / non_VOC     40-49       SGTF   1.512498 0.03022817 Inf 20.703046  3.251550e-95 1.4543973  1.572920 0.19973751 1.372048 1.3333795
# 15 VOC / non_VOC     50-59       SGTF   1.485248 0.03307166 Inf 17.765566  1.306169e-70 1.4218226  1.551502 0.20315030 1.351973 1.3095987
# 16 VOC / non_VOC     60-69       SGTF   1.386901 0.04622484 Inf  9.813253  9.873368e-23 1.2991977  1.480525 0.20784197 1.283675 1.2231359
# 17 VOC / non_VOC     70-79       SGTF   1.423893 0.07677047 Inf  6.554554  5.580841e-11 1.2811030  1.582598 0.21618261 1.304363 1.2077109
# 18 VOC / non_VOC       80+       SGTF   1.173738 0.10154342 Inf  1.851674  6.407261e-02 0.9906753  1.390628 0.19549286 1.135182 0.9924845
# 19 VOC / non_VOC       avg sequencing   1.420378 0.06412250 Inf  7.773301  7.646654e-15 1.3001000  1.551784 0.17459950 1.323254 1.2353699
# 20 VOC / non_VOC       avg       SGTF   1.401612 0.01889904 Inf 25.039168 2.290923e-138 1.3650557  1.439147 0.17367985 1.310221 1.2836675
# RR.UCL
# 1  1.560540
# 2  1.399431
# 3  1.401166
# 4  1.470755
# 5  1.476207
# 6  1.461491
# 7  1.404806
# 8  1.421999
# 9  1.328756
# 10 1.481757
# 11 1.323633
# 12 1.333799
# 13 1.397721
# 14 1.411407
# 15 1.395188
# 16 1.346087
# 17 1.405570
# 18 1.291967
# 19 1.415420
# 20 1.337161
table2csv(fitSA_contrasts, file=".\\multinomial_logistic_fits\\tables\\bGLM_secondary_attack_rates_odds ratios and relative risks.csv")

plot_fitSA_oddsratios_RR = qplot(data=fitSA_contrasts, x=age_group, y=RR, group=data_type, 
                                     geom="blank") +
    facet_wrap(~data_type, ncol=1) +
    # geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL, fill=variant), alpha=I(0.3)) +
    # geom_col(aes(y=prob, colour=NULL, 
    #             fill=variant), position = position_dodge(width = 0.9), alpha=I(1)) +
    geom_point(aes(y=RR), fill=I("red"), colour=I("red"), alpha=I(1)) + # position = position_dodge(width = 0.9)
    geom_linerange(aes(ymin=RR.LCL, ymax=RR.UCL), fill=I("red"), colour=I("red")) + # position = position_dodge(width = 0.9)
    geom_line(data=fitSA_contrasts[fitSA_contrasts$age_group!="avg",],
              aes(x=age_group, y=RR, fill=I("red"), colour=I("red")), alpha=I(1)) +
    geom_hline(yintercept=1, color=alpha("black", 0.5)) +
    ylab("Increased risk to be infected by VOC (relative risk)") +
    xlab("Age of person being infected") +
    theme_hc() + 
#    ggtitle("RELATIVE RISK TO BE INFECTED BY VOC 202012/01 BY AGE") +
    coord_cartesian(# ylim=c(0, 60), 
                    xlim=c(0.5, 10.5), expand=c(0,0)) +
    theme(legend.position = "right") # + # theme(legend.position = "none")
 # labs(tag = "@TWenseleers\ndata PHE SGSS, NHS T&T, COG-UK") +
 # theme(plot.tag.position = "bottomright",
 #       plot.tag = element_text(vjust = 1, hjust = 1, size=8))
plot_fitSA_oddsratios_RR
ggsave(file=".\\multinomial_logistic_fits\\plots\\bGLM_plot secondary_attack_rates_risk ratios VOC_byage.png", width=7, height=5)
saveRDS(plot_fitdSA, file = ".\\multinomial_logistic_fits\\plots\\bGLM_plot secondary_attack_rates_risk ratios VOC_byage.rds")
graph2ppt(file=".\\multinomial_logistic_fits\\plots\\bGLM_plot secondary_attack_rates_risk ratios VOC_byage.pptx", width=6, height=6)


plot_fitSA_multipanel = ggarrange(plot_fitdSA+theme(
    legend.title = element_text(color = "black", size = 8),
    legend.text = element_text(color = "black", size = 7),
    legend.key.size = unit(0.3, "cm")
)+ggtitle("")+xlab("")+ylab("Secondary attack rate (%)")+
    theme(legend.position = c(0.85, 0.135)) + 
    theme(
        legend.title = element_text(color = "black", size = 7),
        legend.text = element_text(color = "black", size = 6),
        legend.key.size = unit(0.3, "cm")), plot_fitSA_oddsratios+ggtitle(""), 
   nrow=2, common.legend=FALSE) 
plot_fitSA_multipanel
saveRDS(plot_fitSA_multipanel, file = ".\\multinomial_logistic_fits\\plots\\bGLM_plot secondary_attack_rates_multipanel.rds")
ggsave(file=".\\multinomial_logistic_fits\\plots\\bGLM_plot secondary_attack_rates_multipanel.png", width=6, height=7)
graph2ppt(file=".\\multinomial_logistic_fits\\plots\\bGLM_plot secondary_attack_rates_multipanel.pptx", width=6, height=7)




# 1.4 CONTRASTS (DIFF IN LOG(ODDS RATIOS) TO BE INFECTED BY VOC VS BY NON-VOC ACROSS DIFFERENT AGE GROUPS ####

# comparison of log(odds ratios) to be infected by VOC with average log(odds ratios) across all age groups
# ie look at whether log(odds ratios) to be infected by VOC are above or below average for each age group
# p values are Sidak adjusted

# = age x variant interaction contrasts testing for above-agerage susceptibility for VOC for each age category
# (p values Sidak corrected for multiple testing)

contrastsLOGORVOC_agegroups_eff = cbind(as.data.frame(contrast(emmeans(fit,  ~ age_group*variant), interaction="eff", adjust="sidak")),
                                  as.data.frame(confint(contrast(emmeans(fit,  ~ age_group*variant), interaction="eff", adjust="sidak")))[,c("asymp.LCL","asymp.UCL")])
contrastsLOGORVOC_agegroups_eff = contrastsLOGORVOC_agegroups_eff[contrastsLOGORVOC_agegroups_eff$variant_eff=="VOC effect",]
contrastsLOGORVOC_agegroups_eff
#     age_group_eff variant_eff     estimate         SE  df    z.ratio    p.value    asymp.LCL  asymp.UCL
# 10   (0-9) effect  VOC effect  0.025165698 0.01376674 Inf  1.8280067 0.71603023 -0.015913518 0.06624491
# 11 (10-19) effect  VOC effect -0.025892748 0.01123982 Inf -2.3036612 0.32054845 -0.059431775 0.00764628
# 12 (20-29) effect  VOC effect -0.009624199 0.01189957 Inf -0.8087853 0.99994245 -0.045131876 0.02588348
# 13 (30-39) effect  VOC effect  0.029406525 0.01151948 Inf  2.5527658 0.17585082 -0.004966969 0.06378002
# 14 (40-49) effect  VOC effect  0.038069832 0.01103908 Inf  3.4486428 0.01009299  0.005129825 0.07100984
# 15 (50-59) effect  VOC effect  0.028979332 0.01185735 Inf  2.4439976 0.23154783 -0.006402356 0.06436102
# 16 (60-69) effect  VOC effect -0.005275627 0.01613089 Inf -0.3270511 1.00000000 -0.053409339 0.04285809
# 17 (70-79) effect  VOC effect  0.007885860 0.02468474 Inf  0.3194630 1.00000000 -0.065772054 0.08154378
# 18   (80+) effect  VOC effect -0.088714674 0.03872377 Inf -2.2909615 0.32953820 -0.204264315 0.02683497
# i.e. after Sidak p value adjustment only 40-49 slightly more susceptible (p=0.01) for VOC than average

# PS estimate=difference in log(odds ratios)
# this significance could be indicated with an asterisk on the odds ratio plot above

table2csv(contrastsLOGORVOC_agegroups_eff, file=".\\multinomial_logistic_fits\\tables\\bGLMfit_secondaryattack_Sidak posthoc tests age variant interaction contrasts_diff log odds ratio.csv")



# PS: this would be all pairwise comparisons between age categories in terms of
# difference in log(odds ratios) to be infected by VOC

contrastsLOGORVOC_agegroups = as.data.frame(contrast(emmeans(fit,  ~ age_group*variant), interaction="revpairwise", adjust="sidak"))
contrastsLOGORVOC_agegroups

# age_group_revpairwise variant_revpairwise     estimate         SE  df     z.ratio      p.value
# 1        (10-19) - (0-9)       VOC - non_VOC -0.102116891 0.03411972 Inf -2.99289915 0.0948191199
# 2        (20-29) - (0-9)       VOC - non_VOC -0.069579793 0.03525165 Inf -1.97380232 0.8323937603
# 3      (20-29) - (10-19)       VOC - non_VOC  0.032537098 0.03029362 Inf  1.07405768 0.9999936449
# 4        (30-39) - (0-9)       VOC - non_VOC  0.008481654 0.03459626 Inf  0.24516103 1.0000000000
# 5      (30-39) - (10-19)       VOC - non_VOC  0.110598545 0.02952840 Inf  3.74549799 0.0064609279
# 6      (30-39) - (20-29)       VOC - non_VOC  0.078061447 0.03082914 Inf  2.53206682 0.3367110933
# 7        (40-49) - (0-9)       VOC - non_VOC  0.025808268 0.03378114 Inf  0.76398461 0.9999999994
# 8      (40-49) - (10-19)       VOC - non_VOC  0.127925159 0.02856902 Inf  4.47775762 0.0002715164
# 9      (40-49) - (20-29)       VOC - non_VOC  0.095388061 0.02991158 Inf  3.18900086 0.0501320321
# 10     (40-49) - (30-39)       VOC - non_VOC  0.017326614 0.02913629 Inf  0.59467458 1.0000000000
# 11       (50-59) - (0-9)       VOC - non_VOC  0.007627269 0.03517840 Inf  0.21681680 1.0000000000
# 12     (50-59) - (10-19)       VOC - non_VOC  0.109744160 0.03020835 Inf  3.63290764 0.0100395106
# 13     (50-59) - (20-29)       VOC - non_VOC  0.077207062 0.03148105 Inf  2.45249309 0.4021347536
# 14     (50-59) - (30-39)       VOC - non_VOC -0.000854385 0.03074536 Inf -0.02778907 1.0000000000
# 15     (50-59) - (40-49)       VOC - non_VOC -0.018180999 0.02982523 Inf -0.60958463 1.0000000000
# 16       (60-69) - (0-9)       VOC - non_VOC -0.060882648 0.04304242 Inf -1.41448027 0.9978831282
# 17     (60-69) - (10-19)       VOC - non_VOC  0.041234242 0.03908545 Inf  1.05497670 0.9999958918
# 18     (60-69) - (20-29)       VOC - non_VOC  0.008697145 0.04007729 Inf  0.21700931 1.0000000000
# 19     (60-69) - (30-39)       VOC - non_VOC -0.069364302 0.03950203 Inf -1.75596814 0.9485063769
# 20     (60-69) - (40-49)       VOC - non_VOC -0.086690916 0.03879015 Inf -2.23486917 0.6043287200
# 21     (60-69) - (50-59)       VOC - non_VOC -0.068509917 0.04001288 Inf -1.71219676 0.9620389580
# 22       (70-79) - (0-9)       VOC - non_VOC -0.034559674 0.06040024 Inf -0.57217780 1.0000000000
# 23     (70-79) - (10-19)       VOC - non_VOC  0.067557216 0.05764739 Inf  1.17190418 0.9999517106
# 24     (70-79) - (20-29)       VOC - non_VOC  0.035020119 0.05832416 Inf  0.60043930 1.0000000000
# 25     (70-79) - (30-39)       VOC - non_VOC -0.043041328 0.05793035 Inf -0.74298405 0.9999999997
# 26     (70-79) - (40-49)       VOC - non_VOC -0.060367942 0.05744735 Inf -1.05083944 0.9999962697
# 27     (70-79) - (50-59)       VOC - non_VOC -0.042186943 0.05827992 Inf -0.72386751 0.9999999999
# 28     (70-79) - (60-69)       VOC - non_VOC  0.026322974 0.06333789 Inf  0.41559600 1.0000000000
# 29         (80+) - (0-9)       VOC - non_VOC -0.227760743 0.09069905 Inf -2.51117004 0.3532669354
# 30       (80+) - (10-19)       VOC - non_VOC -0.125643853 0.08888948 Inf -1.41348391 0.9979094219
# 31       (80+) - (20-29)       VOC - non_VOC -0.158180950 0.08933006 Inf -1.77074714 0.9432470276
# 32       (80+) - (30-39)       VOC - non_VOC -0.236242397 0.08907346 Inf -2.65221982 0.2510096786
# 33       (80+) - (40-49)       VOC - non_VOC -0.253569011 0.08876005 Inf -2.85679214 0.1430632335
# 34       (80+) - (50-59)       VOC - non_VOC -0.235388012 0.08930118 Inf -2.63588911 0.2616788983
# 35       (80+) - (60-69)       VOC - non_VOC -0.166878095 0.09268134 Inf -1.80055758 0.9315210433
# 36       (80+) - (70-79)       VOC - non_VOC -0.193201069 0.10190861 Inf -1.89582683 0.8835564266
