# MULTINOMIAL AND LOGISTIC FITS TO DETERMINE GROWTH RATE & 
# COMPETITIVE ADVANTAGE OF VOC 202012/01 STRAIN
# ACROSS DIFFERENT REGIONS IN THE UK AS WELL AS DENMARK
# T. Wenseleers, 7 Jan. 2021

library(lme4)
library(emmeans)
library(effects)
library(ggplot2)
library(ggthemes)
library(MASS)
library(nlme)
library(mgcv)
library(tidyr)
library(nnet)
library(splines)
library(ggpubr)
# install from https://github.com/tomwenseleers/export
# library(devtools)
# devtools::install_github("tomwenseleers/export")
library(export) 



# FUNCTION TO CALCULATE DIFFERENCE IN R VALUES BETWEEN TWO COMPETING STRAINS
# delta_R FROM THE DIFFERENCE IN THEIR GROWTH RATE delta_r

# GIVEN THAT THE EFFECTIVE REPRODUCTION NUMBER R = 1 + r*g WITH r=GROWTH RATE 
# AND g=GENERATION TIME (EQN 3.1 IN WALLINGA & LIPSITCH PROC B 2007), WE CAN 
# CALCULATE THE DIFFERENCE IN EFFECTIVE REPRODUCTION NUMBER delta_R 
# (ASSUMING IDENTICAL GENERATION TIMES FOR BOTH) AS delta_r*g, OR, IF WE
# HAVE AN EMPIRICAL DISTRIBUTION FOR GENERATION TIME, AS THE WEIGHTED PRODUCT
# OF BOTH.

delta_R.from.delta_r = function (delta_r, mean=5.5, sd=2.1) { 
    require(R0)
    if (sd==0) { delta_R = delta_r*mean } else {
                    time = seq(0, mean+10*sd, by=0.001) # we use a very small time step to calculate weighted avg
                    dens = dgamma(seq(0, mean+10*sd, by=0.001), 
                                  shape = (mean/sd)^2, 
                                  scale = sd^2/mean)/100
                    delta_R = sapply(delta_r, function (delta_r) weighted.mean(x=time*delta_r, w=dens)) }
    return( delta_R ) 
           # weighted average over given gamma generation time distribution
}

# This function calculates delta_R.from.delta_R for two sets of defaults for
# generation time : gamma(mean=5.5d, SD=2.1d) (Ferretti et al. 2020) or
# gamma(mean=3.6d, SD=3.1d) (Abbott et al. 2020, Ganyani et al. 2020).
# It works on an input dataframe df with delta_r values and returns the original
# data frame plus two sets of estimates for delta_R as a dataframe with extra columns
# with column names coln
delta_R.from.delta_r_df = function (df, mean1=5.5, sd1=2.1, mean2=3.6, sd2=3.1,
                                    coln=c("delta_R1","delta_R1.LCL","delta_R1.UCL",
                                           "delta_R2","delta_R2.LCL","delta_R2.UCL")) { 
  df_num = df[,which(unlist(lapply(df, is.numeric))), drop=F]
  df_nonnum = df[,which(!unlist(lapply(df, is.numeric))), drop=F]
  df_out1 = apply(df_num, 2, function (delta_r) delta_R.from.delta_r(delta_r, mean=mean1, sd=sd1))
  if (class(df_out1)[1]=="numeric") df_out1=as.data.frame(t(df_out1), check.names=F)
  df_out2 = apply(df_num, 2, function (delta_r) delta_R.from.delta_r(delta_r, mean=mean2, sd=sd2))
  if (class(df_out2)[1]=="numeric") df_out2=as.data.frame(t(df_out2), check.names=F)
  df_out = data.frame(df_out1, df_out2, check.names=F)
  if (!is.null(coln)) colnames(df_out) = coln
  return( data.frame(df_nonnum, df_num, df_out, check.names=F) )
}





# 1. LOAD & PREPROCESS DATA ####
# desired order of factor levels
levels_nhs_name = c("South East","London","East of England",
                    "South West","Midlands","North East and Yorkshire",
                    "Scotland","North West","Wales")
levels_variant_lineages = c("B","B.1.98","B.40","B.1.1",
                     "B.1.1.257","B.1.1.1","B.1.1.315",
                     "20A.EU1","20B.501Y.V1", # "501Y.WALES", 
                     "other")
# colours to use for variant lineages
n = length(levels_variant_lineages)
lineage_cols = hcl(h = seq(15, 375, length = n + 1), l = 65, c = 200)[1:n]
lineage_cols[which(levels_variant_lineages=="other")] = "grey75"

# read in data
data = read.csv(".\\data\\cog_metadata_microreact_public-2020-12-22-annotated.csv")
data$sample_date = as.Date(data$sample_date)
data$sample_date_num = as.numeric(data$sample_date)
# we do not consider data from Northern Ireland due to low nr of sequences there & absence of new VOC
data = data[data$nhs_name!="Northern Ireland",] 
data$nhs_name = factor(data$nhs_name, levels=levels_nhs_name)
data$lad = as.factor(data$lad) # local authority
# data$B_1_1_7 = (data$lineage=="B.1.1.7")&(data$n501y == "Y") # = 20B/501Y.V1 501Y.V1
unique(data[data$n501y == "Y","lineage"]) # lineages where at least some samples have n501y mutation
# "B.1.1.7"   "B.1.1.70"  "B.1"       "B.1.177"   "B.1.83"    "B.1.1"     "B.1.1.136"
data$week = strftime(data$sample_date, format = "%V") # week number
data$variant_lineage = data$lineage # we slightly recode some variant lineages, e.g. 20B.501Y.V1 = data$lineage=="B.1.1.7"&data$n501y == "Y"
data$variant_lineage[grepl("B.1.177",data$lineage,fixed=T)]="20A.EU1" # B.1.177 & descendant lineages = 20A.EU1 
data$variant_lineage[data$lineage=="B.1.1.7"&data$n501y == "Y"]="20B.501Y.V1" # =VOC 202012/01
data$variant_lineage[data$lineage=="B.1.1.70"&data$n501y == "Y"]="501Y.WALES" # Welsh variant, no official name yet


# aggregate data by week and check which variant lineages reached at least 15% in some week    
data_agbyweek = as.data.frame(table(data$week, data$variant_lineage))
colnames(data_agbyweek) = c("week", "variant_lineage", "count")
data_agbyweek_sum = aggregate(count ~ week, data=data_agbyweek, sum)
data_agbyweek$total = data_agbyweek_sum$count[match(data_agbyweek$week, data_agbyweek_sum$week)]
sum(data_agbyweek[data_agbyweek$variant_lineage=="20B.501Y.V1","total"]) == nrow(data) # correct
data_agbyweek$sample_date = as.Date(paste(2020, data_agbyweek$week, 1, sep="-"), "%Y-%U-%u")
data_agbyweek$variant_lineage = as.factor(data_agbyweek$variant_lineage)
data_agbyweek$prop = data_agbyweek$count/data_agbyweek$total
maxweeklyprop_lineages = aggregate(prop ~ variant_lineage, data=data_agbyweek, max)
# we select lineages that at one point in time reached at least a relative abundance of 15%
selectedlineages = as.character(maxweeklyprop_lineages$variant_lineage[maxweeklyprop_lineages$prop>=0.15])
selectedlineages
# 20A.EU1     20B.501Y.V1 B           B.1.1       B.1.1.1     B.1.1.257   B.1.1.315   B.1.98      B.40
# selectedlineages = as.factor(c(selectedlineages, "501Y.WALES")) # to also include the Welsh 501Y variant
# we recode all other 415 minor lineages as a single category "other"
data$variant_lineage[!data$variant_lineage %in% selectedlineages] = "other"
length(unique(data$lineage[data$variant_lineage=="other"])) # 415 lineages
data$variant_lineage = factor(data$variant_lineage, levels=levels_variant_lineages)

# aggregated data to make Muller plots of raw data
# aggregated by week for selected variant lineages
data_agbyweek = as.data.frame(table(data$week, data$variant_lineage))
colnames(data_agbyweek) = c("week", "variant_lineage", "count")
data_agbyweek_sum = aggregate(count ~ week, data=data_agbyweek, sum)
data_agbyweek$total = data_agbyweek_sum$count[match(data_agbyweek$week, data_agbyweek_sum$week)]
sum(data_agbyweek[data_agbyweek$variant_lineage=="20B.501Y.V1","total"]) == nrow(data) # correct
data_agbyweek$sample_date = as.Date(paste(2020, data_agbyweek$week, 1, sep="-"), "%Y-%U-%u")-3.5
data_agbyweek$variant_lineage = factor(data_agbyweek$variant_lineage, levels=levels_variant_lineages)
data_agbyweek$sample_date_num = as.numeric(data_agbyweek$sample_date)
data_agbyweek$prop = data_agbyweek$count/data_agbyweek$total

# aggregated by day for selected variant lineages
data_agbyday = as.data.frame(table(data$sample_date, data$variant_lineage))
colnames(data_agbyday) = c("sample_date", "variant_lineage", "count")
data_agbyday_sum = aggregate(count ~ sample_date, data=data_agbyday, sum)
data_agbyday$total = data_agbyday_sum$count[match(data_agbyday$sample_date, data_agbyday_sum$sample_date)]
sum(data_agbyday[data_agbyday$variant_lineage=="20B.501Y.V1","total"]) == nrow(data) # correct
data_agbyday$sample_date = as.Date(data_agbyday$sample_date)
data_agbyday$variant_lineage = factor(data_agbyday$variant_lineage, levels=levels_variant_lineages)
data_agbydayregion = data_agbydayregion[data_agbydayregion$total!=0,]
data_agbyday$sample_date_num = as.numeric(data_agbyday$sample_date)
data_agbyday$prop = data_agbyday$count/data_agbyday$total
data_agbyday = data_agbyday[data_agbyday$total!=0,]


# aggregated by week and nhs_name for selected variant lineages
data_agbyweekregion = as.data.frame(table(data$week, data$nhs_name, data$variant_lineage))
colnames(data_agbyweekregion) = c("week", "nhs_name", "variant_lineage", "count")
data_agbyweekregion_sum = aggregate(count ~ week + nhs_name, data=data_agbyweekregion, sum)
data_agbyweekregion$total = data_agbyweekregion_sum$count[match(interaction(data_agbyweekregion$week,data_agbyweekregion$nhs_name), 
                                                          interaction(data_agbyweekregion_sum$week,data_agbyweekregion_sum$nhs_name))]
sum(data_agbyweekregion[data_agbyweekregion$variant_lineage=="20B.501Y.V1","total"]) == nrow(data) # correct
data_agbyweekregion$sample_date = as.Date(paste(2020, data_agbyweekregion$week, 1, sep="-"), "%Y-%U-%u")-3.5
data_agbyweekregion$variant_lineage = factor(data_agbyweekregion$variant_lineage, levels=levels_variant_lineages)
data_agbyweekregion$nhs_name = factor(data_agbyweekregion$nhs_name, levels=levels_nhs_name)
data_agbyweekregion$sample_date_num = as.numeric(data_agbyweekregion$sample_date)
data_agbyweekregion$prop = data_agbyweekregion$count/data_agbyweekregion$total
data_agbyweekregion = data_agbyweekregion[data_agbyweekregion$total!=0,]

# aggregated by day and nhs_name for selected variant lineages
data_agbydayregion = as.data.frame(table(data$sample_date, data$nhs_name, data$variant_lineage))
colnames(data_agbydayregion) = c("sample_date", "nhs_name", "variant_lineage", "count")
data_agbydayregion_sum = aggregate(count ~ sample_date + nhs_name, data=data_agbydayregion, sum)
data_agbydayregion$total = data_agbydayregion_sum$count[match(interaction(data_agbydayregion$sample_date,data_agbydayregion$nhs_name),
                                                              interaction(data_agbydayregion_sum$sample_date,data_agbydayregion_sum$nhs_name))]
sum(data_agbydayregion[data_agbydayregion$variant_lineage=="20B.501Y.V1","total"]) == nrow(data) # correct
data_agbydayregion$sample_date = as.Date(data_agbydayregion$sample_date)
data_agbydayregion$variant_lineage = factor(data_agbydayregion$variant_lineage, levels=levels_variant_lineages)
data_agbydayregion$nhs_name = factor(data_agbydayregion$nhs_name, levels=levels_nhs_name)
data_agbydayregion$sample_date_num = as.numeric(data_agbydayregion$sample_date)
data_agbydayregion$prop = data_agbydayregion$count/data_agbydayregion$total
data_agbydayregion = data_agbydayregion[data_agbydayregion$total!=0,]

# data aggregated by day, nhs_name and lad for selected variant lineages
data_agbydayregionlad = as.data.frame(table(data$sample_date, data$nhs_name, data$lad, data$variant_lineage))
colnames(data_agbydayregionlad) = c("sample_date", "nhs_name", "lad", "variant_lineage", "count")
data_agbydayregionlad_sum = aggregate(count ~ sample_date + nhs_name + lad, data=data_agbydayregionlad, sum)
data_agbydayregionlad$total = data_agbydayregionlad_sum$count[match(interaction(data_agbydayregionlad$sample_date,data_agbydayregionlad$nhs_name,data_agbydayregionlad$lad),
                                                                    interaction(data_agbydayregionlad_sum$sample_date,data_agbydayregionlad_sum$nhs_name,data_agbydayregionlad_sum$lad))]
sum(data_agbydayregionlad[data_agbydayregionlad$variant_lineage=="20B.501Y.V1","total"]) == nrow(data) # correct
data_agbydayregionlad$sample_date = as.Date(data_agbydayregionlad$sample_date)
data_agbydayregionlad$variant_lineage = factor(data_agbydayregionlad$variant_lineage, levels=levels_variant_lineages)
data_agbydayregionlad$nhs_name = factor(data_agbydayregionlad$nhs_name, levels=levels_nhs_name)
data_agbydayregionlad$sample_date_num = as.numeric(data_agbydayregionlad$sample_date)
data_agbydayregionlad$prop = data_agbydayregionlad$count/data_agbydayregionlad$total
data_agbydayregionlad = data_agbydayregionlad[data_agbydayregionlad$total!=0,]
data_agbydayregionlad$obs = factor(1:nrow(data_agbydayregionlad))

# long version of data aggregated by day, nhs_name and lad, including all the zeros here
data_agbydayregionlad$variant_lineage = relevel(data_agbydayregionlad$variant_lineage, ref="other")
data_agbydayregionlad_wide = spread(data_agbydayregionlad, variant_lineage, count)
data_agbydayregionlad_wide[is.na(data_agbydayregionlad_wide)] = 0
sum(colSums(data_agbydayregionlad_wide[,c(8:17)])) == nrow(data) # correct
data_agbydayregionlad_wide$index = factor(1:nrow(data_agbydayregionlad_wide)) # for observation index random factor
data_agbydayregionlad_wide$TOTAL = rowSums(data_agbydayregionlad_wide[,levels_variant_lineages])
data_agbydayregionlad_long = gather(data_agbydayregionlad_wide, variant_lineage, 
                                    count, all_of(levels_variant_lineages), factor_key=TRUE)
data_agbydayregionlad_long$TOTAL = data_agbydayregionlad_wide$TOT[match(data_agbydayregionlad_long$index,
                                                                        data_agbydayregionlad_wide$index)]
sum(data_agbydayregionlad_long$TOTAL) # 1241160 = 10 times bigger than nrow(data), correct
nrow(data) # 1241160
data_agbydayregionlad_long = data_agbydayregionlad_long[data_agbydayregionlad_long$TOT!=0,]
data_agbydayregionlad_long$obs = factor(1:nrow(data_agbydayregionlad_long)) # for observation-level random factor
data_agbydayregionlad_long$variant_lineage = factor(data_agbydayregionlad_long$variant_lineage,
                                                    levels=levels_variant_lineages)
data_agbydayregionlad_long$variant_lineage = relevel(data_agbydayregionlad_long$variant_lineage,
                                                     ref="other")
nrow(data_agbydayregionlad_long) # 270210


# MULLER PLOTS OF RAW DATA ####

# Muller plots showing spread of main variant lineages
n = length(levels(data$variant_lineage))
lineage_cols = hcl(h = seq(15, 375, length = n + 1), l = 65, c = 200)[1:n]
lineage_cols[which(levels(data_agbyday$variant_lineage)=="other")] = "grey75"

# Muller plots using daily & weekly aggregated data (overall across all regions)
ggplot(data=data_agbyday, aes(x=sample_date, y=count, group=variant_lineage)) + 
    geom_area(aes(lwd=I(1.2), colour=NULL, fill=variant_lineage), position="fill") +
    scale_fill_manual("", values=lineage_cols) +
    scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
                       labels=substring(months(as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01"))),1,1),
                       limits=as.Date(c("2020-03-01","2020-12-22")), expand=c(0,0)) +
    guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
    theme_hc() + theme(legend.position="bottom", # c(0.8,0.7)) 
                       axis.title.x=element_blank(),
                       axis.title.y=element_blank()) + 
    labs(title = "MAIN SARS-CoV2 VARIANT LINEAGES IN THE UK") +
    ylab("Relative abundance")
muller_raw0 = ggplot(data=data_agbyweek, aes(x=sample_date, y=count, group=variant_lineage)) + 
    geom_area(aes(lwd=I(1.2), colour=NULL, fill=variant_lineage), position="fill") +
    scale_fill_manual("", values=lineage_cols) +
    scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
                       labels=substring(months(as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01"))),1,1),
                       limits=as.Date(c("2020-03-01","2020-12-22")), expand=c(0,0)) +
    guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
    theme_hc() + theme(legend.position="bottom", # c(0.8,0.7)) 
                       axis.title.x=element_blank()) + 
    labs(title = "MAIN SARS-CoV2 VARIANT LINEAGES IN THE UK") +
    ylab("Relative abundance")
muller_raw0
saveRDS(muller_raw0, file = ".\\multinomial_fits\\ggplot2_muller plot multinom lineage overall_raw data.rds")
graph2ppt(file=".\\multinomial_fits\\muller plot multinom lineage overall_raw data.pptx", width=7, height=5)
ggsave(file=".\\multinomial_fits\\muller plot multinom lineage overall_raw data.png", width=7, height=5)

muller_raw = ggplot(data=data_agbyweekregion, aes(x=sample_date, y=count, group=variant_lineage)) + 
    geom_area(aes(lwd=I(1.2), colour=NULL, fill=variant_lineage), position="fill") +
    facet_wrap(~nhs_name) +
    scale_fill_manual("", values=lineage_cols) +
    scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
                       labels=substring(months(as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01"))),1,1),
                       limits=as.Date(c("2020-03-01","2020-12-22")), expand=c(0,0)) +
    guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
    theme_hc() + theme(legend.position="bottom", # c(0.8,0.7)) 
                       axis.title.x=element_blank()) + 
    labs(title = "MAIN SARS-CoV2 VARIANT LINEAGES IN THE UK") +
    ylab("Relative abundance")
muller_raw
saveRDS(muller_raw, file = ".\\multinomial_fits\\muller plot multinom lineage region_raw data_raw data.rds")
graph2ppt(file=".\\multinomial_fits\\muller plot multinom lineage region_raw data.pptx", width=7, height=5)
ggsave(file=".\\multinomial_fits\\muller plot multinom lineage region_raw data.png", width=7, height=5)




# 2. MULTINOMIAL FITS FOR UK BY REGION ####

# we take the category "other" as reference category
data$variant_lineage = relevel(data$variant_lineage, ref="other") 
# PS if you would like to take a lineage with a low growth rate as a reference, one could use B.40

# we fit both a main effects model and models with date:nhs region interactions effect
set.seed(1)
fit1 = multinom(variant_lineage~sample_date_num+nhs_name, data = data, maxit=1000) # multinomial main effects / ANCOVA model
fit2 = multinom(variant_lineage~sample_date_num*nhs_name, data = data, maxit=10000) # multinomial homogeneity of slopes model
# or with spline terms in function of sample_date (df were tuned based on BIC value)
fit3 = multinom(variant_lineage~ns(sample_date_num, df=6)+nhs_name, data = data, maxit=1000) # with nat cubic splines ifo time with df=6
fit4 = multinom(variant_lineage~ns(sample_date_num, df=6)*nhs_name, data = data, maxit=10000) # with nat cubic splines ifo time with df=6
fit5 = multinom(variant_lineage~ns(sample_date_num, df=5)+nhs_name, data = data, maxit=1000) # with nat cubic splines ifo time with df=6
fit6 = multinom(variant_lineage~ns(sample_date_num, df=5)*nhs_name, data = data, maxit=10000) # with nat cubic splines ifo time with df=6
fit7 = multinom(variant_lineage~ns(sample_date_num, df=4)+nhs_name, data = data, maxit=1000) # with nat cubic splines ifo time with df=6
fit8 = multinom(variant_lineage~ns(sample_date_num, df=4)*nhs_name, data = data, maxit=10000) # with nat cubic splines ifo time with df=6
fit9 = multinom(variant_lineage~ns(sample_date_num, df=3)+nhs_name, data = data, maxit=1000) # with nat cubic splines ifo time with df=6
fit10 = multinom(variant_lineage~ns(sample_date_num, df=3)*nhs_name, data = data, maxit=10000) # with nat cubic splines ifo time with df=6
fit11 = multinom(variant_lineage~ns(sample_date_num, df=2)+nhs_name, data = data, maxit=1000) # with nat cubic splines ifo time with df=6
fit12 = multinom(variant_lineage~ns(sample_date_num, df=2)*nhs_name, data = data, maxit=10000) # with nat cubic splines ifo time with df=6
fit3B = multinom(variant_lineage~ns(sample_date_num, df=7)+nhs_name, data = data, maxit=1000) # with nat cubic splines ifo time with df=6
fit3C = multinom(variant_lineage~ns(sample_date_num, df=8)+nhs_name, data = data, maxit=1000) # with nat cubic splines ifo time with df=6
fit3D = multinom(variant_lineage~ns(sample_date_num, df=9)+nhs_name, data = data, maxit=1000) # with nat cubic splines ifo time with df=6
fit3E = multinom(variant_lineage~ns(sample_date_num, df=10)+nhs_name, data = data, maxit=1000) # with nat cubic splines ifo time with df=6
fit3F = multinom(variant_lineage~ns(sample_date_num, df=12)+nhs_name, data = data, maxit=1000) # with nat cubic splines ifo time with df=6
fit3G = multinom(variant_lineage~ns(sample_date_num, df=14)+nhs_name, data = data, maxit=1000) # with nat cubic splines ifo time with df=6

BIC(fit1, fit2, fit3, fit4, fit5, fit6, fit7, fit8, fit9, fit10, fit11, fit12,
    fit3B, fit3C, fit3D, fit3E, fit3F, fit3G)[order(BIC(fit1, fit2, fit3, fit4, fit5, fit6, 
                                                        fit7, fit8, fit9, fit10, fit11, fit12,fit3B, fit3C, fit3D, fit3E, fit3F, fit3G)$BIC),]
#        df      BIC
# fit10 324 236303.1
# fit8  405 236306.8
# fit6  486 236749.3
# fit12 243 236750.2
# fit3C 153 237050.1
# fit3  135 237053.4
# fit3B 144 237081.2
# fit3E 171 237123.3
# fit5  126 237126.6
# fit3D 162 237142.4
# fit3F 189 237173.3
# fit7  117 237175.5
# fit4  567 237269.0
# fit3G 207 237273.5
# fit9  108 237289.5
# fit11  99 237812.8
# fit1   90 240472.4
# fit2  162 333954.7

# we continue with fit12 as this was a realistic fit & one that also still extrapolated in a stable way
# and allows for differences in growth rates between regions
fit = fit12 # simplest spline model with date*region interaction that extrapolates well (with 2 spline knots)
# PS fit8, fit10 & fit6 all extrapolated in unstable way
# PS2 fit3C had best BIC of main effect spline models that extrapolated OK




# OUTPUT MULTINOMIAL MAIN EFFECTS / ANCOVA MODEL fit1
# summary(fit)


# mean growth rates of strains 20B.501Y.V1, 20A.EU1 & other 
# and pairwise differences in growth rates for the period between Sept 1 2020 and 17 Dec 2020
fit_emtrends = emtrends(fit, revpairwise ~ variant_lineage, var='sample_date_num', 
                         mode="latent", adjust="Tukey", 
                         at=list(sample_date_num=seq(as.numeric(as.Date("2020-09-01")),
                                                     as.numeric(max(data$sample_date))),
                                         variant_lineage=c("20B.501Y.V1","20A.EU1","other"))) 
# trends conditional on date to show temporal evolution of growth rates
fit_emtrends_date = emtrends(fit, revpairwise ~ variant_lineage|sample_date_num, var='sample_date_num', 
                              mode="latent", adjust="Tukey", 
                              at=list(sample_date_num=seq(as.numeric(min(data$sample_date)),
                                                          as.numeric(max(data$sample_date))),
                                      variant_lineage=c("20B.501Y.V1","20A.EU1","other"))) 
fit_emtrends_date = as.data.frame(fit_emtrends_date)
fit_emtrends_date$sample_date = as.Date(as.numeric(fit_emtrends_date$sample_date_num), origin="1970-01-01")
# plot of mean growth rates over time
plot(fit_emtrends_date$sample_date[fit_emtrends_date$variant_lineage=="20B.501Y.V1"],
     fit_emtrends_date$sample_date_num.trend[fit_emtrends_date$variant_lineage=="20B.501Y.V1"],
     type="l", col="red", lwd=2, ylab="Mean growth rate", xlab="Sample date", ylim=c(-0.5,0.5))
lines(fit_emtrends_date$sample_date[fit_emtrends_date$variant_lineage=="20A.EU1"],
      fit_emtrends_date$sample_date_num.trend[fit_emtrends_date$variant_lineage=="20A.EU1"], col="blue", lwd=2)
lines(fit_emtrends_date$sample_date[fit_emtrends_date$variant_lineage=="20A.EU1"],
      fit_emtrends_date$sample_date_num.trend[fit_emtrends_date$variant_lineage=="20A.EU1"], col="blue", lwd=2)
lines(fit_emtrends_date$sample_date[fit_emtrends_date$variant_lineage=="other"],
      fit_emtrends_date$sample_date_num.trend[fit_emtrends_date$variant_lineage=="other"], col="darkgrey", lwd=2)
# pairwise differences in growth rates of strains 20B.501Y.V1 and 20A.EU1 vs "other"
# CHECK IF THERE ARE FITNESS DIFFERENCES RELATED TO PERIOD OF LOCKDOWN 
# 5 NOV - 2 DEC England, Wales XX - XX, Scotland XX - XX
# not sure what to conclude from this plot...
plot(fit_emtrends_date$sample_date[fit_emtrends_date$variant_lineage=="20B.501Y.V1"],
     fit_emtrends_date$sample_date_num.trend[fit_emtrends_date$variant_lineage=="20B.501Y.V1"]-
     fit_emtrends_date$sample_date_num.trend[fit_emtrends_date$variant_lineage=="other"],
     type="l", col="red", lwd=2, ylab="Differences in growth rates with other competing strains", xlab="Sample date", ylim=c(0,0.5))
lines(fit_emtrends_date$sample_date[fit_emtrends_date$variant_lineage=="20A.EU1"],
      fit_emtrends_date$sample_date_num.trend[fit_emtrends_date$variant_lineage=="20A.EU1"]-
        fit_emtrends_date$sample_date_num.trend[fit_emtrends_date$variant_lineage=="other"], 
      col="blue", lwd=2)


# marginal mean growth rates during period from Sept 1 202 until Dec 17 2020
fit_trends = as.data.frame(fit_emtrends$emtrends)
colnames(fit_trends)[which(colnames(fit_trends)=="sample_date_num.trend")] = "r" # column sample_date_num.trend is marginal mean growth rate r
fit_trends
#   variant_lineage          r          SE  df   lower.CL   upper.CL
# 1         20A.EU1 0.06928386 0.003981986 243 0.05970994 0.07885778
# 2     20B.501Y.V1 0.20299107 0.009371493 243 0.18045911 0.22552303
# 3           other 0.04634709 0.003972450 243 0.03679609 0.05589808

# the contrasts (difference) in growth rates between the strains is 
fit_contrasts = data.frame(as.data.frame(fit_emtrends$contrasts),
                          as.data.frame(confint(fit_emtrends$contrasts))[,c("lower.CL","upper.CL")])
colnames(fit_contrasts)[which(colnames(fit_contrasts) %in% c("estimate","lower.CL","upper.CL"))] = 
  c("delta_r","delta_r.lower.CL","delta_r.upper.CL")
# which for our two sets of generation time estimates would result in changes in R
# delta_R1 and delta_R2 of
fit_contrasts = data.frame(fit_contrasts[,c("contrast","SE","t.ratio","p.value")],
  delta_R.from.delta_r_df(fit_contrasts[,c("delta_r","delta_r.lower.CL","delta_r.upper.CL")]))
fit_contrasts
# contrast           SE   t.ratio p.value     delta_r delta_r.lower.CL delta_r.upper.CL   delta_R1
# 1 20B.501Y.V1 - 20A.EU1 0.0095126256  14.05576       0  0.13370721       0.11127486       0.15613956  0.7353896
# 2       other - 20A.EU1 0.0004163371 -55.09184       0 -0.02293678      -0.02391857      -0.02195498 -0.1261523
# 3   other - 20B.501Y.V1 0.0095171128 -16.45919       0 -0.15664398      -0.17908692      -0.13420105 -0.8615419
# delta_R1.LCL delta_R1.UCL    delta_R2 delta_R2.LCL delta_R2.UCL
# 1    0.6120117    0.8587675  0.48131966   0.40056761   0.56207171
# 2   -0.1315521   -0.1207524 -0.08256789  -0.08610214  -0.07903363
# 3   -0.9849780   -0.7381058 -0.56388755  -0.64467769  -0.48309741


# plot model predictions
extrapolate = 60 # nr of days to show extrapolations for
newdat = expand.grid(sample_date_num =
                         seq(min(data$sample_date_num),max(data$sample_date_num)+extrapolate,1),
                     nhs_name = levels_nhs_name)
fit_preds = data.frame(newdat, predict(fit, type="prob", newdata=newdat), check.names=F)
fit_preds = gather(fit_preds, variant_lineage, prob, levels_variant_lineages)
fit_preds$sample_date = as.Date(fit_preds$sample_date_num, origin="1970-01-01")
fit_preds$variant_lineage = factor(fit_preds$variant_lineage, levels=levels_variant_lineages)
fit_preds$nhs_name = factor(fit_preds$nhs_name, levels=levels_nhs_name)
muller_fit = ggplot(data=fit_preds, 
                     aes(x=sample_date, y=prob, group=variant_lineage)) + 
    facet_wrap(~nhs_name) +
    geom_area(aes(lwd=I(1.2), colour=NULL, fill=variant_lineage), position="stack") +
    annotate("rect", xmin=max(data$sample_date)+1, xmax=as.Date("2021-01-31"), ymin=0, ymax=1, alpha=0.4, fill="white") + # extrapolated part
    scale_fill_manual("", values=lineage_cols) +
    scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
                       labels=substring(months(as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01"))),1,1),
                       limits=as.Date(c("2020-03-01","2021-01-31")), expand=c(0,0)) +
    guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
    theme_hc() + theme(legend.position="bottom", 
                       axis.title.x=element_blank()) + 
    labs(title = "MAIN SARS-CoV2 VARIANT LINEAGES IN THE UK") +
    ylab("Relative abundance")
muller_fit
saveRDS(muller_fit, file = ".\\multinomial_fits\\fit12_multinomial fit_muller plot fit_FINAL.rds")
graph2ppt(file=".\\multinomial_fits\\fit12_multinomial fit_muller plot fit_FINAL.pptx", width=7, height=5)
ggsave(".\\multinomial_fits\\fit12_multinomial fit_muller plot fit_FINAL.png", width=7, height=5)
# plot with raw data and fit combined:
ggarrange(muller_raw +
              coord_cartesian(xlim=c(as.Date("2020-03-01"),
                                     as.Date("2021-01-31")), 
                              expand=c(0,0)) + labs(title = "") , 
          muller_fit+labs(title = "")+coord_cartesian(xlim=c(as.Date("2020-03-01"),
                                                             as.Date("2021-01-31")), 
                                                      expand=c(0,0)),
          ncol=1, common.legend=TRUE, legend="bottom")




# plot predictions specifically for strains 20B.501Y.V1 and 20A.EU1

extrapolate = 60 # nr of days to show extrapolations for
# we calculate predictions & 95% confidence intervals using the effects package here (emmeans ran out of memory)
fit_strainpreds.eff = Effect(c("sample_date_num","nhs_name"), 
                              fit, xlevels=list(sample_date_num=seq(min(data$sample_date_num),
                                                                     max(data$sample_date_num)+extrapolate)))
fit_strainpreds = data.frame(fit_strainpreds.eff,check.names=F)
cols = names(fit_strainpreds.eff$x)
fit_strainpreds_probs = fit_strainpreds[,grepl("^prob\\.", colnames(fit_strainpreds))]
colnames(fit_strainpreds_probs) = gsub("prob\\.|prob\\.X","", colnames(fit_strainpreds_probs))
fit_strainpreds_probs = gather(data.frame(fit_strainpreds[,cols], fit_strainpreds_probs, check.names=F), 
                                variant_lineage, prob, levels(data$variant_lineage), factor_key=TRUE)
fit_strainpreds_LCL = fit_strainpreds[,grepl("^L\\.prob\\.", colnames(fit_strainpreds))]
colnames(fit_strainpreds_LCL) = gsub("L\\.prob\\.|L\\.prob\\.X","", colnames(fit_strainpreds_LCL))
fit_strainpreds_LCL = gather(fit_strainpreds_LCL, variant_lineage, LCL, levels(data$variant_lineage), factor_key=TRUE)[,-1,drop=F]
fit_strainpreds_UCL = fit_strainpreds[,grepl("^U\\.prob\\.", colnames(fit_strainpreds))]
colnames(fit_strainpreds_UCL) = gsub("U\\.prob\\.|U\\.prob\\.X","", colnames(fit_strainpreds_UCL))
fit_strainpreds_UCL = gather(fit_strainpreds_UCL, variant_lineage, UCL, levels(data$variant_lineage), factor_key=TRUE)[,-1,drop=F]
fit_strainpreds = data.frame(fit_strainpreds_probs, fit_strainpreds_LCL, fit_strainpreds_UCL, check.names=F)
fit_strainpreds$sample_date = as.Date(fit_strainpreds$sample_date_num, origin="1970-01-01")
fit_strainpreds$nhs_name = factor(fit_strainpreds$nhs_name, levels=levels_nhs_name)
fit_strainpreds$variant_lineage = factor(fit_strainpreds$variant_lineage, levels=levels_variant_lineages)

# It can be seen that the fit for 20B.501Y.V1 corresponds perfectly to a logistical model
# with different dates of introduction for the different regions, but identical rates of spread.
# For 20A.EU1 clonal interference by 20B.501Y.V1 results in the relative abundance of that
# strain starting to decline after an initial prolonged spread from July onwards.
# This is nicely captured by the multinomial fit.
# The raw data plotted on this graph have been aggregated by week.
plotmultinom2strains = qplot(data=fit_strainpreds[fit_strainpreds$variant_lineage %in% c("20A.EU1",
                                                                   "20B.501Y.V1"),], 
      x=sample_date, y=prob, colour=NULL, fill=nhs_name, group=nhs_name, alpha=NULL, geom="blank") +
  facet_wrap(~variant_lineage, ncol=1) +
  geom_ribbon(aes(y=prob, ymin=LCL, ymax=UCL, colour=NULL, fill=nhs_name), alpha=I(0.5)) +
  geom_line(aes(y=prob, colour=nhs_name), alpha=I(0.8)) +
  geom_point(data=data_agbyweekregion[data_agbyweekregion$variant_lineage %in% c("20A.EU1",
                                                                                 "20B.501Y.V1"),], 
             aes(x=sample_date, y=prop, colour=nhs_name, fill=NULL, size=total, alpha=I(0.5))) +
  # labs(tag = "@TWenseleers\ndata COG-UK") +
  # theme(plot.tag.position = "bottomright",
  #       plot.tag = element_text(vjust = 1, hjust = 1, size=8)) +
  ylab("Relative abundance (%)") +
  theme_hc() + xlab("") + 
  # ggtitle("GROWTH OF STRAIN 20A.EU1 AND VOC 20A.EU1") +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01")),
                     labels=c("M","A","M","J","J","A","S","O","N","D","J","F")) +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(# xlim=c(as.Date("2020-09-01"),as.Date("2021-02-01")), 
    xlim=c(as.Date("2020-07-01"),as.Date("2021-02-01")), 
    ylim=c(0.0001,0.99900001), expand=c(0,0)) +
  scale_color_discrete("region", h=c(0, 280), c=150) +
  scale_fill_discrete("region", h=c(0, 280), c=150) +
  scale_size_continuous("total number\nof sequences\nper week", trans="log2", 
                        range=c(0.01, 4), limits=c(100,1600), breaks=c(100,200,400,800,1600)) +
  theme(legend.direction = "vertical", legend.box = "vertical", legend.position="right", 
        legend.key.size=unit(0.45, "cm"))
plotmultinom2strains

saveRDS(plotmultinom2strains, file = ".\\multinomial_fits\\fit12_multinomial fit_growth VOC strain 20B50YV1 and 20AEU1_FINAL.rds")
graph2ppt(file=".\\multinomial_fits\\fit12_multinomial fit_growth VOC strain 20B50YV1 and 20AEU1_FINAL.pptx", width=6, height=6)
ggsave(file=".\\multinomial_fits\\fit12_multinomial fit_growth VOC strain 20B50YV1 and 20AEU1_FINAL.png", width=6, height=6)







# 3. BINOMIAL GLMM FITS TO COMPARE RATES OF SPREAD OF STRAINS 20B.501Y.V1 AND 20A.EU1 ACROSS REGIONS IN UK ####

# 3.1 BINOMIAL MIXED MODEL FIT FOR STRAIN 20B.501Y.V1 ####

# mixed binomial GLMM with nested random effects nhs_name, lad and observation 
# (the latter to take into account overdispersion)
# and with or without fixed effect interaction nhs_name x date

bGLMMfit2 = glmer(cbind(count, total-count) ~  (1|nhs_name/lad) + 
                    nhs_name+scale(sample_date_num), 
                  family=binomial(logit), data=data_agbydayregionlad,
                  subset=data_agbydayregionlad$variant_lineage=="20B.501Y.V1"&
                    data_agbydayregionlad$sample_date>="2020-08-01") 
bGLMMfit2_od = glmer(cbind(count, total-count) ~  (1|nhs_name/lad/obs) + 
                       nhs_name+scale(sample_date_num), 
                     family=binomial(logit), data=data_agbydayregionlad,
                     subset=data_agbydayregionlad$variant_lineage=="20B.501Y.V1"&
                       data_agbydayregionlad$sample_date>="2020-08-01")
bGLMMfit3 = glmer(cbind(count, total-count) ~  (1|nhs_name/lad) + 
                    nhs_name*scale(sample_date_num), 
                  family=binomial(logit), data=data_agbydayregionlad,
                  subset=data_agbydayregionlad$variant_lineage=="20B.501Y.V1"&
                    data_agbydayregionlad$sample_date>="2020-08-01")
bGLMMfit3_od = glmer(cbind(count, total-count) ~  (1|nhs_name/lad/obs) + 
                       nhs_name*scale(sample_date_num), 
                     family=binomial(logit), data=data_agbydayregionlad,
                     subset=data_agbydayregionlad$variant_lineage=="20B.501Y.V1"&
                       data_agbydayregionlad$sample_date>="2020-08-01")
bGLMMfit4 = glmer(cbind(count, total-count) ~  (1|lad) + # 
                    nhs_name*scale(sample_date_num), 
                  family=binomial(logit), data=data_agbydayregionlad,
                  subset=data_agbydayregionlad$variant_lineage=="20B.501Y.V1"&
                    data_agbydayregionlad$sample_date>="2020-08-01")
bGLMMfit4_od = glmer(cbind(count, total-count) ~  (1|lad/obs) + 
                       nhs_name*scale(sample_date_num), 
                     family=binomial(logit), data=data_agbydayregionlad,
                     subset=data_agbydayregionlad$variant_lineage=="20B.501Y.V1"&
                       data_agbydayregionlad$sample_date>="2020-08-01")

# check AIC & BIC values
BIC(bGLMMfit2, bGLMMfit2_od, bGLMMfit3, bGLMMfit3_od, bGLMMfit4, bGLMMfit4_od)
#              df      BIC
# bGLMMfit2    12 3282.132
# bGLMMfit2_od 13 3128.136
# bGLMMfit3    20 3304.187
# bGLMMfit3_od 21 3178.216
# bGLMMfit4    19 3295.137
# bGLMMfit4_od 20 3169.130

AIC(bGLMMfit2, bGLMMfit2_od, bGLMMfit3, bGLMMfit3_od, bGLMMfit4, bGLMMfit4_od)
#              df      AIC
# bGLMMfit2    12 3197.801
# bGLMMfit2_od 13 3036.777
# bGLMMfit3    20 3163.635
# bGLMMfit3_od 21 3030.636
# bGLMMfit4    19 3161.613
# bGLMMfit4_od 20 3028.577

# bGLMMfit2_od has the best BIC (homogenous slopes), but bGLMMfit4_od has the best AIC (heterogeneous slopes)
# since below we specifically would like to test for heterogeneity of slopes across regions we
# will continue with bGLMMfit4_od

# we continue with fit bGLMMfit4_od as that one allows for different rates of spread
# in different regions
fit = bGLMMfit4_od

summary(fit)

# plain effect plot with partial residuals
plot(Effect(c("sample_date_num","nhs_name"), 
            bGLMMfit4_od, xlevels=list(sample_date_num=seq(as.numeric(as.Date("2020-10-01")),
                                                           as.numeric(as.Date("2020-12-31")))),
            x.var="sample_date_num", residuals=TRUE),
     partial.residuals=TRUE, use.spline=FALSE, residuals.pch=16, 
     residuals.color=alpha("steelblue",0.2), ylab="VOC 20B.501Y.V1 (proportion)",
     smooth.residuals=TRUE, span=0.2)
# graph2png(file=".\\multinomial_fits\\fit binomial mixed model growth VOC strain 20B_501Y_V1_effect plot with partial residuals.png", width=8, height=6)


# PLOT MODEL FIT

extrapolate = 60 # 60 nr of days to extrapolate fit into the future
bGLMM_preds = as.data.frame(emmeans(fit, ~ sample_date_num|nhs_name, at=list(sample_date_num=
                                                                    seq(as.numeric(as.Date("2020-09-01")),
                                                                        max(data$sample_date_num[data$sample_date>="2020-09-01"])+extrapolate)), 
                           type="response"))
bGLMM_preds$sample_date = as.Date(bGLMM_preds$sample_date_num, origin="1970-01-01")
bGLMM_preds$nhs_name = factor(bGLMM_preds$nhs_name, 
                              levels=levels_nhs_name)
plot_bGLMMVOC = qplot(data=bGLMM_preds, x=sample_date, y=prob, geom="blank") +
    facet_wrap(~nhs_name) +
    geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL, fill=nhs_name), alpha=I(0.3)) +
    # geom_ribbon(aes(y=prob, ymin=lower.CL, ymax=upper.CL, colour=NULL, fill=nhs_name), alpha=I(0.3)) +
    geom_line(aes(y=prob, colour=nhs_name), alpha=I(0.8)) +
    # labs(tag = "@TWenseleers\ndata COG-UK") +
    # theme(plot.tag.position = "bottomright",
    #       plot.tag = element_text(vjust = 1, hjust = 1, size=8)) +
    # ylab("Relative abundance of VOC strain 202012/01 (%)") +
    ylab("Relative abundance (%)") +
    theme_hc() + xlab("") + 
    ggtitle("GROWTH OF VOC 20B.501Y.V1 BY NHS REGION") +
    scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01")),
                       labels=c("M","A","M","J","J","A","S","O","N","D","J","F")) +
    scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                        labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(# xlim=c(as.Date("2020-09-01"),as.Date("2021-02-01")), 
       xlim=c(as.Date("2020-07-01"),as.Date("2021-01-31")), 
       ylim=c(0.0001,0.999001), expand=c(0,0)) +
   scale_color_discrete("", h=c(0, 280), c=200) +
    scale_fill_discrete("", h=c(0, 280), c=200) +
    geom_point(data=data_agbyweekregion[data_agbyweekregion$variant_lineage=="20B.501Y.V1",], 
                aes(x=sample_date, y=prop, colour=nhs_name, size=total), alpha=I(0.5)) +
    scale_size_continuous("total number\nof sequences\nper week", trans="log2", 
                        range=c(0.01, 4), limits=c(100,1600), breaks=c(100,200,400,800,1600)) +
    guides(fill=FALSE) + guides(colour=FALSE) + theme(legend.position = "right")
    # geom_point(data=data_agbydayregion[data_agbydayregion$variant_lineage=="20B.501Y.V1",], 
    #         aes(x=sample_date, y=prop, colour=nhs_name, size=total), alpha=I(0.4)) +
    # theme(legend.position = "none")
plot_bGLMMVOC
saveRDS(plot_bGLMMVOC, file = ".\\multinomial_fits\\fit binomial mixed model growth VOC strain 20B_501Y_V1.rds")
graph2ppt(file=".\\multinomial_fits\\fit binomial mixed model growth VOC strain 20B_501Y_V1.pptx", width=8, height=6)
ggsave(file=".\\multinomial_fits\\fit binomial mixed model growth VOC strain 20B_501Y_V1.png", width=8, height=6)


#  CALCULATE GROWTH RATES & SELECTIVE ADVANTAGE OF VOC STRAIN IN DIFFERENT REGIONS

# growth rate of strain = slope of logistic regression = selective advantage s / generation time
# where selective advantage s is the prop increase in the R value compared to other circulating strains
# (assuming identical generation times)

bGLMM_VOC_growthrates = as.data.frame(emtrends(fit, ~ nhs_name, var="sample_date_num"))[,-c(3,4)] 
# sample_date_num.trend = 
# logistic growth rate of VOC 202012/01 strain = growth rate VOC strain - growth rate all other strains
colnames(bGLMM_VOC_growthrates)[2] = "logistic_growth_rate"

# for our two sets of default estimated growth rates these logistic growth rates 
# translate to a proportional increase in R value
# compared to all other circulating strains of 
# (assuming that generation times would be the same)
bGLMM_VOC_growthrates = delta_R.from.delta_r_df(bGLMM_VOC_growthrates)
bGLMM_VOC_growthrates
#                   nhs_name logistic_growth_rate  asymp.LCL asymp.UCL  delta_R1 delta_R1.LCL delta_R1.UCL  delta_R2
# 1               South East           0.10279920 0.09275695 0.1128414 0.5653956    0.5101632    0.6206279 0.3700569
# 2                   London           0.11074404 0.09758339 0.1239047 0.6090922    0.5367086    0.6814757 0.3986568
# 3          East of England           0.09328065 0.08276698 0.1037943 0.5130436    0.4552184    0.5708688 0.3357920
# 4               South West           0.11387790 0.08175591 0.1459999 0.6263284    0.4496575    0.8029994 0.4099381
# 5                 Midlands           0.11808105 0.09867610 0.1374860 0.6494458    0.5427185    0.7561730 0.4250686
# 6 North East and Yorkshire           0.15286560 0.12347369 0.1822575 0.8407608    0.6791053    1.0024163 0.5502861
# 7                 Scotland           0.11126736 0.07349657 0.1490382 0.6119705    0.4042311    0.8197099 0.4005406
# 8               North West           0.13183790 0.10057672 0.1630991 0.7251084    0.5531719    0.8970449 0.4745905
# 9                    Wales           0.10330451 0.07238902 0.1342200 0.5681748    0.3981396    0.7382099 0.3718759
#   delta_R2.LCL delta_R2.UCL
# 1    0.3339068    0.4062070
# 2    0.3512810    0.4460325
# 3    0.2979449    0.3736392
# 4    0.2943052    0.5255709
# 5    0.3552146    0.4949226
# 6    0.4444810    0.6560912
# 7    0.2645732    0.5365081
# 8    0.3620564    0.5871246
# 9    0.2605863    0.4831656

# on average across all regions we get
bGLMM_VOC_growthrates_avg = as.data.frame(emtrends(fit, ~ 1, var="sample_date_num"))[,-c(3,4)] 
colnames(bGLMM_VOC_growthrates_avg)[2] = "logistic_growth_rate"

# for our two sets of default estimated growth rates these logistic growth rates 
# translate into a proportional increase in R value
# compared to all other circulating strains of 
# (assuming that generation times would be the same)
bGLMM_VOC_growthrates_avg = delta_R.from.delta_r_df(bGLMM_VOC_growthrates_avg)
bGLMM_VOC_growthrates_avg
# 1         logistic_growth_rate asymp.LCL asymp.UCL  delta_R1 delta_R1.LCL delta_R1.UCL  delta_R2 delta_R2.LCL
# 1 overall            0.1153398 0.1066392 0.1240404 0.6343689    0.5865158     0.682222 0.4152006    0.3838803
# delta_R2.UCL
# 1    0.4465209


# TUKEY POSTHOC TESTS TO TEST FOR DIFFERENCES IN GROWTH RATES IN THE VOC STRAIN ACROSS DIFFERENT NHS REGIONS

tukey_VOC = as.data.frame(emtrends(fit, pairwise ~ nhs_name, var="sample_date_num", adjust="tukey")$contrasts)[,-4]
colnames(tukey_VOC)[2] = "diff_logistic_growth_rate"
tukey_VOC
#                                      contrast diff_logistic_growth_rate          SE     z.ratio     p.value
# 1                         South East - London             -0.0079448409 0.008370080 -0.94919536 0.990050138
# 2                South East - East of England              0.0095185440 0.007287041  1.30622901 0.929859366
# 3                     South East - South West             -0.0110787073 0.017117019 -0.64723346 0.999322699
# 4                       South East - Midlands             -0.0152818546 0.011123865 -1.37379003 0.907861247
# 5       South East - North East and Yorkshire             -0.0500664067 0.015823462 -3.16406149 0.041428077
# 6                       South East - Scotland             -0.0084681686 0.019907061 -0.42538518 0.999971482
# 7                     South East - North West             -0.0290387048 0.016736255 -1.73507778 0.724789811
# 8                          South East - Wales             -0.0005053108 0.016569780 -0.03049593 1.000000000
# 9                    London - East of England              0.0174633849 0.008507569  2.05268810 0.506193123
# 10                        London - South West             -0.0031338664 0.017672156 -0.17733356 0.999999971
# 11                          London - Midlands             -0.0073370137 0.011949768 -0.61398798 0.999540772
# 12          London - North East and Yorkshire             -0.0421215658 0.016421779 -2.56498187 0.201489416
# 13                          London - Scotland             -0.0005233277 0.020392811 -0.02566236 1.000000000
# 14                        London - North West             -0.0210938639 0.017293637 -1.21974709 0.952423812
# 15                             London - Wales              0.0074395301 0.017128478  0.43433691 0.999966520
# 16               East of England - South West             -0.0205972514 0.017184142 -1.19861974 0.957039284
# 17                 East of England - Midlands             -0.0248003986 0.011239010 -2.20663553 0.401144746
# 18 East of England - North East and Yorkshire             -0.0595849508 0.015911595 -3.74475037 0.005652289
# 19                 East of England - Scotland             -0.0179867126 0.019980466 -0.90021488 0.993017193
# 20               East of England - North West             -0.0385572489 0.016809277 -2.29380763 0.345706518
# 21                    East of England - Wales             -0.0100238549 0.016637481 -0.60248633 0.999600828
# 22                      South West - Midlands             -0.0042031473 0.019133650 -0.21967305 0.999999839
# 23      South West - North East and Yorkshire             -0.0389876994 0.022198178 -1.75634684 0.711105381
# 24                      South West - Scotland              0.0026105388 0.025274634  0.10328691 1.000000000
# 25                    South West - North West             -0.0179599975 0.022856595 -0.78576872 0.997266278
# 26                         South West - Wales              0.0105733965 0.022728600  0.46520228 0.999943296
# 27        Midlands - North East and Yorkshire             -0.0347845522 0.017951476 -1.93769873 0.587294857
# 28                        Midlands - Scotland              0.0068136860 0.021643656  0.31481215 0.999997259
# 29                      Midlands - North West             -0.0137568503 0.018769726 -0.73292762 0.998331893
# 30                           Midlands - Wales              0.0147765438 0.018614004  0.79384017 0.997062817
# 31        North East and Yorkshire - Scotland              0.0415982382 0.024365420  1.70726536 0.742309422
# 32      North East and Yorkshire - North West              0.0210277019 0.021890267  0.96059594 0.989234213
# 33           North East and Yorkshire - Wales              0.0495610959 0.021746915  2.27899431 0.354861443
# 34                      Scotland - North West             -0.0205705363 0.025010714 -0.82246897 0.996238694
# 35                           Scotland - Wales              0.0079628577 0.024881848  0.32002679 0.999996882
# 36                         North West - Wales              0.0285333940 0.022427434  1.27225404 0.939449723

tukey_VOC[tukey_VOC$p.value<0.05,]
#                                      contrast diff_logistic_growth_rate         SE   z.ratio     p.value
# 5       South East - North East and Yorkshire               -0.05006641 0.01582346 -3.164061 0.041428077
# 18 East of England - North East and Yorkshire               -0.05958495 0.01591159 -3.744750 0.00565228

# Almost no sign differences in rate of spread of VOC, only 2 small differences 
# (slightly faster spread in North East & Yorkshire vs South East & East of England)



# 3.2 BINOMIAL MIXED MODEL FIT FOR STRAIN 20A.EU1 ####

# equivalent binomial mixed model for earlier strain 20A.EU1

# mixed binomial GLMM with nested random effects nhs_name, lad and observation (the latter to take into account overdispersion)
bGLMMfit4_odB = glmer(cbind(count, total-count) ~  (1|lad/obs) + 
                       nhs_name*scale(sample_date_num), 
                     family=binomial(logit), data=data_agbydayregionlad,
                     subset=data_agbydayregionlad$variant_lineage=="20A.EU1"&
                       data_agbydayregionlad$sample_date>="2020-07-01"&
                       data_agbydayregionlad$sample_date<="2020-09-30")
BIC(bGLMMfit4_odB) # 7689.543

fit = bGLMMfit4_odB

summary(fit)


# PLOT MODEL FIT

extrapolate = 60 # 60 nr of days to extrapolate fit into the future
bGLMM_predsB = as.data.frame(emmeans(fit, ~ sample_date_num|nhs_name, at=list(sample_date_num=
                                                                               seq(as.numeric(as.Date("2020-06-01")),
                                                                                   max(data$sample_date_num[data$sample_date>="2020-09-30"])+extrapolate)), 
                                    type="response"))
bGLMM_predsB$sample_date = as.Date(bGLMM_predsB$sample_date_num, origin="1970-01-01")
bGLMM_predsB$nhs_name = factor(bGLMM_predsB$nhs_name, 
                              levels=levels_nhs_name)
plot_bGLMM20AEU1 = qplot(data=bGLMM_predsB, x=sample_date, y=prob, geom="blank") +
  facet_wrap(~nhs_name) +
  geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL, fill=nhs_name), alpha=I(0.3)) +
  # geom_ribbon(aes(y=prob, ymin=lower.CL, ymax=upper.CL, colour=NULL, fill=nhs_name), alpha=I(0.3)) +
  geom_line(aes(y=prob, colour=nhs_name), alpha=I(0.8)) +
  # labs(tag = "@TWenseleers\ndata COG-UK") +
  # theme(plot.tag.position = "bottomright",
  #       plot.tag = element_text(vjust = 1, hjust = 1, size=8)) +
  # ylab("Relative abundance of VOC strain 202012/01 (%)") +
  ylab("Relative abundance (%)") +
  theme_hc() + xlab("") + 
  ggtitle("GROWTH OF STRAIN 20A.EU1 BY NHS REGION") +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01")),
                     labels=c("M","A","M","J","J","A","S","O","N","D","J","F")) +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(# xlim=c(as.Date("2020-09-01"),as.Date("2021-02-01")), 
    xlim=c(as.Date("2020-06-01"),as.Date("2020-11-30")), 
    ylim=c(0.0001,0.999001), expand=c(0,0)) +
  scale_color_discrete("", h=c(0, 280), c=200) +
  scale_fill_discrete("", h=c(0, 280), c=200) +
  geom_point(data=data_agbyweekregion[data_agbyweekregion$variant_lineage=="20A.EU1"&
                                        data_agbyweekregion$sample_date>="2020-06-01"&
                                        data_agbyweekregion$sample_date<="2020-09-30",], 
             aes(x=sample_date, y=prop, colour=nhs_name, size=total), alpha=I(0.5)) +
  scale_size_continuous("total number\nof sequences\nper week", trans="log2", 
                        range=c(0.01, 4), limits=c(100,1600), breaks=c(100,200,400,800,1600)) +
  guides(fill=FALSE) + guides(colour=FALSE) + theme(legend.position = "right")
  # geom_point(data=data_agbydayregion[data_agbydayregion$variant_lineage=="20B.501Y.V1",], 
  #         aes(x=sample_date, y=prop, colour=nhs_name, size=total), alpha=I(0.4)) +
  # theme(legend.position = "none")
saveRDS(plot_bGLMM20AEU1, file = ".\\multinomial_fits\\fit binomial mixed model growth strain 20A_EU1.rds")
graph2ppt(file=".\\multinomial_fits\\fit binomial mixed model growth strain 20A_EU1.pptx", width=8, height=6)
ggsave(file=".\\multinomial_fits\\fit binomial mixed model growth strain 20A_EU1.png", width=8, height=6)


#  CALCULATE GROWTH RATES & SELECTIVE ADVANTAGE OF VOC STRAIN IN DIFFERENT REGIONS

# growth rate of strain = slope of logistic regression = selective advantage s / generation time
# where selective advantage s is the prop increase in the R value compared to other circulating strains
# (assuming identical generation times)

bGLMM_20AEU1_growthrates = as.data.frame(emtrends(fit, ~ nhs_name, var="sample_date_num"))[,-c(3,4)] 
colnames(bGLMM_20AEU1_growthrates)[2] = "logistic_growth_rate"

# for our two sets of default estimated growth rates these logistic growth rates 
# translate to a proportional increase in R value
# compared to all other circulating strains of 
# (assuming that generation times would be the same)
bGLMM_20AEU1_growthrates = delta_R.from.delta_r_df(bGLMM_20AEU1_growthrates)
bGLMM_20AEU1_growthrates
# nhs_name logistic_growth_rate  asymp.LCL  asymp.UCL  delta_R1 delta_R1.LCL delta_R1.UCL  delta_R2
# 1               South East           0.05691101 0.04619749 0.06762452 0.3130105    0.2540862    0.3719349 0.2048684
# 2                   London           0.03455288 0.02514782 0.04395795 0.1900409    0.1383130    0.2417687 0.1243836
# 3          East of England           0.04350878 0.03312170 0.05389587 0.2392983    0.1821693    0.2964273 0.1566231
# 4               South West           0.04816305 0.03456757 0.06175853 0.2648968    0.1901216    0.3396719 0.1733775
# 5                 Midlands           0.06887487 0.05855373 0.07919601 0.3788118    0.3220455    0.4355780 0.2479360
# 6 North East and Yorkshire           0.05888714 0.05023475 0.06753952 0.3238792    0.2762911    0.3714673 0.2119821
# 7                 Scotland           0.06168384 0.05413392 0.06923375 0.3392611    0.2977366    0.3807856 0.2220497
# 8               North West           0.08083494 0.07073663 0.09093326 0.4445922    0.3890514    0.5001329 0.2909899
# 9                    Wales           0.07098382 0.06084757 0.08112006 0.3904110    0.3346616    0.4461603 0.2555278
# delta_R2.LCL delta_R2.UCL
# 1   0.16630187    0.2434350
# 2   0.09052722    0.1582400
# 3   0.11923160    0.1940145
# 4   0.12443647    0.2223186
# 5   0.21078191    0.2850901
# 6   0.18083524    0.2431290
# 7   0.19487149    0.2492279
# 8   0.25463795    0.3273418
# 9   0.21903930    0.2920163

# on average across all regions we get
bGLMM_20AEU1_growthrates_avg = as.data.frame(emtrends(fit, ~ 1, var="sample_date_num"))[,-c(3,4)] 
colnames(bGLMM_20AEU1_growthrates_avg)[2] = "logistic_growth_rate"

# for our two sets of default estimated growth rates these logistic growth rates 
# translate into a proportional increase in R value
# compared to all other circulating strains of 
# (assuming that generation times would be the same)
bGLMM_20AEU1_growthrates_avg = delta_R.from.delta_r_df(bGLMM_20AEU1_growthrates_avg)
bGLMM_20AEU1_growthrates_avg
# 1         logistic_growth_rate  asymp.LCL  asymp.UCL  delta_R1 delta_R1.LCL delta_R1.UCL  delta_R2 delta_R2.LCL
# 1 overall            0.0582667 0.05467437 0.06185904 0.3204669     0.300709    0.3402247 0.2097487     0.196817
# delta_R2.UCL
# 1    0.2226804


# TUKEY POSTHOC TESTS TO TEST FOR DIFFERENCES IN GROWTH RATES IN THE VOC STRAIN ACROSS DIFFERENT NHS REGIONS

tukey_20AEU1 = as.data.frame(emtrends(fit, pairwise ~ nhs_name, var="sample_date_num", adjust="tukey")$contrasts)[,-4]
colnames(tukey_20AEU1)[2] = "diff_logistic_growth_rate"
tukey_20AEU1
#                                      contrast diff_logistic_growth_rate          SE    z.ratio      p.value
# 1                         South East - London               0.022358121 0.007238036  3.0889762 5.192807e-02
# 2                South East - East of England               0.013402222 0.007582098  1.7676139 7.037628e-01
# 3                     South East - South West               0.008747953 0.008799149  0.9941817 9.865215e-01
# 4                       South East - Midlands              -0.011963862 0.007567615 -1.5809289 8.155317e-01
# 5       South East - North East and Yorkshire              -0.001976130 0.006984826 -0.2829176 9.999988e-01
# 6                       South East - Scotland              -0.004772833 0.006613123 -0.7217215 9.985063e-01
# 7                     South East - North West              -0.023923936 0.007464117 -3.2051930 3.648920e-02
# 8                          South East - Wales              -0.014072812 0.007490843 -1.8786685 6.286171e-01
# 9                    London - East of England              -0.008955899 0.007107364 -1.2600873 9.426522e-01
# 10                        London - South West              -0.013610168 0.008391966 -1.6218092 7.930856e-01
# 11                          London - Midlands              -0.034321983 0.007094580 -4.8377754 4.608480e-05
# 12          London - North East and Yorkshire              -0.024334251 0.006464879 -3.7640688 5.252370e-03
# 13                          London - Scotland              -0.027130954 0.006054266 -4.4812956 2.550603e-04
# 14                        London - North West              -0.046282057 0.006977097 -6.6334257 1.180014e-09
# 15                             London - Wales              -0.036430933 0.007009064 -5.1976889 7.168168e-06
# 16               East of England - South West              -0.004654269 0.008691014 -0.5355266 9.998347e-01
# 17                 East of England - Midlands              -0.025366084 0.007444675 -3.4072791 1.893946e-02
# 18 East of England - North East and Yorkshire              -0.015378352 0.006848936 -2.2453638 3.760630e-01
# 19                 East of England - Scotland              -0.018175055 0.006466088 -2.8108271 1.120541e-01
# 20               East of England - North West              -0.037326159 0.007335013 -5.0887649 1.276551e-05
# 21                    East of England - Wales              -0.027475034 0.007363754 -3.7311181 5.951107e-03
# 22                      South West - Midlands              -0.020711815 0.008681852 -2.3856447 2.916992e-01
# 23      South West - North East and Yorkshire              -0.010724083 0.008173291 -1.3120887 9.281073e-01
# 24                      South West - Scotland              -0.013520786 0.007849238 -1.7225603 7.327295e-01
# 25                    South West - North West              -0.032671889 0.008582856 -3.8066454 4.461116e-03
# 26                         South West - Wales              -0.022820765 0.008610315 -2.6503983 1.660201e-01
# 27        Midlands - North East and Yorkshire               0.009987732 0.006836871  1.4608629 8.735913e-01
# 28                        Midlands - Scotland               0.007191029 0.006462714  1.1126949 9.725284e-01
# 29                      Midlands - North West              -0.011960075 0.007327205 -1.6322834 7.871326e-01
# 30                           Midlands - Wales              -0.002108950 0.007352134 -0.2868487 9.999987e-01
# 31        North East and Yorkshire - Scotland              -0.002796703 0.005740661 -0.4871744 9.999193e-01
# 32      North East and Yorkshire - North West              -0.021947806 0.006710807 -3.2705168 2.969046e-02
# 33           North East and Yorkshire - Wales              -0.012096682 0.006746988 -1.7929011 6.870666e-01
# 34                      Scotland - North West              -0.019151103 0.006301187 -3.0392852 6.004907e-02
# 35                           Scotland - Wales              -0.009299979 0.006356581 -1.4630472 8.726455e-01
# 36                         North West - Wales               0.009851125 0.007238461  1.3609419 9.123505e-01

tukey_20AEU1[tukey_20AEU1$p.value<0.05,]
#                                 contrast diff_logistic_growth_rate          SE   z.ratio      p.value
# 7                South East - North West               -0.02392394 0.007464117 -3.205193 3.648920e-02
# 11                     London - Midlands               -0.03432198 0.007094580 -4.837775 4.608480e-05
# 12     London - North East and Yorkshire               -0.02433425 0.006464879 -3.764069 5.252370e-03
# 13                     London - Scotland               -0.02713095 0.006054266 -4.481296 2.550603e-04
# 14                   London - North West               -0.04628206 0.006977097 -6.633426 1.180014e-09
# 15                        London - Wales               -0.03643093 0.007009064 -5.197689 7.168168e-06
# 17            East of England - Midlands               -0.02536608 0.007444675 -3.407279 1.893946e-02
# 20          East of England - North West               -0.03732616 0.007335013 -5.088765 1.276551e-05
# 21               East of England - Wales               -0.02747503 0.007363754 -3.731118 5.951107e-03
# 25               South West - North West               -0.03267189 0.008582856 -3.806645 4.461116e-03
# 32 North East and Yorkshire - North West               -0.02194781 0.006710807 -3.270517 2.969046e-02

# Unlike the VOC which displaces the other strains at a consistently and very high rate in all regions,
# there were many more differences in logistic growth rates across regions here, implying it was
# not displacing the other strains at a constant, uniform rate.




# 4. COMPARISON OF RATE OF SPREAD OF STRAIN 20B.501Y.V1 IN UK & DENMARK ####

# data from Denmark aggregated by week are provided by the Statens Serum Institut, link
# https://www.ssi.dk/-/media/cdn/files/scenarier_for_udviklingen_i_den_engelske_virusvariant_af_sars-cov-2.pdf?la=da

data_denmark = read.csv(".\\data_denmark.csv")
data_denmark = data_denmark[data_denmark$REGION!="Midtjylland",]
# since this data is aggregated by week, we will compare this data also with 
# by-week aggregated data from the UK
data_agbyweekregion_UK = data_agbyweekregion
data_agbyweekregion_UK$COUNTRY = "UK"
data_agbyweekregion_UK=data_agbyweekregion_UK[data_agbyweekregion_UK$sample_date>="2020-08-01",]
data_agbyweekregion_UK$week = as.numeric(as.character(data_agbyweekregion_UK$week))
data_agbyweekregion_UK = data_agbyweekregion_UK[data_agbyweekregion_UK$variant_lineage=="20B.501Y.V1",
                                            c("week","COUNTRY","nhs_name","count","total")] 
head(data_agbyweekregion_UK)
colnames(data_agbyweekregion_UK) = colnames(data_denmark)
data_denmark_uk = rbind(data_agbyweekregion_UK,data_denmark)
data_denmark_uk$sample_date = as.Date(paste(2020, data_denmark_uk$WEEK, 1, sep="-"), "%Y-%U-%u")-3.5
data_denmark_uk$sample_date_num = as.numeric(data_denmark_uk$sample_date)
data_denmark_uk$obs = factor(1:nrow(data_denmark_uk)) 
data_denmark_uk$prop = data_denmark_uk$B117/data_denmark_uk$TOTAL
# for observation-level random effect to deal with overdispersion
head(data_denmark_uk)
data_denmark = data_denmark_uk[data_denmark_uk$COUNTRY=="DENMARK",] 
levels_region_denmark = c("Sjælland","Syddanmark",
                          "Nordjylland","Hovedstaden")
data_denmark$REGION = factor(data_denmark$REGION, levels=levels_region_denmark)
# we leave out data for Midtjylland (too few data points)


# 4.1 MODEL FOR DENMARK ONLY ####

bGLMMfit_denm = glmer(cbind(B117, TOTAL-B117) ~ (1|obs) +  
                          REGION*scale(sample_date_num), 
                       family=binomial(logit), data=data_denmark) 
bGLMMfit_denm2 = glmer(cbind(B117, TOTAL-B117) ~ (1|obs) + 
                     REGION+scale(sample_date_num), 
                   family=binomial(logit), data=data_denmark) 
BIC(bGLMMfit_denm, bGLMMfit_denm2) # model with different slopes across regions has better BIC
plot(Effect(c("sample_date_num","REGION"), 
            bGLMMfit_denm, xlevels=list(sample_date_num=seq(as.numeric(as.Date("2020-10-01")),
                                                              as.numeric(as.Date("2020-12-31")))),
            x.var="sample_date_num", residuals=TRUE, confint=TRUE, se=TRUE), 
     partial.residuals=TRUE, use.spline=FALSE, residuals.pch=16, 
     residuals.color=alpha("steelblue",0.2), ylab="VOC 20B.501Y.V1 (proportion)",
     smooth.residuals=FALSE, span=0.2)
# graph2png(file="fit binomial GLMM growth VOC strain 20B_501Y_V1_denmark_effect plot with partial residuals.png", width=8, height=6)


# PLOT MODEL FIT

extrapolate = 60 # 60 nr of days to extrapolate fit into the future
fitdenm_preds = as.data.frame(emmeans(bGLMMfit_denm, ~ sample_date_num*REGION, 
                                      at=list(sample_date_num=seq(as.numeric(as.Date("2020-09-01")),
                                                                  max(data$sample_date_num[data$sample_date>="2020-09-01"])+extrapolate)), 
                                      type="response"))
fitdenm_preds$sample_date = as.Date(fitdenm_preds$sample_date_num, origin="1970-01-01")
fitdenm_preds$REGION = factor(fitdenm_preds$REGION, levels=levels_region_denmark)

n = length(levels(data_denmark_uk$REGION))
reg_cols = hcl(h = seq(290, 0, length = n + 1), l = 50, c = 255)[1:n]
reg_cols_denm = reg_cols[1:4]

plot_fitGLMMdenmark = qplot(data=fitdenm_preds, x=sample_date, y=prob, geom="blank") +
  facet_wrap(~REGION) +
  geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL, fill=REGION), alpha=I(0.3)) +
  geom_line(aes(y=prob, colour=REGION), alpha=I(0.8)) +
  # labs(tag = "@TWenseleers\ndata COG-UK") +
  # theme(plot.tag.position = "bottomright",
  #       plot.tag = element_text(vjust = 1, hjust = 1, size=8)) +
  # ylab("Relative abundance of VOC strain 202012/01 (%)") +
  ylab("Relative abundance (%)") +
  theme_hc() + xlab("") + 
  ggtitle("GROWTH OF VOC 20B.501Y.V1 IN DENMARK") +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01")),
                     labels=c("M","A","M","J","J","A","S","O","N","D","J","F")) +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99")) +
  coord_cartesian(# xlim=c(as.Date("2020-09-01"),as.Date("2021-02-01")), 
    xlim=c(as.Date("2020-08-01"),as.Date("2021-01-31")), 
    ylim=c(0.0001,0.99), expand=c(0,0)) +
  scale_color_manual("", values=reg_cols_denm) +
  scale_fill_manual("", values=reg_cols_denm) +
  geom_point(data=data_denmark, 
             aes(x=sample_date, y=prop, colour=REGION, size=TOTAL), alpha=I(0.5)) +
  scale_size_continuous("total number\nof sequences\nper week", trans="log2", 
                      range=c(0.01, 4), limits=c(100,1600), breaks=c(100,200,400,800,1600)) +
  guides(fill=FALSE) + guides(colour=FALSE) + theme(legend.position = "right")
  # theme(legend.position = "none")
plot_fitGLMMdenmark

saveRDS(plot_fitGLMMdenmark, file = ".\\multinomial_fits\\fit binomial GLMM growth VOC strain 20B_501Y_V1 DENMARK.rds")
graph2ppt(file=".\\multinomial_fits\\fit binomial GLMM growth VOC strain 20B_501Y_V1 DENMARK.pptx", width=6, height=6)
ggsave(file=".\\multinomial_fits\\fit binomial GLMM growth VOC strain 20B_501Y_V1 DENMARK.png", width=8, height=6)


#  CALCULATE GROWTH RATES & SELECTIVE ADVANTAGE OF VOC STRAIN IN DIFFERENT REGIONS IN DENMARK

# growth rate of strain = slope of logistic regression = selective advantage s / generation time
# where selective advantage s is the prop increase in the R value compared to other circulating strains
# (assuming identical generation times)

bGLMM_VOC_growthrates_DENM = as.data.frame(emtrends(bGLMMfit_denm, ~ REGION, 
                                                    var="sample_date_num"))[,-c(3,4)] 
# sample_date_num.trend = 
# logistic growth rate of VOC 202012/01 strain = growth rate VOC strain - growth rate all other strains
colnames(bGLMM_VOC_growthrates_DENM)[2] = "logistic_growth_rate"

# for our two sets of default estimated growth rates these logistic growth rates 
# translate to a proportional increase in R value
# compared to all other circulating strains of 
# (assuming that generation times would be the same)
bGLMM_VOC_growthrates_DENM = delta_R.from.delta_r_df(bGLMM_VOC_growthrates_DENM)
bGLMM_VOC_growthrates_DENM
# REGION logistic_growth_rate   asymp.LCL  asymp.UCL  delta_R1 delta_R1.LCL delta_R1.UCL  delta_R2 delta_R2.LCL
# 1    Sjælland           0.31013668  0.07091321 0.54936015 1.7057517   0.39002264    3.0214807 1.1164311   0.25527361
# 2  Syddanmark           0.12005460 -0.01410516 0.25421437 0.6603003  -0.07757837    1.3981790 0.4321730  -0.05077580
# 3 Nordjylland           0.11499083  0.05145337 0.17852830 0.6324496   0.28299353    0.9819056 0.4139444   0.18522202
# 4 Hovedstaden           0.02824482 -0.01271628 0.06920593 0.1553465  -0.06993955    0.3806326 0.1016758  -0.04577612
# delta_R2.UCL
# 1    1.9775885
# 2    0.9151217
# 3    0.6426668
# 4    0.2491277

# on average across all regions we get
bGLMM_VOC_growthrates_DENM_avg = as.data.frame(emtrends(bGLMMfit_denm, ~ 1, var="sample_date_num"))[,-c(3,4)] 
colnames(bGLMM_VOC_growthrates_DENM_avg)[2] = "logistic_growth_rate"

# for our two sets of default estimated growth rates these logistic growth rates 
# translate into a proportional increase in R value
# compared to all other circulating strains of 
# (assuming that generation times would be the same)
bGLMM_VOC_growthrates_DENM_avg = delta_R.from.delta_r_df(bGLMM_VOC_growthrates_DENM_avg)
bGLMM_VOC_growthrates_DENM_avg
# 1           logistic_growth_rate  asymp.LCL asymp.UCL delta_R1 delta_R1.LCL delta_R1.UCL  delta_R2 delta_R2.LCL
# 1   overall            0.1433567 0.06964378 0.2170697 0.788462    0.3830408     1.193883 0.5160561    0.2507039
# delta_R2.UCL
# 1    0.7814082


# TUKEY POSTHOC TESTS TO TEST FOR DIFFERENCES IN GROWTH RATES IN THE VOC STRAIN ACROSS DIFFERENT REGIONS

tukey_VOC_DENM = as.data.frame(emtrends(bGLMMfit_denm, pairwise ~ REGION, 
                                        var="sample_date_num", adjust="tukey")$contrasts)[,-4]
colnames(tukey_VOC_DENM)[2] = "diff_logistic_growth_rate"
tukey_VOC_DENM
#                    contrast diff_logistic_growth_rate         SE    z.ratio   p.value
# 1     Sjælland - Syddanmark                0.19008207 0.13824039 1.37501112 0.5150779
# 2    Sjælland - Nordjylland                0.19514584 0.12369192 1.57767654 0.3914041
# 3    Sjælland - Hovedstaden                0.28189185 0.12365440 2.27967508 0.1027548
# 4  Syddanmark - Nordjylland                0.00506377 0.07345059 0.06894117 0.9998825
# 5  Syddanmark - Hovedstaden                0.09180978 0.07140829 1.28570197 0.5720520
# 6 Nordjylland - Hovedstaden                0.08674601 0.03815799 2.27333804 0.1042759

# no sign differences in slopes across regions
tukey_VOC_DENM[tukey_VOC_DENM$p.value<0.05,]





# 4.2 MODEL FOR DENMARK & UK COMBINED ####

bGLMMfit_denmUK = glmer(cbind(B117, TOTAL-B117) ~ (1|REGION/obs) +  
                     COUNTRY*scale(sample_date_num), 
                   family=binomial(logit), data=data_denmark_uk) 
bGLMMfit_denmUK2 = glmer(cbind(B117, TOTAL-B117) ~ (sample_date_num||REGION/obs) +  
                          COUNTRY*scale(sample_date_num), 
                        family=binomial(logit), data=data_denmark_uk) 
BIC(bGLMMfit_denmUK, bGLMMfit_denmUK2)
#                  df      BIC
# bGLMMfit_denmUK   6 662.3898
# bGLMMfit_denmUK2  8 673.0838

# model with random intercepts in function of REGION fits best, so we continue with that

summary(bGLMMfit_denmUK) # no sign COUNTRY x sample_date interaction effect (p=0.07)
plot(Effect(c("sample_date_num","COUNTRY"), 
            bGLMMfit_denmUK, xlevels=list(sample_date_num=seq(as.numeric(as.Date("2020-10-01")),
                                                           as.numeric(as.Date("2020-12-31")))),
            x.var="sample_date_num", residuals=TRUE, confint=TRUE, se=TRUE), 
     partial.residuals=TRUE, use.spline=FALSE, residuals.pch=16, 
     residuals.color=alpha("steelblue",0.2), ylab="VOC 20B.501Y.V1 (proportion)",
     smooth.residuals=FALSE, span=0.2)
# graph2png(file="fit binomial GLMM growth VOC strain 20B_501Y_V1_denmark UK_effect plot with partial residuals.png", width=8, height=6)


# PLOT MODEL FIT

extrapolate = 60 # 60 nr of days to extrapolate fit into the future
fitdenmUK_preds = as.data.frame(emmeans(bGLMMfit_denmUK, ~ sample_date_num*COUNTRY|REGION, 
                                      at=list(sample_date_num=seq(as.numeric(as.Date("2020-08-01")),
                                                                  max(data$sample_date_num[data$sample_date>="2020-09-01"])+extrapolate)), 
                                      type="response"))
fitdenmUK_preds$sample_date = as.Date(fitdenmUK_preds$sample_date_num, origin="1970-01-01")
fitdenmUK_preds$COUNTRY = factor(fitdenmUK_preds$COUNTRY, levels=c("UK", "DENMARK"))
fitdenmUK_preds$REGION = factor(fitdenmUK_preds$REGION, levels=c(c("Nordjylland",
                                                                   "Sjælland",
                                                                   "Hovedstaden",
                                                                   "Syddanmark"), levels_nhs_name))
data_denmark_uk$REGION = factor(data_denmark_uk$REGION, levels=c(c("Nordjylland",
                                                                   "Sjælland",
                                                                   "Hovedstaden",
                                                                   "Syddanmark"), levels_nhs_name))

n = length(levels(data_denmark_uk$REGION))
reg_cols = hcl(h = seq(290, 0, length = n + 1), l = 50, c = 255)[1:n]
reg_cols[5:n] = rev(reg_cols[5:n])

fitdenmUK_preds_country = as.data.frame(emmeans(bGLMMfit_denmUK, ~ sample_date_num*COUNTRY, 
                                        at=list(sample_date_num=seq(as.numeric(as.Date("2020-08-01")),
                                                                    max(data$sample_date_num[data$sample_date>="2020-09-01"])+extrapolate)), 
                                        type="response"))
fitdenmUK_preds_country$sample_date = as.Date(fitdenmUK_preds_country$sample_date_num, origin="1970-01-01")
fitdenmUK_preds_country$COUNTRY = factor(fitdenmUK_preds_country$COUNTRY, levels=c("UK", "DENMARK"))
plot_denmUK = qplot(data=fitdenmUK_preds, x=sample_date, y=prob, geom="blank") +
  facet_wrap(~COUNTRY, ncol=1) +
  geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL, fill=REGION), alpha=I(0.3)) +
  geom_line(aes(y=prob, colour=REGION), alpha=I(0.8)) +
  # labs(tag = "@TWenseleers\ndata COG-UK & Statens Serum Institut DK") +
  # theme(plot.tag.position = "bottomright",
  #     plot.tag = element_text(vjust = 1, hjust = 1, size=8)) +
  ylab("Relative abundance (%)") +
  theme_hc() + xlab("") + 
  ggtitle("GROWTH OF VOC 20B.501Y.V1 IN DENMARK & UK") +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01")),
                     labels=c("M","A","M","J","J","A","S","O","N","D","J","F")) +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99")) +
  coord_cartesian(# xlim=c(as.Date("2020-09-01"),as.Date("2021-02-01")), 
    xlim=c(as.Date("2020-08-01"),as.Date("2021-02-011")), 
    ylim=c(0.0001,0.99), expand=c(0,0)) +
  scale_color_manual("region", values=reg_cols) +
  scale_fill_manual("region", values=reg_cols) +
  scale_size_continuous("total number\nof sequences\nper week", trans="log2", 
                        range=c(0.01, 4), limits=c(100,1600), breaks=c(100,200,400,800,1600))  +
  geom_point(data=data_denmark_uk, 
             aes(x=sample_date, y=prop, colour=REGION, size=TOTAL), alpha=I(0.4)) +
  theme(legend.direction = "vertical", legend.box = "vertical", legend.position="right", 
        legend.key.size=unit(0.45, "cm"))
plot_denmUK

saveRDS(plot_denmUK, file = ".\\multinomial_fits\\fit binomial GLMM growth VOC strain 20B_501Y_V1 DENMARK UK.rds")
graph2ppt(file=".\\multinomial_fits\\fit binomial GLMM growth VOC strain 20B_501Y_V1 DENMARK UK.pptx", width=6, height=6)
ggsave(file=".\\multinomial_fits\\fit binomial GLMM growth VOC strain 20B_501Y_V1 DENMARK UK.png", width=6, height=6)


#  CALCULATE GROWTH RATES & SELECTIVE ADVANTAGE OF VOC STRAIN IN DIFFERENT REGIONS & COUNTRIES

# TUKEY POSTHOC TESTS TO TEST FOR DIFFERENCES IN GROWTH RATES IN THE VOC STRAIN BETWEEN UK & DENMARK

tukey_VOC_DENMUK = as.data.frame(emtrends(bGLMMfit_denmUK, revpairwise ~ COUNTRY, 
                                          var="sample_date_num", adjust="tukey")$contrasts)[,-4]
colnames(tukey_VOC_DENMUK)[2] = "diff_logistic_growth_rate"
# no sign differnce in slope between Denmark & UK
tukey_VOC_DENMUK
#       contrast diff_logistic_growth_rate        SE  z.ratio    p.value
# 1 UK - DENMARK                 0.0258449 0.0142566 1.812837 0.06985692

# recalculate logistic growth rates to differences in R value in both countries
bGLMM_VOC_growthrates_DENMUK = as.data.frame(emtrends(bGLMMfit_denmUK, ~ COUNTRY, var="sample_date_num"))[,-c(3,4)] 
colnames(bGLMM_VOC_growthrates_DENMUK)[2] = "logistic_growth_rate"

# for our two sets of default estimated growth rates these logistic growth rates 
# translate into a proportional increase in R value of
bGLMM_VOC_growthrates_DENMUK = delta_R.from.delta_r_df(bGLMM_VOC_growthrates_DENMUK)
bGLMM_VOC_growthrates_DENMUK
#   COUNTRY logistic_growth_rate  asymp.LCL asymp.UCL  delta_R1 delta_R1.LCL delta_R1.UCL  delta_R2 delta_R2.LCL
# 1 DENMARK           0.08228956 0.05545675 0.1091224 0.4525926    0.3050121    0.6001730 0.2962262    0.1996334
# 2      UK           0.10813445 0.10004652 0.1162224 0.5947395    0.5502558    0.6392231 0.3892628    0.3601478
# delta_R2.UCL
# 1    0.3928191
# 2    0.4183778
