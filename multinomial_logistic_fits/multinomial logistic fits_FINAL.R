# MULTINOMIAL AND LOGISTIC FITS TO DETERMINE GROWTH RATE & 
# COMPETITIVE ADVANTAGE OF VOC 202012/01 COMPARED TO OTHER CIRCULATING VARIANTS
# ACROSS DIFFERENT REGIONS IN THE UK AS WELL AS DENMARK
# T. Wenseleers, 15 Jan. 2021

library(lme4)
library(emmeans)
library(effects)
library(ggplot2)
library(ggthemes)
library(MASS)
library(car)
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
library(afex)



# FUNCTION TO CALCULATE EXPECTED MULTIPLICATIVE EFFECT M ON R VALUES 
# WHEN COMPARING TWO COMPETING STRAINS FROM THE DIFFERENCE IN 
# THEIR MALTHUSIAN GROWTH RATES delta_r (= r_new-r_old) AND MEAN GENERATION
# TIME g. NOTE THAT THE DIFFERENCE IN MALTHUSIAN GROWTH RATES IS ALSO 
# THE SLOPE OF A FIT ON A LOG-ODDS SCALE OF THE PROPORTION OF THE NEW TYPE.

# GIVEN THAT THE EFFECTIVE REPRODUCTION NUMBER R = (1 + k * r * g)^(1 / k) 
# (Park et al. 2020) WHEN GENERATION TIME IS GAMMA DISTRIBUTED (with mean g and k=(SD/g)^2), 
# WHICH IS APPROX EQUAL TO exp(r*g) WITH r=MALTHUSIAN GROWTH RATE, IT FOLLOWS THAT
# THE EXPECTED MULTIPLICATIVE DIFFERENCE IN THE R VALUE OF TWO COMPETING VARIANTS,
# ASSUMING IDENTICAL GENERATION TIMES, EQUALS exp((r_new-r_old)*g) = exp(delta_r*g),
# WHERE THE DIFFERENCE IN MALTHUSIAN GROWTH RATE delta_r IS SOMETIMES REFERRED TO AS
# THE SELECTION RATE (TRAVISANO & LENSKI 1996) AND delta_r*g IS THE DIMENSIONLESS
# SELECTION COEFFICIENT sT OF CHEVIN (2011).

M.from.delta_r = function (delta_r, g=5.5) { 
    delta_R = exp(delta_r*g)
    return( delta_R ) 
}

M.from.delta_r(0.1, 5.5)

# function to calculate Rt from r if gen time is gamma distributed with mean g and SD sigma
# simplifies to exp(r*g) for k smallish
Rt.from.r = function(r, g=4.7, sigma=2.9) {
  k <- (sigma / g)^2
  Rt <- (1 + k * r * g)^(1 / k)
  return(Rt) }
  

# This function calculates the expected multiplicative effect on R M for two sets 
# of defaults for generation time : gamma(mean=5.5d) (SD=2.1d) (Ferretti et al. 2020) or
# gamma(mean=3.6d) (SD=3.1d) (Abbott et al. 2020, Ganyani et al. 2020).
# It works on an input dataframe df with delta_r values and returns the original
# data frame plus two sets of estimates M1 and M1 as a dataframe with extra columns
# with column names coln
M.from.delta_r_df = function (df, g1=5.5, g2=3.6, 
                                    coln=c("M1","M1.LCL","M1.UCL",
                                           "M2","M2.LCL","M2.UCL")) { 
  df_num = df[,which(unlist(lapply(df, is.numeric))), drop=F]
  df_nonnum = df[,which(!unlist(lapply(df, is.numeric))), drop=F]
  df_out1 = apply(df_num, 2, function (delta_r) M.from.delta_r(delta_r, g1))
  if (class(df_out1)[1]=="numeric") df_out1=as.data.frame(t(df_out1), check.names=F)
  df_out2 = apply(df_num, 2, function (delta_r) M.from.delta_r(delta_r, g2))
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
levels_variants = c("B","B.1.98","B.40","B.1.1",
                     "B.1.1.257","B.1.1.1","B.1.1.315",
                     "B.1.177","VOC 202012/01", # "501Y.WALES", 
                     "minority variants")
# colours to use for variant lineages
n = length(levels_variants)
lineage_cols = hcl(h = seq(15, 375, length = n + 1), l = 65, c = 200)[1:n]
lineage_cols[which(levels_variants=="minority variants")] = "grey75"

# read in data
data = read.csv(".\\data\\cog_metadata_microreact_public-2020-12-22-annotated.csv")
data$sample_date = as.Date(data$sample_date)
data$sample_date_num = as.numeric(data$sample_date)
# we do not consider data from Northern Ireland due to low nr of sequences there & absence of new VOC
data = data[data$nhs_name!="Northern Ireland",] 
data$nhs_name = factor(data$nhs_name, levels=levels_nhs_name) # NHS region
data$lad = as.factor(data$lad) # local authority district
unique(data[data$n501y == "Y","lineage"]) # lineages where at least some samples have n501y mutation
# "B.1.1.7"   "B.1.1.70"  "B.1"       "B.1.177"   "B.1.83"    "B.1.1"     "B.1.1.136"
# table(data$lineage, data$n501y, data$del_21765_6)
data$week = as.numeric(strftime(data$sample_date, format = "%V")) # week number
# we slightly recode some lineages into specific variants
# B.1.177 & descendant lineages = B.1.177 in nextstrain
data$variant[grepl("B.1.177",data$lineage,fixed=T)] = "B.1.177" 
# VOC 202012/01 = VOC 202012/01 = B.1.1.7 + n501y == "Y" + del_21765_6=="del"
data$variant = data$lineage 
sum(data$lineage=="B.1.1.7"&data$n501y =="Y"&data$del_21765_6=="del") # 3262
data$variant[data$lineage=="B.1.1.7"&
             data$n501y =="Y"&
             data$del_21765_6=="del"] = "VOC 202012/01" # VARIANT OF CONCERN, del_21765_6=deletion Δ69/Δ70
table(data[data$lineage=="B.1.1.7","n501y"]=="Y", data[data$lineage=="B.1.1.7","del_21765_6"]=="del")
#        no_del_21765_6 del_21765_6
# FALSE               0          23
# TRUE               84        3262
# data$variant[data$lineage=="B.1.1.70"&data$n501y == "Y"]="501Y.WALES" # Welsh variant, no official name yet


# aggregate data by week and check which variant lineages reached at least 15% in some week    
data_agbyweek = as.data.frame(table(data$week, data$variant))
colnames(data_agbyweek) = c("week", "variant", "count")
data_agbyweek_sum = aggregate(count ~ week, data=data_agbyweek, sum)
data_agbyweek$total = data_agbyweek_sum$count[match(data_agbyweek$week, data_agbyweek_sum$week)]
sum(data_agbyweek[data_agbyweek$variant=="VOC 202012/01","total"]) == nrow(data) # correct
data_agbyweek$sample_date = as.Date(paste(2020, data_agbyweek$week, 1, sep="-"), "%Y-%U-%u")
data_agbyweek$variant = as.factor(data_agbyweek$variant)
data_agbyweek$prop = data_agbyweek$count/data_agbyweek$total
maxweeklyprop_lineages = aggregate(prop ~ variant, data=data_agbyweek, max)
# we select lineages that at one point in time reached at least a relative abundance of 15%
selectedlineages = as.character(maxweeklyprop_lineages$variant[maxweeklyprop_lineages$prop>=0.15])
selectedlineages
# "B"             "B.1.1"         "B.1.1.1"       "B.1.1.257"     "B.1.1.315"     "B.1.177"       "B.1.98"       
# "B.40"          "VOC 202012/01"

# we recode all other 415 minor lineages as a single category "minority variants"
data$variant[!data$variant %in% selectedlineages] = "minority variants"
length(unique(data$lineage[data$variant=="minority variants"])) # 440 lineages
data$variant = factor(data$variant, levels=levels_variants)

# aggregated data to make Muller plots of raw data
# aggregated by week for selected variant lineages
data_agbyweek = as.data.frame(table(data$week, data$variant))
colnames(data_agbyweek) = c("week", "variant", "count")
data_agbyweek_sum = aggregate(count ~ week, data=data_agbyweek, sum)
data_agbyweek$total = data_agbyweek_sum$count[match(data_agbyweek$week, data_agbyweek_sum$week)]
sum(data_agbyweek[data_agbyweek$variant=="VOC 202012/01","total"]) == nrow(data) # correct
data_agbyweek$sample_date = as.Date(paste(2020, data_agbyweek$week, 1, sep="-"), "%Y-%U-%u")-3.5
data_agbyweek$variant = factor(data_agbyweek$variant, levels=levels_variants)
data_agbyweek$sample_date_num = as.numeric(data_agbyweek$sample_date)
data_agbyweek$prop = data_agbyweek$count/data_agbyweek$total

# aggregated by day for selected variant lineages
data_agbyday = as.data.frame(table(data$sample_date, data$variant))
colnames(data_agbyday) = c("sample_date", "variant", "count")
data_agbyday_sum = aggregate(count ~ sample_date, data=data_agbyday, sum)
data_agbyday$total = data_agbyday_sum$count[match(data_agbyday$sample_date, data_agbyday_sum$sample_date)]
sum(data_agbyday[data_agbyday$variant=="VOC 202012/01","total"]) == nrow(data) # correct
data_agbyday$sample_date = as.Date(data_agbyday$sample_date)
data_agbyday$variant = factor(data_agbyday$variant, levels=levels_variants)
data_agbyday = data_agbyday[data_agbyday$total!=0,]
data_agbyday$sample_date_num = as.numeric(data_agbyday$sample_date)
data_agbyday$prop = data_agbyday$count/data_agbyday$total
data_agbyday = data_agbyday[data_agbyday$total!=0,]


# aggregated by week and nhs_name for selected variant lineages
data_agbyweekregion = as.data.frame(table(data$week, data$nhs_name, data$variant))
colnames(data_agbyweekregion) = c("week", "nhs_name", "variant", "count")
data_agbyweekregion_sum = aggregate(count ~ week + nhs_name, data=data_agbyweekregion, sum)
data_agbyweekregion$total = data_agbyweekregion_sum$count[match(interaction(data_agbyweekregion$week,data_agbyweekregion$nhs_name), 
                                                          interaction(data_agbyweekregion_sum$week,data_agbyweekregion_sum$nhs_name))]
sum(data_agbyweekregion[data_agbyweekregion$variant=="VOC 202012/01","total"]) == nrow(data) # correct
data_agbyweekregion$sample_date = as.Date(paste(2020, data_agbyweekregion$week, 1, sep="-"), "%Y-%U-%u")-3.5
data_agbyweekregion$variant = factor(data_agbyweekregion$variant, levels=levels_variants)
data_agbyweekregion$nhs_name = factor(data_agbyweekregion$nhs_name, levels=levels_nhs_name)
data_agbyweekregion$sample_date_num = as.numeric(data_agbyweekregion$sample_date)
data_agbyweekregion$prop = data_agbyweekregion$count/data_agbyweekregion$total
data_agbyweekregion = data_agbyweekregion[data_agbyweekregion$total!=0,]

# aggregated by day and nhs_name for selected variant lineages
data_agbydayregion = as.data.frame(table(data$sample_date, data$nhs_name, data$variant))
colnames(data_agbydayregion) = c("sample_date", "nhs_name", "variant", "count")
data_agbydayregion_sum = aggregate(count ~ sample_date + nhs_name, data=data_agbydayregion, sum)
data_agbydayregion$total = data_agbydayregion_sum$count[match(interaction(data_agbydayregion$sample_date,data_agbydayregion$nhs_name),
                                                              interaction(data_agbydayregion_sum$sample_date,data_agbydayregion_sum$nhs_name))]
sum(data_agbydayregion[data_agbydayregion$variant=="VOC 202012/01","total"]) == nrow(data) # correct
data_agbydayregion$sample_date = as.Date(data_agbydayregion$sample_date)
data_agbydayregion$variant = factor(data_agbydayregion$variant, levels=levels_variants)
data_agbydayregion$nhs_name = factor(data_agbydayregion$nhs_name, levels=levels_nhs_name)
data_agbydayregion$sample_date_num = as.numeric(data_agbydayregion$sample_date)
data_agbydayregion$prop = data_agbydayregion$count/data_agbydayregion$total
data_agbydayregion = data_agbydayregion[data_agbydayregion$total!=0,]

# data aggregated by day, nhs_name and lad for selected variant lineages
data_agbydayregionlad = as.data.frame(table(data$sample_date, data$nhs_name, data$lad, data$variant))
colnames(data_agbydayregionlad) = c("sample_date", "nhs_name", "lad", "variant", "count")
data_agbydayregionlad_sum = aggregate(count ~ sample_date + nhs_name + lad, data=data_agbydayregionlad, sum)
data_agbydayregionlad$total = data_agbydayregionlad_sum$count[match(interaction(data_agbydayregionlad$sample_date,data_agbydayregionlad$nhs_name,data_agbydayregionlad$lad),
                                                                    interaction(data_agbydayregionlad_sum$sample_date,data_agbydayregionlad_sum$nhs_name,data_agbydayregionlad_sum$lad))]
sum(data_agbydayregionlad[data_agbydayregionlad$variant=="VOC 202012/01","total"]) == nrow(data) # correct
data_agbydayregionlad$sample_date = as.Date(data_agbydayregionlad$sample_date)
data_agbydayregionlad$variant = factor(data_agbydayregionlad$variant, levels=levels_variants)
data_agbydayregionlad$nhs_name = factor(data_agbydayregionlad$nhs_name, levels=levels_nhs_name)
data_agbydayregionlad$sample_date_num = as.numeric(data_agbydayregionlad$sample_date)
data_agbydayregionlad$prop = data_agbydayregionlad$count/data_agbydayregionlad$total
data_agbydayregionlad = data_agbydayregionlad[data_agbydayregionlad$total!=0,]
data_agbydayregionlad$obs = factor(1:nrow(data_agbydayregionlad))

# long version of data aggregated by day & nhs_name, including all the zeros here
# (to be able to fit a Poisson mixed multinomial surrogate model)
data_agbydayregion$variant = relevel(data_agbydayregion$variant, ref="minority variants")
data_agbydayregion_wide = spread(data_agbydayregion, variant, count)
data_agbydayregion_wide[is.na(data_agbydayregion_wide)] = 0
sum(colSums(data_agbydayregion_wide[,levels_variants])) == nrow(data) # correct
data_agbydayregion_wide$index = factor(1:nrow(data_agbydayregion_wide)) # for observation index random factor
data_agbydayregion_wide$TOTAL = rowSums(data_agbydayregion_wide[,levels_variants])
data_agbydayregion_long = gather(data_agbydayregion_wide, variant, 
                                 count, all_of(levels_variants), factor_key=TRUE)
data_agbydayregion_long$TOTAL = data_agbydayregion_wide$TOTAL[match(data_agbydayregion_long$index,
                                                                    data_agbydayregion_wide$index)]
sum(data_agbydayregion_long$TOTAL) # 1241160 = 10 times bigger than nrow(data), correct
nrow(data) # 124116
data_agbydayregion_long = data_agbydayregion_long[data_agbydayregion_long$TOTAL!=0,]
data_agbydayregion_long$obs = factor(1:nrow(data_agbydayregion_long)) # for observation-level random factor
data_agbydayregion_long$variant = factor(data_agbydayregion_long$variant,
                                         levels=levels_variants)
data_agbydayregion_long$variant = relevel(data_agbydayregion_long$variant,
                                             ref="minority variants")
nrow(data_agbydayregion_long) # 73980

# long version of data aggregated by day, nhs_name and lad, including all the zeros here
# (to be able to fit a Poisson mixed multinomial surrogate model)
data_agbydayregionlad$variant = relevel(data_agbydayregionlad$variant, ref="minority variants")
data_agbydayregionlad_wide = spread(data_agbydayregionlad, variant, count)
data_agbydayregionlad_wide[is.na(data_agbydayregionlad_wide)] = 0
sum(colSums(data_agbydayregionlad_wide[,levels_variants])) == nrow(data) # correct
data_agbydayregionlad_wide$index = factor(1:nrow(data_agbydayregionlad_wide)) # for observation index random factor
data_agbydayregionlad_wide$TOTAL = rowSums(data_agbydayregionlad_wide[,levels_variants])
data_agbydayregionlad_long = gather(data_agbydayregionlad_wide, variant, 
                                    count, all_of(levels_variants), factor_key=TRUE)
data_agbydayregionlad_long$TOTAL = data_agbydayregionlad_wide$TOTAL[match(data_agbydayregionlad_long$index,
                                                                        data_agbydayregionlad_wide$index)]
sum(data_agbydayregionlad_long$TOTAL) # 1241160 = 10 times bigger than nrow(data), correct
nrow(data) # 124116
data_agbydayregionlad_long = data_agbydayregionlad_long[data_agbydayregionlad_long$TOTAL!=0,]
data_agbydayregionlad_long$obs = factor(1:nrow(data_agbydayregionlad_long)) # for observation-level random factor
data_agbydayregionlad_long$variant = factor(data_agbydayregionlad_long$variant,
                                                    levels=levels_variants)
data_agbydayregionlad_long$variant = relevel(data_agbydayregionlad_long$variant,
                                                     ref="minority variants")
nrow(data_agbydayregionlad_long) # 275300



# MULLER PLOTS OF RAW DATA ####

# Muller plots showing spread of main variant lineages
n = length(levels(data$variant))
lineage_cols = hcl(h = seq(15, 375, length = n + 1), l = 65, c = 200)[1:n]
lineage_cols[which(levels(data_agbyday$variant)=="minority variants")] = "grey75"

# Muller plots using daily & weekly aggregated data (overall across all regions)
ggplot(data=data_agbyday, aes(x=sample_date, y=count, group=variant)) + 
    geom_area(aes(lwd=I(1.2), colour=NULL, fill=variant), position="fill") +
    scale_fill_manual("", values=lineage_cols) +
    scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
                       labels=substring(months(as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01"))),1,1),
                       limits=as.Date(c("2020-03-01","2020-12-22")), expand=c(0,0)) +
    guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
    theme_hc() + theme(legend.position="bottom", # c(0.8,0.7)) 
                       axis.title.x=element_blank(),
                       axis.title.y=element_blank()) + 
    # labs(title = "MAIN SARS-CoV2 VARIANT LINEAGES IN THE UK") +
    ylab("Relative abundance")
muller_raw0 = ggplot(data=data_agbyweek, aes(x=sample_date, y=count, group=variant)) + 
    geom_area(aes(lwd=I(1.2), colour=NULL, fill=variant), position="fill") +
    scale_fill_manual("", values=lineage_cols) +
    scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
                       labels=substring(months(as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01"))),1,1),
                       limits=as.Date(c("2020-03-01","2020-12-22")), expand=c(0,0)) +
    guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
    theme_hc() + theme(legend.position="bottom",  
                       axis.title.x=element_blank()) + 
    # labs(title = "MAIN SARS-CoV2 VARIANT LINEAGES IN THE UK") +
    ylab("Relative abundance")
muller_raw0
saveRDS(muller_raw0, file = ".\\multinomial_logistic_fits\\plots\\muller plot lineages_overall_raw data.rds")
graph2ppt(file=".\\multinomial_logistic_fits\\plots\\muller plot lineages_overall_raw data.pptx", width=7, height=5)
ggsave(file=".\\multinomial_logistic_fits\\plots\\muller plot lineages_overall_raw data.png", width=7, height=5)

muller_raw = ggplot(data=data_agbyweekregion, aes(x=sample_date, y=count, group=variant)) + 
    geom_area(aes(lwd=I(1.2), colour=NULL, fill=variant), position="fill") +
    facet_wrap(~nhs_name) +
    scale_fill_manual("", values=lineage_cols) +
    scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
                       labels=substring(months(as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01"))),1,1),
                       limits=as.Date(c("2020-03-01","2020-12-22")), expand=c(0,0)) +
    guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
    theme_hc() + theme(legend.position="bottom",  
                       axis.title.x=element_blank()) + 
    # labs(title = "MAIN SARS-CoV2 VARIANT LINEAGES IN THE UK") +
    ylab("Relative abundance")
muller_raw
saveRDS(muller_raw, file = ".\\multinomial_logistic_fits\\plots\\muller plot lineages by region_raw data.rds")
graph2ppt(file=".\\multinomial_logistic_fits\\plots\\muller plot lineage by region_raw data.pptx", width=7, height=5)
ggsave(file=".\\multinomial_logistic_fits\\plots\\muller plot lineages by region_raw data.png", width=7, height=5)




# 2. MULTINOMIAL FITS FOR UK BY REGION ####

# we take the category "minority variants" as reference category
data$variant = relevel(data$variant, ref="minority variants") 
# PS if you would like to take a lineage with a low growth rate as a reference, one could use B.40

set_treatment_contrasts()

# to look at the average selective benefit across regions, we fit a model with an additive nhs_name effect
# and we further add a sample_date or a 2-knot spine in function of sample_date
set.seed(1)
mfit1 = multinom(variant~sample_date_num+nhs_name, data = data, maxit=1000) # multinomial common slopes model
mfit2 = multinom(variant~sample_date_num*nhs_name, data = data, maxit=10000) # multinomial separate-slopes model
mfit3 = multinom(variant~ns(sample_date_num, df=2)+nhs_name, data = data, maxit=1000) # with 2-knot nat cubic splines ifo time
mfit4 = multinom(variant~ns(sample_date_num, df=2)*nhs_name, data = data, maxit=10000) # with 2-knot nat cubic splines ifo time

# saveRDS(mfit1, file = ".\\multinomial_logistic_fits\\fits\\mfit1.rds")
# saveRDS(mfit2, file = ".\\multinomial_logistic_fits\\fits\\mfit2.rds")
# saveRDS(mfit3, file = ".\\multinomial_logistic_fits\\fits\\mfit3.rds")
# saveRDS(mfit4, file = ".\\multinomial_logistic_fits\\fits\\mfit4_model1.rds")
# or to directly load previously fitted models
mfit1 = readRDS(file = ".\\multinomial_logistic_fits\\fits\\mfit1.rds")
mfit2 = readRDS(file = ".\\multinomial_logistic_fits\\fits\\mfit2.rds")
mfit3 = readRDS(file = ".\\multinomial_logistic_fits\\fits\\mfit3.rds")
mfit4 = readRDS(file = ".\\multinomial_logistic_fits\\fits\\mfit4_model1.rds")

BIC(mfit1, mfit2, mfit3, mfit4)
# df      BIC
# mfit1  90 255643.1
# mfit2 162 255025.9
# mfit3  99 251968.8
# mfit4 243 250303.3


# we continue with mfit4 as this was a realistic fit, extrapolates in a stable way & 
# has the best BIC of these models

mfit = mfit4 
mfit

# # note about mfit1: 
# mfit_emtrends_date1 = emtrends(mfit1, revpairwise ~ variant|sample_date_num, var='sample_date_num', 
#                               mode="latent", adjust="Tukey", 
#                               at=list(sample_date_num=mean(as.numeric(min(data$sample_date)),
#                                                           as.numeric(max(data$sample_date))),
#                                       variant=c("VOC 202012/01","B.1.177","minority variants"))) 
# mfit_emtrends_date1_df = as.data.frame(mfit_emtrends_date1)
# # the 0.1137356 coef ifo sample data above for the VOC is equal to the difference in growth rate with the reference category other
# mean(mfit_emtrends_date1_df$sample_date_num.trend[mfit_emtrends_date1_df$variant=="VOC 202012/01"]-
#      mfit_emtrends_date1_df$sample_date_num.trend[mfit_emtrends_date1_df$variant=="minority variants"])
# # 0.1137356



# GROWTH RATE CONTRASTS

# the mean difference in growth rates between the variants and corresponding 
# expected multiplicative effects on the R value for the multinomial common slopes model mfit1
mfit_emtrends1 = emtrends(mfit1, revpairwise ~ variant, var="sample_date_num", 
                         mode="latent", adjust="Tukey", 
                         at=list(sample_date_num=as.numeric(max(data$sample_date)), # identical at any sample date with this model, so we just pick a date here
                                 variant=c("VOC 202012/01","B.1.177","minority variants")),
                                 nhs_name=levels_nhs_name[1]) # identical for all regions with this model, so we just pick one
mfit_contrasts1 = data.frame(as.data.frame(mfit_emtrends1$contrasts),
                          as.data.frame(confint(mfit_emtrends1$contrasts))[,c("lower.CL","upper.CL")])
colnames(mfit_contrasts1)[which(colnames(mfit_contrasts1) %in% c("estimate","lower.CL","upper.CL"))] = 
  c("delta_r","delta_r.lower.CL","delta_r.upper.CL")
mfit_contrasts1 = data.frame(mfit_contrasts1[,c("contrast","SE","t.ratio","p.value")],
  M.from.delta_r_df(mfit_contrasts1[,c("delta_r","delta_r.lower.CL","delta_r.upper.CL")]))
mfit_contrasts1
#                             contrast           SE  t.ratio      p.value    delta_r delta_r.lower.CL delta_r.upper.CL       M1   M1.LCL
# 1         B.1.177 - minority variants 4.028407e-07 59949.36 4.353682e-10 0.02415004       0.02414908       0.02415100 1.142050 1.142044
# 2 (VOC 202012/01) - minority variants 1.228263e-06 92598.74 4.353682e-10 0.11373562       0.11373269       0.11373855 1.869266 1.869236
# 3           (VOC 202012/01) - B.1.177 1.215052e-06 73729.86 4.353682e-10 0.08958558       0.08958268       0.08958847 1.636763 1.636737
# M1.UCL       M2   M2.LCL   M2.UCL
# 1 1.142056 1.090831 1.090828 1.090835
# 2 1.869296 1.505987 1.505971 1.506002
# 3 1.636789 1.380586 1.380572 1.380600
table2csv(mfit_contrasts1, file=".\\multinomial_logistic_fits\\tables\\mfit1_multinomial fit_growthrates_UK_homog slopes.csv")



# the mean difference in growth rates between the variants and corresponding 
# expected multiplicative effects on the R value for the multinomial separate slopes model mfit2
# on average across regions, evaluated at the 17th of December 2020 (the last sampling day)
mfit_emtrends2 = emtrends(mfit2, revpairwise ~ variant, var="sample_date_num", 
                          mode="latent", adjust="Tukey", 
                          at=list(sample_date_num=as.numeric(max(data$sample_date)),
                                  variant=c("VOC 202012/01","B.1.177","minority variants"))) 
mfit_contrasts2 = data.frame(as.data.frame(mfit_emtrends2$contrasts),
                             as.data.frame(confint(mfit_emtrends2$contrasts))[,c("lower.CL","upper.CL")])
colnames(mfit_contrasts2)[which(colnames(mfit_contrasts2) %in% c("estimate","lower.CL","upper.CL"))] = 
  c("delta_r","delta_r.lower.CL","delta_r.upper.CL")
mfit_contrasts2 = data.frame(mfit_contrasts2[,c("contrast","SE","t.ratio","p.value")],
                             M.from.delta_r_df(mfit_contrasts2[,c("delta_r","delta_r.lower.CL","delta_r.upper.CL")]))
mfit_contrasts2
#                              contrast           SE  t.ratio  p.value    delta_r delta_r.lower.CL delta_r.upper.CL       M1   M1.LCL   M1.UCL
# 1         B.1.177 - minority variants 4.596212e-07 52930.45 1.08e-13 0.02432796       0.02432687       0.02432905 1.143168 1.143162 1.143175
# 2 (VOC 202012/01) - minority variants 2.641739e-06 45897.02 1.08e-13 0.12124792       0.12124167       0.12125417 1.948118 1.948051 1.948185
# 3           (VOC 202012/01) - B.1.177 2.630108e-06 36850.18 1.08e-13 0.09691996       0.09691374       0.09692618 1.704139 1.704080 1.704197
# M2   M2.LCL   M2.UCL
# 1 1.091530 1.091526 1.091535
# 2 1.547271 1.547236 1.547305
# 3 1.417524 1.417492 1.417556
table2csv(mfit_contrasts2, file=".\\multinomial_logistic_fits\\tables\\mfit2_multinomial fit_growthrates_UK_heter slopes.csv")



# the mean difference in growth rates between the variants and corresponding 
# expected multiplicative effects on the R value for the multinomial separate slopes spline model mfit4
# on average across regions, evaluated for the relevant timeframes where the variants started invading

# the growth contrasts for the VOC - minority variants and for the VOC - B.1.177 comparisons we calculate
# for the period from 1 Nov 2020 until 17th of Dec 2020
mfit_emtrends = emtrends(mfit4, ~ variant, var="sample_date_num", 
                         mode="latent", 
                         at=list(sample_date_num=as.numeric(seq(as.Date("2020-11-01"),
                                                                max(data$sample_date), by=1)),
                                 variant=c("VOC 202012/01","B.1.177","minority variants")))
mfit_contrasts = data.frame(as.data.frame(contrast(mfit_emtrends, method="trt.vs.ctrl", ref=3, reverse=TRUE, adjust="Tukey")),
                            as.data.frame(confint(contrast(mfit_emtrends, method="trt.vs.ctrl", ref=3, reverse=TRUE, adjust="Tukey")))[,c("lower.CL","upper.CL")])
colnames(mfit_contrasts)[which(colnames(mfit_contrasts) %in% c("estimate","lower.CL","upper.CL"))] = 
  c("delta_r","delta_r.lower.CL","delta_r.upper.CL")
mfit_contrasts = data.frame(mfit_contrasts[,c("contrast","SE","t.ratio","p.value")],
                            M.from.delta_r_df(mfit_contrasts[,c("delta_r","delta_r.lower.CL","delta_r.upper.CL")]))
# for the B.1.177 - minority variants we calculate it for the period from 1 July 2020 until 30th of Sept 2020
mfit_emtrends_B117_minor = emtrends(mfit4, revpairwise ~ variant, var="sample_date_num", 
                                    mode="latent", 
                                    at=list(sample_date_num=as.numeric(seq(as.Date("2020-07-01"),
                                                                           as.Date("2020-09-30"), by=1)),
                                            variant=c("B.1.177","minority variants"))) 
mfit_contrasts_B117_minor = data.frame(as.data.frame(mfit_emtrends_B117_minor$contrasts),
                                       as.data.frame(confint(mfit_emtrends_B117_minor$contrasts))[,c("lower.CL","upper.CL")])
colnames(mfit_contrasts_B117_minor)[which(colnames(mfit_contrasts_B117_minor) %in% c("estimate","lower.CL","upper.CL"))] = 
  c("delta_r","delta_r.lower.CL","delta_r.upper.CL")
mfit_contrasts_B117_minor = data.frame(mfit_contrasts_B117_minor[,c("contrast","SE","t.ratio","p.value")],
                                       M.from.delta_r_df(mfit_contrasts_B117_minor[,c("delta_r","delta_r.lower.CL","delta_r.upper.CL")]))
mfit_contrasts_B117_minor

mfit_contrasts = rbind(mfit_contrasts, mfit_contrasts_B117_minor)
mfit_contrasts
# contrast           SE  t.ratio       p.value   delta_r delta_r.lower.CL delta_r.upper.CL       M1   M1.LCL
# 1 (VOC 202012/01) - minority variants 0.0040190615 26.55809  0.000000e+00 0.1067386       0.09769422       0.11578296 1.798697 1.711411
# 2           (VOC 202012/01) - B.1.177 0.0040179010 26.12332  0.000000e+00 0.1049609       0.09591914       0.11400265 1.781196 1.694784
# 3         B.1.177 - minority variants 0.0009922894 45.21956 2.597318e-120 0.0448709       0.04291631       0.04682548 1.279910 1.266225
# M1.UCL       M2   M2.LCL   M2.UCL
# 1 1.890434 1.468526 1.421481 1.517127
# 2 1.872014 1.459158 1.412426 1.507435
# 3 1.293744 1.175314 1.167073 1.183613

table2csv(mfit_contrasts, file=".\\multinomial_logistic_fits\\tables\\model 1_mfit4_multinomial spline fit_growthrates_UK_heter slopes.csv")





# plot model predictions of best fitting multinomial spline model
extrapolate = 60 # nr of days to show extrapolations for
newdat = expand.grid(sample_date_num =
                         seq(min(data$sample_date_num),max(data$sample_date_num)+extrapolate,1),
                     nhs_name = levels_nhs_name)
mfit_preds = data.frame(newdat, predict(mfit4, type="prob", newdata=newdat), check.names=F)
mfit_preds = gather(mfit_preds, variant, prob, levels_variants)
mfit_preds$sample_date = as.Date(mfit_preds$sample_date_num, origin="1970-01-01")
mfit_preds$variant = factor(mfit_preds$variant, levels=levels_variants)
mfit_preds$nhs_name = factor(mfit_preds$nhs_name, levels=levels_nhs_name)
muller_mfit = ggplot(data=mfit_preds, 
                     aes(x=sample_date, y=prob, group=variant)) + 
    facet_wrap(~nhs_name) +
    geom_area(aes(lwd=I(1.2), colour=NULL, fill=variant), position="stack") +
    annotate("rect", xmin=max(data$sample_date)+1, xmax=as.Date("2021-01-31"), ymin=0, ymax=1, alpha=0.4, fill="white") + # extrapolated part
    scale_fill_manual("", values=lineage_cols) +
    scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
                       labels=substring(months(as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01"))),1,1),
                       limits=as.Date(c("2020-03-01","2021-01-31")), expand=c(0,0)) +
    guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
    theme_hc() + theme(legend.position="bottom", 
                       axis.title.x=element_blank()) + 
    # labs(title = "MAIN SARS-CoV2 VARIANT LINEAGES IN THE UK") +
    ylab("Relative abundance")
muller_mfit
saveRDS(muller_mfit, file = ".\\multinomial_logistic_fits\\plots\\model 1_plot multinomial spline fit_muller plot fit.rds")
graph2ppt(file=".\\multinomial_logistic_fits\\plots\\model 1_plot multinomial spline fit_muller plot fit.pptx", width=7, height=5)
ggsave(".\\multinomial_logistic_fits\\plots\\model 1_plot multinomial spline fit_muller plot fit.png", width=7, height=5)
# plot with raw data and fit combined:
ggarrange(muller_raw +
              coord_cartesian(xlim=c(as.Date("2020-03-01"),
                                     as.Date("2021-01-31")), 
                              expand=c(0,0)) + labs(title = "") , 
          muller_mfit+labs(title = "")+coord_cartesian(xlim=c(as.Date("2020-03-01"),
                                                             as.Date("2021-01-31")), 
                                                      expand=c(0,0)),
          ncol=1, common.legend=TRUE, legend="bottom")




# plot predictions specifically for variants VOC 202012/01 and B.1.177

extrapolate = 60 # nr of days to show extrapolations for
# we calculate predictions & 95% confidence intervals using the effects package here (emmeans ran out of memory)
fit_varpreds.eff = Effect(c("sample_date_num","nhs_name"), 
                              mfit, xlevels=list(sample_date_num=seq(min(data$sample_date_num),
                                                                     max(data$sample_date_num)+extrapolate)))
fit_varpreds = data.frame(fit_varpreds.eff,check.names=F)
cols = names(fit_varpreds.eff$x)
fit_varpreds_probs = fit_varpreds[,grepl("^prob\\.", colnames(fit_varpreds))]
colnames(fit_varpreds_probs) = gsub("prob\\.|prob\\.X","", colnames(fit_varpreds_probs))
lvls = colnames(fit_varpreds_probs)
fit_varpreds_probs = gather(data.frame(fit_varpreds[,cols], fit_varpreds_probs, check.names=F), 
                                variant, prob, lvls, factor_key=TRUE)
fit_varpreds_LCL = fit_varpreds[,grepl("^L\\.prob\\.", colnames(fit_varpreds))]
colnames(fit_varpreds_LCL) = gsub("L\\.prob\\.|L\\.prob\\.X","", colnames(fit_varpreds_LCL))
fit_varpreds_LCL = gather(fit_varpreds_LCL, variant, LCL, lvls, factor_key=TRUE)[,-1,drop=F]
fit_varpreds_UCL = fit_varpreds[,grepl("^U\\.prob\\.", colnames(fit_varpreds))]
colnames(fit_varpreds_UCL) = gsub("U\\.prob\\.|U\\.prob\\.X","", colnames(fit_varpreds_UCL))
fit_varpreds_UCL = gather(fit_varpreds_UCL, variant, UCL, lvls, factor_key=TRUE)[,-1,drop=F]
fit_varpreds = data.frame(fit_varpreds_probs, fit_varpreds_LCL, fit_varpreds_UCL, check.names=F)
fit_varpreds$sample_date = as.Date(fit_varpreds$sample_date_num, origin="1970-01-01")
fit_varpreds$nhs_name = factor(fit_varpreds$nhs_name, levels=levels_nhs_name)
fit_varpreds$variant = as.character(fit_varpreds$variant)
fit_varpreds$variant[fit_varpreds$variant=="minority.variants"] = "minority variants"
fit_varpreds$variant[fit_varpreds$variant=="VOC.202012.01"] = "VOC 202012/01"
fit_varpreds$variant = factor(fit_varpreds$variant, levels=levels_variants)

# It can be seen that the fit for VOC 202012/01 corresponds perfectly to a logistical model
# with different dates of introduction for the different regions, but identical rates of spread.
# For B.1.177 clonal interference by VOC 202012/01 results in the relative abundance of that
# variant starting to decline after an initial prolonged spread from July onwards.
# This is nicely captured by the multinomial fit.
# The raw data plotted on this graph have been aggregated by week.
plotmultinom2vars = qplot(data=fit_varpreds[fit_varpreds$variant %in% c("B.1.177",
                                                                   "VOC 202012/01"),], 
      x=sample_date, y=prob, colour=NULL, fill=nhs_name, group=nhs_name, alpha=NULL, geom="blank") +
  facet_wrap(~variant, ncol=1) +
  geom_ribbon(aes(y=prob, ymin=LCL, ymax=UCL, colour=NULL, fill=nhs_name), alpha=I(0.5)) +
  geom_line(aes(y=prob, colour=nhs_name), alpha=I(0.8)) +
  geom_point(data=data_agbyweekregion[data_agbyweekregion$variant %in% c("B.1.177",
                                                                                 "VOC 202012/01"),], 
             aes(x=sample_date, y=prop, colour=nhs_name, fill=NULL, size=total, alpha=I(0.5))) +
  # labs(tag = "@TWenseleers\ndata COG-UK") +
  # theme(plot.tag.position = "bottomright",
  #       plot.tag = element_text(vjust = 1, hjust = 1, size=8)) +
  ylab("Relative abundance (%)") +
  theme_hc() + xlab("") + 
  # ggtitle("GROWTH OF VARIANT B.1.177 AND VOC") +
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
plotmultinom2vars

saveRDS(plotmultinom2vars, file = ".\\multinomial_logistic_fits\\plots\\model 1_plot multinomial spline fit_growth VOC and B1177.rds")
graph2ppt(file=".\\multinomial_logistic_fits\\plots\\model 1_plot multinomial spline fit_growth VOC and B1177.pptx", width=6, height=6)
ggsave(file=".\\multinomial_logistic_fits\\plots\\model 1_plot multinomial spline fit_growth VOC and B1177.png", width=6, height=6)







# 3. BINOMIAL GLMM FITS TO COMPARE RATES OF SPREAD OF STRAINS VOC 202012/01 AND B.1.177 ACROSS REGIONS IN UK ####

# 3.1 BINOMIAL MIXED MODEL FIT FOR VARIANT VOC 202012/01 ####

# 3.1.1 BINOMIAL MIXED MODEL FIT FOR VARIANT VOC 202012/01 vs all other variants ####

# mixed binomial GLMM with random intercept for lad (LTLA) & 
# observation-level random effect to take into account overdispersion
# and fixed effects ns_name and date, with or without interaction effect nhs_name x date

start_date = as.Date("2020-08-01")
stop_date = max(data_agbydayregionlad$sample_date) # 17 Dec 2020
focal_variant = "VOC 202012/01"
reference = "all other"

data_subs = data_agbydayregionlad[(data_agbydayregionlad$variant %in% c(focal_variant, reference))&
                                    data_agbydayregionlad$sample_date>=start_date&
                                    data_agbydayregionlad$sample_date<=stop_date,]

if (reference!="all other") { data_subs2 = data_subs[data_subs$variant==reference,]
data_subs = data_subs[data_subs$variant==focal_variant,]
data_subs$total = data_subs$count + data_subs2$count } else { 
  data_subs = data_subs[data_subs$variant==focal_variant,] } 

set_sum_contrasts()
bGLMMfit1 = glmer(cbind(count, total-count) ~  (1|lad/obs) + 
                       nhs_name+scale(sample_date_num), 
                     family=binomial(logit), data=data_subs)
bGLMMfit2 = glmer(cbind(count, total-count) ~  (1|lad/obs) + 
                       nhs_name*scale(sample_date_num), 
                     family=binomial(logit), data=data_subs)

# saveRDS(bGLMMfit1, file = ".\\multinomial_logistic_fits\\fits\\bGLMMfit_VOCvsall_fit1_homog slopes_model 2a.rds")
# saveRDS(bGLMMfit2, file = ".\\multinomial_logistic_fits\\fits\\bGLMMfit_VOCvsall_fit2_heter slopes_model S1.rds")
# or to directly load previously fitted models
bGLMMfit1 = readRDS(file = ".\\multinomial_logistic_fits\\fits\\bGLMMfit_VOCvsall_fit1_homog slopes_model 2a.rds")
bGLMMfit2 = readRDS(file = ".\\multinomial_logistic_fits\\fits\\bGLMMfit_VOCvsall_fit1_heter slopes_model S1.rds")


# check BIC values
BIC(bGLMMfit1, bGLMMfit2)
#              df      BIC
# bGLMMfit1    12 3057.597
# bGLMMfit2    20 3107.655


# bGLMMfit1 has the best BIC (homogenous slopes across regions), but since below we 
# specifically would like to test for heterogeneity of slopes across regions we
# will then use bGLMMfit2 and use bGLMMfit1 only to estimate
# average selective benefit

# plain effect plot with partial residuals
plot(Effect(c("sample_date_num","nhs_name"), 
            bGLMMfit2, xlevels=list(sample_date_num=seq(as.numeric(as.Date("2020-10-01")),
                                                        as.numeric(as.Date("2020-12-31")))),
            x.var="sample_date_num", residuals=TRUE),
     partial.residuals=TRUE, use.spline=FALSE, residuals.pch=16, 
     residuals.color=alpha("steelblue",0.2), ylab="VOC VOC 202012/01 (proportion)",
     smooth.residuals=TRUE, span=0.2)


# PLOT MODEL FIT

# of model bGLMMfit1 (homogenous slopes, model 2a in table)
extrapolate = 60 # 60 nr of days to extrapolate fit into the future
total.SD = sqrt(sum(sapply(as.data.frame(VarCorr(bGLMMfit1))$sdcor, function (x) x^2))) # see https://cran.r-project.org/web/packages/emmeans/vignettes/transformations.html#bias-adj
bGLMM_preds = as.data.frame(emmeans(bGLMMfit1, ~ sample_date_num, by="nhs_name", at=list(sample_date_num=
                                                                    seq(as.numeric(as.Date("2020-09-01")),
                                                                        max(data$sample_date_num[data$sample_date>="2020-09-01"])+extrapolate)), 
                           type="response"), bias.adjust = TRUE, sigma = total.SD)
bGLMM_preds$sample_date = as.Date(bGLMM_preds$sample_date_num, origin="1970-01-01")
bGLMM_preds$nhs_name = factor(bGLMM_preds$nhs_name, 
                              levels=levels_nhs_name)
plot_bGLMMVOC_hom = qplot(data=bGLMM_preds, x=sample_date, y=prob, geom="blank") +
    facet_wrap(~nhs_name) +
    geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL, fill=nhs_name), alpha=I(0.3)) +
    # geom_ribbon(aes(y=prob, ymin=lower.CL, ymax=upper.CL, colour=NULL, fill=nhs_name), alpha=I(0.3)) +
    geom_line(aes(y=prob, colour=nhs_name), alpha=I(0.8)) +
    # labs(tag = "@TWenseleers\ndata COG-UK") +
    # theme(plot.tag.position = "bottomright",
    #       plot.tag = element_text(vjust = 1, hjust = 1, size=8)) +
    # ylab("Relative abundance of VOC 202012/01 (%)") +
    ylab("Relative abundance (%)") +
    theme_hc() + xlab("") + 
    # ggtitle("GROWTH OF VOC VOC 202012/01 BY NHS REGION") +
    scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01")),
                       labels=c("M","A","M","J","J","A","S","O","N","D","J","F")) +
    scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                        labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(# xlim=c(as.Date("2020-09-01"),as.Date("2021-02-01")), 
       xlim=c(as.Date("2020-07-01"),as.Date("2021-01-31")), 
       ylim=c(0.0001,0.999001), expand=c(0,0)) +
   scale_color_discrete("", h=c(0, 280), c=200) +
    scale_fill_discrete("", h=c(0, 280), c=200) +
    geom_point(data=data_agbyweekregion[data_agbyweekregion$variant=="VOC 202012/01",], 
                aes(x=sample_date, y=prop, colour=nhs_name, size=total), alpha=I(0.5)) +
    scale_size_continuous("total number\nof sequences\nper week", trans="log2", 
                        range=c(0.01, 4), limits=c(100,1600), breaks=c(100,200,400,800,1600)) +
    guides(fill=FALSE) + guides(colour=FALSE) + theme(legend.position = "right")
    # geom_point(data=data_agbydayregion[data_agbydayregion$variant=="VOC 202012/01",], 
    #         aes(x=sample_date, y=prop, colour=nhs_name, size=total), alpha=I(0.4)) +
    # theme(legend.position = "none")
plot_bGLMMVOC_hom

saveRDS(plot_bGLMMVOC_hom, file = ".\\multinomial_logistic_fits\\plots\\model2a_plot VOCvsall_fit bGLMM_homog slopes.rds")
graph2ppt(file=".\\multinomial_logistic_fits\\plots\\model2a_plot VOCvsall_fit bGLMM_homog slopes.pptx", width=8, height=6)
ggsave(file=".\\multinomial_logistic_fits\\plots\\model2a_plot VOCvsall_fit bGLMM_homog slopes.png", width=8, height=6)


# of model bGLMMfit1 (heterogeneous slopes, model S1 in table S1)
extrapolate = 60 # 60 nr of days to extrapolate fit into the future
total.SD = sqrt(sum(sapply(as.data.frame(VarCorr(bGLMMfit2))$sdcor, function (x) x^2))) # see https://cran.r-project.org/web/packages/emmeans/vignettes/transformations.html#bias-adj
bGLMM_preds = as.data.frame(emmeans(bGLMMfit2, ~ sample_date_num, by="nhs_name", at=list(sample_date_num=
                                                                                           seq(as.numeric(as.Date("2020-09-01")),
                                                                                               max(data$sample_date_num[data$sample_date>="2020-09-01"])+extrapolate)), 
                                    type="response"), bias.adjust = TRUE, sigma = total.SD)
bGLMM_preds$sample_date = as.Date(bGLMM_preds$sample_date_num, origin="1970-01-01")
bGLMM_preds$nhs_name = factor(bGLMM_preds$nhs_name, 
                              levels=levels_nhs_name)
plot_bGLMMVOC_het = qplot(data=bGLMM_preds, x=sample_date, y=prob, geom="blank") +
  facet_wrap(~nhs_name) +
  geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL, fill=nhs_name), alpha=I(0.3)) +
  # geom_ribbon(aes(y=prob, ymin=lower.CL, ymax=upper.CL, colour=NULL, fill=nhs_name), alpha=I(0.3)) +
  geom_line(aes(y=prob, colour=nhs_name), alpha=I(0.8)) +
  # labs(tag = "@TWenseleers\ndata COG-UK") +
  # theme(plot.tag.position = "bottomright",
  #       plot.tag = element_text(vjust = 1, hjust = 1, size=8)) +
  # ylab("Relative abundance of VOC 202012/01 (%)") +
  ylab("Relative abundance (%)") +
  theme_hc() + xlab("") + 
  # ggtitle("GROWTH OF VOC VOC 202012/01 BY NHS REGION") +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01")),
                     labels=c("M","A","M","J","J","A","S","O","N","D","J","F")) +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(# xlim=c(as.Date("2020-09-01"),as.Date("2021-02-01")), 
    xlim=c(as.Date("2020-07-01"),as.Date("2021-01-31")), 
    ylim=c(0.0001,0.999001), expand=c(0,0)) +
  scale_color_discrete("", h=c(0, 280), c=200) +
  scale_fill_discrete("", h=c(0, 280), c=200) +
  geom_point(data=data_agbyweekregion[data_agbyweekregion$variant=="VOC 202012/01",], 
             aes(x=sample_date, y=prop, colour=nhs_name, size=total), alpha=I(0.5)) +
  scale_size_continuous("total number\nof sequences\nper week", trans="log2", 
                        range=c(0.01, 4), limits=c(100,1600), breaks=c(100,200,400,800,1600)) +
  guides(fill=FALSE) + guides(colour=FALSE) + theme(legend.position = "right")
# geom_point(data=data_agbydayregion[data_agbydayregion$variant=="VOC 202012/01",], 
#         aes(x=sample_date, y=prop, colour=nhs_name, size=total), alpha=I(0.4)) +
# theme(legend.position = "none")
plot_bGLMMVOC_het

saveRDS(plot_bGLMMVOC_het, file = ".\\multinomial_logistic_fits\\plots\\modelS1_plot VOCvsall_fit bGLMM_heter slopes.rds")
graph2ppt(file=".\\multinomial_logistic_fits\\plots\\modelS1_plot VOCvsall_fit bGLMM_heter slopes.pptx", width=8, height=6)
ggsave(file=".\\multinomial_logistic_fits\\plots\\modelS1_plot VOCvsall_fit bGLMM_heter slopes.png", width=8, height=6)




#  CALCULATE GROWTH RATES & SELECTIVE ADVANTAGE OF VOC IN DIFFERENT REGIONS

# growth rate of variant = slope of logistic regression and corresponding expected multiplicative increase in R value
# as we are calculating this for different regions, we will use model bGLMMfit2 (heterogeneous slopes) here 

bGLMM_VOC_growthrates = as.data.frame(emtrends(bGLMMfit2, ~ nhs_name, var="sample_date_num"))[,-c(3,4)] 
# sample_date_num.trend = 
# logistic growth rate of VOC 202012/01 = Malthusian growth rate VOC - Malthusian growth rate all other variants
colnames(bGLMM_VOC_growthrates)[2] = "logistic_growth_rate"
bGLMM_VOC_growthrates = M.from.delta_r_df(bGLMM_VOC_growthrates)
bGLMM_VOC_growthrates
#                   nhs_name logistic_growth_rate  asymp.LCL asymp.UCL       M1   M1.LCL   M1.UCL       M2   M2.LCL   M2.UCL
# 1               South East           0.10533350 0.09504587 0.1156211 1.784850 1.686663 1.888752 1.461116 1.407993 1.516244
# 2                   London           0.11275275 0.09930020 0.1262053 1.859189 1.726595 2.001965 1.500667 1.429723 1.575132
# 3          East of England           0.09394234 0.08333680 0.1045479 1.676457 1.581466 1.777154 1.402410 1.349876 1.456990
# 4               South West           0.11687015 0.08342033 0.1503200 1.901772 1.582193 2.285900 1.523077 1.350282 1.717985
# 5                 Midlands           0.12110534 0.10114001 0.1410707 1.946590 1.744155 2.172522 1.546477 1.439224 1.661722
# 6 North East and Yorkshire           0.15259907 0.12317527 0.1820229 2.314734 1.968878 2.721344 1.732138 1.558044 1.925686
# 7                 Scotland           0.10087762 0.06064220 0.1411130 1.741639 1.395890 2.173028 1.437865 1.243975 1.661975
# 8               North West           0.13157635 0.10033321 0.1628195 2.061987 1.736432 2.448577 1.605885 1.435050 1.797057
# 9                    Wales           0.10396307 0.07291925 0.1350069 1.771447 1.493401 2.101262 1.453925 1.300189 1.625840
table2csv(bGLMM_VOC_growthrates, file=".\\multinomial_logistic_fits\\tables\\model S1_VOCvsall_bGLMM_VOC_growthrates_UK_by region_heter slopes.csv")

# on average across all regions, now using the most parsimonious model bGLMMfit1_od, we get
bGLMM_VOC_growthrates_avg = as.data.frame(emtrends(bGLMMfit1, ~ 1, var="sample_date_num"))[,-c(3,4)] 
colnames(bGLMM_VOC_growthrates_avg)[2] = "logistic_growth_rate"
bGLMM_VOC_growthrates_avg = M.from.delta_r_df(bGLMM_VOC_growthrates_avg)
bGLMM_VOC_growthrates_avg
# 1         logistic_growth_rate asymp.LCL asymp.UCL       M1   M1.LCL   M1.UCL       M2   M2.LCL   M2.UCL
# 1 overall             0.109303 0.1033037 0.1153023 1.824245 1.765034 1.885443 1.482146 1.450478 1.514504
table2csv(bGLMM_VOC_growthrates_avg, file=".\\multinomial_logistic_fits\\tables\\model 2a_VOCvsall_bGLMM_VOC_avggrowthrate_UK_homog slopes.csv")


# TUKEY POSTHOC TESTS TO TEST FOR DIFFERENCES IN GROWTH RATE OF THE VOC ACROSS DIFFERENT NHS REGIONS

tukey_VOC = as.data.frame(emtrends(bGLMMfit2, pairwise ~ nhs_name, var="sample_date_num", adjust="tukey")$contrasts)[,-4]
colnames(tukey_VOC)[2] = "diff_logistic_growth_rate"
tukey_VOC
#                                      contrast diff_logistic_growth_rate          SE     z.ratio     p.value
# 1                         South East - London              -0.007419253 0.008582327 -0.86448027 0.994697313
# 2                South East - East of England               0.011391163 0.007421248  1.53493899 0.839203805
# 3                     South East - South West              -0.011536649 0.017805261 -0.64793486 0.999317297
# 4                       South East - Midlands              -0.015771834 0.011433777 -1.37940718 0.905852910
# 5       South East - North East and Yorkshire              -0.047265569 0.015880408 -2.97634484 0.071832592
# 6                       South East - Scotland               0.004455886 0.021166979  0.21051118 0.999999885
# 7                     South East - North West              -0.026242851 0.016767999 -1.56505563 0.823896129
# 8                          South East - Wales               0.001370430 0.016662324  0.08224724 1.000000000
# 9                    London - East of England               0.018810416 0.008669573  2.16970503 0.425655310
# 10                        London - South West              -0.004117397 0.018362167 -0.22423260 0.999999811
# 11                          London - Midlands              -0.008352581 0.012270294 -0.68071566 0.999021107
# 12          London - North East and Yorkshire              -0.039846316 0.016499782 -2.41496021 0.275525525
# 13                          London - Scotland               0.011875139 0.021638915  0.54878624 0.999801292
# 14                        London - North West              -0.018823599 0.017345610 -1.08520822 0.976465681
# 15                             London - Wales               0.008789683 0.017248899  0.50957936 0.999886473
# 16               East of England - South West              -0.022927812 0.017844791 -1.28484622 0.936007404
# 17                 East of England - Midlands              -0.027162997 0.011510226 -2.35990118 0.306339546
# 18 East of England - North East and Yorkshire              -0.058656732 0.015943800 -3.67896809 0.007232701
# 19                 East of England - Scotland              -0.006935277 0.021216887 -0.32687533 0.999996320
# 20               East of England - North West              -0.037634014 0.016816163 -2.23796678 0.380800564
# 21                    East of England - Wales              -0.010020733 0.016713706 -0.59955181 0.999615041
# 22                      South West - Midlands              -0.004235185 0.019859512 -0.21325725 0.999999873
# 23      South West - North East and Yorkshire              -0.035728920 0.022714306 -1.57296993 0.819750856
# 24                      South West - Scotland               0.015992535 0.026679959  0.59942129 0.999615663
# 25                    South West - North West              -0.014706202 0.023340870 -0.63006229 0.999444200
# 26                         South West - Wales               0.012907080 0.023264631  0.55479407 0.999784353
# 27        Midlands - North East and Yorkshire              -0.031493735 0.018122742 -1.73780192 0.723050478
# 28                        Midlands - Scotland               0.020227720 0.022897161  0.88341606 0.993853298
# 29                      Midlands - North West              -0.010471017 0.018913775 -0.55361856 0.999787761
# 30                           Midlands - Wales               0.017142264 0.018820714  0.91081904 0.992444198
# 31        North East and Yorkshire - Scotland               0.051721455 0.025385310  2.03745614 0.516893391
# 32      North East and Yorkshire - North West               0.021022718 0.021895169  0.96015328 0.989266846
# 33           North East and Yorkshire - Wales               0.048635999 0.021805038  2.23049370 0.385612855
# 34                      Scotland - North West              -0.030698737 0.025989258 -1.18120866 0.960593616
# 35                           Scotland - Wales              -0.003085456 0.025909048 -0.11908796 0.999999999
# 36                         North West - Wales               0.027613282 0.022466555  1.22908394 0.950275719

tukey_VOC[tukey_VOC$p.value<0.05,]
#                                      contrast diff_logistic_growth_rate        SE   z.ratio     p.value
# 18 East of England - North East and Yorkshire               -0.05865673 0.0159438 -3.678968 0.007232701

table2csv(tukey_VOC, file=".\\multinomial_logistic_fits\\tables\\model S1_VOCvsall_bGLMM_VOC_Tukey contrasts diff growth rates across regions.csv")


# Almost no sign differences in rate of spread of VOC, only 2 small differences 
# (slightly faster spread in North East & Yorkshire vs South East & East of England)




# 3.1.2 BINOMIAL MIXED MODEL FIT FOR VOC 202012/01 vs B.1.177 ####

# mixed binomial GLMM with random intercept for lad (LTLA) & 
# observation-level random effect to take into account overdispersion
# and fixed effects ns_name and date, with or without interaction effect nhs_name x date

start_date = as.Date("2020-08-01")
stop_date = max(data_agbydayregionlad$sample_date) # 17 Dec 2020
focal_variant = "VOC 202012/01"
reference = "B.1.177"

data_subs = data_agbydayregionlad[(data_agbydayregionlad$variant %in% c(focal_variant, reference))&
                                    data_agbydayregionlad$sample_date>=start_date&
                                    data_agbydayregionlad$sample_date<=stop_date,]

if (reference!="all other") { data_subs2 = data_subs[data_subs$variant==reference,]
                              data_subs = data_subs[data_subs$variant==focal_variant,]
                              data_subs$total = data_subs$count + data_subs2$count } else { 
                                data_subs = data_subs[data_subs$variant==focal_variant,] } 

set_sum_contrasts()
bGLMMfit1_vsB1177 = glmer(cbind(count, total-count) ~  (1|lad/obs) + 
                    nhs_name+scale(sample_date_num), 
                  family=binomial(logit), data=data_subs)
bGLMMfit2_vsB1177 = glmer(cbind(count, total-count) ~  (1|lad/obs) + 
                    nhs_name*scale(sample_date_num), 
                  family=binomial(logit), data=data_subs)


# saveRDS(bGLMMfit1_vsB1177, file = ".\\multinomial_logistic_fits\\fits\\bGLMMfit_VOCvsB1177_fit1_homog slopes_model 2b.rds")
# saveRDS(bGLMMfit2_vsB1177, file = ".\\multinomial_logistic_fits\\fits\\bGLMMfit_VOCvsB1177_fit2_heter slopes.rds")
# or to directly load previously fitted models
bGLMMfit1_vsB1177 = readRDS(file = ".\\multinomial_logistic_fits\\fits\\bGLMMfit_VOCvsB1177_fit1_homog slopes_model 2b.rds")
bGLMMfit2_vsB1177 = readRDS(file = ".\\multinomial_logistic_fits\\fits\\bGLMMfit_VOCvsB1177_fit2_heter slopes.rds")


# check BIC values
BIC(bGLMMfit1_vsB1177, bGLMMfit2_vsB1177)
#                   df      BIC
# bGLMMfit1_vsB1177 12 2910.983
# bGLMMfit2_vsB1177 20 2957.068


# bGLMMfit1_vsB1177 has the best BIC (homogenous slopes across regions)


#  GROWTH RATES & SELECTIVE ADVANTAGE OF VOC RELATIVE TO B.1.177

# on average across all regions, using the most parsimonious model bGLMMfit1_vsB1177, we get
bGLMM_VOC_growthrates_avg = as.data.frame(emtrends(bGLMMfit1_vsB1177, ~ 1, var="sample_date_num"))[,-c(3,4)] 
colnames(bGLMM_VOC_growthrates_avg)[2] = "logistic_growth_rate"
bGLMM_VOC_growthrates_avg = M.from.delta_r_df(bGLMM_VOC_growthrates_avg)
bGLMM_VOC_growthrates_avg
# 1         logistic_growth_rate asymp.LCL asymp.UCL       M1   M1.LCL   M1.UCL      M2   M2.LCL   M2.UCL
# 1 overall            0.1059532 0.0992668 0.1126397 1.790944 1.726278 1.858033 1.46438 1.429551 1.500057
table2csv(bGLMM_VOC_growthrates_avg, file=".\\multinomial_logistic_fits\\tables\\model 2b_VOCvsB1177_bGLMM_VOC_growthrates_UKavg_homog slope.csv")



# 3.1.2 BINOMIAL MIXED MODEL FIT FOR VOC 202012/01 vs minor variants ####

# mixed binomial GLMM with random intercept for lad (LTLA) & 
# observation-level random effect to take into account overdispersion
# and fixed effects ns_name and date, with or without interaction effect nhs_name x date

start_date = as.Date("2020-08-01")
stop_date = max(data_agbydayregionlad$sample_date) # 17 Dec 2020
focal_variant = "VOC 202012/01"
reference = "minority variants"

data_subs = data_agbydayregionlad[(data_agbydayregionlad$variant %in% c(focal_variant, reference))&
                                    data_agbydayregionlad$sample_date>=start_date&
                                    data_agbydayregionlad$sample_date<=stop_date,]

if (reference!="all other") { data_subs2 = data_subs[data_subs$variant==reference,]
                              data_subs = data_subs[data_subs$variant==focal_variant,]
                              data_subs$total = data_subs$count + data_subs2$count } else { 
  data_subs = data_subs[data_subs$variant==focal_variant,] } 

set_sum_contrasts()
bGLMMfit1_vsminor = glmer(cbind(count, total-count) ~  (1|lad/obs) + 
                            nhs_name+scale(sample_date_num), 
                          family=binomial(logit), data=data_subs)
bGLMMfit2_vsminor = glmer(cbind(count, total-count) ~  (1|lad/obs) + 
                            nhs_name*scale(sample_date_num), 
                          family=binomial(logit), data=data_subs)

# saveRDS(bGLMMfit1_vsminor, file = ".\\multinomial_logistic_fits\\fits\\bGLMMfit_VOCvsminor_fit1_homog slopes_model 2c.rds")
# saveRDS(bGLMMfit2_vsminor, file = ".\\multinomial_logistic_fits\\fits\\bGLMMfit_VOCvsminor_fit2_heter slopes.rds")
# or to directly load previously fitted models
bGLMMfit1_vsminor = readRDS(file = ".\\multinomial_logistic_fits\\fits\\bGLMMfit_VOCvsminor_fit1_homog slopes_model 2c.rds")
bGLMMfit2_vsminor = readRDS(file = ".\\multinomial_logistic_fits\\fits\\bGLMMfit_VOCvsminor_fit2_heter slopes.rds")


# check BIC values
BIC(bGLMMfit1_vsminor, bGLMMfit2_vsminor)
#                   df      BIC
# bGLMMfit1_vsminor 12 2759.109
# bGLMMfit2_vsminor 20 2818.063




# bGLMMfit1_vsminor has the best BIC (homogenous slopes across regions)


#  GROWTH RATES & SELECTIVE ADVANTAGE OF VOC RELATIVE TO minority variants

# on average across all regions, using the most parsimonious model bGLMMfit1_vsminor, we get
bGLMM_VOC_growthrates_avg = as.data.frame(emtrends(bGLMMfit1_vsminor, ~ 1, var="sample_date_num"))[,-c(3,4)] 
colnames(bGLMM_VOC_growthrates_avg)[2] = "logistic_growth_rate"
bGLMM_VOC_growthrates_avg = M.from.delta_r_df(bGLMM_VOC_growthrates_avg)
bGLMM_VOC_growthrates_avg
# 1         logistic_growth_rate asymp.LCL asymp.UCL       M1   M1.LCL  M1.UCL       M2   M2.LCL   M2.UCL
# 1 overall            0.1173328 0.1110999 0.1235656 1.906617 1.842364 1.97311 1.525616 1.491765 1.560235
table2csv(bGLMM_VOC_growthrates_avg, file=".\\multinomial_logistic_fits\\tables\\model 2c_VOCvsminor_bGLMM_VOC_growthrates_UKavg_homog slope.csv")




# 3.2 BINOMIAL MIXED MODEL FIT FOR VARIANT B.1.177 ####

# 3.2.1 BINOMIAL MIXED MODEL FIT FOR VARIANT B.1.177 vs all other variants ####

# mixed binomial GLMM with random intercept for lad (LTLA) & 
# observation-level random effect to take into account overdispersion
# and fixed effects ns_name and date, with or without interaction effect nhs_name x date

start_date = as.Date("2020-07-01")
stop_date = as.Date("2020-09-30")
focal_variant = "B.1.177"
reference = "all other"

data_subs = data_agbydayregionlad[(data_agbydayregionlad$variant %in% c(focal_variant, reference))&
                                    data_agbydayregionlad$sample_date>=start_date&
                                    data_agbydayregionlad$sample_date<=stop_date,]

if (reference!="all other") { data_subs2 = data_subs[data_subs$variant==reference,]
                              data_subs = data_subs[data_subs$variant==focal_variant,]
                              data_subs$total = data_subs$count + data_subs2$count } else { 
  data_subs = data_subs[data_subs$variant==focal_variant,] } 


set_sum_contrasts()

# mixed binomial GLMM with nested random intercept for lad and observation (the latter to take into account overdispersion)
# and fixed effects nhs_name and date, with or without interaction
bGLMMfit1_B1177 = glmer(cbind(count, total-count) ~  (1|lad/obs) + 
                            nhs_name+scale(sample_date_num), 
                          family=binomial(logit), data=data_subs)
bGLMMfit2_B1177 = glmer(cbind(count, total-count) ~  (1|lad/obs) + 
                            nhs_name*scale(sample_date_num), 
                          family=binomial(logit), data=data_subs)

# saveRDS(bGLMMfit1_B1177, file = ".\\multinomial_logistic_fits\\fits\\bGLMMfit_B1177vsall_fit1_homog slopes_model 2d.rds")
# saveRDS(bGLMMfit2_B1177, file = ".\\multinomial_logistic_fits\\fits\\bGLMMfit_B1177vsall_fit2_heter slopes_model S2.rds")
# or to directly load previously fitted models
bGLMMfit1_B1177 = readRDS(file = ".\\multinomial_logistic_fits\\fits\\bGLMMfit_B1177vsall_fit1_homog slopes_model 2d.rds")
bGLMMfit2_B1177 = readRDS(file = ".\\multinomial_logistic_fits\\fits\\bGLMMfit_B1177vsall_fit2_heter slopes_model S2.rds")

# check BIC values
BIC(bGLMMfit1_B1177, bGLMMfit2_B1177)
#                 df      BIC
# bGLMMfit1_B1177 12 7754.210
# bGLMMfit2_B1177 20 7755.471


# bGLMMfit1_vsB1177 has the best BIC (homogenous slopes across regions)


#  GROWTH RATES & SELECTIVE ADVANTAGE

# on average across all regions, using the most parsimonious model bGLMMfit1_B1177, we get
bGLMM_B1177_growthrates_avg = as.data.frame(emtrends(bGLMMfit1_B1177, ~ 1, var="sample_date_num"))[,-c(3,4)] 
colnames(bGLMM_B1177_growthrates_avg)[2] = "logistic_growth_rate"
bGLMM_B1177_growthrates_avg = M.from.delta_r_df(bGLMM_B1177_growthrates_avg)
bGLMM_B1177_growthrates_avg
# 1         logistic_growth_rate  asymp.LCL  asymp.UCL      M1   M1.LCL   M1.UCL       M2   M2.LCL   M2.UCL
# 1 overall           0.04755705 0.04421853 0.05089557 1.29896 1.275326 1.323031 1.186734 1.172557 1.201083
table2csv(bGLMM_B1177_growthrates_avg, file=".\\multinomial_logistic_fits\\tables\\model 2d_B1177vsall_bGLMM_growthrates_UKavg_homog slope.csv")

# growth rates per region for heterogeneous slope model bGLMMfit2_B1177
bGLMM_B1177_growthrates_region = as.data.frame(emtrends(bGLMMfit2_B1177, ~ nhs_name, var="sample_date_num"))[,-c(3,4)] 
colnames(bGLMM_B1177_growthrates_region)[2] = "logistic_growth_rate"
bGLMM_B1177_growthrates_region = M.from.delta_r_df(bGLMM_B1177_growthrates_region)
bGLMM_B1177_growthrates_region
#                   nhs_name logistic_growth_rate  asymp.LCL  asymp.UCL       M1   M1.LCL   M1.UCL       M2   M2.LCL   M2.UCL
# 1               South East           0.04997777 0.03936576 0.06058979 1.316370 1.241738 1.395488 1.197122 1.152250 1.243740
# 2                   London           0.02891424 0.01948733 0.03834114 1.172371 1.113135 1.234760 1.109702 1.072674 1.148008
# 3          East of England           0.03500830 0.02478004 0.04523655 1.212332 1.146014 1.282487 1.134316 1.093308 1.176862
# 4               South West           0.04091472 0.02759345 0.05423599 1.252362 1.163885 1.347563 1.158693 1.104438 1.215614
# 5                 Midlands           0.05828617 0.04840214 0.06817019 1.377918 1.305011 1.454898 1.233469 1.190350 1.278149
# 6 North East and Yorkshire           0.05369541 0.04507840 0.06231243 1.343563 1.281372 1.408772 1.213251 1.176192 1.251477
# 7                 Scotland           0.03547555 0.02834740 0.04260369 1.215451 1.168722 1.264049 1.136226 1.107440 1.165760
# 8               North West           0.07143431 0.06163817 0.08123045 1.481254 1.403557 1.563251 1.293257 1.248443 1.339678
# 9                    Wales           0.04170383 0.03302300 0.05038465 1.257809 1.199166 1.319319 1.161990 1.126238 1.198876
table2csv(bGLMM_B1177_growthrates_region, 
          file=".\\multinomial_logistic_fits\\tables\\model S2_B1177vsall_bGLMM_growthrates_UK_by region.csv")


# PLOT MODEL FIT

# OF HOMOGENEOUS SLOPE MODEL

extrapolate = 60 # 60 nr of days to extrapolate fit into the future
total.SD = sqrt(sum(sapply(as.data.frame(VarCorr(bGLMMfit1_B1177))$sdcor, function (x) x^2))) # see https://cran.r-project.org/web/packages/emmeans/vignettes/transformations.html#bias-adj
bGLMM_B1177_preds_hom = as.data.frame(emmeans(bGLMMfit1_B1177, ~ sample_date_num|nhs_name, at=list(sample_date_num=
                                                                               seq(as.numeric(as.Date("2020-06-01")),
                                                                                   max(data$sample_date_num[data$sample_date>="2020-09-30"])+extrapolate)), 
                                    type="response"), bias.adjust=TRUE, sigma=total.SD)
bGLMM_B1177_preds_hom$sample_date = as.Date(bGLMM_B1177_preds_hom$sample_date_num, origin="1970-01-01")
bGLMM_B1177_preds_hom$nhs_name = factor(bGLMM_B1177_preds_hom$nhs_name, 
                              levels=levels_nhs_name)
plot_bGLMM_B1177_hom = qplot(data=bGLMM_B1177_preds_hom, x=sample_date, y=prob, geom="blank") +
  facet_wrap(~nhs_name) +
  geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL, fill=nhs_name), alpha=I(0.3)) +
  # geom_ribbon(aes(y=prob, ymin=lower.CL, ymax=upper.CL, colour=NULL, fill=nhs_name), alpha=I(0.3)) +
  geom_line(aes(y=prob, colour=nhs_name), alpha=I(0.8)) +
  # labs(tag = "@TWenseleers\ndata COG-UK") +
  # theme(plot.tag.position = "bottomright",
  #       plot.tag = element_text(vjust = 1, hjust = 1, size=8)) +
  # ylab("Relative abundance of VOC 202012/01 (%)") +
  ylab("Relative abundance (%)") +
  theme_hc() + xlab("") + 
  # ggtitle("GROWTH OF VARIANT B.1.177 BY NHS REGION") +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01")),
                     labels=c("M","A","M","J","J","A","S","O","N","D","J","F")) +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(# xlim=c(as.Date("2020-09-01"),as.Date("2021-02-01")), 
    xlim=c(as.Date("2020-06-01"),as.Date("2020-11-30")), 
    ylim=c(0.0001,0.999001), expand=c(0,0)) +
  scale_color_discrete("", h=c(0, 280), c=200) +
  scale_fill_discrete("", h=c(0, 280), c=200) +
  geom_point(data=data_agbyweekregion[data_agbyweekregion$variant=="B.1.177"&
                                        data_agbyweekregion$sample_date>="2020-06-01"&
                                        data_agbyweekregion$sample_date<="2020-09-30",], 
             aes(x=sample_date, y=prop, colour=nhs_name, size=total), alpha=I(0.5)) +
  scale_size_continuous("total number\nof sequences\nper week", trans="log2", 
                        range=c(0.01, 4), limits=c(100,1600), breaks=c(100,200,400,800,1600)) +
  guides(fill=FALSE) + guides(colour=FALSE) + theme(legend.position = "right")
  # geom_point(data=data_agbydayregion[data_agbydayregion$variant=="VOC 202012/01",], 
  #         aes(x=sample_date, y=prop, colour=nhs_name, size=total), alpha=I(0.4)) +
  # theme(legend.position = "none")
plot_bGLMM_B1177_hom
saveRDS(plot_bGLMM_B1177_hom, file = ".\\multinomial_logistic_fits\\plots\\model 2d_plot B1177vsall_fit bGLMM_homog slopes.rds")
graph2ppt(file=".\\multinomial_logistic_fits\\plots\\model 2d_plot B1177vsall_fit bGLMM_homog slopess.pptx", width=8, height=6)
ggsave(file=".\\multinomial_logistic_fits\\plots\\model 2d_plot B1177vsall_fit bGLMM_homog slopes.png", width=8, height=6)


# OF HETEROGENEOUS SLOPE MODEL
total.SD = sqrt(sum(sapply(as.data.frame(VarCorr(bGLMMfit2_B1177))$sdcor, function (x) x^2))) # see https://cran.r-project.org/web/packages/emmeans/vignettes/transformations.html#bias-adj
bGLMM_B1177_preds_het = as.data.frame(emmeans(bGLMMfit2_B1177, ~ sample_date_num|nhs_name, at=list(sample_date_num=
                                                                                                     seq(as.numeric(as.Date("2020-06-01")),
                                                                                                         max(data$sample_date_num[data$sample_date>="2020-09-30"])+extrapolate)), 
                                              type="response"), bias.adjust=TRUE, sigma=total.SD)
bGLMM_B1177_preds_het$sample_date = as.Date(bGLMM_B1177_preds_het$sample_date_num, origin="1970-01-01")
bGLMM_B1177_preds_het$nhs_name = factor(bGLMM_B1177_preds_het$nhs_name, 
                                        levels=levels_nhs_name)
plot_bGLMM_B1177_het = qplot(data=bGLMM_B1177_preds_het, x=sample_date, y=prob, geom="blank") +
  facet_wrap(~nhs_name) +
  geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL, fill=nhs_name), alpha=I(0.3)) +
  # geom_ribbon(aes(y=prob, ymin=lower.CL, ymax=upper.CL, colour=NULL, fill=nhs_name), alpha=I(0.3)) +
  geom_line(aes(y=prob, colour=nhs_name), alpha=I(0.8)) +
  # labs(tag = "@TWenseleers\ndata COG-UK") +
  # theme(plot.tag.position = "bottomright",
  #       plot.tag = element_text(vjust = 1, hjust = 1, size=8)) +
  # ylab("Relative abundance of VOC 202012/01 (%)") +
  ylab("Relative abundance (%)") +
  theme_hc() + xlab("") + 
  # ggtitle("GROWTH OF VARIANT B.1.177 BY NHS REGION") +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01")),
                     labels=c("M","A","M","J","J","A","S","O","N","D","J","F")) +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(# xlim=c(as.Date("2020-09-01"),as.Date("2021-02-01")), 
    xlim=c(as.Date("2020-06-01"),as.Date("2020-11-30")), 
    ylim=c(0.0001,0.999001), expand=c(0,0)) +
  scale_color_discrete("", h=c(0, 280), c=200) +
  scale_fill_discrete("", h=c(0, 280), c=200) +
  geom_point(data=data_agbyweekregion[data_agbyweekregion$variant=="B.1.177"&
                                        data_agbyweekregion$sample_date>="2020-06-01"&
                                        data_agbyweekregion$sample_date<="2020-09-30",], 
             aes(x=sample_date, y=prop, colour=nhs_name, size=total), alpha=I(0.5)) +
  scale_size_continuous("total number\nof sequences\nper week", trans="log2", 
                        range=c(0.01, 4), limits=c(100,1600), breaks=c(100,200,400,800,1600)) +
  guides(fill=FALSE) + guides(colour=FALSE) + theme(legend.position = "right")
# geom_point(data=data_agbydayregion[data_agbydayregion$variant=="VOC 202012/01",], 
#         aes(x=sample_date, y=prop, colour=nhs_name, size=total), alpha=I(0.4)) +
# theme(legend.position = "none")
plot_bGLMM_B1177_het
saveRDS(plot_bGLMM_B1177_het, file = ".\\multinomial_logistic_fits\\plots\\model S2_plot B1177vsall_fit bGLMM_heter slopes.rds")
graph2ppt(file=".\\multinomial_logistic_fits\\plots\\model S2_plot B1177vsall_fit bGLMM_heter slopess.pptx", width=8, height=6)
ggsave(file=".\\multinomial_logistic_fits\\plots\\model S2_plot B1177vsall_fit bGLMM_heter slopes.png", width=8, height=6)





# 3.2.2 BINOMIAL MIXED MODEL FIT FOR VARIANT B.1.177 vs minority variants ####

# mixed binomial GLMM with random intercept for lad (LTLA) & 
# observation-level random effect to take into account overdispersion
# and fixed effects ns_name and date, with or without interaction effect nhs_name x date

start_date = as.Date("2020-07-01")
stop_date = as.Date("2020-09-30")
focal_variant = "B.1.177"
reference = "minority variants"

data_subs = data_agbydayregionlad[(data_agbydayregionlad$variant %in% c(focal_variant, reference))&
                                    data_agbydayregionlad$sample_date>=start_date&
                                    data_agbydayregionlad$sample_date<=stop_date,]

if (reference!="all other") { data_subs2 = data_subs[data_subs$variant==reference,]
                              data_subs = data_subs[data_subs$variant==focal_variant,]
                              data_subs$total = data_subs$count + data_subs2$count } else { 
  data_subs = data_subs[data_subs$variant==focal_variant,] } 


set_sum_contrasts()

# mixed binomial GLMM with nested random intercept for lad and observation (the latter to take into account overdispersion)
# and fixed effects nhs_name and date, with or without interaction
bGLMMfit1_B1177_vsminority = glmer(cbind(count, total-count) ~  (1|lad/obs) + 
                          nhs_name+scale(sample_date_num), 
                        family=binomial(logit), data=data_subs)
bGLMMfit2_B1177_vsminority = glmer(cbind(count, total-count) ~  (1|lad/obs) + 
                          nhs_name*scale(sample_date_num), 
                        family=binomial(logit), data=data_subs)

# saveRDS(bGLMMfit1_B1177_vsminority, file = ".\\multinomial_logistic_fits\\fits\\bGLMMfit_B1177vsminor_fit1_homog slopes.rds")
# saveRDS(bGLMMfit2_B1177_vsminority, file = ".\\multinomial_logistic_fits\\fits\\bGLMMfit_B1177vsminor_fit2_heter slopes_model 2e_best model.rds")
# or to directly load previously fitted models
bGLMMfit1_B1177_vsminority = readRDS(file = ".\\multinomial_logistic_fits\\fits\\bGLMMfit_B1177vsminor_fit1_homog slopes.rds")
bGLMMfit2_B1177_vsminority = readRDS(file = ".\\multinomial_logistic_fits\\fits\\bGLMMfit_B1177vsminor_fit2_heter slopes_model 2e_best model.rds")

# check BIC values
BIC(bGLMMfit1_B1177_vsminority, bGLMMfit2_B1177_vsminority)
#                            df      BIC
# bGLMMfit1_B1177_vsminority 12 7464.115
# bGLMMfit2_B1177_vsminority 20 7435.799


# bGLMMfit2_B1177_vsminority has the best BIC (heterogeneous slopes across regions)


#  GROWTH RATES & SELECTIVE ADVANTAGE

# on average across all regions, using the most parsimonious heterogeneous slope model bGLMMfit2_B1177_vsminority, we get
bGLMM_B1177_growthrates_avg_vsminority = as.data.frame(emtrends(bGLMMfit2_B1177_vsminority, ~ 1, var="sample_date_num"))[,-c(3,4)] 
colnames(bGLMM_B1177_growthrates_avg_vsminority)[2] = "logistic_growth_rate"
bGLMM_B1177_growthrates_avg_vsminority = M.from.delta_r_df(bGLMM_B1177_growthrates_avg_vsminority)
bGLMM_B1177_growthrates_avg_vsminority
# 1         logistic_growth_rate  asymp.LCL  asymp.UCL      M1   M1.LCL   M1.UCL       M2   M2.LCL  M2.UCL
# 1 overall            0.0419164 0.03838668 0.04544613 1.25928 1.235069 1.283966 1.162879 1.148196 1.17775
table2csv(bGLMM_B1177_growthrates_avg_vsminority, file=".\\multinomial_logistic_fits\\tables\\model 2e_B1177vsminor_bGLMM_growthrates_UKavg_heter slopes.csv")

# growth rates per region for most parsimonious heterogeneous slope model bGLMMfit2_B1177
bGLMM_B1177_growthrates_region = as.data.frame(emtrends(bGLMMfit2_B1177_vsminority, ~ nhs_name, var="sample_date_num"))[,-c(3,4)] 
colnames(bGLMM_B1177_growthrates_region)[2] = "logistic_growth_rate"
bGLMM_B1177_growthrates_region = M.from.delta_r_df(bGLMM_B1177_growthrates_region)
bGLMM_B1177_growthrates_region
#                   nhs_name logistic_growth_rate   asymp.LCL  asymp.UCL       M1   M1.LCL   M1.UCL       M2   M2.LCL   M2.UCL
# 1               South East           0.04954955 0.038935645 0.06016346 1.313273 1.238804 1.392219 1.195278 1.150467 1.241833
# 2                   London           0.02672749 0.017228857 0.03622613 1.158355 1.099394 1.220479 1.101000 1.063988 1.139300
# 3          East of England           0.03233368 0.022177556 0.04248981 1.194628 1.129728 1.263258 1.123447 1.083113 1.165282
# 4               South West           0.04143081 0.027798453 0.05506318 1.255921 1.165199 1.353708 1.160848 1.105253 1.219240
# 5                 Midlands           0.05513953 0.045015582 0.06526347 1.354277 1.280929 1.431824 1.219575 1.175926 1.264844
# 6 North East and Yorkshire           0.05161340 0.042788435 0.06043837 1.328265 1.265334 1.394326 1.204191 1.166536 1.243063
# 7                 Scotland           0.01704629 0.009075004 0.02501758 1.098290 1.051179 1.147513 1.063289 1.033210 1.094244
# 8               North West           0.06788180 0.057967456 0.07779615 1.452593 1.375505 1.534000 1.276822 1.232054 1.323217
# 9                    Wales           0.03552508 0.026884860 0.04416530 1.215783 1.159358 1.274953 1.136428 1.101624 1.172332
table2csv(bGLMM_B1177_growthrates_region, 
          file=".\\multinomial_logistic_fits\\tables\\model 2e_B1177vsminor_bGLMM_growthrates_UK_by region_heter slopes.csv")




# PLOT MODEL FIT

# USING MOST PARSIMONIOUS HETEROGENEOUS SLOPE MODEL
total.SD = sqrt(sum(sapply(as.data.frame(VarCorr(bGLMMfit2_B1177_vsminority))$sdcor, function (x) x^2))) # see https://cran.r-project.org/web/packages/emmeans/vignettes/transformations.html#bias-adj
bGLMM_B1177_preds_vsminority_het = as.data.frame(emmeans(bGLMMfit2_B1177_vsminority, ~ sample_date_num|nhs_name, at=list(sample_date_num=
                                                                                                     seq(as.numeric(as.Date("2020-06-01")),
                                                                                                         max(data$sample_date_num[data$sample_date>="2020-09-30"])+extrapolate)), 
                                              type="response"), bias.adjust=TRUE, sigma=total.SD)
bGLMM_B1177_preds_vsminority_het$sample_date = as.Date(bGLMM_B1177_preds_vsminority_het$sample_date_num, origin="1970-01-01")
bGLMM_B1177_preds_vsminority_het$nhs_name = factor(bGLMM_B1177_preds_vsminority_het$nhs_name, 
                                        levels=levels_nhs_name)
plot_bGLMM_B1177_preds_vsminority_het = qplot(data=bGLMM_B1177_preds_vsminority_het, x=sample_date, y=prob, geom="blank") +
  facet_wrap(~nhs_name) +
  geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL, fill=nhs_name), alpha=I(0.3)) +
  # geom_ribbon(aes(y=prob, ymin=lower.CL, ymax=upper.CL, colour=NULL, fill=nhs_name), alpha=I(0.3)) +
  geom_line(aes(y=prob, colour=nhs_name), alpha=I(0.8)) +
  # labs(tag = "@TWenseleers\ndata COG-UK") +
  # theme(plot.tag.position = "bottomright",
  #       plot.tag = element_text(vjust = 1, hjust = 1, size=8)) +
  # ylab("Relative abundance of VOC 202012/01 (%)") +
  ylab("Relative abundance (%)") +
  theme_hc() + xlab("") + 
  # ggtitle("GROWTH OF VARIANT B.1.177 BY NHS REGION") +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01")),
                     labels=c("M","A","M","J","J","A","S","O","N","D","J","F")) +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(# xlim=c(as.Date("2020-09-01"),as.Date("2021-02-01")), 
    xlim=c(as.Date("2020-06-01"),as.Date("2020-11-30")), 
    ylim=c(0.0001,0.999001), expand=c(0,0)) +
  scale_color_discrete("", h=c(0, 280), c=200) +
  scale_fill_discrete("", h=c(0, 280), c=200) +
  geom_point(data=data_agbyweekregion[data_agbyweekregion$variant=="B.1.177"&
                                        data_agbyweekregion$sample_date>="2020-06-01"&
                                        data_agbyweekregion$sample_date<="2020-09-30",], 
             aes(x=sample_date, y=prop, colour=nhs_name, size=total), alpha=I(0.5)) +
  scale_size_continuous("total number\nof sequences\nper week", trans="log2", 
                        range=c(0.01, 4), limits=c(100,1600), breaks=c(100,200,400,800,1600)) +
  guides(fill=FALSE) + guides(colour=FALSE) + theme(legend.position = "right")
# geom_point(data=data_agbydayregion[data_agbydayregion$variant=="VOC 202012/01",], 
#         aes(x=sample_date, y=prop, colour=nhs_name, size=total), alpha=I(0.4)) +
# theme(legend.position = "none")
plot_bGLMM_B1177_preds_vsminority_het

saveRDS(plot_bGLMM_B1177_preds_vsminority_het, file = ".\\multinomial_logistic_fits\\plots\\model2e_plot B1177vsminor_fit bGLMM_heter slopes.rds")
graph2ppt(file=".\\multinomial_logistic_fits\\plots\\model2e_plot B1177vsminor_fit bGLMM_heter slopes.pptx", width=8, height=6)
ggsave(file=".\\multinomial_logistic_fits\\plots\\model2e_plot B1177vsminor_fit bGLMM_heter slopes.png", width=8, height=6)


# multipanel for suppl Fig. S3
plot_bGLMMVOC_B1177_multipanel = ggarrange(plot_bGLMMVOC_het, 
          plot_bGLMM_B1177_preds_vsminority_het,
          ncol=1, common.legend=TRUE, legend="right")
plot_bGLMMVOC_B1177_multipanel
saveRDS(plot_bGLMMVOC_B1177_multipanel, file = ".\\multinomial_logistic_fits\\plots\\FigS3_modelS1_2e_VOCvsall_B1177vsminor_fit bGLMM_heter slopes.rds")
graph2ppt(file=".\\multinomial_logistic_fits\\plots\\FigS3_modelS1_2e_VOCvsall_B1177vsminor_fit bGLMM_heter slopes.pptx", width=6, height=8)
ggsave(file=".\\multinomial_logistic_fits\\plots\\FigS3_modelS1_2e_VOCvsall_B1177vsminor_fit bGLMM_heter slopes.png", width=6, height=8)




# TUKEY POSTHOC TESTS TO TEST FOR DIFFERENCES IN GROWTH RATES IN THE B1.177 VARIANT ACROSS DIFFERENT NHS REGIONS

tukey_B1177_vsminority = as.data.frame(emtrends(bGLMMfit2_B1177_vsminority, pairwise ~ nhs_name, var="sample_date_num", adjust="tukey")$contrasts)[,-4]
colnames(tukey_B1177_vsminority)[2] = "diff_logistic_growth_rate"
tukey_B1177_vsminority
#                                      contrast diff_logistic_growth_rate          SE    z.ratio      p.value
# 1                         South East - London               0.022822059 0.007232266  3.1555890 4.251376e-02
# 2                South East - East of England               0.017215873 0.007460991  2.3074513 3.373785e-01
# 3                     South East - South West               0.008118740 0.008776733  0.9250298 9.916182e-01
# 4                       South East - Midlands              -0.005589975 0.007457977 -0.7495296 9.980428e-01
# 5       South East - North East and Yorkshire              -0.002063849 0.007001354 -0.2947785 9.999984e-01
# 6                       South East - Scotland               0.032503263 0.006729439  4.8300107 4.790236e-05
# 7                     South East - North West              -0.018332247 0.007358228 -2.4913943 2.359469e-01
# 8                          South East - Wales               0.014024477 0.006934553  2.0224054 5.274930e-01
# 9                    London - East of England              -0.005606186 0.007053248 -0.7948375 9.970368e-01
# 10                        London - South West              -0.014703319 0.008431590 -1.7438370 7.191829e-01
# 11                          London - Midlands              -0.028412034 0.007051627 -4.0291457 1.836597e-03
# 12          London - North East and Yorkshire              -0.024885908 0.006564720 -3.7908557 4.740774e-03
# 13                          London - Scotland               0.009681203 0.006274920  1.5428409 8.352598e-01
# 14                        London - North West              -0.041154306 0.006941850 -5.9284347 1.097556e-07
# 15                             London - Wales              -0.008797583 0.006492597 -1.3550176 9.143717e-01
# 16               East of England - South West              -0.009097133 0.008628480 -1.0543147 9.803650e-01
# 17                 East of England - Midlands              -0.022805848 0.007286052 -3.1300692 4.593221e-02
# 18 East of England - North East and Yorkshire              -0.019279722 0.006815809 -2.8286769 1.070089e-01
# 19                 East of England - Scotland               0.015287390 0.006537176  2.3385310 3.187944e-01
# 20               East of England - North West              -0.035548120 0.007179657 -4.9512280 2.600273e-05
# 21                    East of England - Wales              -0.003191396 0.006746325 -0.4730570 9.999356e-01
# 22                      South West - Midlands              -0.013708714 0.008629779 -1.5885361 8.114523e-01
# 23      South West - North East and Yorkshire              -0.010182589 0.008232281 -1.2369098 9.484231e-01
# 24                      South West - Scotland               0.024384523 0.008002842  3.0469828 5.872551e-02
# 25                    South West - North West              -0.026450987 0.008531976 -3.1002181 5.022539e-02
# 26                         South West - Wales               0.005905737 0.008172771  0.7226113 9.984930e-01
# 27        Midlands - North East and Yorkshire               0.003526126 0.006815406  0.5173758 9.998726e-01
# 28                        Midlands - Scotland               0.038093237 0.006536073  5.8281537 2.009787e-07
# 29                      Midlands - North West              -0.012742272 0.007183177 -1.7739048 6.996365e-01
# 30                           Midlands - Wales               0.019614451 0.006747697  2.9068362 8.699049e-02
# 31        North East and Yorkshire - Scotland               0.034567112 0.006004723  5.7566537 3.074903e-07
# 32      North East and Yorkshire - North West              -0.016268398 0.006697217 -2.4291280 2.679035e-01
# 33           North East and Yorkshire - Wales               0.016088325 0.006231293  2.5818598 1.940965e-01
# 34                      Scotland - North West              -0.050835509 0.006413252 -7.9266354 2.073897e-13
# 35                           Scotland - Wales              -0.018478786 0.005923991 -3.1193137 4.744161e-02
# 36                         North West - Wales               0.032356723 0.006622452  4.8859125 3.620823e-05


tukey_B1177_vsminority[tukey_B1177_vsminority$p.value<0.05,]
#                               contrast diff_logistic_growth_rate          SE   z.ratio      p.value
# 1                  South East - London                0.02282206 0.007232266  3.155589 4.251376e-02
# 6                South East - Scotland                0.03250326 0.006729439  4.830011 4.790236e-05
# 11                   London - Midlands               -0.02841203 0.007051627 -4.029146 1.836597e-03
# 12   London - North East and Yorkshire               -0.02488591 0.006564720 -3.790856 4.740774e-03
# 14                 London - North West               -0.04115431 0.006941850 -5.928435 1.097556e-07
# 17          East of England - Midlands               -0.02280585 0.007286052 -3.130069 4.593221e-02
# 20        East of England - North West               -0.03554812 0.007179657 -4.951228 2.600273e-05
# 28                 Midlands - Scotland                0.03809324 0.006536073  5.828154 2.009787e-07
# 31 North East and Yorkshire - Scotland                0.03456711 0.006004723  5.756654 3.074903e-07
# 34               Scotland - North West               -0.05083551 0.006413252 -7.926635 2.073897e-13
# 35                    Scotland - Wales               -0.01847879 0.005923991 -3.119314 4.744161e-02
# 36                  North West - Wales                0.03235672 0.006622452  4.885912 3.620823e-05

table2csv(tukey_B1177_vsminority, file=".\\multinomial_logistic_fits\\tables\\model 2e_B1177vsminority_bGLMM_Tukey contrasts diff growth rates across regions.csv")


# Unlike the VOC which displaces the other variants at a consistently and very high rate in all regions,
# there were more differences in logistic growth rates across regions here, implying it was
# not displacing the other variants at a constant, uniform rate, but that there was more stochasticity
# going on, e.g. due to varying import from travellers.




# 4. COMPARISON OF RATE OF SPREAD OF VOC 202012/01 IN UK & DENMARK ####

# data from Denmark aggregated by week are provided by the Statens Serum Institut, link
# https://www.covid19genomics.dk/statistics, download on the 14th of January

data_denmark = read.csv(".\\multinomial_logistic_fits\\data\\data_denmark_14jan2021.csv")
# since this data is aggregated by week, we will compare this data also with 
# by-week aggregated data from the UK
data_agbyweekregion_UK = data_agbyweekregion
data_agbyweekregion_UK$COUNTRY = "UK"
data_agbyweekregion_UK=data_agbyweekregion_UK[data_agbyweekregion_UK$sample_date>="2020-08-01",]
data_agbyweekregion_UK$week = as.numeric(as.character(data_agbyweekregion_UK$week))
data_agbyweekregion_UK = data_agbyweekregion_UK[data_agbyweekregion_UK$variant=="VOC 202012/01",
                                            c("week","COUNTRY","nhs_name","count","total")] 
head(data_agbyweekregion_UK)
colnames(data_agbyweekregion_UK) = colnames(data_denmark)[-1]
data_denmark_uk = rbind(data_agbyweekregion_UK,data_denmark[,-1])
data_denmark_uk$sample_date = as.Date(NA)
data_denmark_uk$sample_date[data_denmark_uk$WEEK>=32] = lubridate::ymd( "2020-01-01" ) + 
                                                          lubridate::weeks( data_denmark_uk$WEEK[data_denmark_uk$WEEK>=32] - 1 )
data_denmark_uk$sample_date[data_denmark_uk$WEEK<32] = lubridate::ymd( "2021-01-01" ) + 
  lubridate::weeks( data_denmark_uk$WEEK[data_denmark_uk$WEEK<32] - 1 )
data_denmark_uk$sample_date_num = as.numeric(data_denmark_uk$sample_date)
data_denmark_uk$obs = factor(1:nrow(data_denmark_uk)) 
data_denmark_uk$prop = data_denmark_uk$B117/data_denmark_uk$TOTAL
# for observation-level random effect to deal with overdispersion
head(data_denmark_uk)
data_denmark = data_denmark_uk[data_denmark_uk$COUNTRY=="DENMARK",] 
levels_region_denmark = rev(c("Hovedstaden","Sjælland","Midtjylland","Nordjylland","Syddanmark"))
data_denmark$REGION = factor(data_denmark$REGION, levels=levels_region_denmark)
range(data_denmark_uk$sample_date[data_denmark_uk$COUNTRY=="DENMARK"]) # "2020-09-23" "2021-01-01"
range(data_denmark_uk$sample_date[data_denmark_uk$COUNTRY=="UK"])      # "2020-08-05" "2021-12-16"



# 4.1 MODEL FOR DENMARK ONLY ####

set_sum_contrasts()
bGLMMfit_denm1 = glmer(cbind(B117, TOTAL-B117) ~ (1|obs) +  
                          REGION+scale(sample_date_num), 
                       family=binomial(logit), data=data_denmark) 
bGLMMfit_denm2 = glmer(cbind(B117, TOTAL-B117) ~ (1|obs) + 
                     REGION*scale(sample_date_num), 
                   family=binomial(logit), data=data_denmark) 

# saveRDS(bGLMMfit_denm1, file = ".\\multinomial_logistic_fits\\fits\\bGLMMfit_DK_VOCvsall_fit1_homog slopes_model 3a.rds")
# saveRDS(bGLMMfit_denm2, file = ".\\multinomial_logistic_fits\\fits\\bGLMMfit_DK_VOCvsall_fit2_heter slopes.rds")
# or to directly load previously fitted models
bGLMMfit_denm1 = readRDS(file = ".\\multinomial_logistic_fits\\fits\\bGLMMfit_DK_VOCvsall_fit1_homog slopes_model 3a.rds")
bGLMMfit_denm2 = readRDS(file = ".\\multinomial_logistic_fits\\fits\\bGLMMfit_DK_VOCvsall_fit2_heter slopes.rds")

BIC(bGLMMfit_denm1, bGLMMfit_denm2) # model bGLMMfit_denm1 with homogeneous slopes across regions has best BIC
#                 df      BIC
# bGLMMfit_denm1  7 217.5952
# bGLMMfit_denm2 11 217.7621

fit = bGLMMfit_denm1

plot(Effect(c("sample_date_num","REGION"), 
            bGLMMfit_denm1, xlevels=list(sample_date_num=seq(as.numeric(as.Date("2020-10-01")),
                                                              as.numeric(as.Date("2020-12-31")))),
            x.var="sample_date_num", residuals=TRUE, confint=TRUE, se=TRUE), 
     partial.residuals=TRUE, use.spline=FALSE, residuals.pch=16, 
     residuals.color=alpha("steelblue",0.2), ylab="VOC VOC 202012/01 (proportion)",
     smooth.residuals=FALSE, span=0.2)


# PLOT MODEL FIT

extrapolate = 60 # 60 nr of days to extrapolate fit into the future
total.SD = sqrt(sum(sapply(as.data.frame(VarCorr(fit))$sdcor, function (x) x^2))) # see https://cran.r-project.org/web/packages/emmeans/vignettes/transformations.html#bias-adj
fitdenm_preds = as.data.frame(emmeans(fit, ~ sample_date_num*REGION, 
                                      at=list(sample_date_num=seq(as.numeric(as.Date("2020-09-01")),
                                                                  max(data$sample_date_num[data$sample_date>="2020-09-01"])+extrapolate)), 
                                      type="response"), bias.adjust = TRUE, sigma = total.SD)
fitdenm_preds$sample_date = as.Date(fitdenm_preds$sample_date_num, origin="1970-01-01")
fitdenm_preds$REGION = factor(fitdenm_preds$REGION, levels=levels_region_denmark)

n = length(levels(data_denmark_uk$REGION))
reg_cols = hcl(h = seq(290, 0, length = n + 1), l = 50, c = 255)[1:n]
reg_cols_denm = reg_cols[1:5]

plot_bGLMM_fitdenmark = qplot(data=fitdenm_preds, x=sample_date, y=prob, geom="blank") +
  facet_wrap(~REGION) +
  geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL, fill=REGION), alpha=I(0.3)) +
  geom_line(aes(y=prob, colour=REGION), alpha=I(0.8)) +
  # labs(tag = "@TWenseleers\ndata COG-UK") +
  # theme(plot.tag.position = "bottomright",
  #       plot.tag = element_text(vjust = 1, hjust = 1, size=8)) +
  # ylab("Relative abundance of VOC 202012/01 (%)") +
  ylab("Relative abundance (%)") +
  theme_hc() + xlab("") + 
  # ggtitle("GROWTH OF VOC VOC 202012/01 IN DENMARK") +
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
plot_bGLMM_fitdenmark

saveRDS(plot_bGLMM_fitdenmark, file = ".\\multinomial_logistic_fits\\plots\\model3a_plot VOCvsall_fit bGLMM_DENMARK.rds")
graph2ppt(file=".\\multinomial_logistic_fits\\plots\\model3a_plot VOCvsall_fit bGLMM_DENMARK.pptx", width=8, height=6)
ggsave(file=".\\multinomial_logistic_fits\\plots\\model3a_plot VOCvsall_fit bGLMM_DENMARK.png", width=8, height=6)



#  CALCULATE GROWTH RATES & SELECTIVE ADVANTAGE OF VOC IN DIFFERENT REGIONS IN DENMARK

# growth rate of variant = slope of logistic regression = selective advantage s / generation time
# where selective advantage s is the prop increase in the R value compared to other circulating variants
# (assuming identical generation times)
# for our two sets of default estimated growth rates these logistic growth rates 
# translate to a proportional increase in R value
# compared to all other circulating variants of 
# (assuming that generation times would be the same)

# in function of region (only makes sense for model that allows heterogeneous slopes across regions, so we use bGLMMfit_denm2 here)
bGLMM_VOC_growthrates_DENM = as.data.frame(emtrends(bGLMMfit_denm2, ~ REGION,
                                                    var="sample_date_num"))[,-c(3,4)]
# sample_date_num.trend =
# logistic growth rate of VOC 202012/01 = growth rate VOC - growth rate all other variants
colnames(bGLMM_VOC_growthrates_DENM)[2] = "logistic_growth_rate"

bGLMM_VOC_growthrates_DENM = M.from.delta_r_df(bGLMM_VOC_growthrates_DENM)
bGLMM_VOC_growthrates_DENM
#        REGION logistic_growth_rate  asymp.LCL  asymp.UCL       M1   M1.LCL   M1.UCL       M2   M2.LCL   M2.UCL
# 1  Syddanmark           0.16070682 0.08405484 0.23735879 2.420290 1.587724 3.689435 1.783441 1.353370 2.350179
# 2 Nordjylland           0.07079264 0.03645460 0.10513069 1.476035 1.222014 1.782860 1.290273 1.140237 1.460050
# 3 Midtjylland           0.25095580 0.10589050 0.39602110 3.975923 1.790326 8.829656 2.468081 1.464049 4.160669
# 4    Sjælland           0.11522948 0.03366279 0.19679617 1.884688 1.203393 2.951693 1.514107 1.128835 2.030874
# 5 Hovedstaden           0.05197638 0.02454919 0.07940356 1.330920 1.144560 1.547622 1.205766 1.092400 1.330897


# on average across all regions based on homogenous slope model bGLMMfit_denm1 we get
bGLMM_VOC_growthrates_DENM_avg = as.data.frame(emtrends(bGLMMfit_denm1, ~ 1, var="sample_date_num"))[,-c(3,4)] 
colnames(bGLMM_VOC_growthrates_DENM_avg)[2] = "logistic_growth_rate"
bGLMM_VOC_growthrates_DENM_avg = M.from.delta_r_df(bGLMM_VOC_growthrates_DENM_avg)
bGLMM_VOC_growthrates_DENM_avg
# 1         logistic_growth_rate  asymp.LCL asymp.UCL       M1   M1.LCL   M1.UCL       M2   M2.LCL   M2.UCL
# 1 overall           0.09586908 0.06704701 0.1246912 1.694317 1.445938 1.985362 1.412172 1.272991 1.566569
table2csv(bGLMM_VOC_growthrates_DENM_avg, file=".\\multinomial_logistic_fits\\tables\\model 3a_VOCvsall_bGLMM_growthrates_DENM_hom slope model.csv")



# TUKEY POSTHOC TESTS TO TEST FOR DIFFERENCES IN GROWTH RATES IN THE VOC ACROSS DIFFERENT REGIONS
# (ONLY MAKES SENSE FOR MODEL THAT ALLOWS FOR HETEROGENEOUS SLOPES, SO WE USE MODEL bGLMMfit_denm2 HERE)

tukey_VOC_DENM = as.data.frame(emtrends(bGLMMfit_denm2, pairwise ~ REGION, 
                                        var="sample_date_num", adjust="tukey")$contrasts)[,-4]
colnames(tukey_VOC_DENM)[2] = "diff_logistic_growth_rate"
tukey_VOC_DENM
#                     contrast diff_logistic_growth_rate         SE    z.ratio    p.value
# 1   Syddanmark - Nordjylland                0.08991417 0.04201802  2.1398956 0.20317927
# 2   Syddanmark - Midtjylland               -0.09024899 0.08242389 -1.0949372 0.80926139
# 3      Syddanmark - Sjælland                0.04547734 0.05590539  0.8134697 0.92660082
# 4   Syddanmark - Hovedstaden                0.10873044 0.04131479  2.6317557 0.06472755
# 5  Nordjylland - Midtjylland               -0.18016316 0.07467420 -2.4126559 0.11185169
# 6     Nordjylland - Sjælland               -0.04443683 0.04365923 -1.0178109 0.84728837
# 7  Nordjylland - Hovedstaden                0.01881627 0.02201799  0.8545861 0.91327218
# 8     Midtjylland - Sjælland                0.13572632 0.08252069  1.6447550 0.46863401
# 9  Midtjylland - Hovedstaden                0.19897943 0.07496480  2.6543048 0.06098077
# 10    Sjælland - Hovedstaden                0.06325310 0.04351192  1.4536961 0.59273609

# no sign pairwise differences in slopes across regions
tukey_VOC_DENM[tukey_VOC_DENM$p.value<0.05,]

table2csv(tukey_VOC_DENM, file=".\\multinomial_logistic_fits\\tables\\model3a_VOCvsall_fit_bGLMM_DENMARK_Tukey contrasts diff growth rates across regions.csv")



# 4.2 MODEL FOR DENMARK & UK COMBINED ####

levels_region_denmark = rev(c("Hovedstaden","Sjælland","Midtjylland","Nordjylland","Syddanmark"))
data_denmark_uk$REGION = factor(data_denmark_uk$REGION, levels=c(levels_region_denmark,
                                                                 levels_nhs_name))

set_sum_contrasts()
bGLMMfit_denmUK0 = glmer(cbind(B117, TOTAL-B117) ~ (1|obs) + 
                           REGION+COUNTRY*scale(sample_date_num), 
                         family=binomial(logit), data=data_denmark_uk) 
bGLMMfit_denmUK1 = glmer(cbind(B117, TOTAL-B117) ~ (1|REGION/obs) +  
                           COUNTRY*scale(sample_date_num), 
                         family=binomial(logit), data=data_denmark_uk)
bGLMMfit_denmUK2 = glmer(cbind(B117, TOTAL-B117) ~ (sample_date_num||REGION/obs) +  
                          COUNTRY*scale(sample_date_num), 
                        family=binomial(logit), data=data_denmark_uk) 
bGLMMfit_denmUK3 = glmer(cbind(B117, TOTAL-B117) ~ (1|obs) +  
                           REGION+COUNTRY+REGION:scale(sample_date_num), 
                         family=binomial(logit), data=data_denmark_uk) 
bGLMMfit_denmUK4 = glmer(cbind(B117, TOTAL-B117) ~ (1|COUNTRY/obs) + 
                           REGION*scale(sample_date_num), 
                         family=binomial(logit), data=data_denmark_uk) 

# saveRDS(bGLMMfit_denmUK0, file = ".\\multinomial_logistic_fits\\fits\\bGLMMfit_DKplusUK_VOCvsall_fit0_model 3b.rds")
# saveRDS(bGLMMfit_denmUK1, file = ".\\multinomial_logistic_fits\\fits\\bGLMMfit_DKplusUK_VOCvsall_fit1.rds")
# saveRDS(bGLMMfit_denmUK2, file = ".\\multinomial_logistic_fits\\fits\\bGLMMfit_DKplusUK_VOCvsall_fit2.rds")
# saveRDS(bGLMMfit_denmUK3, file = ".\\multinomial_logistic_fits\\fits\\bGLMMfit_DKplusUK_VOCvsall_fit3.rds")
# saveRDS(bGLMMfit_denmUK4, file = ".\\multinomial_logistic_fits\\fits\\bGLMMfit_DKplusUK_VOCvsall_fit4.rds")
# or to directly load previously fitted models
bGLMMfit_denmUK0 = readRDS(file = ".\\multinomial_logistic_fits\\fits\\bGLMMfit_DKplusUK_VOCvsall_fit0_model 3b.rds")
bGLMMfit_denmUK1 = readRDS(file = ".\\multinomial_logistic_fits\\fits\\bGLMMfit_DKplusUK_VOCvsall_fit1.rds")
bGLMMfit_denmUK2 = readRDS(file = ".\\multinomial_logistic_fits\\fits\\bGLMMfit_DKplusUK_VOCvsall_fit2.rds")
bGLMMfit_denmUK3 = readRDS(file = ".\\multinomial_logistic_fits\\fits\\bGLMMfit_DKplusUK_VOCvsall_fit3.rds")
bGLMMfit_denmUK4 = readRDS(file = ".\\multinomial_logistic_fits\\fits\\bGLMMfit_DKplusUK_VOCvsall_fit4.rds")

BIC(bGLMMfit_denmUK0, bGLMMfit_denmUK1, bGLMMfit_denmUK2, bGLMMfit_denmUK3, bGLMMfit_denmUK4)
#                  df      BIC
# bGLMMfit_denmUK0 17 751.5700
# bGLMMfit_denmUK1  6 754.1519
# bGLMMfit_denmUK2  8 765.2429
# bGLMMfit_denmUK3 29 779.8196
# bGLMMfit_denmUK4 30 785.3719



# model bGLMMfit_denmUK0 fits best, so we continue with that (heterogeneous slopes only across countries but not regions within countries)

summary(bGLMMfit_denmUK0) 
# sign COUNTRY x sample_date interaction effect (slightly slower spread in DK than in UK, z=-2.6, p=0.009)


# PLOT MODEL FIT

summary(bGLMMfit_denmUK0)
total.SD = sqrt(sum(sapply(as.data.frame(VarCorr(bGLMMfit_denmUK0))$sdcor, function (x) x^2))) # see https://cran.r-project.org/web/packages/emmeans/vignettes/transformations.html#bias-adj
extrapolate = 60 # 60 nr of days to extrapolate fit into the future
fitdenmUK_preds = as.data.frame(emmeans(bGLMMfit_denmUK0, ~ sample_date_num, by=c("COUNTRY","REGION"),
                                      at=list(sample_date_num=seq(as.numeric(as.Date("2020-08-01")),
                                                                  max(data$sample_date_num[data$sample_date>="2020-09-01"])+extrapolate)), 
                                      type="response"), bias.adjust = FALSE, sigma = total.SD)
fitdenmUK_preds$sample_date = as.Date(fitdenmUK_preds$sample_date_num, origin="1970-01-01")
fitdenmUK_preds$COUNTRY = factor(fitdenmUK_preds$COUNTRY, levels=c("UK", "DENMARK"))
fitdenmUK_preds$REGION = factor(fitdenmUK_preds$REGION, levels=c(levels_region_denmark, levels_nhs_name))
n = length(levels(data_denmark_uk$REGION))
reg_cols = hcl(h = seq(290, 0, length = n + 1), l = 50, c = 255)[1:n]
reg_cols[6:n] = rev(reg_cols[6:n])
data_denmark_uk$REGION = factor(data_denmark_uk$REGION, levels=c(levels_region_denmark, levels_nhs_name))
plot_bGLMM_denmUK = qplot(data=fitdenmUK_preds, x=sample_date, y=prob, geom="blank") +
  facet_wrap(~COUNTRY, ncol=1) +
  geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL, fill=REGION), 
              alpha=I(0.3)) +
  geom_line(aes(y=prob, colour=REGION), alpha=I(0.8)) +
  # labs(tag = "@TWenseleers\ndata COG-UK & Statens Serum Institut DK") +
  # theme(plot.tag.position = "bottomright",
  #     plot.tag = element_text(vjust = 1, hjust = 1, size=8)) +
  ylab("Relative abundance (%)") +
  theme_hc() + xlab("") + 
  # ggtitle("SPREAD OF VOC 202012/01 IN DENMARK & UK") +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01")),
                     labels=c("M","A","M","J","J","A","S","O","N","D","J","F")) +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99")) +
  coord_cartesian(# xlim=c(as.Date("2020-09-01"),as.Date("2021-02-01")), 
    xlim=c(as.Date("2020-08-01"),as.Date("2021-02-011")), 
    ylim=c(0.0001,0.99), expand=c(0,0)) +
  scale_fill_manual("region", values=reg_cols) +
  scale_size_continuous("total number\nof sequences\nper week", trans="log2", 
                        range=c(0.01, 4), limits=c(100,1600), breaks=c(100,200,400,800,1600))  +
  geom_point(data=data_denmark_uk, 
             aes(x=sample_date, y=prop, colour=REGION, size=TOTAL, fill=NULL), pch=16, alpha=I(0.4)) +
  scale_color_manual("region", values=reg_cols) +
  theme(legend.direction = "vertical", legend.box = "vertical", legend.position="right", 
        legend.key.size=unit(0.45, "cm"))
plot_bGLMM_denmUK
saveRDS(plot_bGLMM_denmUK, file = ".\\multinomial_logistic_fits\\plots\\model3b_plot VOCvsall_fit bGLMM_DENMARKplusUK.rds")
graph2ppt(file=".\\multinomial_logistic_fits\\plots\\model3b_plot VOCvsall_fit bGLMM_DENMARKplusUK.pptx", width=6, height=6)
ggsave(file=".\\multinomial_logistic_fits\\plots\\model3b_plot VOCvsall_fit bGLMM_DENMARKplusUK.png", width=6, height=6)


fitdenmUK_preds_country = as.data.frame(emmeans(bGLMMfit_denmUK0, ~ sample_date_num, by=c("COUNTRY"),
                                        at=list(sample_date_num=seq(as.numeric(as.Date("2020-08-01")),
                                                                    max(data$sample_date_num[data$sample_date>="2020-09-01"])+extrapolate)), 
                                        type="response"), bias.adjust = TRUE, sigma = total.SD)
fitdenmUK_preds_country$sample_date = as.Date(fitdenmUK_preds_country$sample_date_num, origin="1970-01-01")
fitdenmUK_preds_country$COUNTRY = factor(fitdenmUK_preds_country$COUNTRY, levels=c("UK", "DENMARK"))
data_denmark_uk$REGION = factor(data_denmark_uk$REGION, levels=c(levels_region_denmark, levels_nhs_name))
plot_denmUK_country = qplot(data=fitdenmUK_preds_country, x=sample_date, y=prob, geom="blank") +
  facet_wrap(~COUNTRY, ncol=1) +
  geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL), fill=I("darkgrey"), # fill=REGION), 
              alpha=I(0.3)) +
  geom_line(aes(y=prob), colour=I("darkgrey"), alpha=I(0.8)) +
  # labs(tag = "@TWenseleers\ndata COG-UK & Statens Serum Institut DK") +
  # theme(plot.tag.position = "bottomright",
  #     plot.tag = element_text(vjust = 1, hjust = 1, size=8)) +
  ylab("Relative abundance (%)") +
  theme_hc() + xlab("") + 
  ggtitle("SPREAD OF VOC 202012/01 IN DENMARK & UK") +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01")),
                     labels=c("M","A","M","J","J","A","S","O","N","D","J","F")) +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99")) +
  coord_cartesian(# xlim=c(as.Date("2020-09-01"),as.Date("2021-02-01")), 
    xlim=c(as.Date("2020-08-01"),as.Date("2021-02-011")), 
    ylim=c(0.0001,0.99), expand=c(0,0)) +
  scale_fill_manual("region", values=c("red", "blue")) + # reg_cols) +
  scale_size_continuous("total number\nof sequences\nper week", trans="log2", 
                        range=c(0.01, 4), limits=c(100,1600), breaks=c(100,200,400,800,1600))  +
  geom_point(data=data_denmark_uk, 
              aes(x=sample_date, y=prop, colour=REGION, size=TOTAL, fill=NULL), pch=16, alpha=I(0.4)) +
  scale_color_manual("region", values=reg_cols) +
  theme(legend.direction = "vertical", legend.box = "vertical", legend.position="right", 
        legend.key.size=unit(0.45, "cm"))
plot_denmUK_country



#  CALCULATE GROWTH RATES & SELECTIVE ADVANTAGE OF VOC IN UK & DK

# POSTHOC TEST TO TEST FOR DIFFERENCES IN GROWTH RATES IN THE VOC BETWEEN UK & DENMARK

# PER COUNTRY, BASED ON MODEL bGLMMfit_denmUK0 WITH HOMOGENEOUS SLOPES PER REGION (but not country)
tukey_VOC_DENMUK = as.data.frame(contrast(emtrends(bGLMMfit_denmUK0, ~ COUNTRY, 
                                          var="sample_date_num"),method="revpairwise"))[,-4]
colnames(tukey_VOC_DENMUK)[2] = "diff_logistic_growth_rate"
# UK estimated to have slightly faster spread than UK (p=0.03)
tukey_VOC_DENMUK
#       contrast diff_logistic_growth_rate          SE  z.ratio     p.value
# 1 UK - DENMARK                0.02541104 0.009757873 2.604158 0.009210037

# recalculate logistic growth rates to differences in R value in both countries
bGLMM_VOC_growthrates_DENMUK = as.data.frame(emtrends(bGLMMfit_denmUK0, ~ COUNTRY, var="sample_date_num"))[,-c(3,4)] 
colnames(bGLMM_VOC_growthrates_DENMUK)[2] = "logistic_growth_rate"

# for our two sets of default estimated growth rates these logistic growth rates 
# translate into a proportional increase in R value of
bGLMM_VOC_growthrates_DENMUK = M.from.delta_r_df(bGLMM_VOC_growthrates_DENMUK)
bGLMM_VOC_growthrates_DENMUK
#   COUNTRY logistic_growth_rate  asymp.LCL asymp.UCL       M1   M1.LCL   M1.UCL       M2   M2.LCL   M2.UCL
# 1 DENMARK           0.08390644 0.06644913 0.1013637 1.586429 1.441191 1.746302 1.352647 1.270254 1.440384
# 2      UK           0.10931748 0.10131903 0.1173159 1.824391 1.745873 1.906440 1.482223 1.440152 1.525523
table2csv(bGLMM_VOC_growthrates_DENMUK, file=".\\multinomial_logistic_fits\\tables\\model 3b_VOCvsall_bGLMM_growthrates_DENMUK_country but not region specific slopes.csv")

