# MULTINOMIAL AND LOGISTIC FITS TO DETERMINE GROWTH RATE & 
# COMPETITIVE ADVANTAGE OF VOC 202012/01 COMPARED TO OTHER 
# CIRCULATING SARS-CoV2 VARIANTS ACROSS DIFFERENT REGIONS IN THE UK 
# AS WELL AS IN DENMARK, BELGIUM, SWITZERLAND & THE USA 
# AND IMPLIED EXPECTED INCREASE IN R VALUES, 
# ASSUMING UNALTERED GENERATION TIME

# T. Wenseleers, last updated: 14 Febr. 2021

library(lme4)
require(devtools)
install_version("emmeans", version = "1.5.3", repos = "http://cran.us.r-project.org") # we use version 1.5.3, since version 1.5.4 accidentically introduced a bug in its treatment of nnet::multinom models
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
library(mclogit)
library(tidyr)
library(ISOweek)
library(readr)



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





# 1. LOAD & PREPROCESS COG-UK SEQUENCING DATA ####
# desired order of factor levels
levels_nhs_name = c("South East","London","East of England",
                    "South West","Midlands","North East and Yorkshire",
                    "Scotland","North West","Wales")


# read in data
# data = read.csv(".\\data\\cog_metadata_microreact_public-2020-12-22-annotated.csv")
data = read.csv(".\\data\\cog_metadata_microreact_public-2021-01-11-annotated.csv")
data$sample_date = as.Date(data$sample_date)
data$sample_date_num = as.numeric(data$sample_date)
# we do not consider data from Northern Ireland due to low nr of sequences there & absence of new VOC
data = data[data$nhs_name!="Northern Ireland",] 
data$nhs_name = factor(data$nhs_name, levels=levels_nhs_name) # NHS region
data$lad = as.factor(data$lad) # local authority district
unique(data[data$n501y == "Y","lineage"]) # lineages where at least some samples have n501y mutation
# "B.1.1.7"   "B.1.1.70"  "B.1.351"   "B.1.160"   "B.1.1"     "B.1.1.189" "B.1"       "B.1.177"   "B.1.1.136"
# table(data$lineage, data$n501y, data$del_21765_6)
data$isoweek = date2ISOweek(data$sample_date) # isoweek
data$yearweek = sub("[-][^-]+$", "", data$isoweek)
# we slightly recode some lineages into specific variants
data$variant = data$lineage
# B.1.177 & descendant lineages = B.1.177 in nextstrain
data$variant[grepl("B.1.177",data$lineage,fixed=T)] = "B.1.177" 
# VOC 202012/01 = VOC 202012/01 = B.1.1.7 + n501y == "Y" + del_21765_6=="del"
sum(data$lineage=="B.1.1.7"&data$n501y =="Y"&data$del_21765_6=="del") # 13352
data$variant[data$lineage=="B.1.1.7"&
             data$n501y =="Y"&
             data$del_21765_6=="del"] = "VOC 202012/01" # VARIANT OF CONCERN, del_21765_6=deletion Δ69/Δ70
table(data[data$lineage=="B.1.1.7","n501y"]=="Y", data[data$lineage=="B.1.1.7","del_21765_6"]=="del")
#       FALSE  TRUE
# FALSE    10    94
# TRUE    141 13352
# data$variant[data$lineage=="B.1.1.70"&data$n501y == "Y"]="501Y.WALES" # Welsh variant, no official name yet


# aggregate data by year, week and check which variant lineages reached at least 15% in some week    
data_agbyweek = as.data.frame(table(data$yearweek, data$variant))
colnames(data_agbyweek) = c("yearweek", "variant", "count")
data_agbyweek_sum = aggregate(count ~ yearweek, data=data_agbyweek, sum)
data_agbyweek$total = data_agbyweek_sum$count[match(data_agbyweek$yearweek, data_agbyweek_sum$yearweek)]
sum(data_agbyweek[data_agbyweek$variant=="VOC 202012/01","total"]) == nrow(data) # correct
data_agbyweek$sample_date = ISOweek2date(paste0(data_agbyweek$yearweek,"-3"))
data_agbyweek$variant = as.factor(data_agbyweek$variant)
data_agbyweek$prop = data_agbyweek$count/data_agbyweek$total
maxweeklyprop_lineages = aggregate(prop ~ variant, data=data_agbyweek, max)
# we select the 9 most abundant lineages
# (the ones at one point in time reached at least a relative abundance of 13%)
selectedlineages = maxweeklyprop_lineages$variant[order(maxweeklyprop_lineages$prop, decreasing=T)[1:9]] # as.character(maxweeklyprop_lineages$variant[maxweeklyprop_lineages$prop>0.136])
as.character(selectedlineages)
# "B"             "B.1.1"         "B.1.177"       "VOC 202012/01" "B.1.98"        "B.1.1.315"     "B.40"          "B.1.1.1"      
# "B.1" 

levels_variants = c("B","B.1.98","B.40","B.1","B.1.1",
                    "B.1.1.1", "B.1.1.315", "B.1.177",
                    "VOC 202012/01",  
                    "minority variants")
# colours to use for variant lineages
n = length(levels_variants)
lineage_cols = hcl(h = seq(15, 375, length = n + 1), l = 65, c = 200)[1:n]
lineage_cols[which(levels_variants=="minority variants")] = "grey75"


# we recode all other 410 minor lineages as a single category "minority variants"
data$variant[!data$variant %in% selectedlineages] = "minority variants"
length(unique(data$lineage[data$variant=="minority variants"])) # 410 lineages

data$variant = factor(data$variant, levels=levels_variants)

sum(is.na(data$variant)) # 0, OK


# aggregated data to make Muller plots of raw data
# aggregated by week for selected variant lineages
data_agbyweek = as.data.frame(table(data$yearweek, data$variant))
colnames(data_agbyweek) = c("yearweek", "variant", "count")
data_agbyweek_sum = aggregate(count ~ yearweek, data=data_agbyweek, sum)
data_agbyweek$total = data_agbyweek_sum$count[match(data_agbyweek$yearweek, data_agbyweek_sum$yearweek)]
sum(data_agbyweek[data_agbyweek$variant=="VOC 202012/01","total"]) == nrow(data) # correct
data_agbyweek$sample_date = ISOweek2date(paste0(data_agbyweek$yearweek,"-3"))
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
data_agbyweekregion = as.data.frame(table(data$yearweek, data$nhs_name, data$variant))
colnames(data_agbyweekregion) = c("yearweek", "nhs_name", "variant", "count")
data_agbyweekregion_sum = aggregate(count ~ yearweek + nhs_name, data=data_agbyweekregion, sum)
data_agbyweekregion$total = data_agbyweekregion_sum$count[match(interaction(data_agbyweekregion$yearweek,data_agbyweekregion$nhs_name), 
                                                          interaction(data_agbyweekregion_sum$yearweek,data_agbyweekregion_sum$nhs_name))]
sum(data_agbyweekregion[data_agbyweekregion$variant=="VOC 202012/01","total"]) == nrow(data) # correct
data_agbyweekregion$sample_date = ISOweek2date(paste0(data_agbyweekregion$yearweek,"-3"))
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
write.csv(data_agbydayregion, 
          file=".\\multinomial_logistic_fits\\data\\data_agbydayregion.csv", row.names=F)


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
write.csv(data_agbydayregionlad, 
          file=".\\multinomial_logistic_fits\\data\\data_agbydayregionlad.csv", row.names=F)


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
sum(data_agbydayregion_long$TOTAL) # 1472810 = 10 times bigger than nrow(data), correct
nrow(data) # 147281
data_agbydayregion_long = data_agbydayregion_long[data_agbydayregion_long$TOTAL!=0,]
data_agbydayregion_long$obs = factor(1:nrow(data_agbydayregion_long)) # for observation-level random factor
data_agbydayregion_long$variant = factor(data_agbydayregion_long$variant,
                                         levels=levels_variants)
data_agbydayregion_long$variant = relevel(data_agbydayregion_long$variant,
                                             ref="minority variants")
nrow(data_agbydayregion_long) # 87250

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
sum(data_agbydayregionlad_long$TOTAL) # 1472810 = 10 times bigger than nrow(data), correct
nrow(data) # 147281
data_agbydayregionlad_long = data_agbydayregionlad_long[data_agbydayregionlad_long$TOTAL!=0,]
data_agbydayregionlad_long$obs = factor(1:nrow(data_agbydayregionlad_long)) # for observation-level random factor
data_agbydayregionlad_long$variant = factor(data_agbydayregionlad_long$variant,
                                                    levels=levels_variants)
data_agbydayregionlad_long$variant = relevel(data_agbydayregionlad_long$variant,
                                                     ref="minority variants")
nrow(data_agbydayregionlad_long) # 315650



# MULLER PLOTS OF RAW DATA ####

# Muller plots showing spread of main variant lineages
n = length(levels(data$variant))
lineage_cols = hcl(h = seq(15, 375, length = n + 1), l = 65, c = 200)[1:n]
lineage_cols[which(levels(data_agbyday$variant)=="minority variants")] = "grey75"

# Muller plots using daily & weekly aggregated data (overall across all regions)
range(data$sample_date) # "2020-02-05" "2021-01-06"
ggplot(data=data_agbyday, aes(x=sample_date, y=count, group=variant)) + 
    geom_area(aes(lwd=I(1.2), colour=NULL, fill=variant), position="fill") +
    scale_fill_manual("", values=lineage_cols) +
    scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
                       labels=substring(months(as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01"))),1,1),
                       limits=as.Date(c("2020-03-01",NA)), expand=c(0,0)) +
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
                       limits=as.Date(c("2020-03-01",NA)), expand=c(0,0)) +
    guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
    theme_hc() + theme(legend.position="bottom",  
                       axis.title.x=element_blank()) + 
    # labs(title = "MAIN SARS-CoV2 VARIANT LINEAGES IN THE UK") +
    ylab("Relative abundance")
muller_raw0
ggsave(file=".\\multinomial_logistic_fits\\plots\\muller plot lineages_overall_raw data.png", width=7, height=5)
ggsave(file=".\\multinomial_logistic_fits\\plots\\muller plot lineages_overall_raw data.pdf", width=7, height=5)
saveRDS(muller_raw0, file = ".\\multinomial_logistic_fits\\plots\\muller plot lineages_overall_raw data.rds")
graph2ppt(file=".\\multinomial_logistic_fits\\plots\\muller plot lineages_overall_raw data.pptx", width=7, height=5)

muller_raw = ggplot(data=data_agbyweekregion, aes(x=sample_date, y=count, group=variant)) + 
    geom_area(aes(lwd=I(1.2), colour=NULL, fill=variant), position="fill") +
    facet_wrap(~nhs_name) +
    scale_fill_manual("", values=lineage_cols) +
    scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
                       labels=substring(months(as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01"))),1,1),
                       limits=as.Date(c("2020-03-01",NA)), expand=c(0,0)) +
    guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
    theme_hc() + theme(legend.position="bottom",  
                       axis.title.x=element_blank()) + 
    # labs(title = "MAIN SARS-CoV2 VARIANT LINEAGES IN THE UK") +
    ylab("Relative abundance")
muller_raw
ggsave(file=".\\multinomial_logistic_fits\\plots\\FigS2_muller plot lineages by region_raw data.png", width=7, height=5)
ggsave(file=".\\multinomial_logistic_fits\\plots\\FigS2_muller plot lineages by region_raw data.pdf", width=7, height=5)
saveRDS(muller_raw, file = ".\\multinomial_logistic_fits\\plots\\FigS2_muller plot lineages by region_raw data.rds")
graph2ppt(file=".\\multinomial_logistic_fits\\plots\\FigS2_muller plot lineage by region_raw data.pptx", width=7, height=5)




# 2. MULTINOMIAL FITS ON COG-UK SEQUENCING DATA BY REGION ####

# 2.1 SEPARATE-SLOPES MULTINOMIAL SPLINE FIT ####

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
# mfit1  90 460202.4
# mfit2 162 435085.6
# mfit3  99 302278.0
# mfit4 243 299882.3


# we continue with mfit4 as this was a realistic fit, extrapolates in a stable way & 
# has the best BIC of these models

mfit = mfit4 
mfit


# GROWTH RATE CONTRASTS

# the mean difference in growth rates between the variants and corresponding 
# expected multiplicative effects on the R value for the multinomial separate slopes spline model mfit4
# on average across regions, evaluated for the relevant timeframes where the variants started invading

# the growth contrasts for the VOC - minority variants and for the VOC - B.1.177 comparisons we calculate
# for the period from 1 Nov 2020 until the 6th of January 2021
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
#                              contrast           SE  t.ratio       p.value    delta_r delta_r.lower.CL delta_r.upper.CL       M1   M1.LCL
# 1 (VOC 202012/01) - minority variants 0.0018715734 59.88540  0.000000e+00 0.11207992       0.10786819       0.11629165 1.852321 1.809906
# 2           (VOC 202012/01) - B.1.177 0.0018073550 57.68033  0.000000e+00 0.10424884       0.10018163       0.10831606 1.774234 1.734985
# 3         B.1.177 - minority variants 0.0008968844 62.21600 2.767138e-151 0.05580055       0.05403389       0.05756721 1.359209 1.346066
# M1.UCL       M2   M2.LCL   M2.UCL
# 1 1.89573 1.497037 1.474510 1.519908
# 2 1.81437 1.455422 1.434267 1.476889
# 3 1.37248 1.222481 1.214730 1.230280

table2csv(mfit_contrasts, file=".\\multinomial_logistic_fits\\tables\\model 1a_mfit4_multinomial spline fit_growthrates_UK_heter slopes.csv")





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
    annotate("rect", xmin=max(data$sample_date)+1, xmax=max(data$sample_date_num)+extrapolate, ymin=0, ymax=1, alpha=0.4, fill="white") + # extrapolated part
    scale_fill_manual("", values=lineage_cols) +
    scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
                       labels=substring(months(as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01"))),1,1),
                       expand=c(0,0)) +
    guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
    theme_hc() + theme(legend.position="bottom", 
                       axis.title.x=element_blank()) + 
    # labs(title = "MAIN SARS-CoV2 VARIANT LINEAGES IN THE UK") +
    ylab("Relative abundance") + coord_cartesian(xlim=c(as.Date("2020-03-01"),as.Date("2021-02-28")))
muller_mfit
ggsave(".\\multinomial_logistic_fits\\plots\\Fig2C_model 1a_plot multinomial spline fit_muller plot fit.png", width=7, height=5)
ggsave(".\\multinomial_logistic_fits\\plots\\Fig2C_model 1a_plot multinomial spline fit_muller plot fit.pdf", width=7, height=5)
saveRDS(muller_mfit, file = ".\\multinomial_logistic_fits\\plots\\Fig2C_model 1a_plot multinomial spline fit_muller plot fit.rds")
graph2ppt(file=".\\multinomial_logistic_fits\\plots\\Fig2C_model 1a_plot multinomial spline fit_muller plot fit.pptx", width=7, height=5)

# plot with raw data and fit combined:
ggarrange(muller_raw +
              coord_cartesian(xlim=c(as.Date("2020-03-01"),
                                     as.Date("2021-02-28")), 
                              expand=c(0,0)) + labs(title = "") , 
          muller_mfit+labs(title = "")+coord_cartesian(xlim=c(as.Date("2020-03-01"),
                                                             as.Date("2021-02-28")), 
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
  facet_wrap(~variant, ncol=1, scale="free") +
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
  scale_size_continuous("total number\nof sequences\nper week", trans="sqrt", 
                        range=c(0.01, 4), limits=c(1,max(data_agbyweekregion$total)), breaks=c(10,100,1000)) +
  theme(legend.direction = "vertical", legend.box = "vertical", legend.position="right", 
        legend.key.size=unit(0.45, "cm"))
plotmultinom2vars

ggsave(file=".\\multinomial_logistic_fits\\plots\\FigS3_model 1a_plot multinomial spline fit_growth VOC and B1177.png", width=6, height=6)
ggsave(file=".\\multinomial_logistic_fits\\plots\\FigS3_model 1a_plot multinomial spline fit_growth VOC and B1177.pdf", width=6, height=6)
saveRDS(plotmultinom2vars, file = ".\\multinomial_logistic_fits\\plots\\FigS3_model 1a_plot multinomial spline fit_growth VOC and B1177.rds")
graph2ppt(file=".\\multinomial_logistic_fits\\plots\\FigS3_model 1a_plot multinomial spline fit_growth VOC and B1177.pptx", width=6, height=6)



# 2.2 COMMON-SLOPES MULTINOMIAL MIXED MODEL FIT ####
# MULTINOMIAL MIXED MODEL FITS, INCORPORATING RANDOM EFFECT 1|lad

# To take into account spatial dependencies we also fit a mixed baseline category
# multinomial model with a random intercept for lad (LTLA) using the mclogit::mblogit
# and we also use the inbuilt option to take into account overdispersion

# mixed baseline category multinomial model fit using mclogit::mblogit : 
data_refminor = data
data_refminor$variant = relevel(data_refminor$variant, ref="minority variants")
data_refB1177 = data
data_refB1177$variant = relevel(data_refB1177$variant, ref="B.1.177")


# models with lad included as random intercept and using either minority variants or B.1.177 as reference group
# so that we can later interpret the multinomial model coefficients in function of sample_date_num as the contrasts in growth rate
# with this reference group
# a model with common slopes across nhs regions has the best BIC
set_treatment_contrasts()
mblogitfit1_refminor = mblogit(variant ~ nhs_name+sample_date_num, data=data_refminor,
                               random = ~1|lad, dispersion=TRUE)
mblogitfit1_refB1177 = mblogit(variant ~ nhs_name+sample_date_num, data=data_refB1177, 
                               random = ~1|lad, dispersion=TRUE)

# saveRDS(mblogitfit1_refminor, file = ".\\multinomial_logistic_fits\\fits\\mblogitfit1_refminor_model1b.rds")
# saveRDS(mblogitfit1_refB1177, file = ".\\multinomial_logistic_fits\\fits\\mblogitfit1_refB1177_model1b.rds")
# or to directly load previously fitted models
mblogitfit1_refminor = readRDS(file = ".\\multinomial_logistic_fits\\fits\\mblogitfit1_refminor_model1b.rds")
mblogitfit1_refB1177 = readRDS(file = ".\\multinomial_logistic_fits\\fits\\mblogitfit1_refB1177_model1b.rds")


# plot model predictions of best fitting common slopes across regions multinomial mixed model
extrapolate = 60 # nr of days to show extrapolations for
# t = table(data$lad[data$nhs_name==levels_nhs_name[1]])
# sel_lad = names(t[which.max(t)]) # PS with condition=FALSE this selection has no effect
newdat = expand.grid(sample_date_num =
                       seq(min(data$sample_date_num),max(data$sample_date_num)+extrapolate,1),
                     nhs_name = levels_nhs_name, 
                     lad = data$lad[1]) # PS this selection has no effect with conditional=FALSE in predict()
# see https://www.elff.eu/software/mclogit/manual/predict/
preds_mblogit1_wide = data.frame(predict(mblogitfit1_refminor, 
                                         newdata = newdat,
                                         type="response", se.fit=FALSE, conditional=FALSE), check.names=F)
preds_mblogit1_long = gather(preds_mblogit1_wide, variant, prob, levels_variants, factor_key=TRUE)
preds_mblogit1_long$sample_date_num = newdat$sample_date_num
preds_mblogit1_long$sample_date = as.Date(preds_mblogit1_long$sample_date_num, origin="1970-01-01")
preds_mblogit1_long$nhs_name = factor(newdat$nhs_name, levels=levels_nhs_name)
muller_mblogit1 = ggplot(data=preds_mblogit1_long, 
                         aes(x=sample_date, y=prob, group=variant)) + 
  facet_wrap(~nhs_name) +
  geom_area(aes(lwd=I(1.2), colour=NULL, fill=variant), position="fill") +
  annotate("rect", xmin=max(data$sample_date)+1, xmax=max(data$sample_date)+extrapolate, ymin=0, ymax=1, alpha=0.4, fill="white") + # extrapolated part
  scale_fill_manual("", values=lineage_cols) +
  scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
                     labels=substring(months(as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01"))),1,1),
                     limits=as.Date(c("2020-03-01",NA)), expand=c(0,0)) +
  guides(color = guide_legend(reverse=F, nrow=2, byrow=T), fill = guide_legend(reverse=F, nrow=2, byrow=T)) +
  theme_hc() + theme(legend.position="bottom", 
                     axis.title.x=element_blank()) + 
  # labs(title = "MAIN SARS-CoV2 VARIANT LINEAGES IN THE UK") +
  ylab("Relative abundance") +
  coord_cartesian(xlim=c(as.Date("2020-03-01"),
                         as.Date("2021-02-28")), 
                  expand=c(0,0)) + labs(title = "")
muller_mblogit1
ggsave(".\\multinomial_logistic_fits\\plots\\FigS4_model 1b_plot multinomial multinomial mixed model_muller plot fit.png", width=7, height=5)
ggsave(".\\multinomial_logistic_fits\\plots\\FigS4_model 1b_plot multinomial multinomial mixed model_muller plot fit.pdf", width=7, height=5)
saveRDS(muller_mblogit1, file = ".\\multinomial_logistic_fits\\plots\\FigS4_model 1b_plot multinomial mixed model_muller plot fit.rds")
graph2ppt(file=".\\multinomial_logistic_fits\\plots\\FigS4_model 1b_plot multinomial multinomial mixed model_muller plot fit.pptx", width=7, height=5)


# since mblogit was not supported by emmeans at the time of writing, we estimate the mean difference in growth rates
# with the "minority level" reference level and the expected multiplicative effect on the R
# value directly from the model coefficients (which for this model would be identical)
mcoefs = data.frame(coef=coef(mblogitfit1_refminor),confint(mblogitfit1_refminor),check.names=F)
mcoefs = mcoefs[grepl("VOC|B\\.1\\.177",rownames(mcoefs),fixed=F),]
mcoefs = mcoefs[grepl("sample_date_num",rownames(mcoefs)),]
colnames(mcoefs)[1] = "delta_r"
mcoefs = data.frame(contrast = paste0(gsub("~sample_date_num","",rownames(mcoefs), fixed=T)," - minority variants"),
                    mcoefs, check.names=F)
rownames(mcoefs) = NULL
mcoefs2 = data.frame(coef=coef(mblogitfit1_refB1177),confint(mblogitfit1_refB1177),check.names=F)
mcoefs2 = mcoefs2[grepl("VOC",rownames(mcoefs2),fixed=F),]
mcoefs2 = mcoefs2[grepl("sample_date_num",rownames(mcoefs2)),]
colnames(mcoefs2)[1] = "delta_r"
mcoefs2 = data.frame(contrast = paste0(gsub("~sample_date_num","",rownames(mcoefs2), fixed=T)," - B.1.177"),
                     mcoefs2, check.names=F)
rownames(mcoefs2) = NULL
mblogitfit_contrasts = M.from.delta_r_df(rbind(mcoefs2,mcoefs))
mblogitfit_contrasts
#                            contrast    delta_r      2.5 %     97.5 %       M1   M1.LCL   M1.UCL       M2   M2.LCL   M2.UCL
# 1           VOC 202012/01 - B.1.177 0.09306863 0.09117032 0.09496694 1.668421 1.651092 1.685931 1.398006 1.388485 1.407593
# 2       B.1.177 - minority variants 0.03040260 0.02997586 0.03082933 1.182008 1.179237 1.184785 1.115664 1.113951 1.117379
# 3 VOC 202012/01 - minority variants 0.12346639 0.12152546 0.12540732 1.972033 1.951094 1.993198 1.559677 1.548817 1.570614
table2csv(mblogitfit_contrasts, file=".\\multinomial_logistic_fits\\tables\\model1b_mblogitfit_multinomial mixed model fit_growthrates_UK_homog slopes.csv")

# # PS: this doesn't work at the time of writing since mblogit is not supported by the emmeans package
# this may be fixed in a subsequent version
# mblogitfit_emtrends = emtrends(mblogitfit1_refminor, revpairwise ~ variant, var="sample_date_num", 
#                           mode="latent", adjust="Tukey", 
#                           at=list(sample_date_num=as.numeric(max(data$sample_date)),
#                                   variant=c("VOC 202012/01","B.1.177","minority variants"))) 
# mblogitfit_contrasts = data.frame(as.data.frame(mblogitfit_emtrends$contrasts),
#                              as.data.frame(confint(mblogitfit_emtrends$contrasts))[,c("lower.CL","upper.CL")])
# colnames(mblogitfit_contrasts)[which(colnames(mblogitfit_contrasts) %in% c("estimate","lower.CL","upper.CL"))] = 
#     c("delta_r","delta_r.lower.CL","delta_r.upper.CL")
# mblogitfit_contrasts = data.frame(mblogitfit_contrasts[,c("contrast","SE","t.ratio","p.value")],
#                              M.from.delta_r_df(mblogitfit_contrasts[,c("delta_r","delta_r.lower.CL","delta_r.upper.CL")]))
# mblogitfit_contrasts






# 3. BINOMIAL GLMM FITS ON COG-UK SEQUENCING DATA TO COMPARE GROWTH RATE ADVANTAGE OF STRAINS VOC 202012/01 AND B.1.177 ACROSS DIFFERENT REGIONS ####

# 3.1 BINOMIAL MIXED MODEL FIT FOR VARIANT VOC 202012/01 ####

# 3.1.1 BINOMIAL MIXED MODEL FIT FOR VARIANT VOC 202012/01 vs all other variants ####

# mixed binomial GLMM with random intercept for lad (LTLA) & 
# observation-level random effect to take into account overdispersion
# and fixed effects ns_name and date, with or without interaction effect nhs_name x date

start_date = as.Date("2020-08-01")
stop_date = max(data_agbydayregionlad$sample_date) # 6 Jan 2021
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
glmersettings = glmerControl(optimizer="optimx", optCtrl=list(method="nlminb")) # PS : to try all optimizer run all_fit(fit)
glmersettings2 = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=1E4)) # PS : to try all optimizer run all_fit(fit)
bGLMMfit1 = glmer(cbind(count, total-count) ~  (1|lad/obs) + 
                       nhs_name+scale(sample_date_num), 
                     family=binomial(logit), data=data_subs, control=glmersettings2)
bGLMMfit2 = glmer(cbind(count, total-count) ~  (1|lad/obs) + 
                       nhs_name*scale(sample_date_num), 
                     family=binomial(logit), data=data_subs, control=glmersettings2)

# saveRDS(bGLMMfit1, file = ".\\multinomial_logistic_fits\\fits\\bGLMMfit_VOCvsall_fit1_homog slopes.rds")
# saveRDS(bGLMMfit2, file = ".\\multinomial_logistic_fits\\fits\\bGLMMfit_VOCvsall_fit2_heter slopes_model 2a.rds")
# or to directly load previously fitted models
bGLMMfit1 = readRDS(file = ".\\multinomial_logistic_fits\\fits\\bGLMMfit_VOCvsall_fit1_homog slopes.rds")
bGLMMfit2 = readRDS(file = ".\\multinomial_logistic_fits\\fits\\bGLMMfit_VOCvsall_fit1_heter slopes_model 2a.rds")

# check BIC values
BIC(bGLMMfit1, bGLMMfit2)
#           df      BIC
# bGLMMfit1 12 7554.800
# bGLMMfit2 20 7524.907


# bGLMMfit2 has the best BIC (heterogeneous slopes across regions)


# PLOT MODEL FIT

# of model bGLMMfit2 (heterogeneous slopes, model 2 in table S1)
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
    xlim=c(as.Date("2020-07-01"),as.Date("2021-02-28")), 
    ylim=c(0.0001,0.999001), expand=c(0,0)) +
  scale_color_discrete("", h=c(0, 280), c=200) +
  scale_fill_discrete("", h=c(0, 280), c=200) +
  geom_point(data=data_agbyweekregion[data_agbyweekregion$variant=="VOC 202012/01",], 
             aes(x=sample_date, y=prop, colour=nhs_name, size=total), alpha=I(0.5)) +
  scale_size_continuous("total number\nof sequences\nper week", trans="sqrt", 
                        range=c(0.01, 4), limits=c(1,max(data_agbyweekregion$total)), breaks=c(10,100,1000)) +
  guides(fill=FALSE) + guides(colour=FALSE) + theme(legend.position = "right")
# geom_point(data=data_agbydayregion[data_agbydayregion$variant=="VOC 202012/01",], 
#         aes(x=sample_date, y=prop, colour=nhs_name, size=total), alpha=I(0.4)) +
# theme(legend.position = "none")
plot_bGLMMVOC_het

saveRDS(plot_bGLMMVOC_het, file = ".\\multinomial_logistic_fits\\plots\\model2a_plot VOCvsall_fit bGLMM_heter slopes.rds")
graph2ppt(file=".\\multinomial_logistic_fits\\plots\\model2a_plot VOCvsall_fit bGLMM_heter slopes.pptx", width=8, height=6)
ggsave(file=".\\multinomial_logistic_fits\\plots\\model2a_plot VOCvsall_fit bGLMM_heter slopes.png", width=8, height=6)
ggsave(file=".\\multinomial_logistic_fits\\plots\\model2a_plot VOCvsall_fit bGLMM_heter slopes.pdf", width=8, height=6)




#  CALCULATE GROWTH RATE & TRANSMISSION ADVANTAGE OF VOC IN DIFFERENT REGIONS

# growth rate of variant = slope of logistic regression and corresponding expected multiplicative increase in R value
# using best fitting model bGLMMfit2 (heterogeneous slopes) 

bGLMM_VOC_growthrates = as.data.frame(emtrends(bGLMMfit2, ~ nhs_name, var="sample_date_num"))[,-c(3,4)] 
# sample_date_num.trend = 
# logistic growth rate of VOC 202012/01 = Malthusian growth rate VOC - Malthusian growth rate all other variants
colnames(bGLMM_VOC_growthrates)[2] = "logistic_growth_rate"
bGLMM_VOC_growthrates = M.from.delta_r_df(bGLMM_VOC_growthrates)
bGLMM_VOC_growthrates
# nhs_name logistic_growth_rate  asymp.LCL  asymp.UCL       M1   M1.LCL   M1.UCL       M2   M2.LCL   M2.UCL
# 1               South East           0.09807042 0.09253065 0.10361018 1.714956 1.663491 1.768012 1.423407 1.395301 1.452079
# 2                   London           0.10466452 0.09718873 0.11214032 1.778295 1.706660 1.852937 1.457602 1.418896 1.497362
# 3          East of England           0.09223946 0.08660012 0.09787879 1.660829 1.610107 1.713149 1.393839 1.365827 1.422426
# 4               South West           0.10208649 0.08707703 0.11709595 1.753258 1.614336 1.904135 1.444136 1.368174 1.524315
# 5                 Midlands           0.12410097 0.11412366 0.13407829 1.978928 1.873260 2.090557 1.563245 1.508092 1.620414
# 6 North East and Yorkshire           0.12646821 0.11565810 0.13727833 2.004862 1.889136 2.127677 1.576624 1.516446 1.639190
# 7                 Scotland           0.14332825 0.11875351 0.16790300 2.199666 1.921573 2.518004 1.675282 1.533439 1.830247
# 8               North West           0.14496214 0.13031575 0.15960852 2.219522 2.047740 2.405714 1.685165 1.598614 1.776403
# 9                    Wales           0.11395643 0.10224428 0.12566857 1.871538 1.754780 1.996064 1.507184 1.444957 1.572091
table2csv(bGLMM_VOC_growthrates, file=".\\multinomial_logistic_fits\\tables\\model 2a_VOCvsall_bGLMM_VOC_growthrates_UK_by region_heter slopes.csv")

# marginal slopes across all regions
bGLMM_VOC_growthrates_avg = as.data.frame(emtrends(bGLMMfit2, ~ 1, var="sample_date_num"))[,-c(3,4)] 
colnames(bGLMM_VOC_growthrates_avg)[2] = "logistic_growth_rate"
bGLMM_VOC_growthrates_avg = M.from.delta_r_df(bGLMM_VOC_growthrates_avg)
bGLMM_VOC_growthrates_avg
# 1         logistic_growth_rate asymp.LCL asymp.UCL       M1   M1.LCL   M1.UCL       M2   M2.LCL   M2.UCL
# 1 overall             0.116653 0.1122397 0.1210663 1.899501 1.85395 1.946172 1.521887 1.497898 1.546259
table2csv(bGLMM_VOC_growthrates_avg, file=".\\multinomial_logistic_fits\\tables\\model 2a_VOCvsall_bGLMM_VOC_avggrowthrate_UK_heter slopes.csv")


# TUKEY POSTHOC TESTS TO TEST FOR DIFFERENCES IN GROWTH RATE OF THE VOC ACROSS DIFFERENT NHS REGIONS

tukey_VOC = as.data.frame(emtrends(bGLMMfit2, pairwise ~ nhs_name, var="sample_date_num", adjust="tukey")$contrasts)[,-4]
colnames(tukey_VOC)[2] = "diff_logistic_growth_rate"
tukey_VOC
#                                      contrast diff_logistic_growth_rate          SE    z.ratio      p.value
# 1                         South East - London              -0.006594109 0.004702554 -1.4022397 8.974025e-01
# 2                South East - East of England               0.005830960 0.003963607  1.4711249 8.691112e-01
# 3                     South East - South West              -0.004016072 0.008131069 -0.4939169 9.999104e-01
# 4                       South East - Midlands              -0.026030559 0.005786990 -4.4981171 2.359861e-04
# 5       South East - North East and Yorkshire              -0.028397798 0.006147364 -4.6195084 1.334803e-04
# 6                       South East - Scotland              -0.045257839 0.012829825 -3.5275493 1.250741e-02
# 7                     South East - North West              -0.046891719 0.007972765 -5.8814877 1.458678e-07
# 8                          South East - Wales              -0.015886011 0.006575921 -2.4157851 2.750782e-01
# 9                    London - East of England               0.012425069 0.004730412  2.6266357 1.754028e-01
# 10                        London - South West               0.002578037 0.008533622  0.3021034 9.999980e-01
# 11                          London - Midlands              -0.019436450 0.006335142 -3.0680369 5.522961e-02
# 12          London - North East and Yorkshire              -0.021803689 0.006670884 -3.2684857 2.988393e-02
# 13                          London - Scotland              -0.038663730 0.013091185 -2.9534172 7.657275e-02
# 14                        London - North West              -0.040297611 0.008378280 -4.8097715 5.296858e-05
# 15                             London - Wales              -0.009291903 0.007065418 -1.3151243 9.271881e-01
# 16               East of England - South West              -0.009847032 0.008149599 -1.2082843 9.549697e-01
# 17                 East of England - Midlands              -0.031861519 0.005810572 -5.4833703 1.490709e-06
# 18 East of England - North East and Yorkshire              -0.034228758 0.006170381 -5.5472682 1.037504e-06
# 19                 East of England - Scotland              -0.051088799 0.012843050 -3.9779334 2.264178e-03
# 20               East of England - North West              -0.052722680 0.007991042 -6.5977223 1.502050e-09
# 21                    East of England - Wales              -0.021716972 0.006597944 -3.2914754 2.775706e-02
# 22                      South West - Midlands              -0.022014487 0.009174244 -2.3995969 2.839346e-01
# 23      South West - North East and Yorkshire              -0.024381726 0.009403722 -2.5927740 1.894170e-01
# 24                      South West - Scotland              -0.041241767 0.014665319 -2.8121970 1.116605e-01
# 25                    South West - North West              -0.042875647 0.010686293 -4.0122097 1.968868e-03
# 26                         South West - Wales              -0.011869939 0.009688986 -1.2250962 9.512014e-01
# 27        Midlands - North East and Yorkshire              -0.002367239 0.007474135 -0.3167241 9.999971e-01
# 28                        Midlands - Scotland              -0.019227280 0.013515849 -1.4225728 8.894884e-01
# 29                      Midlands - North West              -0.020861160 0.009030644 -2.3100414 3.358090e-01
# 30                           Midlands - Wales               0.010144548 0.007828062  1.2959207 9.328712e-01
# 31        North East and Yorkshire - Scotland              -0.016860041 0.013669219 -1.2334312 9.492525e-01
# 32      North East and Yorkshire - North West              -0.018493921 0.009270447 -1.9949331 5.468852e-01
# 33           North East and Yorkshire - Wales               0.012511786 0.008098672  1.5449184 8.342143e-01
# 34                      Scotland - North West              -0.001633881 0.014583071 -0.1120395 1.000000e+00
# 35                           Scotland - Wales               0.029371827 0.013867866  2.1179774 4.608067e-01
# 36                         North West - Wales               0.031005708 0.009555781  3.2447068 3.223206e-02

tukey_VOC[tukey_VOC$p.value<0.05,]
#                                 contrast diff_logistic_growth_rate          SE   z.ratio      p.value
# 4                       South East - Midlands               -0.02603056 0.005786990 -4.498117 2.359861e-04
# 5       South East - North East and Yorkshire               -0.02839780 0.006147364 -4.619508 1.334803e-04
# 6                       South East - Scotland               -0.04525784 0.012829825 -3.527549 1.250741e-02
# 7                     South East - North West               -0.04689172 0.007972765 -5.881488 1.458678e-07
# 12          London - North East and Yorkshire               -0.02180369 0.006670884 -3.268486 2.988393e-02
# 14                        London - North West               -0.04029761 0.008378280 -4.809772 5.296858e-05
# 17                 East of England - Midlands               -0.03186152 0.005810572 -5.483370 1.490709e-06
# 18 East of England - North East and Yorkshire               -0.03422876 0.006170381 -5.547268 1.037504e-06
# 19                 East of England - Scotland               -0.05108880 0.012843050 -3.977933 2.264178e-03
# 20               East of England - North West               -0.05272268 0.007991042 -6.597722 1.502050e-09
# 21                    East of England - Wales               -0.02171697 0.006597944 -3.291475 2.775706e-02
# 25                    South West - North West               -0.04287565 0.010686293 -4.012210 1.968868e-03
# 36                         North West - Wales                0.03100571 0.009555781  3.244707 3.223206e-02

table2csv(tukey_VOC, file=".\\multinomial_logistic_fits\\tables\\model 2a_VOCvsall_bGLMM_VOC_Tukey contrasts diff growth rates across regions.csv")





# 3.1.2 BINOMIAL MIXED MODEL FIT FOR VOC 202012/01 vs B.1.177 ####

# mixed binomial GLMM with random intercept for lad (LTLA) & 
# observation-level random effect to take into account overdispersion
# and fixed effects ns_name and date, with or without interaction effect nhs_name x date

start_date = as.Date("2020-08-01")
stop_date = max(data_agbydayregionlad$sample_date) # 6 Jan 2021
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
glmersettings = glmerControl(optimizer="optimx", optCtrl=list(method="nlminb")) # PS : to try all optimizer run all_fit(fit)
glmersettings2 = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=1E4)) # PS : to try all optimizer run all_fit(fit)
bGLMMfit1_vsB1177 = glmer(cbind(count, total-count) ~  (1|lad/obs) + 
                    nhs_name+scale(sample_date_num), 
                  family=binomial(logit), data=data_subs, control=glmersettings2)
bGLMMfit2_vsB1177 = glmer(cbind(count, total-count) ~  (1|lad/obs) + 
                    nhs_name*scale(sample_date_num), 
                  family=binomial(logit), data=data_subs, control=glmersettings2)


# saveRDS(bGLMMfit1_vsB1177, file = ".\\multinomial_logistic_fits\\fits\\bGLMMfit_VOCvsB1177_fit1_homog slopes.rds")
# saveRDS(bGLMMfit2_vsB1177, file = ".\\multinomial_logistic_fits\\fits\\bGLMMfit_VOCvsB1177_fit2_heter slopes_model 2b.rds")
# or to directly load previously fitted models
bGLMMfit1_vsB1177 = readRDS(file = ".\\multinomial_logistic_fits\\fits\\bGLMMfit_VOCvsB1177_fit1_homog slopes.rds")
bGLMMfit2_vsB1177 = readRDS(file = ".\\multinomial_logistic_fits\\fits\\bGLMMfit_VOCvsB1177_fit2_heter slopes_model 2b.rds")


# check BIC values
BIC(bGLMMfit1_vsB1177, bGLMMfit2_vsB1177)
#                   df      BIC
# bGLMMfit1_vsB1177 12 7301.352
# bGLMMfit2_vsB1177 20 7271.361


# bGLMMfit2_vsB1177 has the best BIC (heterog slopes across regions)


#  GROWTH RATE & TRANSMISSION ADVANTAGE OF VOC RELATIVE TO B.1.177

# on average across all regions, using the most parsimonious model bGLMMfit2_vsB1177, we get
bGLMM_VOC_growthrates_avg = as.data.frame(emtrends(bGLMMfit2_vsB1177, ~ 1, var="sample_date_num"))[,-c(3,4)] 
colnames(bGLMM_VOC_growthrates_avg)[2] = "logistic_growth_rate"
bGLMM_VOC_growthrates_avg = M.from.delta_r_df(bGLMM_VOC_growthrates_avg)
bGLMM_VOC_growthrates_avg
# 1         logistic_growth_rate asymp.LCL asymp.UCL       M1   M1.LCL   M1.UCL      M2   M2.LCL   M2.UCL
# 1 overall            0.1147566 0.1102424 0.1192707 1.879792 1.833696 1.927047 1.511532 1.487167 1.536296
table2csv(bGLMM_VOC_growthrates_avg, file=".\\multinomial_logistic_fits\\tables\\model 2b_VOCvsB1177_bGLMM_VOC_growthrates_UKavg_heterog slope.csv")



# 3.1.2 BINOMIAL MIXED MODEL FIT FOR VOC 202012/01 vs minor variants ####

# mixed binomial GLMM with random intercept for lad (LTLA) & 
# observation-level random effect to take into account overdispersion
# and fixed effects ns_name and date, with or without interaction effect nhs_name x date

start_date = as.Date("2020-08-01")
stop_date = max(data_agbydayregionlad$sample_date) # 6 Jan 2020
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
glmersettings = glmerControl(optimizer="optimx", optCtrl=list(method="nlminb")) # PS : to try all optimizer run all_fit(fit)
glmersettings2 = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=1E4)) # PS : to try all optimizer run all_fit(fit)
bGLMMfit1_vsminor = glmer(cbind(count, total-count) ~  (1|lad/obs) + 
                            nhs_name+scale(sample_date_num), 
                          family=binomial(logit), data=data_subs, control=glmersettings2)
bGLMMfit2_vsminor = glmer(cbind(count, total-count) ~  (1|lad/obs) + 
                            nhs_name*scale(sample_date_num), 
                          family=binomial(logit), data=data_subs, control=glmersettings2)

# saveRDS(bGLMMfit1_vsminor, file = ".\\multinomial_logistic_fits\\fits\\bGLMMfit_VOCvsminor_fit1_homog slopes.rds")
# saveRDS(bGLMMfit2_vsminor, file = ".\\multinomial_logistic_fits\\fits\\bGLMMfit_VOCvsminor_fit2_heter slopes_model 2c.rds")
# or to directly load previously fitted models
bGLMMfit1_vsminor = readRDS(file = ".\\multinomial_logistic_fits\\fits\\bGLMMfit_VOCvsminor_fit1_homog slopes.rds")
bGLMMfit2_vsminor = readRDS(file = ".\\multinomial_logistic_fits\\fits\\bGLMMfit_VOCvsminor_fit2_heter slopes_model 2c.rds")


# check BIC values
BIC(bGLMMfit1_vsminor, bGLMMfit2_vsminor)
#                   df      BIC
# bGLMMfit1_vsminor 12 5141.523
# bGLMMfit2_vsminor 20 5125.771




# bGLMMfit2_vsminor has the best BIC (heterog slopes across regions)


#  GROWTH RATE & TRANSMISSION ADVANTAGE OF VOC RELATIVE TO minority variants

# on average across all regions, using the most parsimonious model bGLMMfit2_vsminor, we get
bGLMM_VOC_growthrates_avg = as.data.frame(emtrends(bGLMMfit2_vsminor, ~ 1, var="sample_date_num"))[,-c(3,4)] 
colnames(bGLMM_VOC_growthrates_avg)[2] = "logistic_growth_rate"
bGLMM_VOC_growthrates_avg = M.from.delta_r_df(bGLMM_VOC_growthrates_avg)
bGLMM_VOC_growthrates_avg
# 1         logistic_growth_rate asymp.LCL asymp.UCL       M1   M1.LCL  M1.UCL       M2   M2.LCL   M2.UCL
# 1 overall            0.1330159 0.1271237 0.1389081 2.078377 2.012102 2.146835 1.614229 1.580348 1.648835
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
glmersettings = glmerControl(optimizer="optimx", optCtrl=list(method="nlminb")) # PS : to try all optimizer run all_fit(fit)
glmersettings2 = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=1E4)) # PS : to try all optimizer run all_fit(fit)

# mixed binomial GLMM with nested random intercept for lad and observation (the latter to take into account overdispersion)
# and fixed effects nhs_name and date, with or without interaction
bGLMMfit1_B1177 = glmer(cbind(count, total-count) ~  (1|lad/obs) + 
                            nhs_name+scale(sample_date_num), 
                          family=binomial(logit), data=data_subs, control=glmersettings2)
bGLMMfit2_B1177 = glmer(cbind(count, total-count) ~  (1|lad/obs) + 
                            nhs_name*scale(sample_date_num), 
                          family=binomial(logit), data=data_subs, control=glmersettings2)

# saveRDS(bGLMMfit1_B1177, file = ".\\multinomial_logistic_fits\\fits\\bGLMMfit_B1177vsall_fit1_homog slopes_model 2d.rds")
# saveRDS(bGLMMfit2_B1177, file = ".\\multinomial_logistic_fits\\fits\\bGLMMfit_B1177vsall_fit2_heter slopes_model 2e.rds")
# or to directly load previously fitted models
bGLMMfit1_B1177 = readRDS(file = ".\\multinomial_logistic_fits\\fits\\bGLMMfit_B1177vsall_fit1_homog slopes_model 2d.rds")
bGLMMfit2_B1177 = readRDS(file = ".\\multinomial_logistic_fits\\fits\\bGLMMfit_B1177vsall_fit2_heter slopes_model 2e.rds")

# check BIC values
BIC(bGLMMfit1_B1177, bGLMMfit2_B1177)
#                 df      BIC
# bGLMMfit1_B1177 12 7657.740
# bGLMMfit2_B1177 20 7661.333


# bGLMMfit1_vsB1177 has the best BIC (homogenous slopes across regions)


#  GROWTH RATE & TRANSMISSION ADVANTAGE

# on average across all regions, using the most parsimonious model bGLMMfit1_B1177, we get
bGLMM_B1177_growthrates_avg = as.data.frame(emtrends(bGLMMfit1_B1177, ~ 1, var="sample_date_num"))[,-c(3,4)] 
colnames(bGLMM_B1177_growthrates_avg)[2] = "logistic_growth_rate"
bGLMM_B1177_growthrates_avg = M.from.delta_r_df(bGLMM_B1177_growthrates_avg)
bGLMM_B1177_growthrates_avg
# 1         logistic_growth_rate  asymp.LCL  asymp.UCL      M1   M1.LCL   M1.UCL       M2   M2.LCL   M2.UCL
# 1 overall           0.06115157 0.05760982 0.06469332 1.399806 1.372802 1.427341 1.246258 1.230469 1.26225
table2csv(bGLMM_B1177_growthrates_avg, file=".\\multinomial_logistic_fits\\tables\\model 2d_B1177vsall_bGLMM_growthrates_UKavg_homog slope.csv")

# growth rates per region for heterogeneous slope model bGLMMfit2_B1177
bGLMM_B1177_growthrates_region = as.data.frame(emtrends(bGLMMfit2_B1177, ~ nhs_name, var="sample_date_num"))[,-c(3,4)] 
colnames(bGLMM_B1177_growthrates_region)[2] = "logistic_growth_rate"
bGLMM_B1177_growthrates_region = M.from.delta_r_df(bGLMM_B1177_growthrates_region)
bGLMM_B1177_growthrates_region
#                 nhs_name logistic_growth_rate  asymp.LCL  asymp.UCL       M1   M1.LCL   M1.UCL       M2   M2.LCL   M2.UCL
# 1               South East           0.05748963 0.04675292 0.06822633 1.371895 1.293228 1.455348 1.229937 1.183304 1.278407
# 2                   London           0.03491949 0.02548835 0.04435063 1.211740 1.150488 1.276253 1.133953 1.096100 1.173115
# 3          East of England           0.04383788 0.03341639 0.05425936 1.272659 1.201764 1.347736 1.170951 1.127834 1.215717
# 4               South West           0.04849638 0.03491154 0.06208122 1.305688 1.211687 1.406982 1.190754 1.133921 1.250436
# 5                 Midlands           0.06913330 0.05883372 0.07943289 1.462626 1.382074 1.547872 1.282588 1.235902 1.331037
# 6 North East and Yorkshire           0.05884118 0.05020587 0.06747650 1.382131 1.318022 1.449358 1.235936 1.198105 1.274961
# 7                 Scotland           0.06150159 0.05397367 0.06902952 1.402503 1.345620 1.461791 1.247830 1.214467 1.282109
# 8               North West           0.08090281 0.07083684 0.09096877 1.560436 1.476394 1.649263 1.338099 1.290478 1.387478
# 9                    Wales           0.07088988 0.06077399 0.08100577 1.476825 1.396902 1.561320 1.290724 1.244565 1.338595
table2csv(bGLMM_B1177_growthrates_region, 
          file=".\\multinomial_logistic_fits\\tables\\model 2e_B1177vsall_bGLMM_growthrates_UK_by region.csv")


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
    xlim=c(as.Date("2020-06-01"),as.Date("2021-02-28")), 
    ylim=c(0.0001,0.999001), expand=c(0,0)) +
  scale_color_discrete("", h=c(0, 280), c=200) +
  scale_fill_discrete("", h=c(0, 280), c=200) +
  geom_point(data=data_agbyweekregion[data_agbyweekregion$variant=="B.1.177"&
                                        data_agbyweekregion$sample_date>="2020-06-01"&
                                        data_agbyweekregion$sample_date<="2020-09-30",], 
             aes(x=sample_date, y=prop, colour=nhs_name, size=total), alpha=I(0.5)) +
  scale_size_continuous("total number\nof sequences\nper week", trans="sqrt", 
                        range=c(0.01, 4), limits=c(1,max(data_agbyweekregion$total)), breaks=c(10,100,1000)) +
  guides(fill=FALSE) + guides(colour=FALSE) + theme(legend.position = "right")
  # geom_point(data=data_agbydayregion[data_agbydayregion$variant=="VOC 202012/01",], 
  #         aes(x=sample_date, y=prop, colour=nhs_name, size=total), alpha=I(0.4)) +
  # theme(legend.position = "none")
plot_bGLMM_B1177_hom
saveRDS(plot_bGLMM_B1177_hom, file = ".\\multinomial_logistic_fits\\plots\\model 2d_plot B1177vsall_fit bGLMM_homog slopes.rds")
graph2ppt(file=".\\multinomial_logistic_fits\\plots\\model 2d_plot B1177vsall_fit bGLMM_homog slopess.pptx", width=8, height=6)
ggsave(file=".\\multinomial_logistic_fits\\plots\\model 2d_plot B1177vsall_fit bGLMM_homog slopes.png", width=8, height=6)
ggsave(file=".\\multinomial_logistic_fits\\plots\\model 2d_plot B1177vsall_fit bGLMM_homog slopes.pdf", width=8, height=6)


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
    xlim=c(as.Date("2020-06-01"),as.Date("2021-02-28")), 
    ylim=c(0.0001,0.999001), expand=c(0,0)) +
  scale_color_discrete("", h=c(0, 280), c=200) +
  scale_fill_discrete("", h=c(0, 280), c=200) +
  geom_point(data=data_agbyweekregion[data_agbyweekregion$variant=="B.1.177"&
                                        data_agbyweekregion$sample_date>="2020-06-01"&
                                        data_agbyweekregion$sample_date<="2020-09-30",], 
             aes(x=sample_date, y=prop, colour=nhs_name, size=total), alpha=I(0.5)) +
  scale_size_continuous("total number\nof sequences\nper week", trans="sqrt", 
                        range=c(0.01, 4), limits=c(1,max(max(data_agbyweekregion$total))), breaks=c(10,100,1000)) +
  guides(fill=FALSE) + guides(colour=FALSE) + theme(legend.position = "right")
# geom_point(data=data_agbydayregion[data_agbydayregion$variant=="VOC 202012/01",], 
#         aes(x=sample_date, y=prop, colour=nhs_name, size=total), alpha=I(0.4)) +
# theme(legend.position = "none")
plot_bGLMM_B1177_het
saveRDS(plot_bGLMM_B1177_het, file = ".\\multinomial_logistic_fits\\plots\\model 2e_plot B1177vsall_fit bGLMM_heter slopes.rds")
graph2ppt(file=".\\multinomial_logistic_fits\\plots\\model 2e_plot B1177vsall_fit bGLMM_heter slopess.pptx", width=8, height=6)
ggsave(file=".\\multinomial_logistic_fits\\plots\\model 2e_plot B1177vsall_fit bGLMM_heter slopes.png", width=8, height=6)
ggsave(file=".\\multinomial_logistic_fits\\plots\\model 2e_plot B1177vsall_fit bGLMM_heter slopes.pdf", width=8, height=6)




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
glmersettings = glmerControl(optimizer="optimx", optCtrl=list(method="nlminb")) # PS : to try all optimizer run all_fit(fit)
glmersettings2 = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=1E4)) # PS : to try all optimizer run all_fit(fit)

# mixed binomial GLMM with nested random intercept for lad and observation (the latter to take into account overdispersion)
# and fixed effects nhs_name and date, with or without interaction
bGLMMfit1_B1177_vsminority = glmer(cbind(count, total-count) ~  (1|lad/obs) + 
                          nhs_name+scale(sample_date_num), 
                        family=binomial(logit), data=data_subs, control=glmersettings2)
bGLMMfit2_B1177_vsminority = glmer(cbind(count, total-count) ~  (1|lad/obs) + 
                          nhs_name*scale(sample_date_num), 
                        family=binomial(logit), data=data_subs, control=glmersettings2)

# saveRDS(bGLMMfit1_B1177_vsminority, file = ".\\multinomial_logistic_fits\\fits\\bGLMMfit_B1177vsminor_fit1_homog slopes_model 2f_best model.rds")
# saveRDS(bGLMMfit2_B1177_vsminority, file = ".\\multinomial_logistic_fits\\fits\\bGLMMfit_B1177vsminor_fit2_heter slopes_model 2g.rds")
# or to directly load previously fitted models
bGLMMfit1_B1177_vsminority = readRDS(file = ".\\multinomial_logistic_fits\\fits\\bGLMMfit_B1177vsminor_fit1_homog slopes_model 2f_best model.rds")
bGLMMfit2_B1177_vsminority = readRDS(file = ".\\multinomial_logistic_fits\\fits\\bGLMMfit_B1177vsminor_fit2_heter slopes_model 2g.rds")

# check BIC values
BIC(bGLMMfit1_B1177_vsminority, bGLMMfit2_B1177_vsminority)
#                            df      BIC
# bGLMMfit1_B1177_vsminority 12 7022.226
# bGLMMfit2_B1177_vsminority 20 7026.006


# bGLMMfit1_B1177_vsminority has the best BIC (homog slopes across regions)


#  GROWTH RATE & TRANSMISSION ADVANTAGE

# on average across all regions, using the most parsimonious homog slope model bGLMMfit1_B1177_vsminority, we get
bGLMM_B1177_growthrates_avg_vsminority = as.data.frame(emtrends(bGLMMfit1_B1177_vsminority, ~ 1, var="sample_date_num"))[,-c(3,4)] 
colnames(bGLMM_B1177_growthrates_avg_vsminority)[2] = "logistic_growth_rate"
bGLMM_B1177_growthrates_avg_vsminority = M.from.delta_r_df(bGLMM_B1177_growthrates_avg_vsminority)
bGLMM_B1177_growthrates_avg_vsminority
# 1            logistic_growth_rate  asymp.LCL  asymp.UCL      M1   M1.LCL   M1.UCL       M2   M2.LCL  M2.UCL
# 1 overall            0.05722866 0.05365433 0.06080299 1.369927 1.343259 1.397125 1.228782 1.213072 1.244695
table2csv(bGLMM_B1177_growthrates_avg_vsminority, file=".\\multinomial_logistic_fits\\tables\\model 2f_B1177vsminor_bGLMM_growthrates_UKavg_homog slopes.csv")

# growth rates per region for heterogeneous slope model bGLMMfit2_B1177
bGLMM_B1177_growthrates_region = as.data.frame(emtrends(bGLMMfit2_B1177_vsminority, ~ nhs_name, var="sample_date_num"))[,-c(3,4)] 
colnames(bGLMM_B1177_growthrates_region)[2] = "logistic_growth_rate"
bGLMM_B1177_growthrates_region = M.from.delta_r_df(bGLMM_B1177_growthrates_region)
bGLMM_B1177_growthrates_region
#                   nhs_name logistic_growth_rate  asymp.LCL  asymp.UCL       M1   M1.LCL   M1.UCL       M2   M2.LCL   M2.UCL
# 1               South East           0.05580314 0.04521897 0.06638731 1.359228 1.282363 1.440702 1.222492 1.176788 1.269971
# 2                   London           0.03128840 0.02241861 0.04015820 1.187780 1.131226 1.247161 1.119227 1.084053 1.155542
# 3          East of England           0.03766634 0.02753763 0.04779504 1.230185 1.163528 1.300661 1.145222 1.104216 1.187752
# 4               South West           0.04829084 0.03476019 0.06182149 1.304213 1.210679 1.404973 1.189874 1.133303 1.249267
# 5                 Midlands           0.06624827 0.05611741 0.07637914 1.439600 1.361580 1.522091 1.269336 1.223876 1.316484
# 6 North East and Yorkshire           0.05693563 0.04841146 0.06545980 1.367721 1.305078 1.433371 1.227486 1.190390 1.265738
# 7                 Scotland           0.05816930 0.04963971 0.06669889 1.377033 1.313924 1.443173 1.232950 1.195666 1.271397
# 8               North West           0.07574652 0.06592206 0.08557098 1.516805 1.437020 1.601019 1.313490 1.267846 1.360777
# 9                    Wales           0.06234078 0.05264878 0.07203278 1.408992 1.335851 1.486137 1.251605 1.208688 1.296046
table2csv(bGLMM_B1177_growthrates_region, 
          file=".\\multinomial_logistic_fits\\tables\\model 2g_B1177vsminor_bGLMM_growthrates_UK_by region_heter slopes.csv")




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
    xlim=c(as.Date("2020-06-01"),as.Date("2021-02-28")), 
    ylim=c(0.0001,0.999001), expand=c(0,0)) +
  scale_color_discrete("", h=c(0, 280), c=200) +
  scale_fill_discrete("", h=c(0, 280), c=200) +
  geom_point(data=data_agbyweekregion[data_agbyweekregion$variant=="B.1.177"&
                                        data_agbyweekregion$sample_date>="2020-06-01"&
                                        data_agbyweekregion$sample_date<="2020-09-30",], 
             aes(x=sample_date, y=prop, colour=nhs_name, size=total), alpha=I(0.5)) +
  scale_size_continuous("total number\nof sequences\nper week", trans="sqrt", 
                        range=c(0.01, 4), limits=c(1,max(data_agbyweekregion$total)), breaks=c(10,100,1000)) +
  guides(fill=FALSE) + guides(colour=FALSE) + theme(legend.position = "right")
# geom_point(data=data_agbydayregion[data_agbydayregion$variant=="VOC 202012/01",], 
#         aes(x=sample_date, y=prop, colour=nhs_name, size=total), alpha=I(0.4)) +
# theme(legend.position = "none")
plot_bGLMM_B1177_preds_vsminority_het

saveRDS(plot_bGLMM_B1177_preds_vsminority_het, file = ".\\multinomial_logistic_fits\\plots\\model2g_plot B1177vsminor_fit bGLMM_heter slopes.rds")
graph2ppt(file=".\\multinomial_logistic_fits\\plots\\model2g_plot B1177vsminor_fit bGLMM_heter slopes.pptx", width=8, height=6)
ggsave(file=".\\multinomial_logistic_fits\\plots\\model2g_plot B1177vsminor_fit bGLMM_heter slopes.png", width=8, height=6)
ggsave(file=".\\multinomial_logistic_fits\\plots\\model2g_plot B1177vsminor_fit bGLMM_heter slopes.pdf", width=8, height=6)


# multipanel for suppl Fig. S5
plot_bGLMMVOC_B1177_multipanel = ggarrange(plot_bGLMMVOC_het, 
          plot_bGLMM_B1177_preds_vsminority_het,
          ncol=1, common.legend=TRUE, legend="right")
plot_bGLMMVOC_B1177_multipanel
saveRDS(plot_bGLMMVOC_B1177_multipanel, file = ".\\multinomial_logistic_fits\\plots\\FigS5_models_2a_and_2g_VOCvsall_B1177vsminor_fit bGLMM_heter slopes.rds")
graph2ppt(file=".\\multinomial_logistic_fits\\plots\\FigS5_models_2a_and_2g_VOCvsall_B1177vsminor_fit bGLMM_heter slopes.pptx", width=6, height=8)
ggsave(file=".\\multinomial_logistic_fits\\plots\\FigS5_models_2a_and_2g_VOCvsall_B1177vsminor_fit bGLMM_heter slopes.png", width=6, height=8)
ggsave(file=".\\multinomial_logistic_fits\\plots\\FigS5_models_2a_and_2g_VOCvsall_B1177vsminor_fit bGLMM_heter slopes.pdf", width=6, height=8)



# TUKEY POSTHOC TESTS TO TEST FOR DIFFERENCES IN GROWTH RATES IN THE B1.177 VARIANT ACROSS DIFFERENT NHS REGIONS

tukey_B1177_vsminority = as.data.frame(emtrends(bGLMMfit2_B1177_vsminority, pairwise ~ nhs_name, var="sample_date_num", adjust="tukey")$contrasts)[,-4]
colnames(tukey_B1177_vsminority)[2] = "diff_logistic_growth_rate"
tukey_B1177_vsminority
#                                 contrast diff_logistic_growth_rate          SE    z.ratio      p.value
# 1                         South East - London               0.024514735 0.007000001  3.5021046 1.367551e-02
# 2                South East - East of England               0.018136800 0.007438159  2.4383452 2.630139e-01
# 3                     South East - South West               0.007512300 0.008727492  0.8607627 9.948515e-01
# 4                       South East - Midlands              -0.010445137 0.007450908 -1.4018610 8.975464e-01
# 5       South East - North East and Yorkshire              -0.001132496 0.006884202 -0.1645065 1.000000e+00
# 6                       South East - Scotland              -0.002366164 0.006885071 -0.3436658 9.999946e-01
# 7                     South East - North West              -0.019943380 0.007309906 -2.7282678 1.378319e-01
# 8                          South East - Wales              -0.006537645 0.007277191 -0.8983748 9.931130e-01
# 9                    London - East of England              -0.006377935 0.006814495 -0.9359366 9.909368e-01
# 10                        London - South West              -0.017002435 0.008199928 -2.0734860 4.916400e-01
# 11                          London - Midlands              -0.034959872 0.006832972 -5.1163494 1.104234e-05
# 12          London - North East and Yorkshire              -0.025647231 0.006199995 -4.1366534 1.172252e-03
# 13                          London - Scotland              -0.026880898 0.006200578 -4.3352246 4.945730e-04
# 14                        London - North West              -0.044458115 0.006665025 -6.6703597 9.181322e-10
# 15                             London - Wales              -0.031052380 0.006635035 -4.6800627 9.987285e-05
# 16               East of England - South West              -0.010624500 0.008579053 -1.2384234 9.480593e-01
# 17                 East of England - Midlands              -0.028581937 0.007280170 -3.9259988 2.791031e-03
# 18 East of England - North East and Yorkshire              -0.019269296 0.006695039 -2.8781456 9.395714e-02
# 19                 East of England - Scotland              -0.020502964 0.006695790 -3.0620678 5.620250e-02
# 20               East of England - North West              -0.038080180 0.007129803 -5.3409862 3.294271e-06
# 21                    East of England - Wales              -0.024674445 0.007098478 -3.4760192 1.497388e-02
# 22                      South West - Midlands              -0.017957436 0.008594664 -2.0893703 4.805796e-01
# 23      South West - North East and Yorkshire              -0.008644795 0.008100330 -1.0672153 9.788018e-01
# 24                      South West - Scotland              -0.009878463 0.008100774 -1.2194468 9.524918e-01
# 25                    South West - North West              -0.027455679 0.008460782 -3.2450523 3.219683e-02
# 26                         South West - Wales              -0.014049945 0.008437110 -1.6652557 7.678832e-01
# 27        Midlands - North East and Yorkshire               0.009312641 0.006714817  1.3868794 9.031382e-01
# 28                        Midlands - Scotland               0.008078973 0.006715888  1.2029642 9.561176e-01
# 29                      Midlands - North West              -0.009498243 0.007153527 -1.3277707 9.232740e-01
# 30                           Midlands - Wales               0.003907491 0.007117345  0.5490097 9.998007e-01
# 31        North East and Yorkshire - Scotland              -0.001233668 0.006067510 -0.2033235 9.999999e-01
# 32      North East and Yorkshire - North West              -0.018810884 0.006540521 -2.8760528 9.448230e-02
# 33           North East and Yorkshire - Wales              -0.005405149 0.006511457 -0.8300984 9.959899e-01
# 34                      Scotland - North West              -0.017577216 0.006540666 -2.6873743 1.521518e-01
# 35                           Scotland - Wales              -0.004171482 0.006511999 -0.6405840 9.993722e-01
# 36                         North West - Wales               0.013405735 0.006954720  1.9275737 5.944208e-01



tukey_B1177_vsminority[tukey_B1177_vsminority$p.value<0.05,]
#                        contrast diff_logistic_growth_rate          SE   z.ratio      p.value
# 1                South East - London                0.02451473 0.007000001  3.502105 1.367551e-02
# 11                 London - Midlands               -0.03495987 0.006832972 -5.116349 1.104234e-05
# 12 London - North East and Yorkshire               -0.02564723 0.006199995 -4.136653 1.172252e-03
# 13                 London - Scotland               -0.02688090 0.006200578 -4.335225 4.945730e-04
# 14               London - North West               -0.04445811 0.006665025 -6.670360 9.181322e-10
# 15                    London - Wales               -0.03105238 0.006635035 -4.680063 9.987285e-05
# 17        East of England - Midlands               -0.02858194 0.007280170 -3.925999 2.791031e-03
# 20      East of England - North West               -0.03808018 0.007129803 -5.340986 3.294271e-06
# 21           East of England - Wales               -0.02467445 0.007098478 -3.476019 1.497388e-02
# 25           South West - North West               -0.02745568 0.008460782 -3.245052 3.219683e-02

table2csv(tukey_B1177_vsminority, file=".\\multinomial_logistic_fits\\tables\\model 2g_B1177vsminority_bGLMM_Tukey contrasts diff growth rates across regions.csv")




# 4. ANALYSIS OF PILLAR 2 S-GENE TARGET FAILURE DATA UK ####

levels_UKregions = c("South East","London","East of England",
                     "South West","Midlands","North East and Yorkshire",
                     "Scotland","North West","Wales")


sgtfdata_uk = read.csv(".//fitting_data//sgtf-2021-01-18.csv") # Pillar 2 S gene targeted failure data (SGTF) (S dropout)
sgtfdata_uk$other = sgtfdata_uk$other+sgtfdata_uk$sgtf
colnames(sgtfdata_uk) = c("collection_date","REGION","SGTF","TOTAL")
sgtfdata_uk_truepos = read.csv(".//data//sgtfvoc.csv") # modelled proportion of S dropout that was actually the VOC
sgtfdata_uk$TRUEPOS = sgtfdata_uk_truepos$sgtfv[match(interaction(sgtfdata_uk$REGION, sgtfdata_uk$collection_date),
                                                      interaction(sgtfdata_uk_truepos$nhs_name, sgtfdata_uk_truepos$date))] # modelled proportion of S dropout samples that were actually the VOC
sgtfdata_uk$est_n_B117 = sgtfdata_uk$SGTF * sgtfdata_uk$TRUEPOS
sgtfdata_uk$COUNTRY = "UK"
sgtfdata_uk = sgtfdata_uk[,c("collection_date","COUNTRY","REGION","est_n_B117","TOTAL")]
colnames(sgtfdata_uk)[which(colnames(sgtfdata_uk)=="TOTAL")] = "n_pos"
range(sgtfdata_uk$collection_date) # "2020-10-01" "2021-01-17"
sgtfdata_uk$collection_date = as.Date(sgtfdata_uk$collection_date)
sgtfdata_uk$collection_date_num = as.numeric(sgtfdata_uk$collection_date)
sgtfdata_uk$REGION = factor(sgtfdata_uk$REGION, levels=levels_UKregions)
sgtfdata_uk$REGION = droplevels(sgtfdata_uk$REGION)
sgtfdata_uk$obs = factor(1:nrow(sgtfdata_uk))
sgtfdata_uk$propB117 = sgtfdata_uk$est_n_B117 / sgtfdata_uk$n_pos
head(sgtfdata_uk)

set_sum_contrasts()
glmersettings = glmerControl(optimizer="Nelder_Mead", optCtrl=list(maxfun=1e5)) # bobyqa, PS : to try all optimizer run all_fit(fit1)
glmersettings2 = glmerControl(optimizer="optimx", optCtrl=list(method="L-BFGS-B"))
glmersettings3 = glmerControl(optimizer="optimx", optCtrl=list(method="nlminb"))
glmersettings4 = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=1e5))
fit_ukSGTF_1 = glmer(cbind(est_n_B117, n_pos-est_n_B117 ) ~ (1|obs)+scale(collection_date_num)+REGION, family=binomial(logit), 
                     data=sgtfdata_uk, control=glmersettings)  # common slope model
fit_ukSGTF_2 = glmer(cbind(est_n_B117, n_pos-est_n_B117) ~ (1|obs)+scale(collection_date_num)*REGION, family=binomial(logit), 
                     data=sgtfdata_uk, control=glmersettings3) # heter slope model
fit_ukSGTF_3 = glmer(cbind(est_n_B117, n_pos-est_n_B117) ~ (1|obs)+scale(ns(collection_date_num,df=3))+REGION, family=binomial(logit), 
                     data=sgtfdata_uk, control=glmersettings3) # with additive spline term
fit_ukSGTF_4 = glmer(cbind(est_n_B117, n_pos-est_n_B117) ~ (1|obs)+ns(collection_date_num,df=3)*REGION, family=binomial(logit), 
                     data=sgtfdata_uk, control=glmersettings3) # with spline term in interaction with region
BIC(fit_ukSGTF_1, fit_ukSGTF_2, fit_ukSGTF_3, fit_ukSGTF_4) 
# separate-slopes 3 df spline model fit_be_uk2_4 best
# df      BIC
# fit_ukSGTF_1  9 4902.696
# fit_ukSGTF_2 15 4769.405
# fit_ukSGTF_3 11 4905.592
# fit_ukSGTF_4 29 4428.474


# model fit_ukSGTF_4 best

summary(fit_ukSGTF_4)

# avg growth rate advantage for UK (difference in growth rate between B.1.1.7 and old strains):
fit_ukSGTF_4_emtrends = as.data.frame(emtrends(fit_ukSGTF_4, revpairwise ~ 1, 
                                               var="collection_date_num",
                                               mode="link", adjust="Tukey")$emtrends)
fit_ukSGTF_4_emtrends[,c(2,5,6)]
#   collection_date_num.trend asymp.LCL asymp.UCL
# 1                 0.1093804 0.1074621 0.1112988

# with a generation time of 4.7 days this would translate in an increased 
# infectiousness (multiplicative effect on Rt) of
exp(fit_ukSGTF_4_emtrends[,c(2,5,6)]*4.7) 
#    collection_date_num.trend asymp.LCL asymp.UCL
# 1                  1.672113  1.657104  1.687257


# PLOT MODEL FIT

# spline model fit_ukSGTF_4
date.from = as.numeric(as.Date("2020-09-01"))
date.to = as.numeric(as.Date("2021-03-01")) # date to extrapolate to
total.SD = sqrt(sum(sapply(as.data.frame(VarCorr(fit_ukSGTF_4))$sdcor, function (x) x^2))) 
# bias correction for random effects in marginal means, see https://cran.r-project.org/web/packages/emmeans/vignettes/transformations.html#bias-adj
fit_ukSGTF_4_preds = as.data.frame(emmeans(fit_ukSGTF_4, ~ collection_date_num, 
                                           by=c("REGION"), 
                                           at=list(collection_date_num=seq(date.from,
                                                                           date.to)), 
                                           type="response"), bias.adjust = TRUE, sigma = total.SD)
fit_ukSGTF_4_preds$collection_date = as.Date(fit_ukSGTF_4_preds$collection_date_num, origin="1970-01-01")

n = length(levels(fit_ukSGTF_4_preds$REGION))
reg_cols = hcl(h = seq(290, 0, length = n + 1), l = 50, c = 255)[1:n]
# reg_cols[2:n] = rev(reg_cols[2:n])

fit_ukSGTF_4_preds$REGION = factor(fit_ukSGTF_4_preds$REGION, levels=unique(fit_ukSGTF_4_preds$REGION))
sgtfdata_uk$REGION = factor(sgtfdata_uk$REGION, levels=levels(fit_ukSGTF_4_preds$REGION))

# PLOT MODEL FIT (logit scale):
plot_UK_SGTF = qplot(data=fit_ukSGTF_4_preds, x=collection_date, y=prob, geom="blank") +
  # facet_wrap(~COUNTRY) +
  geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL, 
                  fill=REGION
  ), 
  # fill=I("steelblue"), 
  alpha=I(0.3)) +
  geom_line(aes(y=prob, 
                colour=REGION
  ), 
  # colour=I("steelblue"), 
  alpha=I(0.8)) +
  ylab("Relative abundance of B.1.1.7 (%)") +
  theme_hc() + 
  xlab("") + 
  # scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
  #                   labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M")) +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(# xlim=c(as.Date("2020-09-01"),as.Date("2021-02-01")), 
    # xlim=c(as.Date("2020-07-01"),as.Date("2021-01-31")), 
    ylim=c(0.001,99.9), expand=c(0,0)) +
  scale_color_manual("", values=reg_cols) +
  scale_fill_manual("", values=reg_cols) +
  # scale_color_discrete("", h=c(0, 280), c=200) +
  # scale_fill_discrete("", h=c(0, 280), c=200) +
  geom_point(data=sgtfdata_uk, 
             aes(x=collection_date, y=propB117, size=n_pos,
                 colour=REGION
             ), 
             # colour=I("steelblue"), 
             alpha=I(0.5)) +
  scale_size_continuous("number of\npositive tests (Ct<30)", trans="sqrt", 
                        range=c(1, 4), limits=c(1,10000), breaks=c(100,1000,10000)) +
  # guides(fill=FALSE) + 
  # guides(colour=FALSE) + 
  theme(legend.position = "right") +
  xlab("Collection date") +
  ggtitle("UK") +
  theme(plot.title = element_text(hjust = 0.5))
plot_UK_SGTF


# PLOT MODEL FIT (response scale):
plot_UK_SGTF_response = qplot(data=fit_ukSGTF_4_preds, x=collection_date, y=prob*100, geom="blank") +
  # facet_wrap(~COUNTRY) +
  geom_ribbon(aes(y=prob*100, ymin=asymp.LCL*100, ymax=asymp.UCL*100, colour=NULL, 
                  fill=REGION
  ), 
  # fill=I("steelblue"), 
  alpha=I(0.3)) +
  geom_line(aes(y=prob*100, 
                colour=REGION
  ), 
  # colour=I("steelblue"), 
  alpha=I(0.8)) +
  ylab("Relative abundance of B.1.1.7 (%)") +
  theme_hc() + 
  xlab("") + 
  # scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
  #                   labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M")) +
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                    labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(# xlim=c(as.Date("2020-09-01"),as.Date("2021-02-01")), 
    # xlim=c(as.Date("2020-07-01"),as.Date("2021-01-31")), 
    ylim=c(0,100), expand=c(0,0)) +
  scale_color_manual("", values=reg_cols) +
  scale_fill_manual("", values=reg_cols) +
  # scale_color_discrete("", h=c(0, 280), c=200) +
  # scale_fill_discrete("", h=c(0, 280), c=200) +
  geom_point(data=sgtfdata_uk, 
             aes(x=collection_date, y=propB117*100, size=n_pos,
                 colour=REGION
             ), 
             # colour=I("steelblue"), 
             alpha=I(0.5)) +
  scale_size_continuous("number of\npositive tests (Ct<30)", trans="sqrt", 
                        range=c(1, 4), limits=c(1,10000), breaks=c(100,1000,10000)) +
  # guides(fill=FALSE) + 
  # guides(colour=FALSE) + 
  theme(legend.position = "right") +
  xlab("Collection date") +
  ggtitle("UK") +
  theme(plot.title = element_text(hjust = 0.5))
plot_UK_SGTF_response



# 5. INTERNATIONAL COMPARISONS: COMPETITIVE ADVANTAGE OF VOC IN DENMARK, SWITZERLAND & THE USA ####

# 5.1. DATA DENMARK: SEQUENCING DATA ####

# Data source: Danish Covid-19 Genome Consortium & the Statens Serum Institut, https://www.covid19genomics.dk/statistics
# downloaded on the 16th of February 2021

data_denmark = read.csv(".//multinomial_logistic_fits//data/dk//data_denmark_20210216.csv", sep=";", dec=",")
data_denmark$percent = NULL
data_denmark$Region = gsub("SjÃ¦lland","Sjælland",data_denmark$Region)
data_denmark$WEEK = sapply(data_denmark$Week, function(s) as.numeric(strsplit(s, "W")[[1]][[2]]))
data_denmark$date = as.Date(NA)
data_denmark$date[data_denmark$WEEK>=42] = lubridate::ymd( "2020-01-01" ) + 
  lubridate::weeks( data_denmark$WEEK[data_denmark$WEEK>=42] - 1 ) + 1
data_denmark$date[data_denmark$WEEK<42] = lubridate::ymd( "2021-01-01" ) + 
  lubridate::weeks( data_denmark$WEEK[data_denmark$WEEK<42] - 1 ) + 6 
data_denmark$date_num = as.numeric(data_denmark$date)
data_denmark$obs = factor(1:nrow(data_denmark))
colnames(data_denmark)[colnames(data_denmark) %in% c("yes")] = "n_B117"
data_denmark$propB117 = data_denmark$n_B117 / data_denmark$total

data_denmark_whole = data_denmark[data_denmark$Region=="Whole Denmark",]
data_denmark = data_denmark[data_denmark$Region!="Whole Denmark",]
levels_DK = c("Syddanmark","Sjælland","Nordjylland","Hovedstaden","Midtjylland")
data_denmark$Region = factor(data_denmark$Region, levels=levels_DK)
range(data_denmark$date) # "2020-10-15" "2021-01-28"

fit_denmark1 = glmer(cbind(n_B117,total-n_B117) ~ (1|obs) + Region + scale(date_num), family=binomial(logit), data=data_denmark)
fit_denmark2 = glmer(cbind(n_B117,total-n_B117) ~ (1|obs) + Region * scale(date_num), family=binomial(logit), data=data_denmark)
fit_denmark3 = glmer(cbind(n_B117,total-n_B117) ~ (1|Region/obs) + scale(date_num), family=binomial(logit), data=data_denmark)
fit_denmark4 = glmer(cbind(n_B117,total-n_B117) ~ (date_num||Region/obs) + scale(date_num), family=binomial(logit), data=data_denmark)
BIC(fit_denmark1, fit_denmark2, fit_denmark3, fit_denmark4)
# df      BIC
# fit_denmark1  7 381.9215
# fit_denmark2 11 375.6446
# fit_denmark3  4 373.6176
# fit_denmark4  6 382.3817

summary(fit_denmark3)

# common-slope model fit_denmark3 with nested random intercepts fits best

#  GROWTH RATE & TRANSMISSION ADVANTAGE

# on average across all regions, using the most parsimonious model fit_denmark3, we get
dk_growthrates_avg_B117vsallother = as.data.frame(emtrends(fit_denmark3, ~ 1, var="date_num"))[,-c(3,4)] 
colnames(dk_growthrates_avg_B117vsallother)[2] = "logistic_growth_rate"
dk_growthrates_avg_B117vsallother = M.from.delta_r_df(dk_growthrates_avg_B117vsallother)
dk_growthrates_avg_B117vsallother
# 1         logistic_growth_rate  asymp.LCL  asymp.UCL       M1   M1.LCL   M1.UCL       M2   M2.LCL   M2.UCL
# 1 overall            0.0798134 0.06716476 0.09246204 1.551115 1.446875 1.662864 1.332862 1.273531 1.394957
table2csv(dk_growthrates_avg_B117vsallother, file=".\\multinomial_logistic_fits\\tables\\model 3a_B1177vsallother_bGLMM_DK.csv")


# PLOT MODEL FIT
date.from = as.numeric(as.Date("2020-09-01"))
date.to = as.numeric(as.Date("2021-03-01")) # date to extrapolate to
total.SD = sqrt(sum(sapply(as.data.frame(VarCorr(fit_denmark3))$sdcor, function (x) x^2))) 
# bias correction for random effects in marginal means, see https://cran.r-project.org/web/packages/emmeans/vignettes/transformations.html#bias-adj
fit_denmark_preds = as.data.frame(emmeans(fit_denmark3, ~ date_num, 
                                          at=list(date_num=seq(date.from,
                                                               date.to)), 
                                          type="response"), bias.adjust = TRUE, sigma = total.SD)
fit_denmark_preds$date = as.Date(fit_denmark_preds$date_num, origin="1970-01-01")

# n = length(levels(fit_denmark_preds$Region))
# reg_cols = hcl(h = seq(290, 0, length = n + 1), l = 50, c = 255)[1:n]

# PLOT MODEL FIT (response scale)
plot_denmark = qplot(data=fit_denmark_preds, x=date, y=prob, geom="blank") +
  geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL # , 
                  # fill=Region
  ), 
  fill=I("steelblue"), 
  alpha=I(0.3)) +
  geom_line(aes(y=prob# , 
                # colour=Region
  ), 
  colour=I("steelblue"), 
  alpha=I(0.8)) +
  ylab("Relative abundance of B.1.1.7 (%)") +
  theme_hc() + 
  xlab("") + 
  # scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
  #                   labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M")) +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(# xlim=c(as.Date("2020-09-01"),as.Date("2021-02-01")), 
    xlim=c(as.Date("2020-10-01"),as.Date("2021-03-01")), 
    ylim=c(0.001,99.9), expand=c(0,0)) +
  scale_color_manual("", values=reg_cols) +
  scale_fill_manual("", values=reg_cols) +
  # scale_color_discrete("", h=c(0, 280), c=200) +
  # scale_fill_discrete("", h=c(0, 280), c=200) +
  geom_point(data=data_denmark_whole, 
             aes(x=date, y=propB117, size=total,
                 # colour=Region
             ), 
             colour=I("steelblue"), 
             alpha=I(0.5)) +
  scale_size_continuous("number of\nsequences", trans="sqrt", 
                        range=c(1, 4), limits=c(1,max(data_denmark_whole$total)), breaks=c(10,100,1000)) +
  # guides(fill=FALSE) + 
  # guides(colour=FALSE) + 
  theme(legend.position = "right") +
  xlab("Collection date") +
  ggtitle("DENMARK") +
  theme(plot.title = element_text(hjust = 0.5))
plot_denmark


# PLOT MODEL FIT (response scale)
plot_denmark_response = qplot(data=fit_denmark_preds, x=date, y=prob*100, geom="blank") +
  geom_ribbon(aes(y=prob*100, ymin=asymp.LCL*100, ymax=asymp.UCL*100, colour=NULL # , 
                  # fill=Region
  ), 
  fill=I("steelblue"), 
  alpha=I(0.3)) +
  geom_line(aes(y=prob*100 # , 
                # colour=Region
  ), 
  colour=I("steelblue"), 
  alpha=I(0.8)) +
  ylab("Relative abundance of B.1.1.7 (%)") +
  theme_hc() + 
  xlab("") + 
  # scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
  #                   labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M")) +
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                    labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(# xlim=c(as.Date("2020-09-01"),as.Date("2021-02-01")), 
    xlim=c(as.Date("2020-10-01"),as.Date("2021-03-01")), 
    ylim=c(0,100), expand=c(0,0)) +
  # scale_color_manual("", values=reg_cols) +
  # scale_fill_manual("", values=reg_cols) +
  # scale_color_discrete("", h=c(0, 280), c=200) +
  # scale_fill_discrete("", h=c(0, 280), c=200) +
  geom_point(data=data_denmark_whole, 
             aes(x=date, y=propB117*100, size=total # ,
                 # colour=Region
             ), 
             colour=I("steelblue"), 
             alpha=I(0.5)) +
  scale_size_continuous("total n", trans="sqrt", 
                        range=c(1, 4), limits=c(1,max(data_denmark_whole$total)), 
                        breaks=c(10,100,1000)) +
  # guides(fill=FALSE) + 
  # guides(colour=FALSE) + 
  theme(legend.position = "right") +
  xlab("Collection date") +
  ggtitle("DENMARK") +
  theme(plot.title = element_text(hjust = 0.5))
plot_denmark_response



# 5.2. DATA SWITZERLAND ####

# Data source: https://ispmbern.github.io/covid-19/variants (contact: Christian Althaus) 
# & https://github.com/covid-19-Re/variantPlot/raw/master/data/data.csv (https://ibz-shiny.ethz.ch/covidDashboard/variant-plot/index.html, contact: Tanja Stadler)
# data was downloaded on the 17th of February 2021

data_geneva = read.csv("https://ispmbern.github.io/covid-19/variants/data/variants_GE.csv")
data_geneva$date = as.Date(data_geneva$date)
data_geneva$lab = "Geneva"
colnames(data_geneva)[colnames(data_geneva) %in% c("N501Y")] = c("n_B117")
head(data_geneva)
data_zurich = read.csv("https://ispmbern.github.io/covid-19/variants/data/variants_ZH.csv")
data_zurich$date = as.Date(data_zurich$date)
data_zurich$lab = "Zürich"
colnames(data_zurich)[colnames(data_zurich) %in% c("N501Y")] = c("n_B117")
head(data_zurich)
data_bern = read.csv("https://ispmbern.github.io/covid-19/variants/data/variants_BE.csv")
data_bern$date = as.Date(data_bern$date)
data_bern$lab = "Bern"
colnames(data_bern)[colnames(data_bern) %in% c("N501Y")] = c("n_B117")
head(data_bern)

data_viollier_risch = read.csv("https://github.com/covid-19-Re/variantPlot/raw/master/data/data.csv")
data_viollier_risch[is.na(data_viollier_risch)] = 0
data_viollier_risch$date = as.Date(NA)
data_viollier_risch$date[data_viollier_risch$week>=51] = lubridate::ymd( "2020-01-01" ) + 
  lubridate::weeks( data_viollier_risch$week[data_viollier_risch$week>=51] - 1 ) + 1
data_viollier_risch$date[data_viollier_risch$week<51] = lubridate::ymd( "2021-01-01" ) + 
  lubridate::weeks( data_viollier_risch$week[data_viollier_risch$week<51] - 1 ) + 6 # PS dates were made to match the ones given in https://ispmbern.github.io/covid-19/variants/data/variants_CH.csv
colnames(data_viollier_risch)[colnames(data_viollier_risch) %in% c("n","b117")] = c("total","n_B117")
data_viollier_risch = data_viollier_risch[,c("date","total","n_B117","lab")]

data_switzerland = rbind(data_geneva, data_zurich, data_bern, data_viollier_risch)[,c("date","lab","n_B117","total")]
# write_csv(data_switzerland, file=".//multinomial_logistic_fits//data//ch//data_switzerland.csv")

# Details data:
# Viollier data = sequencing of a random subset of all positive cases by ETH/Tanja Stadler (covers large parts of Switzerland, though with a bias towards German speaking Switzerland) - as it is sequencing data is 1-2 weeks later than N501Y screening
# Risch - Taqpath + N501Y re-screening = faster (covers primarily German speaking Switzerland)
# Samples are provided and screened by Labor Risch. Genomic characterization is performed by Labor Risch, the University Hospital Basel (Clinical Mircobiology) and the University Hospitals of Geneva (Group Eckerle and Group Kaiser). 
# Geneva - centre de reference pour infections virales emergentes / university hospital Geneva - N501Y and WGS currently:
# Samples that were sent to the Geneva University Hospitals for primary diagnosis of SARS-CoV-2. All positives were re-screened for 501Y using RT-PCR (mostly B.1.1.7). To cover the period of November and December 2020, we use sequence data from randomly chosen samples from Geneva that were submitted to GISAID by the Swiss Viollier Sequencing Consortium from ETH Zurich.
# Bern: Samples from SARS-CoV-2-positive cases that were re-screened for 501Y using RT-PCR at the Institute for Infectious Diseases, University of Bern.
# Zurich: Samples from SARS-CoV-2-positive cases from the University Hospital Zurich and test centers at Limmattal Hospital in Schlieren (ZH) and Spital Männedorf that were re-screened for 501Y using RT-PCR at the Institute of Medical Virology, University of Zurich. In addition, we use SARS-CoV-2-positive samples from Kantonsspital Winterthur and its walk-in test center that were re-screened for 501Y using RT-PCR.

data_switzerland = read_csv(file=".//multinomial_logistic_fits//data//ch//data_switzerland.csv", col_names=TRUE) # colClasses=c("Date","character","numeric","numeric")
data_switzerland = data.frame(data_switzerland)
data_switzerland$date = as.Date(data_switzerland$date)
data_switzerland$lab = factor(data_switzerland$lab, levels=c("Geneva","Zürich","Bern","Viollier","Risch"),
                              labels=c("Geneva","Zürich","Bern","Switzerland","Switzerland"))
data_switzerland$date_num = as.numeric(data_switzerland$date)
data_switzerland$obs = factor(1:nrow(data_switzerland))
data_switzerland$propB117 = data_switzerland$n_B117 / data_switzerland$total

fit_switerland1 = glmer(cbind(n_B117,total-n_B117) ~ (1|obs) + lab + scale(date_num), family=binomial(logit), data=data_switzerland)
fit_switerland2 = glmer(cbind(n_B117,total-n_B117) ~ (1|obs) + lab * scale(date_num), family=binomial(logit), data=data_switzerland)
fit_switerland3 = glmer(cbind(n_B117,total-n_B117) ~ (1|lab/obs) + scale(date_num), family=binomial(logit), data=data_switzerland)
fit_switerland4 = glmer(cbind(n_B117,total-n_B117) ~ (date_num||lab/obs) + scale(date_num), family=binomial(logit), data=data_switzerland)
BIC(fit_switerland1, fit_switerland2, fit_switerland3, fit_switerland4)
#                 df      BIC
# fit_switerland1  6 483.3556
# fit_switerland2  9 497.0946
# fit_switerland3  4 492.7274
# fit_switerland4  6 502.0293

# fit fit_switerland1 best

summary(fit_switerland1)


#  GROWTH RATE & TRANSMISSION ADVANTAGE

# on average across all regions, using the most parsimonious model fit_switerland1, we get
ch_growthrates_avg_B117vsallother = as.data.frame(emtrends(fit_switerland1, ~ 1, var="date_num"))[,-c(3,4)] 
colnames(ch_growthrates_avg_B117vsallother)[2] = "logistic_growth_rate"
ch_growthrates_avg_B117vsallother = M.from.delta_r_df(ch_growthrates_avg_B117vsallother)
ch_growthrates_avg_B117vsallother
# 1         logistic_growth_rate  asymp.LCL  asymp.UCL     M1   M1.LCL   M1.UCL       M2   M2.LCL   M2.UCL
# 1 overall            0.1006605 0.09215089 0.1091702 1.739561 1.66002 1.822914 1.436742 1.393395 1.481437
table2csv(ch_growthrates_avg_B117vsallother, file=".\\multinomial_logistic_fits\\tables\\model 3b_B1177vsallother_bGLMM_CH.csv")



# PLOT MODEL FIT
date.from = as.numeric(as.Date("2020-09-01"))
date.to = as.numeric(as.Date("2021-03-01")) # date to extrapolate to
total.SD = sqrt(sum(sapply(as.data.frame(VarCorr(fit_switerland1))$sdcor, function (x) x^2))) 
# bias correction for random effects in marginal means, see https://cran.r-project.org/web/packages/emmeans/vignettes/transformations.html#bias-adj
fit_switzerland_preds = as.data.frame(emmeans(fit_switerland1, ~ date_num, 
                                              by=c("lab"), 
                                              at=list(date_num=seq(date.from,
                                                                   date.to)), 
                                              type="response"), bias.adjust = TRUE, sigma = total.SD)
fit_switzerland_preds$date = as.Date(fit_switzerland_preds$date_num, origin="1970-01-01")

n = length(levels(fit_switzerland_preds$lab))
reg_cols = hcl(h = seq(290, 0, length = n + 1), l = 50, c = 255)[1:n]
# reg_cols[2:n] = rev(reg_cols[2:n])

# PLOT MODEL FIT (logit scale):
plot_switzerland = qplot(data=fit_switzerland_preds, x=date, y=prob, geom="blank") +
  geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL, 
                  fill=lab
  ), 
  # fill=I("steelblue"), 
  alpha=I(0.3)) +
  geom_line(aes(y=prob, 
                colour=lab
  ), 
  # colour=I("steelblue"), 
  alpha=I(0.8)) +
  ylab("Relative abundance of B.1.1.7 (%)") +
  theme_hc() + 
  xlab("") + 
  # scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
  #                   labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M")) +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(# xlim=c(as.Date("2020-09-01"),as.Date("2021-02-01")), 
    xlim=c(as.Date("2020-11-01"),as.Date("2021-03-01")), 
    ylim=c(0.001,99.9), expand=c(0,0)) +
  scale_color_manual("", values=reg_cols) +
  scale_fill_manual("", values=reg_cols) +
  # scale_color_discrete("", h=c(0, 280), c=200) +
  # scale_fill_discrete("", h=c(0, 280), c=200) +
  geom_point(data=data_switzerland, 
             aes(x=date, y=propB117, size=total,
                 colour=lab
             ), 
             # colour=I("steelblue"), 
             alpha=I(0.5)) +
  scale_size_continuous("total n", trans="sqrt", 
                        range=c(1, 4), limits=c(1,2000), breaks=c(10,100,1000)) +
  # guides(fill=FALSE) + 
  # guides(colour=FALSE) + 
  theme(legend.position = "right") +
  xlab("Collection date") +
  ggtitle("SWITZERLAND") +
  theme(plot.title = element_text(hjust = 0.5))
plot_switzerland


# PLOT MODEL FIT (response scale):
plot_switzerland_response = qplot(data=fit_switzerland_preds, x=date, y=prob*100, geom="blank") +
  geom_ribbon(aes(y=prob*100, ymin=asymp.LCL*100, ymax=asymp.UCL*100, colour=NULL, 
                  fill=lab
  ), 
  # fill=I("steelblue"), 
  alpha=I(0.3)) +
  geom_line(aes(y=prob*100, 
                colour=lab
  ), 
  # colour=I("steelblue"), 
  alpha=I(0.8)) +
  ylab("Relative abundance of B.1.1.7 (%)") +
  theme_hc() + 
  xlab("") + 
  # scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
  #                   labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M")) +
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                    labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(# xlim=c(as.Date("2020-09-01"),as.Date("2021-02-01")), 
    xlim=c(as.Date("2020-11-01"),as.Date("2021-03-01")), 
    ylim=c(0,100), expand=c(0,0)) +
  scale_color_manual("", values=reg_cols) +
  scale_fill_manual("", values=reg_cols) +
  # scale_color_discrete("", h=c(0, 280), c=200) +
  # scale_fill_discrete("", h=c(0, 280), c=200) +
  geom_point(data=data_switzerland, 
             aes(x=date, y=propB117*100, size=total,
                 colour=lab
             ), 
             # colour=I("steelblue"), 
             alpha=I(0.5)) +
  scale_size_continuous("total n", trans="sqrt", 
                        range=c(1, 4), limits=c(1,2000), breaks=c(500,1000,2000)) +
  # guides(fill=FALSE) + 
  # guides(colour=FALSE) + 
  theme(legend.position = "right") +
  xlab("Collection date") +
  ggtitle("SWITZERLAND") +
  theme(plot.title = element_text(hjust = 0.5))
plot_switzerland_response




# 5.3. S-GENE TARGET FAILURE DATA USA ####

# Data source: Helix® COVID-19 Surveillance, https://github.com/myhelix/helix-covid19db
# see preprint https://www.medrxiv.org/content/10.1101/2021.02.06.21251159v1 & https://github.com/andersen-lab/paper_2021_early-b117-usa/tree/master/b117_frequency/data

# Data were downloaded on the 17th of February 2021 

us_data = read.csv("https://github.com/myhelix/helix-covid19db/raw/master/counts_by_state.csv")
# write.csv(us_data, file=".//multinomial_logistic_fits//data//us//data_us.csv", row.names=F)

us_data = read.csv(file=".//multinomial_logistic_fits//data//us//data_us.csv")

# helix_b117 = read_tsv("https://github.com/andersen-lab/paper_2021_early-b117-usa/raw/master/b117_frequency/data/covid_baseline_for_b117_paper.20210127_update.txt") %>%
#   dplyr::select(state, collection_date, n_b117, n_sgtf_seq) # n_b117/n_sgtf_seq = prop of S dropout samples that are B117
# 
# helix_sgtf = read_tsv("https://github.com/andersen-lab/paper_2021_early-b117-usa/raw/master/b117_frequency/data/covid_baseline_for_b117_paper.20210201_klados20211029_phyloseq.txt") %>%
#   dplyr::select(state, collection_date, n, n_sgtf) # n_sgtf/n = prop of pos tests that have S dropout
# helix_sgtf = helix_sgtf[helix_sgtf$state %in% unique(helix_b117$state),]
# 
# helix_metadata = dplyr::left_join(helix_sgtf, helix_b117, by=c("state", "collection_date"))
# 
# tmp = helix_metadata %>%
#   dplyr::group_by(collection_date) %>%
#   dplyr::summarise(n_sgtf = sum(n_sgtf), n = sum(n)) %>%
#   dplyr::mutate(state = "USA")
# 
# tmp = dplyr::bind_rows(tmp, helix_metadata)
# states_gt_500 = tmp %>% dplyr::group_by(state) %>% dplyr::summarise(n = sum(n), n_sgtf = sum(n_sgtf)) %>% dplyr::filter(n > 500 & n_sgtf > 0) %>% dplyr::select(state) %>% purrr::as_vector()
# states_gt_500
# # state1  state2  state3  state4  state5  state6  state7  state8  state9 state10 state11 state12 state13 state14 state15 state16 state17 
# # "AL"    "AZ"    "CA"    "FL"    "GA"    "IL"    "IN"    "LA"    "MA"    "MI"    "MN"    "NC"    "NJ"    "NY"    "OH"    "PA"    "TX" 
# # state18 
# # "USA" 



us_data$collection_date = as.Date(us_data$collection_date)
us_data$collection_date_num = as.numeric(us_data$collection_date)
us_data$obs = factor(1:nrow(us_data))
# us_data = us_data[us_data$state %in% sel_states,]
us_data$state = factor(us_data$state)


fit_us_propB117amongSGTF = glmer(cbind(B117, sequenced_SGTF-B117) ~ (1|state)+scale(collection_date_num), 
                                                                    family=binomial(logit), data=us_data)

# implied growth rate advantage of B.1.1.7 over other earlier strains showing S dropout:
as.data.frame(emtrends(fit_us_propB117amongSGTF, ~ 1, var="collection_date_num"))[,c(2,5,6)]
#   collection_date_num.trend  asymp.LCL asymp.UCL
# 1                0.09623626 0.07866956  0.113803

# with a generation time of 5.5 days this would translate to a multiplicative effect on Rt
# and estimated increased infectiousness of B.1.1.7 over other strains showing S dropout of
exp(5.5*as.data.frame(emtrends(fit_us_propB117amongSGTF, ~ 1, var="collection_date_num"))[,c(2,5,6)])
#   collection_date_num.trend asymp.LCL asymp.UCL
# 1                   1.697742  1.541387  1.869958


# FIT FOR WHOLE US + PLOT

fitted_truepos = predict(fit_us_propB117amongSGTF, newdat=us_data, type="response") 
# fitted true positive rate, ie prop of S dropout samples that are B.1.1.7 for dates & states in helix_sgtf

us_data$est_n_B117 = us_data$all_SGTF*fitted_truepos # estimated nr of B.1.1.7 samples
us_data$propB117 = us_data$est_n_B117/us_data$positive
fit_us1 = glmer(cbind(est_n_B117, positive-est_n_B117) ~ (1|state/obs)+scale(collection_date_num), 
               family=binomial(logit), data=us_data) # random intercepts by state
fit_us2 = glmer(cbind(est_n_B117, positive-est_n_B117) ~ (collection_date_num||state/obs)+scale(collection_date_num), 
                family=binomial(logit), data=us_data) # random intercepts+slopes by state, with uncorrelated intercepts & slopes
BIC(fit_us1, fit_us2) # random intercept model fit_us1 is best
# df      BIC
# fit_us1  4 1622.360
# fit_us2  6 1638.224
summary(fit_us1)

#  GROWTH RATE & TRANSMISSION ADVANTAGE

# on average across all states, using the most parsimonious model fit_us1, we get
us_growthrates_avg_B117vsallother = as.data.frame(emtrends(fit_us1, ~ 1, var="collection_date_num"))[,-c(3,4)] 
colnames(us_growthrates_avg_B117vsallother)[2] = "logistic_growth_rate"
us_growthrates_avg_B117vsallother = M.from.delta_r_df(us_growthrates_avg_B117vsallother)
us_growthrates_avg_B117vsallother
# 1         logistic_growth_rate  asymp.LCL  asymp.UCL     M1   M1.LCL   M1.UCL       M2   M2.LCL   M2.UCL
# 1 overall            0.08433614 0.08033021 0.08834206 1.590182 1.55553 1.625607 1.354741 1.335344 1.374419
table2csv(us_growthrates_avg_B117vsallother, file=".\\multinomial_logistic_fits\\tables\\model 3c_B1177vsallother_bGLMM_USA.csv")


# plot model fit fit_us

date.to = as.numeric(as.Date("2021-06-01"))
# sel_states = intersect(rownames(ranef(fit_us)$state)[order(ranef(fit_us1)$state[,1], decreasing=T)],states_gt_500)[1:16] # unique(helix_sgtf$state[helix_sgtf$propB117>0.03])
# rem_states = c("NY","NJ","MN","IL","AL","OH","MI") # states with too few data points we don't want to show on plot
# sel_states = setdiff(sel_states,rem_states)


# sel_states = unique(us_data$state)
# we fitted our model on all the available data from all states, but below we will plot just
# the 9 states with the most data
# sel_states=c("FL","NY","CA","NJ","GA","TX","OH","PA","LA","IL","MI","MA","NC","IN","AZ")
# sel_states=c("FL","CA","GA","TX","PA","LA","IL","MI","MA","NC","IN","AZ")
sel_states=c("FL","CA","GA","TX","PA","MA","NC","IN","AZ")
total.SD = sqrt(sum(sapply(as.data.frame(VarCorr(fit_us))$sdcor, function (x) x^2))) 
fit_us_preds = as.data.frame(emmeans(fit_us1, ~ collection_date_num, 
                                     # by="state", 
                                     at=list(collection_date_num=seq(min(us_data$collection_date_num),
                                                                     date.to)), 
                                     type="link"), bias.adjust = TRUE, sigma = total.SD)
fit_us_preds$collection_date = as.Date(fit_us_preds$collection_date_num, origin="1970-01-01")
fit_us_preds2 = do.call(rbind,lapply(unique(us_data$state), function(st) { ranintercs = ranef(fit_us1)$state
raninterc = ranintercs[rownames(ranintercs)==st,]
data.frame(state=st, fit_us_preds, raninterc=raninterc)}))
fit_us_preds2$prob = plogis(fit_us_preds2$emmean+fit_us_preds2$raninterc)
fit_us_preds2$prob.asymp.LCL = plogis(fit_us_preds2$asymp.LCL+fit_us_preds2$raninterc)
fit_us_preds2$prob.asymp.UCL = plogis(fit_us_preds2$asymp.UCL+fit_us_preds2$raninterc)
fit_us_preds2 = fit_us_preds2[as.character(fit_us_preds2$state) %in% sel_states,]
fit_us_preds2$state = droplevels(fit_us_preds2$state)
fit_us_preds2$state = factor(fit_us_preds2$state, # we order states by random intercept, ie date of introduction
                             levels=intersect(rownames(ranef(fit_us1)$state)[order(ranef(fit_us1)$state[,1], decreasing=T)],
                                              sel_states))

# PLOT MODEL FIT (logit scale)
plot_us = qplot(data=fit_us_preds2, x=collection_date, y=prob, geom="blank") +
  facet_wrap(~state, nrow=3) +
  geom_ribbon(aes(y=prob, ymin=prob.asymp.LCL, ymax=prob.asymp.UCL, colour=NULL, 
                  fill=state
  ), 
  # fill=I("steelblue"), 
  alpha=I(0.3)) +
  geom_line(aes(y=prob, 
                colour=state
  ), 
  # colour=I("steelblue"), 
  alpha=I(0.8)) +
  ylab("Relative abundance of B.1.1.7 (%)") +
  theme_hc() + xlab("") + 
  # scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
  #                   labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M")) +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(xlim=c(min(fit_us_preds$collection_date), as.Date("2021-04-01")), 
                  # xlim=c(as.Date("2020-07-01"),as.Date("2021-01-31")), 
                  ylim=c(0.001,0.9990001), expand=c(0,0)) +
  scale_color_discrete("state", h=c(0, 240), c=120, l=50) +
  scale_fill_discrete("state", h=c(0, 240), c=120, l=50) +
  geom_point(data=us_data[us_data$state %in% sel_states,],  
             aes(x=collection_date, y=propB117, size=positive,
                 colour=state
             ), pch=I(16),
             # colour=I("steelblue"), 
             alpha=I(0.3)) +
  scale_size_continuous("number of\npositive tests", trans="sqrt", 
                        range=c(1, 4), limits=c(1,max(us_data$positive)), breaks=c(10,100,1000)) +
  # guides(fill=FALSE) + 
  # guides(colour=FALSE) + 
  theme(legend.position = "right") +
  xlab("Collection date") +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5)) +
  ggtitle("US") +
  theme(plot.title = element_text(hjust = 0.5))
plot_us


# PLOT MODEL FIT (response scale)
plot_us_response = qplot(data=fit_us_preds2, x=collection_date, y=prob*100, geom="blank") +
  facet_wrap(~state) +
  geom_ribbon(aes(y=prob*100, ymin=prob.asymp.LCL*100, ymax=prob.asymp.UCL*100, colour=NULL, 
                  fill=state
  ), 
  # fill=I("steelblue"), 
  alpha=I(0.3)) +
  geom_line(aes(y=prob*100, 
                colour=state
  ), 
  # colour=I("steelblue"), 
  alpha=I(0.8)) +
  ylab("Relative abundance of B.1.1.7 (%)") +
  theme_hc() + xlab("") + 
  # scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
  #                   labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M")) +
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                    labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  coord_cartesian(xlim=c(min(fit_us_preds2$collection_date), as.Date("2021-04-01")), 
                  # xlim=c(as.Date("2020-07-01"),as.Date("2021-01-31")), 
                  ylim=c(0,100), expand=c(0,0)) +
  scale_color_discrete("state", h=c(0, 240), c=120, l=50) +
  scale_fill_discrete("state", h=c(0, 240), c=120, l=50) +
  geom_point(data=us_data[us_data$state %in% sel_states,],  
             aes(x=collection_date, y=propB117*100, size=positive,
                 colour=state
             ), pch=I(16),
             # colour=I("steelblue"), 
             alpha=I(0.3)) +
  scale_size_continuous("number of\npositive tests", trans="log10", 
                        range=c(1, 4), limits=c(1,max(us_data$positive)), breaks=c(10,100,1000)) +
  # guides(fill=FALSE) + 
  # guides(colour=FALSE) + 
  theme(legend.position = "right") +
  xlab("Collection date") +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5)) +
  ggtitle("US") +
  theme(plot.title = element_text(hjust = 0.5))

plot_us_response






# 5.4. MULTIPANEL PLOT INTERNATIONAL COMPARISONS ####

fit_uk_preds2 = fit_ukSGTF_4_preds
fit_uk_preds2$country = "UK"
colnames(fit_uk_preds2)[2] = "REGION"
colnames(fit_uk_preds2)[1] = "date_num"
colnames(fit_uk_preds2)[8] = "date"
fit_switzerland_preds2 = fit_switzerland_preds
fit_switzerland_preds2$country = "Switzerland"
colnames(fit_switzerland_preds2)[2] = "REGION"
colnames(fit_switzerland_preds2)[1] = "date_num"
colnames(fit_switzerland_preds2)[8] = "date"
fit_denmark_preds2 = fit_denmark_preds
fit_denmark_preds2$country = "Denmark"
fit_denmark_preds2$REGION = "Denmark"
colnames(fit_denmark_preds2)[1] = "date_num"
colnames(fit_denmark_preds2)[7] = "date"
fit_us_preds3 = fit_us_preds2
fit_us_preds3$country = "USA"
fit_us_preds3 = fit_us_preds3[,-which(colnames(fit_us_preds3) %in% c("asymp.LCL","asymp.UCL"))]
colnames(fit_us_preds3)[1] = "REGION"
colnames(fit_us_preds3)[2] = "date_num"
colnames(fit_us_preds3)[6] = "date"
colnames(fit_us_preds3)[9] = "asymp.LCL"
colnames(fit_us_preds3)[10] = "asymp.UCL"
fit_us_preds3 = fit_us_preds3[fit_us_preds3$REGION %in% c("FL","CA"),]
fit_us_preds3$REGION = factor(fit_us_preds3$REGION, levels=c("FL","CA"), labels=c("Florida","California"))
fit_us_preds3 = fit_us_preds3[,c("date_num","REGION","prob","SE","df","asymp.LCL","asymp.UCL","date","country")]

fits_international = rbind(fit_uk_preds2,fit_denmark_preds2,fit_switzerland_preds2,fit_us_preds3)
fits_international$country = factor(fits_international$country, levels=c("UK","Denmark","Switzerland","USA"))

sgtfdata_uk2 = sgtfdata_uk
sgtfdata_uk2$country = "UK"
colnames(sgtfdata_uk2)[colnames(sgtfdata_uk2) %in% c("collection_date","n_pos")] = c("date","total")
sgtfdata_uk2 = sgtfdata_uk2[,c("date","country","REGION","propB117","total")]

data_switzerland2 = data_switzerland
data_switzerland2$country = "Switzerland"
colnames(data_switzerland2)[colnames(data_switzerland2) %in% c("lab")] = c("REGION")
data_switzerland2 = data_switzerland2[,c("date","country","REGION","propB117","total")]

data_denmark2 = data_denmark_whole
data_denmark2$country = "Denmark"
data_denmark2$REGION = "Denmark"
data_denmark2 = data_denmark2[,c("date","country","REGION","propB117","total")]

data_us2 = data.frame(us_data)
data_us2$country = "USA"
colnames(data_us2)[1] = "REGION"
colnames(data_us2)[2] = "date"
colnames(data_us2)[3] = "total"
data_us2 = data_us2[,c("date","country","REGION","propB117","total")]
data_us2 = data_us2[data_us2$REGION %in% c("FL","CA"),]
data_us2$REGION = factor(data_us2$REGION, levels=c("FL","CA"), labels=c("Florida","California"))

data_international = rbind(sgtfdata_uk2, data_denmark2, data_switzerland2, data_us2)
data_international$country = factor(data_international$country, levels=c("UK","Denmark","Switzerland","USA"))

# n1 = length(levels(fit_uk_preds2$REGION))
# n2 = length(levels(fit_switzerland_preds2$REGION))
# n3 = length(levels(fit_denmark_preds2$REGION))
# reg_cols = c(hcl(h = seq(290, 0, length = n1), l = 50, c = 255),
#              muted(hcl(h = seq(290, 0, length = n2+n3), l = 50, c = 255), c=200, l=40))

# ymin = 0.001
ymax = 0.999
data_international$propB117[data_international$propB117>ymax] = ymax
fits_international$prob[fits_international$prob>ymax] = ymax
fits_international$asymp.LCL[fits_international$asymp.LCL>ymax] = ymax
fits_international$asymp.UCL[fits_international$asymp.UCL>ymax] = ymax

fits_international$REGION = factor(fits_international$REGION, levels=levels(fits_international$REGION))
data_international$REGION = factor(data_international$REGION, levels=levels(fits_international$REGION))

# PLOT MODEL FITS (response scale)
plot_international = qplot(data=fits_international, x=date, y=prob, geom="blank") +
  facet_wrap(~country, ncol=1, scales="fixed") +
  geom_ribbon(aes(y=prob, ymin=asymp.LCL, ymax=asymp.UCL, colour=NULL, 
                  fill=REGION
  ), 
  # fill=I("steelblue"), 
  alpha=I(0.3)) +
  geom_line(aes(y=prob, 
                colour=REGION
  ), 
  # colour=I("steelblue"), 
  alpha=I(0.8)) +
  ylab("Relative abundance of B.1.1.7 (%)") +
  theme_hc() + 
  xlab("") + 
  # scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
  #                   labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M")) +
  scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
                      labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9") # ,
                      # limits = c(ymin,ymax+1E-7)
                      ) +
  # scale_color_manual("", values=reg_cols) +
  # scale_fill_manual("", values=reg_cols) +
  scale_color_discrete("region", h=c(0, 280), c=155, l=50) +
  scale_fill_discrete("region", h=c(0, 280), c=155, l=50) +
  # scale_color_discrete("", h=c(0, 280), c=200) +
  # scale_fill_discrete("", h=c(0, 280), c=200) +
  geom_point(data=data_international, 
             aes(x=date, y=propB117, size=total, shape=country,
                 colour=REGION, fill=REGION
             ), 
             # colour=I("steelblue"), 
             alpha=I(0.5)) +
  scale_size_continuous("total n", trans="sqrt", 
                        range=c(1, 2), limits=c(1,max(data_international$total)), breaks=c(10,100,1000,10000)) +
  scale_shape_manual(values=21:25) +
  # guides(fill=FALSE) + 
  # guides(colour=FALSE) + 
  theme(legend.position = "right") +
  xlab("Collection date") +
  guides(
    shape = guide_legend(order = 1),
    color = guide_legend(order = 2),
    fill = guide_legend(order = 2),
    size = guide_legend(order = 3)
  ) + 
  coord_cartesian( 
    xlim=c(as.Date("2020-09-01"),as.Date("2021-03-01")),
    ylim=c(ymin,ymax+1E-7), 
    expand=FALSE)
# ggtitle("INTERNATIONAL SPREAD OF SARS-CoV2 VARIANT B.1.1.7") +
# theme(plot.title = element_text(hjust = 0.5))
plot_international

saveRDS(plot_international, file = ".\\multinomial_logistic_fits\\plots\\FigS6_binomGLMM_B117vsall_fits_UK_DK_CH_USA.rds")
graph2ppt(file = ".\\multinomial_logistic_fits\\plots\\FigS6_binomGLMM_B117vsall_fits_UK_DK_CH_USA.pptx", width=7, height=8)
ggsave(file = ".\\multinomial_logistic_fits\\plots\\FigS6_binomGLMM_B117vsall_fits_UK_DK_CH_USA.png", width=7, height=8)
ggsave(file = ".\\multinomial_logistic_fits\\plots\\FigS6_binomGLMM_B117vsall_fits_UK_DK_CH_USA.pdf", width=7, height=8)




# PLOT MODEL FITS (response scale)
plot_international_response = qplot(data=fits_international, x=date, y=prob*100, geom="blank") +
  facet_wrap(~country, ncol=1) +
  geom_ribbon(aes(y=prob*100, ymin=asymp.LCL*100, ymax=asymp.UCL*100, colour=NULL, 
                  fill=REGION
  ), 
  # fill=I("steelblue"), 
  alpha=I(0.3)) +
  geom_line(aes(y=prob*100, 
                colour=REGION
  ), 
  # colour=I("steelblue"), 
  alpha=I(0.8)) +
  ylab("Relative abundance of B.1.1.7 (%)") +
  theme_hc() + 
  xlab("") + 
  # scale_x_continuous(breaks=as.Date(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01","2021-02-01","2021-03-01")),
  #                   labels=c("M","A","M","J","J","A","S","O","N","D","J","F","M")) +
  # scale_y_continuous( trans="logit", breaks=c(10^seq(-5,0),0.5,0.9,0.99,0.999),
  #                    labels = c("0.001","0.01","0.1","1","10","100","50","90","99","99.9")) +
  scale_color_discrete("region", h=c(0, 280), c=155, l=50) +
  scale_fill_discrete("region", h=c(0, 280), c=155, l=50) +
  #   scale_color_manual("", values=reg_cols) +
  #  scale_fill_manual("", values=reg_cols) +
  # scale_color_discrete("", h=c(0, 280), c=200) +
  # scale_fill_discrete("", h=c(0, 280), c=200) +
  geom_point(data=data_international, 
             aes(x=date, y=propB117*100, size=total, shape=country,
                 colour=REGION, fill=REGION
             ), 
             # colour=I("steelblue"), 
             alpha=I(0.5)) +
  scale_size_continuous("total n", trans="sqrt", 
                        range=c(1, 2), limits=c(1,max(data_international$total)), breaks=c(10,100,1000,10000)) +
  scale_shape_manual(values=21:25) +
  # guides(fill=FALSE) + 
  # guides(colour=FALSE) + 
  theme(legend.position = "right") +
  xlab("Collection date") +
  guides(
    shape = guide_legend(order = 1),
    color = guide_legend(order = 2),
    fill = guide_legend(order = 2),
    size = guide_legend(order = 3)
  ) +
  coord_cartesian( 
    xlim=c(as.Date("2020-09-01"),as.Date("2021-03-01")),
    ylim=c(0,100), expand=c(0,0))
# +
# ggtitle("INTERNATIONAL SPREAD OF SARS-CoV2 VARIANT B.1.1.7") +
# theme(plot.title = element_text(hjust = 0.5))
plot_international_response

saveRDS(plot_international_response, file = paste0(".\\plots\\",dat,"\\binomGLMM_B117vsall_fits_UK_DK_CH_USA_response.rds"))
graph2ppt(file=paste0(".\\plots\\",dat,"\\binomGLMM_B117vsall_fits_UK_DK_CH_USA_response.pptx"), width=7, height=8)
ggsave(file=paste0(".\\plots\\",dat,"\\binomGLMM_B117vsall_fits_UK_DK_CH_USA_response.png"), width=7, height=8)
ggsave(file=paste0(".\\plots\\",dat,"\\binomGLMM_B117vsall_fits_UK_DK_CH_USA_response.pdf"), width=7, height=8)


plot_us2 = plot_us + coord_cartesian(xlim=c(as.Date("2020-11-01"), as.Date("2021-03-31")),
                                     ylim=c(0.001,99.9), expand=c(0,0)) # + ggtitle("SPREAD OF VARIANT B.1.1.7 IN THE US")
plot_us2

saveRDS(plot_us2, file = ".\\multinomial_logistic_fits\\plots\\FigS7_model3c_plot B1177vsall_fit bGLMM_US.rds")
graph2ppt(file = ".\\multinomial_logistic_fits\\plots\\FigS7_model3c_plot B1177vsall_fit bGLMM_US.pptx", width=8, height=6)
ggsave(file = ".\\multinomial_logistic_fits\\plots\\FigS7_model3c_plot B1177vsall_fit bGLMM_US.png", width=8, height=6)
ggsave(file = ".\\multinomial_logistic_fits\\plots\\FigS7_model3c_plot B1177vsall_fit bGLMM_US.pdf", width=8, height=6)
