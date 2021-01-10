library(data.table)
library(magrittr)
library(parallel)
library(KFAS)
library(zoo)
library(ggplot2)
library(ggpubr)

ncores=7


# r(t) estimation function ---------------------------------------------
nbss <- function(N){
    nb_ss_model <- function(N,dispersion){
        model <- SSModel(N~SSMtrend(2, Q=list(0, NA),
                                    P1=diag(c(10, 1)),
                                    a1=c(0, 0),
                                    state_names=c("abundance", "growth_rate"))+
                             SSMseasonal(7),
                         u=rep(exp(dispersion),length(N)),distribution='negative binomial')
        fit <- fitSSM(model,inits = 0,method='L-BFGS-B',control=list(maxit=200))
    }
    
    nb_log_likelihood <- function(N, disp){
        fit <- nb_ss_model(N, disp)
        ll <- logLik(fit$model, marginal = TRUE)
        return(-ll)
    }
    
    
    ### maximum marginal likelihood estimation
    mldisp <- optim(c(-1), function(disp) nb_log_likelihood(N, disp),
                    method="Brent", lower=-2, upper=2) 
    fit <- nb_ss_model(N, mldisp$par)
    
    
    estimates <- KFS(fit$model, smoothing="state")
    out <- data.table(q2.5_abundance= c(qnorm(0.025, estimates$alphahat[,1], (sqrt(estimates$V[1,1,])))), 
                      q97.5_abundance= c(qnorm(0.975, estimates$alphahat[,1],(sqrt(estimates$V[1,1,])))), 
                      abundance = (c(estimates$alphahat[,1])),
                      q2.5_growth_rate = c(qnorm(0.025, estimates$alphahat[,2], (sqrt(estimates$V[2,2,])))), 
                      q97.5_growth_rate = c(qnorm(0.975, estimates$alphahat[,2], (sqrt(estimates$V[2,2,])))),
                      q25_growth_rate = c(qnorm(0.25, estimates$alphahat[,2], (sqrt(estimates$V[2,2,])))), 
                      q75_growth_rate = c(qnorm(0.75, estimates$alphahat[,2], (sqrt(estimates$V[2,2,])))), 
                      growth_rate = (c(estimates$alphahat[,2])),
                      dispersion=mldisp$par[1],
                      sd_growth_rate = sqrt(estimates$V[2,2,]),
                      N=N)
    return(out)
}


# Data load + prep --------------------------------------------------------

COG <- read.csv('data/cog_metadata_microreact_public-2020-12-22-annotated.csv') %>%
    as.data.table()
COG[,sample_date:=as.Date(sample_date)]

COG[lineage == 'B.1.1.7' & n501y == 'Y' & del_21765_6 == 'del',lineage:='VOC']

###pooling data across all UK
UK_pooled <- COG[,list(N=.N),by=c('lineage','sample_date')]
UK_pooled[,n_obs:=.N,by=lineage]

##remove lineages observed fewer than 10 times
UK_pooled <- UK_pooled[n_obs>=10]

### filling dates

lineage_dates <- expand.grid('lineage'=unique(COG$lineage),
                             'sample_date'=unique(COG$sample_date)) %>%
    as.data.table

setkey(lineage_dates,lineage,sample_date)
setkey(UK_pooled,lineage,sample_date)

UK_pooled <- UK_pooled[lineage_dates]
UK_pooled[is.na(N),N:=0]

### remove leading zeros - start lineage when first observed three days in a week
UK_pooled[,weekly_observations:=rollapply(N,FUN=function(x) sum(x!=0),
                                          width=7,align='left',fill=NA),
          by=lineage]
UK_pooled[,start_date:=min(sample_date[weekly_observations>=3],na.rm=T),by=lineage]
UK_pooled <- UK_pooled[!is.infinite(start_date)]

### remove trailing zeroes after 2 weeks of no observation and no subsequent observation up to date of last sample
UK_pooled[,rev_cum_tot:=rev(cumsum(rev(N))),by=lineage]
UK_pooled[,last_date:=max(sample_date[rev_cum_tot!=0]),by=lineage]

UK_pooled <- UK_pooled[sample_date >= start_date & 
                           sample_date<=(last_date+14),c('lineage','sample_date','N')]


# estimating abundance & growth rates -------------------------------------

lineage_nbs <- function(lin,UK_pooled.=UK_pooled){
    out <- tryCatch(nbss(UK_pooled[lineage==lin,N]),error=function(e) NULL)
    if (!is.null(out)){
        out$lineage <- lin
        out[,sample_date:=UK_pooled[lineage==lin]$sample_date]
    }
    return(out)
}


lineages <- unique(UK_pooled$lineage)

cl <- makeCluster(ncores)
clusterEvalQ(cl,{library(data.table)
    library(KFAS)})
clusterExport(cl,varlist = c('nbss','UK_pooled'))

UK <- parLapply(cl,lineages,lineage_nbs) %>% rbindlist

stopCluster(cl)
rm('cl')
write.csv(UK,file = 'data/UK_lineage_rt_estimates_VOC_within_B117.csv')

# relativized growth rates ------------------------------------------------

UK[,relativized_rt:=(growth_rate-mean(growth_rate,na.rm=T))/sd(growth_rate,na.rm=T),by=sample_date]

IPOs <- UK[,list(initial_rt=mean(relativized_rt[1:31],na.rm=T),
                 sd_rt=sd(relativized_rt[1:31],na.rm=T),
                 start_date=min(sample_date)),by=lineage]
top_lineages <- c('B.1.177','VOC','B.1.1.1','B.1.1.301','B.1.1',"B.1.159",'B.1.1.7')
g_rel <- ggplot(IPOs,aes(start_date,initial_rt))+
    geom_point()+
    geom_point(data=IPOs[lineage %in% top_lineages],aes(col=lineage),cex=6)+
    theme_bw(base_size=15)+
    geom_hline(yintercept = 0)+
    scale_y_continuous('<rel(t)>')+
    geom_smooth(data=IPOs[!lineage %in% top_lineages],col='black',fill='black',alpha=0.2)+
    ggtitle("Average relativized Growth Rates 1 month post-IPO")+
    theme(legend.position = c(0.6,0.3))


# IPO analysis ------------------------------------------------------------

UK[,days_since_IPO:=1:.N,by=lineage]
g_t <- ggplot(UK,aes(days_since_IPO,relativized_rt,group=lineage))+
    geom_line(lwd=2,alpha=0.02)+
    geom_line(data=UK[lineage %in% top_lineages],aes(color=lineage),lwd=2)+
    theme_bw(base_size=15)+
    geom_hline(yintercept = 0)+
    scale_y_continuous('rel(t)')+
    theme(legend.position='none')+
    ggtitle('rel(t) dynamics post-IPO')

ggarrange(g_rel,g_t,ncol=2)


# Rt conversion -----------------------------------------------------------

growth_to_R <- function(r, gamma_mean, gamma_sd) {
    k <- (gamma_sd / gamma_mean)^2
    R <- (1 + k * r * gamma_mean)^(1 / k)
    return(R)
}