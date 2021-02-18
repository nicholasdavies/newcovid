library(data.table)
library(magrittr)
library(parallel)
library(KFAS)
library(zoo)
library(ggplot2)
library(ggpubr)

ncores <- 7
n_lagging_zeroes <- 7
min_obs_total <- 10
min_obs_of_IPO_week <- 3
max_date <- as.Date('2020-12-12')
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

COG <- read.csv('data/cog_metadata_microreact_public-2021-01-11-annotated.csv') %>%
    as.data.table()
COG[,sample_date:=as.Date(sample_date)]

COG[lineage == 'B.1.1.7' & n501y == 'Y' & del_21765_6 == 'del',lineage:='VOC']

###pooling data across all UK
UK_pooled <- COG[,list(N=.N),by=c('lineage','sample_date')]
UK_pooled[,n_obs:=.N,by=lineage]

##remove lineages observed fewer than 10 times
UK_pooled <- UK_pooled[n_obs>=min_obs_total]

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
UK_pooled[,start_date:=min(sample_date[weekly_observations>=min_obs_of_IPO_week],na.rm=T),by=lineage]
UK_pooled <- UK_pooled[!is.infinite(start_date)]

### remove trailing zeroes after 2 weeks of no observation and no subsequent observation up to date of last sample
UK_pooled[,rev_cum_tot:=rev(cumsum(rev(N))),by=lineage]
UK_pooled[,last_date:=max(sample_date[rev_cum_tot!=0]),by=lineage]

UK_pooled <- UK_pooled[sample_date >= start_date & 
                           sample_date<=(last_date+n_lagging_zeroes),c('lineage','sample_date','N')]


# estimating abundance & growth rates -------------------------------------

UK_pooled <- UK_pooled[sample_date<=max_date]

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
write.csv(UK,file = 'data/UK_lineage_rt_estimates_VOC_within_B117_Dec_12_stop_7_lagging_zeros.csv')

# relativized growth rates ------------------------------------------------
## we remove lineage B 1.1.7 because this lineage contains very few counts after our definition of the VOC subset in line 67 above.
UK[lineage!='B.1.1.7',relativized_rt:=(growth_rate-mean(growth_rate,na.rm=T))/sd(growth_rate,na.rm=T),by=sample_date]


IPO_length <- 31
IPOs <- UK[,list(initial_rt=mean(relativized_rt[1:IPO_length],na.rm=T),
                 sd_rt=sd(relativized_rt[1:IPO_length],na.rm=T),
                 start_date=min(sample_date)),by=lineage]
top_lineages <- c('B.1.177','B.1.177.7','VOC','B.1.1.1','B.1.1.301','B.1.1')
g_rel <- ggplot(IPOs[start_date<as.Date('2020-11-13')],aes(start_date,initial_rt))+
    geom_ribbon(aes(ymin=qnorm(0.025),ymax=qnorm(0.975)),fill='grey',col='black',alpha=0.2)+
    geom_point()+
    geom_point(data=IPOs[lineage %in% top_lineages],aes(col=lineage),cex=5)+
    geom_point(data=IPOs[lineage=='VOC'],aes(col=lineage),cex=8)+
    theme_bw(base_size=15)+
    geom_hline(yintercept = 0)+
    scale_x_date('IPO date',breaks=as.Date(c('2020-03-01','2020-05-01','2020-07-01','2020-09-01','2020-11-01')),
                       labels=c('March','May','July','Sept','Nov'))+
    scale_y_continuous(expression(paste('<',rho,'(t)>',sep='')),limits = c(-3,3))+
    geom_smooth(data=IPOs[!lineage %in% top_lineages & start_date<as.Date('2020-11-13')],col='black',fill='black',alpha=0.2)+
    ggtitle("Average relativized growth rates 1 month post-IPO")+
    guides(col=guide_legend(ncol=6))+
    theme(legend.position = 'bottom')

saveRDS(g_rel,file='avg_relativized_growth_rates_post_IPO.Rd')

# IPO analysis ------------------------------------------------------------

UK[,days_since_IPO:=1:.N,by=lineage]
g_t <- ggplot(UK,aes(days_since_IPO,relativized_rt))+
    geom_line(lwd=2,alpha=0.02,aes(group=lineage))+
    geom_line(data=UK[lineage %in% top_lineages],aes(color=lineage),lwd=2)+
    theme_bw(base_size=15)+
    geom_hline(yintercept = 0)+
    geom_smooth()+
    scale_y_continuous(expression(paste(rho,'(t)',sep='')),limits=c(-5,5))+
    scale_x_continuous('Days post-IPO')+
    theme(legend.position='none')+
    geom_ribbon(aes(ymin=-1.96,ymax=1.96),fill='grey',col='black',alpha=0.1)+
    # geom_vline(xintercept = 31,lty=2)+
    ggtitle('Relativized growth rate dynamics post-IPO')

saveRDS(g_t,'relativized_growth_rates_over_days_since_IPO.Rd')

save(list=ls(),file = 'relativized_growth_rate_workspace.Rd')

ggarrange(g_rel,g_t,nrow=2,labels = c("A","B"))
ggsave('~/COVID/COG_B117/figures/relativized_fitness_analysis_VOC_20201212.png',height=11,width=10,units='in')

# lineage_plot -----------------------------------------------------------

plot_lineage <- function(lin,UK.=UK){
    dum <- UK[lineage==lin]
    
    avg <- UK[lineage!=lin & sample_date %in% unique(dum$sample_date) & N!=0,
              list(growth_rate=mean(growth_rate,na.rm=T)),by=sample_date]
    g_N <- ggplot(dum,aes(sample_date,N))+
        geom_line(aes(sample_date,exp(abundance)),lwd=2)+
        geom_ribbon(aes(ymin=exp(q2.5_abundance),ymax=exp(q97.5_abundance)),fill='grey',alpha=0.4)+
        geom_point(cex=2)+
        theme_bw()+
        scale_y_continuous('Sequences')+
        ggtitle(paste(lin,'Abundance'))
    g_comp <- ggplot(dum,aes(sample_date,N/N_seq))+
        geom_point(cex=2)+
        geom_smooth()+
        theme_bw()+
        scale_y_continuous('Relative Abundance',limits=c(0,1))+
        ggtitle(paste(lin,'Relative Abundance'))
    g_r <- ggplot(dum,aes(sample_date,growth_rate))+
        geom_line(lwd=2)+
        geom_ribbon(aes(ymin=q2.5_growth_rate,ymax=q97.5_growth_rate),fill='grey',alpha=0.4)+
        geom_line(data=avg,col='steelblue',lwd=2)+
        theme_bw()+
        geom_hline(yintercept = 0)+
        scale_y_continuous('r(t)')+
        ggtitle(paste(lin,'Growth Rate'))
    ggarrange(g_N,g_comp,g_r,nrow=3,align='v')
}

plot_lineage('VOC')

plot_list <- lapply(top_lineages,plot_lineage)
ggarrange(plotlist = plot_list,ncol=7)


# lineage-richess dynamics ------------------------------------------------


g_N <- COG[,list(N=.N),by=sample_date] %>% 
    ggplot(aes(sample_date,N))+
    geom_bar(stat='identity',fill='grey',col='grey')+
    theme_bw(base_size=15)+
    geom_vline(xintercept = IPOs[lineage=='VOC',start_date],aes(col=lineage),lwd=2)+
    geom_vline(xintercept = as.Date('2020-12-12'),lwd=2,lty=2)+
    ggtitle('COG Sequences')+
    scale_x_continuous('')
g_S <- COG[,list(n_lineages=length(unique(lineage))),by=sample_date] %>% 
    ggplot(aes(sample_date,n_lineages))+
    theme_bw(base_size = 15)+
    geom_bar(stat='identity',fill='steelblue',col='steelblue',alpha=0.9)+
    ggtitle('COG Lineage Richness')+
    scale_y_continuous('Number of Lineages')+
    geom_vline(xintercept = IPOs[lineage=='VOC',start_date],aes(col=lineage),lwd=2)+
    geom_vline(xintercept = as.Date('2020-12-12'),lwd=2,lty=2)

ggarrange(g_N,g_S,nrow=2,align='v')

