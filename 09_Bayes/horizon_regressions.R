# read in data ####
load('horizon_regressions_output.Rdata')

# previous data processing ####
# load('horizon_preprocessed_data.Rdata')
# rm(quest_data)
# rm(trial_data)
# rm(trial_data_pre)

# trial_data_pre=read.csv(paste0(data_path,'horizon_trial_data_2021-04-13.csv'),
#                    stringsAsFactors=F,header=T)
# load(paste0(data_path,'scripts/preprocessing/excluded_subjects.Rdata'))
# horizon_exclude_IDs=unique(c(flag_horizon_IDs,flag_questionnaire_consistency,
#   demo_flag_IDs,dass_flag_IDs,oasis_flag_IDs,flag_grid_skips))
# trial_data=trial_data_pre[!(trial_data_pre$ID %in% horizon_exclude_IDs),]
# trial_data$valence=ifelse(trial_data$outcome_valenced>0,1,0)
# trial_data$valence=factor(trial_data$valence,labels=c('loss','gain'))
# trial_data$valence=relevel(trial_data$valence,ref="gain")
# trial_data$correct=ifelse(trial_data$valence=='gain',
#   ifelse(trial_data$Rvalue>trial_data$Lvalue,
#     ifelse(trial_data$Rchosen==1,1,0),
#     ifelse(trial_data$Rchosen==0,1,0)),
#   ifelse(trial_data$Rvalue>trial_data$Lvalue,
#     ifelse(trial_data$Rchosen==0,1,0),
#     ifelse(trial_data$Rchosen==1,1,0)))
# trial_data$free_trial=ifelse(trial_data$trial_block<5,NA,trial_data$trial_block-4)
# trial_data$times_Rchosen=NA
# trial_data$uncR=-1*trial_data$uncertainty #coded 1/-1 if R/L chosen more, should be -1/1
# for (t in 1:dim(trial_data)[1]) {
#   if (!is.na(trial_data$free_trial[t]) & trial_data$free_trial[t]==1) {
#     trial_data$times_Rchosen[t]=trial_data$Rchosen[t]
#   } else if (!is.na(trial_data$free_trial[t]) & trial_data$free_trial[t]>1) {
#     trial_data$times_Rchosen[t]=trial_data$times_Rchosen[t-1]+
#       trial_data$Rchosen[t]
#   }
# }
# trial_data$diff_Rchosen=ifelse(is.na(trial_data$free_trial),NA,
# 2*trial_data$times_Rchosen-trial_data$free_trial) # times R - (total - times R)
# trial_data$mean_diff=trial_data$Rvalue-trial_data$Lvalue
# trial_data$mean_diff_valenced=ifelse(trial_data$valence=='gain',trial_data$mean_diff,
#                                      -1*trial_data$mean_diff)
# trial_data$horizon=as.factor(trial_data$horizon)
# quest_data=read.csv(file=paste0(data_path,'merged_questionnaire_data_2021-04-13.csv'))
# trial_quest_data=merge(trial_data,quest_data,by.x='ID',by.y='subject')
# trial_quest_data$med_anx=as.factor(ifelse(trial_quest_data$dass_anx<14,0,1))
# trial_quest_data_ff=trial_quest_data[trial_quest_data$free_trial==1,]
# save.image('horizon_preprocessed_data.Rdata')

# plot data ####
library(ggplot2)
trial_quest_data_ff_avg=aggregate(Rchosen~horizon*valence*uncR,
                                 data=trial_quest_data_ff,FUN=mean)
ggplot(trial_quest_data_ff_avg,aes(y=Rchosen,x=as.factor(uncR),
                                  color=interaction(valence,horizon)))+
  geom_point(size=20,shape = '_')+
  labs(x='Relative Uncertainty of Right-sided Option',
       y='Proportion Right-sided Option Chosen')+
  theme_classic()+geom_hline(yintercept=0.5)

# frequentist analyses ####
library(lme4)
lm3=glmer(Rchosen~horizon*valence*(uncR+scale(mean_diff))+
  (horizon*valence*(uncR+scale(mean_diff))|ID),data=trial_quest_data_ff,
  family='binomial',
  control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
summary(lm3)

# bayesian analyses ####
library(brms)
library(bayesplot)
library(rstan)
library(matrixStats)

#stan doesn't like NA in data
trial_quest_data_ff=trial_quest_data_ff[!is.na(trial_quest_data_ff$Rchosen),]

#start with intercept-only model

#get default priors: a bit wide but should work
get_prior(data=trial_quest_data_ff,family=bernoulli,
  Rchosen~1+(1|ID))

#run intercept-only model
intercept_only=brm(data=trial_quest_data_ff,family=bernoulli,
  Rchosen~1+(1|ID),
  iter=3000,warmup=500,chains=3,cores=3,seed=123)

#diagnostics: 
check_hmc_diagnostics(intercept_only$fit)
traceplot(intercept_only$fit,pars=c('b_Intercept','sd_1'))

#plot model parameters: quick and dirty way
fit_io=extract(intercept_only$fit)
plot(density(fit_io$b_Intercept),xlim=c(-0.1,0.35))
abline(v=0)
plot(density(fit_io$sd_1))
plot(seq(1,190,by=1),colMedians(fit_io$r_1_1),pch=19,col='blue',
     ylab='Intercept Value',xlab='Participant Number')
abline(h=median(fit_io$b_Intercept),lwd=3)

#get raw stan code
io_code=make_stancode(data=trial_quest_data_ff,family=bernoulli,
  Rchosen~1+(1|ID),save_model='intercept_only.stan')

#another way of running this - helpful for understanding steps/debugging
# io_code=make_stancode(data=trial_quest_data_ff,family=bernoulli,
#   Rchosen~1+(1|ID))
# io_data=make_standata(data=trial_quest_data_ff,family=bernoulli,
#   Rchosen~1+(1|ID))
# io_model=stan_model(model_code=io_code)
# io_fit=sampling(io_model,data=io_data,
#   iter=3000,warmup=500,chains=3,cores=1,seed=123)

#prior sensitivity check - narrow prior for intercept
intercept_only_narrow=brm(data=trial_quest_data_ff,family=bernoulli,
  Rchosen~1+(1|ID),
  prior=c(prior(student_t(3,0,1),class=Intercept)),
  iter=3000,warmup=500,chains=3,cores=3,seed=123)
fit_io_narrow=extract(intercept_only_narrow$fit)

#compare output with different priors to see if it changes things
plot(density(fit_io_narrow$b_Intercept),xlim=c(-0.1,0.35),ylim=c(0,16))
abline(v=0)
par(new=TRUE)
plot(density(fit_io$b_Intercept),xlim=c(-0.1,0.35),col='red',ylim=c(0,16))

plot(density(fit_io_narrow$sd_1),xlim=c(0.15,0.5),ylim=c(0,16))
par(new=TRUE)
plot(density(fit_io$sd_1),col='red',xlim=c(0.15,0.5),ylim=c(0,16))

plot(seq(1,190,by=1),colMedians(fit_io_narrow$r_1_1),pch=19,col='blue',
     ylab='Intercept Value',xlab='Participant Number')
abline(h=median(fit_io_narrow$b_Intercept),lwd=3)
points(seq(1,190,by=1),colMedians(fit_io$r_1_1),pch=19,col=rgb(1,0,0,0.3))

#build up model - abbreviated version
#first, test whether allowing effect of uncertainty to vary by person improves
#fit
uncertainty_only=brm(data=trial_quest_data_ff,family=bernoulli,
  Rchosen~uncR+(1|ID),
  iter=3000,warmup=500,chains=3,cores=3,seed=123)
print(summary(uncertainty_only$fit,pars=c('b_Intercept','b'))$summary)

uncertainty_ind=brm(data=trial_quest_data_ff,family=bernoulli,
  Rchosen~uncR+(1+uncR|ID),
  iter=3000,warmup=500,chains=3,cores=3,seed=123)
print(summary(uncertainty_ind$fit,pars=c('b_Intercept','b'))$summary)

#example of comparing fit: loo_compare sorts so best fitting model is on top
uncertainty_only=add_criterion(uncertainty_only,"waic")
uncertainty_ind=add_criterion(uncertainty_ind,"waic")
loo_compare(uncertainty_only,uncertainty_ind,criterion="waic")

#better to start with simple model and build up in case deubugging is needed. but,
# for this example, skip to the full model
full_code=make_stancode(data=trial_quest_data_ff,family=bernoulli,
  Rchosen~1+horizon*valence*(uncR+scale(mean_diff))+
    (1+horizon*valence*(uncR+scale(mean_diff))|ID),
  save_model='full_model.stan')
full_model=brm(data=trial_quest_data_ff,family=bernoulli,
  Rchosen~1+horizon*valence*(uncR+scale(mean_diff))+
    (1+horizon*valence*(uncR+scale(mean_diff))|ID),
  iter=3000,warmup=500,chains=3,cores=3,seed=123)
check_hmc_diagnostics(full_model$fit)
print(summary(full_model$fit,pars=c('b_Intercept','b'))$summary)

#fancier plotting with bayesplot::mcmc_intervals
mcmc_intervals(full_model,pars=c('b_Intercept','b_horizon10','b_valenceloss',
'b_uncR','b_scalemean_diff','b_horizon10:valenceloss','b_horizon10:uncR',
'b_horizon10:scalemean_diff','b_valenceloss:uncR',
'b_valenceloss:scalemean_diff','b_horizon10:valenceloss:uncR',
'b_horizon10:valenceloss:scalemean_diff'),prob_outer = 0.95)
mcmc_intervals(full_model,pars=c('b_uncR','b_horizon10:uncR',
'b_valenceloss:uncR','b_horizon10:valenceloss:uncR'),prob_outer = 0.95)

#how strong of evidence? calculate proportion of samples in a direction
# total # samples = #chains*(total iterations-warmup iterations) = 3*(3000-500)=7500
fit_full=extract(full_model$fit,pars=c('b_Intercept','b_horizon10','b_valenceloss',
'b_uncR','b_scalemean_diff','b_horizon10:valenceloss','b_horizon10:uncR',
'b_horizon10:scalemean_diff','b_valenceloss:uncR',
'b_valenceloss:scalemean_diff','b_horizon10:valenceloss:uncR',
'b_horizon10:valenceloss:scalemean_diff'))
sum(fit_full$`b_horizon10:uncR`>0)/7500 #100%
sum(fit_full$b_uncR>0)/7500 #98%

#look at correlations among effects
fit_full_corr=extract(full_model$fit,pars='cor_1')$cor_1
hist(rowMedians(fit_full_corr),main='median correlation among participant-level predictors')
fit_full_popcorr=vcov(full_model,correlation=TRUE)
hist(fit_full_popcorr,main='median correlation among population-level predictors')

#predict reponses from fitted model & compare to empirical data plotted above
pp_full=posterior_predict(full_model)
pp_full_median=colMedians(pp_full)
trial_quest_data_ff_pred=trial_quest_data_ff
trial_quest_data_ff_pred$Rchosen=pp_full_median
trial_quest_data_ff_pred$predict=1
trial_quest_data_ff$predict=0
trial_quest_data_ff_comb=rbind(trial_quest_data_ff,trial_quest_data_ff_pred)
trial_quest_data_ff_comb$predict=as.factor(trial_quest_data_ff_comb$predict)
levels(trial_quest_data_ff_comb$predict)=c("empirical","predicted")
trial_quest_data_ff_avgpred=aggregate(Rchosen~horizon*valence*uncR*predict,
                                 data=trial_quest_data_ff_comb,FUN=mean)
ggplot(trial_quest_data_ff_avgpred,aes(y=Rchosen,x=as.factor(uncR),
  color=interaction(valence,horizon)))+facet_wrap(~predict)+
  geom_point(size=20,shape = '_')+
  labs(x='Relative Uncertainty of Right-sided Option',
       y='Proportion Right-sided Option Chosen')+
  theme_classic()+geom_hline(yintercept=0.5)

# extra: add L1 predictor
full_anx_code=make_stancode(data=trial_quest_data_ff,family=bernoulli,
  Rchosen~1+scale(dass_anx)*horizon*valence*(uncR+scale(mean_diff))+
    (1+horizon*valence*(uncR+scale(mean_diff))|ID),
  save_model='full_anx_model.stan')
full_anx_model=brm(data=trial_quest_data_ff,family=bernoulli,
  Rchosen~1+scale(dass_anx)*horizon*valence*(uncR+scale(mean_diff))+
    (1+horizon*valence*(uncR+scale(mean_diff))|ID),
  iter=3000,warmup=500,chains=3,cores=3,seed=123)
check_hmc_diagnostics(full_anx_model$fit)
print(summary(full_anx_model$fit,pars=c('b_Intercept','b'))$summary)
mcmc_intervals(full_anx_model,pars=c('b_uncR','b_horizon10:uncR',
'b_valenceloss:uncR','b_horizon10:valenceloss:uncR','b_scaledass_anx:uncR',
'b_scaledass_anx:horizon10:uncR','b_scaledass_anx:valenceloss:uncR',
'b_horizon10:valenceloss:uncR','b_scaledass_anx:horizon10:valenceloss:uncR'),
prob_outer = 0.95)

save.image('horizon_regressions_output.Rdata')


