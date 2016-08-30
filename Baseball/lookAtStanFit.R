# library(shinystan)
library(coda)
library(ggplot2)
library(ggmcmc)
library(CircStats)
library(reshape2)
library(plyr)
library(grid)

setwd("/Users/ghislaindentremont/Documents/TOJ/Baseball")
load("toj_color_post_June16th2016" )  # Rstan results
load("toj_trials.Rdata")  # actual toj data
load("color_trials.Rdata")  # actual color data
source("../EndogenousVisualPriorEntry-BayesianHierarchicalModel/functions.R")



############################################################################################
####                                        Diagnostics                                 ####
############################################################################################
# convert stanfit sample to dataframe table 
gg_toj_color_post = ggs(toj_color_post)

# look
gg_toj_color_post

# look at structure
str(gg_toj_color_post)

# THIS LIST IS ADAPTED FOR FUTURE ANALYSES
# list of parameters to examine
param_list = c("population_logit_rho_attention_effect_mean"
               , "population_logit_rho_intercept_mean"
               , "population_logit_rho_convention_effect_mean"
               , "population_logit_rho_attention_convention_interaction_effect_mean"
               , "population_log_kappa_attention_effect_mean"
               , "population_log_kappa_convention_effect_mean"
               , "population_log_kappa_attention_convention_interaction_effect_mean"
               , "population_log_kappa_intercept_mean"
               , "population_log_jnd_effect_mean"
               , "population_log_jnd_convention_effect_mean"
               , "population_log_jnd_intercept_mean"
               , "population_log_jnd_convention_interaction_effect_mean"
               , "population_pss_effect_mean"
               , "population_pss_convention_effect_mean" 
               , "population_pss_intercept_mean"
               , "population_pss_convention_interaction_effect_mean"
               , "zpopulation_logit_rho_effect_sd" 
               , "zpopulation_logit_rho_intercept_sd"
               , "zpopulation_log_kappa_effect_sd"
               , "zpopulation_log_kappa_intercept_sd"
               , "zpopulation_log_jnd_effect_sd"
               , "zpopulation_log_jnd_intercept_sd"
               , "zpopulation_pss_effect_sd"
              , "zpopulation_pss_intercept_sd"
              )
# THIS LIST IS ADAPTED FOR FUTURE ANALYSES

# look at posteriors
for (param in param_list) {
  ptm = proc.time()
  print( ggs_histogram(gg_toj_color_post, family = param) )
  print( proc.time() - ptm )
}

# posterior by chain
for (param in param_list) {
  ptm = proc.time()
  print( ggs_density(gg_toj_color_post, family = param) )
  print( proc.time() - ptm )
}

# traceplots
for (param in param_list) {
  ptm = proc.time()
  print( ggs_traceplot(gg_toj_color_post, family = param) )
  print( proc.time() - ptm )
}

# running means
for (param in param_list) {
  ptm = proc.time()
  print( ggs_running(gg_toj_color_post, family = param) )
  print( proc.time() - ptm )
}

# compare complete and last part of chains
for (param in param_list) {
  ptm = proc.time()
  print( ggs_compare_partial(gg_toj_color_post, family = param) )
  print( proc.time() - ptm )
}

# autocorrelation
for (param in param_list) {
  ptm = proc.time()
  print( ggs_autocorrelation(gg_toj_color_post, family = param) )
  print( proc.time() - ptm )
}
# NOTE: autocorrelation is not indicative of lack of convergence per se, but is indicative of misbehavior perhaps
# solution to autocorrelation is thinining



############################################################################################
####                         Posterior Predictive Checks                                ####
############################################################################################


#-------------------------------------- TOJ Actual Data -----------------------------------#
real_toj = aggregate(safe ~ soa2 + base_probe_dist + know_tie_goes_runner, data = toj_trials, FUN = mean)
real_toj[,2] = as.character(real_toj[,2])
real_toj[,3] = as.character(real_toj[,3])
#-------------------------------------- TOJ Actual Data -----------------------------------#


#-------------------------------------- TOJ Simulated Data --------------------------------#
### Get PSS Parameters
pss_intercept_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "population_pss_intercept_mean",]$value
pss_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "population_pss_effect_mean",]$value
pss_convention_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "population_pss_convention_effect_mean",]$value
pss_convention_interaction_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "population_pss_convention_interaction_effect_mean",]$value

pss_glove_know_mean_reps = get_condition_mean_sample(
  pss_intercept_mean - pss_convention_effect_mean/2
  , (pss_effect_mean - pss_convention_interaction_effect_mean )
  , FALSE
  , "null"
)

pss_glove_dontknow_mean_reps = get_condition_mean_sample(
  pss_intercept_mean + pss_convention_effect_mean/2
  , (pss_effect_mean + pss_convention_interaction_effect_mean )
  , FALSE
  , "null"
)

pss_base_know_mean_reps = get_condition_mean_sample(
  pss_intercept_mean - pss_convention_effect_mean/2
  , (pss_effect_mean - pss_convention_interaction_effect_mean )
  , TRUE
  , "null"
)

pss_base_dontknow_mean_reps = get_condition_mean_sample(
  pss_intercept_mean + pss_convention_effect_mean/2
  , (pss_effect_mean + pss_convention_interaction_effect_mean )
  , TRUE
  , "null"
)


### Get JND Parameters
jnd_intercept_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "population_logjnd_intercept_mean",]$value
jnd_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "population_logjnd_effect_mean",]$value
jnd_convention_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "population_logjnd_convention_effect_mean",]$value
jnd_convention_interaction_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "population_logjnd_convention_interaction_effect_mean",]$value

jnd_glove_know_mean_reps = get_condition_mean_sample(
  jnd_intercept_mean - jnd_convention_effect_mean/2
  , (jnd_effect_mean - jnd_convention_interaction_effect_mean )
  , FALSE
  , "log"
)

jnd_glove_dontknow_mean_reps = get_condition_mean_sample(
  jnd_intercept_mean + jnd_convention_effect_mean/2
  , (jnd_effect_mean + jnd_convention_interaction_effect_mean )
  , FALSE
  , "log"
)

jnd_base_know_mean_reps = get_condition_mean_sample(
  jnd_intercept_mean - jnd_convention_effect_mean/2
  , (jnd_effect_mean - jnd_convention_interaction_effect_mean )
  , TRUE
  , "log"
)

jnd_base_dontknow_mean_reps = get_condition_mean_sample(
  jnd_intercept_mean + jnd_convention_effect_mean/2
  , (jnd_effect_mean + jnd_convention_interaction_effect_mean )
  , TRUE
  , "log"
)
#-------------------------------------- TOJ Simulated Data --------------------------------#


#-------------------------------------- Do TOJ PPC ----------------------------------------#
SOAs = c(-250, -150, -100, -50, -17, 17, 50, 100, 150, 250)

do_toj_ppc(
  pss_glove_know_mean_reps
  , jnd_glove_know_mean_reps
  , "know convention and attend glove"
  , c("know_tie_goes_runner", "base_probe_dist")
  , c("TRUE", "0.2")
  , "safe proportion"
)

do_toj_ppc(
  pss_glove_dontknow_mean_reps
  , jnd_glove_dontknow_mean_reps
  , "don't know convention and attend glove"
  , c("know_tie_goes_runner", "base_probe_dist")
  , c("FALSE", "0.2")
  , "safe proportion"
)

do_toj_ppc(
  pss_base_know_mean_reps
  , jnd_base_know_mean_reps
  , "know convention and attend base"
  , c("know_tie_goes_runner", "base_probe_dist")
  , c("TRUE", "0.8")
  , "safe proportion"
)

do_toj_ppc(
  pss_base_dontknow_mean_reps
  , jnd_base_dontknow_mean_reps
  , "don't know convention and attend base"
  , c("know_tie_goes_runner", "base_probe_dist")
  , c("FALSE", "0.8")
  , "safe proportion"
)
#-------------------------------------- Do TOJ PPC ----------------------------------------#


#---------------------------------------- TOJ SDs -----------------------------------------#
id_all = unique(toj_trials$id)
get_pss_jnd = function(id_list) {
  toj_by_condition = ddply(
    .data = toj_trials[toj_trials$id %in% id_list,]
    , .variables = .(id, base_probe_dist, know_tie_goes_runner)
    , .fun = function(x){
      fit = glm(
        formula = safe~soa2
        , data = x
        # NOTE: The results change trivially depending on the link used
        , family = binomial(link = "probit")  # default is logit, but I used probit (normal cdf) for Bayesian analysis
      )
      to_return = data.frame(
        id = x$id[1]
        , pss = -coef(fit)[1]/coef(fit)[2]
        , slope = coef(fit)[2]
        , jnd = qnorm(0.84)/coef(fit)[2]  # mathematically the same as half the difference between .16 and .84 points 
      )
      return(to_return)
    }
  )
  return(toj_by_condition)
}
toj_by_condition = get_pss_jnd(id_all)

toj_pss_SD =  sd(aggregate(pss ~ id, data = toj_by_condition, FUN = mean)$pss)

get_violin( 
  tan(gg_toj_color_post[gg_toj_color_post$Parameter == "zpopulation_pss_intercept_sd",]$value) * 250
  , "PSS Intercept sd"
  , y_lab = "PSS (ms)"
  , hline = F
)

# JND
toj_jnd_SD =  sd(aggregate(jnd ~ id, data = toj_by_condition, FUN = mean)$jnd)

get_violin( 
  (
  exp(  # ROUGH ESTIMATE
    gg_toj_color_post[gg_toj_color_post$Parameter == "population_logjnd_intercept_mean",]$value 
     + tan(gg_toj_color_post[gg_toj_color_post$Parameter == "zpopulation_logjnd_intercept_sd",]$value) 
    ) - 
      exp(
        gg_toj_color_post[gg_toj_color_post$Parameter == "population_logjnd_intercept_mean",]$value 
      ) 
  ) * 250
  , "JND Intercept sd"
  , y_lab = "JND (ms)"
  , hline = F
)
#---------------------------------------- TOJ SDs -----------------------------------------#


#-------------------------------------- Color Actual Data ---------------------------------#
hist(color_trials[color_trials$attended == TRUE,]$color_diff_radians, breaks = 30, freq = F, col = rgb(.1,.1,.1,.5))
hist(color_trials[color_trials$attended == FALSE,]$color_diff_radians, breaks = 30, freq = F, col = rgb(.9,.9,.9,.5), add = T)
#-------------------------------------- Color Actual Data ---------------------------------#


#-------------------------------------- Color Simulated Data ------------------------------#
### Get Rho Parameters 
rho_intercept_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "population_logit_rho_intercept_mean",]$value
rho_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "population_logit_rho_attention_effect_mean",]$value
rho_convention_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "population_logit_rho_convention_effect_mean",]$value
rho_convention_interaction_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "population_logit_rho_attention_convention_interaction_effect_mean",]$value

rho_attend_know_reps = get_condition_mean_sample(
  rho_intercept_mean - rho_convention_effect_mean/2
  , (rho_effect_mean - rho_convention_interaction_effect_mean)
  , TRUE
  , "logit"
  )

rho_attend_dontknow_reps = get_condition_mean_sample(
  rho_intercept_mean + rho_convention_effect_mean/2
  , (rho_effect_mean + rho_convention_interaction_effect_mean)
  , TRUE
  , "logit"
)

rho_unattend_know_reps = get_condition_mean_sample(
  rho_intercept_mean - rho_convention_effect_mean/2
  , (rho_effect_mean - rho_convention_interaction_effect_mean)
  , FALSE
  , "logit"
)

rho_unattend_dontknow_reps = get_condition_mean_sample(
  rho_intercept_mean + rho_convention_effect_mean/2
  , (rho_effect_mean + rho_convention_interaction_effect_mean)
  , FALSE
  , "logit"
)

### Get Kappa Parameters
kappa_intercept_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "population_log_kappa_intercept_mean",]$value
kappa_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "population_log_kappa_attention_effect_mean",]$value
kappa_convention_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "population_log_kappa_convention_effect_mean",]$value
kappa_convention_interaction_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "population_log_kappa_attention_convention_interaction_effect_mean",]$value

kappa_attend_know_reps = get_condition_mean_sample(
  kappa_intercept_mean - kappa_convention_effect_mean/2
  , (kappa_effect_mean - kappa_convention_interaction_effect_mean)
  , TRUE
  , "log_free"
)

kappa_attend_dontknow_reps = get_condition_mean_sample(
  kappa_intercept_mean + kappa_convention_effect_mean/2
  , (kappa_effect_mean + kappa_convention_interaction_effect_mean)
  , TRUE
  , "log_free"
)

kappa_unattend_know_reps = get_condition_mean_sample(
  kappa_intercept_mean - kappa_convention_effect_mean/2
  , (kappa_effect_mean - kappa_convention_interaction_effect_mean)
  , FALSE
  , "log_free"
)

kappa_unattend_dontknow_reps = get_condition_mean_sample(
  kappa_intercept_mean + kappa_convention_effect_mean/2
  , (kappa_effect_mean + kappa_convention_interaction_effect_mean)
  , FALSE
  , "log_free"
)
#-------------------------------------- Color Simulated Data ------------------------------#


#-------------------------------------- Do Color PPC --------------------------------------#
do_color_ppc(
  rho_attend_know_reps
  , kappa_attend_know_reps
  , "know convention and attend"
  , c("know_tie_goes_runner", "attended")
  , c("TRUE", "TRUE")
)

do_color_ppc(
  rho_attend_dontknow_reps
  , kappa_attend_dontknow_reps
  , "don't know convention and attend"
  , c("know_tie_goes_runner", "attended")
  , c("FALSE", "TRUE")
)

do_color_ppc(
  rho_unattend_know_reps
  , kappa_unattend_know_reps
  , "know convention and unattend"
  , c("know_tie_goes_runner", "attended")
  , c("TRUE", "FALSE")
)

do_color_ppc(
  rho_unattend_dontknow_reps
  , kappa_unattend_dontknow_reps
  , "don't know convention and unattend"
  , c("know_tie_goes_runner", "attended")
  , c("FALSE", "FALSE")
)
#-------------------------------------- Do Color PPC --------------------------------------#



############################################################################################
####                                       Analysis                                     ####
############################################################################################

#------------------------------------------------------------------------------------------#
#--------------------------------- Correlations -------------------------------------------#
#------------------------------------------------------------------------------------------#
## ORDER:
# (1) population_pss_intercept_mean      
# (2) population_pss_effect_mean          
# (3) population_log_jnd_intercept_mean    
# (4) population_log_jnd_effect_mean     
# (5) population_logit_rho_intercept_mean                         
# (6) population_log_kappa_intercept_mean                        
# (7) population_logit_rho_attention_effect_mean                 
# (8) population_log_kappa_attention_effect_mean                   

# for quick look
ggs_caterpillar(gg_toj_color_post, family = "cor", thick_ci = c(0.25, 0.75) ) + geom_vline(xintercept = 0, col = "red")


#-------------------------------------- Get Betas -----------------------------------------#
# extract samples
detach('package:rstan', unload = T)  # to ensure 
library(rstan)

# FOR ANALYSIS USED FOR SUBMISSION
# version so nothing gets flipped twice by accident 
ex_toj_color_post2 = extract(toj_color_post)

# flip pss and jnd effects to make positive values indicate predicted results
ex_toj_color_post = ex_toj_color_post2
ex_toj_color_post$population_pss_effect_mean = -ex_toj_color_post2$population_pss_effect_mean
ex_toj_color_post$population_pss_convention_interaction_effect_mean = -ex_toj_color_post2$population_pss_convention_interaction_effect_mean
ex_toj_color_post$population_logjnd_effect_mean = -ex_toj_color_post2$population_logjnd_effect_mean
ex_toj_color_post$population_logjnd_convention_interaction_effect_mean = -ex_toj_color_post2$population_logjnd_convention_interaction_effect_mean
# FOR ANALYSIS USED FOR SUBMISSION

# get correlation posteriors 
pos_corr2 = data.frame(value = ex_toj_color_post$cor)
pos_corr2$id = rownames(pos_corr2)
pos_corr = melt( pos_corr2 )
names(pos_corr)[2] = c("parameter")
# WARNING: CORRELATIONS CONTAINING 2s OR (EXCLUSIVE) 4s WILL BE IN THE WRONG DIRECTION

# (1) population_pss_intercept_mean      
# (2) population_pss_effect_mean          
# (3) population_log_jnd_intercept_mean    
# (4) population_log_jnd_effect_mean     
# (5) population_logit_rho_intercept_mean                         
# (6) population_log_kappa_intercept_mean                        
# (7) population_logit_rho_attention_effect_mean                 
# (8) population_log_kappa_attention_effect_mean 
betas2 = data.frame(value = ex_toj_color_post$beta)
betas2$iteration = rownames(betas2)

# melt
ptm = proc.time()
betas3 = melt( betas2 )  # this takes a while
proc.time() - ptm

betas3$parameter = rep( c(
  "population_pss_intercept_mean"      
  , "population_pss_effect_mean"          
  , "population_logjnd_intercept_mean"    
  , "population_logjnd_effect_mean"     
  , "population_logit_rho_intercept_mean"                       
  , "population_log_kappa_intercept_mean"                       
  , "population_logit_rho_attention_effect_mean"                 
  , "population_log_kappa_attention_effect_mean"
)
, times = 1
, each = nrow(betas2)*length(unique(betas3$variable))/8  # 8 is number of parameters 
)  
betas3$participant = rep(c(1:length(unique(toj_trials$id))), times = 8, each = nrow(betas2))

# flip betas 
betas = betas3
betas[betas$parameter=="population_pss_effect_mean",]$value = -betas3[betas3$parameter=="population_pss_effect_mean" ,]$value
betas[betas$parameter=="population_logjnd_effect_mean",]$value = -betas3[betas3$parameter=="population_logjnd_effect_mean",]$value
#-------------------------------------- Get Betas -----------------------------------------#


#--------------------------- PSS Intercept Shrinkage --------------------------------------#
pssinterceptmean = extract_samples("population_pss_intercept_mean")

pssinterceptsd = extract_samples("zpopulation_pss_intercept_sd", TRUE)

pssconventioneffect = extract_samples("population_pss_convention_effect_mean")

conventionfactor = ifelse(aggregate(know_tie_goes_runner~id,data = toj_trials, FUN =unique)$know_tie_goes_runner, -1, 1)

pssintercept_ids = ddply(
  .data = betas
  , .variables = .(participant)
  , .fun = function(x){
    i = unique(x$participant)
    x_use = x[x$parameter ==  "population_pss_intercept_mean",]$value
    pssintercept = median(
      pssinterceptmean + pssinterceptsd*x_use 
    + pssconventioneffect*conventionfactor[i]/2
    )
    df = data.frame(pssintercept*250, conventionfactor[i])
    names(df) = c("pssintercept", "conventionfactor")
    return(df)
  }
)

realpssintercept_ids = aggregate(pss ~ id, data = toj_by_condition, FUN = mean)
  
plot(as.numeric(realpssintercept_ids$id), realpssintercept_ids$pss, col = 'red')
abline(h = mean(realpssintercept_ids$pss), col = 'red')
points(pssintercept_ids$pssintercept, col = 'blue')
abline(h = mean(pssintercept_ids$pssintercept), col = 'blue')

# compare sd values computed in this way
sd(realpssintercept_ids$pss)
sd(pssintercept_ids$pssintercept)
#--------------------------- PSS Intercept Shrinkage --------------------------------------#


#--------------------------- JND Intercept Shrinkage --------------------------------------#
jndinterceptmean = extract_samples("population_logjnd_intercept_mean")

jndinterceptsd = extract_samples("zpopulation_logjnd_intercept_sd", TRUE)

jndconventioneffect = extract_samples("population_logjnd_convention_effect_mean")

conventionfactor = ifelse(aggregate(know_tie_goes_runner~id,data = toj_trials, FUN =unique)$know_tie_goes_runner, -1, 1)

jndintercept_ids = ddply(
  .data = betas
  , .variables = .(participant)
  , .fun = function(x){
    i = unique(x$participant)
    x_use = x[x$parameter ==  "population_logjnd_intercept_mean",]$value
    jndintercept = median(
      exp(
        jndinterceptmean + jndinterceptsd*x_use 
        + jndconventioneffect*conventionfactor[i]/2
      )
    ) 
    df = data.frame(jndintercept*250, conventionfactor[i])
    names(df) = c("jndintercept", "conventionfactor")
    return(df)
  }
)

realjndintercept_ids = aggregate(jnd ~ id, data = toj_by_condition, FUN = mean)

plot(as.numeric(realjndintercept_ids$id), realjndintercept_ids$jnd, col = 'red')
abline(h = mean(realjndintercept_ids$jnd), col = 'red')
points(jndintercept_ids$jndintercept, col = 'blue')
abline(h = mean(jndintercept_ids$jndintercept), col = 'blue')

# compare sd values computed in this way
sd(realjndintercept_ids$jnd)
sd(jndintercept_ids$jndintercept)
#--------------------------- JND Intercept Shrinkage --------------------------------------#



#------------------------------------------------------------------------------------------#
#--------------------------------- PSS Scale ----------------------------------------------#
#------------------------------------------------------------------------------------------#

#---------------------------------- PSS Conditions-----------------------------------------#

#---------------------------------- PSS Conditions-----------------------------------------#


#---------------------------------- PSS Intercept -----------------------------------------#
get_violin(
  ex_toj_color_post$population_pss_intercept_mean * 250  
  , "PSS Intercept Mean" 
  , y_lab = "SOA (ms)"
)
#---------------------------------- PSS Intercept -----------------------------------------#


#---------------------------------- PSS Attention Effects ---------------------------------#
# effect of attention on PSS and JND
get_violin(
  ex_toj_color_post$population_pss_effect_mean * 250
  , "PSS Attention Effect Mean"
  , y_lab = "SOA (Glove - Base; ms)"
)
#---------------------------------- PSS Attention Effects ---------------------------------#


#--------------------------------- PSS Convention Effects ---------------------------------#
# effect of convention knowlege on PSS 
get_violin(
  ex_toj_color_post$population_pss_convention_effect_mean * 250
  , "PSS Convention\nEffect Mean" 
  , y_lab = "SOA (Don't Know - Know; ms)"
)
# PSS: DK 
get_95_HDI(
  (ex_toj_color_post$population_pss_intercept_mean + ex_toj_color_post$population_pss_convention_effect_mean/2) * 250
)
# PSS: K
get_95_HDI(
( ex_toj_color_post$population_pss_intercept_mean - ex_toj_color_post$population_pss_convention_effect_mean/2)  * 250
)

# effect of interaction between convention knowledge and attention on PSS
get_violin(
  # multiply by two
  2 * ex_toj_color_post$population_pss_convention_interaction_effect_mean * 250 
  , "PSS Convention Knowledge\nInteraction Effect Mean"
  , y_lab = "SOA (ms)"
)

# effect of attention on PSS by knowledge of convention
get_violin(
  c(
  ( (ex_toj_color_post$population_pss_intercept_mean + ex_toj_color_post$population_pss_convention_effect_mean/2 + (ex_toj_color_post$population_pss_effect_mean + ex_toj_color_post$population_pss_convention_interaction_effect_mean)/2) 
    - (ex_toj_color_post$population_pss_intercept_mean + ex_toj_color_post$population_pss_convention_effect_mean/2 - (ex_toj_color_post$population_pss_effect_mean + ex_toj_color_post$population_pss_convention_interaction_effect_mean)/2 ) ) * 250
  , ( (ex_toj_color_post$population_pss_intercept_mean -  ex_toj_color_post$population_pss_convention_effect_mean/2 + (ex_toj_color_post$population_pss_effect_mean - ex_toj_color_post$population_pss_convention_interaction_effect_mean)/2) 
      - (ex_toj_color_post$population_pss_intercept_mean -  ex_toj_color_post$population_pss_convention_effect_mean/2 - (ex_toj_color_post$population_pss_effect_mean - ex_toj_color_post$population_pss_convention_interaction_effect_mean)/2 ) ) * 250
  )
  , c("PSS Attention Effect\nGiven Don't Know", "PSS Attention Effect\nGiven Know")
  , y_lab = "SOA (Glove - Base; ms)"
)
# DK + A
get_95_HDI(
  (ex_toj_color_post$population_pss_intercept_mean + ex_toj_color_post$population_pss_convention_effect_mean/2 + (ex_toj_color_post$population_pss_effect_mean + ex_toj_color_post$population_pss_convention_interaction_effect_mean)/2) * 250 
)
# DK + U
get_95_HDI(
  (ex_toj_color_post$population_pss_intercept_mean + ex_toj_color_post$population_pss_convention_effect_mean/2 - (ex_toj_color_post$population_pss_effect_mean + ex_toj_color_post$population_pss_convention_interaction_effect_mean)/2  ) * 250
)
# K + A
get_95_HDI(
  (ex_toj_color_post$population_pss_intercept_mean -  ex_toj_color_post$population_pss_convention_effect_mean/2 + (ex_toj_color_post$population_pss_effect_mean - ex_toj_color_post$population_pss_convention_interaction_effect_mean)/2) * 250 
)
# K + U
get_95_HDI(
  (ex_toj_color_post$population_pss_intercept_mean -  ex_toj_color_post$population_pss_convention_effect_mean/2 - (ex_toj_color_post$population_pss_effect_mean - ex_toj_color_post$population_pss_convention_interaction_effect_mean)/2 )  * 250
)
#--------------------------------- PSS Convention Effects ---------------------------------#



#------------------------------------------------------------------------------------------#
#--------------------------------- JND Scale ----------------------------------------------#
#------------------------------------------------------------------------------------------#

#---------------------------------- JND Conditions-----------------------------------------#
# DK + A
jnd_pos_attn_pos_conv =
  exp(ex_toj_color_post$population_logjnd_intercept_mean + ex_toj_color_post$population_logjnd_convention_effect_mean/2 
      + (ex_toj_color_post$population_logjnd_effect_mean + ex_toj_color_post$population_logjnd_convention_interaction_effect_mean)/2)
# DK + U
jnd_neg_attn_pos_conv =
  exp(ex_toj_color_post$population_logjnd_intercept_mean + ex_toj_color_post$population_logjnd_convention_effect_mean/2 
      - (ex_toj_color_post$population_logjnd_effect_mean + ex_toj_color_post$population_logjnd_convention_interaction_effect_mean)/2)  
# K + A
jnd_pos_attn_neg_conv =
  exp(ex_toj_color_post$population_logjnd_intercept_mean -  ex_toj_color_post$population_logjnd_convention_effect_mean/2 
      + (ex_toj_color_post$population_logjnd_effect_mean - ex_toj_color_post$population_logjnd_convention_interaction_effect_mean)/2)
# K + U
jnd_neg_attn_neg_conv =
  exp(ex_toj_color_post$population_logjnd_intercept_mean -  ex_toj_color_post$population_logjnd_convention_effect_mean/2 
      - (ex_toj_color_post$population_logjnd_effect_mean - ex_toj_color_post$population_logjnd_convention_interaction_effect_mean)/2)
#---------------------------------- JND Conditions-----------------------------------------#


#---------------------------------- JND Intercept -----------------------------------------#
# BAD:
get_violin(
  exp( ex_toj_color_post$population_logjnd_intercept_mean ) * 250 
  , "JND Intercept Mean"
  , y_lab = "SOA (ms)"
  , hline = FALSE
)

# GOOD:
get_violin(
  ( jnd_neg_attn_neg_conv
  + jnd_pos_attn_neg_conv
  + jnd_neg_attn_pos_conv
  + jnd_pos_attn_pos_conv
  )/4 * 250 
  , "JND Intercept Mean"
  , y_lab = "SOA (ms)"
  , hline = FALSE
)

# BETTER:
get_violin(
  ex_toj_color_post$population_logjnd_intercept_mean
  , "JND Intercept Mean"
  , y_lab = "Log of SOA (ms)"
  , hline = FALSE
)
#---------------------------------- JND Intercept -----------------------------------------#


#---------------------------------- JND Attention Effects ---------------------------------#
# BAD:
get_violin(
  ( exp( ex_toj_color_post$population_logjnd_intercept_mean + ex_toj_color_post$population_logjnd_effect_mean/2 )
    - exp( ex_toj_color_post$population_logjnd_intercept_mean - ex_toj_color_post$population_logjnd_effect_mean/2  ) ) * 250 
  , "JND Attention Effect Mean"
  , y_lab = "SOA (Glove - Base; ms)"
)

# GOOD:
get_violin(
  ( jnd_pos_attn_neg_conv - jnd_neg_attn_neg_conv
  + jnd_pos_attn_pos_conv - jnd_neg_attn_pos_conv
  )/2 * 250 
  , "JND Attention Effect Mean"
  , y_lab = "SOA (Glove - Base; ms)"
)

# BETTER:
get_violin(
  ex_toj_color_post$population_logjnd_effect_mean
  , "JND Attention Effect Mean"
  , y_lab = "Log of SOA (Glove - Base; ms)"
)
#---------------------------------- JND Attention Effects ---------------------------------#


#--------------------------------- JND Convention Effects ---------------------------------#
### BAD:
get_violin(
  ( exp( ex_toj_color_post$population_logjnd_intercept_mean + ex_toj_color_post$population_logjnd_convention_effect_mean/2 )
        - exp( ex_toj_color_post$population_logjnd_intercept_mean - ex_toj_color_post$population_logjnd_convention_effect_mean/2  ) ) * 250 
  ,  "JND Convention\nEffect Mean"
  , y_lab = "SOA (Don't Know - Know; ms)"
)
# DK 
get_95_HDI(
  exp(ex_toj_color_post$population_log_jnd_intercept_mean + ex_toj_color_post$population_log_jnd_convention_effect_mean/2) * 250
)

# K
get_95_HDI(
  exp( ex_toj_color_post$population_log_jnd_intercept_mean - ex_toj_color_post$population_log_jnd_convention_effect_mean/2)  * 250
)

# GOOD:
get_violin(
  ( jnd_neg_attn_pos_conv - jnd_neg_attn_neg_conv
  + jnd_pos_attn_pos_conv - jnd_pos_attn_neg_conv
  )/2 * 250 
  ,  "JND Convention\nEffect Mean"
  , y_lab = "SOA (Don't Know - Know; ms)"
)
# DK 
get_95_HDI(
  (jnd_neg_attn_pos_conv + jnd_pos_attn_pos_conv)/2 * 250
  )
# K
get_95_HDI(
  (jnd_neg_attn_neg_conv + jnd_pos_attn_neg_conv)/2 * 250
)

# BETTER:
get_violin(
  ex_toj_color_post$population_logjnd_convention_effect_mean
  ,  "JND Convention\nEffect Mean"
  , y_lab = "Log of SOA (Don't Know - Know; ms)"
)
# DK 
get_95_HDI(
  ex_toj_color_post$population_logjnd_intercept_mean + ex_toj_color_post$population_logjnd_convention_effect_mean/2
)
# K
get_95_HDI(
  ex_toj_color_post$population_logjnd_intercept_mean - ex_toj_color_post$population_logjnd_convention_effect_mean/2
)


### BAD: effect of interaction between convention knowledge and attention on JND
get_violin(
  ( exp(ex_toj_color_post$population_log_jnd_intercept_mean + ex_toj_color_post$population_log_jnd_convention_interaction_effect_mean/2) 
    - exp(ex_toj_color_post$population_log_jnd_intercept_mean - ex_toj_color_post$population_log_jnd_convention_interaction_effect_mean/2) ) * 250
  , "JND Convention Knowledge\nInteraction Effect Mean"
  , y_lab = "SOA (ms)"
)

# GOOD:
get_violin(
  (
    (jnd_pos_attn_pos_conv - jnd_neg_attn_pos_conv) 
    - (jnd_pos_attn_neg_conv - jnd_neg_attn_neg_conv ) 
  ) * 250 
  , "JND Convention Knowledge\nInteraction Effect Mean"
  , y_lab = "SOA (ms)"
)

# BETTER: effect of interaction between convention knowledge and attention on JND
get_violin(
  2 * ex_toj_color_post$population_logjnd_convention_interaction_effect_mean
  , "JND Convention Knowledge\nInteraction Effect Mean"
  , y_lab = "Log of SOA (ms)"
)


### GOOD:  effect of attention on JND by knowledge of convention
get_violin(
  c(
    exp(ex_toj_color_post$population_logjnd_intercept_mean + ex_toj_color_post$population_logjnd_convention_effect_mean/2 + (ex_toj_color_post$population_logjnd_effect_mean + ex_toj_color_post$population_logjnd_convention_interaction_effect_mean)/2) *250
      - exp(ex_toj_color_post$population_logjnd_intercept_mean + ex_toj_color_post$population_logjnd_convention_effect_mean/2 - (ex_toj_color_post$population_logjnd_effect_mean + ex_toj_color_post$population_logjnd_convention_interaction_effect_mean)/2 ) * 250
    , exp(ex_toj_color_post$population_logjnd_intercept_mean -  ex_toj_color_post$population_logjnd_convention_effect_mean/2 + (ex_toj_color_post$population_logjnd_effect_mean - ex_toj_color_post$population_logjnd_convention_interaction_effect_mean)/2) * 250
        - exp(ex_toj_color_post$population_logjnd_intercept_mean -  ex_toj_color_post$population_logjnd_convention_effect_mean/2 - (ex_toj_color_post$population_logjnd_effect_mean - ex_toj_color_post$population_logjnd_convention_interaction_effect_mean)/2 )  * 250
  )
  , c("JND Attention Effect\nGiven Don't Know"  , "JND Attention Effect\nGiven Know")
  , y_lab = "SOA (Glove - Base; ms)"
)
get_95_HDI(
  exp(ex_toj_color_post$population_logjnd_intercept_mean -  ex_toj_color_post$population_logjnd_convention_effect_mean/2 + (ex_toj_color_post$population_logjnd_effect_mean - ex_toj_color_post$population_logjnd_convention_interaction_effect_mean)/2) * 250
  - exp(ex_toj_color_post$population_logjnd_intercept_mean -  ex_toj_color_post$population_logjnd_convention_effect_mean/2 - (ex_toj_color_post$population_logjnd_effect_mean - ex_toj_color_post$population_logjnd_convention_interaction_effect_mean)/2 )  * 250
)

# DK + A
get_95_HDI(
  exp(ex_toj_color_post$population_log_jnd_intercept_mean + ex_toj_color_post$population_log_jnd_convention_effect_mean/2 + (ex_toj_color_post$population_log_jnd_effect_mean + ex_toj_color_post$population_log_jnd_convention_interaction_effect_mean/2)/2) * 250
)
# DK + U
get_95_HDI(
  exp(ex_toj_color_post$population_log_jnd_intercept_mean + ex_toj_color_post$population_log_jnd_convention_effect_mean/2 - (ex_toj_color_post$population_log_jnd_effect_mean + ex_toj_color_post$population_log_jnd_convention_interaction_effect_mean/2)/2  ) * 250
)
# K + A
get_95_HDI(
  exp(ex_toj_color_post$population_log_jnd_intercept_mean -  ex_toj_color_post$population_log_jnd_convention_effect_mean/2 + (ex_toj_color_post$population_log_jnd_effect_mean - ex_toj_color_post$population_log_jnd_convention_interaction_effect_mean/2)/2) * 250 
)
# K + U
get_95_HDI(
  exp(ex_toj_color_post$population_log_jnd_intercept_mean -  ex_toj_color_post$population_log_jnd_convention_effect_mean/2 - (ex_toj_color_post$population_log_jnd_effect_mean - ex_toj_color_post$population_log_jnd_convention_interaction_effect_mean/2)/2 )  * 250
)

# BETTER:  effect of attention on JND by knowledge of convention
get_violin(
  c(
    (ex_toj_color_post$population_logjnd_intercept_mean + ex_toj_color_post$population_logjnd_convention_effect_mean/2 + (ex_toj_color_post$population_logjnd_effect_mean + ex_toj_color_post$population_logjnd_convention_interaction_effect_mean)/2) 
    - (ex_toj_color_post$population_logjnd_intercept_mean + ex_toj_color_post$population_logjnd_convention_effect_mean/2 - (ex_toj_color_post$population_logjnd_effect_mean + ex_toj_color_post$population_logjnd_convention_interaction_effect_mean)/2 ) 
    , (ex_toj_color_post$population_logjnd_intercept_mean -  ex_toj_color_post$population_logjnd_convention_effect_mean/2 + (ex_toj_color_post$population_logjnd_effect_mean - ex_toj_color_post$population_logjnd_convention_interaction_effect_mean)/2) 
    - (ex_toj_color_post$population_logjnd_intercept_mean -  ex_toj_color_post$population_logjnd_convention_effect_mean/2 - (ex_toj_color_post$population_logjnd_effect_mean - ex_toj_color_post$population_logjnd_convention_interaction_effect_mean)/2 )
  )
  , c("JND Attention Effect\nGiven Don't Know"  , "JND Attention Effect\nGiven Know")
  , y_lab = "Log of SOA (Glove - Base; ms)"
)
# DK + A
get_95_HDI(
  (ex_toj_color_post$population_logjnd_intercept_mean + ex_toj_color_post$population_logjnd_convention_effect_mean/2 + (ex_toj_color_post$population_logjnd_effect_mean + ex_toj_color_post$population_logjnd_convention_interaction_effect_mean)/2) 
)
# DK + U
get_95_HDI(
  (ex_toj_color_post$population_logjnd_intercept_mean + ex_toj_color_post$population_logjnd_convention_effect_mean/2 - (ex_toj_color_post$population_logjnd_effect_mean + ex_toj_color_post$population_logjnd_convention_interaction_effect_mean)/2  ) 
)
# K + A
get_95_HDI(
  (ex_toj_color_post$population_logjnd_intercept_mean -  ex_toj_color_post$population_logjnd_convention_effect_mean/2 + (ex_toj_color_post$population_logjnd_effect_mean - ex_toj_color_post$population_logjnd_convention_interaction_effect_mean)/2) 
)
# K + U
get_95_HDI(
  (ex_toj_color_post$population_logjnd_intercept_mean -  ex_toj_color_post$population_logjnd_convention_effect_mean/2 - (ex_toj_color_post$population_logjnd_effect_mean - ex_toj_color_post$population_logjnd_convention_interaction_effect_mean)/2 )
)

#--------------------------------- SOA Convention Effects ---------------------------------#



#------------------------------------------------------------------------------------------#
#--------------------------------- Rho Scale ----------------------------------------------#
#------------------------------------------------------------------------------------------#

#---------------------------------- Rho Conditions-----------------------------------------#
# DK + A
rho_pos_attn_pos_conv =
  plogis(ex_toj_color_post$logitRhoMean + ex_toj_color_post$logitRhoConventionEffectMean/2 
      + (ex_toj_color_post$logitRhoEffectMean + ex_toj_color_post$logitRhoConventionInteractionEffectMean)/2)
# DK + U
rho_neg_attn_pos_conv =
  plogis(ex_toj_color_post$logitRhoMean + ex_toj_color_post$logitRhoConventionEffectMean/2 
      - (ex_toj_color_post$logitRhoEffectMean + ex_toj_color_post$logitRhoConventionInteractionEffectMean)/2)  
# K + A
rho_pos_attn_neg_conv =
  plogis(ex_toj_color_post$logitRhoMean -  ex_toj_color_post$logitRhoConventionEffectMean/2 
      + (ex_toj_color_post$logitRhoEffectMean - ex_toj_color_post$logitRhoConventionInteractionEffectMean)/2)
# K + U
rho_neg_attn_neg_conv =
  plogis(ex_toj_color_post$logitRhoMean -  ex_toj_color_post$logitRhoConventionEffectMean/2 
      - (ex_toj_color_post$logitRhoEffectMean - ex_toj_color_post$logitRhoConventionInteractionEffectMean)/2)
#---------------------------------- Rho Conditions-----------------------------------------#



#---------------------------------- Rho Intercept -----------------------------------------#
# BAD:
get_violin(
  plogis(ex_toj_color_post$population_logit_rho_intercept_mean)
  , "Probability of Encoding Intercept Mean"
  , y_lab = "\u03C1"
  , hline = FALSE
)

# GOOD:
get_violin(
  ( rho_neg_attn_neg_conv
  + rho_neg_attn_pos_conv
  + rho_pos_attn_neg_conv
  + rho_pos_attn_pos_conv
  )/4
  , "Probability of Encoding Intercept Mean"
  , y_lab = "\u03C1"
  , hline = FALSE
)

# BETTER: 
get_violin(
  ex_toj_color_post$logitRhoMean
  , "Probability of Encoding Intercept Mean"
  , y_lab = "Log-odds of \u03C1"
  , hline = FALSE
)
#---------------------------------- Rho Intercept -----------------------------------------#


#-------------------------------- Rho Attention Effect ------------------------------------#
# BAD:
get_violin(
  ( plogis(ex_toj_color_post$population_logit_rho_intercept_mean + ex_toj_color_post$population_logit_rho_attention_effect_mean/2 )
    - plogis(ex_toj_color_post$population_logit_rho_intercept_mean - ex_toj_color_post$population_logit_rho_attention_effect_mean/2 ) )
  , "Probability of Encoding Effect Mean"
  , y_lab = "\u03C1 (Attended - Unattended)"
)

# GOOD:
get_violin(
  ( rho_pos_attn_neg_conv - rho_neg_attn_neg_conv
  + rho_pos_attn_pos_conv - rho_neg_attn_pos_conv
  )/2
  , "Probability of Encoding Effect Mean"
  , y_lab = "\u03C1 (Attended - Unattended)"
)

# BETTER:
get_violin(
  ex_toj_color_post$logitRhoEffectMean
  , "Probability of Encoding Effect Mean"
  , y_lab = "Log-odds of \u03C1 (Attended - Unattended)"
)
#-------------------------------- Rho Attention Effect ------------------------------------#


#-------------------------------- Rho Convention Effects ----------------------------------#
### BAD:
get_violin(
  ( plogis(ex_toj_color_post$population_logit_rho_intercept_mean + ex_toj_color_post$population_logit_rho_convention_effect_mean/2 )
    - plogis(ex_toj_color_post$population_logit_rho_intercept_mean - ex_toj_color_post$population_logit_rho_convention_effect_mean/2 ) )
  , "Probability of Encoding\nConvention Effect Mean"
  , y_lab = "\u03C1 (Don't Know - Know)"
)
# DK
get_95_HDI(
  plogis(ex_toj_color_post$population_logit_rho_intercept_mean + ex_toj_color_post$population_logit_rho_convention_effect_mean/2)
  )
# K
get_95_HDI(
  plogis(ex_toj_color_post$population_logit_rho_intercept_mean - ex_toj_color_post$population_logit_rho_convention_effect_mean/2 )
)

# GOOD:
get_violin(
  ( rho_pos_attn_pos_conv - rho_pos_attn_neg_conv 
  + rho_neg_attn_pos_conv - rho_neg_attn_neg_conv
  )/2
  , "Probability of Encoding\nConvention Effect Mean"
  , y_lab = "\u03C1 (Don't Know - Know)"
)
# DK
get_95_HDI(
  (rho_neg_attn_pos_conv + rho_pos_attn_pos_conv)/2
)
# K
get_95_HDI(
  (rho_neg_attn_neg_conv + rho_pos_attn_neg_conv)/2
)

# BETTER:
get_violin(
  ex_toj_color_post$logitRhoConventionEffectMean
  , "Probability of Encoding\nConvention Effect Mean"
  , y_lab = "Log-odds of \u03C1 (Don't Know - Know)"
)
# DK
get_95_HDI(
  (ex_toj_color_post$logitRhoMean + ex_toj_color_post$logitRhoConventionEffectMean/2)
)
# K
get_95_HDI(
  (ex_toj_color_post$logitRhoMean - ex_toj_color_post$logitRhoConventionEffectMean/2 )
)


### BAD: interaction
get_violin(
  ( plogis(ex_toj_color_post$population_logit_rho_intercept_mean + ex_toj_color_post$population_logit_rho_attention_convention_interaction_effect_mean/2 )
    - plogis(ex_toj_color_post$population_logit_rho_intercept_mean - ex_toj_color_post$population_logit_rho_attention_convention_interaction_effect_mean/2 ) )
  , "Probability of Encoding\nConvention Knowledge\nInteraction Effect Mean"
  , y_lab = "\u03C1"
)

# GOOD: interaction
get_violin(
  ( rho_pos_attn_pos_conv - rho_neg_attn_pos_conv )
  - ( rho_pos_attn_neg_conv - rho_neg_attn_neg_conv )
  , "Probability of Encoding\nConvention Knowledge\nInteraction Effect Mean"
  , y_lab = "\u03C1"
)

# BETTER: interaction
get_violin(
  2 * ex_toj_color_post$logitRhoConventionInteractionEffectMean
  , "Probability of Encoding\nConvention Knowledge\nInteraction Effect Mean"
  , y_lab = "Log-odds of \u03C1"
)

### GOOD: effect of attention on rho for each convention knowledge level
get_violin(
  c(
  ( plogis(ex_toj_color_post$population_logit_rho_intercept_mean + ex_toj_color_post$population_logit_rho_convention_effect_mean/2
           + (ex_toj_color_post$population_logit_rho_attention_effect_mean + ex_toj_color_post$population_logit_rho_attention_convention_interaction_effect_mean/2)/2 )
    - plogis(ex_toj_color_post$population_logit_rho_intercept_mean + ex_toj_color_post$population_logit_rho_convention_effect_mean/2
             - (ex_toj_color_post$population_logit_rho_attention_effect_mean +  ex_toj_color_post$population_logit_rho_attention_convention_interaction_effect_mean/2)/2 ) )
  , ( plogis(ex_toj_color_post$population_logit_rho_intercept_mean - ex_toj_color_post$population_logit_rho_convention_effect_mean/2
             + (ex_toj_color_post$population_logit_rho_attention_effect_mean - ex_toj_color_post$population_logit_rho_attention_convention_interaction_effect_mean/2)/2 )
      - plogis(ex_toj_color_post$population_logit_rho_intercept_mean - ex_toj_color_post$population_logit_rho_convention_effect_mean/2
               - (ex_toj_color_post$population_logit_rho_attention_effect_mean -  ex_toj_color_post$population_logit_rho_attention_convention_interaction_effect_mean/2)/2 ) )
  )
  , c("Probability of Encoding\nAttention Effect\nGiven Don't Know"  , "Probability of Encoding\nAttention Effect\nGiven Know")
  , y_lab = "\u03C1 (Attended - Unattended)"
)
# Know condition estimate 
get_95_HDI(
  ( plogis(ex_toj_color_post$population_logit_rho_intercept_mean - ex_toj_color_post$population_logit_rho_convention_effect_mean/2
           + (ex_toj_color_post$population_logit_rho_attention_effect_mean - ex_toj_color_post$population_logit_rho_attention_convention_interaction_effect_mean/2)/2 )
    - plogis(ex_toj_color_post$population_logit_rho_intercept_mean - ex_toj_color_post$population_logit_rho_convention_effect_mean/2
             - (ex_toj_color_post$population_logit_rho_attention_effect_mean -  ex_toj_color_post$population_logit_rho_attention_convention_interaction_effect_mean/2)/2 ) )
)

# A + DK
get_95_HDI(
  plogis(ex_toj_color_post$population_logit_rho_intercept_mean + ex_toj_color_post$population_logit_rho_convention_effect_mean/2
                     + (ex_toj_color_post$population_logit_rho_attention_effect_mean + ex_toj_color_post$population_logit_rho_attention_convention_interaction_effect_mean/2)/2 )
)
# U + DK
get_95_HDI(
  plogis(ex_toj_color_post$population_logit_rho_intercept_mean + ex_toj_color_post$population_logit_rho_convention_effect_mean/2
         - (ex_toj_color_post$population_logit_rho_attention_effect_mean +  ex_toj_color_post$population_logit_rho_attention_convention_interaction_effect_mean/2)/2 ) 
)
# A + K
get_95_HDI(
  plogis(ex_toj_color_post$population_logit_rho_intercept_mean - ex_toj_color_post$population_logit_rho_convention_effect_mean/2
         + (ex_toj_color_post$population_logit_rho_attention_effect_mean - ex_toj_color_post$population_logit_rho_attention_convention_interaction_effect_mean/2)/2 )
)
# U + K
get_95_HDI(
  plogis(ex_toj_color_post$population_logit_rho_intercept_mean - ex_toj_color_post$population_logit_rho_convention_effect_mean/2
         - (ex_toj_color_post$population_logit_rho_attention_effect_mean -  ex_toj_color_post$population_logit_rho_attention_convention_interaction_effect_mean/2)/2 ) 
  )

# BETTER:
get_violin(
  c(
    ( (ex_toj_color_post$logitRhoMean + ex_toj_color_post$logitRhoConventionEffectMean/2
             + (ex_toj_color_post$logitRhoEffectMean + ex_toj_color_post$logitRhoConventionInteractionEffectMean)/2 )
      - (ex_toj_color_post$logitRhoMean + ex_toj_color_post$logitRhoConventionEffectMean/2
               - (ex_toj_color_post$logitRhoEffectMean +  ex_toj_color_post$logitRhoConventionInteractionEffectMean)/2 ) )
    , ( (ex_toj_color_post$logitRhoMean - ex_toj_color_post$logitRhoConventionEffectMean/2
               + (ex_toj_color_post$logitRhoEffectMean - ex_toj_color_post$logitRhoConventionInteractionEffectMean)/2 )
        - (ex_toj_color_post$logitRhoMean - ex_toj_color_post$logitRhoConventionEffectMean/2
                 - (ex_toj_color_post$logitRhoEffectMean -  ex_toj_color_post$logitRhoConventionInteractionEffectMean)/2 ) )
  )
  , c("Probability of Encoding\nAttention Effect\nGiven Don't Know"  , "Probability of Encoding\nAttention Effect\nGiven Know")
  , y_lab = "Log-odds of \u03C1 (Attended - Unattended)"
)
# Know condition estimate 
get_95_HDI(
  ( (ex_toj_color_post$logitRhoMean - ex_toj_color_post$logitRhoConventionEffectMean/2
           + (ex_toj_color_post$logitRhoEffectMean - ex_toj_color_post$logitRhoConventionInteractionEffectMean)/2 )
    - (ex_toj_color_post$logitRhoMean - ex_toj_color_post$logitRhoConventionEffectMean/2
             - (ex_toj_color_post$logitRhoEffectMean -  ex_toj_color_post$logitRhoConventionInteractionEffectMean)/2 ) )
)

# A + DK
get_95_HDI(
  (ex_toj_color_post$logitRhoMean + ex_toj_color_post$logitRhoConventionEffectMean/2
         + (ex_toj_color_post$logitRhoEffectMean + ex_toj_color_post$logitRhoConventionInteractionEffectMean)/2 )
)
# U + DK
get_95_HDI(
  (ex_toj_color_post$logitRhoMean + ex_toj_color_post$logitRhoConventionEffectMean/2
         - (ex_toj_color_post$logitRhoEffectMean +  ex_toj_color_post$logitRhoConventionInteractionEffectMean)/2 ) 
)
# A + K
get_95_HDI(
  (ex_toj_color_post$logitRhoMean - ex_toj_color_post$logitRhoConventionEffectMean/2
         + (ex_toj_color_post$logitRhoEffectMean - ex_toj_color_post$logitRhoConventionInteractionEffectMean)/2 )
)
# U + K
get_95_HDI(
  (ex_toj_color_post$logitRhoMean - ex_toj_color_post$logitRhoConventionEffectMean/2
         - (ex_toj_color_post$logitRhoEffectMean -  ex_toj_color_post$logitRhoConventionInteractionEffectMean)/2 ) 
)

#-------------------------------- Rho Convention Effects ----------------------------------#



#------------------------------------------------------------------------------------------#
#--------------------------------- Kappa Scale ----------------------------------------------#
#------------------------------------------------------------------------------------------#

#---------------------------------- Kappa Conditions-----------------------------------------#
# DK + A
Kappa_pos_attn_pos_conv =
  exp(ex_toj_color_post$logKappaMean + ex_toj_color_post$logKappaConventionEffectMean/2 
         + (ex_toj_color_post$logKappaEffectMean + ex_toj_color_post$logKappaConventionInteractionEffectMean)/2)
# DK + U
Kappa_neg_attn_pos_conv =
  exp(ex_toj_color_post$logKappaMean + ex_toj_color_post$logKappaConventionEffectMean/2 
         - (ex_toj_color_post$logKappaEffectMean + ex_toj_color_post$logKappaConventionInteractionEffectMean)/2)  
# K + A
Kappa_pos_attn_neg_conv =
  exp(ex_toj_color_post$logKappaMean -  ex_toj_color_post$logKappaConventionEffectMean/2 
         + (ex_toj_color_post$logKappaEffectMean - ex_toj_color_post$logKappaConventionInteractionEffectMean)/2)
# K + U
Kappa_neg_attn_neg_conv =
  exp(ex_toj_color_post$logKappaMean -  ex_toj_color_post$logKappaConventionEffectMean/2 
         - (ex_toj_color_post$logKappaEffectMean - ex_toj_color_post$logKappaConventionInteractionEffectMean)/2)
#---------------------------------- Kappa Conditions-----------------------------------------#



#---------------------------------- Kappa Intercept -----------------------------------------#
# BAD:
get_violin(
  exp(ex_toj_color_post$logKappaMean)
  , "Fidelity of Encoding Intercept Mean"
  , y_lab = "\u03BA"
  , hline = FALSE
)

# GOOD:
get_violin(
  ( Kappa_neg_attn_neg_conv
    + Kappa_neg_attn_pos_conv
    + Kappa_pos_attn_neg_conv
    + Kappa_pos_attn_pos_conv
  )/4
  , "Fidelity of Encoding Intercept Mean"
  , y_lab = "\u03BA"
  , hline = FALSE
)

# BETTER:
get_violin(
  ex_toj_color_post$logKappaMean
  , "Fidelity of Encoding Intercept Mean"
  , y_lab = "Log of \u03BA"
  , hline = FALSE
)
#---------------------------------- Kappa Intercept -----------------------------------------#


#-------------------------------- Kappa Attention Effect ------------------------------------#
# BAD:
get_violin(
  ( exp(ex_toj_color_post$logKappaMean + ex_toj_color_post$logKappaEffectMean/2 )
    - exp(ex_toj_color_post$logKappaMean - ex_toj_color_post$logKappaEffectMean/2 ) )
  , "Fidelity of Encoding Effect Mean"
  , y_lab = "\u03BA (Attended - Unattended)"
)

# GOOD:
get_violin(
  ( Kappa_pos_attn_neg_conv - Kappa_neg_attn_neg_conv
    + Kappa_pos_attn_pos_conv - Kappa_neg_attn_pos_conv
  )/2
  , "Fidelity of Encoding Effect Mean"
  , y_lab = "\u03BA (Attended - Unattended)"
)

# BETTER:
get_violin(
  ex_toj_color_post$logKappaEffectMean
  , "Fidelity of Encoding Effect Mean"
  , y_lab = "Log of \u03BA (Attended - Unattended)"
)
#-------------------------------- Kappa Attention Effect ------------------------------------#


#-------------------------------- Kappa Convention Effects ----------------------------------#
### BAD:
get_violin(
  ( exp(ex_toj_color_post$population_log_kappa_intercept_mean + ex_toj_color_post$population_log_kappa_convention_effect_mean/2 )
    - exp(ex_toj_color_post$population_log_kappa_intercept_mean - ex_toj_color_post$population_log_kappa_convention_effect_mean/2 ) )
  , "Fidelity of Encoding\nConvention Effect Mean"
  , y_lab = "\u03BA (Don't Know - Know)"
)
# DK
get_95_HDI(
  exp(ex_toj_color_post$logKappaMean + ex_toj_color_post$logKappaConventionEffectMean/2)
)
# K
get_95_HDI(
  exp(ex_toj_color_post$logKappaMean - ex_toj_color_post$logKappaConventionEffectMean/2 )
)

# GOOD:
get_violin(
  ( Kappa_pos_attn_pos_conv - Kappa_pos_attn_neg_conv 
    + Kappa_neg_attn_pos_conv - Kappa_neg_attn_neg_conv
  )/2
  , "Fidelity of Encoding\nConvention Effect Mean"
  , y_lab = "\u03BA (Don't Know - Know)"
)
# DK
get_95_HDI(
  (Kappa_neg_attn_pos_conv + Kappa_pos_attn_pos_conv)/2
)
# K
get_95_HDI(
  (Kappa_neg_attn_neg_conv + Kappa_pos_attn_neg_conv)/2
)

# BETTER:
get_violin(
  ( (ex_toj_color_post$logKappaMean + ex_toj_color_post$logKappaConventionEffectMean/2 )
    - (ex_toj_color_post$logKappaMean - ex_toj_color_post$logKappaConventionEffectMean/2 ) )
  , "Fidelity of Encoding\nConvention Effect Mean"
  , y_lab = "Log of \u03BA (Don't Know - Know)"
)
# DK
get_95_HDI(
  (ex_toj_color_post$logKappaMean + ex_toj_color_post$logKappaConventionEffectMean/2)
)
# K
get_95_HDI(
  (ex_toj_color_post$logKappaMean - ex_toj_color_post$logKappaConventionEffectMean/2 )
)


### BAD: interaction
get_violin(
  ( exp(ex_toj_color_post$population_log_kappa_intercept_mean + ex_toj_color_post$population_log_kappa_attention_convention_interaction_effect_mean/2 )
    - exp(ex_toj_color_post$population_log_kappa_intercept_mean - ex_toj_color_post$population_log_kappa_attention_convention_interaction_effect_mean/2 ) )
  , "Fidelity of Encoding\nConvention Knowledge\nInteraction Effect Mean"
  , y_lab = "\u03BA"
)

# GOOD: interaction
get_violin(
  ( Kappa_pos_attn_pos_conv - Kappa_neg_attn_pos_conv )
  - ( Kappa_pos_attn_neg_conv - Kappa_neg_attn_neg_conv )
  , "Fidelity of Encoding\nConvention Knowledge\nInteraction Effect Mean"
  , y_lab = "\u03BA"
)

# BETTER: interaction
get_violin(
  2 * ex_toj_color_post$logKappaConventionInteractionEffectMean
  , "Fidelity of Encoding\nConvention Knowledge\nInteraction Effect Mean"
  , y_lab = "Log of \u03BA"
)


### GOOD: effect of attention on Kappa for each convention knowledge level
get_violin(
  c(
    ( exp(ex_toj_color_post$logKappaMean + ex_toj_color_post$logKappaConventionEffectMean/2
             + (ex_toj_color_post$logKappaEffectMean + ex_toj_color_post$logKappaConventionInteractionEffectMean)/2 )
      - exp(ex_toj_color_post$logKappaMean + ex_toj_color_post$logKappaConventionEffectMean/2
               - (ex_toj_color_post$logKappaEffectMean +  ex_toj_color_post$logKappaConventionInteractionEffectMean)/2 ) )
    , ( exp(ex_toj_color_post$logKappaMean - ex_toj_color_post$logKappaConventionEffectMean/2
               + (ex_toj_color_post$logKappaEffectMean - ex_toj_color_post$logKappaConventionInteractionEffectMean)/2 )
        - exp(ex_toj_color_post$logKappaMean - ex_toj_color_post$logKappaConventionEffectMean/2
                 - (ex_toj_color_post$logKappaEffectMean -  ex_toj_color_post$logKappaConventionInteractionEffectMean)/2 ) )
  )
  , c("Fidelity of Encoding\nAttention Effect\nGiven Don't Know"  , "Fidelity of Encoding\nAttention Effect\nGiven Know")
  , y_lab = "\u03BA (Attended - Unattended)"
)
# Know condition estimate 
get_95_HDI(
  ( exp(ex_toj_color_post$logKappaMean - ex_toj_color_post$logKappaConventionEffectMean/2
           + (ex_toj_color_post$logKappaEffectMean - ex_toj_color_post$logKappaConventionInteractionEffectMean)/2 )
    - exp(ex_toj_color_post$logKappaMean - ex_toj_color_post$logKappaConventionEffectMean/2
             - (ex_toj_color_post$logKappaEffectMean -  ex_toj_color_post$logKappaConventionInteractionEffectMean)/2 ) )
)

# A + DK
get_95_HDI(
  exp(ex_toj_color_post$logKappaMean + ex_toj_color_post$logKappaConventionEffectMean/2
         + (ex_toj_color_post$logKappaEffectMean + ex_toj_color_post$logKappaConventionInteractionEffectMean)/2 )
)
# U + DK
get_95_HDI(
  exp(ex_toj_color_post$logKappaMean + ex_toj_color_post$logKappaConventionEffectMean/2
         - (ex_toj_color_post$logKappaEffectMean +  ex_toj_color_post$logKappaConventionInteractionEffectMean)/2 ) 
)
# A + K
get_95_HDI(
  exp(ex_toj_color_post$logKappaMean - ex_toj_color_post$logKappaConventionEffectMean/2
         + (ex_toj_color_post$logKappaEffectMean - ex_toj_color_post$logKappaConventionInteractionEffectMean)/2 )
)
# U + K
get_95_HDI(
  exp(ex_toj_color_post$logKappaMean - ex_toj_color_post$logKappaConventionEffectMean/2
         - (ex_toj_color_post$logKappaEffectMean -  ex_toj_color_post$logKappaConventionInteractionEffectMean)/2 ) 
)

# BETTER: effect of attention on Kappa for each convention knowledge level
get_violin(
  c(
    ( (ex_toj_color_post$logKappaMean + ex_toj_color_post$logKappaConventionEffectMean/2
          + (ex_toj_color_post$logKappaEffectMean + ex_toj_color_post$logKappaConventionInteractionEffectMean)/2 )
      - (ex_toj_color_post$logKappaMean + ex_toj_color_post$logKappaConventionEffectMean/2
            - (ex_toj_color_post$logKappaEffectMean +  ex_toj_color_post$logKappaConventionInteractionEffectMean)/2 ) )
    , ( (ex_toj_color_post$logKappaMean - ex_toj_color_post$logKappaConventionEffectMean/2
            + (ex_toj_color_post$logKappaEffectMean - ex_toj_color_post$logKappaConventionInteractionEffectMean)/2 )
        - (ex_toj_color_post$logKappaMean - ex_toj_color_post$logKappaConventionEffectMean/2
              - (ex_toj_color_post$logKappaEffectMean -  ex_toj_color_post$logKappaConventionInteractionEffectMean)/2 ) )
  )
  , c("Fidelity of Encoding\nAttention Effect\nGiven Don't Know"  , "Fidelity of Encoding\nAttention Effect\nGiven Know")
  , y_lab = "Log of \u03BA (Attended - Unattended)"
)
# Know condition estimate 
get_95_HDI(
  ( (ex_toj_color_post$logKappaMean - ex_toj_color_post$logKappaConventionEffectMean/2
        + (ex_toj_color_post$logKappaEffectMean - ex_toj_color_post$logKappaConventionInteractionEffectMean)/2 )
    - (ex_toj_color_post$logKappaMean - ex_toj_color_post$logKappaConventionEffectMean/2
          - (ex_toj_color_post$logKappaEffectMean -  ex_toj_color_post$logKappaConventionInteractionEffectMean)/2 ) )
)

# A + DK
get_95_HDI(
  (ex_toj_color_post$logKappaMean + ex_toj_color_post$logKappaConventionEffectMean/2
      + (ex_toj_color_post$logKappaEffectMean + ex_toj_color_post$logKappaConventionInteractionEffectMean)/2 )
)
# U + DK
get_95_HDI(
  (ex_toj_color_post$logKappaMean + ex_toj_color_post$logKappaConventionEffectMean/2
      - (ex_toj_color_post$logKappaEffectMean +  ex_toj_color_post$logKappaConventionInteractionEffectMean)/2 ) 
)
# A + K
get_95_HDI(
  (ex_toj_color_post$logKappaMean - ex_toj_color_post$logKappaConventionEffectMean/2
      + (ex_toj_color_post$logKappaEffectMean - ex_toj_color_post$logKappaConventionInteractionEffectMean)/2 )
)
# U + K
get_95_HDI(
  (ex_toj_color_post$logKappaMean - ex_toj_color_post$logKappaConventionEffectMean/2
      - (ex_toj_color_post$logKappaEffectMean -  ex_toj_color_post$logKappaConventionInteractionEffectMean)/2 ) 
)
#-------------------------------- Kappa Convention Effects ----------------------------------#



#------------------------------------------------------------------------------------------#
#--------------------------------- Graphs -------------------------------------------------#
#------------------------------------------------------------------------------------------#

#---------------------------------- NCFs --------------------------------------------------#
### BAD:
yGloveDontKnow = pnorm(
  -250:250
  , mean = ( 
    median(ex_toj_color_post$population_pss_intercept_mean) + median(ex_toj_color_post$population_pss_convention_effect_mean)/2 
    + (
     median(ex_toj_color_post$population_pss_effect_mean) 
     + median(ex_toj_color_post$population_pss_convention_interaction_effect_mean/2)
     )/2 
  )* 250
  , sd = ( exp( 
    median(ex_toj_color_post$population_log_jnd_intercept_mean) + median(ex_toj_color_post$population_log_jnd_convention_effect_mean)/2 
    + (
      median(ex_toj_color_post$population_log_jnd_effect_mean) 
      + median(ex_toj_color_post$population_log_jnd_convention_interaction_effect_mean/2)
      ) /2
  ) ) * 250
)

yGloveKnow = pnorm(
  -250:250
  , mean = ( 
    median(ex_toj_color_post$population_pss_intercept_mean) - median(ex_toj_color_post$population_pss_convention_effect_mean)/2 
    + (
      median(ex_toj_color_post$population_pss_effect_mean) 
      - median(ex_toj_color_post$population_pss_convention_interaction_effect_mean/2)
    )/2 
  )* 250
  , sd = ( exp( 
    median(ex_toj_color_post$population_log_jnd_intercept_mean) - median(ex_toj_color_post$population_log_jnd_convention_effect_mean)/2 
    + (
      median(ex_toj_color_post$population_log_jnd_effect_mean) 
      - median(ex_toj_color_post$population_log_jnd_convention_interaction_effect_mean/2)
    ) /2
  ) ) * 250
)

yBaseDontKnow = pnorm(
  -250:250
  , mean = ( 
    median(ex_toj_color_post$population_pss_intercept_mean) + median(ex_toj_color_post$population_pss_convention_effect_mean)/2 
    - (
      median(ex_toj_color_post$population_pss_effect_mean) 
      + median(ex_toj_color_post$population_pss_convention_interaction_effect_mean/2)
    )/2 
  )* 250
  , sd = ( exp( 
    median(ex_toj_color_post$population_log_jnd_intercept_mean) + median(ex_toj_color_post$population_log_jnd_convention_effect_mean)/2 
    - (
      median(ex_toj_color_post$population_log_jnd_effect_mean) 
      + median(ex_toj_color_post$population_log_jnd_convention_interaction_effect_mean/2)
    ) /2
  ) ) * 250
)

yBaseKnow = pnorm(
  -250:250
  , mean = ( 
    median(ex_toj_color_post$population_pss_intercept_mean) - median(ex_toj_color_post$population_pss_convention_effect_mean)/2 
    - (
      median(ex_toj_color_post$population_pss_effect_mean) 
      - median(ex_toj_color_post$population_pss_convention_interaction_effect_mean/2)
    )/2 
  )* 250
  , sd = ( exp( 
    median(ex_toj_color_post$population_log_jnd_intercept_mean) - median(ex_toj_color_post$population_log_jnd_convention_effect_mean)/2 
    - (
      median(ex_toj_color_post$population_log_jnd_effect_mean) 
      - median(ex_toj_color_post$population_log_jnd_convention_interaction_effect_mean/2)
    ) /2
  ) ) * 250
)

# GOOD:
yGloveDontKnow = pnorm(
  -250:250
  , mean = median(
    ( ex_toj_color_post$population_pss_intercept_mean + ex_toj_color_post$population_pss_convention_effect_mean/2 
    + (
      ex_toj_color_post$population_pss_effect_mean
      + ex_toj_color_post$population_pss_convention_interaction_effect_mean
    )/2 ) * 250
  )
  , sd = median( exp( 
    ex_toj_color_post$population_logjnd_intercept_mean + ex_toj_color_post$population_logjnd_convention_effect_mean/2 
    + (
     ex_toj_color_post$population_logjnd_effect_mean
      + ex_toj_color_post$population_logjnd_convention_interaction_effect_mean
    ) /2 ) * 250 
  )
)

yGloveKnow = pnorm(
  -250:250
  , mean = median(
    ( ex_toj_color_post$population_pss_intercept_mean - ex_toj_color_post$population_pss_convention_effect_mean/2 
      + (
        ex_toj_color_post$population_pss_effect_mean
        - ex_toj_color_post$population_pss_convention_interaction_effect_mean
      )/2 ) * 250
  )
  , sd = median( exp( 
    ex_toj_color_post$population_logjnd_intercept_mean - ex_toj_color_post$population_logjnd_convention_effect_mean/2 
    + (
      ex_toj_color_post$population_logjnd_effect_mean
      - ex_toj_color_post$population_logjnd_convention_interaction_effect_mean
    ) /2 ) * 250 
  )
)

yBaseDontKnow = pnorm(
  -250:250
  , mean = median(
    ( ex_toj_color_post$population_pss_intercept_mean + ex_toj_color_post$population_pss_convention_effect_mean/2 
      - (
        ex_toj_color_post$population_pss_effect_mean
        + ex_toj_color_post$population_pss_convention_interaction_effect_mean
      )/2 ) * 250
  )
  , sd = median( exp( 
    ex_toj_color_post$population_logjnd_intercept_mean + ex_toj_color_post$population_logjnd_convention_effect_mean/2 
    - (
      ex_toj_color_post$population_logjnd_effect_mean
      + ex_toj_color_post$population_logjnd_convention_interaction_effect_mean
    ) /2 ) * 250 
  )
)

yBaseKnow = pnorm(
  -250:250
  , mean = median(
    ( ex_toj_color_post$population_pss_intercept_mean - ex_toj_color_post$population_pss_convention_effect_mean/2 
      - (
        ex_toj_color_post$population_pss_effect_mean
        - ex_toj_color_post$population_pss_convention_interaction_effect_mean
      )/2 ) * 250
  )
  , sd = median( exp( 
    ex_toj_color_post$population_logjnd_intercept_mean - ex_toj_color_post$population_logjnd_convention_effect_mean/2 
    - (
      ex_toj_color_post$population_logjnd_effect_mean
      - ex_toj_color_post$population_logjnd_convention_interaction_effect_mean
    ) /2 ) * 250 
  )
)

df = data.frame(SOA = -250:250
                , Prop = c(yGloveDontKnow, yGloveKnow, yBaseDontKnow, yBaseKnow)
                , Attend = c(rep("Glove",1002), rep("Base", 1002))
                , Convention = c(rep( c(rep("Don't Know", 501), rep("Know", 501)), 2) )
)

toj_means_by_id_by_condition2 = ddply(
  .data = toj_trials
  , .variables = .(base_probe_dist, know_tie_goes_runner, soa2)
  , .fun = function(x){
    to_return = data.frame(
      Prop = mean(x$safe)
      , Attend = paste( unique(x$base_probe_dist) )
      , Convention = paste( unique(x$know_tie_goes_runner) )
      , SOA = unique(x$soa2)
    )
    return(to_return)
  }
)

levels(toj_means_by_id_by_condition2$Attend) = c("Glove", "Base")
levels(toj_means_by_id_by_condition2$Convention) = c("Don't Know", "Know")


gg = ggplot(data = df, aes(y = Prop, x = SOA, colour = Attend))+
  geom_line(size = 1.25)+
  # scale_color_manual("Attend", values = c("red", "blue"))+
  scale_color_hue("Attend", l = c(60, 15), c = c(100, 50), h = c(240, 360) ) +
  labs(x = "SOA (ms)", y = "Proportion of 'Safe' Responses")+
  geom_point(data = toj_means_by_id_by_condition2, aes(y = Prop, x = SOA, colour = Attend), size = 4)+
  geom_point(data = toj_means_by_id_by_condition2, aes(y = Prop, x = SOA), size = 2.5, colour = "grey90")+
  facet_grid(Convention~.)+
  theme_gray(base_size = 30)+
  theme(panel.grid.major = element_line(size = 1.5)
        ,panel.grid.minor = element_line(size = 1))
  # define text to add
  Text1 = textGrob(label = paste("Out"), gp = gpar(fontsize= 30))
  Text2 = textGrob(label = paste("Safe"), gp = gpar(fontsize= 30)) 
  gg = gg+
  annotation_custom(grob = Text1,  xmin = -200, xmax = -200, ymin = -0.25, ymax = -0.25)+
  annotation_custom(grob = Text2,  xmin = 200, xmax = 200, ymin = -0.25, ymax = -0.25)
  # Code to override clipping
  gg2 <- ggplot_gtable(ggplot_build(gg))
  gg2$layout$clip[gg2$layout$name=="panel"] <- "off"
  grid.draw(gg2)
#---------------------------------- NCFs --------------------------------------------------#

