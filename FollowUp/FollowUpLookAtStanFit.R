# library(shinystan)
library(coda)
library(ggplot2)
library(ggmcmc)
library(CircStats)
library(grid)
library(sprintfr)
library(plyr)
library(gtools)

setwd("~/Documents/TOJ/Follow-Up")
# load("FollowUptoj_color_post_June21th2016")
load("FollowUp_color_trials.Rdata")
load("FollowUp_toj_trials.Rdata")
source("../EndogenousVisualPriorEntry-BayesianHierarchicalModel/functions.R")


############################################################################################
####                                        Diagnostics                                 ####
############################################################################################
# convert stanfit sample to dataframe table 
gg_toj_color_post = ggs(toj_color_post)

# list of parameters to examine
param_list = c("logitRhoAttentionEffectMean"
               , "logitRhoMean"
               , "logitRhoJudgementTypeEffectMean"
               , "logitRhoJudgementTypeInteractionEffectMean"
               , "logitRhoInitialBiasEffectMean"
               , "logitRhoInitialBiasInteractionEffectMean"
               , "logKappaAttentionEffectMean"
               , "logKappaMean"
               , "logKappaJudgementTypeEffectMean"
               , "logKappaJudgementTypeInteractionEffectMean"
               , "logKappaInitialBiasEffectMean"
               , "logKappaInitialBiasInteractionEffectMean"
               , "population_logjnd_intercept_mean"
               , "population_logjnd_effect_mean"
               , "population_logjnd_initial_bias_effect_mean"
               , "population_logjnd_initial_bias_interaction_effect_mean"
               , "population_logjnd_judgement_type_effect_mean"
               , "population_logjnd_judgement_type_interaction_effect_mean"
               , "population_pss_intercept_mean"
               , "population_pss_effect_mean"
               , "population_pss_initial_bias_effect_mean"
               , "population_pss_initial_bias_interaction_effect_mean"
               , "population_pss_judgement_type_effect_mean"
               , "population_pss_judgement_type_interaction_effect_mean"
               , "zlogitRhoEffectSD" 
               , "zlogitRhoSD"
               , "zlogKappaEffectSD"
               , "zlogKappaSD"
               , "zpopulation_logjnd_effect_sd"
               , "zpopulation_logjnd_intercept_sd"
               , "zpopulation_pss_effect_sd"
               , "zpopulation_pss_intercept_sd")

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
real_toj = aggregate(left_first_TF ~ soa2 + block_bias + toj_judgement_type + probe_initial_bias, data = toj_trials, FUN = mean)
#-------------------------------------- TOJ Actual Data -----------------------------------#


#-------------------------------------- TOJ Simulated Data --------------------------------#

### Get PSS Parameters
pss_intercept_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "population_pss_intercept_mean",]$value
pss_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "population_pss_effect_mean",]$value

pss_judgement_type_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "population_pss_judgement_type_effect_mean",]$value
pss_judgement_type_interaction_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "population_pss_judgement_type_interaction_effect_mean",]$value

pss_initial_bias_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "population_pss_initial_bias_effect_mean",]$value
pss_initial_bias_interaction_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "population_pss_initial_bias_interaction_effect_mean",]$value

# judgement type
pss_right_first_right_mean_reps = get_condition_mean_sample(
  ( pss_intercept_mean 
    - pss_judgement_type_effect_mean/2 
    - pss_initial_bias_effect_mean/2 )
  , ( pss_effect_mean 
      - pss_judgement_type_interaction_effect_mean 
      - pss_initial_bias_interaction_effect_mean)
  , TRUE
  , "null"
)

### Get JND Parameters
logjnd_intercept_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "population_logjnd_intercept_mean",]$value
logjnd_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "population_logjnd_effect_mean",]$value

logjnd_judgement_type_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "population_logjnd_judgement_type_effect_mean",]$value
logjnd_judgement_type_interaction_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "population_logjnd_judgement_type_interaction_effect_mean",]$value

logjnd_initial_bias_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "population_logjnd_initial_bias_effect_mean",]$value
logjnd_initial_bias_interaction_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "population_logjnd_initial_bias_interaction_effect_mean",]$value

# judgement type
logjnd_right_first_right_short_mean_reps = get_condition_mean_sample(
  ( logjnd_intercept_mean 
  - logjnd_judgement_type_effect_mean/2 
  - logjnd_initial_bias_effect_mean/2 )
  , ( logjnd_effect_mean 
     - logjnd_judgement_type_interaction_effect_mean 
     - logjnd_initial_bias_interaction_effect_mean )
  , TRUE
  , "log"
)
#-------------------------------------- TOJ Simulated Data --------------------------------#


#-------------------------------------- Do TOJ PPC ----------------------------------------#
SOAs = c(-250, -150, -100, -50, -17, 17, 50, 100, 150, 250)

# judgement type
do_toj_ppc(
  pss_right_first_right_short_mean_reps
  , logjnd_right_first_right_short_mean_reps
  , "'which first?' & initial right & short duration & attend right"
  , c("toj_judgement_type", "probe_initial_bias","onehundredms" , "block_bias")
  , c("first", "RIGHT", "TRUE","RIGHT")
  , "left proportion"
  , real = real_toj
)
#-------------------------------------- Do TOJ PPC ----------------------------------------#


#-------------------------------------- Color Actual Data ---------------------------------#
hist(color_trials[color_trials$attended == TRUE,]$color_diff_radians, breaks = 30, freq = F, col = rgb(.1,.1,.1,.5))
hist(color_trials[color_trials$attended == FALSE,]$color_diff_radians, breaks = 30, freq = F, col = rgb(.9,.9,.9,.5), add = T)
#-------------------------------------- Color Actual Data ---------------------------------#


#-------------------------------------- Color Simulated Data ------------------------------#
### Get Rho Parameters
rho_intercept_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "logitRhoMean",]$value
rho_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "logitRhoEffectMean",]$value

rho_judgement_type_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "logitRhoJudgementTypeEffectMean",]$value
rho_judgement_type_interaction_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "logitRhoJudgementTypeInteractionEffectMean",]$value

rho_initial_bias_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "logitRhoInitialBiasEffectMean",]$value
rho_initial_bias_interaction_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "logitRhoInitialBiasInteractionEffectMean",]$value

rho_probe_duration_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "logitRhoProbeEffectMean",]$value
rho_probe_duration_interaction_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "logitRhoProbeInteractionEffectMean",]$value

rho_initial_bias_judgement_type_interaction_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "logitRhoInitialBiasJudgementTypeInteractionEffectMean",]$value
rho_initial_bias_probe_interaction_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "logitRhoInitialBiasProbeInteractionEffectMean",]$value
rho_judgement_type_probe_interaction_effect_mean =  gg_toj_color_post[gg_toj_color_post$Parameter == "logitRhoJudgementTypeProbeInteractionEffectMean",]$value
rho_four_way_interaction_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "logitRhoFourWayInteractionEffectMean",]$value

# get condition combination
condition_toj_judgement = -1
condition_initial_bias = -1
condition_probe = -1
rho_right_first_right_short_mean_reps = get_condition_mean_sample(
  ( rho_intercept_mean 
    - rho_judgement_type_effect_mean/2 
    - rho_initial_bias_effect_mean/2 
    - rho_probe_duration_effect_mean/2 )
  , ( rho_effect_mean 
      - rho_judgement_type_interaction_effect_mean 
      - rho_initial_bias_interaction_effect_mean
      - rho_probe_duration_interaction_effect_mean
      + condition_toj_judgement * condition_initial_bias * rho_initial_bias_judgement_type_interaction_mean
      + condition_toj_judgement * condition_probe * rho_judgement_type_probe_interaction_effect_mean
      + condition_initial_bias * condition_probe * rho_initial_bias_probe_interaction_mean
      + condition_initial_bias * condition_probe * condition_toj_judgement * rho_four_way_interaction_effect_mean)
  , TRUE
  , "logit"
)

### Get Kappa Parameters
kappa_intercept_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "logKappaMean",]$value
kappa_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "logKappaEffectMean",]$value

kappa_judgement_type_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "logKappaJudgementTypeEffectMean",]$value
kappa_judgement_type_interaction_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "logKappaJudgementTypeInteractionEffectMean",]$value

kappa_initial_bias_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "logKappaInitialBiasEffectMean",]$value
kappa_initial_bias_interaction_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "logKappaInitialBiasInteractionEffectMean",]$value

kappa_probe_duration_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "logKappaProbeEffectMean",]$value
kappa_probe_duration_interaction_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "logKappaProbeInteractionEffectMean",]$value

kappa_initial_bias_judgement_type_interaction_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "logKappaInitialBiasJudgementTypeInteractionEffectMean",]$value
kappa_initial_bias_probe_interaction_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "logKappaInitialBiasProbeInteractionEffectMean",]$value
kappa_judgement_type_probe_interaction_effect_mean =  gg_toj_color_post[gg_toj_color_post$Parameter == "logKappaJudgementTypeProbeInteractionEffectMean",]$value
kappa_four_way_interaction_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "logKappaFourWayInteractionEffectMean",]$value

# get condition combination
condition_toj_judgement = -1
condition_initial_bias = -1
condition_probe = -1
kappa_right_first_right_short_mean_reps = get_condition_mean_sample(
  ( kappa_intercept_mean 
    - kappa_judgement_type_effect_mean/2 
    - kappa_initial_bias_effect_mean/2 
    - kappa_probe_duration_effect_mean/2 )
  , ( kappa_effect_mean 
      - kappa_judgement_type_interaction_effect_mean 
      - kappa_initial_bias_interaction_effect_mean
      - kappa_probe_duration_interaction_effect_mean
      + condition_toj_judgement * condition_initial_bias * kappa_initial_bias_judgement_type_interaction_mean
      + condition_toj_judgement * condition_probe * kappa_judgement_type_probe_interaction_effect_mean
      + condition_initial_bias * condition_probe * kappa_initial_bias_probe_interaction_mean
      + condition_initial_bias * condition_probe * condition_toj_judgement * kappa_four_way_interaction_effect_mean)
  , TRUE
  , "log_free"
)
#-------------------------------------- Color Simulated Data ------------------------------#


#-------------------------------------- Do Color PPC --------------------------------------#
do_color_ppc(
  rho_right_first_right_short_mean_reps
  , kappa_right_first_right_short_mean_reps
  , "'which first?' & initial right & short duration & attend right"
  , c("toj_judgement_type", "probe_initial_bias","onehundredms" , "block_bias")
  , c("first", "RIGHT", "TRUE","RIGHT")
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
# (3) population_logjnd_intercept_mean    
# (4) population_logjnd_effect_mean     
# (5) logitRhoMean                         
# (6) logKappaMean                        
# (7) logitRhoEffectMean                 
# (8) logKappaEffectMean     

# NOTE: for quick look
# not necessarily HDI
ggs_caterpillar(gg_toj_color_post, family = "cor", thick_ci = c(0.25, 0.75) ) + geom_vline(xintercept = 0, col = "red")


#-------------------------------------- Get Betas -----------------------------------------#
# extract samples
detach('package:rstan', unload = T)  # to ensure 
library(rstan)
ex_toj_color_post = extract(toj_color_post)

# for violin plots later
pos_corr2 = data.frame(value = ex_toj_color_post$cor)
pos_corr2$id = rownames(pos_corr2)
pos_corr = melt( pos_corr2 )
names(pos_corr)[2] = c("parameter")

# (1) population_pss_intercept_mean      
# (2) population_pss_effect_mean          
# (3) population_logjnd_intercept_mean    
# (4) population_logjnd_effect_mean     
# (5) logitRhoMean                         
# (6) logKappaMean                        
# (7) logitRhoEffectMean                 
# (8) logKappaEffectMean 
# JND and PSS intercepts
library(reshape)
betas2 = data.frame(value = ex_toj_color_post$beta)
betas2$iteration = rownames(betas2)
betas = melt( betas2 )
betas$parameter = rep( c(
  "population_pss_intercept_mean"      
  , "population_pss_effect_mean"          
  , "population_logjnd_intercept_mean"    
  , "population_logjnd_effect_mean"     
  , "logitRhoMean"                       
  , "logKappaMean"                       
  , "logitRhoEffectMean"                 
  , "logKappaEffectMean"
)
, times = 1
, each = nrow(betas2)*length(unique(betas$variable))/8  # 8 is number of parameters 
)  
betas$participant = rep(c(1:length(unique(toj_trials$id))), times = 8, each = nrow(betas2))
#-------------------------------------- Get Betas -----------------------------------------#


# #---------------------------- Rho vs. PSS Effects -----------------------------------------#
# psseffect = extract_samples("population_pss_effect_mean")
# 
# psseffectsd = extract_samples("zpopulation_pss_effect_sd", TRUE)
# 
# pssinteraction = extract_samples("population_pss_probe_interaction_effect_mean")
# 
# pssjudgementinteraction = extract_samples("population_pss_judgement_type_interaction_effect_mean")
# 
# pssinitialbiasinteraction = extract_samples("population_pss_initial_bias_interaction_effect_mean")
# 
# probefactor = ifelse(aggregate(onehundredms ~ id, data = color_trials, unique)$onehundred, -1, 1)
# 
# judgementfactor = ifelse(aggregate(toj_judgement_type ~ id, data = toj_trials, FUN = unique)$toj_judgement_type == "first", -1, 1)
# 
# initialbiasfactor = ifelse(aggregate(probe_initial_bias ~ id, data = toj_trials, FUN = unique)$probe_initial_bias == "RIGHT", -1, 1)
# 
# psseffect_ids = ddply(
#   .data = betas
#   , .variables = .(participant)
#   , .fun = function(x){
#     i = unique(x$participant)
#     x_use = x[x$parameter ==  "population_pss_effect_mean",]$value
#     psseffect_use = median(psseffect)  + median(psseffectsd)*median(x_use) + median(pssinteraction)*probefactor[i]+ median(pssjudgementinteraction)*judgementfactor[i]+ median(pssinitialbiasinteraction)*initialbiasfactor[i]
#     df = data.frame(psseffect_use*250, probefactor[i], judgementfactor[i],initialbiasfactor[i])
#     names(df) = c("psseffect", "probefactor", "judgementfactor", "initialbiasfactor")
#     return(df)
#   }
# )
# 
# logitrhoeffect = extract_samples("logitRhoEffectMean")
# 
# logitrhoeffectsd = extract_samples("zlogitRhoEffectSD", TRUE)
# 
# logitrhointeractioneffect = extract_samples("logitRhoProbeInteractionEffectMean")
# 
# logitrhojudgementinteractioneffect = extract_samples("logitRhoJudgementTypeInteractionEffectMean")
# 
# logitrhoinitialbiasinteractioneffect = extract_samples("logitRhoInitialBiasInteractionEffectMean")
# 
# rhoeffect_ids = ddply(
#   .data = betas
#   , .variables = .(participant)
#   , .fun = function(x){
#     i = unique(x$participant)
#     x_use = x[x$parameter == "logitRhoEffectMean",]$value
#     logitrhoeffect_use =  median(logitrhoeffect) + median(logitrhoeffectsd)*median(x_use)+ median(logitrhointeractioneffect)*probefactor[i] + median(logitrhojudgementinteractioneffect)*judgementfactor[i]+ median(logitrhoinitialbiasinteractioneffect)*initialbiasfactor[i]
#     df = data.frame(logitrhoeffect_use, probefactor[i], judgementfactor[i], initialbiasfactor[i])
#     names(df) = c("value", "probefactor", "judgementfactor", "initialbiasfactor")
#     return(df)
#   }
# )
# 
# psseffect_v_rhoeffect = merge(rhoeffect_ids, psseffect_ids)
# 
# ggplot(data = psseffect_v_rhoeffect, aes(y = psseffect, x = value, colour = factor(probefactor), shape = factor(probefactor)))+ #, fill = factor(judgementfactor)))+
#   scale_y_continuous(name = "PSS Effect Mean")+
#   scale_x_continuous(name = "Logit \u03C1 Effect Mean")+
#   geom_vline(xintercept = 0, linetype = 2, size = 1)+
#   geom_hline(yintercept = 0, linetype = 2, size = 1)+
#   geom_point(size = 4)+
#   scale_shape_manual(name = "Probe\nDuration", labels = c("Short", "Long") , values = c(16,17) )+
#   # scale_fill_manual(name = "Judgement\nType", labels = c("Second", "First"), values = c("white", "black"))+
#   # scale_fill_manual(name = "Initial\nBias", labels = c("Left", "Right"), values = c("white", "black"))+
#   scale_colour_manual(name = "Probe\nDuration", labels =c("Short", "Long"), values = c("red", "blue") )+
#   theme_gray(base_size = 30)+
#   theme(panel.grid.major = element_line(size = 1.5)
#         ,panel.grid.minor = element_line(size = 1))
# 
# get_corr(
#   "value.2.7"
#   , "Logit \u03C1 vs. PSS Effect Means"
# )
# #---------------------------- Rho vs. PSS Effects -----------------------------------------#
# 
# 
# #---------------------------- Kappa vs. PSS Effects ---------------------------------------#
# logkappaeffect = extract_samples("logKappaEffectMean")
# 
# logkappaeffectsd = extract_samples("zlogKappaEffectSD", TRUE)
# 
# logkappainteractioneffect = extract_samples("logKappaProbeInteractionEffectMean")
# 
# logkappajudgementtypeinteractioneffect = extract_samples("logKappaJudgementTypeInteractionEffectMean")
# 
# logkappainitialbiasinteractioneffect = extract_samples("logKappaInitialBiasInteractionEffectMean")
# 
# kappaeffect_ids = ddply(
#   .data = betas
#   , .variables = .(participant)
#   , .fun = function(x){
#     i = unique(x$participant)
#     x_use = x[x$parameter == "logKappaEffectMean",]$value
#     logkappaeffect_use =  median(logkappaeffect) + median(logkappaeffectsd)*median(x_use)  +median(logkappainteractioneffect)*probefactor[i] + median(logkappajudgementtypeinteractioneffect)*judgementfactor[i] + median(logkappainitialbiasinteractioneffect)*initialbiasfactor[i]
#     df = data.frame(logkappaeffect_use, probefactor[i], judgementfactor[i], initialbiasfactor[i])
#     names(df) = c("value","probefactor", "judgementfactor", "initialbiasfactor")
#     return(df)
#   }
# )
# 
# psseffect_v_kappaeffect = merge(kappaeffect_ids, psseffect_ids)
# 
# ggplot(data = psseffect_v_kappaeffect, aes(y = psseffect, x = value, colour = factor(probefactor), shape = factor(probefactor)))+ #, fill = factor(initialbiasfactor)))+
#   scale_y_continuous(name = "PSS Effect Mean")+
#   scale_x_continuous(name = "Log \u03BA Effect Mean")+
#   geom_vline(xintercept = 0, linetype = 2, size = 1)+
#   geom_hline(yintercept = 0, linetype = 2, size = 1)+
#   geom_point(size = 4)+
#   scale_shape_manual(name = "Probe\nDuration", labels = c("Short", "Long") , values = c(16,17) )+
#   # scale_fill_manual(name = "Judgement\nType", labels = c("Second", "First"), values = c("white", "black"))+
#   # scale_fill_manual(name = "Initial\nBias", labels = c("Left", "Right"), values = c("white", "black"))+
#   scale_colour_manual(name = "Probe\nDuration", labels =c("Short", "Long"), values = c("red", "blue") )+
#   theme_gray(base_size = 30)+
#   theme(panel.grid.major = element_line(size = 1.5)
#         ,panel.grid.minor = element_line(size = 1))
# 
# get_corr(
#   "value.2.8"
#   , "Log \u03BA vs. PSS Effect Means"
# )        
# #---------------------------- Kappa vs. PSS Effects ---------------------------------------#
# 
# 
# #------------------------- Rho vs. JND Intercept Means ------------------------------------#
# logitrhomean = extract_samples("logitRhoMean")
# 
# logitrhosd = extract_samples("zlogitRhoSD")
# 
# logitrhojudgementeffect = extract_samples("logitRhoJudgementTypeEffectMean")
# 
# logitrhoprobeeffect = extract_samples("logitRhoProbeEffectMean")
# 
# logitrhoinitialbiaseffect = extract_samples("logitRhoInitialBiasEffectMean")
# 
# rhointercept_ids = ddply(
#   .data = betas
#   , .variables = .(participant)
#   , .fun = function(x){
#     i = unique(x$participant)
#     x_use = x[x$parameter == "logitRhoMean",]$value
#     logitrhointercept_use =  median(logitrhomean) + median(logitrhosd)*median(x_use)+ median(logitrhoprobeeffect)*probefactor[i] + median(logitrhojudgementeffect)*judgementfactor[i]+ median(logitrhoinitialbiaseffect)*initialbiasfactor[i]
#     df = data.frame(logitrhointercept_use, probefactor[i], judgementfactor[i], initialbiasfactor[i])
#     names(df) = c("logitrhointercepts", "probefactor", "judgementfactor", "initialbiasfactor")
#     return(df)
#   }
# )
# 
# logjndmean = extract_samples("population_logjnd_intercept_mean")
# 
# logjndsd = extract_samples("zpopulation_logjnd_intercept_sd")
# 
# logjndjudgementeffect = extract_samples("population_logjnd_judgement_type_effect_mean")
# 
# logjndprobeeffect = extract_samples("population_logjnd_probe_effect_mean")
# 
# logjndinitialbiaseffect = extract_samples("population_logjnd_initial_bias_effect_mean")
# 
# jndintercept_ids = ddply(
#   .data = betas
#   , .variables = .(participant)
#   , .fun = function(x){
#     i = unique(x$participant)
#     x_use = x[x$parameter == "population_logjnd_intercept_mean",]$value
#     logjndintercept_use =  median(logjndmean) + median(logjndsd)*median(x_use)+ median(logjndprobeeffect)*probefactor[i] + median(logjndjudgementeffect)*judgementfactor[i]+ median(logjndinitialbiaseffect)*initialbiasfactor[i]
#     df = data.frame(logjndintercept_use, probefactor[i], judgementfactor[i], initialbiasfactor[i])
#     names(df) = c("logjndintercepts", "probefactor", "judgementfactor", "initialbiasfactor")
#     return(df)
#   }
# )
# 
# jndintercept_v_rhointercept = merge(rhointercept_ids, jndintercept_ids)
# 
# ggplot(data = jndintercept_v_rhointercept, aes(y = logjndintercepts, x = logitrhointercepts, colour = factor(probefactor), shape = factor(probefactor)))+ #, fill = factor(judgementfactor)) )+
#   scale_y_continuous(name = "Log JND Intercept Mean")+
#   scale_x_continuous(name = "Logit \u03C1 Intercept Mean")+
#   geom_vline(xintercept = 0, linetype = 2, size = 1)+
#   geom_hline(yintercept = 0, linetype = 2, size = 1)+
#   geom_point(size = 4)+
#   scale_shape_manual(name = "Probe\nDuration", labels = c("Short", "Long") , values = c(21,22) )+
#   # scale_fill_manual(name = "Judgement\nType", labels = c("Second", "First"), values = c("white", "black"))+
#   # scale_fill_manual(name = "Initial\nBias", labels = c("Left", "Right"), values = c("white", "black"))+
#   scale_colour_manual(name = "Probe\nDuration", labels =c("Short", "Long"), values = c("red", "blue") )+
#   theme_gray(base_size = 30)+
#   theme(panel.grid.major = element_line(size = 1.5)
#         ,panel.grid.minor = element_line(size = 1))
# 
# get_corr(
#   "value.3.5"
#   , "Logit \u03C1 vs. JND Intercept Means"
# )  
# #------------------------- Rho vs. JND Intercept Means ------------------------------------#


#------------------------------------------------------------------------------------------#
#--------------------------------- SOA Scale ----------------------------------------------#
#------------------------------------------------------------------------------------------#

#---------------------------------- SOA Intercepts ----------------------------------------#
get_violin(
  c(
    ex_toj_color_post$population_pss_intercept_mean * 250
    , exp( ex_toj_color_post$population_logjnd_intercept_mean ) * 250
    )
  , c("PSS Intercept Mean", "JND Intercept Mean")
  , y_lab = "SOA (ms)"
  , hline = FALSE
  , facet = TRUE
)
#---------------------------------- SOA Intercepts ----------------------------------------#


#---------------------------------- Main Effects ---------------------=--------------------#
# effect of attention on PSS and JND
get_violin(
  c(
  ( (ex_toj_color_post$population_pss_intercept_mean + ex_toj_color_post$population_pss_effect_mean/2) 
    - (ex_toj_color_post$population_pss_intercept_mean - ex_toj_color_post$population_pss_effect_mean/2) ) * 250
  , ( exp( ex_toj_color_post$population_logjnd_intercept_mean + ex_toj_color_post$population_logjnd_effect_mean/2 )
  - exp( ex_toj_color_post$population_logjnd_intercept_mean - ex_toj_color_post$population_logjnd_effect_mean/2  ) ) * 250 
  )
  , c("PSS Effect Mean"  , "JND Effect Mean")
  , y_lab = "SOA (Right - Left; ms)"
)

# effect of judgement type (Q) on PSS and JND
get_violin(
  c(
  ( ex_toj_color_post$population_pss_judgement_type_effect_mean ) * 250 
  , ( exp( ex_toj_color_post$population_logjnd_intercept_mean + ex_toj_color_post$population_logjnd_judgement_type_effect_mean/2 )
      - exp( ex_toj_color_post$population_logjnd_intercept_mean - ex_toj_color_post$population_logjnd_judgement_type_effect_mean/2  ) ) * 250
  )
  , c("PSS Judgement\nType Effect Mean"  , "JND Judgement\nType Effect Mean")
  , y_lab = "SOA (Second - First; ms)"
)

# effect of initial bias on PSS and JND
get_violin(
  c(
  ( ex_toj_color_post$population_pss_initial_bias_effect_mean ) * 250
  , ( exp( ex_toj_color_post$population_logjnd_intercept_mean + ex_toj_color_post$population_logjnd_initial_bias_effect_mean/2 )
      - exp( ex_toj_color_post$population_logjnd_intercept_mean - ex_toj_color_post$population_logjnd_initial_bias_effect_mean/2  ) ) * 250 
  )
  , c("PSS Initial Probe\nBias Effect Mean"  , "JND Initial Probe\nBias Effect Mean")
  , y_lab = "SOA (Left - Right; ms)"
)

# effect of probe duration on PSS and JND
get_violin(
  c(
  ( ex_toj_color_post$population_pss_probe_effect_mean ) * 250
  , ( exp( ex_toj_color_post$population_logjnd_intercept_mean + ex_toj_color_post$population_logjnd_probe_effect_mean/2 )
      - exp( ex_toj_color_post$population_logjnd_intercept_mean - ex_toj_color_post$population_logjnd_probe_effect_mean/2  ) ) * 250 
  )
  , c("PSS Probe Duration\nBias Effect Mean" , "JND Probe Duration\nBias Effect Mean")
  , y_lab = "SOA (Long - Short; ms)"
)
#---------------------------------- Main Effects ---------------------=--------------------#


#------------------------------- Two-way Interactions -------------------------------------#
#  effect of interaction between judgement type and attention on PSS 
get_violin(
  c(
  (ex_toj_color_post$population_pss_judgement_type_interaction_effect_mean) * 250
  ,   ( exp( ex_toj_color_post$population_logjnd_intercept_mean + ex_toj_color_post$population_logjnd_judgement_type_interaction_effect_mean/2 )
        - exp( ex_toj_color_post$population_logjnd_intercept_mean - ex_toj_color_post$population_logjnd_judgement_type_interaction_effect_mean/2  ) ) * 250 
  )
  , c("PSS Attention\n& Judgement Type\nInteraction Effect Mean"  , "JND Attention\n& Judgement Type\nInteraction Effect Mean")
  , y_lab = "SOA (ms)"
)

#  effect of interaction between initial bias and attention on PSS 
get_violin(
  c(
  (ex_toj_color_post$population_pss_initial_bias_interaction_effect_mean) * 250
  ,  ( exp( ex_toj_color_post$population_logjnd_intercept_mean + ex_toj_color_post$population_logjnd_initial_bias_interaction_effect_mean/2 )
       - exp( ex_toj_color_post$population_logjnd_intercept_mean - ex_toj_color_post$population_logjnd_initial_bias_interaction_effect_mean/2  ) ) * 250 
  )
  , c("PSS Attention\n& Initial Probe Bias\nInteraction Effect Mean", "JND Attention\n& Initial Probe Bias\nInteraction Effect Mean")
  , y_lab = "SOA (ms)"
)

# effect of interaction between probe duration and attention on PSS
get_violin(
  c(
  (ex_toj_color_post$population_pss_probe_interaction_effect_mean) * 250
  ,  ( exp( ex_toj_color_post$population_logjnd_intercept_mean + ex_toj_color_post$population_logjnd_probe_interaction_effect_mean/2 )
       - exp( ex_toj_color_post$population_logjnd_intercept_mean - ex_toj_color_post$population_logjnd_probe_interaction_effect_mean/2  ) ) * 250 
  )
  , c("PSS Attention\n& Probe Duration\nInteraction Effect Mean"  , "JND Attention\n& Probe Duration\nInteraction Effect Mean")
  , y_lab = "SOA (ms)"
)
#------------------------------- Two-way Interactions -------------------------------------#


#------------------------------- Three-way Interactions -----------------------------------#
#  effect of interaction between judgement type and initial bias and attention on PSS 
get_violin(
  c(
    (ex_toj_color_post$population_pss_initial_bias_judgement_type_interaction_effect_mean) * 250
    ,   ( exp( ex_toj_color_post$population_logjnd_intercept_mean + ex_toj_color_post$population_logjnd_initial_bias_judgement_type_interaction_effect_mean/2 )
          - exp( ex_toj_color_post$population_logjnd_intercept_mean - ex_toj_color_post$population_logjnd_initial_bias_judgement_type_interaction_effect_mean/2  ) ) * 250 
  )
  , c("PSS Attention\n& Judgement Type\n& Initial Bias\nInteraction Effect Mean"  , "JND Attention\n& Judgement Type\n& Initial Bias\nInteraction Effect Mean")
  , y_lab = "SOA (ms)"
)

#  effect of interaction between initial bias and probe duration and attention on PSS 
get_violin(
  c(
    (ex_toj_color_post$population_pss_initial_bias_probe_interaction_effect_mean) * 250
    ,   ( exp( ex_toj_color_post$population_logjnd_intercept_mean + ex_toj_color_post$population_logjnd_initial_bias_probe_interaction_effect_mean/2 )
          - exp( ex_toj_color_post$population_logjnd_intercept_mean - ex_toj_color_post$population_logjnd_initial_bias_probe_interaction_effect_mean/2  ) ) * 250 
  )
  , c("PSS Attention\n& Initial Bias\n& Probe Duration\nInteraction Effect Mean"  , "JND Attention\n& Initial Bias\n& Probe Duration\nInteraction Effect Mean")
  , y_lab = "SOA (ms)"
)

#  effect of interaction between judgement type and probe duration and attention on PSS 
get_violin(
  c(
    (ex_toj_color_post$population_pss_judgement_type_probe_interaction_effect_mean) * 250
    ,   ( exp( ex_toj_color_post$population_logjnd_intercept_mean + ex_toj_color_post$population_logjnd_judgement_type_probe_interaction_effect_mean/2 )
          - exp( ex_toj_color_post$population_logjnd_intercept_mean - ex_toj_color_post$population_logjnd_judgement_type_probe_interaction_effect_mean/2  ) ) * 250 
  )
  , c("PSS Attention\n& Judgement Type\n& Probe Duration\nInteraction Effect Mean"  , "JND Attention\n& Judgement Type\n& Probe Duration\nInteraction Effect Mean")
  , y_lab = "SOA (ms)"
)
#------------------------------- Three-way Interactions -----------------------------------#


#------------------------------- Four-way Interactions -----------------------------------#
get_violin(
  c(
    (ex_toj_color_post$population_pss_four_way_interaction_effect_mean) * 250
    ,   ( exp( ex_toj_color_post$population_logjnd_intercept_mean + ex_toj_color_post$population_logjnd_four_way_interaction_effect_mean/2 )
          - exp( ex_toj_color_post$population_logjnd_intercept_mean - ex_toj_color_post$population_logjnd_four_way_interaction_effect_mean/2  ) ) * 250 
  )
  , c("PSS Four-Way\nInteraction Effect Mean"  , "JND Four-Way\nInteraction Effect Mean")
  , y_lab = "SOA (ms)"
)
#------------------------------- Four-way Interactions -----------------------------------#



#------------------------------------------------------------------------------------------#
#--------------------------------- Rho Scale ----------------------------------------------#
#------------------------------------------------------------------------------------------#

#---------------------------------- Rho Intercept -----------------------------------------#
get_violin(
  plogis(ex_toj_color_post$logitRhoMean)
  , "Probability of Memory Intercept Mean"
  , y_lab = "\u03C1"
  , hline = FALSE
)
#---------------------------------- Rho Intercept -----------------------------------------#


#-------------------------------- Main Effects --------------------------------------------#
get_violin(
  ( plogis(ex_toj_color_post$logitRhoMean + ex_toj_color_post$logitRhoEffectMean/2 )
    - plogis(ex_toj_color_post$logitRhoMean - ex_toj_color_post$logitRhoEffectMean/2 ) )
  , "Probability of Memory Effect Mean"
  , y_lab = "\u03C1 (Attended - Unattended)"
)

get_violin(
  ( plogis(ex_toj_color_post$logitRhoMean + ex_toj_color_post$logitRhoJudgementTypeEffectMean/2 )
    - plogis(ex_toj_color_post$logitRhoMean - ex_toj_color_post$logitRhoJudgementTypeEffectMean/2 ) )
  , "Probability of Memory\nJudgement Type Effect Mean"
  , y_lab = "\u03C1 (Second - First)"
)

get_violin(
  ( plogis(ex_toj_color_post$logitRhoMean + ex_toj_color_post$logitRhoInitialBiasEffectMean/2 )
    - plogis(ex_toj_color_post$logitRhoMean - ex_toj_color_post$logitRhoInitialBiasEffectMean/2 ) )
  , "Probability of Memory\nInitial Probe Bias Effect Mean"
  , y_lab = "\u03C1 (Left - Right)"
)

get_violin(
  ( plogis(ex_toj_color_post$logitRhoMean + ex_toj_color_post$logitRhoProbeEffectMean/2 )
    - plogis(ex_toj_color_post$logitRhoMean - ex_toj_color_post$logitRhoProbeEffectMean/2 ) )
  , "Probability of Memory\nProbe Duration Effect Mean"
  , y_lab = "\u03C1 (Long - Short)"
)
#-------------------------------- Main Effects --------------------------------------------#


#------------------------------- Two-way Interactions -------------------------------------#
get_violin(
  ( plogis(ex_toj_color_post$logitRhoMean + ex_toj_color_post$logitRhoJudgementTypeInteractionEffectMean/2 )
    - plogis(ex_toj_color_post$logitRhoMean - ex_toj_color_post$logitRhoJudgementTypeInteractionEffectMean/2 ) )
  , "Probability of Memory\nAttention\n& Judgement Type\nInteraction Effect Mean"
  , y_lab = "\u03C1"
)

get_violin(
  ( plogis(ex_toj_color_post$logitRhoMean + ex_toj_color_post$logitRhoInitialBiasInteractionEffectMean/2 )
    - plogis(ex_toj_color_post$logitRhoMean - ex_toj_color_post$logitRhoInitialBiasInteractionEffectMean/2 ) )
  , "Probability of Memory\nAttention\nInitial Probe Bias\nInteraction Effect Mean"
  , y_lab = "\u03C1"
)

get_violin(
  ( plogis(ex_toj_color_post$logitRhoMean + ex_toj_color_post$logitRhoProbeInteractionEffectMean/2 )
    - plogis(ex_toj_color_post$logitRhoMean - ex_toj_color_post$logitRhoProbeInteractionEffectMean/2 ) )
  , "Probability of Memory\nAttention\nProbe Duration Interaction Effect Mean"
  , y_lab = "\u03C1"
)
#------------------------------- Two-way Interactions -------------------------------------#


#-------------------------------- Three-way Interactions ----------------------------------#
get_violin(
  ( plogis(ex_toj_color_post$logitRhoMean + ex_toj_color_post$logitRhoInitialBiasJudgementTypeInteractionEffectMean/2 )
    - plogis(ex_toj_color_post$logitRhoMean - ex_toj_color_post$logitRhoInitialBiasJudgementTypeInteractionEffectMean/2 ) )
  , "Probability of Memory\nAttention\n& Initial Bias\n& Judgement Type\nInteraction Effect Mean"
  , y_lab = "\u03C1"
)

get_violin(
  ( plogis(ex_toj_color_post$logitRhoMean + ex_toj_color_post$logitRhoInitialBiasProbeInteractionEffectMean/2 )
    - plogis(ex_toj_color_post$logitRhoMean - ex_toj_color_post$logitRhoInitialBiasProbeInteractionEffectMean/2 ) )
  , "Probability of Memory\nAttention\n& Initial Bias\n& Probe Duration\nInteraction Effect Mean"
  , y_lab = "\u03C1"
)

get_violin(
  ( plogis(ex_toj_color_post$logitRhoMean + ex_toj_color_post$logitRhoJudgementTypeProbeInteractionEffectMean/2 )
    - plogis(ex_toj_color_post$logitRhoMean - ex_toj_color_post$logitRhoJudgementTypeProbeInteractionEffectMean/2 ) )
  , "Probability of Memory\nAttention\n& Probe Duration\n& Judgement Type\nInteraction Effect Mean"
  , y_lab = "\u03C1"
)
#-------------------------------- Three-way Interactions ----------------------------------#


#-------------------------------- Four-way Interactions -----------------------------------#
get_violin(
  ( plogis(ex_toj_color_post$logitRhoMean + ex_toj_color_post$logitRhoFourWayInteractionEffectMean/2 )
    - plogis(ex_toj_color_post$logitRhoMean - ex_toj_color_post$logitRhoFourWayInteractionEffectMean/2 ) )
  , "Probability of Memory\nFour-Way\nInteraction Effect Mean"
  , y_lab = "\u03C1"
)

get_violin(
  c(
    plogis(ex_toj_color_post$logitRhoMean 
           + ex_toj_color_post$logitRhoInitialBiasEffectMean/2
           + ex_toj_color_post$logitRhoJudgementTypeEffectMean/2
           + ex_toj_color_post$logitRhoProbeEffectMean/2
           + (
             ex_toj_color_post$logitRhoInitialBiasInteractionEffectMean
             + ex_toj_color_post$logitRhoJudgementTypeInteractionEffectMean
             + ex_toj_color_post$logitRhoProbeInteractionEffectMean
             
             + ex_toj_color_post$logitRhoInitialBiasJudgementTypeInteractionEffectMean
             + ex_toj_color_post$logitRhoInitialBiasProbeInteractionEffectMean
             + ex_toj_color_post$logitRhoJudgementTypeProbeInteractionEffectMean
             + ex_toj_color_post$logitRhoFourWayInteractionEffectMean
           )/2
    )
    -     plogis(ex_toj_color_post$logitRhoMean 
                 + ex_toj_color_post$logitRhoInitialBiasEffectMean/2
                 + ex_toj_color_post$logitRhoJudgementTypeEffectMean/2
                 + ex_toj_color_post$logitRhoProbeEffectMean/2
                 - (
                   ex_toj_color_post$logitRhoInitialBiasInteractionEffectMean
                   + ex_toj_color_post$logitRhoJudgementTypeInteractionEffectMean
                   + ex_toj_color_post$logitRhoProbeInteractionEffectMean
                   
                   + ex_toj_color_post$logitRhoInitialBiasJudgementTypeInteractionEffectMean
                   + ex_toj_color_post$logitRhoInitialBiasProbeInteractionEffectMean
                   + ex_toj_color_post$logitRhoJudgementTypeProbeInteractionEffectMean
                   + ex_toj_color_post$logitRhoFourWayInteractionEffectMean
                 )/2
    )
    ,     plogis(ex_toj_color_post$logitRhoMean 
                 + ex_toj_color_post$logitRhoInitialBiasEffectMean/2
                 + ex_toj_color_post$logitRhoJudgementTypeEffectMean/2
                 - ex_toj_color_post$logitRhoProbeEffectMean/2
                 + (
                   ex_toj_color_post$logitRhoInitialBiasInteractionEffectMean
                   + ex_toj_color_post$logitRhoJudgementTypeInteractionEffectMean
                   - ex_toj_color_post$logitRhoProbeInteractionEffectMean
                   
                   + ex_toj_color_post$logitRhoInitialBiasJudgementTypeInteractionEffectMean
                   - ex_toj_color_post$logitRhoInitialBiasProbeInteractionEffectMean
                   - ex_toj_color_post$logitRhoJudgementTypeProbeInteractionEffectMean
                   - ex_toj_color_post$logitRhoFourWayInteractionEffectMean
                 )/2
    )
    -     plogis(ex_toj_color_post$logitRhoMean 
                 + ex_toj_color_post$logitRhoInitialBiasEffectMean/2
                 + ex_toj_color_post$logitRhoJudgementTypeEffectMean/2
                 - ex_toj_color_post$logitRhoProbeEffectMean/2
                 - (
                   ex_toj_color_post$logitRhoInitialBiasInteractionEffectMean
                   + ex_toj_color_post$logitRhoJudgementTypeInteractionEffectMean
                   - ex_toj_color_post$logitRhoProbeInteractionEffectMean
                   
                   + ex_toj_color_post$logitRhoInitialBiasJudgementTypeInteractionEffectMean
                   - ex_toj_color_post$logitRhoInitialBiasProbeInteractionEffectMean
                   - ex_toj_color_post$logitRhoJudgementTypeProbeInteractionEffectMean
                   - ex_toj_color_post$logitRhoFourWayInteractionEffectMean
                 )/2
    )
    ,     plogis(ex_toj_color_post$logitRhoMean 
                 + ex_toj_color_post$logitRhoInitialBiasEffectMean/2
                 - ex_toj_color_post$logitRhoJudgementTypeEffectMean/2
                 - ex_toj_color_post$logitRhoProbeEffectMean/2
                 + (
                   ex_toj_color_post$logitRhoInitialBiasInteractionEffectMean
                   - ex_toj_color_post$logitRhoJudgementTypeInteractionEffectMean
                   - ex_toj_color_post$logitRhoProbeInteractionEffectMean
                   
                   - ex_toj_color_post$logitRhoInitialBiasJudgementTypeInteractionEffectMean
                   - ex_toj_color_post$logitRhoInitialBiasProbeInteractionEffectMean
                   + ex_toj_color_post$logitRhoJudgementTypeProbeInteractionEffectMean
                   + ex_toj_color_post$logitRhoFourWayInteractionEffectMean
                 )/2
    )
    -     plogis(ex_toj_color_post$logitRhoMean 
                 + ex_toj_color_post$logitRhoInitialBiasEffectMean/2
                 - ex_toj_color_post$logitRhoJudgementTypeEffectMean/2
                 - ex_toj_color_post$logitRhoProbeEffectMean/2
                 - (
                   ex_toj_color_post$logitRhoInitialBiasInteractionEffectMean
                   - ex_toj_color_post$logitRhoJudgementTypeInteractionEffectMean
                   - ex_toj_color_post$logitRhoProbeInteractionEffectMean
                   
                   - ex_toj_color_post$logitRhoInitialBiasJudgementTypeInteractionEffectMean
                   - ex_toj_color_post$logitRhoInitialBiasProbeInteractionEffectMean
                   + ex_toj_color_post$logitRhoJudgementTypeProbeInteractionEffectMean
                   + ex_toj_color_post$logitRhoFourWayInteractionEffectMean
                 )/2
    )
    ,     plogis(ex_toj_color_post$logitRhoMean 
                 + ex_toj_color_post$logitRhoInitialBiasEffectMean/2
                 - ex_toj_color_post$logitRhoJudgementTypeEffectMean/2
                 + ex_toj_color_post$logitRhoProbeEffectMean/2
                 + (
                   ex_toj_color_post$logitRhoInitialBiasInteractionEffectMean
                   - ex_toj_color_post$logitRhoJudgementTypeInteractionEffectMean
                   + ex_toj_color_post$logitRhoProbeInteractionEffectMean
                   
                   - ex_toj_color_post$logitRhoInitialBiasJudgementTypeInteractionEffectMean
                   + ex_toj_color_post$logitRhoInitialBiasProbeInteractionEffectMean
                   - ex_toj_color_post$logitRhoJudgementTypeProbeInteractionEffectMean
                   - ex_toj_color_post$logitRhoFourWayInteractionEffectMean
                 )/2
    )
    -     plogis(ex_toj_color_post$logitRhoMean 
                 + ex_toj_color_post$logitRhoInitialBiasEffectMean/2
                 - ex_toj_color_post$logitRhoJudgementTypeEffectMean/2
                 + ex_toj_color_post$logitRhoProbeEffectMean/2
                 - (
                   ex_toj_color_post$logitRhoInitialBiasInteractionEffectMean
                   - ex_toj_color_post$logitRhoJudgementTypeInteractionEffectMean
                   + ex_toj_color_post$logitRhoProbeInteractionEffectMean
                   
                   - ex_toj_color_post$logitRhoInitialBiasJudgementTypeInteractionEffectMean
                   + ex_toj_color_post$logitRhoInitialBiasProbeInteractionEffectMean
                   - ex_toj_color_post$logitRhoJudgementTypeProbeInteractionEffectMean
                   - ex_toj_color_post$logitRhoFourWayInteractionEffectMean
                 )/2
    )
    
    
    
    , plogis(ex_toj_color_post$logitRhoMean 
           - ex_toj_color_post$logitRhoInitialBiasEffectMean/2
           + ex_toj_color_post$logitRhoJudgementTypeEffectMean/2
           + ex_toj_color_post$logitRhoProbeEffectMean/2
           + (
             -ex_toj_color_post$logitRhoInitialBiasInteractionEffectMean
             + ex_toj_color_post$logitRhoJudgementTypeInteractionEffectMean
             + ex_toj_color_post$logitRhoProbeInteractionEffectMean
             
             - ex_toj_color_post$logitRhoInitialBiasJudgementTypeInteractionEffectMean
             - ex_toj_color_post$logitRhoInitialBiasProbeInteractionEffectMean
             + ex_toj_color_post$logitRhoJudgementTypeProbeInteractionEffectMean
             - ex_toj_color_post$logitRhoFourWayInteractionEffectMean
           )/2
    )
    -     plogis(ex_toj_color_post$logitRhoMean 
                 - ex_toj_color_post$logitRhoInitialBiasEffectMean/2
                 + ex_toj_color_post$logitRhoJudgementTypeEffectMean/2
                 + ex_toj_color_post$logitRhoProbeEffectMean/2
                 - (
                   -ex_toj_color_post$logitRhoInitialBiasInteractionEffectMean
                   + ex_toj_color_post$logitRhoJudgementTypeInteractionEffectMean
                   + ex_toj_color_post$logitRhoProbeInteractionEffectMean
                   
                   - ex_toj_color_post$logitRhoInitialBiasJudgementTypeInteractionEffectMean
                   - ex_toj_color_post$logitRhoInitialBiasProbeInteractionEffectMean
                   + ex_toj_color_post$logitRhoJudgementTypeProbeInteractionEffectMean
                   - ex_toj_color_post$logitRhoFourWayInteractionEffectMean
                 )/2
    )
    ,     plogis(ex_toj_color_post$logitRhoMean 
                 - ex_toj_color_post$logitRhoInitialBiasEffectMean/2
                 + ex_toj_color_post$logitRhoJudgementTypeEffectMean/2
                 - ex_toj_color_post$logitRhoProbeEffectMean/2
                 + (
                   - ex_toj_color_post$logitRhoInitialBiasInteractionEffectMean
                   + ex_toj_color_post$logitRhoJudgementTypeInteractionEffectMean
                   - ex_toj_color_post$logitRhoProbeInteractionEffectMean
                   
                   - ex_toj_color_post$logitRhoInitialBiasJudgementTypeInteractionEffectMean
                   + ex_toj_color_post$logitRhoInitialBiasProbeInteractionEffectMean
                   - ex_toj_color_post$logitRhoJudgementTypeProbeInteractionEffectMean
                   + ex_toj_color_post$logitRhoFourWayInteractionEffectMean
                 )/2
    )
    -     plogis(ex_toj_color_post$logitRhoMean 
                 - ex_toj_color_post$logitRhoInitialBiasEffectMean/2
                 + ex_toj_color_post$logitRhoJudgementTypeEffectMean/2
                 - ex_toj_color_post$logitRhoProbeEffectMean/2
                 - (
                   - ex_toj_color_post$logitRhoInitialBiasInteractionEffectMean
                   + ex_toj_color_post$logitRhoJudgementTypeInteractionEffectMean
                   - ex_toj_color_post$logitRhoProbeInteractionEffectMean
                   
                   - ex_toj_color_post$logitRhoInitialBiasJudgementTypeInteractionEffectMean
                   + ex_toj_color_post$logitRhoInitialBiasProbeInteractionEffectMean
                   - ex_toj_color_post$logitRhoJudgementTypeProbeInteractionEffectMean
                   + ex_toj_color_post$logitRhoFourWayInteractionEffectMean
                 )/2
    )
    ,     plogis(ex_toj_color_post$logitRhoMean 
                 - ex_toj_color_post$logitRhoInitialBiasEffectMean/2
                 - ex_toj_color_post$logitRhoJudgementTypeEffectMean/2
                 - ex_toj_color_post$logitRhoProbeEffectMean/2
                 + (
                   - ex_toj_color_post$logitRhoInitialBiasInteractionEffectMean
                   - ex_toj_color_post$logitRhoJudgementTypeInteractionEffectMean
                   - ex_toj_color_post$logitRhoProbeInteractionEffectMean
                   
                   + ex_toj_color_post$logitRhoInitialBiasJudgementTypeInteractionEffectMean
                   + ex_toj_color_post$logitRhoInitialBiasProbeInteractionEffectMean
                   + ex_toj_color_post$logitRhoJudgementTypeProbeInteractionEffectMean
                   - ex_toj_color_post$logitRhoFourWayInteractionEffectMean
                 )/2
    )
    -     plogis(ex_toj_color_post$logitRhoMean 
                 - ex_toj_color_post$logitRhoInitialBiasEffectMean/2
                 - ex_toj_color_post$logitRhoJudgementTypeEffectMean/2
                 - ex_toj_color_post$logitRhoProbeEffectMean/2
                 - (
                   - ex_toj_color_post$logitRhoInitialBiasInteractionEffectMean
                   - ex_toj_color_post$logitRhoJudgementTypeInteractionEffectMean
                   - ex_toj_color_post$logitRhoProbeInteractionEffectMean
                   
                   + ex_toj_color_post$logitRhoInitialBiasJudgementTypeInteractionEffectMean
                   + ex_toj_color_post$logitRhoInitialBiasProbeInteractionEffectMean
                   + ex_toj_color_post$logitRhoJudgementTypeProbeInteractionEffectMean
                   - ex_toj_color_post$logitRhoFourWayInteractionEffectMean
                 )/2
    )
    ,     plogis(ex_toj_color_post$logitRhoMean 
                 - ex_toj_color_post$logitRhoInitialBiasEffectMean/2
                 - ex_toj_color_post$logitRhoJudgementTypeEffectMean/2
                 + ex_toj_color_post$logitRhoProbeEffectMean/2
                 + (
                   - ex_toj_color_post$logitRhoInitialBiasInteractionEffectMean
                   - ex_toj_color_post$logitRhoJudgementTypeInteractionEffectMean
                   + ex_toj_color_post$logitRhoProbeInteractionEffectMean
                   
                   + ex_toj_color_post$logitRhoInitialBiasJudgementTypeInteractionEffectMean
                   - ex_toj_color_post$logitRhoInitialBiasProbeInteractionEffectMean
                   - ex_toj_color_post$logitRhoJudgementTypeProbeInteractionEffectMean
                   + ex_toj_color_post$logitRhoFourWayInteractionEffectMean
                 )/2
    )
    -     plogis(ex_toj_color_post$logitRhoMean 
                 - ex_toj_color_post$logitRhoInitialBiasEffectMean/2
                 - ex_toj_color_post$logitRhoJudgementTypeEffectMean/2
                 + ex_toj_color_post$logitRhoProbeEffectMean/2
                 - (
                   - ex_toj_color_post$logitRhoInitialBiasInteractionEffectMean
                   - ex_toj_color_post$logitRhoJudgementTypeInteractionEffectMean
                   + ex_toj_color_post$logitRhoProbeInteractionEffectMean
                   
                   + ex_toj_color_post$logitRhoInitialBiasJudgementTypeInteractionEffectMean
                   - ex_toj_color_post$logitRhoInitialBiasProbeInteractionEffectMean
                   - ex_toj_color_post$logitRhoJudgementTypeProbeInteractionEffectMean
                   + ex_toj_color_post$logitRhoFourWayInteractionEffectMean
                 )/2
    )  
    
  )
  , c("Left\nSecond\nLong", "Left\nSecond\nShort", "Left\nFirst\nShort", "Left\nFirst\nLong"
      , "Right\nSecond\nLong", "Right\nSecond\nShort", "Right\nFirst\nShort", "Right\nFirst\nLong")
  , y_lab = "\u03C1 (Attended - Unattended)"
)
#-------------------------------- Four-way Interactions -----------------------------------#



#------------------------------------------------------------------------------------------#
#--------------------------------- Kappa Scale --------------------------------------------#
#------------------------------------------------------------------------------------------#
 
#---------------------------------- Kappa Intercept ---------------------------------------#
get_violin(
  exp(ex_toj_color_post$logKappaMean)
  , "Fidelity of Memory Intercept Mean"
  , y_lab = "\u03BA"
  , hline = FALSE
)
#---------------------------------- Kappa Intercept ---------------------------------------#


#-------------------------------- Main Effects --------------------------------------------#
get_violin(
  ( exp(ex_toj_color_post$logKappaMean + ex_toj_color_post$logKappaEffectMean/2 )
    - exp(ex_toj_color_post$logKappaMean - ex_toj_color_post$logKappaEffectMean/2 ) )
  , "Fidelity of Memory Effect Mean"
  , y_lab = "\u03BA (Attended - Unattended)"
)

get_violin(
  ( exp(ex_toj_color_post$logKappaMean + ex_toj_color_post$logKappaJudgementTypeEffectMean/2 )
    - exp(ex_toj_color_post$logKappaMean - ex_toj_color_post$logKappaJudgementTypeEffectMean/2 ) )
  , "Fidelity of Memory\nJudgement Type Effect Mean"
  , y_lab = "\u03BA (Second - First)"
)

get_violin(
  ( exp(ex_toj_color_post$logKappaMean + ex_toj_color_post$logKappaInitialBiasEffectMean/2 )
    - exp(ex_toj_color_post$logKappaMean - ex_toj_color_post$logKappaInitialBiasEffectMean/2 ) )
  , "Fidelity of Memory\nInitial Probe Bias Effect Mean"
  , y_lab = "\u03BA (Left - Right)"
)

get_violin(
  ( exp(ex_toj_color_post$logKappaMean + ex_toj_color_post$logKappaProbeEffectMean/2 )
    - exp(ex_toj_color_post$logKappaMean - ex_toj_color_post$logKappaProbeEffectMean/2 ) )
  , "Fidelity of Memory\nProbe Duration Effect Mean"
  , y_lab = "\u03BA (Long - Short)"
)
#-------------------------------- Main Effects --------------------------------------------#


#------------------------------- Two-way Interactions -------------------------------------#
get_violin(
  ( exp(ex_toj_color_post$logKappaMean + ex_toj_color_post$logKappaJudgementTypeInteractionEffectMean/2 )
    - exp(ex_toj_color_post$logKappaMean - ex_toj_color_post$logKappaJudgementTypeInteractionEffectMean/2 ) )
  , "Fidelity of Memory\nAttention\n& Judgement Type\nInteraction Effect Mean"
  , y_lab = "\u03BA"
)

get_violin(
  ( exp(ex_toj_color_post$logKappaMean + ex_toj_color_post$logKappaInitialBiasInteractionEffectMean/2 )
    - exp(ex_toj_color_post$logKappaMean - ex_toj_color_post$logKappaInitialBiasInteractionEffectMean/2 ) )
  , "Fidelity of Memory\nAttention\nInitial Probe Bias\nInteraction Effect Mean"
  , y_lab = "\u03BA"
)

get_violin(
  ( exp(ex_toj_color_post$logKappaMean + ex_toj_color_post$logKappaProbeInteractionEffectMean/2 )
    - exp(ex_toj_color_post$logKappaMean - ex_toj_color_post$logKappaProbeInteractionEffectMean/2 ) )
  , "Fidelity of Memory\nAttention\nProbe Duration Interaction Effect Mean"
  , y_lab = "\u03BA"
)
#------------------------------- Two-way Interactions -------------------------------------#


#-------------------------------- Three-way Interactions ----------------------------------#
get_violin(
  ( exp(ex_toj_color_post$logKappaMean + ex_toj_color_post$logKappaInitialBiasJudgementTypeInteractionEffectMean/2 )
    - exp(ex_toj_color_post$logKappaMean - ex_toj_color_post$logKappaInitialBiasJudgementTypeInteractionEffectMean/2 ) )
  , "Fidelity of Memory\nAttention\n& Initial Bias\n& Judgement Type\nInteraction Effect Mean"
  , y_lab = "\u03BA"
)

get_violin(
  ( exp(ex_toj_color_post$logKappaMean + ex_toj_color_post$logKappaInitialBiasProbeInteractionEffectMean/2 )
    - exp(ex_toj_color_post$logKappaMean - ex_toj_color_post$logKappaInitialBiasProbeInteractionEffectMean/2 ) )
  , "Fidelity of Memory\nAttention\n& Initial Bias\n& Probe Duration\nInteraction Effect Mean"
  , y_lab = "\u03BA"
)

get_violin(
  ( exp(ex_toj_color_post$logKappaMean + ex_toj_color_post$logKappaJudgementTypeProbeInteractionEffectMean/2 )
    - exp(ex_toj_color_post$logKappaMean - ex_toj_color_post$logKappaJudgementTypeProbeInteractionEffectMean/2 ) )
  , "Fidelity of Memory\nAttention\n& Probe Duration\n& Judgement Type\nInteraction Effect Mean"
  , y_lab = "\u03BA"
)
#-------------------------------- Three-way Interactions ----------------------------------#


#-------------------------------- Four-way Interactions -----------------------------------#
get_violin(
  ( exp(ex_toj_color_post$logKappaMean + ex_toj_color_post$logKappaFourWayInteractionEffectMean/2 )
    - exp(ex_toj_color_post$logKappaMean - ex_toj_color_post$logKappaFourWayInteractionEffectMean/2 ) )
  , "Fidelity of Memory\nFour-Way\nInteraction Effect Mean"
  , y_lab = "\u03BA"
)

get_violin(
  c(
    exp(ex_toj_color_post$logKappaMean 
           + ex_toj_color_post$logKappaInitialBiasEffectMean/2
           + ex_toj_color_post$logKappaJudgementTypeEffectMean/2
           + ex_toj_color_post$logKappaProbeEffectMean/2
           + (
             ex_toj_color_post$logKappaInitialBiasInteractionEffectMean
             + ex_toj_color_post$logKappaJudgementTypeInteractionEffectMean
             + ex_toj_color_post$logKappaProbeInteractionEffectMean
             
             + ex_toj_color_post$logKappaInitialBiasJudgementTypeInteractionEffectMean
             + ex_toj_color_post$logKappaInitialBiasProbeInteractionEffectMean
             + ex_toj_color_post$logKappaJudgementTypeProbeInteractionEffectMean
             + ex_toj_color_post$logKappaFourWayInteractionEffectMean
           )/2
    )
    -     exp(ex_toj_color_post$logKappaMean 
                 + ex_toj_color_post$logKappaInitialBiasEffectMean/2
                 + ex_toj_color_post$logKappaJudgementTypeEffectMean/2
                 + ex_toj_color_post$logKappaProbeEffectMean/2
                 - (
                   ex_toj_color_post$logKappaInitialBiasInteractionEffectMean
                   + ex_toj_color_post$logKappaJudgementTypeInteractionEffectMean
                   + ex_toj_color_post$logKappaProbeInteractionEffectMean
                   
                   + ex_toj_color_post$logKappaInitialBiasJudgementTypeInteractionEffectMean
                   + ex_toj_color_post$logKappaInitialBiasProbeInteractionEffectMean
                   + ex_toj_color_post$logKappaJudgementTypeProbeInteractionEffectMean
                   + ex_toj_color_post$logKappaFourWayInteractionEffectMean
                 )/2
    )
    ,     exp(ex_toj_color_post$logKappaMean 
                 + ex_toj_color_post$logKappaInitialBiasEffectMean/2
                 + ex_toj_color_post$logKappaJudgementTypeEffectMean/2
                 - ex_toj_color_post$logKappaProbeEffectMean/2
                 + (
                   ex_toj_color_post$logKappaInitialBiasInteractionEffectMean
                   + ex_toj_color_post$logKappaJudgementTypeInteractionEffectMean
                   - ex_toj_color_post$logKappaProbeInteractionEffectMean
                   
                   + ex_toj_color_post$logKappaInitialBiasJudgementTypeInteractionEffectMean
                   - ex_toj_color_post$logKappaInitialBiasProbeInteractionEffectMean
                   - ex_toj_color_post$logKappaJudgementTypeProbeInteractionEffectMean
                   - ex_toj_color_post$logKappaFourWayInteractionEffectMean
                 )/2
    )
    -     exp(ex_toj_color_post$logKappaMean 
                 + ex_toj_color_post$logKappaInitialBiasEffectMean/2
                 + ex_toj_color_post$logKappaJudgementTypeEffectMean/2
                 - ex_toj_color_post$logKappaProbeEffectMean/2
                 - (
                   ex_toj_color_post$logKappaInitialBiasInteractionEffectMean
                   + ex_toj_color_post$logKappaJudgementTypeInteractionEffectMean
                   - ex_toj_color_post$logKappaProbeInteractionEffectMean
                   
                   + ex_toj_color_post$logKappaInitialBiasJudgementTypeInteractionEffectMean
                   - ex_toj_color_post$logKappaInitialBiasProbeInteractionEffectMean
                   - ex_toj_color_post$logKappaJudgementTypeProbeInteractionEffectMean
                   - ex_toj_color_post$logKappaFourWayInteractionEffectMean
                 )/2
    )
    ,     exp(ex_toj_color_post$logKappaMean 
                 + ex_toj_color_post$logKappaInitialBiasEffectMean/2
                 - ex_toj_color_post$logKappaJudgementTypeEffectMean/2
                 - ex_toj_color_post$logKappaProbeEffectMean/2
                 + (
                   ex_toj_color_post$logKappaInitialBiasInteractionEffectMean
                   - ex_toj_color_post$logKappaJudgementTypeInteractionEffectMean
                   - ex_toj_color_post$logKappaProbeInteractionEffectMean
                   
                   - ex_toj_color_post$logKappaInitialBiasJudgementTypeInteractionEffectMean
                   - ex_toj_color_post$logKappaInitialBiasProbeInteractionEffectMean
                   + ex_toj_color_post$logKappaJudgementTypeProbeInteractionEffectMean
                   + ex_toj_color_post$logKappaFourWayInteractionEffectMean
                 )/2
    )
    -     exp(ex_toj_color_post$logKappaMean 
                 + ex_toj_color_post$logKappaInitialBiasEffectMean/2
                 - ex_toj_color_post$logKappaJudgementTypeEffectMean/2
                 - ex_toj_color_post$logKappaProbeEffectMean/2
                 - (
                   ex_toj_color_post$logKappaInitialBiasInteractionEffectMean
                   - ex_toj_color_post$logKappaJudgementTypeInteractionEffectMean
                   - ex_toj_color_post$logKappaProbeInteractionEffectMean
                   
                   - ex_toj_color_post$logKappaInitialBiasJudgementTypeInteractionEffectMean
                   - ex_toj_color_post$logKappaInitialBiasProbeInteractionEffectMean
                   + ex_toj_color_post$logKappaJudgementTypeProbeInteractionEffectMean
                   + ex_toj_color_post$logKappaFourWayInteractionEffectMean
                 )/2
    )
    ,     exp(ex_toj_color_post$logKappaMean 
                 + ex_toj_color_post$logKappaInitialBiasEffectMean/2
                 - ex_toj_color_post$logKappaJudgementTypeEffectMean/2
                 + ex_toj_color_post$logKappaProbeEffectMean/2
                 + (
                   ex_toj_color_post$logKappaInitialBiasInteractionEffectMean
                   - ex_toj_color_post$logKappaJudgementTypeInteractionEffectMean
                   + ex_toj_color_post$logKappaProbeInteractionEffectMean
                   
                   - ex_toj_color_post$logKappaInitialBiasJudgementTypeInteractionEffectMean
                   + ex_toj_color_post$logKappaInitialBiasProbeInteractionEffectMean
                   - ex_toj_color_post$logKappaJudgementTypeProbeInteractionEffectMean
                   - ex_toj_color_post$logKappaFourWayInteractionEffectMean
                 )/2
    )
    -     exp(ex_toj_color_post$logKappaMean 
                 + ex_toj_color_post$logKappaInitialBiasEffectMean/2
                 - ex_toj_color_post$logKappaJudgementTypeEffectMean/2
                 + ex_toj_color_post$logKappaProbeEffectMean/2
                 - (
                   ex_toj_color_post$logKappaInitialBiasInteractionEffectMean
                   - ex_toj_color_post$logKappaJudgementTypeInteractionEffectMean
                   + ex_toj_color_post$logKappaProbeInteractionEffectMean
                   
                   - ex_toj_color_post$logKappaInitialBiasJudgementTypeInteractionEffectMean
                   + ex_toj_color_post$logKappaInitialBiasProbeInteractionEffectMean
                   - ex_toj_color_post$logKappaJudgementTypeProbeInteractionEffectMean
                   - ex_toj_color_post$logKappaFourWayInteractionEffectMean
                 )/2
    )
    
    
    
    , exp(ex_toj_color_post$logKappaMean 
             - ex_toj_color_post$logKappaInitialBiasEffectMean/2
             + ex_toj_color_post$logKappaJudgementTypeEffectMean/2
             + ex_toj_color_post$logKappaProbeEffectMean/2
             + (
               -ex_toj_color_post$logKappaInitialBiasInteractionEffectMean
               + ex_toj_color_post$logKappaJudgementTypeInteractionEffectMean
               + ex_toj_color_post$logKappaProbeInteractionEffectMean
               
               - ex_toj_color_post$logKappaInitialBiasJudgementTypeInteractionEffectMean
               - ex_toj_color_post$logKappaInitialBiasProbeInteractionEffectMean
               + ex_toj_color_post$logKappaJudgementTypeProbeInteractionEffectMean
               - ex_toj_color_post$logKappaFourWayInteractionEffectMean
             )/2
    )
    -     exp(ex_toj_color_post$logKappaMean 
                 - ex_toj_color_post$logKappaInitialBiasEffectMean/2
                 + ex_toj_color_post$logKappaJudgementTypeEffectMean/2
                 + ex_toj_color_post$logKappaProbeEffectMean/2
                 - (
                   -ex_toj_color_post$logKappaInitialBiasInteractionEffectMean
                   + ex_toj_color_post$logKappaJudgementTypeInteractionEffectMean
                   + ex_toj_color_post$logKappaProbeInteractionEffectMean
                   
                   - ex_toj_color_post$logKappaInitialBiasJudgementTypeInteractionEffectMean
                   - ex_toj_color_post$logKappaInitialBiasProbeInteractionEffectMean
                   + ex_toj_color_post$logKappaJudgementTypeProbeInteractionEffectMean
                   - ex_toj_color_post$logKappaFourWayInteractionEffectMean
                 )/2
    )
    ,     exp(ex_toj_color_post$logKappaMean 
                 - ex_toj_color_post$logKappaInitialBiasEffectMean/2
                 + ex_toj_color_post$logKappaJudgementTypeEffectMean/2
                 - ex_toj_color_post$logKappaProbeEffectMean/2
                 + (
                   - ex_toj_color_post$logKappaInitialBiasInteractionEffectMean
                   + ex_toj_color_post$logKappaJudgementTypeInteractionEffectMean
                   - ex_toj_color_post$logKappaProbeInteractionEffectMean
                   
                   - ex_toj_color_post$logKappaInitialBiasJudgementTypeInteractionEffectMean
                   + ex_toj_color_post$logKappaInitialBiasProbeInteractionEffectMean
                   - ex_toj_color_post$logKappaJudgementTypeProbeInteractionEffectMean
                   + ex_toj_color_post$logKappaFourWayInteractionEffectMean
                 )/2
    )
    -     exp(ex_toj_color_post$logKappaMean 
                 - ex_toj_color_post$logKappaInitialBiasEffectMean/2
                 + ex_toj_color_post$logKappaJudgementTypeEffectMean/2
                 - ex_toj_color_post$logKappaProbeEffectMean/2
                 - (
                   - ex_toj_color_post$logKappaInitialBiasInteractionEffectMean
                   + ex_toj_color_post$logKappaJudgementTypeInteractionEffectMean
                   - ex_toj_color_post$logKappaProbeInteractionEffectMean
                   
                   - ex_toj_color_post$logKappaInitialBiasJudgementTypeInteractionEffectMean
                   + ex_toj_color_post$logKappaInitialBiasProbeInteractionEffectMean
                   - ex_toj_color_post$logKappaJudgementTypeProbeInteractionEffectMean
                   + ex_toj_color_post$logKappaFourWayInteractionEffectMean
                 )/2
    )
    ,     exp(ex_toj_color_post$logKappaMean 
                 - ex_toj_color_post$logKappaInitialBiasEffectMean/2
                 - ex_toj_color_post$logKappaJudgementTypeEffectMean/2
                 - ex_toj_color_post$logKappaProbeEffectMean/2
                 + (
                   - ex_toj_color_post$logKappaInitialBiasInteractionEffectMean
                   - ex_toj_color_post$logKappaJudgementTypeInteractionEffectMean
                   - ex_toj_color_post$logKappaProbeInteractionEffectMean
                   
                   + ex_toj_color_post$logKappaInitialBiasJudgementTypeInteractionEffectMean
                   + ex_toj_color_post$logKappaInitialBiasProbeInteractionEffectMean
                   + ex_toj_color_post$logKappaJudgementTypeProbeInteractionEffectMean
                   - ex_toj_color_post$logKappaFourWayInteractionEffectMean
                 )/2
    )
    -     exp(ex_toj_color_post$logKappaMean 
                 - ex_toj_color_post$logKappaInitialBiasEffectMean/2
                 - ex_toj_color_post$logKappaJudgementTypeEffectMean/2
                 - ex_toj_color_post$logKappaProbeEffectMean/2
                 - (
                   - ex_toj_color_post$logKappaInitialBiasInteractionEffectMean
                   - ex_toj_color_post$logKappaJudgementTypeInteractionEffectMean
                   - ex_toj_color_post$logKappaProbeInteractionEffectMean
                   
                   + ex_toj_color_post$logKappaInitialBiasJudgementTypeInteractionEffectMean
                   + ex_toj_color_post$logKappaInitialBiasProbeInteractionEffectMean
                   + ex_toj_color_post$logKappaJudgementTypeProbeInteractionEffectMean
                   - ex_toj_color_post$logKappaFourWayInteractionEffectMean
                 )/2
    )
    ,     exp(ex_toj_color_post$logKappaMean 
                 - ex_toj_color_post$logKappaInitialBiasEffectMean/2
                 - ex_toj_color_post$logKappaJudgementTypeEffectMean/2
                 + ex_toj_color_post$logKappaProbeEffectMean/2
                 + (
                   - ex_toj_color_post$logKappaInitialBiasInteractionEffectMean
                   - ex_toj_color_post$logKappaJudgementTypeInteractionEffectMean
                   + ex_toj_color_post$logKappaProbeInteractionEffectMean
                   
                   + ex_toj_color_post$logKappaInitialBiasJudgementTypeInteractionEffectMean
                   - ex_toj_color_post$logKappaInitialBiasProbeInteractionEffectMean
                   - ex_toj_color_post$logKappaJudgementTypeProbeInteractionEffectMean
                   + ex_toj_color_post$logKappaFourWayInteractionEffectMean
                 )/2
    )
    -     exp(ex_toj_color_post$logKappaMean 
                 - ex_toj_color_post$logKappaInitialBiasEffectMean/2
                 - ex_toj_color_post$logKappaJudgementTypeEffectMean/2
                 + ex_toj_color_post$logKappaProbeEffectMean/2
                 - (
                   - ex_toj_color_post$logKappaInitialBiasInteractionEffectMean
                   - ex_toj_color_post$logKappaJudgementTypeInteractionEffectMean
                   + ex_toj_color_post$logKappaProbeInteractionEffectMean
                   
                   + ex_toj_color_post$logKappaInitialBiasJudgementTypeInteractionEffectMean
                   - ex_toj_color_post$logKappaInitialBiasProbeInteractionEffectMean
                   - ex_toj_color_post$logKappaJudgementTypeProbeInteractionEffectMean
                   + ex_toj_color_post$logKappaFourWayInteractionEffectMean
                 )/2
    )  
    
  )
  , c("Left\nSecond\nLong", "Left\nSecond\nShort", "Left\nFirst\nShort", "Left\nFirst\nLong"
      , "Right\nSecond\nLong", "Right\nSecond\nShort", "Right\nFirst\nShort", "Right\nFirst\nLong")
  , y_lab = "\u03BA (Attended - Unattended)"
)
#-------------------------------- Four-way Interactions -----------------------------------#



#------------------------------------------------------------------------------------------#
#--------------------------------- Graphs -------------------------------------------------#
#------------------------------------------------------------------------------------------#

#---------------------------------- NCFs --------------------------------------------------#
# including effects 
yLeft = pnorm(
  -250:250
  , mean = ( median(ex_toj_color_post$population_pss_intercept_mean) - median(ex_toj_color_post$population_pss_effect_mean)/2 ) * 250
  , sd = ( exp( median(ex_toj_color_post$population_logjnd_intercept_mean) - median(ex_toj_color_post$population_logjnd_effect_mean)/2 )   ) * 250
)
yRight= pnorm(
  -250:250
  , mean = ( median(ex_toj_color_post$population_pss_intercept_mean) + median(ex_toj_color_post$population_pss_effect_mean)/2 ) * 250
  , sd = ( exp( median(ex_toj_color_post$population_logjnd_intercept_mean) + median(ex_toj_color_post$population_logjnd_effect_mean)/2 )   ) * 250
  
)
df = data.frame(SOA = -250:250, Prop = c(yRight, yLeft), Attend = c(rep("Right",501), rep("Left", 501)))


gg = ggplot(data = df, aes(y = Prop, x = SOA, colour = Attend))+
  geom_line(size = 1.25)+
  # scale_color_manual("Attend", values = c("red", "blue"))+
  scale_color_hue("Attend", l = c(60, 15), c = c(100, 50), h = c(240, 360) ) +
  labs(x = "SOA (ms)", y = "Proportion of 'Left' Responses")+
  theme_gray(base_size = 24)+
  theme(panel.grid.major = element_line(size = 1.5)
        ,panel.grid.minor = element_line(size = 1))
# define text to add
Text1 = textGrob(label = paste("Right"), gp = gpar(fontsize= 24))
Text2 = textGrob(label = paste("Left"), gp = gpar(fontsize= 24)) 
gg = gg+
  annotation_custom(grob = Text1,  xmin = -200, xmax = -200, ymin = -0.115, ymax = -0.115)+
  annotation_custom(grob = Text2,  xmin = 200, xmax = 200, ymin = -0.115, ymax = -0.115)
# Code to override clipping
gg2 <- ggplot_gtable(ggplot_build(gg))
gg2$layout$clip[gg2$layout$name=="panel"] <- "off"
grid.draw(gg2)
#---------------------------------- NCFs --------------------------------------------------#


# #---------------------------------- Posteriors --------------------------------------------#
# ### Rho
# pos_rhoMean_WithEffect = data.frame(
#   c( 
#     plogis(ex_toj_color_post$logitRhoMean + ex_toj_color_post$logitRhoProbeEffectMean/2
#            + (ex_toj_color_post$logitRhoEffectMean +  ex_toj_color_post$logitRhoProbeInteractionEffectMean)/2)
#     ,  plogis(ex_toj_color_post$logitRhoMean + ex_toj_color_post$logitRhoProbeEffectMean/2
#               - (ex_toj_color_post$logitRhoEffectMean +  ex_toj_color_post$logitRhoProbeInteractionEffectMean)/2)
#     ,  plogis(ex_toj_color_post$logitRhoMean - ex_toj_color_post$logitRhoProbeEffectMean/2
#               + (ex_toj_color_post$logitRhoEffectMean -  ex_toj_color_post$logitRhoProbeInteractionEffectMean)/2)
#     ,  plogis(ex_toj_color_post$logitRhoMean - ex_toj_color_post$logitRhoProbeEffectMean/2
#               - (ex_toj_color_post$logitRhoEffectMean -  ex_toj_color_post$logitRhoProbeInteractionEffectMean)/2)
#   ) 
#   , c(rep("AttendedLong",80000), rep("UnattendedLong",80000), rep("AttendedShort",80000), rep("UnattendedShort",80000))
# )
# names(pos_rhoMean_WithEffect) = c("rhoMean", "Effect")
# 
# 
# # overlapping
# ggplot(pos_rhoMean_WithEffect, aes(x = rhoMean, ..density.., fill = Effect))+
#   geom_density(data = pos_rhoMean_WithEffect[pos_rhoMean_WithEffect$Effect == "AttendedLong",],alpha = 0.5)+
#   geom_density(data = pos_rhoMean_WithEffect[pos_rhoMean_WithEffect$Effect == "UnattendedLong",],alpha = 0.5)+
#   geom_density(data = pos_rhoMean_WithEffect[pos_rhoMean_WithEffect$Effect == "AttendedShort",],alpha = 0.5)+
#   geom_density(data = pos_rhoMean_WithEffect[pos_rhoMean_WithEffect$Effect == "UnattendedShort",],alpha = 0.5)+
#   scale_fill_hue("Effect", l = c(90, 45, 70, 30) , c = c(100, 50, 100, 50) ) +
#   labs(x = "Probability of Memory Population Mean", y = "Density", colour = "")+
#   theme_gray(base_size = 24)+
#   theme(panel.grid.major = element_line(size = 1.5)
#         ,panel.grid.minor = element_line(size = 1))
# 
# ### Kappa
# pos_kappaMean_WithEffect = data.frame(
#   c( 
#     plogis(ex_toj_color_post$logKappaMean + ex_toj_color_post$logKappaProbeEffectMean/2
#            + (ex_toj_color_post$logKappaEffectMean +  ex_toj_color_post$logKappaProbeInteractionEffectMean)/2)
#     ,  plogis(ex_toj_color_post$logKappaMean + ex_toj_color_post$logKappaProbeEffectMean/2
#               - (ex_toj_color_post$logKappaEffectMean +  ex_toj_color_post$logKappaProbeInteractionEffectMean)/2)
#     ,  plogis(ex_toj_color_post$logKappaMean - ex_toj_color_post$logKappaProbeEffectMean/2
#               + (ex_toj_color_post$logKappaEffectMean -  ex_toj_color_post$logKappaProbeInteractionEffectMean)/2)
#     ,  plogis(ex_toj_color_post$logKappaMean - ex_toj_color_post$logKappaProbeEffectMean/2
#               - (ex_toj_color_post$logKappaEffectMean -  ex_toj_color_post$logKappaProbeInteractionEffectMean)/2)
#   ) 
#   , c(rep("AttendedLong",80000), rep("UnattendedLong",80000), rep("AttendedShort",80000), rep("UnattendedShort",80000))
# )
# names(pos_kappaMean_WithEffect) = c("kappaMean", "Effect")
# 
# 
# # overlapping
# ggplot(pos_rhoMean_WithEffect, aes(x = kappaMean, ..density.., fill = Effect))+
#   geom_density(data = pos_kappaMean_WithEffect[pos_kappaMean_WithEffect$Effect == "AttendedLong",],alpha = 0.5)+
#   geom_density(data = pos_kappaMean_WithEffect[pos_kappaMean_WithEffect$Effect == "UnattendedLong",],alpha = 0.5)+
#   geom_density(data = pos_kappaMean_WithEffect[pos_kappaMean_WithEffect$Effect == "AttendedShort",],alpha = 0.5)+
#   geom_density(data = pos_kappaMean_WithEffect[pos_kappaMean_WithEffect$Effect == "UnattendedShort",],alpha = 0.5)+
#   scale_fill_hue("Effect", l = c(90, 45, 70, 30) , c = c(100, 50, 100, 50) ) +
#   labs(x = "Fidelity of Memory Population Mean", y = "Density", colour = "")+
#   theme_gray(base_size = 24)+
#   theme(panel.grid.major = element_line(size = 1.5)
#         ,panel.grid.minor = element_line(size = 1))
# #---------------------------------- Posteriors --------------------------------------------#





