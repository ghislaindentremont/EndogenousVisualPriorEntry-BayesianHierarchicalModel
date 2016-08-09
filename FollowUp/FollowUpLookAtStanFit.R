# library(shinystan)
library(plyr)
library(coda)
library(ggplot2)
library(ggmcmc)
library(CircStats)
library(grid)
library(sprintfr)
library(gtools)

setwd("~/Documents/TOJ/Follow-Up")
load("FollowUptoj_color_post_Aug8th2016")
load("FollowUp_color_trials.Rdata")
load("FollowUp_toj_trials.Rdata")
source("../EndogenousVisualPriorEntry-BayesianHierarchicalModel/functions.R")


############################################################################################
####                                        Diagnostics                                 ####
############################################################################################
# convert stanfit sample to dataframe table 
gg_toj_color_post = ggs(toj_color_post)

# list of parameters to examine
param_list = c(
               "population_pss_intercept_mean"
               , "population_pss_attention_effect_mean"
               # , "population_pss_initial_bias_effect_mean"
               , "population_pss_judgement_type_effect_mean"
               , "population_pss_probe_duration_effect_mean"
               # , "population_pss_attention_initial_bias_interaction_effect_mean"
               , "population_pss_attention_judgement_type_interaction_effect_mean"	
               , "population_pss_attention_probe_duration_interaction_effect_mean"
               , "population_log_jnd_intercept_mean"
               , "population_log_jnd_attention_effect_mean"
               # , "population_log_jnd_initial_bias_effect_mean"
               , "population_log_jnd_judgement_type_effect_mean"
               , "population_log_jnd_probe_duration_effect_mean"
               # , "population_log_jnd_attention_initial_bias_interaction_effect_mean"
               , "population_log_jnd_attention_judgement_type_interaction_effect_mean"	
               , "population_log_jnd_attention_probe_duration_interaction_effect_mean"
               , "population_logit_rho_intercept_mean"
               , "population_logit_rho_attention_effect_mean"
               , "population_logit_rho_probe_duration_effect_mean"
               , "population_logit_rho_attention_probe_duration_interaction_effect_mean"
               # , "population_logit_rho_initial_bias_effect_mean"
               # , "population_logit_rho_attention_initial_bias_interaction_effect_mean"
               , "population_log_kappa_intercept_mean"
               , "population_log_kappa_attention_effect_mean"
               , "population_log_kappa_probe_duration_effect_mean"
               , "population_log_kappa_attention_probe_duration_interaction_effect_mean"
               # , "population_log_kappa_initial_bias_effect_mean"
               # , "population_log_kappa_attention_initial_bias_interaction_effect_mean"
               , "population_pss_intercept_sd"
               , "population_pss_effect_sd" 
               , "population_log_jnd_intercept_sd"
               , "population_log_jnd_effect_sd" 
               , "population_logit_rho_intercept_sd"
               , "population_log_kappa_intercept_sd"
               , "population_logit_rho_effect_sd"
               , "population_log_kappa_effect_sd"
)

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
real_toj = aggregate(left_first_TF ~ soa2 + block_bias + toj_judgement_type + longprobe, data = toj_trials, FUN = mean)
#-------------------------------------- TOJ Actual Data -----------------------------------#


#-------------------------------------- TOJ Simulated Data --------------------------------#

### Get PSS Parameters
pss_intercept_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "population_pss_intercept_mean",]$value
pss_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "population_pss_attention_effect_mean",]$value

pss_judgement_type_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "population_pss_judgement_type_effect_mean",]$value
pss_judgement_type_interaction_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "population_pss_attention_judgement_type_interaction_effect_mean",]$value

pss_probe_duration_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "population_pss_probe_duration_effect_mean",]$value
pss_probe_duration_interaction_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "population_pss_attention_probe_duration_interaction_effect_mean",]$value

# pss_initial_bias_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "population_pss_initial_bias_effect_mean",]$value
# pss_initial_bias_interaction_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "population_pss_attention_initial_bias_interaction_effect_mean",]$value

# judgement type
pss_right_second_long_mean_reps = get_condition_mean_sample(
  ( pss_intercept_mean 
    + pss_judgement_type_effect_mean/2
    + pss_probe_duration_effect_mean/2
    )
   # - pss_initial_bias_effect_mean/2 )
  , ( pss_effect_mean 
      + pss_judgement_type_interaction_effect_mean
      + pss_probe_duration_interaction_effect_mean 
      )
      # - pss_initial_bias_interaction_effect_mean)
  , TRUE
  , "null"
)

### Get JND Parameters
log_jnd_intercept_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "population_log_jnd_intercept_mean",]$value
log_jnd_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "population_log_jnd_attention_effect_mean",]$value

log_jnd_judgement_type_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "population_log_jnd_judgement_type_effect_mean",]$value
log_jnd_judgement_type_interaction_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "population_log_jnd_attention_judgement_type_interaction_effect_mean",]$value

log_jnd_probe_duration_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "population_log_jnd_probe_duration_effect_mean",]$value
log_jnd_probe_duration_interaction_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "population_log_jnd_attention_probe_duration_interaction_effect_mean",]$value

# log_jnd_initial_bias_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "population_log_jnd_initial_bias_effect_mean",]$value
# log_jnd_initial_bias_interaction_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "population_log_jnd_attention_initial_bias_interaction_effect_mean",]$value

# judgement type
log_jnd_right_second_long_mean_reps = get_condition_mean_sample(
  ( log_jnd_intercept_mean 
  + log_jnd_judgement_type_effect_mean/2
  + log_jnd_probe_duration_effect_mean/2 )
  # - log_jnd_initial_bias_effect_mean/2 )
  , ( log_jnd_effect_mean 
     + log_jnd_judgement_type_interaction_effect_mean
     + log_jnd_probe_duration_interaction_effect_mean )
     # - log_jnd_initial_bias_interaction_effect_mean )
  , TRUE
  , "log"
)
#-------------------------------------- TOJ Simulated Data --------------------------------#


#-------------------------------------- Do TOJ PPC ----------------------------------------#
SOAs = c(-250, -150, -100, -50, -17, 17, 50, 100, 150, 250)

# judgement type
do_toj_ppc(
  pss_right_second_long_mean_reps
  , log_jnd_right_second_long_mean_reps
  , "'which second?' & attend right & long probe duration"
  , c("toj_judgement_type", "block_bias", "longprobe")
  , c("second", "RIGHT", "TRUE")
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
rho_intercept_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "population_logit_rho_intercept_mean",]$value
rho_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "population_logit_rho_attention_effect_mean",]$value

# rho_judgement_type_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "population_logit_rho_judgement_type_effect_mean",]$value
# rho_judgement_type_interaction_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "population_logit_rho_attention_judgement_type_interaction_effect_mean",]$value

rho_probe_duration_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "population_logit_rho_probe_duration_effect_mean",]$value
rho_probe_duration_interaction_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "population_logit_rho_attention_probe_duration_interaction_effect_mean",]$value

# rho_initial_bias_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "logitRhoInitialBiasEffectMean",]$value
# rho_initial_bias_interaction_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "population_rho_attention_probe_duration_interaction_effect_mean",]$value

# get condition combination
rho_right_long_mean_reps = get_condition_mean_sample(
  ( rho_intercept_mean 
    # + rho_judgement_type_effect_mean/2 
    + rho_probe_duration_effect_mean/2 )
#     - rho_initial_bias_effect_mean/2 )
  , ( rho_effect_mean 
      # + rho_judgement_type_interaction_effect_mean 
      + rho_probe_duration_interaction_effect_mean )
#       - rho_initial_bias_interaction_effect_mean )
  , TRUE
  , "logit"
)

### Get Kappa Parameters
kappa_intercept_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "population_log_kappa_intercept_mean",]$value
kappa_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "population_log_kappa_attention_effect_mean",]$value

# kappa_judgement_type_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "logKappaJudgementTypeEffectMean",]$value
# kappa_judgement_type_interaction_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "population_log_kappa_attention_judgement_type_interaction_effect_mean",]$value

kappa_probe_duration_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "population_log_kappa_probe_duration_effect_mean",]$value
kappa_probe_duration_interaction_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "population_log_kappa_attention_probe_duration_interaction_effect_mean",]$value

# kappa_initial_bias_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "logKappaInitialBiasEffectMean",]$value
# kappa_initial_bias_interaction_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "logKappaAttentionInitialBiasEffectMean",]$value

# get condition combination
kappa_right_long_mean_reps = get_condition_mean_sample(
  ( kappa_intercept_mean 
    + kappa_probe_duration_effect_mean/2)
#     - kappa_judgement_type_effect_mean/2 
#     - kappa_initial_bias_effect_mean/2 )
  , ( kappa_effect_mean 
        + kappa_probe_duration_interaction_effect_mean )
#       - kappa_judgement_type_interaction_effect_mean 
#       - kappa_initial_bias_interaction_effect_mean )
  , TRUE
  , "log_free"
)
#-------------------------------------- Color Simulated Data ------------------------------#


#-------------------------------------- Do Color PPC --------------------------------------#
do_color_ppc(
  rho_right_long_mean_reps
  , kappa_right_long_mean_reps
  , "attend right & probe duration long"
  , c("block_bias", "longprobe")
  , c("RIGHT", "TRUE")
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
# (2) population_pss_attention_effect_mean          
# (3) population_log_jnd_intercept_mean    
# (4) population_log_jnd_attention_effect_mean     
# (5) population_logit_rho_intercept_mean                         
# (6) population_log_kappa_intercept_mean                        
# (7) population_logit_rho_attention_effect_mean                 
# (8) population_log_kappa_attention_effect_mean     

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
# (2) population_pss_attention_effect_mean          
# (3) population_log_jnd_intercept_mean    
# (4) population_log_jnd_attention_effect_mean     
# (5) population_logit_rho_intercept_mean                         
# (6) population_log_kappa_intercept_mean                        
# (7) population_logit_rho_attention_effect_mean                 
# (8) population_log_kappa_attention_effect_mean 
# JND and PSS intercepts
library(reshape)
betas2 = data.frame(value = ex_toj_color_post$beta)
betas2$iteration = rownames(betas2)
betas = melt( betas2 )
betas$parameter = rep( c(
  "population_pss_intercept_mean"      
  , "population_pss_attention_effect_mean"          
  , "population_log_jnd_intercept_mean"    
  , "population_log_jnd_attention_effect_mean"     
  , "population_logit_rho_intercept_mean"                       
  , "population_log_kappa_intercept_mean"                       
  , "population_logit_rho_attention_effect_mean"                 
  , "population_log_kappa_attention_effect_mean"
)
, times = 1
, each = nrow(betas2)*length(unique(betas$variable))/8  # 8 is number of parameters 
)  
betas$participant = rep(c(1:length(unique(toj_trials$id))), times = 8, each = nrow(betas2))
#-------------------------------------- Get Betas -----------------------------------------#


#---------------------------- Rho vs. PSS Effects -----------------------------------------#
psseffect = extract_samples("population_pss_attention_effect_mean")

psseffectsd = extract_samples("population_pss_effect_sd")

pssjudgementinteraction = extract_samples("population_pss_attention_judgement_type_interaction_effect_mean")

pssinitialbiasinteraction = extract_samples("population_pss_attention_initial_bias_interaction_effect_mean")

pssprobeinteraction = extract_samples("population_pss_attention_probe_duration_interaction_effect_mean")

judgementfactor = ifelse(aggregate(toj_judgement_type ~ id, data = toj_trials, FUN = unique)$toj_judgement_type == "first", -1, 1)

probefactor = ifelse(aggregate(longprobe ~ id, data = color_trials, FUN = unique)$longprobe == "FALSE", -1, 1)

psseffect_ids = ddply(
  .data = betas
  , .variables = .(participant)
  , .fun = function(x){
    i = unique(x$participant)
    x_use = x[x$parameter ==  "population_pss_attention_effect_mean",]$value
    psseffect_use = median(psseffect)  + median(psseffectsd)*median(x_use)  + median(pssjudgementinteraction)*judgementfactor[i]+ median(pssprobeinteraction)*probefactor[i]
    df = data.frame(psseffect_use*250, judgementfactor[i], probefactor[i])
    names(df) = c("psseffect",  "judgementfactor", "probefactor")
    return(df)
  }
)

logitrhoeffect = extract_samples("population_logit_rho_attention_effect_mean")

logitrhoeffectsd = extract_samples("population_logit_rho_effect_sd")

# logitrhojudgementinteractioneffect = extract_samples("population_logit_rho_attention_judgement_type_interaction_effect_mean")

logitrhoprobeinteractioneffect = extract_samples("population_logit_rho_attention_probe_duration_interaction_effect_mean")

rhoeffect_ids = ddply(
  .data = betas
  , .variables = .(participant)
  , .fun = function(x){
    i = unique(x$participant)
    x_use = x[x$parameter == "population_logit_rho_attention_effect_mean",]$value
    logitrhoeffect_use =  median(logitrhoeffect) + median(logitrhoeffectsd)*median(x_use)+ median(logitrhoprobeinteractioneffect)*probefactor[i] #+ median(logitrhojudgementinteractioneffect)*judgementfactor[i]
    df = data.frame(logitrhoeffect_use, probefactor[i])
    names(df) = c("logitrhoeffect", "probefactor")
    return(df)
  }
)

psseffect_v_rhoeffect = merge(rhoeffect_ids, psseffect_ids)

ggplot(data = psseffect_v_rhoeffect, aes(y = psseffect, x = logitrhoeffect, colour = factor(judgementfactor), shape = factor(judgementfactor), fill = factor(probefactor)))+
  scale_y_continuous(name = "PSS Effect Mean")+
  scale_x_continuous(name = "Logit \u03C1 Effect Mean")+
  geom_vline(xintercept = 0, linetype = 2, size = 1)+
  geom_hline(yintercept = 0, linetype = 2, size = 1)+
  geom_point(size = 4)+
  scale_shape_manual(name = "Judgement\nType", labels = c("First", "Second") , values = c(21,22) )+
  scale_colour_manual(name = "Judgement\nType", labels =c("First", "Second"), values = c("red", "blue") )+
  scale_fill_manual(name = "Probe\nDuration", labels = c("Short", "Long"), values = c("black", "white")) +
  theme_gray(base_size = 30)+
  theme(panel.grid.major = element_line(size = 1.5)
        ,panel.grid.minor = element_line(size = 1))

# NOTE: RECALL THAT COR BELOW IS FOR ATTENTION EFFECT ONLY
# get_corr(
#   "value.2.7"
#   , "Logit \u03C1 vs. PSS Effect Means"
# )
#---------------------------- Rho vs. PSS Effects -----------------------------------------#


#---------------------------- Kappa vs. PSS Effects ---------------------------------------#
logkappaeffect = extract_samples("population_log_kappa_attention_effect_mean")

logkappaeffectsd = extract_samples("population_log_kappa_effect_sd")

# logkappajudgementinteractioneffect = extract_samples("population_log_kappa_attention_judgement_type_interaction_effect_mean")

logkappaprobeinteractioneffect = extract_samples("population_log_kappa_attention_probe_duration_interaction_effect_mean")

kappaeffect_ids = ddply(
  .data = betas
  , .variables = .(participant)
  , .fun = function(x){
    i = unique(x$participant)
    x_use = x[x$parameter == "population_log_kappa_attention_effect_mean",]$value
    logkappaeffect_use =  median(logkappaeffect) + median(logkappaeffectsd)*median(x_use) + median(logkappaprobeinteractioneffect)*probefactor[i] # + median(logkappajudgementinteractioneffect)*judgementfactor[i]
    df = data.frame(logkappaeffect_use, probefactor[i])
    names(df) = c("kappaeffect", "probefactor")
    return(df)
  }
)

psseffect_v_kappaeffect = merge(kappaeffect_ids, psseffect_ids)

ggplot(data = psseffect_v_kappaeffect, aes(y = psseffect, x = kappaeffect, colour = factor(judgementfactor), shape = factor(judgementfactor), fill = factor(probefactor)))+ 
  scale_y_continuous(name = "PSS Effect Mean")+
  scale_x_continuous(name = "Log \u03BA Effect Mean")+
  geom_vline(xintercept = 0, linetype = 2, size = 1)+
  geom_hline(yintercept = 0, linetype = 2, size = 1)+
  geom_point(size = 4)+
  scale_shape_manual(name = "Judgement\nType", labels = c("First", "Second") , values = c(21,22) )+
  scale_colour_manual(name = "Judgement\nType", labels =c("First", "Second"), values = c("red", "blue") )+
  scale_fill_manual(name = "Probe\nDuration", labels = c("Short", "Long"), values = c("black", "white")) +
  theme_gray(base_size = 30)+
  theme(panel.grid.major = element_line(size = 1.5)
        ,panel.grid.minor = element_line(size = 1))

# get_corr(
#   "value.2.8"
#   , "Log \u03BA vs. PSS Effect Means"
# )        
#---------------------------- Kappa vs. PSS Effects ---------------------------------------#



#------------------------------------------------------------------------------------------#
#--------------------------------- SOA Scale ----------------------------------------------#
#------------------------------------------------------------------------------------------#

#---------------------------------- SOA Intercepts ----------------------------------------#
get_violin(
  c(
    ex_toj_color_post$population_pss_intercept_mean * 250
    , exp( ex_toj_color_post$population_log_jnd_intercept_mean ) * 250
    )
  , c("PSS Intercept Mean", "JND Intercept Mean")
  , y_lab = "SOA (ms)"
  , hline = FALSE
  , facet = TRUE
)
get_95_HDI(exp( ex_toj_color_post$population_log_jnd_intercept_mean ) * 250)
#---------------------------------- SOA Intercepts ----------------------------------------#


#---------------------------------- Main Effects ------------------------------------------#
# effect of attention on PSS and JND
get_violin(
  c(
  ( (ex_toj_color_post$population_pss_intercept_mean + ex_toj_color_post$population_pss_attention_effect_mean/2) 
    - (ex_toj_color_post$population_pss_intercept_mean - ex_toj_color_post$population_pss_attention_effect_mean/2) ) * 250
  , ( exp( ex_toj_color_post$population_log_jnd_intercept_mean + ex_toj_color_post$population_log_jnd_attention_effect_mean/2 )
  - exp( ex_toj_color_post$population_log_jnd_intercept_mean - ex_toj_color_post$population_log_jnd_attention_effect_mean/2  ) ) * 250 
  )
  , c("PSS\nAttention Effect Mean"  , "JND\nAttention Effect Mean")
  , y_lab = "SOA (Right - Left; ms)"
)
get_95_HDI( ( exp( ex_toj_color_post$population_log_jnd_intercept_mean + ex_toj_color_post$population_log_jnd_attention_effect_mean/2 )
              - exp( ex_toj_color_post$population_log_jnd_intercept_mean - ex_toj_color_post$population_log_jnd_attention_effect_mean/2  ) ) * 250 
)

# effect of judgement type (Q) on PSS and JND
get_violin(
  c(
  ( ex_toj_color_post$population_pss_judgement_type_effect_mean ) * 250 
  , ( exp( ex_toj_color_post$population_log_jnd_intercept_mean + ex_toj_color_post$population_log_jnd_judgement_type_effect_mean/2 )
      - exp( ex_toj_color_post$population_log_jnd_intercept_mean - ex_toj_color_post$population_log_jnd_judgement_type_effect_mean/2  ) ) * 250
  )
  , c("PSS Judgement\nType Effect Mean"  , "JND Judgement\nType Effect Mean")
  , y_lab = "SOA (Second - First; ms)"
)
get_95_HDI(
  ( exp( ex_toj_color_post$population_log_jnd_intercept_mean + ex_toj_color_post$population_log_jnd_judgement_type_effect_mean/2 )
             - exp( ex_toj_color_post$population_log_jnd_intercept_mean - ex_toj_color_post$population_log_jnd_judgement_type_effect_mean/2  ) ) * 250
)

# # effect of initial bias on PSS and JND
# get_violin(
#   c(
#   ( ex_toj_color_post$population_pss_initial_bias_effect_mean ) * 250
#   , ( exp( ex_toj_color_post$population_log_jnd_intercept_mean + ex_toj_color_post$population_log_jnd_initial_bias_effect_mean/2 )
#       - exp( ex_toj_color_post$population_log_jnd_intercept_mean - ex_toj_color_post$population_log_jnd_initial_bias_effect_mean/2  ) ) * 250 
#   )
#   , c("PSS Initial Probe\nBias Effect Mean"  , "JND Initial Probe\nBias Effect Mean")
#   , y_lab = "SOA (Left - Right; ms)"
# )
# get_95_HDI(
#   ( exp( ex_toj_color_post$population_log_jnd_intercept_mean + ex_toj_color_post$population_log_jnd_initial_bias_effect_mean/2 )
#     - exp( ex_toj_color_post$population_log_jnd_intercept_mean - ex_toj_color_post$population_log_jnd_initial_bias_effect_mean/2  ) ) * 250 
# )
#---------------------------------- Main Effects ---------------------=--------------------#


#------------------------------- Two-way Interactions -------------------------------------#
#  effect of interaction between judgement type and attention on PSS 
get_violin(
  2* c(  # MULTIPLY INTERACTION BY 2 TO GET DIFFERENCE IN DIFFERENCES (Meaningful)
  (ex_toj_color_post$population_pss_attention_judgement_type_interaction_effect_mean) * 250
  ,   ( exp( ex_toj_color_post$population_log_jnd_intercept_mean + ex_toj_color_post$population_log_jnd_attention_judgement_type_interaction_effect_mean/2 )
        - exp( ex_toj_color_post$population_log_jnd_intercept_mean - ex_toj_color_post$population_log_jnd_attention_judgement_type_interaction_effect_mean/2  ) ) * 250 
  )
  , c("PSS Attention\n& Judgement Type\nInteraction Effect Mean"  , "JND Attention\n& Judgement Type\nInteraction Effect Mean")
  , y_lab = "SOA (ms)"
)
get_95_HDI(
  2* (  # MULTIPLY INTERACTION BY 2 TO GET DIFFERENCE IN DIFFERENCES (Meaningful)
    exp( ex_toj_color_post$population_log_jnd_intercept_mean + ex_toj_color_post$population_log_jnd_attention_judgement_type_interaction_effect_mean/2 )
    - exp( ex_toj_color_post$population_log_jnd_intercept_mean - ex_toj_color_post$population_log_jnd_attention_judgement_type_interaction_effect_mean/2  ) ) * 250 
)

get_violin(
  c(
  (ex_toj_color_post$population_pss_intercept_mean + ex_toj_color_post$population_pss_judgement_type_effect_mean/2
  + (ex_toj_color_post$population_pss_attention_effect_mean + ex_toj_color_post$population_pss_attention_judgement_type_interaction_effect_mean)/2 ) *250
  - (ex_toj_color_post$population_pss_intercept_mean + ex_toj_color_post$population_pss_judgement_type_effect_mean/2
  - (ex_toj_color_post$population_pss_attention_effect_mean + ex_toj_color_post$population_pss_attention_judgement_type_interaction_effect_mean)/2 ) * 250
  ,   (ex_toj_color_post$population_pss_intercept_mean - ex_toj_color_post$population_pss_judgement_type_effect_mean/2
       + (ex_toj_color_post$population_pss_attention_effect_mean - ex_toj_color_post$population_pss_attention_judgement_type_interaction_effect_mean)/2 ) *250
  - (ex_toj_color_post$population_pss_intercept_mean - ex_toj_color_post$population_pss_judgement_type_effect_mean/2
     - (ex_toj_color_post$population_pss_attention_effect_mean - ex_toj_color_post$population_pss_attention_judgement_type_interaction_effect_mean)/2 ) *250
  )
  , c("PSS Attention Effect\nGiven Second", "PSS Attention Effect\nGiven First")
  , y_lab = "SOA (Right - Left; ms)"
)
get_95_HDI(
  (ex_toj_color_post$population_pss_intercept_mean - ex_toj_color_post$population_pss_judgement_type_effect_mean/2
   + (ex_toj_color_post$population_pss_attention_effect_mean - ex_toj_color_post$population_pss_attention_judgement_type_interaction_effect_mean)/2 ) *250
  - (ex_toj_color_post$population_pss_intercept_mean - ex_toj_color_post$population_pss_judgement_type_effect_mean/2
     - (ex_toj_color_post$population_pss_attention_effect_mean - ex_toj_color_post$population_pss_attention_judgement_type_interaction_effect_mean)/2 ) *250
)

#  effect of interaction between probe duration and attention on PSS 
get_violin(
  2* c(  # MULTIPLY INTERACTION BY 2 TO GET DIFFERENCE IN DIFFERENCES (Meaningful)
    (ex_toj_color_post$population_pss_attention_probe_duration_interaction_effect_mean) * 250
    ,   ( exp( ex_toj_color_post$population_log_jnd_intercept_mean + ex_toj_color_post$population_log_jnd_attention_probe_duration_interaction_effect_mean/2 )
          - exp( ex_toj_color_post$population_log_jnd_intercept_mean - ex_toj_color_post$population_log_jnd_attention_probe_duration_interaction_effect_mean/2  ) ) * 250 
  )
  , c("PSS Attention\n& Probe Duration\nInteraction Effect Mean"  , "JND Attention\n& Probe Duration\nInteraction Effect Mean")
  , y_lab = "SOA (ms)"
)
get_95_HDI(
  2*(  # MULTIPLY INTERACTION BY 2 TO GET DIFFERENCE IN DIFFERENCES (Meaningful) 
    exp( ex_toj_color_post$population_log_jnd_intercept_mean + ex_toj_color_post$population_log_jnd_attention_probe_duration_interaction_effect_mean/2 )
    - exp( ex_toj_color_post$population_log_jnd_intercept_mean - ex_toj_color_post$population_log_jnd_attention_probe_duration_interaction_effect_mean/2  ) ) * 250 
)

get_violin(
  c(
    (ex_toj_color_post$population_pss_intercept_mean + ex_toj_color_post$population_pss_probe_duration_effect_mean/2
     + (ex_toj_color_post$population_pss_attention_effect_mean + ex_toj_color_post$population_pss_attention_probe_duration_interaction_effect_mean)/2 ) *250
    - (ex_toj_color_post$population_pss_intercept_mean + ex_toj_color_post$population_pss_probe_duration_effect_mean/2
       - (ex_toj_color_post$population_pss_attention_effect_mean + ex_toj_color_post$population_pss_attention_probe_duration_interaction_effect_mean)/2 ) * 250
    ,   (ex_toj_color_post$population_pss_intercept_mean - ex_toj_color_post$population_pss_probe_duration_effect_mean/2
         + (ex_toj_color_post$population_pss_attention_effect_mean - ex_toj_color_post$population_pss_attention_probe_duration_interaction_effect_mean)/2 ) *250
    - (ex_toj_color_post$population_pss_intercept_mean - ex_toj_color_post$population_pss_probe_duration_effect_mean/2
       - (ex_toj_color_post$population_pss_attention_effect_mean - ex_toj_color_post$population_pss_attention_probe_duration_interaction_effect_mean)/2 ) *250
  )
  , c("PSS Attention Effect\nGiven Long", "PSS Attention Effect\nGiven Short")
  , y_lab = "SOA (Right - Left; ms)"
)
get_95_HDI(
  (ex_toj_color_post$population_pss_intercept_mean - ex_toj_color_post$population_pss_probe_duration_effect_mean/2
   + (ex_toj_color_post$population_pss_attention_effect_mean - ex_toj_color_post$population_pss_attention_probe_duration_interaction_effect_mean)/2 ) *250
  - (ex_toj_color_post$population_pss_intercept_mean - ex_toj_color_post$population_pss_probe_duration_effect_mean/2
     - (ex_toj_color_post$population_pss_attention_effect_mean - ex_toj_color_post$population_pss_attention_probe_duration_interaction_effect_mean)/2 ) *250
)
#------------------------------- Two-way Interactions -------------------------------------#


#------------------------------------------------------------------------------------------#
#--------------------------------- Rho Scale ----------------------------------------------#
#------------------------------------------------------------------------------------------#

#---------------------------------- Rho Intercept -----------------------------------------#
get_violin(
  plogis(ex_toj_color_post$population_logit_rho_intercept_mean)
  , "Probability of Encoding Intercept Mean"
  , y_lab = "\u03C1"
  , hline = FALSE
)
#---------------------------------- Rho Intercept -----------------------------------------#


#-------------------------------- Main Effects --------------------------------------------#
get_violin(
  ( plogis(ex_toj_color_post$population_logit_rho_intercept_mean + ex_toj_color_post$population_logit_rho_attention_effect_mean/2 )
    - plogis(ex_toj_color_post$population_logit_rho_intercept_mean - ex_toj_color_post$population_logit_rho_attention_effect_mean/2 ) )
  , "Probability of Encoding Attention Effect Mean"
  , y_lab = "\u03C1 (Attended - Unattended)"
)

get_violin(
    ( plogis(ex_toj_color_post$population_logit_rho_intercept_mean + ex_toj_color_post$population_logit_rho_probe_duration_effect_mean/2 )
    - plogis(ex_toj_color_post$population_logit_rho_intercept_mean - ex_toj_color_post$population_logit_rho_probe_duration_effect_mean/2 ) )
  , "Probability of Encoding\nProbe Duration Effect Mean"
  , y_lab = "\u03C1 (Long - Short)"
)

# get_violin(
#   ( plogis(ex_toj_color_post$population_logit_rho_intercept_mean + ex_toj_color_post$population_logit_rho_judgement_type_effect_mean/2 )
#     - plogis(ex_toj_color_post$population_logit_rho_intercept_mean - ex_toj_color_post$population_logit_rho_judgement_type_effect_mean/2 ) )
#   , "Probability of Memory\nJudgement Type Effect Mean"
#   , y_lab = "\u03C1 (Second - First)"
# )
# 
# get_violin(
#   ( plogis(ex_toj_color_post$population_logit_rho_intercept_mean + ex_toj_color_post$logitRhoInitialBiasEffectMean/2 )
#     - plogis(ex_toj_color_post$population_logit_rho_intercept_mean - ex_toj_color_post$logitRhoInitialBiasEffectMean/2 ) )
#   , "Probability of Memory\nInitial Probe Bias Effect Mean"
#   , y_lab = "\u03C1 (Left - Right)"
# )
#-------------------------------- Main Effects --------------------------------------------#


#------------------------------- Two-way Interactions -------------------------------------#
get_violin(
  2* (  # MULTIPLY INTERACTION BY 2 TO GET DIFFERENCE IN DIFFERENCES (Meaningful)
    plogis(ex_toj_color_post$population_logit_rho_intercept_mean + ex_toj_color_post$population_logit_rho_attention_probe_duration_interaction_effect_mean/2 )
    - plogis(ex_toj_color_post$population_logit_rho_intercept_mean - ex_toj_color_post$population_logit_rho_attention_probe_duration_interaction_effect_mean/2 ) )
  , "Probability of Encoding\nProbe Duration Attention\nInteraction Effect Mean"
  , y_lab = "\u03C1"
)

get_violin(
  c(
  plogis(ex_toj_color_post$population_logit_rho_intercept_mean + ex_toj_color_post$population_logit_rho_probe_duration_effect_mean/2 
    + (ex_toj_color_post$population_logit_rho_attention_effect_mean + ex_toj_color_post$population_logit_rho_attention_probe_duration_interaction_effect_mean )/2 )
  - plogis(ex_toj_color_post$population_logit_rho_intercept_mean + ex_toj_color_post$population_logit_rho_probe_duration_effect_mean/2 
    - (ex_toj_color_post$population_logit_rho_attention_effect_mean + ex_toj_color_post$population_logit_rho_attention_probe_duration_interaction_effect_mean )/2 )
  , plogis(ex_toj_color_post$population_logit_rho_intercept_mean - ex_toj_color_post$population_logit_rho_probe_duration_effect_mean/2 
     + (ex_toj_color_post$population_logit_rho_attention_effect_mean - ex_toj_color_post$population_logit_rho_attention_probe_duration_interaction_effect_mean )/2 )
    - plogis(ex_toj_color_post$population_logit_rho_intercept_mean - ex_toj_color_post$population_logit_rho_probe_duration_effect_mean/2 
     - (ex_toj_color_post$population_logit_rho_attention_effect_mean - ex_toj_color_post$population_logit_rho_attention_probe_duration_interaction_effect_mean )/2 ) 
  )
  , c("Probability of Encoding\nAttention Effect \nGiven Long Probe Duration","Probability of Encoding\nAttention Effect\nGiven Short Probe Duration")
  , y_lab = "\u03C1 (Attended - Unattended)"
)
#------------------------------- Two-way Interactions -------------------------------------#



#------------------------------------------------------------------------------------------#
#--------------------------------- Kappa Scale --------------------------------------------#
#------------------------------------------------------------------------------------------#
 
#---------------------------------- Kappa Intercept ---------------------------------------#
get_violin(
  exp(ex_toj_color_post$population_log_kappa_intercept_mean)
  , "Fidelity of Memory Intercept Mean"
  , y_lab = "\u03BA"
  , hline = FALSE
)
#---------------------------------- Kappa Intercept ---------------------------------------#


#-------------------------------- Main Effects --------------------------------------------#
get_violin(
  ( exp(ex_toj_color_post$population_log_kappa_intercept_mean + ex_toj_color_post$population_log_kappa_attention_effect_mean/2 )
    - exp(ex_toj_color_post$population_log_kappa_intercept_mean - ex_toj_color_post$population_log_kappa_attention_effect_mean/2 ) )
  , "Fidelity of Encoding Attention Effect Mean"
  , y_lab = "\u03BA (Attended - Unattended)"
)

get_violin(
  ( exp(ex_toj_color_post$population_log_kappa_intercept_mean + ex_toj_color_post$population_log_kappa_probe_duration_effect_mean/2 )
    - exp(ex_toj_color_post$population_log_kappa_intercept_mean - ex_toj_color_post$population_log_kappa_probe_duration_effect_mean/2 ) )
  , "Fidelity of Encoding Probe Duration Effect Mean"
  , y_lab = "\u03BA (Long - Short)"
)

# get_violin(
#   ( exp(ex_toj_color_post$population_log_kappa_intercept_mean + ex_toj_color_post$logKappaJudgementTypeEffectMean/2 )
#     - exp(ex_toj_color_post$population_log_kappa_intercept_mean - ex_toj_color_post$logKappaJudgementTypeEffectMean/2 ) )
#   , "Fidelity of Memory\nJudgement Type Effect Mean"
#   , y_lab = "\u03BA (Second - First)"
# )
# 
# get_violin(
#   ( exp(ex_toj_color_post$population_log_kappa_intercept_mean + ex_toj_color_post$logKappaInitialBiasEffectMean/2 )
#     - exp(ex_toj_color_post$population_log_kappa_intercept_mean - ex_toj_color_post$logKappaInitialBiasEffectMean/2 ) )
#   , "Fidelity of Memory\nInitial Probe Bias Effect Mean"
#   , y_lab = "\u03BA (Left - Right)"
# )
#-------------------------------- Main Effects --------------------------------------------#


#------------------------------- Two-way Interactions -------------------------------------#
get_violin(
  2*  # MULTIPLY INTERACTION BY 2 TO GET DIFFERENCE IN DIFFERENCES (Meaningful)
  ( exp(ex_toj_color_post$population_log_kappa_intercept_mean + ex_toj_color_post$population_log_kappa_attention_probe_duration_interaction_effect_mean/2 )
    - exp(ex_toj_color_post$population_log_kappa_intercept_mean - ex_toj_color_post$population_log_kappa_attention_probe_duration_interaction_effect_mean/2 ) )
  , "Fidelity of Encoding\nProbe Duration Attention\nInteraction Effect Mean"
  , y_lab = "\u03BA"
)

get_violin(
  c(
    exp(ex_toj_color_post$population_log_kappa_intercept_mean + ex_toj_color_post$population_log_kappa_probe_duration_effect_mean/2 
           + (ex_toj_color_post$population_log_kappa_attention_effect_mean + ex_toj_color_post$population_log_kappa_attention_probe_duration_interaction_effect_mean )/2 )
    - exp(ex_toj_color_post$population_log_kappa_intercept_mean + ex_toj_color_post$population_log_kappa_probe_duration_effect_mean/2 
             - (ex_toj_color_post$population_log_kappa_attention_effect_mean + ex_toj_color_post$population_log_kappa_attention_probe_duration_interaction_effect_mean )/2 )
    , exp(ex_toj_color_post$population_log_kappa_intercept_mean - ex_toj_color_post$population_log_kappa_probe_duration_effect_mean/2 
             + (ex_toj_color_post$population_log_kappa_attention_effect_mean - ex_toj_color_post$population_log_kappa_attention_probe_duration_interaction_effect_mean )/2 )
    - exp(ex_toj_color_post$population_log_kappa_intercept_mean - ex_toj_color_post$population_log_kappa_probe_duration_effect_mean/2 
             - (ex_toj_color_post$population_log_kappa_attention_effect_mean - ex_toj_color_post$population_log_kappa_attention_probe_duration_interaction_effect_mean )/2 ) 
  )
  , c("Fidelity of Encoding\nAttention Effect \nGiven Long Probe Duration","Fidelity of Encoding\nAttention Effect\nGiven Short Probe Duration")
  , y_lab = "\u03BA"
)
#------------------------------- Two-way Interactions -------------------------------------#





#------------------------------------------------------------------------------------------#
#--------------------------------- Graphs -------------------------------------------------#
#------------------------------------------------------------------------------------------#

#---------------------------------- NCFs --------------------------------------------------#
# including effects 
yLeftFirst = pnorm(
  -250:250
  , mean = ( 
    median(ex_toj_color_post$population_pss_intercept_mean) - median(ex_toj_color_post$population_pss_judgement_type_effect_mean)/2 
   - ( 
     median(ex_toj_color_post$population_pss_attention_effect_mean) 
     - median(ex_toj_color_post$population_pss_attention_judgement_type_interaction_effect_mean)
    ) /2 
   ) * 250
  , sd = ( exp( 
    median(ex_toj_color_post$population_log_jnd_intercept_mean) - median(ex_toj_color_post$population_log_jnd_judgement_type_effect_mean)/2 
    - ( 
      median(ex_toj_color_post$population_log_jnd_attention_effect_mean) 
      - median(ex_toj_color_post$population_log_jnd_attention_judgement_type_interaction_effect_mean)
    ) /2 
    ) ) * 250
)

yLeftSecond = pnorm(
  -250:250
  , mean = ( 
    median(ex_toj_color_post$population_pss_intercept_mean) + median(ex_toj_color_post$population_pss_judgement_type_effect_mean)/2 
    - ( 
      median(ex_toj_color_post$population_pss_attention_effect_mean) 
      + median(ex_toj_color_post$population_pss_attention_judgement_type_interaction_effect_mean)
    ) /2 
  ) * 250
  , sd = ( exp( 
    median(ex_toj_color_post$population_log_jnd_intercept_mean) + median(ex_toj_color_post$population_log_jnd_judgement_type_effect_mean)/2 
    - ( 
      median(ex_toj_color_post$population_log_jnd_attention_effect_mean) 
      + median(ex_toj_color_post$population_log_jnd_attention_judgement_type_interaction_effect_mean)
    ) /2 
  ) ) * 250
)

yRightFirst = pnorm(
  -250:250
  , mean = ( 
    median(ex_toj_color_post$population_pss_intercept_mean) - median(ex_toj_color_post$population_pss_judgement_type_effect_mean)/2 
    + ( 
      median(ex_toj_color_post$population_pss_attention_effect_mean) 
      - median(ex_toj_color_post$population_pss_attention_judgement_type_interaction_effect_mean)
    ) /2 
  ) * 250
  , sd = ( exp( 
    median(ex_toj_color_post$population_log_jnd_intercept_mean) - median(ex_toj_color_post$population_log_jnd_judgement_type_effect_mean)/2 
    + ( 
      median(ex_toj_color_post$population_log_jnd_attention_effect_mean) 
      - median(ex_toj_color_post$population_log_jnd_attention_judgement_type_interaction_effect_mean)
    ) /2 
  ) ) * 250
)

yRightSecond = pnorm(
  -250:250
  , mean = ( 
    median(ex_toj_color_post$population_pss_intercept_mean) + median(ex_toj_color_post$population_pss_judgement_type_effect_mean)/2 
    + ( 
      median(ex_toj_color_post$population_pss_attention_effect_mean) 
      + median(ex_toj_color_post$population_pss_attention_judgement_type_interaction_effect_mean)
    ) /2 
  ) * 250
  , sd = ( exp( 
    median(ex_toj_color_post$population_log_jnd_intercept_mean) - median(ex_toj_color_post$population_log_jnd_judgement_type_effect_mean)/2 
    + ( 
      median(ex_toj_color_post$population_log_jnd_attention_effect_mean) 
      - median(ex_toj_color_post$population_log_jnd_attention_judgement_type_interaction_effect_mean)
    ) /2 
  ) ) * 250
)

df = data.frame(SOA = -250:250
                , Prop = c(yRightFirst, yRightSecond, yLeftFirst, yLeftSecond)
                , Attend = c(rep("Right",1002), rep("Left", 1002))
                , Judgement = c(rep( c(rep("Which First", 501), rep("Which Second", 501)), 2) )
                )

toj_means_by_id_by_condition2 = ddply(
  .data = toj_trials
  , .variables = .(block_bias, toj_judgement_type, soa2)
  , .fun = function(x){
    to_return = data.frame(
      Prop = mean(x$left_first_TF)
      , Judgement = paste(unique(x$toj_judgement_type))
      , Attend = paste( unique(x$block_bias) )
      , SOA = unique(x$soa2)
    )
    return(to_return)
  }
)

levels(toj_means_by_id_by_condition2$Judgement) = c("Which First", "Which Second")
levels(toj_means_by_id_by_condition2$Attend) = c("Left", "Right")

gg = ggplot(data = df, aes(y = Prop, x = SOA, colour = Attend))+
  geom_line(size = 1.25)+
  # scale_color_manual("Attend", values = c("red", "blue"))+
  scale_color_hue("Attend", l = c(60, 15), c = c(100, 50), h = c(240, 360) ) +
  labs(x = "SOA (ms)", y = "Proportion of 'Left First' Responses")+
  geom_point(data = toj_means_by_id_by_condition2, aes(y = Prop, x = SOA, colour = Attend), size = 4)+
  geom_point(data = toj_means_by_id_by_condition2, aes(y = Prop, x = SOA), size = 2.5, colour = "grey90")+
  facet_grid(Judgement~.)+
  theme_gray(base_size = 30)+
  theme(panel.grid.major = element_line(size = 1.5)
        ,panel.grid.minor = element_line(size = 1))
# define text to add
Text1 = textGrob(label = paste("Right First"), gp = gpar(fontsize= 30))
Text2 = textGrob(label = paste("Left First"), gp = gpar(fontsize= 30)) 
gg = gg+
  annotation_custom(grob = Text1,  xmin = -200, xmax = -200, ymin = -0.25, ymax = -0.25)+
  annotation_custom(grob = Text2,  xmin = 200, xmax = 200, ymin = -0.25, ymax = -0.25)
# Code to override clipping
gg2 <- ggplot_gtable(ggplot_build(gg))
gg2$layout$clip[gg2$layout$name=="panel"] <- "off"
grid.draw(gg2)
#---------------------------------- NCFs --------------------------------------------------#






