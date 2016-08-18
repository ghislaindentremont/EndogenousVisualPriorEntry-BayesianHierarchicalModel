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
               , "population_pss_judgement_type_effect_mean"
               , "population_pss_probe_duration_effect_mean"
               , "population_pss_attention_judgement_type_interaction_effect_mean"	
               , "population_pss_attention_probe_duration_interaction_effect_mean"
               , "population_log_jnd_intercept_mean"
               , "population_log_jnd_attention_effect_mean"
               , "population_log_jnd_judgement_type_effect_mean"
               , "population_log_jnd_probe_duration_effect_mean"
               , "population_log_jnd_attention_judgement_type_interaction_effect_mean"	
               , "population_log_jnd_attention_probe_duration_interaction_effect_mean"
               , "population_logit_rho_intercept_mean"
               , "population_logit_rho_attention_effect_mean"
               , "population_logit_rho_probe_duration_effect_mean"
               , "population_logit_rho_attention_probe_duration_interaction_effect_mean"
               , "population_log_kappa_intercept_mean"
               , "population_log_kappa_attention_effect_mean"
               , "population_log_kappa_probe_duration_effect_mean"
               , "population_log_kappa_attention_probe_duration_interaction_effect_mean"
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

# judgement type
log_jnd_right_second_long_mean_reps = get_condition_mean_sample(
  ( log_jnd_intercept_mean 
  + log_jnd_judgement_type_effect_mean/2
  + log_jnd_probe_duration_effect_mean/2 )
  , ( log_jnd_effect_mean 
     + log_jnd_judgement_type_interaction_effect_mean
     + log_jnd_probe_duration_interaction_effect_mean )
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

rho_probe_duration_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "population_logit_rho_probe_duration_effect_mean",]$value
rho_probe_duration_interaction_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "population_logit_rho_attention_probe_duration_interaction_effect_mean",]$value

# get condition combination
rho_right_long_mean_reps = get_condition_mean_sample(
  ( rho_intercept_mean 
    + rho_probe_duration_effect_mean/2 )
  , ( rho_effect_mean 
      + rho_probe_duration_interaction_effect_mean )
  , TRUE
  , "logit"
)

### Get Kappa Parameters
kappa_intercept_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "population_log_kappa_intercept_mean",]$value
kappa_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "population_log_kappa_attention_effect_mean",]$value

kappa_probe_duration_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "population_log_kappa_probe_duration_effect_mean",]$value
kappa_probe_duration_interaction_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "population_log_kappa_attention_probe_duration_interaction_effect_mean",]$value

# get condition combination
kappa_right_long_mean_reps = get_condition_mean_sample(
  ( kappa_intercept_mean 
    + kappa_probe_duration_effect_mean/2)
  , ( kappa_effect_mean 
        + kappa_probe_duration_interaction_effect_mean )
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

# (1) population_pss_intercept_mean     *
# (2) population_pss_attention_effect_mean          
# (3) population_log_jnd_intercept_mean * ^
# (4) population_log_jnd_attention_effect_mean     
# (5) population_logit_rho_intercept_mean ^                  
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


# #---------------------------- Rho vs. PSS Effects -----------------------------------------#
# psseffect = extract_samples("population_pss_attention_effect_mean")
# 
# psseffectsd = extract_samples("population_pss_effect_sd")
# 
# pssjudgementinteraction = extract_samples("population_pss_attention_judgement_type_interaction_effect_mean")
# 
# pssinitialbiasinteraction = extract_samples("population_pss_attention_initial_bias_interaction_effect_mean")
# 
# pssprobeinteraction = extract_samples("population_pss_attention_probe_duration_interaction_effect_mean")
# 
# judgementfactor = ifelse(aggregate(toj_judgement_type ~ id, data = toj_trials, FUN = unique)$toj_judgement_type == "first", -1, 1)
# 
# probefactor = ifelse(aggregate(longprobe ~ id, data = color_trials, FUN = unique)$longprobe == "FALSE", -1, 1)
# 
# psseffect_ids = ddply(
#   .data = betas
#   , .variables = .(participant)
#   , .fun = function(x){
#     i = unique(x$participant)
#     x_use = x[x$parameter ==  "population_pss_attention_effect_mean",]$value
#     psseffect_use = median(psseffect)  + median(psseffectsd)*median(x_use)  + median(pssjudgementinteraction)*judgementfactor[i]+ median(pssprobeinteraction)*probefactor[i]
#     df = data.frame(psseffect_use*250, judgementfactor[i], probefactor[i])
#     names(df) = c("psseffect",  "judgementfactor", "probefactor")
#     return(df)
#   }
# )
# 
# logitrhoeffect = extract_samples("population_logit_rho_attention_effect_mean")
# 
# logitrhoeffectsd = extract_samples("population_logit_rho_effect_sd")
# 
# # logitrhojudgementinteractioneffect = extract_samples("population_logit_rho_attention_judgement_type_interaction_effect_mean")
# 
# logitrhoprobeinteractioneffect = extract_samples("population_logit_rho_attention_probe_duration_interaction_effect_mean")
# 
# rhoeffect_ids = ddply(
#   .data = betas
#   , .variables = .(participant)
#   , .fun = function(x){
#     i = unique(x$participant)
#     x_use = x[x$parameter == "population_logit_rho_attention_effect_mean",]$value
#     logitrhoeffect_use =  median(logitrhoeffect) + median(logitrhoeffectsd)*median(x_use)+ median(logitrhoprobeinteractioneffect)*probefactor[i] #+ median(logitrhojudgementinteractioneffect)*judgementfactor[i]
#     df = data.frame(logitrhoeffect_use, probefactor[i])
#     names(df) = c("logitrhoeffect", "probefactor")
#     return(df)
#   }
# )
# 
# psseffect_v_rhoeffect = merge(rhoeffect_ids, psseffect_ids)
# 
# ggplot(data = psseffect_v_rhoeffect, aes(y = psseffect, x = logitrhoeffect, colour = factor(judgementfactor), shape = factor(judgementfactor), fill = factor(probefactor)))+
#   scale_y_continuous(name = "PSS Effect Mean")+
#   scale_x_continuous(name = "Logit \u03C1 Effect Mean")+
#   geom_vline(xintercept = 0, linetype = 2, size = 1)+
#   geom_hline(yintercept = 0, linetype = 2, size = 1)+
#   geom_point(size = 4)+
#   scale_shape_manual(name = "Judgement\nType", labels = c("First", "Second") , values = c(21,22) )+
#   scale_colour_manual(name = "Judgement\nType", labels =c("First", "Second"), values = c("red", "blue") )+
#   scale_fill_manual(name = "Probe\nDuration", labels = c("Short", "Long"), values = c("black", "white")) +
#   theme_gray(base_size = 30)+
#   theme(panel.grid.major = element_line(size = 1.5)
#         ,panel.grid.minor = element_line(size = 1))
# 
# # NOTE: RECALL THAT COR BELOW IS FOR ATTENTION EFFECT ONLY
# # get_corr(
# #   "value.2.7"
# #   , "Logit \u03C1 vs. PSS Effect Means"
# # )
# #---------------------------- Rho vs. PSS Effects -----------------------------------------#
# 
# 
# #---------------------------- Kappa vs. PSS Effects ---------------------------------------#
# logkappaeffect = extract_samples("population_log_kappa_attention_effect_mean")
# 
# logkappaeffectsd = extract_samples("population_log_kappa_effect_sd")
# 
# # logkappajudgementinteractioneffect = extract_samples("population_log_kappa_attention_judgement_type_interaction_effect_mean")
# 
# logkappaprobeinteractioneffect = extract_samples("population_log_kappa_attention_probe_duration_interaction_effect_mean")
# 
# kappaeffect_ids = ddply(
#   .data = betas
#   , .variables = .(participant)
#   , .fun = function(x){
#     i = unique(x$participant)
#     x_use = x[x$parameter == "population_log_kappa_attention_effect_mean",]$value
#     logkappaeffect_use =  median(logkappaeffect) + median(logkappaeffectsd)*median(x_use) + median(logkappaprobeinteractioneffect)*probefactor[i] # + median(logkappajudgementinteractioneffect)*judgementfactor[i]
#     df = data.frame(logkappaeffect_use, probefactor[i])
#     names(df) = c("kappaeffect", "probefactor")
#     return(df)
#   }
# )
# 
# psseffect_v_kappaeffect = merge(kappaeffect_ids, psseffect_ids)
# 
# ggplot(data = psseffect_v_kappaeffect, aes(y = psseffect, x = kappaeffect, colour = factor(judgementfactor), shape = factor(judgementfactor), fill = factor(probefactor)))+ 
#   scale_y_continuous(name = "PSS Effect Mean")+
#   scale_x_continuous(name = "Log \u03BA Effect Mean")+
#   geom_vline(xintercept = 0, linetype = 2, size = 1)+
#   geom_hline(yintercept = 0, linetype = 2, size = 1)+
#   geom_point(size = 4)+
#   scale_shape_manual(name = "Judgement\nType", labels = c("First", "Second") , values = c(21,22) )+
#   scale_colour_manual(name = "Judgement\nType", labels =c("First", "Second"), values = c("red", "blue") )+
#   scale_fill_manual(name = "Probe\nDuration", labels = c("Short", "Long"), values = c("black", "white")) +
#   theme_gray(base_size = 30)+
#   theme(panel.grid.major = element_line(size = 1.5)
#         ,panel.grid.minor = element_line(size = 1))
# 
# # get_corr(
# #   "value.2.8"
# #   , "Log \u03BA vs. PSS Effect Means"
# # )        
# #---------------------------- Kappa vs. PSS Effects ---------------------------------------#



#------------------------------------------------------------------------------------------#
#--------------------------------- PSS Scale ----------------------------------------------#
#------------------------------------------------------------------------------------------#

#---------------------------------- PSS Intercept -----------------------------------------#
get_violin(
    ex_toj_color_post$population_pss_intercept_mean * 250
  , "PSS Intercept Mean"
  , y_lab = "SOA (ms)"
  , hline = TRUE
)
#---------------------------------- PSS Intercept -----------------------------------------#


#---------------------------------- Main Effects ------------------------------------------#
# effect of attention on PSS
get_violin(
  ( ex_toj_color_post$population_pss_attention_effect_mean ) * 250
  , c("PSS\nAttention Effect Mean")
  , y_lab = "SOA (Right - Left; ms)"
)

# effect of judgement type (Q) on PSS 
get_violin(
  ( ex_toj_color_post$population_pss_judgement_type_effect_mean ) * 250 
  , "PSS Judgement\nType Effect Mean"
  , y_lab = "SOA (Second - First; ms)"
)

# effect of probe duration on PSS 
get_violin(
  ( ex_toj_color_post$population_pss_probe_duration_effect_mean ) * 250 
  , "PSS Probe\nDuration Effect Mean" 
  , y_lab = "SOA (Long - Short; ms)"
)
#---------------------------------- Main Effects ---------------------=--------------------#


#------------------------------- Two-way Interactions -------------------------------------#
#  effect of interaction between judgement type and attention on PSS 
get_violin(
  # MULTIPLY BY TWO because not properly calculated in .stan 
  2 * ( ex_toj_color_post$population_pss_attention_judgement_type_interaction_effect_mean) * 250
  , "PSS Attention\n& Judgement Type\nInteraction Effect Mean"
  , y_lab = "SOA (ms)"
)

# NOTE: probe duration effects will cancel out (because no back-transformation)
get_violin(
  c(
  (ex_toj_color_post$population_pss_intercept_mean + ex_toj_color_post$population_pss_judgement_type_effect_mean/2
    + (ex_toj_color_post$population_pss_attention_effect_mean + ex_toj_color_post$population_pss_attention_judgement_type_interaction_effect_mean)/2 ) *250
  - (ex_toj_color_post$population_pss_intercept_mean + ex_toj_color_post$population_pss_judgement_type_effect_mean/2
    - (ex_toj_color_post$population_pss_attention_effect_mean + ex_toj_color_post$population_pss_attention_judgement_type_interaction_effect_mean)/2 ) * 250
  , (ex_toj_color_post$population_pss_intercept_mean - ex_toj_color_post$population_pss_judgement_type_effect_mean/2
       + (ex_toj_color_post$population_pss_attention_effect_mean - ex_toj_color_post$population_pss_attention_judgement_type_interaction_effect_mean)/2 ) *250
  - (ex_toj_color_post$population_pss_intercept_mean - ex_toj_color_post$population_pss_judgement_type_effect_mean/2
     - (ex_toj_color_post$population_pss_attention_effect_mean - ex_toj_color_post$population_pss_attention_judgement_type_interaction_effect_mean)/2 ) *250
  )
  , c("PSS Attention Effect\nGiven Which Second", "PSS Attention Effect\nGiven Which First")
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
  # MULTIPLY BY TWO because not properly calculated in .stan
  2 * (ex_toj_color_post$population_pss_attention_probe_duration_interaction_effect_mean) * 250
  , "PSS Attention\n& Probe Duration\nInteraction Effect Mean"
  , y_lab = "SOA (ms)"
)

# NOTE: judgement_type effects will cancel out (because no back-transformation)
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
  , c("PSS Attention Effect\nGiven Long\nProbe Duration", "PSS Attention Effect\nGiven Short\nProbe Duration")
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
#--------------------------------- JND Scale ----------------------------------------------#
#------------------------------------------------------------------------------------------#


#---------------------------------- JND Conditions ----------------------------------------#
pos_attn_pos_judge_pos_probe = 
  exp(ex_toj_color_post$population_log_jnd_intercept_mean + ex_toj_color_post$population_log_jnd_judgement_type_effect_mean/2 + ex_toj_color_post$population_log_jnd_probe_duration_effect_mean/2 
      + (ex_toj_color_post$population_log_jnd_attention_effect_mean + ex_toj_color_post$population_log_jnd_attention_judgement_type_interaction_effect_mean + ex_toj_color_post$population_log_jnd_attention_probe_duration_interaction_effect_mean)/2 ) 

pos_attn_pos_judge_neg_probe = 
  exp(ex_toj_color_post$population_log_jnd_intercept_mean + ex_toj_color_post$population_log_jnd_judgement_type_effect_mean/2 - ex_toj_color_post$population_log_jnd_probe_duration_effect_mean/2 
      + (ex_toj_color_post$population_log_jnd_attention_effect_mean + ex_toj_color_post$population_log_jnd_attention_judgement_type_interaction_effect_mean - ex_toj_color_post$population_log_jnd_attention_probe_duration_interaction_effect_mean)/2 ) 

pos_attn_neg_judge_neg_probe = 
  exp(ex_toj_color_post$population_log_jnd_intercept_mean - ex_toj_color_post$population_log_jnd_judgement_type_effect_mean/2 - ex_toj_color_post$population_log_jnd_probe_duration_effect_mean/2 
      + (ex_toj_color_post$population_log_jnd_attention_effect_mean - ex_toj_color_post$population_log_jnd_attention_judgement_type_interaction_effect_mean - ex_toj_color_post$population_log_jnd_attention_probe_duration_interaction_effect_mean)/2 ) 

pos_attn_neg_judge_pos_probe = 
  exp(ex_toj_color_post$population_log_jnd_intercept_mean - ex_toj_color_post$population_log_jnd_judgement_type_effect_mean/2 + ex_toj_color_post$population_log_jnd_probe_duration_effect_mean/2 
      + (ex_toj_color_post$population_log_jnd_attention_effect_mean - ex_toj_color_post$population_log_jnd_attention_judgement_type_interaction_effect_mean + ex_toj_color_post$population_log_jnd_attention_probe_duration_interaction_effect_mean)/2 ) 

neg_attn_pos_judge_pos_probe =
  exp(ex_toj_color_post$population_log_jnd_intercept_mean + ex_toj_color_post$population_log_jnd_judgement_type_effect_mean/2 + ex_toj_color_post$population_log_jnd_probe_duration_effect_mean/2 
      - (ex_toj_color_post$population_log_jnd_attention_effect_mean + ex_toj_color_post$population_log_jnd_attention_judgement_type_interaction_effect_mean + ex_toj_color_post$population_log_jnd_attention_probe_duration_interaction_effect_mean)/2 ) 

neg_attn_pos_judge_neg_probe = 
  exp(ex_toj_color_post$population_log_jnd_intercept_mean + ex_toj_color_post$population_log_jnd_judgement_type_effect_mean/2 - ex_toj_color_post$population_log_jnd_probe_duration_effect_mean/2 
      - (ex_toj_color_post$population_log_jnd_attention_effect_mean + ex_toj_color_post$population_log_jnd_attention_judgement_type_interaction_effect_mean - ex_toj_color_post$population_log_jnd_attention_probe_duration_interaction_effect_mean)/2 ) 

neg_attn_neg_judge_neg_probe =
  exp(ex_toj_color_post$population_log_jnd_intercept_mean - ex_toj_color_post$population_log_jnd_judgement_type_effect_mean/2 - ex_toj_color_post$population_log_jnd_probe_duration_effect_mean/2 
      - (ex_toj_color_post$population_log_jnd_attention_effect_mean - ex_toj_color_post$population_log_jnd_attention_judgement_type_interaction_effect_mean - ex_toj_color_post$population_log_jnd_attention_probe_duration_interaction_effect_mean)/2 ) 

neg_attn_neg_judge_pos_probe =
  exp(ex_toj_color_post$population_log_jnd_intercept_mean - ex_toj_color_post$population_log_jnd_judgement_type_effect_mean/2 + ex_toj_color_post$population_log_jnd_probe_duration_effect_mean/2 
      - (ex_toj_color_post$population_log_jnd_attention_effect_mean - ex_toj_color_post$population_log_jnd_attention_judgement_type_interaction_effect_mean + ex_toj_color_post$population_log_jnd_attention_probe_duration_interaction_effect_mean)/2 ) 

#---------------------------------- JND Conditions ----------------------------------------#


#---------------------------------- JND Intercept -----------------------------------------#
### BAD:
get_violin(
  exp( ex_toj_color_post$population_log_jnd_intercept_mean ) * 250
  , "JND Intercept Mean"
  , y_lab = "SOA (ms)"
  , hline = FALSE
)

# GOOD:
get_violin(
   ( pos_attn_neg_judge_pos_probe
   + pos_attn_neg_judge_neg_probe
   + pos_attn_pos_judge_neg_probe
   + pos_attn_pos_judge_pos_probe
   + neg_attn_neg_judge_pos_probe
   + neg_attn_neg_judge_neg_probe
   + neg_attn_pos_judge_neg_probe
   + neg_attn_pos_judge_pos_probe
   )/8 * 250
  , "JND Intercept Mean"
  , y_lab = "SOA (ms)"
  , hline = FALSE
)

# BETTER:
get_violin(
  ex_toj_color_post$population_log_jnd_intercept_mean
  , "JND Intercept Mean"
  , y_lab = "Log of SOA (ms)"
  , hline = FALSE
)
#---------------------------------- JND Intercept -----------------------------------------#


#---------------------------------- Main Effects ------------------------------------------#
### BAD: effect of attention on JND
get_violin(
  ( exp( ex_toj_color_post$population_log_jnd_intercept_mean + ex_toj_color_post$population_log_jnd_attention_effect_mean/2 )
    - exp( ex_toj_color_post$population_log_jnd_intercept_mean - ex_toj_color_post$population_log_jnd_attention_effect_mean/2  ) ) * 250 
  , "JND\nAttention Effect Mean"
  , y_lab = "SOA (Right - Left; ms)"
)

# GOOD: effect of attention on JND
get_violin(
  ( 
    ( pos_attn_neg_judge_pos_probe
      + pos_attn_neg_judge_neg_probe
      + pos_attn_pos_judge_neg_probe
      + pos_attn_pos_judge_pos_probe
    )/4 
    - 
      ( neg_attn_neg_judge_pos_probe
        + neg_attn_neg_judge_neg_probe
        + neg_attn_pos_judge_neg_probe
        + neg_attn_pos_judge_pos_probe
      )/4 
  ) *250
  , "JND\nAttention Effect Mean"
  , y_lab = "SOA (Right - Left; ms)"
)

# BETTER: effect of attention on JND
get_violin(
  ex_toj_color_post$population_log_jnd_attention_effect_mean
  , "JND\nAttention Effect Mean"
  , y_lab = "Log of SOA (Right - Left; ms)"
)


### BAD: effect of judgement type (Q) on JND
get_violin(
  ( exp( ex_toj_color_post$population_log_jnd_intercept_mean + ex_toj_color_post$population_log_jnd_judgement_type_effect_mean/2 )
    - exp( ex_toj_color_post$population_log_jnd_intercept_mean - ex_toj_color_post$population_log_jnd_judgement_type_effect_mean/2  ) ) * 250
  , "JND Judgement\nType Effect Mean"
  , y_lab = "SOA (Second - First; ms)"
)

# GOOD: effect of judgement type (Q) on JND
get_violin(
  ( 
    ( pos_attn_pos_judge_neg_probe
      + pos_attn_pos_judge_pos_probe
      + neg_attn_pos_judge_neg_probe
      + neg_attn_pos_judge_pos_probe
    )/4 
    - 
      ( pos_attn_neg_judge_neg_probe
        + pos_attn_neg_judge_pos_probe
        + neg_attn_neg_judge_neg_probe
        + neg_attn_neg_judge_pos_probe
      )/4 
  ) *250  
  , "JND Judgement\nType Effect Mean"
  , y_lab = "SOA (Second - First; ms)"
)

# BETTER: effect of judgement type (Q) on JND
get_violin(
  ex_toj_color_post$population_log_jnd_judgement_type_effect_mean
  , "JND Judgement\nType Effect Mean"
  , y_lab = "Log of SOA (Second - First; ms)"
)


### BAD: effect of probe duration on JND
get_violin(
  ( exp( ex_toj_color_post$population_log_jnd_intercept_mean + ex_toj_color_post$population_log_jnd_probe_duration_effect_mean/2 )
    - exp( ex_toj_color_post$population_log_jnd_intercept_mean - ex_toj_color_post$population_log_jnd_probe_duration_effect_mean/2  ) ) * 250
  , "JND Probe\nDuration Effect Mean"
  , y_lab = "SOA (Long - Short; ms)"
)

# GOOD: effect of probe duration on JND
get_violin(
  ( 
    ( pos_attn_pos_judge_pos_probe
      + neg_attn_neg_judge_pos_probe
      + neg_attn_pos_judge_pos_probe
      + pos_attn_neg_judge_pos_probe
    )/4 
    - 
      ( pos_attn_pos_judge_neg_probe
        + pos_attn_neg_judge_neg_probe
        + neg_attn_neg_judge_neg_probe
        + neg_attn_pos_judge_neg_probe
      )/4 
  ) *250  
  , "JND Probe\nDuration Effect Mean"
  , y_lab = "SOA (Long - Short; ms)"
)

# BETTER: effect of probe duration on JND
get_violin(
  ex_toj_color_post$population_log_jnd_probe_duration_effect_mean
  , "JND Probe\nDuration Effect Mean"
  , y_lab = "Log of SOA (Long - Short; ms)"
)
#---------------------------------- Main Effects ---------------------=--------------------#


#------------------------------- Two-way Interactions -------------------------------------#
###  BAD: effect of interaction between judgement type and attention on JND 
get_violin(
  # don't devide effects by two...
  ( 
    exp( ex_toj_color_post$population_log_jnd_intercept_mean + ex_toj_color_post$population_log_jnd_attention_judgement_type_interaction_effect_mean)
    - exp(  ex_toj_color_post$population_log_jnd_intercept_mean - ex_toj_color_post$population_log_jnd_attention_judgement_type_interaction_effect_mean)
  ) *250 
  , "JND Attention\n& Judgement Type\nInteraction Effect Mean"
  , y_lab = "SOA (ms)"
)

#  GOOD: effect of interaction between judgement type and attention on JND 
get_violin(
  ( 
    # average attention effect for positive judgement across probe conditions 
    ( pos_attn_pos_judge_pos_probe - neg_attn_pos_judge_pos_probe
      + pos_attn_pos_judge_neg_probe - neg_attn_pos_judge_neg_probe
    )/2 
    -
      # average attention effect for negative judgement across probe conditions 
      ( pos_attn_neg_judge_pos_probe - neg_attn_neg_judge_pos_probe
        + pos_attn_neg_judge_neg_probe - neg_attn_neg_judge_neg_probe
      )/2
  ) *250 
  , "JND Attention\n& Judgement Type\nInteraction Effect Mean"
  , y_lab = "SOA (ms)"
)

#  BETTER: effect of interaction between judgement type and attention on JND 
get_violin(
  2 * ex_toj_color_post$population_log_jnd_attention_judgement_type_interaction_effect_mean
  , "JND Attention\n& Judgement Type\nInteraction Effect Mean"
  , y_lab = "Log of SOA (ms)"  # awkward to label difference in difference direction in parentheses 
  # so explain in figure caption
)


###  BAD: effect of interaction between probe duration and attention on JND 
get_violin(
  # don't devide effects by two...
  ( 
    exp( ex_toj_color_post$population_log_jnd_intercept_mean + ex_toj_color_post$population_log_jnd_attention_probe_duration_interaction_effect_mean)
    - exp(  ex_toj_color_post$population_log_jnd_intercept_mean - ex_toj_color_post$population_log_jnd_attention_probe_duration_interaction_effect_mean)
  ) *250 
  , "JND Attention\n& Probe Duration\nInteraction Effect Mean"
  , y_lab = "SOA (ms)"
)

#  GOOD: effect of interaction between probe duration and attention on JND
get_violin(
  ( 
    # average attention effect for positive probe across judgement conditions 
    ( pos_attn_pos_judge_pos_probe - neg_attn_pos_judge_pos_probe
      + pos_attn_neg_judge_pos_probe - neg_attn_neg_judge_pos_probe 
    )/2 
    -
      # average attention effect for negative probe across judgement conditions 
      ( pos_attn_pos_judge_neg_probe - neg_attn_pos_judge_neg_probe
        + pos_attn_neg_judge_neg_probe - neg_attn_neg_judge_neg_probe
      )/2
  ) *250 
  , "JND Attention\n& Probe Duration\nInteraction Effect Mean"
  , y_lab = "SOA (ms)"
)

#  BETTER: effect of interaction between probe duration and attention on JND 
get_violin(
  2 * ex_toj_color_post$population_log_jnd_attention_probe_duration_interaction_effect_mean
  , "JND Attention\n& Probe Duration\nInteraction Effect Mean"
  , y_lab = "Log of SOA (ms)"
)
#------------------------------- Two-way Interactions -------------------------------------#



#------------------------------------------------------------------------------------------#
#--------------------------------- Rho Scale ----------------------------------------------#
#------------------------------------------------------------------------------------------#

#--------------------------------- Rho Conditions -----------------------------------------#
rho_pos_attn_pos_probe = 
  plogis(
    ex_toj_color_post$population_logit_rho_intercept_mean + ex_toj_color_post$population_logit_rho_probe_duration_effect_mean/2
    + (ex_toj_color_post$population_logit_rho_attention_effect_mean + ex_toj_color_post$population_logit_rho_attention_probe_duration_interaction_effect_mean)/2
  )

rho_pos_attn_neg_probe = 
  plogis(
    ex_toj_color_post$population_logit_rho_intercept_mean - ex_toj_color_post$population_logit_rho_probe_duration_effect_mean/2
    + (ex_toj_color_post$population_logit_rho_attention_effect_mean - ex_toj_color_post$population_logit_rho_attention_probe_duration_interaction_effect_mean)/2
  )

rho_neg_attn_neg_probe = 
  plogis(
    ex_toj_color_post$population_logit_rho_intercept_mean - ex_toj_color_post$population_logit_rho_probe_duration_effect_mean/2
    - (ex_toj_color_post$population_logit_rho_attention_effect_mean - ex_toj_color_post$population_logit_rho_attention_probe_duration_interaction_effect_mean)/2
  )

rho_neg_attn_pos_probe = 
  plogis(
    ex_toj_color_post$population_logit_rho_intercept_mean + ex_toj_color_post$population_logit_rho_probe_duration_effect_mean/2
    - (ex_toj_color_post$population_logit_rho_attention_effect_mean + ex_toj_color_post$population_logit_rho_attention_probe_duration_interaction_effect_mean)/2
  )
#--------------------------------- Rho Conditions -----------------------------------------#


#---------------------------------- Rho Intercept -----------------------------------------#
### BAD:
get_violin(
  plogis(ex_toj_color_post$population_logit_rho_intercept_mean)
  , "Probability of Encoding Intercept Mean"
  , y_lab = "\u03C1"
  , hline = FALSE
)

# GOOD:
get_violin(
  ( rho_neg_attn_pos_probe
    + rho_neg_attn_neg_probe
    + rho_pos_attn_neg_probe
    + rho_pos_attn_pos_probe
  )/4
  , "Probability of Encoding Intercept Mean"
  , y_lab = "\u03C1"
  , hline = FALSE
)

# BETTER:
get_violin(
  ex_toj_color_post$population_logit_rho_intercept_mean
  , "Probability of Encoding Intercept Mean"
  , y_lab = "Log-odds of \u03C1"
  , hline = FALSE
)
#---------------------------------- Rho Intercept -----------------------------------------#


#-------------------------------- Main Effects --------------------------------------------#
### BAD: attention effect
get_violin(
  ( plogis(ex_toj_color_post$population_logit_rho_intercept_mean + ex_toj_color_post$population_logit_rho_attention_effect_mean/2 )
    - plogis(ex_toj_color_post$population_logit_rho_intercept_mean - ex_toj_color_post$population_logit_rho_attention_effect_mean/2 ) )
  , "Probability of Encoding Attention Effect Mean"
  , y_lab = "\u03C1 (Attended - Unattended)"
)

# GOOD: attention effect
get_violin(
  ( rho_pos_attn_pos_probe - rho_neg_attn_pos_probe
  + rho_pos_attn_neg_probe - rho_neg_attn_neg_probe 
  )/2
  , "Probability of Encoding Attention Effect Mean"
  , y_lab = "\u03C1 (Attended - Unattended)"
)

# BETTER: attention effect
get_violin(
  ex_toj_color_post$population_logit_rho_attention_effect_mean
  , "Probability of Encoding Attention Effect Mean"
  , y_lab = "Log-odds of \u03C1 (Attended - Unattended)"
)


###  BAD: probe effect
get_violin(
    ( plogis(ex_toj_color_post$population_logit_rho_intercept_mean + ex_toj_color_post$population_logit_rho_probe_duration_effect_mean/2 )
    - plogis(ex_toj_color_post$population_logit_rho_intercept_mean - ex_toj_color_post$population_logit_rho_probe_duration_effect_mean/2 ) )
  , "Probability of Encoding\nProbe Duration Effect Mean"
  , y_lab = "\u03C1 (Long - Short)"
)

# GOOD: probe effect
get_violin(
  ( rho_pos_attn_pos_probe - rho_pos_attn_neg_probe
  + rho_neg_attn_pos_probe - rho_neg_attn_neg_probe
  )/2
  , "Probability of Encoding\nProbe Duration Effect Mean"
  , y_lab = "\u03C1 (Long - Short)"
)

#  BETTER: probe effect
get_violin(
  ex_toj_color_post$population_logit_rho_probe_duration_effect_mean
  , "Probability of Encoding\nProbe Duration Effect Mean"
  , y_lab = "Log-odds of \u03C1 (Long - Short)"
)
#-------------------------------- Main Effects --------------------------------------------#


#------------------------------- Two-way Interactions -------------------------------------#
### BAD:
get_violin(
  # don't devide effects by two...
  plogis(ex_toj_color_post$population_logit_rho_intercept_mean + ex_toj_color_post$population_logit_rho_attention_probe_duration_interaction_effect_mean )
  - plogis(ex_toj_color_post$population_logit_rho_intercept_mean - ex_toj_color_post$population_logit_rho_attention_probe_duration_interaction_effect_mean ) 
  , "Probability of Encoding\nProbe Duration Attention\nInteraction Effect Mean"
  , y_lab = "\u03C1"
)

# GOOD:
get_violin(
  (rho_pos_attn_pos_probe - rho_neg_attn_pos_probe) - (rho_pos_attn_neg_probe - rho_neg_attn_neg_probe)
  , "Probability of Encoding\nProbe Duration Attention\nInteraction Effect Mean"
  , y_lab = "\u03C1"
)

# BETTER:
get_violin(
  2 * ex_toj_color_post$population_logit_rho_attention_probe_duration_interaction_effect_mean
  , "Probability of Encoding\nProbe Duration Attention\nInteraction Effect Mean"
  , y_lab = "Log-odds of \u03C1"
)

# get by condition
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
  , c("Probability of Encoding\nAttention Effect \nGiven Long\nProbe Duration","Probability of Encoding\nAttention Effect\nGiven Short\nProbe Duration")
  , y_lab = "\u03C1 (Attended - Unattended)"
)
get_95_HDI(
  plogis(ex_toj_color_post$population_logit_rho_intercept_mean - ex_toj_color_post$population_logit_rho_probe_duration_effect_mean/2 
         + (ex_toj_color_post$population_logit_rho_attention_effect_mean - ex_toj_color_post$population_logit_rho_attention_probe_duration_interaction_effect_mean )/2 )
  - plogis(ex_toj_color_post$population_logit_rho_intercept_mean - ex_toj_color_post$population_logit_rho_probe_duration_effect_mean/2 
           - (ex_toj_color_post$population_logit_rho_attention_effect_mean - ex_toj_color_post$population_logit_rho_attention_probe_duration_interaction_effect_mean )/2 ) 
)

# CHECK simpler 
get_violin(
  c(
    rho_pos_attn_pos_probe - rho_neg_attn_pos_probe
    , rho_pos_attn_neg_probe - rho_neg_attn_neg_probe )
  , c("Probability of Encoding\nAttention Effect \nGiven Long\nProbe Duration","Probability of Encoding\nAttention Effect\nGiven Short\nProbe Duration")
  , y_lab = "\u03C1 (Attended - Unattended)"
)

# BETTER: 
get_violin(
  c(
    ( ex_toj_color_post$population_logit_rho_intercept_mean + ex_toj_color_post$population_logit_rho_probe_duration_effect_mean/2 
           + (ex_toj_color_post$population_logit_rho_attention_effect_mean + ex_toj_color_post$population_logit_rho_attention_probe_duration_interaction_effect_mean )/2 )
    - ( ex_toj_color_post$population_logit_rho_intercept_mean + ex_toj_color_post$population_logit_rho_probe_duration_effect_mean/2 
             - (ex_toj_color_post$population_logit_rho_attention_effect_mean + ex_toj_color_post$population_logit_rho_attention_probe_duration_interaction_effect_mean )/2 )
    , ( ex_toj_color_post$population_logit_rho_intercept_mean - ex_toj_color_post$population_logit_rho_probe_duration_effect_mean/2 
             + (ex_toj_color_post$population_logit_rho_attention_effect_mean - ex_toj_color_post$population_logit_rho_attention_probe_duration_interaction_effect_mean )/2 )
    - ( ex_toj_color_post$population_logit_rho_intercept_mean - ex_toj_color_post$population_logit_rho_probe_duration_effect_mean/2 
             - (ex_toj_color_post$population_logit_rho_attention_effect_mean - ex_toj_color_post$population_logit_rho_attention_probe_duration_interaction_effect_mean )/2 )
  )
  , c("Probability of Encoding\nAttention Effect \nGiven Long\nProbe Duration","Probability of Encoding\nAttention Effect\nGiven Short\nProbe Duration")
  , y_lab = "Log-odds of \u03C1 (Attended - Unattended)"
)
get_95_HDI(
  ( ex_toj_color_post$population_logit_rho_intercept_mean - ex_toj_color_post$population_logit_rho_probe_duration_effect_mean/2 
    + (ex_toj_color_post$population_logit_rho_attention_effect_mean - ex_toj_color_post$population_logit_rho_attention_probe_duration_interaction_effect_mean )/2 )
  - ( ex_toj_color_post$population_logit_rho_intercept_mean - ex_toj_color_post$population_logit_rho_probe_duration_effect_mean/2 
      - (ex_toj_color_post$population_logit_rho_attention_effect_mean - ex_toj_color_post$population_logit_rho_attention_probe_duration_interaction_effect_mean )/2 )
)
#------------------------------- Two-way Interactions -------------------------------------#



#------------------------------------------------------------------------------------------#
#--------------------------------- Kappa Scale --------------------------------------------#
#------------------------------------------------------------------------------------------#
 
#---------------------------------- Kappa Conditions --------------------------------------#
kappa_pos_attn_pos_probe = 
  exp(
    ex_toj_color_post$population_log_kappa_intercept_mean + ex_toj_color_post$population_log_kappa_probe_duration_effect_mean/2
    + (ex_toj_color_post$population_log_kappa_attention_effect_mean + ex_toj_color_post$population_log_kappa_attention_probe_duration_interaction_effect_mean)/2
  )

kappa_pos_attn_neg_probe = 
  exp(
    ex_toj_color_post$population_log_kappa_intercept_mean - ex_toj_color_post$population_log_kappa_probe_duration_effect_mean/2
    + (ex_toj_color_post$population_log_kappa_attention_effect_mean - ex_toj_color_post$population_log_kappa_attention_probe_duration_interaction_effect_mean)/2
  )

kappa_neg_attn_neg_probe = 
  exp(
    ex_toj_color_post$population_log_kappa_intercept_mean - ex_toj_color_post$population_log_kappa_probe_duration_effect_mean/2
    - (ex_toj_color_post$population_log_kappa_attention_effect_mean - ex_toj_color_post$population_log_kappa_attention_probe_duration_interaction_effect_mean)/2
  )

kappa_neg_attn_pos_probe = 
  exp(
    ex_toj_color_post$population_log_kappa_intercept_mean + ex_toj_color_post$population_log_kappa_probe_duration_effect_mean/2
    - (ex_toj_color_post$population_log_kappa_attention_effect_mean + ex_toj_color_post$population_log_kappa_attention_probe_duration_interaction_effect_mean)/2
  )
#---------------------------------- Kappa Conditions --------------------------------------#


#---------------------------------- Kappa Intercept ---------------------------------------#
### BAD:
get_violin(
  exp(ex_toj_color_post$population_log_kappa_intercept_mean)
  , "Fidelity of Memory Intercept Mean"
  , y_lab = "\u03BA"
  , hline = FALSE
)

# GOOD:
get_violin(
  ( kappa_neg_attn_pos_probe
  + kappa_neg_attn_neg_probe
  + kappa_pos_attn_neg_probe
  + kappa_pos_attn_pos_probe
  )/4
  , "Fidelity of Memory Intercept Mean"
  , y_lab = "\u03BA"
  , hline = FALSE
)

# BETTER::
get_violin(
  ex_toj_color_post$population_log_kappa_intercept_mean
  , "Fidelity of Memory Intercept Mean"
  , y_lab = "Log of \u03BA"
  , hline = FALSE
)
#---------------------------------- Kappa Intercept ---------------------------------------#


#-------------------------------- Main Effects --------------------------------------------#
### BAD: attention effect
get_violin(
  ( exp(ex_toj_color_post$population_log_kappa_intercept_mean + ex_toj_color_post$population_log_kappa_attention_effect_mean/2 )
    - exp(ex_toj_color_post$population_log_kappa_intercept_mean - ex_toj_color_post$population_log_kappa_attention_effect_mean/2 ) )
  , "Fidelity of Encoding Attention Effect Mean"
  , y_lab = "\u03BA (Attended - Unattended)"
)

# GOOD: attention effect
get_violin(
  ( kappa_pos_attn_pos_probe - kappa_neg_attn_pos_probe
  + kappa_pos_attn_neg_probe - kappa_neg_attn_neg_probe 
  )/2
  , "Fidelity of Encoding Attention Effect Mean"
  , y_lab = "\u03BA (Attended - Unattended)"
)

# BETTER: attention effect
get_violin(
  ex_toj_color_post$population_logit_rho_attention_effect_mean
  , "Fidelity of Encoding Attention Effect Mean"
  , y_lab = "Log of \u03BA (Attended - Unattended)"
)


### BAD: probe effect
get_violin(
  ( exp(ex_toj_color_post$population_log_kappa_intercept_mean + ex_toj_color_post$population_log_kappa_probe_duration_effect_mean/2 )
    - exp(ex_toj_color_post$population_log_kappa_intercept_mean - ex_toj_color_post$population_log_kappa_probe_duration_effect_mean/2 ) )
  , "Fidelity of Encoding Probe Duration Effect Mean"
  , y_lab = "\u03BA (Long - Short)"
)

# GOOD: probe effect
get_violin(
  ( kappa_pos_attn_pos_probe - kappa_pos_attn_neg_probe
  + kappa_neg_attn_pos_probe - kappa_neg_attn_neg_probe
  )/2
  , "Fidelity of Encoding Probe Duration Effect Mean"
  , y_lab = "\u03BA (Long - Short)"
)

# BETTER: probe effect
get_violin(
  ex_toj_color_post$population_log_kappa_probe_duration_effect_mean
  , "Fidelity of Encoding Probe Duration Effect Mean"
  , y_lab = "Log of \u03BA (Long - Short)"
)
#-------------------------------- Main Effects --------------------------------------------#


#------------------------------- Two-way Interactions -------------------------------------#
### BAD:
get_violin(
  # don't devide effect by two...
  ( exp(ex_toj_color_post$population_log_kappa_intercept_mean + ex_toj_color_post$population_log_kappa_attention_probe_duration_interaction_effect_mean)
    - exp(ex_toj_color_post$population_log_kappa_intercept_mean - ex_toj_color_post$population_log_kappa_attention_probe_duration_interaction_effect_mean) )
  , "Fidelity of Encoding\nProbe Duration Attention\nInteraction Effect Mean"
  , y_lab = "\u03BA"
)

# GOOD:
get_violin(
  (kappa_pos_attn_pos_probe - kappa_neg_attn_pos_probe) - (kappa_pos_attn_neg_probe - kappa_neg_attn_neg_probe)
  , "Fidelity of Encoding\nProbe Duration Attention\nInteraction Effect Mean"
  , y_lab = "\u03BA"
)

# BETTER:
get_violin(
  2 * ex_toj_color_post$population_log_kappa_attention_probe_duration_interaction_effect_mean
  , "Fidelity of Encoding\nProbe Duration Attention\nInteraction Effect Mean"
  , y_lab = "Log of \u03BA"
)

# by condition
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
  , c("Fidelity of Encoding\nAttention Effect \nGiven Long\nProbe Duration","Fidelity of Encoding\nAttention Effect\nGiven Short\nProbe Duration")
  , y_lab = "\u03BA (Attended - Unattended)"
)
get_95_HDI(
  exp(ex_toj_color_post$population_log_kappa_intercept_mean - ex_toj_color_post$population_log_kappa_probe_duration_effect_mean/2 
      + (ex_toj_color_post$population_log_kappa_attention_effect_mean - ex_toj_color_post$population_log_kappa_attention_probe_duration_interaction_effect_mean )/2 )
  - exp(ex_toj_color_post$population_log_kappa_intercept_mean - ex_toj_color_post$population_log_kappa_probe_duration_effect_mean/2 
        - (ex_toj_color_post$population_log_kappa_attention_effect_mean - ex_toj_color_post$population_log_kappa_attention_probe_duration_interaction_effect_mean )/2 ) 
)

# simpler
get_violin(
  c(kappa_pos_attn_pos_probe - kappa_neg_attn_pos_probe
    , kappa_pos_attn_neg_probe - kappa_neg_attn_neg_probe)
  , c("Fidelity of Encoding\nAttention Effect \nGiven Long\nProbe Duration","Fidelity of Encoding\nAttention Effect\nGiven Short\nProbe Duration")
  , y_lab = "\u03BA (Attended - Unattended)"
)

# BETTER:
get_violin(
  c(
    (ex_toj_color_post$population_log_kappa_intercept_mean + ex_toj_color_post$population_log_kappa_probe_duration_effect_mean/2 
        + (ex_toj_color_post$population_log_kappa_attention_effect_mean + ex_toj_color_post$population_log_kappa_attention_probe_duration_interaction_effect_mean )/2 )
    - (ex_toj_color_post$population_log_kappa_intercept_mean + ex_toj_color_post$population_log_kappa_probe_duration_effect_mean/2 
          - (ex_toj_color_post$population_log_kappa_attention_effect_mean + ex_toj_color_post$population_log_kappa_attention_probe_duration_interaction_effect_mean )/2 )
    , (ex_toj_color_post$population_log_kappa_intercept_mean - ex_toj_color_post$population_log_kappa_probe_duration_effect_mean/2 
          + (ex_toj_color_post$population_log_kappa_attention_effect_mean - ex_toj_color_post$population_log_kappa_attention_probe_duration_interaction_effect_mean )/2 )
    - (ex_toj_color_post$population_log_kappa_intercept_mean - ex_toj_color_post$population_log_kappa_probe_duration_effect_mean/2 
          - (ex_toj_color_post$population_log_kappa_attention_effect_mean - ex_toj_color_post$population_log_kappa_attention_probe_duration_interaction_effect_mean )/2 ) 
  )
  , c("Fidelity of Encoding\nAttention Effect \nGiven Long\nProbe Duration","Fidelity of Encoding\nAttention Effect\nGiven Short\nProbe Duration")
  , y_lab = "Log of \u03BA (Attended - Unattended)"
)
get_95_HDI(
  (ex_toj_color_post$population_log_kappa_intercept_mean - ex_toj_color_post$population_log_kappa_probe_duration_effect_mean/2 
      + (ex_toj_color_post$population_log_kappa_attention_effect_mean - ex_toj_color_post$population_log_kappa_attention_probe_duration_interaction_effect_mean )/2 )
  - (ex_toj_color_post$population_log_kappa_intercept_mean - ex_toj_color_post$population_log_kappa_probe_duration_effect_mean/2 
        - (ex_toj_color_post$population_log_kappa_attention_effect_mean - ex_toj_color_post$population_log_kappa_attention_probe_duration_interaction_effect_mean )/2 ) 
)
#------------------------------- Two-way Interactions -------------------------------------#





#------------------------------------------------------------------------------------------#
#--------------------------------- Graphs -------------------------------------------------#
#------------------------------------------------------------------------------------------#

#---------------------------------- NCFs --------------------------------------------------#
###BADS:
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


### GOODs:
yLeftFirst = pnorm(
  -250:250
  , mean = ( 
      median(
      (
        ( ex_toj_color_post$population_pss_intercept_mean - ex_toj_color_post$population_pss_judgement_type_effect_mean/2 + ex_toj_color_post$population_pss_probe_duration_effect_mean/2
        - ( ex_toj_color_post$population_pss_attention_effect_mean - ex_toj_color_post$population_pss_attention_judgement_type_interaction_effect_mean + ex_toj_color_post$population_pss_attention_probe_duration_interaction_effect_mean)/2)
      + ( ex_toj_color_post$population_pss_intercept_mean - ex_toj_color_post$population_pss_judgement_type_effect_mean/2 - ex_toj_color_post$population_pss_probe_duration_effect_mean/2
        - ( ex_toj_color_post$population_pss_attention_effect_mean - ex_toj_color_post$population_pss_attention_judgement_type_interaction_effect_mean - ex_toj_color_post$population_pss_attention_probe_duration_interaction_effect_mean)/2)
      )/2 * 250 
    )
  )
  , sd = median(
     (neg_attn_neg_judge_neg_probe + neg_attn_neg_judge_pos_probe)/2 * 250
    )
)

yLeftSecond = pnorm(
  -250:250
  , mean = ( 
      median(
      (
        ( ex_toj_color_post$population_pss_intercept_mean + ex_toj_color_post$population_pss_judgement_type_effect_mean/2 + ex_toj_color_post$population_pss_probe_duration_effect_mean/2
          - ( ex_toj_color_post$population_pss_attention_effect_mean + ex_toj_color_post$population_pss_attention_judgement_type_interaction_effect_mean + ex_toj_color_post$population_pss_attention_probe_duration_interaction_effect_mean)/2)
        + ( ex_toj_color_post$population_pss_intercept_mean + ex_toj_color_post$population_pss_judgement_type_effect_mean/2 - ex_toj_color_post$population_pss_probe_duration_effect_mean/2
            - ( ex_toj_color_post$population_pss_attention_effect_mean + ex_toj_color_post$population_pss_attention_judgement_type_interaction_effect_mean - ex_toj_color_post$population_pss_attention_probe_duration_interaction_effect_mean)/2)
      )/2 * 250 
    )
  )
  , sd = median(
    (neg_attn_pos_judge_neg_probe + neg_attn_pos_judge_pos_probe)/2 * 250
    ) 
)

yRightFirst = pnorm(
  -250:250
  , mean = ( 
      median(
      (
        ( ex_toj_color_post$population_pss_intercept_mean - ex_toj_color_post$population_pss_judgement_type_effect_mean/2 + ex_toj_color_post$population_pss_probe_duration_effect_mean/2
          + ( ex_toj_color_post$population_pss_attention_effect_mean - ex_toj_color_post$population_pss_attention_judgement_type_interaction_effect_mean + ex_toj_color_post$population_pss_attention_probe_duration_interaction_effect_mean)/2)
        + ( ex_toj_color_post$population_pss_intercept_mean - ex_toj_color_post$population_pss_judgement_type_effect_mean/2 - ex_toj_color_post$population_pss_probe_duration_effect_mean/2
            + ( ex_toj_color_post$population_pss_attention_effect_mean - ex_toj_color_post$population_pss_attention_judgement_type_interaction_effect_mean - ex_toj_color_post$population_pss_attention_probe_duration_interaction_effect_mean)/2)
      )/2 * 250 
    )
  )
  , sd = median(
    (pos_attn_neg_judge_neg_probe + pos_attn_neg_judge_pos_probe)/2 * 250
  )
)

yRightSecond = pnorm(
  -250:250
  , mean = ( 
    median(
      (
        ( ex_toj_color_post$population_pss_intercept_mean + ex_toj_color_post$population_pss_judgement_type_effect_mean/2 + ex_toj_color_post$population_pss_probe_duration_effect_mean/2
          + ( ex_toj_color_post$population_pss_attention_effect_mean + ex_toj_color_post$population_pss_attention_judgement_type_interaction_effect_mean + ex_toj_color_post$population_pss_attention_probe_duration_interaction_effect_mean)/2)
        + ( ex_toj_color_post$population_pss_intercept_mean + ex_toj_color_post$population_pss_judgement_type_effect_mean/2 - ex_toj_color_post$population_pss_probe_duration_effect_mean/2
            + ( ex_toj_color_post$population_pss_attention_effect_mean + ex_toj_color_post$population_pss_attention_judgement_type_interaction_effect_mean - ex_toj_color_post$population_pss_attention_probe_duration_interaction_effect_mean)/2)
      )/2 * 250 
    )
  )
  , sd =  median( 
    (pos_attn_pos_judge_neg_probe + pos_attn_pos_judge_pos_probe)/2 * 250
  )
)

### Execute
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






