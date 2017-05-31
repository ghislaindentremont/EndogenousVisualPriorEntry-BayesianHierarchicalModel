# library(shinystan)
library(plyr)
library(coda)
library(ggplot2)
library(ggmcmc)
library(CircStats)
library(grid)
library(reshape2)
# library(sprintfr)
library(gtools)

setwd("~/Documents/Experiments/TOJ/Follow-Up")
load("FollowUptoj_color_post_Aug8th2016")
load("FollowUp_color_trials.Rdata")
load("FollowUp_toj_trials.Rdata")
# 10,000 * 6 / 2 = 30,000 iters 
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


#-------------------------------------- Get SDs -------------------------------------------#
get_violin(
  (
  (ex_toj_color_post$population_pss_intercept_mean + ex_toj_color_post$population_pss_intercept_sd/2)
  -  (ex_toj_color_post$population_pss_intercept_mean - ex_toj_color_post$population_pss_intercept_sd/2)
  ) * 250
  , labels = "PSS SD"
  , y_lab = "PSS (ms)"
  )

get_violin(
  (
    exp(ex_toj_color_post$population_log_jnd_intercept_mean + ex_toj_color_post$population_log_jnd_intercept_sd/2)
    -  exp(ex_toj_color_post$population_log_jnd_intercept_mean - ex_toj_color_post$population_log_jnd_intercept_sd/2)
  ) * 250
  , labels = "JND SD"
  , y_lab = "JND (ms)"
)

get_violin(
  (
    plogis(ex_toj_color_post$population_logit_rho_intercept_mean + ex_toj_color_post$population_logit_rho_intercept_sd/2)
    -  plogis(ex_toj_color_post$population_logit_rho_intercept_mean - ex_toj_color_post$population_logit_rho_intercept_sd/2)
  ) 
  , labels = "Probability SD"
  , y_lab = "Probability"
)

get_violin(
  (
    exp(ex_toj_color_post$population_log_kappa_intercept_mean + ex_toj_color_post$population_log_kappa_intercept_sd/2)
    -  exp(ex_toj_color_post$population_log_kappa_intercept_mean - ex_toj_color_post$population_log_kappa_intercept_sd/2)
  ) 
  , labels = "Fidelity SD"
  , y_lab = "Fidelity"
)
#-------------------------------------- Get SDs -------------------------------------------#


#-------------------------------------- Get Corr ------------------------------------------#
pos_corr_use2 = data.frame(
  "PSS Intercept v PSS Effect" = pos_corr[pos_corr$parameter == "value.1.2",]$value
  , "PSS Intercept v JND Intercept" = pos_corr[pos_corr$parameter == "value.1.3",]$value
  , "PSS Intercept v JND Effect" =  pos_corr[pos_corr$parameter == "value.1.4",]$value
  , "PSS Intercept v Probability Intercept" =  pos_corr[pos_corr$parameter == "value.1.5",]$value
  , "PSS Intercept v Fidelity Intercept" = pos_corr[pos_corr$parameter == "value.1.6",]$value
  , "PSS Intercept v Probability Effect" = pos_corr[pos_corr$parameter == "value.1.7",]$value
  , "PSS Intercept v Fidelity Effect" =  pos_corr[pos_corr$parameter == "value.1.8",]$value
  
  , "PSS Effect v JND Intercept" = pos_corr[pos_corr$parameter == "value.2.3",]$value
  , "PSS Effect v JND Effect" = pos_corr[pos_corr$parameter == "value.2.4",]$value
  , "PSS Effect v Probability Intercept" = pos_corr[pos_corr$parameter == "value.2.5",]$value
  , "PSS Effect v Fidelity Intercept" = pos_corr[pos_corr$parameter == "value.2.6",]$value
  , "PSS Effect v Probability Effect" = pos_corr[pos_corr$parameter == "value.2.7",]$value
  , "PSS Effect v Fidelity Effect" = pos_corr[pos_corr$parameter == "value.2.8",]$value
  
  , "JND Intercept v JND Effect" = pos_corr[pos_corr$parameter == "value.3.4",]$value
  , "JND Intercept v Probability Intercept" = pos_corr[pos_corr$parameter == "value.3.5",]$value
  , "JND Intercept v Fidelity Intercept" = pos_corr[pos_corr$parameter == "value.3.6",]$value
  , "JND Intercept v Probability Effect" = pos_corr[pos_corr$parameter == "value.3.7",]$value
  , "JND Intercept v Fidelity Effect" = pos_corr[pos_corr$parameter == "value.3.8",]$value
  
  , "JND Effect v Probability Intercept" = pos_corr[pos_corr$parameter == "value.4.5",]$value
  , "JND Effect v Fidelity Intercept" = pos_corr[pos_corr$parameter == "value.4.6",]$value
  , "JND Effect v Probability Effect" = pos_corr[pos_corr$parameter == "value.4.7",]$value
  , "JND Effect v Fidelity Effect" = pos_corr[pos_corr$parameter == "value.4.8",]$value
  
  , "Probability Intercept v Fidelity Intercept" = pos_corr[pos_corr$parameter == "value.5.6",]$value
  , "Probability Intercept v Probability Effect" = pos_corr[pos_corr$parameter == "value.5.7",]$value
  , "Probability Intercept v Fidelity Effect" = pos_corr[pos_corr$parameter == "value.5.8",]$value
  
  , "Fidelity Intercept v Probability Effect" = pos_corr[pos_corr$parameter == "value.6.7",]$value
  , "Fidelity Intercept v Fidelity Effect" = pos_corr[pos_corr$parameter == "value.6.8",]$value
  
  , "Probability Effect v Fidelity Effect" = pos_corr[pos_corr$parameter == "value.7.8",]$value
  
  , check.names = FALSE
)

pos_corr_use = melt(pos_corr_use2)
names(pos_corr_use)[1] = "parameter" 

ggplot(
  data = pos_corr_use
  , aes(x = reorder(parameter, value, FUN = median), y = value)
)+
  labs(x = "", y = "Correlation Coefficient (r)")+
  stat_summary(fun.data = get_95_HDI, size = 0.7)+
  stat_summary(fun.data = get_50_HDI, size = 2.5)+  
  geom_hline(yintercept = 0, linetype = 2, size = 1)+
  coord_flip()+
  theme_gray(base_size = 30)+
  theme(panel.grid.major = element_line(size = 1.5)
        ,panel.grid.minor = element_line(size = 1)
        , axis.ticks.x = element_blank()) 
#-------------------------------------- Get Corr ------------------------------------------#


#--------------------------------- Notable Corrs ------------------------------------------#
get_corr(cor_val = "value.2.6", cor_lab = "PSS Effect v Fidelity Intercept")
get_corr(cor_val = "value.5.7", cor_lab = "Probability Intercept v Probability Effect")
get_corr(cor_val = "value.3.5", cor_lab = "JND Intercept v Probability Intercept")
get_corr(cor_val = "value.1.3", cor_lab = "PSS Intercept v JND Intercept")
#--------------------------------- Notable Corrs ------------------------------------------#


#-------------------------------------- Raw Data ------------------------------------------#
source("/Users/ghislaindentremont/Documents/TOJ/EndogenousVisualPriorEntry-BayesianHierarchicalModel/FollowUp/conventional_analysis/followup_analysis.R")
#-------------------------------------- Raw Data ------------------------------------------#




#----------------------- Prior Entry vs. Probability Differences --------------------------#
pssattentioneffectmean = extract_samples("population_pss_attention_effect_mean")
pssattentioneffectsd = extract_samples("population_pss_effect_sd")
pssprobeinteractioneffectmean = extract_samples("population_pss_attention_probe_duration_interaction_effect_mean")
pssjudgementinteractioneffectmean = extract_samples("population_pss_attention_judgement_type_interaction_effect_mean")

# conditions
probedurationcondition = ifelse(aggregate(longprobe ~ id, data = color_trials, FUN = unique)$longprobe == "FALSE", -1, 1)
judgementtypecondition = ifelse(aggregate(toj_judgement_type ~ id, data = color_trials, FUN = unique)$toj_judgement_type == "first", -1, 1) 

pssdiff_ids = ddply(
  .data = betas
  , .variables = .(participant)
  , .fun = function(x){
    i = unique(x$participant)
    x_use = x[x$parameter ==  "population_pss_attention_effect_mean",]$value
    pssdiff = median(
      pssattentioneffectmean + pssattentioneffectsd*x_use 
      + probedurationcondition[i]*pssprobeinteractioneffectmean
      + judgementtypecondition[i]*pssjudgementinteractioneffectmean
    ) 
    df = data.frame(pssdiff*250, probedurationcondition[i], judgementtypecondition[i])
    names(df) = c("pssdiff", "probecondition", "judgementcondition")
    return(df)
  }
)

rhoattentioneffectmean = extract_samples("population_logit_rho_attention_effect_mean")
rhoattentioneffectsd = extract_samples("population_logit_rho_effect_sd")
rhointeractioneffectmean = extract_samples("population_logit_rho_attention_probe_duration_interaction_effect_mean")

rhodiff_ids = ddply(
  .data = betas
  , .variables = .(participant)
  , .fun = function(x){
    i = unique(x$participant)
    x_use = x[x$parameter ==  "population_logit_rho_attention_effect_mean",]$value
    rhodiff = median(
      rhoattentioneffectmean + rhoattentioneffectsd*x_use 
      + probedurationcondition[i]*rhointeractioneffectmean
    ) 
    df = data.frame(rhodiff, probedurationcondition[i])
    names(df) = c("rhodiff", "probecondition")
    return(df)
  }
)

# rid redundant Cs
rhodiff_vs_pssdiff2 = cbind(pssdiff_ids, rhodiff_ids[-c(1,3)])

# groupings = ddply(
#   .data = rhodiff_vs_pssdiff
#   , .variables = .(participant)
#   , .fun = function(x){
#     i = unique(x$participant)
#     if (x$probecondition == 1 & x$judgementcondition == 1) {
#       x$grouping = "long and second"
#     } else if (x$probecondition == 1 & x$judgementcondition == -1) {
#       x$grouping = "long and first"
#     } else if (x$probecondition == -1 & x$judgementcondition == -1) {
#       x$grouping = "short and first"
#     } else if (x$probecondition == -1 & x$judgementcondition == 1) {
#       x$grouping = "short and second"
#     } else {
#       print("ERROR")
#     }
#   }
# )
# 
# # RECONSIDER 'merge' argument
# rhodiff_vs_pssdiff = merge(rhodiff_vs_pssdiff, groupings)

rhodiff_vs_pssdiff = cbind(rhodiff_vs_pssdiff2, ids[,c("PSS.Difference", "Probability.Difference")] )

# shrinkage for PSS differences
plot(rhodiff_vs_pssdiff$participant, rhodiff_vs_pssdiff$PSS.Difference, col = "red")
abline(h = mean(rhodiff_vs_pssdiff$PSS.Difference), col = 'red')
points(rhodiff_vs_pssdiff$pssdiff, col = "blue")
abline(h = mean(rhodiff_vs_pssdiff$pssdiff), col = "blue")
# correlation between Bayesian estimates and actual values
cor(rhodiff_vs_pssdiff$pssdiff, rhodiff_vs_pssdiff$PSS.Difference)

# shrinkage for probability Differences
plot(rhodiff_vs_pssdiff$participant, rhodiff_vs_pssdiff$Probability.Difference, col = "red")
abline(h = mean(rhodiff_vs_pssdiff$Probability.Difference), col = 'red')
points(rhodiff_vs_pssdiff$rhodiff, col = "blue")
abline(h = mean(rhodiff_vs_pssdiff$rhodiff), col = "blue")
# correlation between Bayesian estimates and actual values
cor(rhodiff_vs_pssdiff$rhodiff, rhodiff_vs_pssdiff$Probability.Difference)

# see 2-D
ggplot(data = rhodiff_vs_pssdiff, aes(x = rhodiff, y = pssdiff))+
  geom_point(size = 2, aes(color = "Bayesian"))+
  scale_x_continuous(name = "Log-Odds of Probability Differences")+
  scale_y_continuous(name = "PSS Differences") +
  geom_vline(xintercept = 0, linetype = 2, size = 1)+
  geom_hline(yintercept = 0, linetype = 2, size = 1)+
  geom_point(data = ids, aes(x = Probability.Difference, y = PSS.Difference,  color = "Raw"), size = 2)+
  scale_color_discrete(name = "")+
  theme_gray(base_size = 30)+
  theme(
    panel.grid.major = element_line(size = 1.5)
    , panel.grid.minor = element_line(size = 1)
    , strip.background = element_blank()
    , strip.text.x = element_blank() 
  ) 

# plot with conditions
ggplot(data = rhodiff_vs_pssdiff, aes(x = rhodiff, y = pssdiff, color = factor(probedurationcondition), shape = factor(judgementtypecondition)))+ # color = factor(V1)))+
  geom_point(size = 4)+
  scale_x_continuous(name = "Log-Odds of Probability Differences")+
  scale_y_continuous(name = "PSS Differences") +
  scale_color_discrete(name = "Probe Duration", labels = c("Short", "Long")) +
  scale_shape_discrete(name = "Judgement Type", labels = c("First", "Second"))+
  geom_vline(xintercept = 0, linetype = 2, size = 1)+
  geom_hline(yintercept = 0, linetype = 2, size = 1)+
# scale_color_discrete(name = "Between-Subject Conditions", labels = c("'Long' and 'First'", "'Long' and 'Second'", "'Short' and 'First'", "'Short' and 'Second'"))+
  theme_gray(base_size = 30)+
  theme(
    panel.grid.major = element_line(size = 1.5)
    , panel.grid.minor = element_line(size = 1)
    , strip.background = element_blank()
    , strip.text.x = element_blank() 
  ) 
#----------------------- Prior Entry vs. Probability Differences --------------------------#


#----------------------- Prior Entry vs. Probability EFFECTS ------------------------------#
priorentry_ids = ddply(
  .data = betas
  , .variables = .(participant)
  , .fun = function(x){
    i = unique(x$participant)
    x_use = x[x$parameter ==  "population_pss_attention_effect_mean",]$value
    priorentry = median(
      pssattentioneffectmean + pssattentioneffectsd*x_use 
    ) 
    df = data.frame(priorentry*250)
    names(df) = c("priorentry")
    return(df)
  }
)

rhoeffect_ids = ddply(
  .data = betas
  , .variables = .(participant)
  , .fun = function(x){
    i = unique(x$participant)
    x_use = x[x$parameter ==  "population_logit_rho_attention_effect_mean",]$value
    rhoeffect = median(
      rhoattentioneffectmean + rhoattentioneffectsd*x_use 
    ) 
    df = data.frame(rhoeffect)
    names(df) = c("rhoeffect")
    return(df)
  }
)

# rid redundant Cs
rhoeffect_vs_priorentry2 = cbind(rhoeffect_ids, priorentry_ids[-1])

# shrinkage for prior entry effects
rhoeffect_vs_priorentry = cbind(rhoeffect_vs_priorentry2, ids[,c("PSS.Effect", "Probability.Effect")] )

plot(rhoeffect_vs_priorentry$participant, rhoeffect_vs_priorentry$PSS.Effect, col = "red")
abline(h = mean(rhoeffect_vs_priorentry$PSS.Effect), col = 'red')
points(rhoeffect_vs_priorentry$priorentry, col = "blue")
abline(h = mean(rhoeffect_vs_priorentry$priorentry), col = "blue")
# correlation between Bayesian estimate and raw value
cor(rhoeffect_vs_priorentry$priorentry, rhoeffect_vs_priorentry$PSS.Effect)

# shrinkage for probability effects
plot(rhoeffect_vs_priorentry$participant, rhoeffect_vs_priorentry$Probability.Effect, col = "red")
abline(h = mean(rhoeffect_vs_priorentry$Probability.Effect), col = 'red')
points(rhoeffect_vs_priorentry$rhoeffect, col = "blue")
abline(h = mean(rhoeffect_vs_priorentry$rhoeffect), col = "blue")
# correlation between Bayesian estimate and raw value
cor(rhoeffect_vs_priorentry$rhoeffect, rhoeffect_vs_priorentry$Probability.Effect)

# see 2-D
ggplot(data = rhoeffect_vs_priorentry, aes(x = rhoeffect, y = priorentry))+
  geom_point(size = 2, aes(color = "Bayesian"))+
  scale_x_continuous(name = "Log-Odds of Probability Effects")+
  scale_y_continuous(name = "Prior Entry Effects") +
  geom_vline(xintercept = 0, linetype = 2, size = 1)+
  geom_hline(yintercept = 0, linetype = 2, size = 1)+
  geom_point(data = ids, aes(x = Probability.Effect, y = PSS.Effect,  color = "Raw"), size = 2)+
  scale_color_discrete(name = "")+
  theme_gray(base_size = 30)+
  theme(
    panel.grid.major = element_line(size = 1.5)
    , panel.grid.minor = element_line(size = 1)
    , strip.background = element_blank()
    , strip.text.x = element_blank() 
  ) 

# # plot with conditions
# ggplot(data = rhoeffect_vs_priorentry, aes(x = rhoeffect, y = priorentry, color = factor(probedurationcondition), shape = factor(judgementtypecondition)))+ # color = factor(V1)))+
#   geom_point(size = 4)+
#   scale_x_continuous(name = "Log-Odds of Probability Effects")+
#   scale_y_continuous(name = "Prior Entry Effects") +
#   scale_color_discrete(name = "Probe Duration", labels = c("Short", "Long")) +
#   scale_shape_discrete(name = "Judgement Type", labels = c("First", "Second"))+
#   geom_vline(xintercept = 0, linetype = 2, size = 1)+
#   geom_hline(yintercept = 0, linetype = 2, size = 1)+
#   # scale_color_discrete(name = "Between-Subject Conditions", labels = c("'Long' and 'First'", "'Long' and 'Second'", "'Short' and 'First'", "'Short' and 'Second'"))+
#   theme_gray(base_size = 30)+
#   theme(
#     panel.grid.major = element_line(size = 1.5)
#     , panel.grid.minor = element_line(size = 1)
#     , strip.background = element_blank()
#     , strip.text.x = element_blank() 
#   ) 
# #----------------------- Prior Entry vs. Probability EFFECTS ------------------------------#


#------------------ Probability Intercept vs. Probability Effects--------------------------#
rhointerceptmean =  extract_samples("population_logit_rho_intercept_mean")
rhointerceptsd = extract_samples("population_logit_rho_intercept_sd")

rhointercept_ids = ddply(
  .data = betas
  , .variables = .(participant)
  , .fun = function(x){
    i = unique(x$participant)
    x_use = x[x$parameter ==  "population_logit_rho_intercept_mean",]$value
    rhointercept = median(
      rhointerceptmean + rhointerceptsd*x_use 
    ) 
    df = data.frame(rhointercept)
    names(df) = c("rhointercept")
    return(df)
  }
)

# rid redundant Cs
rhoeffect_vs_rhointercept2 = cbind(rhointercept_ids, rhoeffect_ids[-1])

rhoeffect_vs_rhointercept = cbind(rhoeffect_vs_rhointercept2, ids[,c("Probability.Intercept", "Probability.Effect")] )

# shrinkage for probability intercepts
plot(rhoeffect_vs_rhointercept$participant, rhoeffect_vs_rhointercept$Probability.Intercept, col = "red")
abline(h = mean(rhoeffect_vs_rhointercept$Probability.Intercept), col = 'red')
points(rhoeffect_vs_rhointercept$rhointercept, col = "blue")
abline(h = mean(rhoeffect_vs_rhointercept$rhointercept), col = "blue")
# correlation between Bayesian estimate and raw value
cor(rhoeffect_vs_rhointercept$rhointercept, rhoeffect_vs_rhointercept$Probability.Intercept)

# see 2-D
ggplot(data = rhoeffect_vs_rhointercept, aes(x = rhointercept, y = rhoeffect))+
  geom_point(size = 2, aes(color = "Bayesian"))+
  scale_y_continuous(name = "Log-Odds of Probability Effects")+
  scale_x_continuous(name = "Log-Odds of Probability Intercepts") +
  geom_vline(xintercept = 0, linetype = 2, size = 1)+
  geom_hline(yintercept = 0, linetype = 2, size = 1)+
  geom_point(data = ids, aes(x = Probability.Intercept, y = Probability.Effect,  color = "Raw"), size = 2)+
  scale_color_discrete(name = "")+
  theme_gray(base_size = 30)+
  theme(
    panel.grid.major = element_line(size = 1.5)
    , panel.grid.minor = element_line(size = 1)
    , strip.background = element_blank()
    , strip.text.x = element_blank() 
  ) 

# Bayesian correlation
cor(rhoeffect_vs_rhointercept$rhointercept, rhoeffect_vs_rhointercept$rhoeffect)

# Raw correlation
cor(ids$Probability.Intercept, ids$Probability.Effect)

# # plot with conditions
# ggplot(data = rhoeffect_vs_rhointercept, aes(x = rhointercept, y = rhoeffect, color = factor(probedurationcondition), shape = factor(judgementtypecondition)))+ # color = factor(V1)))+
#   geom_point(size = 4)+
#   scale_y_continuous(name = "Log-Odds of Probability Effects")+
#   scale_x_continuous(name = "Log-Odds of Probability Intercepts") +
#   scale_color_discrete(name = "Probe Duration", labels = c("Short", "Long")) +
#   scale_shape_discrete(name = "Judgement Type", labels = c("First", "Second"))+
#   geom_vline(xintercept = 0, linetype = 2, size = 1)+
#   geom_hline(yintercept = 0, linetype = 2, size = 1)+
#   theme_gray(base_size = 30)+
#   theme(
#     panel.grid.major = element_line(size = 1.5)
#     , panel.grid.minor = element_line(size = 1)
#     , strip.background = element_blank()
#     , strip.text.x = element_blank() 
#   )
#------------------ Probability Intercept vs. Probability Effects--------------------------#


#------------------------- Prior Entry vs. Fidelity Differences ---------------------------#
kappaattentioneffectmean = extract_samples("population_log_kappa_attention_effect_mean")
kappaattentioneffectsd = extract_samples("population_log_kappa_effect_sd")
kappainteractioneffectmean = extract_samples("population_log_kappa_attention_probe_duration_interaction_effect_mean")

kappadiff_ids = ddply(
  .data = betas
  , .variables = .(participant)
  , .fun = function(x){
    i = unique(x$participant)
    x_use = x[x$parameter ==  "population_log_kappa_attention_effect_mean",]$value
    kappadiff = median(
      kappaattentioneffectmean + kappaattentioneffectsd*x_use 
      + probedurationcondition[i]*kappainteractioneffectmean
    ) 
    df = data.frame(kappadiff, probedurationcondition[i])
    names(df) = c("kappadiff", "probecondition")
    return(df)
  }
)

kappadiff_vs_pssdiff = cbind(kappadiff_ids[-c(1,3)], pssdiff_ids)

# shrinkage for prior entry effects
kappadiff_vs_pssdiff = cbind(kappadiff_vs_pssdiff, ids[,c("PSS.Difference", "Fidelity.Difference")] )

# shrinkage for Fidelity effects
plot(kappadiff_vs_pssdiff$participant, kappadiff_vs_pssdiff$Fidelity.Difference, col = "red")
abline(h = mean(kappadiff_vs_pssdiff$Fidelity.Difference), col = 'red')
points(kappadiff_vs_pssdiff$kappadiff, col = "blue")
abline(h = mean(kappadiff_vs_pssdiff$kappadiff), col = "blue")

# see 2-D
ggplot(data = kappadiff_vs_pssdiff, aes(x = kappadiff, y = pssdiff))+
  geom_point(size = 2, aes(color = "Bayesian"))+
  scale_x_continuous(name = "Log of Fidelity Differences")+
  scale_y_continuous(name = "Prior Entry Differences") +
  geom_vline(xintercept = 0, linetype = 2, size = 1)+
  geom_hline(yintercept = 0, linetype = 2, size = 1)+
  geom_point(data = ids, aes(x = Fidelity.Difference, y = PSS.Difference,  color = "Raw"), size = 2)+
  scale_color_discrete(name = "")+
  theme_gray(base_size = 30)+
  theme(
    panel.grid.major = element_line(size = 1.5)
    , panel.grid.minor = element_line(size = 1)
    , strip.background = element_blank()
    , strip.text.x = element_blank() 
  ) 

# plot with conditions
ggplot(data = kappadiff_vs_pssdiff, aes(x = kappadiff, y = pssdiff, color = factor(probedurationcondition), shape = factor(judgementtypecondition)))+ # color = factor(V1)))+
  geom_point(size = 4)+
  scale_x_continuous(name = "Log of Fidelity Differences")+
  scale_y_continuous(name = "Prior Entry Differences") +
  scale_color_discrete(name = "Probe Duration", labels = c("Short", "Long")) +
  scale_shape_discrete(name = "Judgement Type", labels = c("First", "Second"))+
  geom_vline(xintercept = 0, linetype = 2, size = 1)+
  geom_hline(yintercept = 0, linetype = 2, size = 1)+
  theme_gray(base_size = 30)+
  theme(
    panel.grid.major = element_line(size = 1.5)
    , panel.grid.minor = element_line(size = 1)
    , strip.background = element_blank()
    , strip.text.x = element_blank() 
  ) 
#------------------------- Prior Entry vs. Fidelity Differences ---------------------------#


#----------------------- Prior Entry vs. Fidelity EFFECTS ---------------------------------#
kappaeffect_ids = ddply(
  .data = betas
  , .variables = .(participant)
  , .fun = function(x){
    i = unique(x$participant)
    x_use = x[x$parameter ==  "population_log_kappa_attention_effect_mean",]$value
    kappaeffect = median(
      kappaattentioneffectmean + kappaattentioneffectsd*x_use 
    ) 
    df = data.frame(kappaeffect)
    names(df) = c("kappaeffect")
    return(df)
  }
)

# rid redundant Cs
kappaeffect_vs_priorentry2 = cbind(kappaeffect_ids, priorentry_ids[-1])

# shrinkage for prior entry effects
kappaeffect_vs_priorentry = cbind(kappaeffect_vs_priorentry2, ids[,c("PSS.Effect", "Fidelity.Effect")] )

# shrinkage for fidelity effects
plot(kappaeffect_vs_priorentry$participant, kappaeffect_vs_priorentry$Fidelity.Effect, col = "red")
abline(h = mean(kappaeffect_vs_priorentry$Fidelity.Effect), col = 'red')
points(kappaeffect_vs_priorentry$kappaeffect, col = "blue")
abline(h = mean(kappaeffect_vs_priorentry$kappaeffect), col = "blue")
# correlation between Bayesian estimate and raw value
cor(kappaeffect_vs_priorentry$kappaeffect, kappaeffect_vs_priorentry$Fidelity.Effect)

# see 2-D
ggplot(data = kappaeffect_vs_priorentry, aes(x = kappaeffect, y = priorentry))+
  geom_point(size = 2, aes(color = "Bayesian"))+
  scale_x_continuous(name = "Log of Fidelity Effects")+
  scale_y_continuous(name = "Prior Entry Effects") +
  geom_vline(xintercept = 0, linetype = 2, size = 1)+
  geom_hline(yintercept = 0, linetype = 2, size = 1)+
  geom_point(data = ids, aes(x = Fidelity.Effect, y = PSS.Effect,  color = "Raw"), size = 2)+
  scale_color_discrete(name = "")+
  theme_gray(base_size = 30)+
  theme(
    panel.grid.major = element_line(size = 1.5)
    , panel.grid.minor = element_line(size = 1)
    , strip.background = element_blank()
    , strip.text.x = element_blank() 
  ) 

# # plot with conditions
# ggplot(data = kappaeffect_vs_priorentry, aes(x = kappaeffect, y = priorentry, color = factor(probedurationcondition), shape = factor(judgementtypecondition)))+ # color = factor(V1)))+
#   geom_point(size = 4)+
#   scale_x_continuous(name = "Log of Fidelity Effects")+
#   scale_y_continuous(name = "Prior Entry Effects") +
#   scale_color_discrete(name = "Probe Duration", labels = c("Short", "Long")) +
#   scale_shape_discrete(name = "Judgement Type", labels = c("First", "Second"))+
#   geom_vline(xintercept = 0, linetype = 2, size = 1)+
#   geom_hline(yintercept = 0, linetype = 2, size = 1)+
#   theme_gray(base_size = 30)+
#   theme(
#     panel.grid.major = element_line(size = 1.5)
#     , panel.grid.minor = element_line(size = 1)
#     , strip.background = element_blank()
#     , strip.text.x = element_blank() 
#   ) 
# #----------------------- Prior Entry vs. Fidelity EFFECTS ---------------------------------#



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
  , c("Prior Entry Effect Mean")
  , y_lab = "PSS (Right - Left; ms)"
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

# NOTE: judgement_type effects will cancel out (because no back-transformation)
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

# NOTE: probe duration effects will cancel out (because no back-transformation)
get_violin(
  c(
      (ex_toj_color_post$population_pss_intercept_mean - ex_toj_color_post$population_pss_probe_duration_effect_mean/2
         + (ex_toj_color_post$population_pss_attention_effect_mean - ex_toj_color_post$population_pss_attention_probe_duration_interaction_effect_mean)/2 ) *250
    - (ex_toj_color_post$population_pss_intercept_mean - ex_toj_color_post$population_pss_probe_duration_effect_mean/2
       - (ex_toj_color_post$population_pss_attention_effect_mean - ex_toj_color_post$population_pss_attention_probe_duration_interaction_effect_mean)/2 ) *250
    
    , (ex_toj_color_post$population_pss_intercept_mean + ex_toj_color_post$population_pss_probe_duration_effect_mean/2
     + (ex_toj_color_post$population_pss_attention_effect_mean + ex_toj_color_post$population_pss_attention_probe_duration_interaction_effect_mean)/2 ) *250
    - (ex_toj_color_post$population_pss_intercept_mean + ex_toj_color_post$population_pss_probe_duration_effect_mean/2
       - (ex_toj_color_post$population_pss_attention_effect_mean + ex_toj_color_post$population_pss_attention_probe_duration_interaction_effect_mean)/2 ) * 250
  )
  , c("Short", "Long")
  , y_lab = "PSS (Right - Left; ms)"
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

# Back-transformed:
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

# Transformed
get_violin(
  ex_toj_color_post$population_log_jnd_intercept_mean
  , "JND Intercept Mean"
  , y_lab = "Log of SOA (ms)"
  , hline = FALSE
)
#---------------------------------- JND Intercept -----------------------------------------#


#---------------------------------- Main Effects ------------------------------------------#
# Back-transformed: effect of attention on JND
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

# Transformed effect of attention on JND
get_violin(
  ex_toj_color_post$population_log_jnd_attention_effect_mean
  , "JND\nAttention Effect Mean"
  , y_lab = "Log of SOA (Right - Left; ms)"
)

# Back-transformed: effect of judgement type (Q) on JND
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

# Transformed effect of judgement type (Q) on JND
get_violin(
  ex_toj_color_post$population_log_jnd_judgement_type_effect_mean
  , "JND Judgement\nType Effect Mean"
  , y_lab = "Log of SOA (Second - First; ms)"
)

# Back-transformed: effect of probe duration on JND
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

# Transformed effect of probe duration on JND
get_violin(
  ex_toj_color_post$population_log_jnd_probe_duration_effect_mean
  , "JND Probe\nDuration Effect Mean"
  , y_lab = "Log of SOA (Long - Short; ms)"
)
#---------------------------------- Main Effects ---------------------=--------------------#


#------------------------------- Two-way Interactions -------------------------------------#
#  Back-transformed: effect of interaction between judgement type and attention on JND 
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

#  Transformed effect of interaction between judgement type and attention on JND 
get_violin(
  2 * ex_toj_color_post$population_log_jnd_attention_judgement_type_interaction_effect_mean
  , "JND Attention\n& Judgement Type\nInteraction Effect Mean"
  , y_lab = "Log of SOA (ms)"  # awkward to label difference in difference direction in parentheses 
  # so explain in figure caption
)

#  Back-transformed: effect of interaction between probe duration and attention on JND
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

#  Transformed effect of interaction between probe duration and attention on JND 
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

# Back-transformed:
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

# with long condition only 
get_violin(
  ( rho_neg_attn_pos_probe
    + rho_pos_attn_pos_probe
  )/2
  , "Probability of Encoding\nLong Probe Duration\nIntercept Mean"
  , y_lab = "\u03C1"
  , hline = FALSE
)

# Transformed
get_violin(
  ex_toj_color_post$population_logit_rho_intercept_mean
  , "Probability of Encoding Intercept Mean"
  , y_lab = "Log-odds of \u03C1"
  , hline = FALSE
)
#---------------------------------- Rho Intercept -----------------------------------------#


#-------------------------------- Main Effects --------------------------------------------#
# Back-transformed: attention effect
get_violin(
  ( rho_pos_attn_pos_probe - rho_neg_attn_pos_probe
  + rho_pos_attn_neg_probe - rho_neg_attn_neg_probe 
  )/2
  , "Probability of Encoding Attention Effect Mean"
  , y_lab = "\u03C1 (Attended - Unattended)"
)

# Transformed attention effect
get_violin(
  ex_toj_color_post$population_logit_rho_attention_effect_mean
  , "Probability of Encoding Attention Effect Mean"
  , y_lab = "Log-odds of \u03C1 (Attended - Unattended)"
)

# Back-transformed: probe effect
get_violin(
  ( rho_pos_attn_pos_probe - rho_pos_attn_neg_probe
  + rho_neg_attn_pos_probe - rho_neg_attn_neg_probe
  )/2
  , "Probability of Encoding\nProbe Duration Effect Mean"
  , y_lab = "\u03C1 (Long - Short)"
)

#  Transformed probe effect
get_violin(
  ex_toj_color_post$population_logit_rho_probe_duration_effect_mean
  , "Probability of Encoding\nProbe Duration Effect Mean"
  , y_lab = "Log-odds of \u03C1 (Long - Short)"
)
#-------------------------------- Main Effects --------------------------------------------#


#------------------------------- Two-way Interactions -------------------------------------#
# Back-transformed:
get_violin(
  (rho_pos_attn_pos_probe - rho_neg_attn_pos_probe) - (rho_pos_attn_neg_probe - rho_neg_attn_neg_probe)
  , "Probability of Encoding\nProbe Duration Attention\nInteraction Effect Mean"
  , y_lab = "\u03C1"
)

# Transformed
get_violin(
  2 * ex_toj_color_post$population_logit_rho_attention_probe_duration_interaction_effect_mean
  , "Probability of Encoding\nProbe Duration Attention\nInteraction Effect Mean"
  , y_lab = "Log-odds of \u03C1"
)

# Back-transformed:
values = c( rho_pos_attn_neg_probe - rho_neg_attn_neg_probe
    , rho_pos_attn_pos_probe - rho_neg_attn_pos_probe )
labels = c("Short", "Long")
y_lab = "\u03C1 (Attended - Unattended)"
samps = iters
hline = T

label_vec = NULL
for (i in 1:length(labels)) {
  label = rep(labels[i], samps)
  label_vec = c(label_vec, label)
}

df = data.frame(  
  value = values
  , parameter = label_vec
)

df$parameter = factor(df$parameter, labels)

gg = ggplot(data = df)+
  geom_violin(aes(x = parameter, y = value))+
  labs(x = "", y = y_lab)+
  stat_summary(aes(x = parameter, y = value), fun.data = get_95_HDI, size = 0.7)+  
  stat_summary(aes(x = parameter, y = value), fun.data = get_50_HDI, size = 1.7)

if (hline) {
  gg = gg + geom_hline(yintercept = 0, linetype = 2, size = 1)
}

if (facet) {
  gg = gg + facet_wrap(~parameter, scales = "free")
} 

gg = gg + theme_gray(base_size = 30)+
  theme(
    panel.grid.major = element_line(size = 1.5)
    , panel.grid.minor = element_line(size = 1)
    , strip.background = element_blank()
    , strip.text.x = element_blank() 
    , axis.ticks.x = element_blank() 
  ) 

print(gg)


# Transformed 
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

# Back-transformed:
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

# Transformed:
get_violin(
  ex_toj_color_post$population_log_kappa_intercept_mean
  , "Fidelity of Memory Intercept Mean"
  , y_lab = "Log of \u03BA"
  , hline = FALSE
)
#---------------------------------- Kappa Intercept ---------------------------------------#


#-------------------------------- Main Effects --------------------------------------------#
# Back-transformed: attention effect
get_violin(
  ( kappa_pos_attn_pos_probe - kappa_neg_attn_pos_probe
  + kappa_pos_attn_neg_probe - kappa_neg_attn_neg_probe 
  )/2
  , "Fidelity of Encoding Attention Effect Mean"
  , y_lab = "\u03BA (Attended - Unattended)"
)

# Transformed attention effect
get_violin(
  ex_toj_color_post$population_logit_rho_attention_effect_mean
  , "Fidelity of Encoding Attention Effect Mean"
  , y_lab = "Log of \u03BA (Attended - Unattended)"
)

# Back-transformed: probe effect
get_violin(
  ( kappa_pos_attn_pos_probe - kappa_pos_attn_neg_probe
  + kappa_neg_attn_pos_probe - kappa_neg_attn_neg_probe
  )/2
  , "Fidelity of Encoding Probe Duration Effect Mean"
  , y_lab = "\u03BA (Long - Short)"
)

# Transformed probe effect
get_violin(
  ex_toj_color_post$population_log_kappa_probe_duration_effect_mean
  , "Fidelity of Encoding Probe Duration Effect Mean"
  , y_lab = "Log of \u03BA (Long - Short)"
)
#-------------------------------- Main Effects --------------------------------------------#


#------------------------------- Two-way Interactions -------------------------------------#
# Back-transformed:
get_violin(
  (kappa_pos_attn_pos_probe - kappa_neg_attn_pos_probe) - (kappa_pos_attn_neg_probe - kappa_neg_attn_neg_probe)
  , "Fidelity of Encoding\nProbe Duration Attention\nInteraction Effect Mean"
  , y_lab = "\u03BA"
)

# Transformed
get_violin(
  2 * ex_toj_color_post$population_log_kappa_attention_probe_duration_interaction_effect_mean
  , "Fidelity of Encoding\nProbe Duration Attention\nInteraction Effect Mean"
  , y_lab = "Log of \u03BA"
)

# simpler
get_violin(
  c(kappa_pos_attn_neg_probe - kappa_neg_attn_neg_probe
    , kappa_pos_attn_pos_probe - kappa_neg_attn_pos_probe)
  , c("Short", "Long")
  , y_lab = "\u03BA (Attended - Unattended)"
)

# Transformed
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
yLeftFirstLong = pnorm(
  -250:250
  , mean = ( 
      median(
      (
        ( ex_toj_color_post$population_pss_intercept_mean - ex_toj_color_post$population_pss_judgement_type_effect_mean/2 + ex_toj_color_post$population_pss_probe_duration_effect_mean/2
        - ( ex_toj_color_post$population_pss_attention_effect_mean - ex_toj_color_post$population_pss_attention_judgement_type_interaction_effect_mean + ex_toj_color_post$population_pss_attention_probe_duration_interaction_effect_mean)/2)
      ) * 250 
    )
  )
  , sd = median(
     (neg_attn_neg_judge_pos_probe) * 250
    )
)

yLeftFirstShort = pnorm(
  -250:250
  , mean = ( 
    median(
      (
        ( ex_toj_color_post$population_pss_intercept_mean - ex_toj_color_post$population_pss_judgement_type_effect_mean/2 - ex_toj_color_post$population_pss_probe_duration_effect_mean/2
            - ( ex_toj_color_post$population_pss_attention_effect_mean - ex_toj_color_post$population_pss_attention_judgement_type_interaction_effect_mean - ex_toj_color_post$population_pss_attention_probe_duration_interaction_effect_mean)/2)
      ) * 250 
    )
  )
  , sd = median(
    (neg_attn_neg_judge_neg_probe) * 250
  )
)

yLeftFirst = (yLeftFirstShort + yLeftFirstLong)/2

yLeftSecondLong = pnorm(
  -250:250
  , mean = ( 
      median(
      (
        ( ex_toj_color_post$population_pss_intercept_mean + ex_toj_color_post$population_pss_judgement_type_effect_mean/2 + ex_toj_color_post$population_pss_probe_duration_effect_mean/2
          - ( ex_toj_color_post$population_pss_attention_effect_mean + ex_toj_color_post$population_pss_attention_judgement_type_interaction_effect_mean + ex_toj_color_post$population_pss_attention_probe_duration_interaction_effect_mean)/2)
      ) * 250 
    )
  )
  , sd = median(
    (neg_attn_pos_judge_pos_probe) * 250
    ) 
)

yLeftSecondShort = pnorm(
  -250:250
  , mean = ( 
    median(
      (
        ( ex_toj_color_post$population_pss_intercept_mean + ex_toj_color_post$population_pss_judgement_type_effect_mean/2 - ex_toj_color_post$population_pss_probe_duration_effect_mean/2
            - ( ex_toj_color_post$population_pss_attention_effect_mean + ex_toj_color_post$population_pss_attention_judgement_type_interaction_effect_mean - ex_toj_color_post$population_pss_attention_probe_duration_interaction_effect_mean)/2)
      )* 250 
    )
  )
  , sd = median(
    (neg_attn_pos_judge_neg_probe) * 250
  ) 
)

yLeftSecond = (yLeftSecondShort + yLeftSecondLong)/2

yRightFirstLong = pnorm(
  -250:250
  , mean = ( 
      median(
      (
        ( ex_toj_color_post$population_pss_intercept_mean - ex_toj_color_post$population_pss_judgement_type_effect_mean/2 + ex_toj_color_post$population_pss_probe_duration_effect_mean/2
          + ( ex_toj_color_post$population_pss_attention_effect_mean - ex_toj_color_post$population_pss_attention_judgement_type_interaction_effect_mean + ex_toj_color_post$population_pss_attention_probe_duration_interaction_effect_mean)/2)
       ) * 250 
    )
  )
  , sd = median(
    (pos_attn_neg_judge_pos_probe) * 250
  )
)

yRightFirstShort = pnorm(
  -250:250
  , mean = ( 
    median(
      (
        ( ex_toj_color_post$population_pss_intercept_mean - ex_toj_color_post$population_pss_judgement_type_effect_mean/2 - ex_toj_color_post$population_pss_probe_duration_effect_mean/2
            + ( ex_toj_color_post$population_pss_attention_effect_mean - ex_toj_color_post$population_pss_attention_judgement_type_interaction_effect_mean - ex_toj_color_post$population_pss_attention_probe_duration_interaction_effect_mean)/2)
      ) * 250 
    )
  )
  , sd = median(
    (pos_attn_neg_judge_neg_probe) * 250
  )
)

yRightFirst = (yRightFirstShort + yRightFirstLong)/2

yRightSecondLong = pnorm(
  -250:250
  , mean = ( 
    median(
      (
        ( ex_toj_color_post$population_pss_intercept_mean + ex_toj_color_post$population_pss_judgement_type_effect_mean/2 + ex_toj_color_post$population_pss_probe_duration_effect_mean/2
          + ( ex_toj_color_post$population_pss_attention_effect_mean + ex_toj_color_post$population_pss_attention_judgement_type_interaction_effect_mean + ex_toj_color_post$population_pss_attention_probe_duration_interaction_effect_mean)/2)
      ) * 250 
    )
  )
  , sd =  median( 
    (pos_attn_pos_judge_pos_probe) * 250
  )
)

yRightSecondShort = pnorm(
  -250:250
  , mean = ( 
    median(
      (
       ( ex_toj_color_post$population_pss_intercept_mean + ex_toj_color_post$population_pss_judgement_type_effect_mean/2 - ex_toj_color_post$population_pss_probe_duration_effect_mean/2
            + ( ex_toj_color_post$population_pss_attention_effect_mean + ex_toj_color_post$population_pss_attention_judgement_type_interaction_effect_mean - ex_toj_color_post$population_pss_attention_probe_duration_interaction_effect_mean)/2)
      ) * 250 
    )
  )
  , sd =  median( 
    (pos_attn_pos_judge_neg_probe) * 250
  )
)

yRightSecond = (yRightSecondShort + yRightSecondLong)/2

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


### Just attention conditions
yRight = (yRightFirst + yRightSecond)/2
yLeft = (yLeftFirst + yLeftSecond)/2

df = data.frame(SOA = -250:250
                , Prop = c(yRight, yLeft)
                , Attend = c(rep("Right",501), rep("Left", 501))
)

toj_means_by_id_by_condition3 = ddply(
  .data = toj_trials
  , .variables = .(block_bias, soa2)
  , .fun = function(x){
    to_return = data.frame(
      Prop = mean(x$left_first_TF)
      , Attend = paste( unique(x$block_bias) )
      , SOA = unique(x$soa2)
    )
    return(to_return)
  }
)

levels(toj_means_by_id_by_condition3$Attend) = c("Left", "Right")

gg = ggplot(data = df, aes(y = Prop, x = SOA, colour = Attend))+
  geom_line(size = 1.25)+
  # scale_color_manual("Attend", values = c("red", "blue"))+
  scale_color_hue("Attend", l = c(60, 15), c = c(100, 50), h = c(240, 360) ) +
  labs(x = "SOA (ms)", y = "Proportion of 'Left First' Responses")+
  geom_point(data = toj_means_by_id_by_condition3, aes(y = Prop, x = SOA, colour = Attend), size = 4)+
  geom_point(data = toj_means_by_id_by_condition3, aes(y = Prop, x = SOA), size = 2.5, colour = "grey90")+
  theme_gray(base_size = 30)+
  theme(panel.grid.major = element_line(size = 1.5)
        ,panel.grid.minor = element_line(size = 1))
# define text to add
Text1 = textGrob(label = paste("Right First"), gp = gpar(fontsize= 30))
Text2 = textGrob(label = paste("Left First"), gp = gpar(fontsize= 30)) 
gg = gg+
  annotation_custom(grob = Text1,  xmin = -225, xmax = -225, ymin = -0.14, ymax = -0.14)+
  annotation_custom(grob = Text2,  xmin = 225, xmax = 225, ymin = -0.14, ymax = -0.14)
# Code to override clipping
gg2 <- ggplot_gtable(ggplot_build(gg))
gg2$layout$clip[gg2$layout$name=="panel"] <- "off"
grid.draw(gg2)
#---------------------------------- NCFs --------------------------------------------------#






