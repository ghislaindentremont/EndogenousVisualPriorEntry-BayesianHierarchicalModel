#### Libraries ####
library(plyr)
library(ggplot2)
library(grid)
library(rstan)
library(corrplot)

setwd("/Users/ghislaindentremont/Documents/TOJ/Follow-Up")

# read in cleaned-up data structures 
load("FollowUp_color_trials.Rdata")
load("FollowUp_toj_trials.Rdata")


###############################################################################
####                        Statistical Analyses                           ####
###############################################################################

color_means = aggregate(abs_color_diff ~ attended + longprobe + id , data = color_trials, FUN = mean)

#-----------------------------------------------------------------------------#
#                            Absolute Color Diff                              #
#-----------------------------------------------------------------------------#


#----------------------------------- Degrees ---------------------------------#
### with convention knowledge 
color_know = aggregate(abs_color_diff ~ attended + longprobe + id, data=color_trials, FUN = mean)
color_know$longprobe = as.factor(as.character(color_know$longprobe))

# descriptive 
color_know_means = aggregate(abs_color_diff ~ attended + longprobe, data = color_know, FUN = mean)
color_know_SD =  aggregate(abs_color_diff ~ attended + longprobe, data = color_know, FUN = sd)
#----------------------------------- Degrees ---------------------------------#


#----------------------------------- Log Degrees -----------------------------#
color_means$log_abs_color_diff = log(color_means$abs_color_diff)
# hist(color_means$log_abs_color_diff)

color_log_know = aggregate(log_abs_color_diff ~ attended + longprobe + id, data=color_means, FUN = mean)
color_log_know$longprobe = as.factor(as.character(color_log_know$longprobe))

# descriptive 
color_log_know_means = aggregate(log_abs_color_diff ~ attended + longprobe, data = color_log_know, FUN = mean)
color_log_know_SD =  aggregate(log_abs_color_diff ~ attended + longprobe, data = color_log_know, FUN = sd)
#----------------------------------- Log Degrees -----------------------------#


#-----------------------------------------------------------------------------#
#                             Mixture Model                                   #
#-----------------------------------------------------------------------------#

source("/Users/ghislaindentremont/Documents/TOJ/EndogenousVisualPriorEntry-BayesianHierarchicalModel/Baseball/conventional_analysis/fit_uvm.R")

fitted = ddply(
    .data = color_trials
    , .variables = .(id, attended, longprobe)
    , .fun = function(piece_of_df){
      fit = fit_uvm(piece_of_df$color_diff_radians, do_mu = TRUE)
      if (fit$rho == 1) {
        n = length(piece_of_df$color_diff_radians)
        fit$rho = (2*n-1)/(2*n)
      }
      to_return = data.frame(
        kappa_prime = fit$kappa_prime
        , kappa = exp(fit$kappa_prime)
        , rho = fit$rho
        , logit_rho = qlogis(fit$rho)
      )
      return(to_return)
    }
  , .progress = 'time'
) 

# just by id  
fitted_id = ddply(
  .data = color_trials
  , .variables = .(id)
  , .fun = function(piece_of_df){
    fit = fit_uvm(piece_of_df$color_diff_radians, do_mu = TRUE)
    if (fit$rho == 1) {
      n = length(piece_of_df$color_diff_radians)
      fit$rho = (2*n-1)/(2*n)
    }
    to_return = data.frame(
      kappa_prime = fit$kappa_prime
      , kappa = exp(fit$kappa_prime)
      , rho = fit$rho
      , logit_rho = qlogis(fit$rho)
    )
    return(to_return)
  }
  , .progress = 'time'
)


#----------------------- Fidelity (kappa prime) ------------------------------#
color_kappa_prime = aggregate(kappa_prime ~ attended + longprobe + id, data=fitted, FUN = mean)
color_kappa_prime$longprobe = as.factor(as.character(color_kappa_prime$longprobe))

# descriptve 
color_kappa_prime_means = aggregate(kappa_prime ~ attended + longprobe, data = color_kappa_prime, FUN = mean)
color_kappa_prime_SD =  aggregate(kappa_prime ~ attended + longprobe, data = color_kappa_prime, FUN = sd)
#----------------------- Fidelity (kappa prime) ------------------------------#



#------------------------- Fidelity (kappa) ----------------------------------#
color_kappa = aggregate(kappa ~ attended + longprobe + id, data=fitted, FUN = mean)
color_kappa$longprobe = as.factor(as.character(color_kappa$longprobe))

# descriptve 
color_kappa_means = aggregate(kappa ~ attended + longprobe, data = color_kappa, FUN = mean)
color_kappa_SD =  aggregate(kappa ~ attended + longprobe, data = color_kappa, FUN = sd)
#------------------------- Fidelity (kappa) ----------------------------------#


#------------------------ Probability (rho) ----------------------------------# 
color_rho = aggregate(rho ~ attended + longprobe + id, data=fitted, FUN = mean)
color_rho$longprobe = as.factor(as.character(color_rho$longprobe))

# descriptive 
color_rho_means = aggregate(rho ~ attended + longprobe, data = color_rho, FUN = mean)
color_rho_SD =  aggregate(rho ~ attended + longprobe, data = color_rho, FUN = sd)
#------------------------ Probability (rho) ----------------------------------# 


#------------------------ Prob. (logit rho) ----------------------------------#
color_logit_rho = aggregate(logit_rho ~ attended + longprobe + id, data=fitted, FUN = mean)
color_logit_rho$longprobe = as.factor(as.character(color_logit_rho$longprobe))

# descriptive 
color_logit_rho_means = aggregate(logit_rho ~ attended + longprobe, data = color_logit_rho, FUN = mean)
color_logit_rho_SD =  aggregate(logit_rho ~ attended + longprobe, data = color_logit_rho, FUN = sd)
#------------------------ Prob. (logit rho) ----------------------------------# 



#-----------------------------------------------------------------------------#
#                                 TOJ                                         #
#-----------------------------------------------------------------------------#
# get pss and jnd
id_all = unique(toj_trials$id)
get_pss_jnd = function(id_list, norma = F) {
  if (norma) {
    toj_trials$soa_use = toj_trials$soa3
  } else {
    toj_trials$soa_use = toj_trials$soa2
  }
  toj_by_condition = ddply(
    .data = toj_trials[toj_trials$id %in% id_list,]
    , .variables = .(id, block_bias, longprobe, toj_judgement_type)
    , .fun = function(x){
      fit = glm(
        formula = left_first_TF~soa_use
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

toj_by_condition = get_pss_jnd(id_all, norma = T)

toj_by_condition$pss = toj_by_condition$pss * 250
toj_by_condition$log_jnd = log(toj_by_condition$jnd)
toj_by_condition$jnd = toj_by_condition$jnd * 250
 

#------------------------------- PSS -----------------------------------------# 
toj_pss = aggregate(pss ~ block_bias + longprobe + toj_judgement_type + id, data=toj_by_condition, FUN = mean)
toj_pss$longprobe = as.factor(as.character(toj_pss$longprobe))

# descriptive 
toj_pss_means = aggregate(pss ~ block_bias + longprobe + toj_judgement_type, data = toj_pss, FUN = mean)
toj_pss_SD =  aggregate(pss ~ block_bias + longprobe + toj_judgement_type, data = toj_pss, FUN = sd)
#------------------------------- PSS -----------------------------------------# 


#-------------------------------- JND ----------------------------------------# 
toj_jnd = aggregate(jnd ~ block_bias + longprobe + toj_judgement_type + id, data=toj_by_condition, FUN = mean)
toj_jnd$longprobe = as.factor(as.character(toj_jnd$longprobe))
toj_jnd$toj_judgement_type = as.factor(as.character(toj_jnd$toj_judgement_type))

# descriptive 
toj_jnd_means = aggregate(jnd ~ block_bias + longprobe + toj_judgement_type, data = toj_jnd, FUN = mean)
toj_jnd_SD =  aggregate(jnd ~ block_bias + longprobe + toj_judgement_type, data = toj_jnd, FUN = sd)
#-------------------------------- JND ----------------------------------------# 


#-------------------------------- log JND ------------------------------------# 
toj_log_jnd = aggregate(log_jnd ~ block_bias + longprobe + toj_judgement_type + id, data=toj_by_condition, FUN = mean)
toj_log_jnd$longprobe = as.factor(as.character(toj_log_jnd$longprobe))

# descriptive 
toj_log_jnd_means = aggregate(log_jnd ~ block_bias + longprobe + toj_judgement_type, data = toj_log_jnd, FUN = mean)
toj_log_jnd_SD =  aggregate(log_jnd ~ block_bias + longprobe + toj_judgement_type, data = toj_log_jnd, FUN = sd)
#-------------------------------- log JND ------------------------------------# 


#----------------------- Correlation Matrix ----------------------------------# 
ids_temp = data.frame(
  id = aggregate(pss ~id, data = toj_by_condition, FUN = mean)$id
  , "PSS Mean" = aggregate(pss ~ id, data = toj_by_condition, FUN = mean)$pss
  , "JND Mean" = aggregate(log_jnd ~ id, data = toj_by_condition, FUN = mean)$log_jnd
  , "Probability Mean" = fitted_id$logit_rho
  , "Fidelity Mean" = fitted_id$kappa_prime
  # natural right - left
  , "PSS Difference" = aggregate(pss ~ id, data = toj_by_condition, FUN = diff)$pss
  , "JND Difference" = aggregate(log_jnd ~ id, data = toj_by_condition, FUN = diff)$log_jnd
  # this is attended - unattended, naturally
  , "Probability Difference" = aggregate(logit_rho ~ id, data = fitted, FUN = diff)$logit_rho
  , "Fidelity Difference" = aggregate(kappa_prime ~ id, data = fitted, FUN = diff)$kappa_prime
)

### get matrix
corrs = ids_temp[,-1]
M = cor(corrs)

# plot
corrplot(M, method = "number", type = "lower")
#----------------------- Correlation Matrix ----------------------------------# 


#----------------------- Isolate Interactions --------------------------------# 
cdns = aggregate(pss ~longprobe + toj_judgement_type + id, data = toj_by_condition, FUN = mean)[,c("id", "longprobe", "toj_judgement_type")]

cdns$longprobe = ifelse(aggregate(longprobe ~ id, data = cdns, FUN = unique)$longprobe == "FALSE", -1, 1)   # if longprobe T, then +1
cdns$toj_judgement_type = ifelse(aggregate(toj_judgement_type ~ id, data = cdns, FUN = unique)$toj_judgement_type == "first", -1, 1) 

ids2 = cbind(cdns, ids_temp)
ids = ids2[-4]  # double checked allignment, now remove extra id column





### PSS intercepts 
# probe duration
pss_probe_effect_temp = aggregate(pss ~ longprobe, data = toj_pss_means, FUN = mean)
# 'Long' minus 'Short'
pss_probe_effect = pss_probe_effect_temp$pss[2] - pss_probe_effect_temp$pss[1]

# judgement type
pss_judgement_effect_temp = aggregate(pss ~toj_judgement_type, data = toj_pss_means, FUN = mean)
# 'second' minus 'first'
pss_judgement_effect = pss_judgement_effect_temp$pss[2] - pss_judgement_effect_temp$pss[1]

# need to devide effect by two before subtracting out
# THINK ABOUT IT: you can get at full difference in differences because you only have one between subject group per P
ids$PSS.Intercept = ids$PSS.Mean - (pss_probe_effect*ids$longprobe + pss_judgement_effect*ids$toj_judgement_type)/2


### JND intercepts
# probe duration
log_jnd_probe_effect_temp = aggregate(log_jnd ~ longprobe, data = toj_log_jnd_means, FUN = mean)
# 'Long' minus 'Short'
log_jnd_probe_effect = log_jnd_probe_effect_temp$log_jnd[2] - log_jnd_probe_effect_temp$log_jnd[1]

# judgement type
log_jnd_judgement_effect_temp = aggregate(log_jnd ~toj_judgement_type, data = toj_log_jnd_means, FUN = mean)
# 'second' minus 'first'
log_jnd_judgement_effect = log_jnd_judgement_effect_temp$log_jnd[2] - log_jnd_judgement_effect_temp$log_jnd[1]

# need to devide effect by two before subtracting out
# THINK ABOUT IT: you can get at full difference in differences because you only have one between subject group per P
ids$JND.Intercept = ids$JND.Mean - (log_jnd_probe_effect*ids$longprobe + log_jnd_judgement_effect*ids$toj_judgement_type)/2


### Probability intercepts
rho_probe_effect_temp = aggregate(logit_rho ~ longprobe, data = fitted, FUN = mean)
# 'Long' minus 'Short'
rho_probe_effect = rho_probe_effect_temp$logit_rho[2] - rho_probe_effect_temp$logit_rho[1]

# need to devide effect by two before subtracting out
ids$Probability.Intercept = ids$Probability.Mean - (rho_probe_effect*ids$longprobe)/2


### Fidelity intercepts
kappa_probe_effect_temp = aggregate(kappa_prime ~ longprobe, data = fitted, FUN = mean)
# 'Long' minus 'Short'
kappa_probe_effect = kappa_probe_effect_temp$kappa_prime[2] - kappa_probe_effect_temp$kappa_prime[1]

# need to devide effect by two before subtracting out
ids$Fidelity.Intercept = ids$Fidelity.Mean- (kappa_probe_effect*ids$longprobe)/2



### PSS effects
for_pss_interactions = aggregate(pss ~ toj_judgement_type + longprobe, data = toj_pss_means, FUN = diff)

# probe duration
pss_probe_interaction_temp = aggregate(pss ~ longprobe, data = for_pss_interactions, FUN = mean)
# 'Long' minus 'Short'
pss_probe_interaction = pss_probe_interaction_temp$pss[2] - pss_probe_interaction_temp$pss[1]

# judgement type
pss_judgement_interaction_temp = aggregate(pss ~toj_judgement_type, data = for_pss_interactions, FUN = mean)
# 'second' minus 'first'
pss_judgement_interaction = pss_judgement_interaction_temp$pss[2] - pss_judgement_interaction_temp$pss[1]

# need to devide interaction by two before subtracting out
# THINK ABOUT IT: you can get at full difference in differences because you only have one between subject group per P
ids$PSS.Effect = ids$PSS.Difference - (pss_probe_interaction*ids$longprobe + pss_judgement_interaction*ids$toj_judgement_type)/2


### JND effects
for_log_jnd_interactions = aggregate(log_jnd ~ toj_judgement_type + longprobe, data = toj_log_jnd_means, FUN = diff)

# probe duration
log_jnd_probe_interaction_temp = aggregate(log_jnd ~ longprobe, data = for_log_jnd_interactions, FUN = mean)
# 'Long' minus 'Short'
log_jnd_probe_interaction = log_jnd_probe_interaction_temp$log_jnd[2] - log_jnd_probe_interaction_temp$log_jnd[1]

# judgement type
log_jnd_judgement_interaction_temp = aggregate(log_jnd ~toj_judgement_type, data = for_log_jnd_interactions, FUN = mean)
# 'second' minus 'first'
log_jnd_judgement_interaction = log_jnd_judgement_interaction_temp$log_jnd[2] - log_jnd_judgement_interaction_temp$log_jnd[1]

# need to devide interaction by two before subtracting out
# THINK ABOUT IT: you can get at full difference in differences because you only have one between subject group per P
ids$JND.Effect = ids$JND.Difference - (log_jnd_probe_interaction*ids$longprobe + log_jnd_judgement_interaction*ids$toj_judgement_type)/2


### Probability effects
for_rho_interactions2 = aggregate(logit_rho ~ longprobe + attended, data = fitted, FUN = mean)
rho_probe_interaction_temp = aggregate(logit_rho ~ longprobe, data = for_rho_interactions2, FUN = diff)
# 'Long' minus 'Short'
rho_probe_interaction = rho_probe_interaction_temp$logit_rho[2] - rho_probe_interaction_temp$logit_rho[1]

# need to devide interaction by two before subtracting out
ids$Probability.Effect = ids$Probability.Difference - (rho_probe_interaction*ids$longprobe)/2


### Fidelity effects
for_kappa_interactions2 = aggregate(kappa_prime ~ longprobe + attended, data = fitted, FUN = mean)
kappa_probe_interaction_temp = aggregate(kappa_prime ~ longprobe, data = for_kappa_interactions2, FUN = diff)
# 'Long' minus 'Short'
kappa_probe_interaction = kappa_probe_interaction_temp$kappa_prime[2] - kappa_probe_interaction_temp$kappa_prime[1]

# need to devide interaction by two before subtracting out
ids$Fidelity.Effect = ids$Fidelity.Difference - (kappa_probe_interaction*ids$longprobe)/2


### Corr Matrix Again
corrse = ids[c(12:19)]
Me = cor(corrse)

# plot
corrplot(Me, method = "number", type = "lower")
#----------------------- Isolate Interactions --------------------------------# 