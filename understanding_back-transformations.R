# I likely don't need all of these, but for the functions...
library(coda)
library(ggplot2)
library(ggmcmc)
library(CircStats)
library(reshape2)
library(plyr)
library(grid)

# to get get_violin function
source('/Users/ghislaindentremont/Documents/TOJ/EndogenousVisualPriorEntry-BayesianHierarchicalModel/functions.R')


#### See curves ####
# curve(plogis(x), -5 ,5)
# Rho
curve(qlogis(x), 0,1)
# Kappa
curve(log(x), 0, 20)
# JND
curve(log(x), 0, 1)



#### back-transformation of posterior from logit space (rho) ####
# input intercept and effect in prob space
logit_rho_effect_backtransformation = function(intercept, effect, quasi_SD){
  logit_intercept = qlogis(intercept)
  logit_effect = qlogis(intercept+effect/2) - qlogis(intercept - effect/2)
  logit_SD = ( qlogis(intercept + quasi_SD/2) - qlogis(intercept - quasi_SD/2) )/2
  
  v1 = rnorm(80000, logit_intercept+logit_effect/2, logit_SD)
  v2= rnorm(80000, logit_intercept-logit_effect/2, logit_SD)
  
  get_violin(
    c(
      v1 - v2
      , plogis(v1) - plogis(v2)
    )
    , c("logit space", "prob space")
    , sprintf("rho (%s)", intercept)
    , facet = TRUE
  )
  
}



#### from normal log JND distribution to JND scale ####
log_par_dist_to_par_dist = function(parmean, parsd, is_jnd = T) {
  if (is_jnd) {
    logparmean = log(parmean/250)
    logparsd =  log((parmean+parsd)/250) - logparmean
  } else if (!is_jnd) {
    logparmean = log(parmean)
    logparsd =  log(parmean+parsd) - logparmean
  }
  logpar = rnorm(10000, logparmean, logparsd)
  hist(logpar)
  
  if (is_jnd) {
    par = exp(logpar)*250 
  } else if (!is_jnd) {
    par = exp(logpar)
  }
  hist(par) 
}

log_par_dist_to_par_dist(100, 20)
log_par_dist_to_par_dist(60, 10)



#### from normal log kappa distribution to kappa scale ####
log_par_dist_to_par_dist(12, 2, is_jnd = F)
log_par_dist_to_par_dist(6, 2, is_jnd = F)



#### from normal logit rho distribution to prob scale ####
logit_rho_dist_to_rho_dist = function(rhomean, rhosd) {
  logitrhomean = qlogis(rhomean)
  logitrhosd =  qlogis(rhomean+rhosd) - logitrhomean
  logitrho = rnorm(10000, logitrhomean, logitrhosd)
  hist(logitrho)
  
  rho = plogis(logitrho)
  hist(rho) 
}

logit_rho_dist_to_rho_dist(.95, .02)
logit_rho_dist_to_rho_dist(.8, .04)
