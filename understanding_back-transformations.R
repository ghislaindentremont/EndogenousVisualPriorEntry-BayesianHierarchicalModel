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

#### logit space (rho) ####
curve(plogis(x), -5 ,5)

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


