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

# look structure
str(gg_toj_color_post)

# list of parameters to examine
param_list = c("logitRhoEffectMean"
               , "logitRhoMean"
               , "logitRhoConventionEffectMean"
               , "logitRhoConventionInteractionEffectMean"
               , "logKappaEffectMean"
               , "logKappaConventionEffectMean"
               , "logKappaConventionInteractionEffectMean"
               , "logKappaMean"
               , "population_logjnd_effect_mean"
               , "population_logjnd_convention_effect_mean"
               , "population_logjnd_intercept_mean"
               , "population_logjnd_convention_interaction_effect_mean"
               , "population_pss_effect_mean"
               , "population_pss_convention_effect_mean" 
               , "population_pss_intercept_mean"
               , "population_pss_convention_interaction_effect_mean"
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


#-------------------------------------- Color Actual Data ---------------------------------#
hist(color_trials[color_trials$attended == TRUE,]$color_diff_radians, breaks = 30, freq = F, col = rgb(.1,.1,.1,.5))
hist(color_trials[color_trials$attended == FALSE,]$color_diff_radians, breaks = 30, freq = F, col = rgb(.9,.9,.9,.5), add = T)
#-------------------------------------- Color Actual Data ---------------------------------#


#-------------------------------------- Color Simulated Data ------------------------------#
### Get Rho Parameters 
rho_intercept_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "logitRhoMean",]$value
rho_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "logitRhoEffectMean",]$value
rho_convention_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "logitRhoConventionEffectMean",]$value
rho_convention_interaction_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "logitRhoConventionInteractionEffectMean",]$value

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
kappa_intercept_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "logKappaMean",]$value
kappa_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "logKappaEffectMean",]$value
kappa_convention_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "logKappaConventionEffectMean",]$value
kappa_convention_interaction_effect_mean = gg_toj_color_post[gg_toj_color_post$Parameter == "logKappaConventionInteractionEffectMean",]$value

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

# version so nothing gets flipped twice by accident 
ex_toj_color_post2 = extract(toj_color_post)

# flip pss and jnd effects to make positive values indicate predicted results
ex_toj_color_post = ex_toj_color_post2

ex_toj_color_post$population_pss_effect_mean = -ex_toj_color_post2$population_pss_effect_mean
ex_toj_color_post$population_pss_convention_interaction_effect_mean = -ex_toj_color_post2$population_pss_convention_interaction_effect_mean

ex_toj_color_post$population_logjnd_effect_mean = -ex_toj_color_post2$population_logjnd_effect_mean
ex_toj_color_post$population_logjnd_convention_interaction_effect_mean = -ex_toj_color_post2$population_logjnd_convention_interaction_effect_mean


# get correlation posteriors 
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
  , "logitRhoMean"                       
  , "logKappaMean"                       
  , "logitRhoEffectMean"                 
  , "logKappaEffectMean"
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


# #---------------------------- Rho vs. PSS Effects -----------------------------------------#
# psseffect = extract_samples("population_pss_effect_mean")
# 
# psseffectsd = extract_samples("zpopulation_pss_effect_sd", TRUE)
# 
# pssinteractioneffect = extract_samples("population_pss_convention_interaction_effect_mean")
# 
# conventionfactor = ifelse(aggregate(know_tie_goes_runner~id,data = toj_trials, FUN =unique)$know_tie_goes_runner, -1, 1)
# 
# psseffect_ids = ddply(
#   .data = betas
#   , .variables = .(participant)
#   , .fun = function(x){
#     i = unique(x$participant)
#     x_use = x[x$parameter ==  "population_pss_effect_mean",]$value
#     psseffect = median(psseffect)  + median(psseffectsd)*median(x_use) 
#     + median(pssinteractioneffect)*conventionfactor[i]
#     df = data.frame(psseffect*250, conventionfactor[i])
#     names(df) = c("psseffect", "conventionfactor")
#     return(df)
#   }
# )
# 
# logitrhoeffect = extract_samples("logitRhoEffectMean")
# 
# logitrhoeffectsd = extract_samples("zlogitRhoEffectSD", TRUE)
# 
# logitrhointeractioneffect = extract_samples("logitRhoConventionInteractionEffectMean")
# 
# rhoeffect_ids = ddply(
#   .data = betas
#   , .variables = .(participant)
#   , .fun = function(x){
#     i = unique(x$participant)
#     x_use = x[x$parameter == "logitRhoEffectMean",]$value
#     logitrhoeffect_use = median(logitrhoeffect) + median(logitrhoeffectsd)*median(x_use)+ median(logitrhointeractioneffect)*conventionfactor[i]
#     df = data.frame(logitrhoeffect_use, conventionfactor[i])
#     names(df) = c("logitrhoeffect", "conventionfactor")
#     return(df)
#   }
# )
# 
# psseffect_v_rhoeffect = merge(rhoeffect_ids, psseffect_ids)
# 
# ggplot(data = psseffect_v_rhoeffect, aes(x = logitrhoeffect, y =psseffect, colour = factor(conventionfactor), shape = factor(conventionfactor)))+
#   scale_y_continuous(name = "PSS Effect Mean")+
#   scale_x_continuous(name = "Logit \u03C1 Effect Mean")+
#   geom_hline(yintercept = 0, linetype = 2, size = 1)+
#   geom_vline(xintercept = 0, linetype = 2, size = 1)+
#   geom_point(size = 4)+
#   scale_shape_manual(name = "Know\nConvention", labels = c("True", "False"), values = c(16,17) )+
#   scale_colour_manual(name = "Know\nConvention", labels =c("True", "False") , values = c("red", "blue") )+
#   theme_gray(base_size = 30)+
#   theme(panel.grid.major = element_line(size = 1.5)
#         ,panel.grid.minor = element_line(size = 1))
#  
# ### Violin
# # must take negative because I flipped pss effect sign
# pos_corr_f = pos_corr
# pos_corr_f[pos_corr_f$parameter == "value.2.7",]$value = -pos_corr[pos_corr$parameter == "value.2.7",]$value    
# # plot
# ggplot(
#   data = pos_corr_f[pos_corr_f$parameter == "value.2.7",]  
#   , aes(x = parameter, y = value)
# )+
#   geom_violin()+
#   labs(x = "Logit \u03C1 vs. PSS Effect Means", y = "Correlation Coefficient (r)")+
#   stat_summary(fun.data = get_95_HDI, size = 0.7)+
#   stat_summary(fun.data = get_50_HDI, size = 2.5)+  
#   geom_hline(yintercept = 0, linetype = 2, size = 1)+
#   theme_gray(base_size = 30)+
#   theme(panel.grid.major = element_line(size = 1.5)
#         ,panel.grid.minor = element_line(size = 1)
#         , axis.text.x = element_blank()
#         , axis.ticks.x = element_blank()) 
# 
# # get 95% HDI
# print("NOTE: make sure to flip sign (see above code)")
# get_95_HDI(pos_corr_f[pos_corr_f$parameter == "value.2.7",]$value)
# #---------------------------- Rho vs. PSS Effects -----------------------------------------#


# #---------------------------- Kappa vs. PSS Effects ---------------------------------------#
# logkappaeffect = extract_samples("logKappaEffectMean")
# 
# logkappaeffectsd = extract_samples("zlogKappaEffectSD", TRUE)
# 
# logkappainteractioneffect = extract_samples("logKappaConventionInteractionEffectMean")
# 
# kappaeffect_ids = ddply(
#   .data = betas
#   , .variables = .(participant)
#   , .fun = function(x){
#     i = unique(x$participant)
#     x_use = x[x$parameter == "logKappaEffectMean",]$value
#     logkappaeffect_use = median(logkappaeffect) + median(logkappaeffectsd)*median(x_use) + median(logkappainteractioneffect)*conventionfactor[i]
#     df = data.frame(logkappaeffect_use, conventionfactor[i])
#     names(df) = c("logkappaeffect", "conventionfactor")
#     return(df)
#   }
# )
# 
# psseffect_v_kappaeffect = merge(kappaeffect_ids, psseffect_ids)
# 
# ggplot(data = psseffect_v_kappaeffect, aes(x = logkappaeffect, y =psseffect, colour = factor(conventionfactor), shape = factor(conventionfactor)))+
#   scale_y_continuous(name = "PSS Effect Mean")+
#   scale_x_continuous(name = "Log \u03BA Effect Mean")+
#   geom_hline(yintercept = 0, linetype = 2, size = 1)+
#   geom_vline(xintercept = 0, linetype = 2, size = 1)+
#   geom_point(size = 4)+
#   scale_shape_manual(name = "Know\nConvention", labels = c("True", "False"), values = c(16,17) )+
#   scale_colour_manual(name = "Know\nConvention", labels =c("True", "False") , values = c("red", "blue") )+
#   theme_gray(base_size = 30)+
#   theme(panel.grid.major = element_line(size = 1.5)
#         ,panel.grid.minor = element_line(size = 1))
# 
# ### Violin
# # must take negative because I flipped pss effect sign
# pos_corr_f = pos_corr
# pos_corr_f[pos_corr_f$parameter == "value.2.8",]$value = -pos_corr[pos_corr$parameter == "value.2.8",]$value    
# # plot
# ggplot(
#   data = pos_corr_f[pos_corr_f$parameter == "value.2.8",]  
#   , aes(x = parameter, y = value)
# )+
#   geom_violin()+
#   labs(x = "Log \u03BA vs. PSS Effect Means", y = "Correlation Coefficient (r)")+
#   stat_summary(fun.data = get_95_HDI, size = 0.7)+
#   stat_summary(fun.data = get_50_HDI, size = 2.5)+  
#   geom_hline(yintercept = 0, linetype = 2, size = 1)+
#   theme_gray(base_size = 30)+
#   theme(panel.grid.major = element_line(size = 1.5)
#         ,panel.grid.minor = element_line(size = 1)
#         , axis.text.x = element_blank()
#         , axis.ticks.x = element_blank()) 
# 
# # get 95% HDI
# print("NOTE: make sure to flip sign (see above code)")
# get_95_HDI(pos_corr_f[pos_corr_f$parameter == "value.2.8",]$value)
# #---------------------------- Kappa vs. PSS Effects ---------------------------------------#



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
  exp(ex_toj_color_post$population_logjnd_intercept_mean + ex_toj_color_post$population_logjnd_convention_effect_mean/2) * 250
)
# K
get_95_HDI(
  exp( ex_toj_color_post$population_logjnd_intercept_mean - ex_toj_color_post$population_logjnd_convention_effect_mean/2)  * 250
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
  ( exp(ex_toj_color_post$population_logjnd_intercept_mean + ex_toj_color_post$population_logjnd_convention_interaction_effect_mean/2) 
    - exp(ex_toj_color_post$population_logjnd_intercept_mean - ex_toj_color_post$population_logjnd_convention_interaction_effect_mean/2) ) * 250
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
  exp(ex_toj_color_post$population_logjnd_intercept_mean + ex_toj_color_post$population_logjnd_convention_effect_mean/2 + (ex_toj_color_post$population_logjnd_effect_mean + ex_toj_color_post$population_logjnd_convention_interaction_effect_mean)/2) * 250
)
# DK + U
get_95_HDI(
  exp(ex_toj_color_post$population_logjnd_intercept_mean + ex_toj_color_post$population_logjnd_convention_effect_mean/2 - (ex_toj_color_post$population_logjnd_effect_mean + ex_toj_color_post$population_logjnd_convention_interaction_effect_mean)/2  ) * 250
)
# K + A
get_95_HDI(
  exp(ex_toj_color_post$population_logjnd_intercept_mean -  ex_toj_color_post$population_logjnd_convention_effect_mean/2 + (ex_toj_color_post$population_logjnd_effect_mean - ex_toj_color_post$population_logjnd_convention_interaction_effect_mean)/2) * 250 
)
# K + U
get_95_HDI(
  exp(ex_toj_color_post$population_logjnd_intercept_mean -  ex_toj_color_post$population_logjnd_convention_effect_mean/2 - (ex_toj_color_post$population_logjnd_effect_mean - ex_toj_color_post$population_logjnd_convention_interaction_effect_mean)/2 )  * 250
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
  plogis(ex_toj_color_post$logitRhoMean)
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
  ( plogis(ex_toj_color_post$logitRhoMean + ex_toj_color_post$logitRhoEffectMean/2 )
    - plogis(ex_toj_color_post$logitRhoMean - ex_toj_color_post$logitRhoEffectMean/2 ) )
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
  ( plogis(ex_toj_color_post$logitRhoMean + ex_toj_color_post$logitRhoConventionEffectMean/2 )
    - plogis(ex_toj_color_post$logitRhoMean - ex_toj_color_post$logitRhoConventionEffectMean/2 ) )
  , "Probability of Encoding\nConvention Effect Mean"
  , y_lab = "\u03C1 (Don't Know - Know)"
)
# DK
get_95_HDI(
  plogis(ex_toj_color_post$logitRhoMean + ex_toj_color_post$logitRhoConventionEffectMean/2)
  )
# K
get_95_HDI(
  plogis(ex_toj_color_post$logitRhoMean - ex_toj_color_post$logitRhoConventionEffectMean/2 )
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
  ( plogis(ex_toj_color_post$logitRhoMean + ex_toj_color_post$logitRhoConventionInteractionEffectMean/2 )
    - plogis(ex_toj_color_post$logitRhoMean - ex_toj_color_post$logitRhoConventionInteractionEffectMean/2 ) )
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
  ( plogis(ex_toj_color_post$logitRhoMean + ex_toj_color_post$logitRhoConventionEffectMean/2
           + (ex_toj_color_post$logitRhoEffectMean + ex_toj_color_post$logitRhoConventionInteractionEffectMean)/2 )
    - plogis(ex_toj_color_post$logitRhoMean + ex_toj_color_post$logitRhoConventionEffectMean/2
             - (ex_toj_color_post$logitRhoEffectMean +  ex_toj_color_post$logitRhoConventionInteractionEffectMean)/2 ) )
  , ( plogis(ex_toj_color_post$logitRhoMean - ex_toj_color_post$logitRhoConventionEffectMean/2
             + (ex_toj_color_post$logitRhoEffectMean - ex_toj_color_post$logitRhoConventionInteractionEffectMean)/2 )
      - plogis(ex_toj_color_post$logitRhoMean - ex_toj_color_post$logitRhoConventionEffectMean/2
               - (ex_toj_color_post$logitRhoEffectMean -  ex_toj_color_post$logitRhoConventionInteractionEffectMean)/2 ) )
  )
  , c("Probability of Encoding\nAttention Effect\nGiven Don't Know"  , "Probability of Encoding\nAttention Effect\nGiven Know")
  , y_lab = "\u03C1 (Attended - Unattended)"
)
# Know condition estimate 
get_95_HDI(
  ( plogis(ex_toj_color_post$logitRhoMean - ex_toj_color_post$logitRhoConventionEffectMean/2
           + (ex_toj_color_post$logitRhoEffectMean - ex_toj_color_post$logitRhoConventionInteractionEffectMean)/2 )
    - plogis(ex_toj_color_post$logitRhoMean - ex_toj_color_post$logitRhoConventionEffectMean/2
             - (ex_toj_color_post$logitRhoEffectMean -  ex_toj_color_post$logitRhoConventionInteractionEffectMean)/2 ) )
)

# A + DK
get_95_HDI(
  plogis(ex_toj_color_post$logitRhoMean + ex_toj_color_post$logitRhoConventionEffectMean/2
                     + (ex_toj_color_post$logitRhoEffectMean + ex_toj_color_post$logitRhoConventionInteractionEffectMean)/2 )
)
# U + DK
get_95_HDI(
  plogis(ex_toj_color_post$logitRhoMean + ex_toj_color_post$logitRhoConventionEffectMean/2
         - (ex_toj_color_post$logitRhoEffectMean +  ex_toj_color_post$logitRhoConventionInteractionEffectMean)/2 ) 
)
# A + K
get_95_HDI(
  plogis(ex_toj_color_post$logitRhoMean - ex_toj_color_post$logitRhoConventionEffectMean/2
         + (ex_toj_color_post$logitRhoEffectMean - ex_toj_color_post$logitRhoConventionInteractionEffectMean)/2 )
)
# U + K
get_95_HDI(
  plogis(ex_toj_color_post$logitRhoMean - ex_toj_color_post$logitRhoConventionEffectMean/2
         - (ex_toj_color_post$logitRhoEffectMean -  ex_toj_color_post$logitRhoConventionInteractionEffectMean)/2 ) 
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
  ( exp(ex_toj_color_post$logKappaMean + ex_toj_color_post$logKappaConventionEffectMean/2 )
    - exp(ex_toj_color_post$logKappaMean - ex_toj_color_post$logKappaConventionEffectMean/2 ) )
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
  ( exp(ex_toj_color_post$logKappaMean + ex_toj_color_post$logKappaConventionInteractionEffectMean/2 )
    - exp(ex_toj_color_post$logKappaMean - ex_toj_color_post$logKappaConventionInteractionEffectMean/2 ) )
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
     + median(ex_toj_color_post$population_pss_convention_interaction_effect_mean)
     )/2 
  )* 250
  , sd = ( exp( 
    median(ex_toj_color_post$population_logjnd_intercept_mean) + median(ex_toj_color_post$population_logjnd_convention_effect_mean)/2 
    + (
      median(ex_toj_color_post$population_logjnd_effect_mean) 
      + median(ex_toj_color_post$population_logjnd_convention_interaction_effect_mean)
      ) /2
  ) ) * 250
)

yGloveKnow = pnorm(
  -250:250
  , mean = ( 
    median(ex_toj_color_post$population_pss_intercept_mean) - median(ex_toj_color_post$population_pss_convention_effect_mean)/2 
    + (
      median(ex_toj_color_post$population_pss_effect_mean) 
      - median(ex_toj_color_post$population_pss_convention_interaction_effect_mean)
    )/2 
  )* 250
  , sd = ( exp( 
    median(ex_toj_color_post$population_logjnd_intercept_mean) - median(ex_toj_color_post$population_logjnd_convention_effect_mean)/2 
    + (
      median(ex_toj_color_post$population_logjnd_effect_mean) 
      - median(ex_toj_color_post$population_logjnd_convention_interaction_effect_mean)
    ) /2
  ) ) * 250
)

yBaseDontKnow = pnorm(
  -250:250
  , mean = ( 
    median(ex_toj_color_post$population_pss_intercept_mean) + median(ex_toj_color_post$population_pss_convention_effect_mean)/2 
    - (
      median(ex_toj_color_post$population_pss_effect_mean) 
      + median(ex_toj_color_post$population_pss_convention_interaction_effect_mean)
    )/2 
  )* 250
  , sd = ( exp( 
    median(ex_toj_color_post$population_logjnd_intercept_mean) + median(ex_toj_color_post$population_logjnd_convention_effect_mean)/2 
    - (
      median(ex_toj_color_post$population_logjnd_effect_mean) 
      + median(ex_toj_color_post$population_logjnd_convention_interaction_effect_mean)
    ) /2
  ) ) * 250
)

yBaseKnow = pnorm(
  -250:250
  , mean = ( 
    median(ex_toj_color_post$population_pss_intercept_mean) - median(ex_toj_color_post$population_pss_convention_effect_mean)/2 
    - (
      median(ex_toj_color_post$population_pss_effect_mean) 
      - median(ex_toj_color_post$population_pss_convention_interaction_effect_mean)
    )/2 
  )* 250
  , sd = ( exp( 
    median(ex_toj_color_post$population_logjnd_intercept_mean) - median(ex_toj_color_post$population_logjnd_convention_effect_mean)/2 
    - (
      median(ex_toj_color_post$population_logjnd_effect_mean) 
      - median(ex_toj_color_post$population_logjnd_convention_interaction_effect_mean)
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

