#### Libraries ####
library(plyr)
library(ggplot2)
library(grid)
library(rstan)
library(ez)

setwd("~/Documents/TOJ/Baseball/baseballtojdata")



###############################################################################
####                         Data Import                                   ####
###############################################################################

check_before_read = function(file){
  temp = tools::md5sum(file)
  if(temp %in% checksums){
    return(NULL)
  }else{
    checksums <<- c(checksums,temp)
    to_return = read.table(
      file
      , header = T
      , sep = '\t'
    )
    return(to_return)
  }
}

checksums = NULL

a = ldply(
  .data = list.files(
    pattern = ".txt"
    , full.names = T
    , path = './Data_alpha'
  )
  , .fun = check_before_read
  , .progress = 'text'
)
a$id = paste('alpha',a$participant_id,a$created)

b = ldply(
  .data = list.files(
    pattern = ".txt"
    , full.names = T
    , path = './Data_delta'
  )
  , .fun = check_before_read
  , .progress = 'text'
)
b$id = paste('delta',b$participant_id,b$created)

a = rbind(a,b)
length(unique(a$id))

length(checksums)

#count trials per id
temp = data.frame(table(id=a$id))

#discard those with too few trials
keep = temp$id[temp$Freq==480]
a = a[a$id %in% keep,]

length(unique(a$id))



#toss baseball-naive Ss
a = a[substr(a$id,1,8)!='alpha 12',]
a = a[substr(a$id,1,7)!='delta 7',]
a = a[substr(a$id,1,8)!='delta 20',]
length(unique(a$id))  # count files


#### Type ####
#As factor 
a$gender = as.factor(a$gender)
a$age = as.factor(a$age)
a$handedness = as.factor(a$handedness)
a$created = as.factor(a$created)
a$block_num = as.factor(a$block_num)
a$trial_num = as.factor(a$trial_num)
#a$soa = as.factor(a$soa)
a$baserun_offset = as.factor(a$baserun_offset)
a$first_arrival = as.factor(a$first_arrival)
a$probed_trial = as.factor(a$probed_trial)
a$glove_probe_dist = as.factor(a$glove_probe_dist)
a$base_probe_dist = as.factor(a$base_probe_dist)
a$probe_location = as.factor(a$probe_location)
a$toj_response = as.factor(a$toj_response)

#As numeric
#a$probe_color = as.numeric(a$probe_color)
#a$color_response = as.numeric(a$color_response)
a$color_diff = as.numeric(a$color_diff)
a$response_time = as.numeric(a$response_time)



#### Summary ####
summary(a)

# double check 20/80
table(a$probe_location, a$base_probe_dist)  # by base
table(a$probe_location, a$glove_probe_dist)  # by glove. Should be symetric 
toss = NULL
for(i in unique(a$id)){
  temp = with(
    a[a$id==i,]
    ,data.frame(table(
      probe_location, glove_probe_dist
    ))
  )
  if(
    (!all(temp$Freq==c(64,16,0,0,16,64,0,0))) &
    (!all(temp$Freq==c(0,0,64,16,0,0,16,64)))
  ){
    toss = c(toss,i)
    print(i)
    print(temp)
  }
}

a = a[!(a$id %in% toss),]

length(unique(a$id))

### add at 'tie goes to runner' summary
a$know_tie_goes_runner = FALSE
a[a$id == "alpha 4 2015-02-04 13:52:56"
  | a$id == "alpha 8 2015-02-09 13:38:03"
  | a$id == "alpha 15 2015-02-10 17:07:37"
  | a$id == "delta 3 2014-12-03 09:03:59"
  | a$id == "delta 4 2014-12-03 10:11:45"
  | a$id == "delta 5 2014-12-04 10:10:39"
  | a$id == "delta 6 2014-12-04 11:07:22"
  | a$id == "delta 11 2014-12-05 09:58:42"
  | a$id == "delta 14 2014-12-05 13:04:02"
  | a$id == "delta 26 2014-12-10 15:18:43"
  | a$id == "delta 28 2015-02-02 13:41:26"
  | a$id == "delta 29 2015-02-02 14:42:42",]$know_tie_goes_runner = TRUE

a$use_tie_goes_runner = FALSE
a[a$id == "alpha 4 2015-02-04 13:52:56"
  | a$id == "alpha 15 2015-02-10 17:07:37"
  | a$id == "delta 3 2014-12-03 09:03:59"
  | a$id == "delta 4 2014-12-03 10:11:45"
  | a$id == "delta 5 2014-12-04 10:10:39"
  | a$id == "delta 6 2014-12-04 11:07:22"
  | a$id == "delta 14 2014-12-05 13:04:02"
  | a$id == "delta 28 2015-02-02 13:41:26"
  | a$id == "delta 29 2015-02-02 14:42:42",]$use_tie_goes_runner = TRUE

# get rid of delta 8 because flat psychometric function
# REASON TO BELIEVE PARTICIPANT WAS NOT DOING THE TASK PROPERLY
a= a[a$id != "delta 8 2014-12-04 14:23:37",]

length(unique(a$id))



###############################################################################
####                              TOJ                                      ####
###############################################################################

toj_trials = a
toj_trials = toj_trials[!is.na(toj_trials$toj_response), ]
toj_trials$safe = FALSE
toj_trials$safe[toj_trials$toj_response == "safe"] = TRUE

### SOAs 
toj_trials$soa2 = toj_trials$soa
# correct soas 
toj_trials[toj_trials$soa2 == "15",]$soa2 = 17
toj_trials[toj_trials$soa2 == "45",]$soa2 = 50
toj_trials[toj_trials$soa2 == "90",]$soa2 = 100
toj_trials[toj_trials$soa2 == "135",]$soa2 = 150
toj_trials[toj_trials$soa2 == "240",]$soa2 = 250
# Negative SOAs means Ball first 
toj_trials$soa2[toj_trials$first_arrival == "ball"] = -toj_trials$soa2[toj_trials$first_arrival == "ball"]

toj_means_by_id_by_condition = ddply(
  .data = toj_trials
  , .variables = .(id,base_probe_dist, soa2)
  , .fun = function(x){
    to_return = data.frame(
      value = mean(x$safe)
    )
    return(to_return)
  }
)
toj_means_by_id_by_condition$soa2 = as.numeric(as.character(toj_means_by_id_by_condition$soa2))

ggplot(
  data = toj_means_by_id_by_condition
  , mapping = aes(
    x = soa2
    , y =  value
    , shape = base_probe_dist
    , linetype = base_probe_dist
    , group = base_probe_dist
  )
)+
  facet_wrap(
    ~ id
  )+
#   geom_line()+
#   geom_point()+
  geom_smooth(
    method = "glm"
    , method.args = list(family = "binomial")
    , formula = y ~ splines::ns(x,3)
     , se = FALSE
#     , level = 0.95  # 95% confidence interval
    )

# function
medsd = function(x){
  sqrt(   median((x-median(x))^2)    )
}

#Assesing the bias of particpants (accross conditions).
toj_means_by_id = ddply(
  .data = toj_trials
  , .variables = .(id)
  , .fun = function(x){
    fit = glm(
      formula = safe~soa2
      , data = x
      , family = binomial
    )
    to_return = data.frame(
      id = x$id[1]
      , pss = -coef(fit)[1]/coef(fit)[2]
      , slope = coef(fit)[2]
    )
    return(to_return)
  }
)

hist(toj_means_by_id$pss,br=100)

# #### PSS Cutoff ####
# pss_cuttoff = medsd(toj_means_by_id$pss)*5
# abline(v=-pss_cuttoff)
# abline(v=pss_cuttoff)
# 
# toj_means_by_id$toss = abs(toj_means_by_id$pss)>pss_cuttoff
# unique(length(toj_means_by_id$id[!toj_means_by_id$toss]))
# 
# hist(toj_means_by_id$slope,br=100)
# slope_cutoff = with(
#   toj_means_by_id[!toj_means_by_id$toss,]
# #### SLOPE Cutoff ####
#   , median(slope)-medsd(slope)*5
# )
# abline(v=slope_cutoff)
# toj_means_by_id$toss[!toj_means_by_id$toss] = toj_means_by_id$slope[!toj_means_by_id$toss]<slope_cutoff
# unique(length(toj_means_by_id$id[!toj_means_by_id$toss]))
# 
# ggplot(
#   data = toj_means_by_id
#   , mapping = aes(
#     x = pss
#     , y = slope
#     , label = id
#     , colour = toss
#   )
# )+
#   geom_text()
# 
# toj_trials$toss = toj_trials$id %in% toj_means_by_id$id[toj_means_by_id$toss]
# print("TOJ Tossed Count:")
# length(unique(toj_trials[toj_trials$toss == TRUE,]$id))



###############################################################################
####                              Color                                    ####
###############################################################################
short_angle = function(x, y)
{
  return(((x - y + 180) %% 360) - 180)
}

degree_to_rad = function(x)
{
  return(x*pi / 180)
}


rad_to_degrees = function(x)
{
  return(x*180 / pi)
}

color_trials = a
color_trials$probe_location[color_trials$probe_location == "[1088, 896]"] = "base"
color_trials$probe_location[color_trials$probe_location == "[1328, 581]"] = "glove"
color_trials = color_trials[!is.na(color_trials$color_diff), ]
hist(color_trials$color_diff,br=100)
color_trials[color_trials$color_diff > 180,]$color_diff = - (360 - color_trials[color_trials$color_diff > 180,]$color_diff)
color_trials[color_trials$color_diff < (-180),]$color_diff = color_trials[color_trials$color_diff < (-180),]$color_diff + 360
hist(color_trials$color_diff,br=100)

color_trials$abs_color_diff = abs(color_trials$color_diff)
hist(color_trials$abs_color_diff, br = 100)

color_trials$attended = FALSE
color_trials$attended[ (color_trials$base_probe_dist == 0.8 & color_trials$probe_location == "base") | (color_trials$base_probe_dist == 0.2 & color_trials$probe_location == "glove")] = TRUE

# color_trials$toj_tossed = color_trials$id %in% unique(toj_trials$id[toj_trials$toss])
# 
# ggplot(
#   data = color_trials
#   , mapping = aes(
#     x = color_diff
#     , fill = toj_tossed
#   )
# )+
#   facet_wrap(
#     ~ id
#   )+
#   geom_histogram()
# 
# ggplot(
# 	data = color_trials
# 	, mapping = aes(
# 		x = abs_color_diff
# 		, fill = toj_tossed
# 	)
# )+
# 	facet_wrap(
# 		~ id
# 	)+
# 	geom_histogram()
# 
# # toss outlier participants 
# color_trials = color_trials[!color_trials$toj_tossed,]
# toj_trials = toj_trials[!toj_trials$toss,]
# 
# # make sure we exclude the same participants for each task data set 
# if (length(unique(color_trials$id)) == length(unique(toj_trials$id)) ){
#   print("Final Count")
#   length(unique(color_trials$id))  # what is the participant count?
# } else {
#   print("Error: uneven final color wheel and toj trial counts")
# }

 

###############################################################################
####                        Statistical Analyses                           ####
###############################################################################

color_means = aggregate(abs_color_diff ~ attended + know_tie_goes_runner + id , data = color_trials, FUN = mean)

#-----------------------------------------------------------------------------#
#                            Absolute Color Diff                              #
#-----------------------------------------------------------------------------#


#----------------------------------- Degrees ---------------------------------#
# descriptive statistics 
color_tot_means = aggregate(abs_color_diff ~ attended, data=color_means, FUN=mean)
color_tot_SD = aggregate(abs_color_diff ~ attended, data=color_means, FUN=sd)

# # ASSUMPTIONS
# # (1) Normality
# m = lm(abs_color_diff ~ attended, data = color_means)
# plot(m) # look at second plot: residuals vs. theoreticcal quantiles 
# # (2) Equal variance
# color_vars =  aggregate(abs_color_diff ~ attended, data = color_means, FUN = var)
# # these are pretty damn equal 
# 
# # t-test
# t.test(
#   formula = abs_color_diff ~ attended
#   , dat = color_means
#   , paired = T
#   , var.equal = T 
# )
# 
# ### with convention knowledge 
# color_know = aggregate(abs_color_diff ~ attended + know_tie_goes_runner + id, data=color_trials, FUN = mean)
# color_know$know_tie_goes_runner = as.factor(as.character(color_know$know_tie_goes_runner))

# descriptive 
color_know_means = aggregate(abs_color_diff ~ attended + know_tie_goes_runner, data = color_know, FUN = mean)
color_know_SD =  aggregate(abs_color_diff ~ attended + know_tie_goes_runner, data = color_know, FUN = sd)

# # ASSUMPTIONS
# # (1) Normality 
# # NA NA 
# # res_list = NULL
# # res_list = c(m$id$residuals, m$`id:know_tie_goes_runner`$residuals, m$`id:attended`$residuals, m$`id:attended:know_tie_goes_runner`$residuals)
# # qqnorm(res_list)
# # qqline(res_list)
# # # (2) Sphericity - cannot check because only two levels by factor 
# 
# # model
# m = aov(
#   formula = abs_color_diff ~ attended*know_tie_goes_runner + Error(id/(attended)) # not devided by between subject..
#   , dat = color_know
# )
# 
# # anova
# summary(m)
#----------------------------------- Degrees ---------------------------------#


#----------------------------------- Log Degrees -----------------------------#
color_means$log_abs_color_diff = log(color_means$abs_color_diff)

# descriptive statistics 
color_log_tot_means = aggregate(log_abs_color_diff ~ attended, data=color_means, FUN=mean)
color_log_tot_SD = aggregate(log_abs_color_diff ~ attended, data=color_means, FUN=sd)

# ASSUMPTIONS
# (1) Normality
m = lm(log_abs_color_diff ~ attended, data = color_means)
plot(m) # look at second plot: residuals vs. theoreticcal quantiles 
# STILL OUTLIERS AT FAR RIGHT 
# TAKING LOG A SECOND TIME FIXES THIS BUT PERHAPS BAD IDEA
# (2) Equal variance
color_vars =  aggregate(log_abs_color_diff ~ attended, data = color_means, FUN = var)
# these are pretty damn equal 

# t-test
t.test(
  formula = log_abs_color_diff ~ attended
  , dat = color_means
  , paired = T
  , var.equal = T 
)

### with convention knowledge 
color_log_know = aggregate(log_abs_color_diff ~ attended + know_tie_goes_runner + id, data=color_means, FUN = mean)
color_log_know$know_tie_goes_runner = as.factor(as.character(color_log_know$know_tie_goes_runner))

# descriptive 
color_log_know_means = aggregate(log_abs_color_diff ~ attended + know_tie_goes_runner, data = color_log_know, FUN = mean)
color_log_know_SD =  aggregate(log_abs_color_diff ~ attended + know_tie_goes_runner, data = color_log_know, FUN = sd)

# ASSUMPTIONS
# (1) Normality 
# NA NA 
# res_list = NULL
# res_list = c(m$id$residuals, m$`id:know_tie_goes_runner`$residuals, m$`id:attended`$residuals, m$`id:attended:know_tie_goes_runner`$residuals)
# qqnorm(res_list)
# qqline(res_list)
# # (2) Sphericity - cannot check because only two levels by factor 

# model
m = aov(
  formula = log_abs_color_diff ~ attended*know_tie_goes_runner + Error(id/(attended)) # not devided by between subject..
  , dat = color_log_know
)

# anova
summary(m)
#----------------------------------- Log Degrees -----------------------------#



# #### Median Split Subgroup ####
# color_means_diff = aggregate(abs_color_diff ~ id, data = color_means, FUN = diff)
# color_means_diff$abs_color_diff = -color_means_diff$abs_color_diff  # make unattended - attended 
# med = median(color_means_diff$abs_color_diff)
# id_keep = color_means_diff[color_means_diff$abs_color_diff > med,]$id



#-----------------------------------------------------------------------------#
#                             Mixture Model                                   #
#-----------------------------------------------------------------------------#

source("/Users/ghislaindentremont/Documents/TOJ/EndogenousVisualPriorEntry-BayesianHierarchicalModel/Baseball/conventional_analysis/fit_uvm.R")

color_trials$color_diff_radians = color_trials$color_diff*pi/180

fitted = ddply(
    .data = color_trials
    , .variables = .(id, attended, know_tie_goes_runner)
    , .fun = function(piece_of_df){
      fit = fit_uvm(piece_of_df$color_diff_radians, do_mu = TRUE)
      to_return = data.frame(
        kappa_prime = fit$kappa_prime
        , kappa = exp(fit$kappa_prime)
        , rho = fit$rho
        , logit_rho = qlogis(fit$rho - 0.000001) # This is PROBLEMATIC
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
    to_return = data.frame(
      kappa_prime = fit$kappa_prime
      , kappa = exp(fit$kappa_prime)
      , rho = fit$rho
      , logit_rho = qlogis(fit$rho - 0.000001) # This is PROBLEMATIC 
    )
    return(to_return)
  }
  , .progress = 'time'
)


#----------------------- Fidelity (kappa prime) ------------------------------#
kappa_prime_tot_means = aggregate(kappa_prime ~ attended, data = fitted, FUN = mean)
kappa_prime_tot_SD = aggregate(kappa_prime ~ attended, data = fitted, FUN = sd)

# ASSUMPTIONS
# (1) Normality
m = lm(kappa_prime ~ attended, data = fitted)
plot(m) # look at second plot: residuals vs. theoreticcal quantiles 
# (2) Equal variance
kappa_vars =  aggregate(kappa_prime ~ attended, data = fitted, FUN = var)
diffreal = kappa_vars$kappa_prime[1] - kappa_vars$kappa_prime[2]
# these seem very close 
# but how close?
# randomization test:
kap = fitted$kappa_prime
varz = NULL
for (i in 1:10000) {
  samp = sample(kap)
  len = length(kap)
  one = samp[1:(len/2)]
  two = samp[(len/2 + 1):len]
  temp = var(one) - var(two)
  varz = c(varz, temp)
}
hist(varz, br = 100)
abline(v = diffreal, col = 'red') # definately within range of what is possible with population equal variances

# test 
t.test(
  formula = kappa_prime ~ attended
  , data = fitted
  , paired = TRUE
  , var.equal = T  # equal variances doesn't make difference
)

### with convention
color_kappa_prime = aggregate(kappa_prime ~ attended + know_tie_goes_runner + id, data=fitted, FUN = mean)
color_kappa_prime$know_tie_goes_runner = as.factor(as.character(color_kappa_prime$know_tie_goes_runner))

# descriptve 
color_kappa_prime_means = aggregate(kappa_prime ~ attended + know_tie_goes_runner, data = color_kappa, FUN = mean)
color_kappa_prime_SD =  aggregate(kappa_prime ~ attended + know_tie_goes_runner, data = color_kappa, FUN = sd)

# ASSUMPTIONS
# (1) Normality 
# NOT SURE HOW TO GET RESIDUALS
# res_list = NULL
# res_list = c(m$id$residuals, m$`id:know_tie_goes_runner`$residuals, m$`id:attended`$residuals, m$`id:attended:know_tie_goes_runner`$residuals)
# qqnorm(res_list)
# qqline(res_list)
# (2) Sphericity - cannot check because only two levels by factor 

# model
m = aov(
  formula = kappa_prime ~ attended*know_tie_goes_runner + Error(id/(attended))
  , dat = color_kappa
)

# anova
summary(m)
#----------------------- Fidelity (kappa prime) ------------------------------#



#------------------------- Fidelity (kappa) ----------------------------------#
kappa_tot_means = aggregate(kappa ~ attended, data = fitted_attn, FUN = mean)
kappa_tot_SD = aggregate(kappa ~ attended, data = fitted, FUN = sd)

# # ASSUMPTIONS
# # (1) Normality
# m = lm(kappa ~ attended, data = fitted)
# plot(m) # look at second plot: residuals vs. theoreticcal quantiles 
# # (2) Equal variance
# kappa_vars =  aggregate(kappa ~ attended, data = fitted, FUN = var)
# diffreal = kappa_vars$kappa[1] - kappa_vars$kappa[2]
# # these seem very close 
# # but how close?
# # randomization test:
# kap = fitted$kappa
# varz = NULL
# for (i in 1:10000) {
#   samp = sample(kap)
#   len = length(kap)
#   one = samp[1:(len/2)]
#   two = samp[(len/2 + 1):len]
#   temp = var(one) - var(two)
#   varz = c(varz, temp)
# }
# hist(varz, br = 100)
# abline(v = diffreal, col = 'red') # definately within range of what is possible with population equal variances
# 
# # test 
# t.test(
#   formula = kappa ~ attended
#   , data = fitted
#   , paired = TRUE
#   , var.equal = T  # equal variances doesn't make difference
# )
# 
# ### with convention knowledge 
# color_kappa = aggregate(kappa ~ attended + know_tie_goes_runner + id, data=fitted, FUN = mean)
# color_kappa$know_tie_goes_runner = as.factor(as.character(color_kappa$know_tie_goes_runner))

# descriptve 
color_kappa_means = aggregate(kappa ~ attended + know_tie_goes_runner, data = color_kappa, FUN = mean)
color_kappa_SD =  aggregate(kappa ~ attended + know_tie_goes_runner, data = color_kappa, FUN = sd)

# # ASSUMPTIONS
# # (1) Normality 
# # NOT SURE HOW TO GET RESIDUALS
# # res_list = NULL
# # res_list = c(m$id$residuals, m$`id:know_tie_goes_runner`$residuals, m$`id:attended`$residuals, m$`id:attended:know_tie_goes_runner`$residuals)
# # qqnorm(res_list)
# # qqline(res_list)
# # (2) Sphericity - cannot check because only two levels by factor 
# 
# # model
# m = aov(
#   formula = kappa ~ attended*know_tie_goes_runner + Error(id/(attended))
#   , dat = color_kappa
# )
# 
# # anova
# summary(m)
#------------------------- Fidelity (kappa) ----------------------------------#


#------------------------ Probability (rho) ----------------------------------# 
rho_tot_means = aggregate(rho ~ attended, data = fitted, FUN = mean)
rho_tot_SD = aggregate(rho ~ attended, data = fitted, FUN = sd)

# ASSUMPTIONS
# (1) Normality
m = lm(rho ~ attended, data = fitted)
plot(m) # look at second plot: residuals vs. theoreticcal quantiles 
# NOT NORMALLY DISTRIBUTED! 
# (2) Equal variance
rho_vars =  aggregate(rho ~ attended, data = fitted, FUN = var)
# close enough, perhaps...

# how many participants are perfect across the board?
perfection_count = sum(fitted_all$rho == 1)
perfection_rate = perfection_count/nrow(fitted_all)

ggplot(
  data = fitted_all
  , mapping = aes(rho)  #, fill = attended)
)+ 
  geom_histogram(bins = 30)+
  labs(x = "Probability of Memory", y = "Number of Occurences")+
  theme_gray(base_size = 18)

ggplot(
  data = fitted
  , mapping = aes(rho, fill = attended)
)+
  geom_histogram(bins = 20)+
  # geom_histogram(alpha = 0.5, position = "identity")+
  scale_fill_discrete(name = "", labels = c("Unattended", "Attended"))+
  labs(x = "Probability of Memory", y = "Number of Occurences")+
  theme_gray(base_size = 18)+
#   geom_vline(xintercept = 0.952, colour = "red")+
#   geom_vline(xintercept= 0.967, colour = "blue") 
  facet_grid(attended~.) +
  theme(strip.background = element_blank(), strip.text = element_blank()) 

# test 
t.test(
  formula = rho ~ attended
  , data = fitted
  , paired = TRUE
  , var.equal = T
)

### with convention knowledge 
color_rho = aggregate(rho ~ attended + know_tie_goes_runner + id, data=fitted, FUN = mean)
color_rho$know_tie_goes_runner = as.factor(as.character(color_rho$know_tie_goes_runner))

# descriptive 
color_rho_means = aggregate(rho ~ attended + know_tie_goes_runner, data = color_rho, FUN = mean)
color_rho_SD =  aggregate(rho ~ attended + know_tie_goes_runner, data = color_rho, FUN = sd)

# ASSUMPTIONS
# (1) Normality 
# NOT SURE HOW TO GET RESIDUALS
# res_list = NULL
# res_list = c(m$id$residuals, m$`id:know_tie_goes_runner`$residuals, m$`id:attended`$residuals, m$`id:attended:know_tie_goes_runner`$residuals)
# qqnorm(res_list)
# qqline(res_list)
# (2) Sphericity - cannot check because only two levels by factor 

# model
m = aov(
  formula = rho ~ attended*know_tie_goes_runner + Error(id/(attended))
  , dat = color_rho
)

# anova
summary(m)

# Count how many P = 1 by condition
fitted$new_rho = fitted$rho == 1
fitted$new_rho
table(fitted$attended, fitted$new_rho)  # pretty clear effect here

# look at effects: attended - unattended 
# positive means attention increases rho
rho_effects = aggregate(rho ~ id, data = fitted, FUN = diff) 
positive_effects_ratio = sum(rho_effects$rho > 0)/nrow(rho_effects) 
negative_effects_ratio= sum(rho_effects$rho < 0)/nrow(rho_effects) 
null_effects_ratio = sum(rho_effects$rho == 0)/nrow(rho_effects)

ggplot(
  data = rho_effects
  , mapping = aes(rho)  #, fill = attended)
)+ 
  geom_histogram(bins = 15)+
  labs(x = "Probability of Memory Differences: Attended - Unattended", y = "Number of Occurences")+
  theme_gray(base_size = 18)+
  scale_x_continuous(breaks = seq(-0.1, 0.2, 0.05))
#------------------------ Probability (rho) ----------------------------------# 


# NOTE: What to do with ceilings as per EM method? 
# #------------------------ Prob. (logit rho) ----------------------------------#
# logit_rho_tot_means = aggregate(logit_rho ~ attended, data = fitted, FUN = mean)
# logit_rho_tot_SD = aggregate(logit_rho ~ attended, data = fitted, FUN = sd)
# 
# # ASSUMPTIONS
# # (1) Normality
# m = lm(logit_rho ~ attended, data = fitted)
# plot(m) # look at second plot: residuals vs. theoreticcal quantiles 
# 
# # test 
# t.test(
#   formula = logit_rho ~ attended
#   , data = fitted
#   , paired = TRUE
#   , var.equal = T
# )
# 
# ### with convention knowledge 
# color_logit_rho = aggregate(logit_rho ~ attended + know_tie_goes_runner + id, data=fitted, FUN = mean)
# color_logit_rho$know_tie_goes_runner = as.factor(as.character(color_logit_rho$know_tie_goes_runner))
# 
# # descriptive 
# color_logit_rho_means = aggregate(logit_rho ~ attended + know_tie_goes_runner, data = color_logit_rho, FUN = mean)
# color_logit_rho_SD =  aggregate(logit_rho ~ attended + know_tie_goes_runner, data = color_logit_rho, FUN = sd)
# 
# # ASSUMPTIONS
# # (1) Normality 
# # NA NA NA
# # res_list = NULL
# # res_list = c(m$id$residuals, m$`id:know_tie_goes_runner`$residuals, m$`id:attended`$residuals, m$`id:attended:know_tie_goes_runner`$residuals)
# # qqnorm(res_list)
# # qqline(res_list)
# # (2) Sphericity - cannot check because only two levels by factor 
# 
# # model
# m = aov(
#   formula = logit_rho ~ attended*know_tie_goes_runner + Error(id/(attended*know_tie_goes_runner))
#   , dat = color_logit_rho
# )
# 
# # anova
# summary(m)
# #------------------------ Prob. (logit rho) ----------------------------------# 



#-----------------------------------------------------------------------------#
#                                 TOJ                                         #
#-----------------------------------------------------------------------------#

# define function to graph mean TOJs
get_psycho_function = function(id_list) {
  toj_means = ddply(
    .data = toj_trials[toj_trials$id %in% id_list,]
    , .variables = .(base_probe_dist, soa2)  # ,room)
    , .fun = function(x){
      to_return = data.frame(
        value = mean(x$safe)
      )
      return(to_return)
    }
  )
  toj_means$soa2 = as.numeric(as.character(toj_means$soa2))
  
  gg = ggplot(
    data = toj_means
    , mapping = aes(
      x = soa2
      , y = value
      , shape = base_probe_dist
      , linetype = base_probe_dist
      , group = base_probe_dist
      #   , colour = room
    )
  ) +
    # facet_grid(.~room)+
    # theme(panel.margin = unit(2, "lines"))+
    scale_x_discrete(breaks = unique(toj_means$soa2), limits = unique(toj_means$soa2) )+
#     geom_line()+ 
#     geom_point()+
    geom_smooth(
      method = "glm"
      , method.args = list(family = "binomial")
      , formula = y ~ splines::ns(x,3)
      , se = FALSE
      # , level = 0.95  # 95% confidence interval 
    )+
    labs(x = "Stimulus Onset Asychrony (SOA)", y = "Proportion of Safe Responses")+
    scale_shape_discrete(name = "Attend"
                         , labels = c("Glove", "Base") )+
    scale_linetype_discrete(name = "Attend"
                            , labels = c("Glove", "Base") )+
    guides(colour = FALSE)
  
  Text1 = textGrob(label = paste("Out"))
  Text2 = textGrob(label = paste("Safe"))
  gg = gg + 
    annotation_custom(grob = Text1,  xmin = -200, xmax = -200, ymin = -0.1, ymax = -0.1)+
    annotation_custom(grob = Text2,  xmin = 200, xmax = 200, ymin = -0.1, ymax = -0.1) +
    theme_gray(base_size = 18)
  
  # Code to override clipping
  gg2 <- ggplot_gtable(ggplot_build(gg))
  gg2$layout$clip[gg2$layout$name=="panel"] <- "off"
  grid.draw(gg2)
  
  # return()   
}

# use all data first 
id_all = unique(toj_trials$id)
get_psycho_function(id_all)

# function for pss and slope extraction
get_pss_slope = function(id_list) {
  toj_by_condition = ddply(
    .data = toj_trials[toj_trials$id %in% id_list,]
    , .variables = .(id, base_probe_dist, know_tie_goes_runner)
    , .fun = function(x){
      fit = glm(
        formula = safe~soa2
        , data = x
        , family = binomial
      )
      to_return = data.frame(
        id = x$id[1]
        , pss = -coef(fit)[1]/coef(fit)[2]
        , slope = coef(fit)[2]
        , jnd = qlogis(0.84)/coef(fit)[2]  # mathematicall the same as half the difference between .16 and .84 points 
      )
      return(to_return)
    }
  )
  return(toj_by_condition)
}

 
#------------------------------- PSS -----------------------------------------# 
toj_by_condition = get_pss_slope(id_all)

# descriptive 
TOJ_pss_means = aggregate(pss ~ base_probe_dist, data = toj_by_condition, FUN = mean)
TOJ_pss_SD = aggregate(pss ~ base_probe_dist, data = toj_by_condition, FUN = sd)

# ASSUMPTIONS
# (1) Normality
m = lm(pss ~ base_probe_dist, data = toj_by_condition)
plot(m) # look at second plot: residuals vs. theoreticcal quantiles 
# (2) Equal variance
pss_varz = aggregate(pss ~ base_probe_dist, data = toj_by_condition, FUN = var)
diffreal = pss_varz$pss[1] - pss_varz$pss[2]
# these seem quite close 
# but how close?
# randomization test:
pss = toj_by_condition$pss
varz = NULL
for (i in 1:10000) {
  samp = sample(pss)
  len = length(pss)
  one = samp[1:(len/2)]
  two = samp[(len/2 + 1):len]
  temp = var(one) - var(two)
  varz = c(varz, temp)
}
hist(varz, br = 100)
abline(v = diffreal, col = 'red') # definately within range of what is possible with population equal variances

# pss
t.test(
  formula = pss ~ base_probe_dist
  , dat = toj_by_condition
  , paired = T
  , var.equal = T
)

### with convention knowledge 
toj_pss = aggregate(pss ~ base_probe_dist + know_tie_goes_runner + id, data=toj_by_condition, FUN = mean)
toj_pss$know_tie_goes_runner = as.factor(as.character(toj_pss$know_tie_goes_runner))

# descriptive 
toj_pss_means = aggregate(pss ~ base_probe_dist + know_tie_goes_runner, data = toj_pss, FUN = mean)
toj_pss_SD =  aggregate(pss ~ base_probe_dist + know_tie_goes_runner, data = toj_pss, FUN = sd)

# model
m = aov(
  formula = pss ~ base_probe_dist*know_tie_goes_runner + Error(id/(base_probe_dist))
  , dat = toj_pss
)

# anova
summary(m)
#------------------------------- PSS -----------------------------------------# 


#------------------------------- slope ---------------------------------------# 
# descriptive 
TOJ_slope_means = aggregate(slope ~ base_probe_dist, data = toj_by_condition, FUN = mean)
TOJ_slope_SD = aggregate(slope ~ base_probe_dist, data = toj_by_condition, FUN = sd)

# ASSUMPTIONS
# (1) Normality
m = lm(slope ~ base_probe_dist, data = toj_by_condition)
plot(m) # look at second plot: residuals vs. theoreticcal quantiles 
# (2) Equal variance
slope_varz = aggregate(slope ~ base_probe_dist, data = toj_by_condition, FUN = var)
diffreal = slope_varz$slope[1] - slope_varz$slope[2]
# these seem quite close 
# but how close?
# randomization test:
slope = toj_by_condition$slope
varz = NULL
for (i in 1:10000) {
  samp = sample(slope)
  len = length(slope)
  one = samp[1:(len/2)]
  two = samp[(len/2 + 1):len]
  temp = var(one) - var(two)
  varz = c(varz, temp)
}
hist(varz, br = 100)
abline(v = diffreal, col = 'red') # definately within range of what is possible with population equal variances

# slope
t.test(
  formula = slope ~ base_probe_dist
  , dat = toj_by_condition
  , paired = T
  , var.equal = T
)

### with convention knowledge 
toj_jnd = aggregate(slope ~ base_probe_dist + know_tie_goes_runner + id, data=toj_by_condition, FUN = mean)
toj_jnd$know_tie_goes_runner = as.factor(as.character(toj_jnd$know_tie_goes_runner))

# descriptive 
toj_jnd_means = aggregate(slope ~ base_probe_dist + know_tie_goes_runner, data = toj_jnd, FUN = mean)
toj_jnd_SD =  aggregate(slope ~ base_probe_dist + know_tie_goes_runner, data = toj_jnd, FUN = sd)

# model
m = aov(
  formula = slope ~ base_probe_dist*know_tie_goes_runner + Error(id/(base_probe_dist))
  , dat = toj_jnd
)

# anova
summary(m)
#------------------------------- slope ---------------------------------------# 


#-------------------------------- JND ----------------------------------------# 
# descriptive 
TOJ_jnd_means = aggregate(jnd ~ base_probe_dist, data = toj_by_condition, FUN = mean)
TOJ_jnd_SD = aggregate(jnd ~ base_probe_dist, data = toj_by_condition, FUN = sd)

# # ASSUMPTIONS
# # (1) Normality
# m = lm(jnd ~ base_probe_dist, data = toj_by_condition)
# plot(m) # look at second plot: residuals vs. theoreticcal quantiles 
# # NOT NORMAL AT ALL
# # (2) Equal variance
# jnd_varz = aggregate(jnd ~ base_probe_dist, data = toj_by_condition, FUN = var)
# diffreal = jnd_varz$jnd[1] - jnd_varz$jnd[2]
# # these seem quite close 
# # but how close?
# # randomization test:
# jnd = toj_by_condition$jnd
# varz = NULL
# for (i in 1:10000) {
#   samp = sample(jnd)
#   len = length(jnd)
#   one = samp[1:(len/2)]
#   two = samp[(len/2 + 1):len]
#   temp = var(one) - var(two)
#   varz = c(varz, temp)
# }
# hist(varz, br = 100)
# abline(v = diffreal, col = 'red') # definately within range of what is possible with population equal variances
# 
# # jnd
# t.test(
#   formula = jnd ~ base_probe_dist
#   , dat = toj_by_condition
#   , paired = T
#   , var.equal = T
# )
# 
# ### with convention knowledge 
# toj_jnd = aggregate(jnd ~ base_probe_dist + know_tie_goes_runner + id, data=toj_by_condition, FUN = mean)
# toj_jnd$know_tie_goes_runner = as.factor(as.character(toj_jnd$know_tie_goes_runner))
# 
# descriptive 
toj_jnd_means = aggregate(jnd ~ base_probe_dist + know_tie_goes_runner, data = toj_jnd, FUN = mean)
toj_jnd_SD =  aggregate(jnd ~ base_probe_dist + know_tie_goes_runner, data = toj_jnd, FUN = sd)
# 
# # model
# m = aov(
#   formula = jnd ~ base_probe_dist*know_tie_goes_runner + Error(id/(base_probe_dist))
#   , dat = toj_jnd
# )
# 
# # anova
# summary(m)
#-------------------------------- JND ----------------------------------------# 


#-------------------------------- log JND ------------------------------------# 
toj_by_condition$log_jnd = log(toj_by_condition$jnd)

# descriptive 
TOJ_log_jnd_means = aggregate(log_jnd ~ base_probe_dist, data = toj_by_condition, FUN = mean)
TOJ_log_jnd_SD = aggregate(log_jnd ~ base_probe_dist, data = toj_by_condition, FUN = sd)

# ASSUMPTIONS
# (1) Normality
m = lm(log_jnd ~ base_probe_dist, data = toj_by_condition)
plot(m) # look at second plot: residuals vs. theoreticcal quantiles 
# BETTER
# (2) Equal variance
log_jnd_varz = aggregate(log_jnd ~ base_probe_dist, data = toj_by_condition, FUN = var)
diffreal = log_jnd_varz$log_jnd[1] - log_jnd_varz$log_jnd[2]
# these seem quite close 
# but how close?
# randomization test:
log_jnd = toj_by_condition$log_jnd
varz = NULL
for (i in 1:10000) {
  samp = sample(log_jnd)
  len = length(log_jnd)
  one = samp[1:(len/2)]
  two = samp[(len/2 + 1):len]
  temp = var(one) - var(two)
  varz = c(varz, temp)
}
hist(varz, br = 100)
abline(v = diffreal, col = 'red') # definately within range of what is possible with population equal variances

# log_jnd
t.test(
  formula = log_jnd ~ base_probe_dist
  , dat = toj_by_condition
  , paired = T
  , var.equal = T
)

### with convention knowledge 
toj_log_jnd = aggregate(log_jnd ~ base_probe_dist + know_tie_goes_runner + id, data=toj_by_condition, FUN = mean)
toj_log_jnd$know_tie_goes_runner = as.factor(as.character(toj_log_jnd$know_tie_goes_runner))

# descriptive 
toj_log_jnd_means = aggregate(log_jnd ~ base_probe_dist + know_tie_goes_runner, data = toj_log_jnd, FUN = mean)
toj_log_jnd_SD =  aggregate(log_jnd ~ base_probe_dist + know_tie_goes_runner, data = toj_log_jnd, FUN = sd)

# model
m = aov(
  formula = log_jnd ~ base_probe_dist*know_tie_goes_runner + Error(id/(base_probe_dist))
  , dat = toj_log_jnd
)

# anova
summary(m)
#-------------------------------- log JND ------------------------------------# 