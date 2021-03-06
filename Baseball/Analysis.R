#### Libraries ####
library(plyr)
library(ggplot2)
library(grid)
library(rstan)

setwd("~/Documents/TOJ/Baseball/baseballtojdata")

##########################################
####          Data Import             #### 
##########################################

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



##########################################
####             TOJ                  #### 
##########################################
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

# ### save toj trials for posterior predictive check
save(toj_trials, file = "../toj_trials.Rdata")


### Graph participant-wise psychometric functions 
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
  geom_smooth(
    method = "glm"
    , method.args = list(family = binomial(link="probit"))
    , se = FALSE
  ) +
  geom_point()


#### Group functions for Ray ####
toj_means_by_id_by_condition2 = ddply(
  .data = toj_trials
  , .variables = .(base_probe_dist, know_tie_goes_runner, soa2)
  , .fun = function(x){
    to_return = data.frame(
      value = mean(x$safe)
      , factor = paste(unique(x$base_probe_dist), unique(x$know_tie_goes_runner))
    )
    return(to_return)
  }
)
toj_means_by_id_by_condition2$soa2 = as.numeric(as.character(toj_means_by_id_by_condition2$soa2))

ggplot(
  data = toj_means_by_id_by_condition2
  , mapping = aes(
    x = soa2
    , y =  value
    , shape = factor
    , linetype = base_probe_dist
    , group = factor
    , colour = know_tie_goes_runner
  )
)+
  geom_smooth(
    method = "glm"
    , method.args = list(family = binomial(link="probit"))
    , se = FALSE
  ) +
  geom_point()



##########################################
####             Color                #### 
##########################################
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

color_trials$attended = FALSE
color_trials$attended[ (color_trials$base_probe_dist == 0.8 & color_trials$probe_location == "base") | (color_trials$base_probe_dist == 0.2 & color_trials$probe_location == "glove")] = TRUE

color_trials$color_diff_radians = color_trials$color_diff*pi/180

### save color trials for posterior predicitve check 
save(color_trials, file = "../color_trials.Rdata")



##########################################
####             Stan                 #### 
##########################################
toj_color_data_for_stan = list(
  N_toj = length(unique(toj_trials$id))
  , L_toj = nrow(toj_trials)
  , y_toj = as.numeric(toj_trials$safe)  
  , x_toj = (as.numeric(toj_trials$soa2))/250  # we normalize soas, and therefore pss
  , id_toj = as.numeric(factor(toj_trials$id))
  , condition_toj = ifelse(toj_trials$glove_probe_dist==.2,-1,1)  # base is -1 and glove is +1
  , N_color = length(unique(color_trials$id))
  , L_color = nrow(color_trials)
  , id_color = as.numeric(factor(color_trials$id))
  , condition_color = as.numeric(as.factor(color_trials$attended)) # TRUE is 2, FALSE is 1 
  , condition_convention = ifelse(aggregate(know_tie_goes_runner~id,data = toj_trials, FUN =unique)$know_tie_goes_runner, -1, 1) # know ('safe' bias) is -1...
  , y_color = pi+degree_to_rad(color_trials$color_diff)  # want from 0 to 360 instead of -180 to 180
)
  
toj_color_model = stan_model(
  file = '../../EndogenousVisualPriorEntry-BayesianHierarchicalModel/Baseball/toj_color.stan'
)

toj_color_post = sampling(
  object = toj_color_model
  , data = toj_color_data_for_stan
#   , iter = 1e4*2
#   , chains = 8
#   , cores = 8
  , iter = 100
  , chains = 4
  , cores = 4
  , pars = c('trial_prob', 'p')
  , include = FALSE
)
print(toj_color_post)
# save(toj_color_post, file = "toj_color_post_June16th2016")
save(toj_color_post, file = "TESTCLEANCODE_toj_color_post_Aug162016")





  
  