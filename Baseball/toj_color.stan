data{
  // NOTE: N_toj = N_color
	int N_toj ; // number of participants for toj task
	int L_toj ; // number of observations for toj task
	int y_toj[L_toj] ; // 'Left First' judgement for each observation
	real x_toj[L_toj] ; // SOA for each observation
	int id_toj[L_toj] ; // participant ID for each observation for toj task 
	int condition_toj[L_toj] ; // attention condition ('Attend Left' or 'Attend Right') for each observation
	int<lower=0> N_color; // number of participants for color task                   
  int<lower=0> L_color; // number of observations for color task                       
  int<lower=0> id_color[L_color]; // participant ID for each observation for color task               
  int<lower=1,upper=2> condition_color[L_color]; // attention condition ('Attended' or 'Unattended') for each observation 
  int condition_convention[N_toj]; // convention condition ('Know' or 'Don't Know') for each participant
  real<lower=0,upper=2*pi()> y_color[L_color]; // color wheel responses (in radians) for each observation
}
transformed data{
  int<lower=1> K;
	vector[8] zeros ;
	real neglog2pi;
  neglog2pi = -log(2.0 * pi()); // log-probability of uniform component (i.e. data invariant)
	for(i in 1:8){
		zeros[i] = 0 ;
	}
	K = max(condition_color);
}
parameters{
	// Population Means
	real population_pss_intercept_mean;
	real population_pss_effect_mean; 
	real population_pss_convention_effect_mean;
	real population_pss_attention_convention_interaction_effect_mean;
	real population_log_jnd_intercept_mean;
	real population_log_jnd_effect_mean;
	real population_log_jnd_convention_effect_mean;
	real population_log_jnd_attention_convention_interaction_effect_mean; 
  real population_logit_rho_intercept_mean;
  real population_logit_rho_attention_effect_mean;
  real population_logit_rho_convention_effect_mean;
  real population_logit_rho_attention_convention_interaction_effect_mean;
  real population_log_kappa_intercept_mean;
  real population_log_kappa_convention_effect_mean;
  real population_log_kappa_attention_effect_mean;
  real population_log_kappa_attention_convention_interaction_effect_mean;
  
  // Population SDs
  // 'Tan Trick' --> equivalent to half-cauchy with lower bound at zero
	real<lower=0,upper=pi()/2> zpopulation_pss_intercept_sd;
	real<lower=0,upper=pi()/2> zpopulation_pss_effect_sd; 
	real<lower=0,upper=pi()/2> zpopulation_log_jnd_intercept_sd;
	real<lower=0,upper=pi()/2> zpopulation_log_jnd_effect_sd; 
  real<lower=0,upper=pi()/2> zpopulation_logit_rho_intercept_sd;
  real<lower=0,upper=pi()/2> zpopulation_logit_rho_effect_sd;
  real<lower=0,upper=pi()/2> zpopulation_log_kappa_intercept_sd;
  real<lower=0,upper=pi()/2> zpopulation_log_kappa_effect_sd;
  
	corr_matrix[8] cor ; 
	// dummy variable for 'Matt Trick'
	vector[8] beta[N_toj];
}

transformed parameters{
	real trial_prob[L_toj];
  vector[L_color] p; // for storing log-probabilities
	{
	  // Local Inits for TOJ
		real id_pss_intercept[N_toj]; 
		real id_pss_effect[N_toj]; 
		real id_log_jnd_intercept[N_toj]; 
		real id_log_jnd_effect[N_toj]; 
		real trial_pss[L_toj]; 
		real trial_log_jnd[L_toj]; 
		real population_pss_intercept_sd;
		real population_pss_effect_sd; 
		real population_log_jnd_intercept_sd;
		real population_log_jnd_effect_sd; 
		
    // Local Inits for Color Wheel
    real population_logit_rho_intercept_sd ;
    real population_log_kappa_intercept_sd ;
    real population_logit_rho_effect_sd ;
    real population_log_kappa_effect_sd ;
    vector[N_color] id_logit_rho_intercept ;
    vector[N_color] id_logit_rho_effect ;
    vector[N_color] id_logit_rho[K] ;
    vector[N_color] id_log_kappa_intercept ;
    vector[N_color] id_log_kappa_effect ;
    vector[N_color] id_log_kappa[K] ;
    // useful transformations
    vector[N_color] id_kappa[K] ;
    vector[N_color] id_inv_sqrt_kappa[K] ;
    vector[N_color] id_rho[K] ;
		
		// Computations for TOJ
		population_pss_intercept_sd = tan(zpopulation_pss_intercept_sd) ;
		population_pss_effect_sd = tan(zpopulation_pss_effect_sd) ;
		population_log_jnd_intercept_sd = tan(zpopulation_log_jnd_intercept_sd);
		population_log_jnd_effect_sd = tan(zpopulation_log_jnd_effect_sd) ;
		// compute unit-level parameters
		for(n in 1:N_toj){
			id_pss_intercept[n] = 
			beta[n,1]*population_pss_intercept_sd + population_pss_intercept_mean 
			+ population_pss_convention_effect_mean*condition_convention[n]/2 ;
			
			id_pss_effect[n] = (
			beta[n,2]*population_pss_effect_sd + population_pss_effect_mean
			+ population_pss_attention_convention_interaction_effect_mean*condition_convention[n]/2
			)/2 ;
			
			id_log_jnd_intercept[n] = 
			beta[n,3]*population_log_jnd_intercept_sd + population_log_jnd_intercept_mean 
			+ population_log_jnd_convention_effect_mean*condition_convention[n]/2 ;
			
			id_log_jnd_effect[n] = (
			beta[n,4]*population_log_jnd_effect_sd + population_log_jnd_effect_mean
			+ population_log_jnd_attention_convention_interaction_effect_mean*condition_convention[n]/2 
			)/2 ;
		}
		// compute trial-level parameters
		for(l in 1:L_toj){
			trial_pss[l] = id_pss_intercept[id_toj[l]] + id_pss_effect[id_toj[l]]*condition_toj[l] ;  
			trial_log_jnd[l] = id_log_jnd_intercept[id_toj[l]] + id_log_jnd_effect[id_toj[l]]*condition_toj[l]; 
			trial_prob[l] = Phi_approx((x_toj[l]-trial_pss[l])/exp(trial_log_jnd[l])) ;
		}
		
  	// Computations for Color Wheel
    population_logit_rho_intercept_sd = tan(zpopulation_logit_rho_intercept_sd) ;
    population_log_kappa_intercept_sd = tan(zpopulation_log_kappa_intercept_sd) ;
    population_logit_rho_effect_sd = tan(zpopulation_logit_rho_effect_sd) ;
    population_log_kappa_effect_sd = tan(zpopulation_log_kappa_effect_sd) ;
    // compute unit-level parameters
    for(n in 1:N_color){
      id_logit_rho_intercept[n] = 
      beta[n,5]*population_logit_rho_intercept_sd + population_logit_rho_intercept_mean
      + population_logit_rho_convention_effect_mean*condition_convention[n]/2 ;
      
      id_log_kappa_intercept[n] = 
      beta[n,6]*population_log_kappa_intercept_sd + population_log_kappa_intercept_mean 
      + population_log_kappa_convention_effect_mean*condition_convention[n]/2 ;
      
      id_logit_rho_effect[n] = 
      beta[n,7]*population_logit_rho_effect_sd + population_logit_rho_attention_effect_mean
      + population_logit_rho_attention_convention_interaction_effect_mean*condition_convention[n]/2 ;
      
      id_log_kappa_effect[n] = 
      beta[n,8]*population_log_kappa_effect_sd + population_log_kappa_attention_effect_mean
      + population_log_kappa_attention_convention_interaction_effect_mean*condition_convention[n]/2 ;
      
      id_logit_rho[1][n] = id_logit_rho_intercept[n] + id_logit_rho_effect[n]/2 ;
      id_logit_rho[2][n] = id_logit_rho_intercept[n] - id_logit_rho_effect[n]/2 ;
      id_log_kappa[1][n] = id_log_kappa_intercept[n] + id_log_kappa_effect[n]/2 ;
      id_log_kappa[2][n] = id_log_kappa_intercept[n] - id_log_kappa_effect[n]/2 ;
    }
    // compute the transforms
    for(i in 1:K){
      id_kappa[i] = exp(id_log_kappa[i]) ;
      for(n in 1:N_color){
        id_inv_sqrt_kappa[i][n] = sqrt(1/id_kappa[i][n]) ;
        id_rho[i][n] = inv_logit(id_logit_rho[i][n]) ;
      }
    }
  	// compute trial-level parameters
    for(l in 1:L_color){
        p[l] = log_mix( //mixture for angle
                      id_rho[condition_color[l]][id_color[l]]
                      , von_mises_lpdf(
                                       y_color[l] | pi()
                                       , id_kappa[condition_color[l]][id_color[l]]
                                       )
                      , neglog2pi
                      ) ;
    }
	}
}
model{
  // Population Means
	population_pss_intercept_mean ~ normal(0,1) ;
	population_pss_effect_mean ~ normal(0,1) ;
	population_pss_convention_effect_mean ~ normal(0,1) ;
	population_pss_attention_convention_interaction_effect_mean ~ normal(0,1) ; 
	population_log_jnd_intercept_mean ~ normal(-1,.5) ;
	population_log_jnd_effect_mean ~ normal(0,1) ;
	population_log_jnd_convention_effect_mean ~ normal(0,1) ;
	population_log_jnd_attention_convention_interaction_effect_mean ~ normal(0,1) ;
  population_logit_rho_intercept_mean ~ normal(3,3);
  population_logit_rho_attention_effect_mean ~ normal(0,3) ;
  population_logit_rho_convention_effect_mean ~ normal(0,3);
  population_logit_rho_attention_convention_interaction_effect_mean ~ normal(0,3);
  population_log_kappa_intercept_mean ~ normal(3,3);
  population_log_kappa_attention_effect_mean ~ normal(0,3) ;
  population_log_kappa_convention_effect_mean ~ normal(0,3);
  population_log_kappa_attention_convention_interaction_effect_mean ~ normal(0,3);

  cor ~ lkj_corr(4) ;  // prior on correlation matrix 
	
	#sample the betas from standard multivariate normal
	for(n in 1:N_toj){
	  beta[n] ~ multi_student_t(1,zeros,cor) ;
	}
	
	y_toj ~ bernoulli(trial_prob) ;
	//update the log-probability from p (defined in transformed parameters)
  target += p ;
}
