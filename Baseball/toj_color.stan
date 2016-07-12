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
	vector[8] zeros ;
	real neglog2pi;
  neglog2pi <- -log(2.0 * pi()); // log-probability of uniform component (i.e. data invariant)
	for(i in 1:8){
		zeros[i] <- 0 ;
	}
}
parameters{
	// Population Means
	real population_pss_intercept_mean;
	real population_pss_effect_mean; 
	real population_pss_convention_effect_mean;
	real population_pss_attention_convention_interaction_effect_mean;
	real population_logjnd_intercept_mean;
	real population_logjnd_effect_mean;
	real population_logjnd_convention_effect_mean;
	real population_logjnd_attention_convention_interaction_effect_mean; 
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
	real<lower=0,upper=pi()/2> zpopulation_logjnd_intercept_sd;
	real<lower=0,upper=pi()/2> zpopulation_logjnd_effect_sd; 
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
		real pss_intercept_per_id[N_toj]; 
		real pss_effect_per_id[N_toj]; 
		real logjnd_intercept_per_id[N_toj]; 
		real logjnd_effect_per_id[N_toj]; 
		real trial_pss[L_toj]; 
		real trial_logjnd[L_toj]; 
		real population_pss_intercept_sd;
		real population_pss_effect_sd; 
		real population_logjnd_intercept_sd;
		real population_logjnd_effect_sd; 
		
    // Local Inits for Color Wheel
    real population_logit_rho_intercept_sd ;
    real population_log_kappa_intercept_sd ;
    real population_logit_rho_effect_sd ;
    real population_log_kappa_effect_sd ;
    vector[N_color] logRho[2]; // subject-level log-rho for each condition
    vector[N_color] log1mrho_neglog2pi[2]; // subject-level log(1-rho) for each condition
    vector[N_color] kappa[2]; // eventually used in the model block
    vector[N_color] rho[2]; // intermediate value
    vector[N_color] id_logit_rho_intercept; // intermediate value
    vector[N_color] id_log_kappa_intercept; // intermediate value
    vector[N_color] id_logit_rho_effect; // intermediate value
    vector[N_color] id_log_kappa_effect; // intermediate value
		
		// Computations for TOJ
		population_pss_intercept_sd <- tan(zpopulation_pss_intercept_sd) ;
		population_pss_effect_sd <- tan(zpopulation_pss_effect_sd) ;
		population_logjnd_intercept_sd <- tan(zpopulation_logjnd_intercept_sd);
		population_logjnd_effect_sd <- tan(zpopulation_logjnd_effect_sd) ;
		// compute unit-level parameters
		for(n in 1:N_toj){
			pss_intercept_per_id[n] <- 
			beta[n,1]*population_pss_intercept_sd + population_pss_intercept_mean 
			+ population_pss_convention_effect_mean*condition_convention[n]/2 ;
			
			pss_effect_per_id[n] <- (
			beta[n,2]*population_pss_effect_sd + population_pss_effect_mean
			+ population_pss_attention_convention_interaction_effect_mean*condition_convention[n]
			)/2 ;
			
			logjnd_intercept_per_id[n] <- 
			beta[n,3]*population_logjnd_intercept_sd + population_logjnd_intercept_mean 
			+ population_logjnd_convention_effect_mean*condition_convention[n]/2 ;
			
			logjnd_effect_per_id[n] <- (
			beta[n,4]*population_logjnd_effect_sd + population_logjnd_effect_mean
			+ population_logjnd_attention_convention_interaction_effect_mean*condition_convention[n]
			)/2 ;
		}
		// compute trial-level parameters
		for(l in 1:L_toj){
			trial_pss[l] <- pss_intercept_per_id[id_toj[l]] + pss_effect_per_id[id_toj[l]]*condition_toj[l] ;  
			trial_logjnd[l] <- logjnd_intercept_per_id[id_toj[l]] + logjnd_effect_per_id[id_toj[l]]*condition_toj[l]; 
			trial_prob[l] <- Phi_approx((x_toj[l]-trial_pss[l])/exp(trial_logjnd[l])) ;
		}
		
  	// Computations for Color Wheel
    population_logit_rho_intercept_sd <- tan(zpopulation_logit_rho_intercept_sd) ;
    population_log_kappa_intercept_sd <- tan(zpopulation_log_kappa_intercept_sd) ;
    population_logit_rho_effect_sd <- tan(zpopulation_logit_rho_effect_sd) ;
    population_log_kappa_effect_sd <- tan(zpopulation_log_kappa_effect_sd) ;
    // compute unit-level parameters
    for(n in 1:N_color){
      id_logit_rho_intercept[n] <- 
      beta[n,5]*population_logit_rho_intercept_sd + population_logit_rho_intercept_mean
      + population_logit_rho_convention_effect_mean*condition_convention[n]/2 ;
      
      id_log_kappa_intercept[n] <- 
      beta[n,6]*population_log_kappa_intercept_sd + population_log_kappa_intercept_mean 
      + population_log_kappa_convention_effect_mean*condition_convention[n]/2 ;
      
      id_logit_rho_effect[n] <- 
      beta[n,7]*population_logit_rho_effect_sd + population_logit_rho_attention_effect_mean
      + population_logit_rho_attention_convention_interaction_effect_mean*condition_convention[n] ;
      
      id_log_kappa_effect[n] <- 
      beta[n,8]*population_log_kappa_effect_sd + population_log_kappa_attention_effect_mean
      + population_log_kappa_attention_convention_interaction_effect_mean*condition_convention[n] ;
      
      // compute the transforms 
      rho[1,n] <- inv_logit( id_logit_rho_intercept[n] - id_logit_rho_effect[n]/2 );  
      rho[2,n] <- inv_logit( id_logit_rho_intercept[n] + id_logit_rho_effect[n]/2 );  
      kappa[1,n] <- exp( id_log_kappa_intercept[n] - id_log_kappa_effect[n]/2 );
      kappa[2,n] <- exp( id_log_kappa_intercept[n] + id_log_kappa_effect[n]/2 );
      logRho[1,n] <- log(rho[1,n]);
      logRho[2,n] <- log(rho[2,n]);
      log1mrho_neglog2pi[1,n] <- log1m(rho[1,n]) + neglog2pi;
      log1mrho_neglog2pi[2,n] <- log1m(rho[2,n]) + neglog2pi;
    }
  	// compute trial-level parameters
    for (l in 1:L_color){
      p[l] <- log_sum_exp( 
                            logRho[condition_color[l],id_color[l]] 
                            + von_mises_log(
                            y_color[l]
                            ,pi()
                            ,kappa[condition_color[l],id_color[l]]
                            )   
                            , log1mrho_neglog2pi[condition_color[l],id_color[l]]
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
	population_logjnd_intercept_mean ~ normal(-1,.5) ;
	population_logjnd_effect_mean ~ normal(0,1) ;
	population_logjnd_convention_effect_mean ~ normal(0,1) ;
	population_logjnd_attention_convention_interaction_effect_mean ~ normal(0,1) ;
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
  increment_log_prob(p) ;
}
