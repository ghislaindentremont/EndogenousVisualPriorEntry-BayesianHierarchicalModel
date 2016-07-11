data{
	int<lower=0> N_toj ;
	int<lower=0> L_toj ;
	int<lower=0,upper=1> y_toj[L_toj] ;
	real x_toj[L_toj] ;
	int id_toj[L_toj] ;
	int<lower=-1,upper=1> condition_toj[L_toj] ;
	int<lower=0> N_color;                      
  int<lower=0> L_color;                         
  int<lower=0> id_color[L_color];                 
  int<lower=1,upper=2> condition_color[L_color];   
  // int<lower=-1,upper=1> condition_initial_bias[N_toj];
  int<lower=-1,upper=1> condition_probe[N_toj];
  int<lower=-1,upper=1> condition_judgement_type[N_toj];
  real<lower=0,upper=2*pi()> y_color[L_color]; 
}
transformed data{
  int<lower=1> K ;
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
	real population_pss_intercept_mean ;
	real population_pss_attention_effect_mean ;
	// real population_pss_initial_bias_effect_mean;
	real population_pss_judgement_type_effect_mean;
	real population_pss_probe_duration_effect_mean;
	// real population_pss_attention_initial_bias_interaction_effect_mean;
	real population_pss_attention_judgement_type_interaction_effect_mean;	
	real population_pss_attention_probe_duration_interaction_effect_mean;
	real population_log_jnd_intercept_mean ;
	real population_log_jnd_attention_effect_mean ;
	// real population_log_jnd_initial_bias_effect_mean;
	real population_log_jnd_judgement_type_effect_mean;
	real population_log_jnd_probe_duration_effect_mean;
	// real population_log_jnd_attention_initial_bias_interaction_effect_mean;
	real population_log_jnd_attention_judgement_type_interaction_effect_mean;	
	real population_log_jnd_attention_probe_duration_interaction_effect_mean;
	real population_logit_rho_intercept_mean;
  real population_logit_rho_attention_effect_mean;
  real population_logit_rho_probe_duration_effect_mean;
  real population_logit_rho_attention_probe_duration_interaction_effect_mean;
  // real population_logit_rho_initial_bias_effect_mean;
  // real population_logit_rho_attention_initial_bias_interaction_effect_mean;
  real population_log_kappa_intercept_mean;
  real population_log_kappa_attention_effect_mean;
  real population_log_kappa_probe_duration_effect_mean;
  real population_log_kappa_attention_probe_duration_interaction_effect_mean;
  // real population_log_kappa_initial_bias_effect_mean;
  // real population_log_kappa_attention_initial_bias_interaction_effect_mean;

	// Population SDs
	real<lower=0> population_pss_intercept_sd ;
	real<lower=0> population_pss_effect_sd ; 
	real<lower=0> population_log_jnd_intercept_sd ;
	real<lower=0> population_log_jnd_effect_sd ; 
  real<lower=0> population_logit_rho_intercept_sd;
  real<lower=0> population_log_kappa_intercept_sd;
  real<lower=0> population_logit_rho_effect_sd;
  real<lower=0> population_log_kappa_effect_sd;
  
	// correlation
	corr_matrix[8] cor ; 
	// dummy variable for matt trick
	vector[8] beta[N_toj]; // N_toj = N_color
}

transformed parameters{
	real trial_prob[L_toj] ;
  vector[L_color] p ; // for storing log-probabilities
	{
	  // local inits for TOJ
		real id_pss_intercept[N_toj] ; 
		real id_pss_effect[N_toj] ; 
		real id_log_jnd_intercept[N_toj] ; 
		real id_log_jnd_effect[N_toj] ; 
		real trial_pss[L_toj] ; 
		real trial_log_jnd[L_toj] ; 
		
    // local inits for color wheel
    vector[N_color] id_logit_rho_intercept ;
    vector[N_color] id_logit_rho_effect ;
    vector[N_color] id_logit_rho[K] ;
    vector[N_color] id_log_kappa_intercept ;
    vector[N_color] id_log_kappa_effect ;
    vector[N_color] id_log_kappa[K] ;
    //useful transformations
    vector[N_color] id_kappa[K] ;
    vector[N_color] id_inv_sqrt_kappa[K] ;
    vector[N_color] id_rho[K] ;
		
		// Computations for TOJ
		for(n in 1:N_toj){
			id_pss_intercept[n] = 
			beta[n,1]*population_pss_intercept_sd + population_pss_intercept_mean
			// + population_pss_initial_bias_effect_mean* condition_initial_bias[n]/2
			+ population_pss_judgement_type_effect_mean * condition_judgement_type[n]/2
			+ population_pss_probe_duration_effect_mean * condition_probe[n]/2 ;

			id_pss_effect[n] = ( 
			 beta[n,2]*population_pss_effect_sd + population_pss_attention_effect_mean
			// + population_pss_attention_initial_bias_interaction_effect_mean * condition_initial_bias[n] 
			+ population_pss_attention_judgement_type_interaction_effect_mean * condition_judgement_type[n]
			+ population_pss_attention_probe_duration_interaction_effect_mean * condition_probe[n]
			)/2 ;

		  id_log_jnd_intercept[n] = 
		  beta[n,3]*population_log_jnd_intercept_sd + population_log_jnd_intercept_mean
			// + population_log_jnd_initial_bias_effect_mean* condition_initial_bias[n]/2 
			+ population_log_jnd_judgement_type_effect_mean * condition_judgement_type[n]/2
			+ population_log_jnd_probe_duration_effect_mean * condition_probe[n]/2 ;
			
			id_log_jnd_effect[n] = ( 
			beta[n,4]*population_log_jnd_effect_sd + population_log_jnd_attention_effect_mean
			// + population_log_jnd_attention_initial_bias_interaction_effect_mean * condition_initial_bias[n] 
			+ population_log_jnd_attention_judgement_type_interaction_effect_mean * condition_judgement_type[n] 
			+ population_log_jnd_attention_probe_duration_interaction_effect_mean * condition_probe[n]
			)/2 ;
		}
		// compute trial-level parameters
		for(l in 1:L_toj){
			trial_pss[l] = id_pss_intercept[id_toj[l]] + id_pss_effect[id_toj[l]]*condition_toj[l];  
			trial_log_jnd[l] = id_log_jnd_intercept[id_toj[l]] + id_log_jnd_effect[id_toj[l]]*condition_toj[l]; 
			trial_prob[l] = Phi_approx((x_toj[l]-trial_pss[l])/exp(trial_log_jnd[l])) ;
		}
		
  	// Computations for Color wheel
    for(n in 1:N_color){
      id_logit_rho_intercept[n] = 
      beta[n,5]*population_logit_rho_intercept_sd + population_logit_rho_intercept_mean
      + population_logit_rho_probe_duration_effect_mean * condition_probe[n]/2 ;
      // + population_logit_rho_initial_bias_effect_mean * condition_initial_bias[n]/2 ;
      
      id_logit_rho_effect[n] = 
      beta[n,7]*population_logit_rho_effect_sd + population_logit_rho_attention_effect_mean 
      + population_logit_rho_attention_probe_duration_interaction_effect_mean * condition_probe[n] ;
      // + population_logit_rho_attention_initial_bias_interaction_effect_mean * condition_initial_bias[n] ;
      
      id_log_kappa_intercept[n] = 
      beta[n,6]*population_log_kappa_intercept_sd + population_log_kappa_intercept_mean 
      + population_log_kappa_probe_duration_effect_mean * condition_probe[n]/2 ;
      // + population_log_kappa_initial_bias_effect_mean * condition_initial_bias[n]/2 ;
      
      id_log_kappa_effect[n] = 
      beta[n,8]*population_log_kappa_effect_sd + population_log_kappa_attention_effect_mean
      + population_log_kappa_attention_probe_duration_interaction_effect_mean * condition_probe[n] ;
      // + population_log_kappa_attention_initial_bias_interaction_effect_mean * condition_initial_bias[n] ;
      
      // compute unit-level parameters
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
	population_pss_intercept_mean ~ student_t(4,0,1);
	population_pss_attention_effect_mean ~ student_t(4,0,1);
	// population_pss_initial_bias_effect_mean ~ student_t(4,0,1);
	population_pss_judgement_type_effect_mean ~ student_t(4,0,1);
	population_pss_probe_duration_effect_mean ~ student_t(4,0,1);
	// population_pss_attention_initial_bias_interaction_effect_mean ~ student_t(4,0,1);
	population_pss_attention_judgement_type_interaction_effect_mean ~ student_t(4,0,1);	
	population_pss_attention_probe_duration_interaction_effect_mean ~ student_t(4,0,1);
	population_log_jnd_intercept_mean ~ student_t(4,-1,.5);
	population_log_jnd_attention_effect_mean ~ student_t(4,0,1);
	// population_log_jnd_initial_bias_effect_mean ~ student_t(4,0,1);
	population_log_jnd_judgement_type_effect_mean ~ student_t(4,0,1);
	population_log_jnd_probe_duration_effect_mean ~ student_t(4,0,1);
	// population_log_jnd_attention_initial_bias_interaction_effect_mean ~ student_t(4,0,1);
	population_log_jnd_attention_judgement_type_interaction_effect_mean ~ student_t(4,0,1);	
	population_log_jnd_attention_probe_duration_interaction_effect_mean ~ student_t(4,0,1);
	population_logit_rho_intercept_mean ~ student_t(4,3,3);
  population_logit_rho_attention_effect_mean ~ student_t(4,0,2);
  population_logit_rho_probe_duration_effect_mean ~ student_t(4,0,2);
  population_logit_rho_attention_probe_duration_interaction_effect_mean ~ student_t(4,0,2);
  // population_logit_rho_initial_bias_effect_mean ~ student_t(4,0,2);
  // population_logit_rho_attention_initial_bias_interaction_effect_mean ~ student_t(4,0,2);
  population_log_kappa_intercept_mean ~ student_t(4,3,3);
  population_log_kappa_attention_effect_mean ~ student_t(4,0,2);
  population_log_kappa_probe_duration_effect_mean ~ student_t(4,0,2);
  population_log_kappa_attention_probe_duration_interaction_effect_mean ~ student_t(4,0,2);
  // population_log_kappa_initial_bias_effect_mean ~ student_t(4,0,2);
  // population_log_kappa_attention_initial_bias_interaction_effect_mean ~ student_t(4,0,2);

	// Population SDs
	population_pss_intercept_sd ~ student_t(4,0,1);
	population_pss_effect_sd ~ student_t(4,0,1); 
	population_log_jnd_intercept_sd ~ student_t(4,0,1);
	population_log_jnd_effect_sd ~ student_t(4,0,1); 
  population_logit_rho_intercept_sd ~ student_t(4,0,2);
  population_log_kappa_intercept_sd ~ student_t(4,0,2);
  population_logit_rho_effect_sd ~ student_t(4,0,2);
  population_log_kappa_effect_sd ~ student_t(4,0,2);

  // prior on correlation matrix 
  cor ~ lkj_corr(4) ;

	// sample the betas from standard multivariate normal
	for(n in 1:N_toj){
	  beta[n] ~ multi_student_t(4,zeros,cor) ;
	}
	
	y_toj ~ bernoulli(trial_prob) ;
  //update the log-probability from p (defined in transformed parameters)
  target += p ;//used to be: increment_log_prob(p) ;
}
