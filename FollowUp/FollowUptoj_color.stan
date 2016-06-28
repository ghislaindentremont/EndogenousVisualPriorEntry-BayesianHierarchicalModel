data{
	int<lower=0> N_toj ;
	int<lower=0> L_toj ;
	int<lower=0,upper=1> y_toj[L_toj] ;
	real x_toj[L_toj] ;
	int id_toj[L_toj] ;
	int<lower=-1,upper=1> condition_toj[L_toj] ;
	int<lower=0> N_color;                      
  int<lower=0> L_color;                         
  int<lower=0> unit_color[L_color];                 
  int<lower=1,upper=2> condition_color[L_color];   
  // int<lower=-1,upper=1> condition_initial_bias[N_toj];
  int<lower=-1,upper=1> condition_judgement_type[N_toj];
  real<lower=0,upper=2*pi()> y_color[L_color]; 
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
	// pss population means
	real population_pss_intercept_mean ;
	
	real population_logjnd_intercept_mean ;
	
	// main effects 
	real population_pss_attention_effect_mean ;
	// real population_pss_initial_bias_effect_mean;
	real population_pss_judgement_type_effect_mean;

	real population_logjnd_attention_effect_mean ;
	// real population_logjnd_initial_bias_effect_mean;
	real population_logjnd_judgement_type_effect_mean;
	
	// two-way interactions (within vs. between)
	// real population_pss_attention_initial_bias_interaction_effect_mean;
	real population_pss_attention_judgement_type_interaction_effect_mean;	
	
	// real population_logjnd_attention_initial_bias_interaction_effect_mean;
	real population_logjnd_attention_judgement_type_interaction_effect_mean;	

	//population sds
	real<lower=0,upper=pi()/2> zpopulation_pss_intercept_sd ;
	real<lower=0,upper=pi()/2> zpopulation_pss_effect_sd ; 
	real<lower=0,upper=pi()/2> zpopulation_logjnd_intercept_sd ;
	real<lower=0,upper=pi()/2> zpopulation_logjnd_effect_sd ; 
	
	//population parameters
  real logitRhoMean;
  
  real logKappaMean;
  
  // main effects
  real logitRhoAttentionEffectMean;
  // real logitRhoInitialBiasEffectMean;
  // real logitRhoJudgementTypeEffectMean;
  
  real logKappaAttentionEffectMean;
  // real logKappaInitialBiasEffectMean;
  // real logKappaJudgementTypeEffectMean;
  
  // two-way interations
  // real logitRhoAttentionInitialBiasEffectMean;
  // real logitRhoAttentionJudgementTypeEffectMean;

  // real logKappaAttentionInitialBiasEffectMean;
  // real logKappaAttentionJudgementTypeEffectMean;
 
  // SDs for population parameters
  real<lower=0,upper=pi()/2> zlogitRhoSD;
  real<lower=0,upper=pi()/2> zlogKappaSD;
  real<lower=0,upper=pi()/2> zlogitRhoEffectSD;
  real<lower=0,upper=pi()/2> zlogKappaEffectSD;
  
	// correlation
	corr_matrix[8] cor ; 
	//dummy variable for matt trick
	vector[8] beta[N_toj]; // N_toj = N_color
}

transformed parameters{
	real trial_prob[L_toj] ;
  vector[L_color] p ; // for storing log-probabilities
	{
	  // local inits for TOJ
		real pss_intercept_per_id[N_toj] ; 
		real pss_effect_per_id[N_toj] ; 
		real logjnd_intercept_per_id[N_toj] ; 
		real logjnd_effect_per_id[N_toj] ; 
		real trial_pss[L_toj] ; 
		real trial_logjnd[L_toj] ; 
		real population_pss_intercept_sd ;
		real population_pss_effect_sd ; 
		real population_logjnd_intercept_sd ;
		real population_logjnd_effect_sd ; 
		
    //local inits for color wheel
    real logitRhoSD ;
    real logKappaSD ;
    real logitRhoEffectSD ;
    real logKappaEffectSD ;
    vector[N_color] logRho[2]; //subject-level log-rho for each condition
    vector[N_color] log1mrho_neglog2pi[2]; //subject-level log(1-rho) for each condition
    vector[N_color] kappa[2]; //eventually used in the model block
    vector[N_color] rho[2]; //intermediate value
    vector[N_color] logitRho; //intermediate value
    vector[N_color] logKappa; //intermediate value
    vector[N_color] logitRhoEffect; //intermediate value
    vector[N_color] logKappaEffect; //intermediate value
		
		// computations for TOJ
		population_pss_intercept_sd <- tan(zpopulation_pss_intercept_sd) ;
		population_pss_effect_sd <- tan(zpopulation_pss_effect_sd) ;
		population_logjnd_intercept_sd <- tan(zpopulation_logjnd_intercept_sd);
		population_logjnd_effect_sd <- tan(zpopulation_logjnd_effect_sd) ;
		for(this_id in 1:N_toj){
			pss_intercept_per_id[this_id] <- beta[this_id,1]*population_pss_intercept_sd + population_pss_intercept_mean
			// + population_pss_initial_bias_effect_mean* condition_initial_bias[this_id]/2
			+ population_pss_judgement_type_effect_mean * condition_judgement_type[this_id]/2 ;

			pss_effect_per_id[this_id] <- ( beta[this_id,2]*population_pss_effect_sd + population_pss_attention_effect_mean
			// + population_pss_attention_initial_bias_interaction_effect_mean * condition_initial_bias[this_id] 
			+ population_pss_attention_judgement_type_interaction_effect_mean * condition_judgement_type[this_id] )/2 ;

		  logjnd_intercept_per_id[this_id] <- beta[this_id,3]*population_logjnd_intercept_sd + population_logjnd_intercept_mean
			// + population_logjnd_initial_bias_effect_mean* condition_initial_bias[this_id]/2 
			+ population_logjnd_judgement_type_effect_mean * condition_judgement_type[this_id]/2 ;
			
			logjnd_effect_per_id[this_id] <- ( beta[this_id,4]*population_logjnd_effect_sd + population_logjnd_attention_effect_mean
			// + population_logjnd_attention_initial_bias_interaction_effect_mean * condition_initial_bias[this_id] 
			+ population_logjnd_attention_judgement_type_interaction_effect_mean * condition_judgement_type[this_id] )/2 ;
		}
		for(this_obs in 1:L_toj){
			trial_pss[this_obs] <- pss_intercept_per_id[id_toj[this_obs]] + pss_effect_per_id[id_toj[this_obs]]*condition_toj[this_obs];  
			trial_logjnd[this_obs] <- logjnd_intercept_per_id[id_toj[this_obs]] + logjnd_effect_per_id[id_toj[this_obs]]*condition_toj[this_obs]; 
			trial_prob[this_obs] <- Phi_approx((x_toj[this_obs]-trial_pss[this_obs])/exp(trial_logjnd[this_obs])) ;
		}
		
  	// computations for color wheel
    logitRhoSD <- tan(zlogitRhoSD) ;
    logKappaSD <- tan(zlogKappaSD) ;
    logitRhoEffectSD <- tan(zlogitRhoEffectSD) ;
    logKappaEffectSD <- tan(zlogKappaEffectSD) ;
    // compute unit-level parameters
    for(n in 1:N_color){
      logitRho[n] <- logitRhoMean + beta[n,5]*logitRhoSD ;
      // + logitRhoInitialBiasEffectMean * condition_initial_bias[n]/2 ;
      // + logitRhoJudgementTypeEffectMean * condition_initial_bias[n]/2 ;
      
      logKappa[n] <- logKappaMean + beta[n,6]*logKappaSD ;
      // + logKappaInitialBiasEffectMean * condition_initial_bias[n]/2
      // + logKappaJudgementTypeEffectMean * condition_initial_bias[n]/2 ;
      
      logitRhoEffect[n] <- logitRhoAttentionEffectMean + beta[n,7]*logitRhoEffectSD ;
      // + logitRhoAttentionInitialBiasEffectMean * condition_initial_bias[n]
      // + logitRhoAttentionJudgementTypeEffectMean * condition_judgement_type[n] ;
      
      logKappaEffect[n] <- logKappaAttentionEffectMean + beta[n,8]*logKappaEffectSD ;
      // + logKappaAttentionInitialBiasEffectMean * condition_initial_bias[n]
      // + logKappaAttentionJudgementTypeEffectMean * condition_judgement_type[n] ;
      
      rho[1,n] <- inv_logit( logitRho[n] - logitRhoEffect[n]/2 ) ; 
      rho[2,n] <- inv_logit( logitRho[n] + logitRhoEffect[n]/2 ) ; 
      kappa[1,n] <- exp( logKappa[n] - logKappaEffect[n]/2  ) ;
      kappa[2,n] <- exp( logKappa[n] + logKappaEffect[n]/2 );
      logRho[1,n] <- log(rho[1,n]);
      logRho[2,n] <- log(rho[2,n]);
      log1mrho_neglog2pi[1,n] <- log1m(rho[1,n]) + neglog2pi;
      log1mrho_neglog2pi[2,n] <- log1m(rho[2,n]) + neglog2pi;
    }
  	// compute trial-level parameters
    for (i in 1:L_color){
      if(kappa[condition_color[i],unit_color[i]]<10){
        p[i] <- log_sum_exp(  
                            logRho[condition_color[i],unit_color[i]] + von_mises_log(y_color[i],pi(),kappa[condition_color[i],unit_color[i]])   
                            , log1mrho_neglog2pi[condition_color[i],unit_color[i]]
        ) ;
      }else{
        p[i] <- log_sum_exp(  
                            logRho[condition_color[i],unit_color[i]] + normal_log(y_color[i],pi(),1/sqrt(kappa[condition_color[i],unit_color[i]]))   
                            , log1mrho_neglog2pi[condition_color[i],unit_color[i]]
        ) ;
      }
    }
	}
}
model{
  //population means
	population_pss_intercept_mean ~ normal(0,1);
	
	population_logjnd_intercept_mean ~ normal(-1,.5);
	
	// main effects 
	population_pss_attention_effect_mean ~ normal(0,1);
	// population_pss_initial_bias_effect_mean~ normal(0,1);
	population_pss_judgement_type_effect_mean~ normal(0,1);

	population_logjnd_attention_effect_mean ~ normal(0,1);
	// population_logjnd_initial_bias_effect_mean~ normal(0,1);
	population_logjnd_judgement_type_effect_mean~ normal(0,1);
	
	// two-way interactions
	// population_pss_attention_initial_bias_interaction_effect_mean~ normal(0,1);
	population_pss_attention_judgement_type_interaction_effect_mean~ normal(0,1);	

	// population_logjnd_attention_initial_bias_interaction_effect_mean~ normal(0,1);
	population_logjnd_attention_judgement_type_interaction_effect_mean~ normal(0,1);	

	//population parameters
  logitRhoMean~ normal(3,3);
  
  logKappaMean~ normal(3,3);
  
  // main effects
  logitRhoAttentionEffectMean~ normal(0,3);
  // logitRhoInitialBiasEffectMean~ normal(0,3);
  // logitRhoJudgementTypeEffectMean~ normal(0,3);
  
  logKappaAttentionEffectMean~ normal(0,3);
  // logKappaInitialBiasEffectMean~ normal(0,3);
  // logKappaJudgementTypeEffectMean~ normal(0,3);
  
  // two-way interations
  // logitRhoAttentionInitialBiasEffectMean~ normal(0,3);
  // logitRhoAttentionJudgementTypeEffectMean~ normal(0,3);
  
  // logKappaAttentionInitialBiasEffectMean~ normal(0,3);
  // logKappaAttentionJudgementTypeEffectMean~ normal(0,3);
  
  cor ~ lkj_corr(4) ;  // prior on correlation matrix 
	
	#sample the betas from standard multivariate normal
	for(this_id in 1:N_toj){
	  beta[this_id] ~ multi_student_t(4,zeros,cor) ;
	}
	
	y_toj ~ bernoulli(trial_prob) ;
	//update the log-probability from p (defined in transformed parameters)
  increment_log_prob(p) ;
}
