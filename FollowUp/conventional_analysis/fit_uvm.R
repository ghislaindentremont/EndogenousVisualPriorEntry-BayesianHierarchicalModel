#Load the CicStats package
library(CircStats)

#Define a customization of circ.mean{CircStats} that permits weighting observations
circ.weighted.mean = function (x,rho){
	sinr = sum(rho*sin(x))
	cosr = sum(rho*cos(x))
	circmean = atan2(sinr, cosr)
	return(circmean)
}

#Define a Customization of est.kappa{CircStats} that permits assuming a mean of zero
est_kappa = function (x, rho, mu ){
	kappa = A1inv(sum(rho*cos(x-mu))/sum(rho))
	return(kappa)
}

#Define a function to perform Expectation-Maximization
#   estimation of a uniform+VonMises mixture model
em_uvm = function( x , rho_start , do_mu , max_steps = 1e4 , max_reset = 1e3 , rel_tol=1.e-03 , trace = TRUE){
	rho = rho_start
	if(do_mu){
		mu = circ.weighted.mean(x,rho)
	}else{
		mu = 0
	}
	kappa = est_kappa(x,rho,mu)
	eps = Inf
	last_NSLL = Inf
	steps = 0
	reset = 0
	min_eps = eps
	unif = dvm(0,0,0)
	suppressWarnings(vm <- dvm(x,mu,kappa))
	while (eps > rel_tol) {
		if(steps>max_steps){
			rel_tol = min_eps
			kappa = kappa_start
			rho = rho_start
			eps = 1.0
			steps = 0
		}else{
			rho = rho*vm/(rho*vm+(1-rho)*unif)
			if(any(!is.finite(rho))){
				if(reset<max_reset){
					#try to reset the fit by re-starting at a random location
					kappa = runif(1,0,exp(8))
					mu = runif(1,0,2*pi)
					rho = runif(1,0,1)
					suppressWarnings(vm <- dvm(x,mu,kappa))
					reset = reset+1
					steps = 0
				}else{
					eps = -1
					kappa = NA
					rho =  NA
				}
			}else{
				if(all(rho==0)){
					rho = 0
					kappa = NA
					eps = -1
				}else{
					if(do_mu){
						mu = circ.weighted.mean(x,rho)
					}else{
						mu = 0
					}
					kappa = est_kappa(x,rho,mu)					
					rho = mean(rho)
					if(kappa<=0){
						rho = 0
						kappa = NA
						eps = -1
					}else{
						steps = steps + 1
						suppressWarnings(vm <- dvm(x,mu,kappa))
						this_NSLL = -sum(log(rho*vm+(1-rho)*unif))
						eps = last_NSLL-this_NSLL
						min_eps = ifelse(eps<min_eps,eps,min_eps)
						last_NSLL = this_NSLL
					}
				}
			}
		}
		if(trace){
			cat("trace:",steps,reset,rho,kappa,eps,"\n")
		}
	}
	return(
		list(
			rho = rho
			, kappa_prime = log(kappa)
			, mu = ifelse(do_mu,mu,NA)
			, rel_tol = rel_tol
			, steps = steps
			, reset = reset
		)
	)
}


#define some functions to help obtain an r-squared measure of fit
get_mix_density = function(x,rho,mu,kappa_prime){
	suppressWarnings(dmixedvm(x,mu,0,exp(kappa_prime),0,rho))
}
get_mix_cdf = function(x,rho,mu,kappa_prime){
	value = integrate(
		get_mix_density
		, lower = 0
		, upper = x
		, rho = rho
		, mu = mu
		, kappa_prime = kappa_prime
	)$value
	return(value)
}
get_rsq = function(x,rho,mu,kappa_prime){
	x = sort(x)
	n = length(x)
	obs_p = (1:n)/n
	exp_p = rep(NA,n)
	for(i in 1:n){
		exp_p[i]=get_mix_cdf(x[i],rho,mu,kappa_prime)
	}
	rsq = cor(obs_p,exp_p)^2
	return(rsq)
}

#Define a function to run em_uvm2() across a range of 
#   starting rho values 
fit_uvm = function(data,do_mu){
	rho_start = seq(.1,1,.1)
	b = expand.grid(
		rho_start = rho_start
		, rho = NA
		, kappa_prime = NA
		, mu = NA
		, rel_tol = NA
		, steps = NA
		, reset = NA
		, SLL = NA
	)
	for(i in 1:nrow(b)){
		fit=NA
		try(fit<-em_uvm(
				data
				, rho_start = b$rho_start[i]
				, do_mu = do_mu
				, trace = FALSE
			)
			, silent = TRUE
		)
		if(!is.list(fit)){
			fit=list(
				rho = NA
				, kappa_prime = NA
				, mu = NA
				, rel_tol = NA
				, steps = NA
				, reset = NA
			)
		}
		b[i,2:7] = as.numeric(fit)
		if(!is.na(fit$kappa_prime)){
			suppressWarnings(b$SLL[i] <- sum(log(dmixedvm(data,0,0,exp(fit$kappa_prime),0,fit$rho))))
		}else{
			if(!is.na(fit$rho)){
				b$SLL[i] = sum(log(dunif(data,0,2*pi)))
			}
		}
	}
	if(all(is.na(b$SLL))){
		fit=list(
			rho = NA
			, kappa_prime = NA
			, mu = NA
			, rel_tol = NA
			, steps = NA
			, reset = NA
			, SLL = NA
			, rsq = NA
		)
	}else{
		fit = b[which.max(b$SLL),2:ncol(b)]
		if(do_mu){
			if(!is.na(fit$kappa_prime)){
				fit$rsq = get_rsq(data,fit$rho,fit$mu,fit$kappa_prime)
			}else{
				fit$rsq = get_rsq(data,fit$rho,fit$mu,0)
			}
		}else{
			if(!is.na(fit$kappa_prime)){
				fit$rsq = get_rsq(data,fit$rho,0,fit$kappa_prime)
			}else{
				fit$rsq = get_rsq(data,fit$rho,0,0)
			}
		}
	}
	return(fit)
}

