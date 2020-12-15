data {
	int<lower=1> n;		// number of data points
	real y[n];		// observations
	int nRIX;		// number of RIX groups
	int indRIX[n]; 		// indicator of RIX group
	int nK; 
}

parameters{
	real mu_a; 
	real<lower=0> sigma; 			// scale of mixtures
	simplex[nK] p[nRIX];			// mixing proportion
}		


transformed parameters { 	
	//vector[nK] log_p;
	//for(i in 1:K)
	//	log_p[i] = log( p[i] ); 
	
	//for(i in 1:n){
	//	vector[nK] lps = log(p[ indRIX[i] ]);
	//	lps[ 1 ] += normal_lpdf( y[i] | -mu_a,  sigma);
	//	lps[ 2 ] += normal_lpdf( y[i] | 0,  sigma);
	//	lps[ 3 ] += normal_lpdf( y[i] | mu_a,  sigma);
	//	target += log_sum_exp(lps);
	//}
	// muOfRIX = rep_vector( mu_a, nRIX ); 
	// vector[nRIX] muOfRIX;			// locations of RIX groups
	real<lower=0.0001> mu_a_sq = mu_a * mu_a;		// squared mu_a
	real<lower=0.0001> mu_abs  = sqrt( mu_a_sq );	
}

model {
  	mu_a ~ normal(0, 10);	
	for(r in 1:nRIX)
		p[r,] ~ dirichlet(rep_vector(1, nK));
	sigma ~ gamma(0.01, 0.01);
	
	//for(k in 1:nK){
	//	target += log_sum_exp(lp[k]);
	//}

	for(i in 1:n){
		vector[nK] lps = log(p [ indRIX[i], ]);
		lps[ 1 ] += normal_lpdf( y[i] | -mu_abs,  sigma);
		lps[ 2 ] += normal_lpdf( y[i] | 0,  	sigma);
		lps[ 3 ] += normal_lpdf( y[i] | mu_abs, sigma);
		target += log_sum_exp(lps);
	}
}


###########################

data {
	int<lower=1> n;		// number of data points
	real y[n];		// observations
	int nRIX;		// number of RIX groups
	int indRIX[n]; 		// indicator of RIX group
}

parameters{
	//vector[nRIX] muOfRIX;			// locations of RIX groups
	// vector[nRIX] muOfRIXsq;		// squared mu's RIX-length vector
	real<lower=0> sigma; 			// scale of mixtures
	simplex[2] p;				// mixing proportion
	real mu_a;             		// locations of mixture component
}		


transformed parameters {
	p[2] = p[3]; 	
	vector[3] lp[nRIX]; 
	lp = log(p);	

	//real mu_a_sq = pow( mu_a, 2 );		// squared mu_a
	for(i in 1:n){
		lp[1] += normal_lpdf( y[i] | 0 , sigma);
		lp[2] += normal_lpdf( y[i] | mu_a , sigma);
	}
}

model {
  	mu_a ~ normal(0, 10);
	p ~ beta(1,1);
	sigma ~ gamma(0.01, 0.01);
	target += log_sum_exp(lp);
}

		//lp[i] = lp[i] + normal_lpdf(y[i] | muOfRIXsq[ indRIX[i] ] = mu_a ? muOfRIX[ indRIX[i] ] : 0, sigma);
		//lp[1] += normal_lpdf(y[i] | muOfRIX[ indRIX[i] ], sigma);		// for mu_i=mu_a
		//lp[2] += normal_lpdf(y[i] | 0, sigma);				// for mu_i=0

  	//mu_a ~ normal(0, 10);
	//p ~ beta(1,1);


	//sigma ~ gamma(0.01, 0.01);