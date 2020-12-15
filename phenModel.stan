data {
	int<lower=1> n;		// number of data points
	real y[n];		// observations
	int nRIX;		// number of RIX groups
	int indRIX[n]; 		// indicator of RIX group
}

parameters{
	row_vector[nRIX] muOfRIX;		// locations of RIX groups
	simplex[2] p; 				// mixing proportion
	real mu_a;				// location of mixtures
	real<lower=0> sigma; 			// scale of mixtures
}		


transformed parameters {
	real mu_sq = pow( mu_a, 2 );		// location of mixtures
	real muOfRIXsq[nRIX];
	//vector[n] log_lik;
	vector[2] lp;			

	for(r in 1:nRIX)
		muOfRIXsq[r] = pow( muOfRIX[r], 2 ); 

  	for (i in 1:n){
		//log_lik[i] = log_lik[i] + normal_lpdf(y[i] | muOfRIXsq[ indRIX[i] ] = mu_a ? muOfRIX[ indRIX[i] ] : 0, sigma);
		lp = log(p); 
		if(muOfRIXsq[ indRIX[i] ] == mu_sq){
			lp[1] += normal_lpdf(y[i] | muOfRIX[ indRIX[i] ], sigma);
			//lp[3] += normal_lpdf(y[i] | -mu_a, sigma);
		}
		else {
			lp[2] += normal_lpdf(y[i] | 0, sigma);
		}

	}
}

model {
  	//target += log_sum_exp(lp);
	//target += uniform_lpdf(p | 0, 1);
 	//target += normal_lpdf(mu_a | 0, 10);
	//target += gamma_lpdf(sigma | 0.01, 0.01);

	mu_a ~ normal(0, 10);
	p ~ beta(1,1);
	sigma ~ gamma(0.01, 0.01);
	target += log_sum_exp(lp);
}



