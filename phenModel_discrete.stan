data {
	int<lower=1> n;			// number of data points
	real y[n];			// observations
	int nRIX;			// number of RIX groups
	int indRIX[n]; 			// indicator of RIX group
	int nK; 
}

parameters{
	real mu_a; 
	real<lower=0> sigma; 		// scale of mixtures
	simplex[nK] p[nRIX];		// mixing proportion
}		


transformed parameters { 	
	real<lower=0.0001> mu_a_sq = mu_a * mu_a;		// squared mu_a
	real<lower=0.0001> mu_abs  = sqrt( mu_a_sq );		// absolute val mu_a
}

model {
  	mu_a ~ normal(0, 10);	
	for(r in 1:nRIX)
		p[r,] ~ dirichlet(rep_vector(1, nK));
	sigma ~ gamma(0.01, 0.01);
	
	for(i in 1:n){
		vector[nK] lps = log(p [ indRIX[i], ]);
		lps[ 1 ] += normal_lpdf( y[i] | -mu_abs,  sigma);
		lps[ 2 ] += normal_lpdf( y[i] | 0,  	sigma);
		lps[ 3 ] += normal_lpdf( y[i] | mu_abs, sigma);
		target += log_sum_exp(lps);
	}
}