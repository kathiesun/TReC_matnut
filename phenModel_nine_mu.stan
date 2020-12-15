data {
	int<lower=1> n;		// number of data points
	real y[n];		// observations
	int nRIX;		// number of RIX groups
	int indRIX[n]; 		// indicator of RIX group
}

parameters{
	vector[nRIX] muOfRIX;			// locations of RIX groups	
	real<lower=0> sigma; 			// scale of mixtures
}		

model {
  	muOfRIX ~ normal(0, 10);
	sigma ~ gamma(0.01, 0.01);
	y ~ normal( muOfRIX[ indRIX ], sigma); 
}

