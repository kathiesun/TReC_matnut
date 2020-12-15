data {
   int<lower=0> 	N;		// Number of observations
   vector[N] 		y;      	// outcome vector

   int<lower=0> 	K_d;   		// number of diet predictor levels
   matrix[N, K_d] 	x_d;   		// diet contrast matrix
   int<lower=0> 	K_s;   		// number of strain predictors
   matrix[N, K_s] 	x_s;   		// strain contrast matrix
   int<lower=0> 	K_sd;   	// number of strain:diet predictor levels
   matrix[N, K_sd] 	x_sd;   	// strain:diet contrast matrix

   vector[N] 		SPO; 		// strain x PO matrix


}
 
parameters {
   real 		beta_0;			// intercept
   vector[K_s]		beta_SPO_raw; 		// coefficients for SPO (unique for each rix)
  
   vector[K_d] 		beta_d_raw;     	// coefficients for diet effect predictors
   vector[K_s] 		beta_s_raw;     	// coefficients for strain effect predictors
   vector[K_sd] 	beta_sd_raw;    	// coefficients for strain:diet effect predictors
   vector[K_sd] 	beta_sdp_raw;   	// coefficients for PO:strain:diet effect predictors

   real<lower=0> 	sigma;  		// error on values
   real<lower=0> 	lambda;  		// scale for double exponential (laplace) prior on beta_SPO
   // real<lower=0, upper=1> lambda;  		// indicator variable spike-or-slab
   // real<lower=0, upper=1> pi;  		// prior proportion on spike-or-slab
   // real 		c;  			// scale for non-zero beta_SPO
   real<lower=0> 	sigma_beta_s;  		// error on beta values
   real<lower=0> 	sigma_beta_sd;  	// error on beta values
   real<lower=0> 	sigma_beta_sdp;		// error on beta values
}	

transformed parameters {
	// sum-to-zero constraints
   vector[K_s+1]  beta_SPO = append_row(beta_SPO_raw, -sum(beta_SPO_raw));
   vector[K_d+1]  beta_d   = append_row(beta_d_raw, -sum(beta_d_raw));
   vector[K_s+1]  beta_s   = append_row(beta_s_raw, -sum(beta_s_raw));
   vector[K_sd+1] beta_sd  = append_row(beta_sd_raw, -sum(beta_sd_raw));
   vector[K_sd+1] beta_sdp = append_row(beta_sdp_raw, -sum(beta_sdp_raw));
}

model {
   beta_0  		~ normal(0, 10);          		// estimation of intercept
   beta_SPO_raw  	~ double_exponential(0, lambda);	// estimation of SPO betas
   // beta_SPO_raw  	~ normal(0, pow(c,2) * lambda);	// estimation of SPO betas
   beta_d_raw    	~ normal(0, 10); 			// estimation of fixed diet effects
   beta_s_raw  		~ normal(0, sigma_beta_s);  		// estimation of strain effect betas
   beta_sd_raw  	~ normal(0, sigma_beta_sd);  		// estimation of strain:diet effect betas
   beta_sdp_raw  	~ normal(0, sigma_beta_sdp);  		// estimation of PO:strain:diet effect betas
   
   sigma		~ normal(0, 10);	 		// prior on error
   lambda		~ gamma(2, 0.1); 			// prior on SPO lambda
   // lambda 		~ beta(1,1);				// prior on spike-or-slab pi
   // c 		~ normal(0, 10); 			// prior on non-zero SPO beta
   sigma_beta_s		~ normal(0, 10); 			// prior on strain random effect betas
   sigma_beta_sd	~ normal(0, 10); 			// prior on strain:diet random effect betas
   sigma_beta_sdp	~ normal(0, 10); 			// prior on PO:strain:diet random effect betas
  
   y ~ normal(beta_0 + (SPO .* (x_s * beta_SPO_raw)) + x_d * beta_d_raw + x_s * beta_s_raw + x_sd * beta_sd_raw + (SPO .* (x_sd * beta_sdp_raw)) , sigma);	// likelihood
}
 
 
