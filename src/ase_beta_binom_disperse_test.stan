/*
*Simple beta-binomial example
*/

data {
  int N; 		//the number of observations
  int y[N]; 		//the response
  int W[N]; 		//the number of trials per observations
  matrix[N,1] X; 	//the model matrix
}

parameters {
  vector[1] betas;	//the regression parameters
  real phi; 		//the overdispersion parameter
}
transformed parameters {
  vector[N] mu; 	//the linear predictor
  vector[N] A; 		//the first shape parameter for the beta distribution
  vector[N] B; 		//the second shape parameter for the beta distribution
  
  for(n in 1:N)
    mu[n] = inv_logit(X[n,]*betas); 	//using logit link

  A = mu * phi;
  B = (1-mu) * phi;
}
model {  
  betas ~ cauchy(0,10); 		//prior for the intercept following Gelman 2008
  y 	~ beta_binomial(W,A,B);
}

generated quantities {
  vector[N] y_rep;
  for(n in 1:N){
    y_rep[n] = beta_binomial_rng(W[n],A[n],B[n]); //posterior draws to get posterior predictive checks
  }
}