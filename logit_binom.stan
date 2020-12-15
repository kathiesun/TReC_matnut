data {
   // Define variables in data
   // Number of observations (an integer)
   int<lower=0> N;
   // Number of kmers
   int<lower=0> K;
   // Number of genes
   int<lower=0> G; 
   // Number of predictors
   int<lower=0> P;
   // Variables
   int y_gk[N];
   int kmer[N];
   int n_gk[N];
   int gene[N];
}
 
parameters {
   // Define parameters to estimate
   real beta[P];
   real kmerEf[K];
   real geneEf[G];
   real<lower=0,upper=10> tau_k;
   real<lower=0,upper=10> tau_g;  
}

transformed parameters  {
   // Probability trasformation from linear predictor
   // real<lower=0> odds[N];
   real<lower=0, upper=1> prob[N];
 
   for (i in 1:N) {
   	// odds[i] = exp(beta[1] + geneEf[gene[i]] + kmerEf[kmer[i]]);
   	prob[i] = exp(beta[1] + geneEf[gene[i]] + kmerEf[kmer[i]]) /
		(exp(beta[1] + geneEf[gene[i]] + kmerEf[kmer[i]]) + 1);
     	// prob[i] = odds[i] / (odds[i] + 1);
   }
}

model {
  kmerEf ~ normal(0,tau_k);
  geneEf ~ normal(0,tau_g);
  beta 	~ normal(0, 100);

  for(i in 1:N) {
	y_gk[i] ~ binomial(n_gk[i], prob[i] );
  }
}
 
 
