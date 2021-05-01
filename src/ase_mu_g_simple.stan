data {
   int<lower=0> 		N; 		// total number of k-mer counts
   int<lower=0> 		nG; 		// total number of genes
   int				N_gk[N];	// maternal counts of k-mer
   int				y_gk[N];   	// summed counts of k-mer (mat + pat)
   matrix[N, nG] 		indG;   	// indicator vector of genes
   matrix[N, nRIX] 		indRIX;   	// indicator vector of genes
   int<lower=0> 		mu_g; 		// estimated mean for gene g
}
 
parameters {
   vector<lower=0.001, upper=0.999>[nG]		nu_g;		// mean skewing proportion for each gene
   real<lower=0> 				alpha;     	// variance on mu_g
   real<lower=0>[N] 				e_i;
   mu_g = b0 + bRIX * indRIX
}
	

model {
   alpha	~ normal( 0, 10 ); 
   nu_g 	~ beta( mu_g * phi, (1-mu_g) * phi );
   e_i 		~ normal( 0, v_g );
   v_g 		~ gamma(a, b); 	
   y_gk 	~ binomial( N_gk, nu_g) ;	

}

// likelihood   mu_g[mu_g_ind[i]] 