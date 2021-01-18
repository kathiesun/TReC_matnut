data {
   int<lower=0> 		N; 		// total number of k-mer counts
   int<lower=0> 		nG; 		// total number of genes
   int				N_gk[N];	// maternal counts of k-mer
   int				y_gk[N];   	// summed counts of k-mer (mat + pat)
   matrix[N, nG] 		indG;   	// indicator vector of genes
}
 
parameters {
   vector<lower=0.001, upper=0.999>[nG]		mu_g;		// mean skewing proportion for each gene
   real<lower=0> 				alpha;     	// variance on mu_g
}	

model {
   alpha	~ normal( 0, 10 ); 
   mu_g 	~ beta( 1, 1 );
  
   y_gk ~ binomial( N_gk, indG*mu_g );	
}

// likelihood   mu_g[mu_g_ind[i]] 