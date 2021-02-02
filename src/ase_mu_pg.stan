data {
   int<lower=0> 		N; 		// total number of k-mer counts
   int<lower=0> 		N_gk[N];	// maternal counts of k-mer
   int<lower=0> 		y_gk[N];   	// summed counts of k-mer (mat + pat)
   int<lower=0> 		nP;   		// number of samples
   matrix[N, nP]  		indP;   	// indicator vector of samples
   int<lower=0> 		nG;   		// number of genes
   matrix[N, nG] 		indG;   	// indicator vector of genes
   int<lower=0>			nGP; 		// number of gene-sample combinations
   matrix[N, nGP] 		indGP;   	// indicator vector of gene-sample combinations
   matrix[nG, nP] 		map_g_p; 
   matrix[nGP, nP] 		map_gp_p; 
   matrix[nGP, nG] 		map_gp_g; 
   vector[nGP]	  		indSPO;   	// strain:PO matrix predictor level
}
 
parameters {
   vector<lower=0.001, upper=0.999>[nG]		mu_g;		// mean skewing proportion for each gene
   real<lower=0> 				alpha0;     	// variance on mu_g
   real						beta_PORIX;	// intercept for the RIX
}	

transformed parameters {
   vector[nGP] 					eta_gp; 
   vector<lower=0.001, upper=0.999>[nGP]	mu_gp;		// mean skewing proportion for each sample  
   vector[N]		 			mu_gp_ind; 
   //vector[nGP]		 			mu_po_ind; 

   eta_gp 	= logit( to_vector(map_gp_g*mu_g ) ); // + mu_po_ind; 
   mu_gp 	= inv_logit( eta_gp ); 
   mu_gp_ind 	= to_vector( indGP * mu_gp ); 				// NxnGP * nGPx1 = Nx1
   //mu_po_ind 	= indSPO * beta_PORIX;   		// nGPxN * Nx1 = nGPx1 * real = nGPx1
}

model {
   //beta_PORIX ~ normal( 0, 10 );
   alpha0	~ normal( 0, 10 );
   mu_g 	~ beta( 1, 1 );
  
   y_gk ~ binomial( N_gk, mu_gp_ind );			// likelihood   
}



