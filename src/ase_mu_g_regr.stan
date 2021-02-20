data {
   int<lower=0> 		N; 		// total number of k-mer counts
   int				N_gk[N];	// maternal counts of k-mer
   int				y_gk[N];   	// summed counts of k-mer (mat + pat)
   int<lower=0> 		nP;   		// number of samples
   matrix[N, nP]  		indP;   	// indicator vector of samples
   int<lower=0> 		nG;   		// number of genes
   matrix[N, nG] 		indG;   	// indicator vector of genes
   int<lower=0>			nGP; 		// number of gene-sample combinations
   matrix[N, nGP] 		indGP;   	// indicator vector of genes x samples
   matrix[nGP, nP] 		map_gp_p; 
   matrix[nG, nGP] 		map_g_gp; 
   vector[nG] 			weight_g; 	// number of samples that contributed to the gene_mu
   row_vector[nP]	  	indSPO;   	// strain:PO matrix predictor level
  
   //int<lower=0>  		nDiet;   	// number of diets
   //matrix[nP, nDiet-1]	XC;		// STZ matrix indicating diet group per sample
   //matrix[nDiet, nDiet-1]	C_diet;		// STZ matrix for diet
}
 
parameters {
   real<lower=0.001, upper=0.999>		mu_r;		// mean skewing proportion for each RIX
   vector<lower=0.001, upper=0.999>[nG]		mu_g;		// mean skewing proportion for each gene
   vector<lower=0.001, upper=0.999>[nGP]	mu_gp;		// mean skewing proportion for each gene-sample
   real<lower=0> 				alpha_gp;     	// variance on mu_gp
   real<lower=0> 				alpha_g;     	// variance on mu_g
   real						beta_PORIX;     // estimated effect for PO in this RIX
   
   //matrix[nDiet-1, 1]				a_d;
   //matrix[nDiet-1, 1]				a_sdp;
   //vector[nDiet]				beta_d;
   //vector[nDiet]				beta_sdp;
}	

transformed parameters {
   real						beta_0;		// intercept for the RIX
   vector[nP]					eta_p; 
   vector<lower=0.001, upper=0.999>[nP]		mu_p;		// mean skewing proportion for each sample 
   vector<lower=0.001, upper=0.999>[nG]		weightMu_g;    
   vector<lower=0>[nG]				mu_g_a; 
   vector<lower=0>[nG]				mu_g_b;  
   vector<lower=0>[nGP]				mu_gp_a; 
   vector<lower=0>[nGP]				mu_gp_b;    
   vector[N]		 			mu_g_ind; 
  
   beta_0 = logit( mu_r );
   for (i in 1:nP ){ 
	eta_p[i]  =  beta_0 + indSPO[i] * beta_PORIX;	
   }		
   mu_p   = inv_logit( eta_p );
   weightMu_g   = (map_g_gp * mu_gp) .* weight_g; 
  
   mu_g_a  = weightMu_g 		* alpha_g  + 1;
   mu_g_b  = (1  - weightMu_g ) 	* alpha_g  + 1;
   mu_gp_a = (map_gp_p * mu_p) 		* alpha_gp + 1;
   mu_gp_b = (1 - (map_gp_p * mu_p)) 	* alpha_gp + 1; 
   mu_g_ind = indG * mu_g; 
}

model {
   beta_PORIX   ~ normal( 0, 10 );
   beta_0 	~ normal( 0, 10 );
   alpha_g	~ normal( 0, 10 ); 
   alpha_gp	~ normal( 0, 10 ); 
   mu_r  	~ beta( 1, 1 );	   	     		// flat prior on RIX-wide proportion
   mu_gp 	~ beta( mu_gp_a, mu_gp_b );       
   mu_g 	~ beta( mu_g_a , mu_g_b );
  
   y_gk ~ binomial( N_gk, mu_g_ind );	// likelihood   mu_g[mu_g_ind[i]] 
   
}

   // XC %*% a_d + ( SPO * (XC %*% a_sdp) )
   // beta_d   = aC_diet * a_d;
   // beta_sdp = C_diet * a_sdp; 

   // a_d	~ normal(0, 10);			// estimation of diet effects
   // a_sdp  	~ normal(0, 10);			// estimation of diet:strain:PO effects
   // mu_g 	~ beta( 1, 1 );    

   // mu_g_a = (map_g_gp * mu_gp)		* alpha  + 1;
   // mu_g_b = (1 - (map_g_gp * mu_gp)) 	* alpha  + 1;