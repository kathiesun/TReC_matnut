#########
# Model #
#########


############
modelFull = "model
{
  # hyper parameters on sigma (error) and tau (rix variance)
  #tauInv2  <- pow(tau, -2)
  #tau      ~ dunif(0,1000)
  #sigInv2  <- pow(sig, -2)
  #sig      ~ dunif(0,1000)
  #tauPEInv2  <- pow(tauPE, -2)
  #tauPE      ~ dunif(0,1000)
  #tauDRInv2  <- pow(tauDR, -2)
  #tauDR      ~ dunif(0,1000)
  #tauPDRInv2  <- pow(tauPDR, -2)
  #tauPDR      ~ dunif(0,1000)
  
  # other hyper parameters on sigma (error) and tau (rix variance)
  tauInv2      ~ dgamma(0.001, 0.001)
  sigInv2      ~ dgamma(0.001, 0.001)
  tauPEInv2    ~ dgamma(0.001, 0.001)
  tauDRInv2    ~ dgamma(0.001, 0.001)
  tauPDRInv2   ~ dgamma(0.001, 0.001)
  
  for (k in 1:ndiet){
  # prior on fixed diet effect
  Diet[k] ~ dnorm(0,1/10000^2)
  }
  
  for (j in 1:nrix){
  # prior on random rix effect
  RIX[j] ~ dnorm(0, tauInv2)
  }
  
  for (j in 1:nrix){
  # prior on rix-by-poe effect
  PORIX[j] ~ dnorm(0,tauPEInv2)
  }
  
  for (a in 1:ndietrix){
  # prior on diet-by-rix effect
  DietRIX[a] ~ dnorm(0, tauDRInv2)
  }
  
  for (a in 1:ndietrix){
  # prior on po-diet-by-rix effect
  PODietRIX[a] ~ dnorm(0, tauPDRInv2)
  }
  
  for (i in 1:N){
  # Likelihood
    mu[i,1] <- RIX[rix[i]] + Diet[diet[i]] + DietRIX[dietrix[i]] + poe[i]*PORIX[rix[i]] + 
              poe[i]*PODietRIX[dietrix[i]]
    y[i]  ~ dnorm(mu[i,1], sigInv2)
  }
}"

############  
modelBreedBatch = "model
{
  # hyper parameters on sigma (error) and tau (rix variance)
  tauInv2      ~ dgamma(0.001, 0.001)
  sigInv2      ~ dgamma(0.001, 0.001)
  tauPEInv2    ~ dgamma(0.001, 0.001)
  tauDRInv2    ~ dgamma(0.001, 0.001)
  tauPDRInv2   ~ dgamma(0.001, 0.001)
  tauBBInv2    ~ dgamma(0.001, 0.001)
  tauDamInv2   ~ dgamma(0.001, 0.001)
  
  for (k in 1:ndiet){
  # prior on fixed diet effect
  Diet[k] ~ dnorm(0,1/10000^2)
  }
  
  for (m in 1:ndam){
  # prior on random dam effect
  DamID[m] ~ dnorm(0, tauDamInv2)
  }
  
  for (m in 1:nbatch){
  # prior on random batch effect
  BreedingBatch[m] ~ dnorm(0, tauBBInv2)
  }
  
  for (j in 1:nrix){
  # prior on random rix effect
  RIX[j] ~ dnorm(0, tauInv2)
  }
  
  for (j in 1:nrix){
  # prior on rix-by-poe effect
  PORIX[j] ~ dnorm(0,tauPEInv2)
  }
  
  for (a in 1:ndietrix){
  # prior on diet-by-rix effect
  DietRIX[a] ~ dnorm(0, tauDRInv2)
  }
  
  for (a in 1:ndietrix){
  # prior on po-diet-by-rix effect
  PODietRIX[a] ~ dnorm(0, tauPDRInv2)
  }
  
  for (i in 1:N){
  # Likelihood
  mu[i,1] <- BreedingBatch[batch[i]] + DamID[dam[i]] + RIX[rix[i]] + Diet[diet[i]] + 
  DietRIX[dietrix[i]] + poe[i]*PORIX[rix[i]] + poe[i]*PODietRIX[dietrix[i]] 
  y[i]  ~ dnorm(mu[i,1], sigInv2)
  }
  }"
  
  
############  
  modelDam = "model
{
  # hyper parameters on sigma (error) and tau (rix variance)
  tauInv2      ~ dgamma(0.001, 0.001)
  sigInv2      ~ dgamma(0.001, 0.001)
  tauPEInv2    ~ dgamma(0.001, 0.001)
  tauDRInv2    ~ dgamma(0.001, 0.001)
  tauPDRInv2   ~ dgamma(0.001, 0.001)
  #tauBBInv2    ~ dgamma(0.001, 0.001)
  tauDamInv2   ~ dgamma(0.001, 0.001)
  
  for (k in 1:ndiet){
  # prior on fixed diet effect
  Diet[k] ~ dnorm(0,1/10000^2)
  }

  for (m in 1:ndam){
  # prior on random rix effect
  DamID[m] ~ dnorm(0, tauDamInv2)
  }

  #for (m in 1:nbatch){
  # prior on random batch effect
  #BreedingBatch[m] ~ dnorm(0, tauBBInv2)
  #}

  for (j in 1:nrix){
  # prior on random rix effect
  RIX[j] ~ dnorm(0, tauInv2)
  }
  
  for (j in 1:nrix){
  # prior on rix-by-poe effect
  PORIX[j] ~ dnorm(0,tauPEInv2)
  }
  
  for (a in 1:ndietrix){
  # prior on diet-by-rix effect
  DietRIX[a] ~ dnorm(0, tauDRInv2)
  }
  
  for (a in 1:ndietrix){
  # prior on po-diet-by-rix effect
  PODietRIX[a] ~ dnorm(0, tauPDRInv2)
  }
  
  for (i in 1:N){
  # Likelihood
  mu[i,1] <- DamID[dam[i]] + RIX[rix[i]] + Diet[diet[i]] + DietRIX[dietrix[i]] + poe[i]*PORIX[rix[i]] + 
  poe[i]*PODietRIX[dietrix[i]] #+ BreedingBatch[batch[i]]
  y[i]  ~ dnorm(mu[i,1], sigInv2)
  }
  }"
  
  ############  
  modelCage = "model
  {
  # hyper parameters on sigma (error) and tau (rix variance)
  tauInv2      ~ dgamma(0.001, 0.001)
  sigInv2      ~ dgamma(0.001, 0.001)
  tauPEInv2    ~ dgamma(0.001, 0.001)
  tauDRInv2    ~ dgamma(0.001, 0.001)
  tauPDRInv2   ~ dgamma(0.001, 0.001)
  #tauBBInv2    ~ dgamma(0.001, 0.001)
  tauDamInv2   ~ dgamma(0.001, 0.001)
  tauCageInv2  ~ dgamma(0.001, 0.001)
  
  for (k in 1:ndiet){
  # prior on fixed diet effect
  Diet[k] ~ dnorm(0,1/10000^2)
  }
  
  for (m in 1:ndam){
  # prior on random dam effect
  DamID[m] ~ dnorm(0, tauDamInv2)
  }

  for (m in 1:ncage){
  # prior on random cage effect
  Cage[m] ~ dnorm(0, tauCageInv2)
  }
  
  #for (m in 1:nbatch){
  # prior on random batch effect
  #BreedingBatch[m] ~ dnorm(0, tauBBInv2)
  #}
  
  for (j in 1:nrix){
  # prior on random rix effect
  RIX[j] ~ dnorm(0, tauInv2)
  }
  
  for (j in 1:nrix){
  # prior on rix-by-poe effect
  PORIX[j] ~ dnorm(0,tauPEInv2)
  }
  
  for (a in 1:ndietrix){
  # prior on diet-by-rix effect
  DietRIX[a] ~ dnorm(0, tauDRInv2)
  }
  
  for (a in 1:ndietrix){
  # prior on po-diet-by-rix effect
  PODietRIX[a] ~ dnorm(0, tauPDRInv2)
  }
  
  for (i in 1:N){
  # Likelihood
  mu[i,1] <- Cage[cage[i]] + DamID[dam[i]] + RIX[rix[i]] + Diet[diet[i]] + 
  DietRIX[dietrix[i]] + poe[i]*PORIX[rix[i]] + poe[i]*PODietRIX[dietrix[i]] #+ BreedingBatch[batch[i]] 
  y[i]  ~ dnorm(mu[i,1], sigInv2) 
  }
  }"
###############  
  modelNOPODietRIX = "model
{
  # hyper parameters on sigma (error) and tau (rix variance)
  #tauInv2  <- pow(tau, -2)
  #tau      ~ dunif(0,1000)
  #sigInv2  <- pow(sig, -2)
  #sig      ~ dunif(0,1000)
  
  # other hyper parameters on sigma (error) and tau (rix variance)
  tauInv2      ~ dgamma(0.01, 0.01)
  sigInv2      ~ dgamma(0.01, 0.01)
  tauPEInv2    ~ dgamma(0.01, 0.01)
  tauDRInv2    ~ dgamma(0.01, 0.01)

  for (k in 1:ndiet){
  # prior on fixed diet effect
  Diet[k] ~ dnorm(0,1/10000^2)
  }
  
  for (j in 1:nrix){
  # prior on random rix effect
  RIX[j] ~ dnorm(0, tauInv2)
  }
  
  for (j in 1:nrix){
  # prior on rix-by-poe effect
  PORIX[j] ~ dnorm(0,tauPEInv2)
  }
  
  for (a in 1:ndietrix){
  # prior on diet-by-rix effect
  DietRIX[a] ~ dnorm(0, tauDRInv2)
  }

  for (i in 1:N){
  # Likelihood
  mu[i,1] <- RIX[rix[i]] + Diet[diet[i]] + DietRIX[dietrix[i]] + poe[i]*PORIX[rix[i]]
  y[i,1]  ~ dnorm(mu[i,1], sigInv2)
  }
  }"
  
  modelNOPO = "model
{
  # hyper parameters on sigma (error) and tau (rix variance)
  #tauInv2  <- pow(tau, -2)
  #tau      ~ dunif(0,1000)
  #sigInv2  <- pow(sig, -2)
  #sig      ~ dunif(0,1000)
  
  # other hyper parameters on sigma (error) and tau (rix variance)
  tauInv2      ~ dgamma(0.01, 0.01)
  sigInv2      ~ dgamma(0.01, 0.01)
  tauDRInv2    ~ dgamma(0.01, 0.01)
  
  for (k in 1:ndiet){
  # prior on fixed diet effect
  Diet[k] ~ dnorm(0,1/10000^2)
  }
  
  for (j in 1:nrix){
  # prior on random rix effect
  RIX[j] ~ dnorm(0, tauInv2)
  }
  
  for (a in 1:ndietrix){
  # prior on diet-by-rix effect
  DietRIX[a] ~ dnorm(0, tauDRInv2)
  }
  
  for (i in 1:N){
  # Likelihood
  mu[i,1] <- RIX[rix[i]] + Diet[diet[i]] + DietRIX[dietrix[i]]
  y[i,1]  ~ dnorm(mu[i,1], sigInv2)
  }
  }"
  
  modelNOrand = "model
{
  # hyper parameters on sigma (error) and tau (rix variance)
  #tauInv2  <- pow(tau, -2)
  #tau      ~ dunif(0,1000)
  #sigInv2  <- pow(sig, -2)
  #sig      ~ dunif(0,1000)
  
  # other hyper parameters on sigma (error) and tau (rix variance)
  tauInv2      ~ dgamma(0.01, 0.01)
  sigInv2      ~ dgamma(0.01, 0.01)

  for (k in 1:ndiet){
  # prior on fixed diet effect
  Diet[k] ~ dnorm(0,1/10000^2)
  }
  
  for (j in 1:nrix){
  # prior on random rix effect
  RIX[j] ~ dnorm(0, tauInv2)
  }

  for (i in 1:N){
  # Likelihood
  mu[i,1] <- RIX[rix[i]] + Diet[diet[i]]
  y[i,1]  ~ dnorm(mu[i,1], sigInv2)
  }
  }"