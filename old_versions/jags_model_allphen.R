
lmerForm <- "~ -1 + (1 | DamID) + (1 | SireID) + (1 | BehaviorBatch) + OFBox + Diet + (1+PO|RIX) + (1+PO|Diet:RIX)" 
lmerForm <- "~ -1 + (1 | DamID) + (1 | SireID) + (1 | BehaviorBatch) + SIHOrder + Diet + (1+PO|RIX) + (1+PO|Diet:RIX)" 
lmerForm <- "~ -1 + (1 | DamID) + (1 | SireID) + (1 | BehaviorBatch) + FSTChamber + Diet + (1+PO|RIX) + (1+PO|Diet:RIX)" 
lmerForm <- "~ -1 + (1 | DamID) + (1 | SireID) + (1 | BehaviorBatch) + RestraintExperimenter + RestraintOrder + Diet + (1+PO|RIX) + (1+PO|Diet:RIX)" 
lmerForm <- "~ -1 + (1 | DamID) + (1 | SireID) + (1 | BehaviorBatch) + Diet + (1+PO|RIX) + (1+PO|Diet:RIX)" 


##################
modelMatch = "model
{
  # hyper parameters on sigma (error) and tau (rix variance)
  tauInv2      ~ dgamma(0.001, 0.001)     #residual error
  sigInv2      ~ dgamma(0.001, 0.001)     #random RIX effect
  tauDRInv2    ~ dgamma(0.001, 0.001)     #Diet-RIX effect
  
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
    mu[i,1] <- RIX[rix[i]] + DietRIX[dietrix[i]] 
    y[i]  ~ dnorm(mu[i,1], sigInv2) 
  }
}"

#################


modelWeight = "model
  {
  # hyper parameters on sigma (error) and tau (rix variance)
  tauInv2      ~ dgamma(0.001, 0.001)     #residual error
  sigInv2      ~ dgamma(0.001, 0.001)     #random RIX effect
  
  tauPEInv2    ~ dgamma(0.001, 0.001)     #POE effect
  tauDRInv2    ~ dgamma(0.001, 0.001)     #Diet-RIX effect
  tauPDRInv2   ~ dgamma(0.001, 0.001)     #PO-diet-RIX effect
  tauBEInv2    ~ dgamma(0.001, 0.001)     #Behav batch effect
  tauDamInv2   ~ dgamma(0.001, 0.001)     #Dam effect
  tauSirInv2   ~ dgamma(0.001, 0.001)     #Sire effect

  for (k in 1:ndiet){
  # prior on fixed diet effect
  Diet[k] ~ dnorm(0,1/10000^2)
  }
  
  for (m in 1:ndam){
  # prior on random dam effect
  DamID[m] ~ dnorm(0, tauDamInv2)
  }
  
  for (m in 1:nsir){
  # prior on random sire effect
  SireID[m] ~ dnorm(0, tauSirInv2)
  }
  
  for (m in 1:nbatch){
  # prior on random behav batch effect
  BehaviorBatch[m] ~ dnorm(0, tauBEInv2)
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
  mu[i,1] <- SireID[sire[i]] + DamID[dam[i]] + BehaviorBatch[be[i]] + RIX[rix[i]] + Diet[diet[i]] + DietRIX[dietrix[i]] + 
  poe[i]*PORIX[rix[i]] + poe[i]*PODietRIX[dietrix[i]] 
  
  y[i]  ~ dnorm(mu[i,1], sigInv2) 
  }
}"


##########

modelOF = "model
  {
# hyper parameters on sigma (error) and tau (rix variance)
  tauInv2      ~ dgamma(0.001, 0.001)     #residual error
  sigInv2      ~ dgamma(0.001, 0.001)     #random RIX effect
  
  tauPEInv2    ~ dgamma(0.001, 0.001)     #POE effect
  tauDRInv2    ~ dgamma(0.001, 0.001)     #Diet-RIX effect
  tauPDRInv2   ~ dgamma(0.001, 0.001)     #PO-diet-RIX effect
  tauBEInv2    ~ dgamma(0.001, 0.001)     #Behav batch effect
  tauDamInv2   ~ dgamma(0.001, 0.001)     #Dam effect
  tauSirInv2   ~ dgamma(0.001, 0.001)     #Sire effect

  OFBox[1] = 0
  
  for (k in 1:ndiet){
  # prior on fixed diet effect
  Diet[k] ~ dnorm(0,1/10000^2)
  }

  for (m in 2:nbox){
  # prior on fixed test chamber effect
  OFBox[m] ~ dnorm(0, 1/10000^2)
  }

  for (m in 1:ndam){
  # prior on random dam effect
  DamID[m] ~ dnorm(0, tauDamInv2)
  }
  
  for (m in 1:nsir){
  # prior on random sire effect
  SireID[m] ~ dnorm(0, tauSirInv2)
  }
  
  for (m in 1:nbatch){
  # prior on random behav batch effect
  BehaviorBatch[m] ~ dnorm(0, tauBEInv2)
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
  mu[i,1] <- SireID[sire[i]] + DamID[dam[i]] + BehaviorBatch[be[i]] + OFBox[of[i]] + RIX[rix[i]] + Diet[diet[i]] + 
  DietRIX[dietrix[i]] + poe[i]*PORIX[rix[i]] + poe[i]*PODietRIX[dietrix[i]] 
  y[i]  ~ dnorm(mu[i,1], sigInv2) 
  }
}"

##########

modelLD = "model
{
  # hyper parameters on sigma (error) and tau (rix variance)
  tauInv2      ~ dgamma(0.001, 0.001)     #residual error
  sigInv2      ~ dgamma(0.001, 0.001)     #random RIX effect
  
  tauPEInv2    ~ dgamma(0.001, 0.001)     #POE effect
  tauDRInv2    ~ dgamma(0.001, 0.001)     #Diet-RIX effect
  tauPDRInv2   ~ dgamma(0.001, 0.001)     #PO-diet-RIX effect
  tauBEInv2    ~ dgamma(0.001, 0.001)     #Behav batch effect
  tauDamInv2   ~ dgamma(0.001, 0.001)     #Dam effect
  tauSirInv2   ~ dgamma(0.001, 0.001)     #Sire effect

  LDChamber[1] = 0
  
  for (k in 1:ndiet){
  # prior on fixed diet effect
  Diet[k] ~ dnorm(0,1/10000^2)
  }

  for (m in 2:nbox){
  # prior on fixed chamber effect
  LDChamber[m] ~ dnorm(0,1/10000^2)
  }

  for (m in 1:ndam){
  # prior on random dam effect
  DamID[m] ~ dnorm(0, tauDamInv2)
  }
  
  for (m in 1:nsir){
  # prior on random sire effect
  SireID[m] ~ dnorm(0, tauSirInv2)
  }
  
  for (m in 1:nbatch){
  # prior on random test chamber
  BehaviorBatch[m] ~ dnorm(0, tauBEInv2)
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
  mu[i,1] <- SireID[sire[i]] + DamID[dam[i]] + BehaviorBatch[be[i]] + LDChamber[ld[i]] + RIX[rix[i]] + Diet[diet[i]] + 
  DietRIX[dietrix[i]] + poe[i]*PORIX[rix[i]] + poe[i]*PODietRIX[dietrix[i]] 
  y[i]  ~ dnorm(mu[i,1], sigInv2) 
  }
}"
  
##########
  
modelSIH = "model
  {
  # hyper parameters on sigma (error) and tau (rix variance)
  tauInv2      ~ dgamma(0.001, 0.001)     #residual error
  sigInv2      ~ dgamma(0.001, 0.001)     #random RIX effect
  
  tauPEInv2    ~ dgamma(0.001, 0.001)     #POE effect
  tauDRInv2    ~ dgamma(0.001, 0.001)     #Diet-RIX effect
  tauPDRInv2   ~ dgamma(0.001, 0.001)     #PO-diet-RIX effect
  tauBEInv2    ~ dgamma(0.001, 0.001)     #Behav batch effect
  tauDamInv2   ~ dgamma(0.001, 0.001)     #Dam effect
  tauSirInv2   ~ dgamma(0.001, 0.001)     #Sire effect

  SIHOrder[1] = 0
  
  for (k in 1:ndiet){
  # prior on fixed diet effect
  Diet[k] ~ dnorm(0,1/10000^2)
  }

  for (m in 2:norder){
  # prior on fixed test order
  SIHOrder[m] ~ dnorm(0, 1/10000^2)
  }

  for (m in 1:ndam){
  # prior on random dam effect
  DamID[m] ~ dnorm(0, tauDamInv2)
  }
  
  for (m in 1:nsir){
  # prior on random sire effect
  SireID[m] ~ dnorm(0, tauSirInv2)
  }
  
  for (m in 1:nbatch){
  # prior on random behav batch effect
  BehaviorBatch[m] ~ dnorm(0, tauBEInv2)
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
  mu[i,1] <- SireID[sire[i]] + DamID[dam[i]] + BehaviorBatch[be[i]] + SIHOrder[sih[i]] + RIX[rix[i]] + Diet[diet[i]] + 
  DietRIX[dietrix[i]] + poe[i]*PORIX[rix[i]] + poe[i]*PODietRIX[dietrix[i]] 
  y[i]  ~ dnorm(mu[i,1], sigInv2) 
  }
}"  

##########

modelFST = "model
{
  # hyper parameters on sigma (error) and tau (rix variance)
  tauInv2      ~ dgamma(0.001, 0.001)     #residual error
  sigInv2      ~ dgamma(0.001, 0.001)     #random RIX effect
  
  tauPEInv2    ~ dgamma(0.001, 0.001)     #POE effect
  tauDRInv2    ~ dgamma(0.001, 0.001)     #Diet-RIX effect
  tauPDRInv2   ~ dgamma(0.001, 0.001)     #PO-diet-RIX effect
  tauBEInv2    ~ dgamma(0.001, 0.001)     #Behav batch effect
  tauDamInv2   ~ dgamma(0.001, 0.001)     #Dam effect
  tauSirInv2   ~ dgamma(0.001, 0.001)     #Sire effect

  FSTChamber[1] = 0
  
  for (k in 1:ndiet){
  # prior on fixed diet effect
  Diet[k] ~ dnorm(0,1/10000^2)
  }

  for (m in 2:nbox){
  # prior on fixed test chamber effect
  FSTChamber[m] ~ dnorm(0, 1/10000^2)
  }
  
  for (m in 1:ndam){
  # prior on random dam effect
  DamID[m] ~ dnorm(0, tauDamInv2)
  }
  
  for (m in 1:nsir){
  # prior on random sire effect
  SireID[m] ~ dnorm(0, tauSirInv2)
  }
  
  for (m in 1:nbatch){
  # prior on random behav batch effect
  BehaviorBatch[m] ~ dnorm(0, tauBEInv2)
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
  mu[i,1] <- SireID[sire[i]] + DamID[dam[i]] + BehaviorBatch[be[i]] + FSTChamber[fst[i]] + RIX[rix[i]] + Diet[diet[i]] + 
  DietRIX[dietrix[i]] + poe[i]*PORIX[rix[i]] + poe[i]*PODietRIX[dietrix[i]] 
  y[i]  ~ dnorm(mu[i,1], sigInv2) 
  }
}"
  
##########
  
modelStress = "model
  {
  # hyper parameters on sigma (error) and tau (rix variance)
  tauInv2      ~ dgamma(0.001, 0.001)     #residual error
  sigInv2      ~ dgamma(0.001, 0.001)     #random RIX effect
  
  tauPEInv2    ~ dgamma(0.001, 0.001)     #POE effect
  tauDRInv2    ~ dgamma(0.001, 0.001)     #Diet-RIX effect
  tauPDRInv2   ~ dgamma(0.001, 0.001)     #PO-diet-RIX effect
  tauBEInv2    ~ dgamma(0.001, 0.001)     #Behav batch effect
  tauDamInv2   ~ dgamma(0.001, 0.001)     #Dam effect
  tauSirInv2   ~ dgamma(0.001, 0.001)     #Sire effect

  RestraintOrder[1]          = 0
  RestraintExperimenter[1]   = 0

  for (k in 1:ndiet){
  # prior on fixed diet effect
  Diet[k] ~ dnorm(0, 1/10000^2)
  }

  for (m in 2:norder){
  # prior on fixed test order effect
  RestraintOrder[m] ~ dnorm(0, 1/10000^2)      
  }

  for (m in 2:nex){
  # prior on fixed test experimenter effect
  RestraintExperimenter[m] ~ dnorm(0, 1/10000^2) 
  }
  
  for (m in 1:ndam){
  # prior on random dam effect
  DamID[m] ~ dnorm(0, tauDamInv2)
  }
  
  for (m in 1:nsir){
  # prior on random sire effect
  SireID[m] ~ dnorm(0, tauSirInv2)
  }
  
  for (m in 1:nbatch){
  # prior on random behav batch effect
  BehaviorBatch[m] ~ dnorm(0, tauBEInv2)
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
  mu[i,1] <- SireID[sire[i]] + DamID[dam[i]] + BehaviorBatch[be[i]] + RestraintOrder[restord[i]] + 
  RestraintExperimenter[restex[i]] + RIX[rix[i]] + Diet[diet[i]] + DietRIX[dietrix[i]] + 
  poe[i]*PORIX[rix[i]] + poe[i]*PODietRIX[dietrix[i]] 
  y[i]  ~ dnorm(mu[i,1], sigInv2) 
  }
}"