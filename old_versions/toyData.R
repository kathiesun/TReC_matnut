setwd("~/matnut/src")
library(rstan)
library(tidyverse)

# ---------------------------------
# Generate toy data
# ---------------------------------
set.seed(1)

nmice=2
ngenes=3
nkmer=4

total_counts_per_kmer = pk_matrix = data = list()
pg                    = rbeta(ngenes, 1, 1)
total_counts_gene     = floor(runif(ngenes,0,1) * 10000)


for(j in 1:nmice){
  pk_matrix[[j]]              = t(sapply(1:ngenes, function(x) rbeta(nkmer,total_counts_gene[x]*pg[x] + 1,(1-pg[x])*total_counts_gene[x]+1)))
  total_counts_per_kmer       = do.call("rbind", lapply(total_counts_gene, function(x) rpois(nkmer,x/nkmer)))
  counts_per_kmer_a1          = sapply(1:length(pk_matrix[[j]]), function(x) rbinom(1, t(total_counts_per_kmer)[x], t(pk_matrix[[j]])[x]))
  temp_dat1  =  cbind(allele  = 1, 
                      count   = counts_per_kmer_a1, 
                      kmer    = seq(1:nkmer), 
                      gene    = rep(1:ngenes, each=nkmer),
                      mouse   = j)
  counts_per_kmer_a2          = sapply(1:length(pk_matrix[[j]]), function(x) rbinom(1, t(total_counts_per_kmer)[x], t(1-pk_matrix[[j]])[x]))
  temp_dat2   = cbind(allele  = 2, 
                      count   = counts_per_kmer_a2, 
                      kmer    = seq(1:nkmer), 
                      gene    = rep(1:ngenes, each=nkmer),
                      mouse   = j)
  data[[j]] = rbind(temp_dat1, temp_dat2)
}
data <- do.call("rbind", data)

as.tibble(data) %>% 
  filter(mouse == 1) %>%
  group_by(gene, kmer) %>%
  summarise(y1 = count[allele == 1],
            sum = (count[allele == 1] + count[allele == 2]), 
            ratio = count[allele == 1] / (count[allele == 1] + count[allele == 2])) %>% 
  filter(gene==1) -> 
  data_model

confInt %>% filter(gene_name == "Reps2", Pup.ID == 379) %>% 
  group_by(pos) %>% dplyr::slice(1) -> data_model



## Create Stan data
standat <-  list(N    = nrow(remNanRatios),
                 K    = length(unique(remNanRatios$pos)),
                 G    = length(unique(remNanRatios$gene_name)),
                 P    = 1,
                 y_gk = remNanRatios$count,
                 n_gk = remNanRatios$total,
                 kmer = as.numeric(as.factor(remNanRatios$pos)),
                 gene = as.numeric(as.factor(remNanRatios$gene_name)))

fileName <- "matnut/logit_binom.stan"
stan_code <- readChar(fileName, file.info(fileName)$size)
cat(stan_code)

chains=2
iter=10000
warmup=round((floor((iter*0.25)/(iter/10))*(iter/10)))
thin=50
resStan <- stan(model_code = stan_code, data = standat,
                chains = 2, iter = iter, warmup = warmup, thin = thin)
stanmcmc<- As.mcmc.list(resStan)
summcmc <- summary(stanmcmc)
traceplot(stanmcmc[,1,drop=F])
tset <- sapply(1:(ncol(stanmcmc[[1]])-1), function(x) HPDinterval(stanmcmc[,x,drop=T]), simplify=F)

sigTest <- rep(0, length = length(tset))
sigTest[grep("odds", colnames(stanmcmc[[1]]))] <- 1 
sigTest[grep("prob", colnames(stanmcmc[[1]]))] <- 0.5
        
isSig <- sapply(1:length(tset), function(x) sapply(1:chains, function(y) ifelse(((tset[[x]][[y]][1] < sigTest[x] && tset[[x]][[y]][2] < sigTest[x]) 
                                          || (tset[[x]][[y]][1] > sigTest[x] && tset[[x]][[y]][2] > sigTest[x])), T, F)))
trueSig <- intersect(which(isSig[1,] == T), which(isSig[2,] == T))

for(i in 1:10){
  ind = trueSig[i]
  traceplot(stanmcmc[,ind,drop=F])
}

######################## IGNORE ############################
##########
#  STAN  #
##########
stan_file <- file.path("../../..","Dropbox", "doubleGLM.stan")
#dissector.sm <- stan_model(file=stan_file)
tsstan
an_fit <- list()
vceg <- list()

for(j in 1:2){
  covar <- ifelse(j==1, "Sex", "Diet")
  
  for(i in 1:7){ 
    miss <- which(is.na(compl_phen[[j]][,allvars[[j]][i]]))
    mouseID_miss <- compl_phen[[j]]$MouseID[miss]
    remove <- which(as.numeric(unlist(strsplit(colnames(compl_kin[[j]]),"[.]"))[c(F,T,F)]) %in% mouseID_miss)
    if (length(mouseID_miss) > 0){
      temp_phen <- compl_phen[[j]][-which(compl_phen[[j]]$MouseID %in% mouseID_miss),]
      temp_kin <- compl_kin[[j]][-remove, -remove]
    } else {
      temp_phen <- data.frame(compl_phen[[j]])
      temp_kin <- compl_kin[[j]]
    }
    
    vceg$x <- temp_phen[,covar]
    vceg$R <- temp_kin
    vceg$Rinv <- solve(vceg$R)
    vceg$y <- as.vector(scale(as.numeric(temp_phen[,allvars[[j]][i]])))
    stan_fit <- stan(file = stan_file, control = list(adapt_delta = 0.8), 
                     data = list(num_cov = length(unique(vceg$x)),
                                 N = nrow(temp_kin),
                                 phenotype = vceg$y,
                                 cov = vceg$x,
                                 R = temp_kin), chains=3, iter = 5000, 
                     thin=10, warmup = 1000)
  }
}

#h2_mcmc <- as.mcmc(stan_fit@sim$samples[[1]]$h2)
#mcmc_trace(regex_pars = 'sigma')
rstan::traceplot(stan_fit, pars = 'h2', inc_warmup=F)
rstan::traceplot(stan_fit, pars = 'grand_sig2', inc_warmup=F)
rstan::traceplot(stan_fit, pars = 'cov_vef', inc_warmup=F)


# ---------------------------------
# Prior parameters
# ---------------------------------

n=10  
a <- 1; b <-1;
S <- 20000
tau0 <- 200
sig0 <- 1000
nu0 <- 10
mu0=120

X <- matrix(NA, nrow=n, ncol=S)
sig.1 <- numeric(S)
sig.2 <- numeric(S)
theta.1 <- numeric(S)
theta.2 <- numeric(S)
p <- numeric(S)
theta.min <- numeric(S)
theta.max <- numeric(S)


# ---------------------------------
# Random initialization to groups
# ---------------------------------

p[1] <- rbeta(1, a,b)
X[,1]<- rbinom(n, 1, p[1])
X[,1]<- ifelse(X[,1]==1, 2, 1)
theta.1[1]<- mean(Y) 
theta.2[1]<- mean(Y) 
sig.1[1]<-  var(Y)
sig.2[1]<-  var(Y)



# ---------------------------------
# Gibbs sampling algorithm
# ---------------------------------

for (i in 2:S){
  
  # Update for p
  n1 <- sum(X[,i-1]==1)
  n2 <- sum(X[,i-1]==2)
  p[i] <- rbeta(1, a + n1, b + n2)
  
  mean.y1 <- mean(Y[which(X[,i-1] == 1)])
  mean.y2 <- mean(Y[which(X[,i-1] == 2)])
  var.y1 <- var(Y[which(X[,i-1] == 1)])
  var.y2 <- var(Y[which(X[,i-1] == 2)])
  
  # Update parameters
  sig.1[i] <- 1/rgamma(1, (n1/2)+(nu0/2), 0.5*((n1-1)*var.y1 + n1*(mean.y1-theta.1[i-1])^2))
  
  sig.2[i] <- 1/rgamma(1, (n2/2)+(nu0/2), 0.5*((n2-1)*var.y2 + n1*(mean.y2-theta.2[i-1])^2))
  
  a1 <-n1*(1/sig.1[i-1])+(1/tau0)
  b1 <- n1*(1/sig.1[i-1])*mean.y1+(1/tau0)*mu0
  theta.1[i] <- rnorm(1, b1/a1, sqrt(1/a1))
  
  a2 <-n2*(1/sig.2[i-1])+(1/tau0)
  b2 <- n2*(1/sig.2[i-1])*mean.y2+(1/tau0)*mu0
  theta.2[i] <- rnorm(1, b2/a2, sqrt(1/a2))
  
  # calculate theta statistics
  theta.min[i] <- min(theta.1[i], theta.2[i])
  theta.max[i] <- max(theta.1[i], theta.2[i])
  
  
  # Update for X's 
  for (j in 1:n){
    w1 <- p[i-1]*dnorm(Y[j], theta.1[i], sqrt(sig.1[i]))
    w0 <- (1-p[i-1])*dnorm(Y[j], theta.2[i], sqrt(sig.2[i]))
    w  <- w1/(w1+w0)
    X[j,i]<- sample(1:0, 1, prob=c(w,1-w))
    X[j,i]<- ifelse(X[j,i]==0, 2, 1)
  }
}

hist(theta.1, freq = FALSE, xlim = c(90, 190), ylim = c(0, 0.17),
     xlab = expression(theta),
     main = expression(paste("Approx. of p(",theta,"|y)", sep = "")))
lines(density(theta.1), col = "blue")
lines(density(theta.2), col = "red")
hist(theta.2, freq = FALSE, 
     add = T)

# ---------------------------------
# Part c
# ---------------------------------

par(mfrow=c(1,2))
plot(acf(theta.min))
plot(acf(theta.max))

print(acf(theta.min, plot=F))
print(acf(theta.max, plot=F))


c(effectiveSize(mcmc(theta.min)),
  effectiveSize(mcmc(theta.max)))
