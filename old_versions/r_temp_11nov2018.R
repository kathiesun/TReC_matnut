#######IGNORE#######

"model {\n 
  sig2inv ~ dgamma(0.001,0.001)\n
  tauRIX ~ dgamma(0.001,0.001)\n
  for (i in 1:nRIX){ RIX[i] ~ dnorm(0, tauRIX) }\n
  tauDietRIX ~ dgamma(0.001,0.001)\n
  for (i in 1:nDietRIX){ DietRIX[i] ~ dnorm(0, tauDietRIX) }
  \nfor (i in 1:1){ \n 
    grand_mu ~ dunif(0, 1000) } \n 
  for (i in 1:N){ \n 
    mu[i,1] <- RIX[indRIX[i]]+DietRIX[indDietRIX[i]]+grand_mu \n 
    y[i]  ~ dnorm(mu[i,1], sig2inv) \n 
  }\n
}"




#$WeightPND60$modelFin
"model {\n 
  sig2inv ~ dgamma(0.001,0.001)\n
  tauDamID ~ dgamma(0.001,0.001)\n
  for (i in 1:nDamID){ \n 
    betaDamID[i] ~ dnorm(0, tauDamID) }\n
    tauRIX ~ dgamma(0.001,0.001)\n
  for (i in 1:nRIX){ \n 
    betaRIX[i] ~ dnorm(0, tauRIX) }\n
    tauDietRIX ~ dgamma(0.001,0.001)\n
  for (i in 1:nDietRIX){ \n 
    betaDietRIX[i] ~ dnorm(0, tauDietRIX) }\n
    tauPORIX ~ dgamma(0.001,0.001)\n
  for (i in 1:nPORIX){ \n
    betaPORIX[i] ~ dnorm(0, tauPORIX) }\n
    tauPODietRIX ~ dgamma(0.001,0.001)\n
  for (i in 1:nPODietRIX){ \n 
    betaPODietRIX[i] ~ dnorm(0, tauPODietRIX) }\n
  for (i in 1:nDiet){ betaDiet[i] ~ dnorm(0, 1e-08) }  \n
  for (i in 1:N){ \n
    mu[i,1] <- betaDamID[indDamID[i]]+betaRIX[indRIX[i]]+
    betaDietRIX[indDietRIX[i]]+PO[i]*betaPORIX[indPORIX[i]]+PO[i]*betaPODietRIX[indPODietRIX[i]]+betaDiet[indDiet[i]] \n 
  y[i]  ~ dnorm(mu[i,1], sig2inv) \n
  }\n
}"

############################


setwd("~/rna_seq/kmerSearch/")
library(tidyverse)
library(rjags)
library(coda)


dataSource <- "C:/DB Mount/Dropbox\ (ValdarLab)/outputs/matnut_outputs"
dataOrig <- readRDS(file.path(dataSource, "plotMyCounts_23oct2018.rds"))

data_genes <- dataOrig[[2]]
data_genes$fin1 <- data_genes$sum1
data_genes$fin2 <- data_genes$sum2
data_genes[which(data_genes$Xist),"fin1"] <- data_genes[which(data_genes$Xist),"sum2"]
data_genes[which(data_genes$Xist),"fin2"] <- data_genes[which(data_genes$Xist),"sum1"]

data_kmers <- dataOrig[[1]]
data_kmers$RRIX <- as.numeric(gsub("[a-z]","", data_kmers$RIX))
data_kmers$fin1 <- data_kmers$CC_1
data_kmers$fin2 <- data_kmers$CC_2
data_kmers[which(data_kmers$Xist),"fin1"] <- data_kmers[which(data_kmers$Xist),"CC_2"]
data_kmers[which(data_kmers$Xist),"fin2"] <- data_kmers[which(data_kmers$Xist),"CC_1"]

############ DON'T DO #################
data %>% select("Pup.ID") %>%
  distinct() -> merge_indP
merge_indP$indP <- seq(1:nrow(merge_indP))


data_kmers$pup.gene <- paste0(data_kmers$Pup.ID,".",data_kmers$seq.Gene)

data %>% select(one_of("seq.Gene","Pup.ID")) %>%
  distinct() %>%
  group_by(Pup.ID) %>%
  mutate(indG=seq(n())) -> merge_indG


############# encoding stuff #################
getEncoding <- function(df, terms){
  encoded <- data.frame()
  for(i in 1:length(terms)){
    
    var <- paste(terms[i])
    vec <- as.vector(t(df[,var]))
    len <- length(levels(as.factor(vec)))
    tempdf <- data.frame(Level = as.character(levels(as.factor(vec))), 
                         Index = 1:len, 
                         Variable = rep(var,len))
    encoded <- rbind(encoded,tempdf)
  }
  return(encoded)
}

##########################################

pups <- unique(data_kmers$Pup.ID)  #[which(data_genes$CCs == "CC001/CC011")])

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

jags.run <- function(data, pups, mu_g0=0.5, niter=1000, nchains=4,
                     P0=1, Q0=1, mu0=0.5, alpha0=100,
                     model=c("beta_beta.jags", "beta_gamma.jags")){
  res <- list()
  if(length(model) > 1) model = model[1]
  
  gamma = F
  if(length(grep("gamma", model)) > 0) gamma = T
  
  for(i in 1:length(pups)){
    dataPup <- data %>% filter(Pup.ID == pups[i])
    
    y_g <- dataPup$fin1   
    #as.vector(t(dataPup %>% filter(Pup.ID == pups[i]) %>% group_by(fin1) %>% dplyr::select("fin1")))
    N_g <- dataPup$sum    
    #as.vector(t(data_genes %>% filter(Pup.ID == pups[i]) %>% group_by(sum) %>% dplyr::select("sum")))
    nG <- length(y_g)
    
    #jags.data <- list("y_g", "N_g", "nG")
    jags.params <- c("mu", "alpha", "mu_g")
    if(gamma) jags.params <- c("P", "Q", "mu_g")
    jags.init <- list("mu"=mu0, "alpha"=alpha0, "mu_g"=rep(mu_g0, nG))
    if(gamma) jags.init <- list("P"=P0, "Q"=Q0, "mu_g"=rep(mu_g0, nG))
    
    
    jags <- jags.model(model,
                       data = list("y_g"=y_g, "N_g"=N_g, "nG"=nG),
                       inits = jags.init, 
                       n.chains = nchains,
                       n.adapt = niter*0.1)
    
    update(jags, niter)
    
    res[[i]] <- jags.samples(jags,
                             jags.params,
                             niter)
  }
  if(gamma){
    mu.est <- lapply(res, function(x) as.mcmc.list(x$P/(x$P + x$Q)))
  } else {
    mu.est <- lapply(res, function(x) x$mu)
  }
  if(gamma){
    alpha.est <- lapply(res, function(x) as.mcmc.list(x$P + x$Q))
  } else {
    alpha.est <- lapply(res, function(x) x$alpha)
  }
  return(list(res=res, mu.est=mu.est, alpha.est=alpha.est))
}




jags.genes.run <- function(data, mu_g0=0.5, niter=1000, nchains=4,
                     P0=1, Q0=1, mu0=0.5, alpha0=100,
                     model=c("beta_beta_reg.jags", "beta_gamma_kmers.jags")){

  if(length(model) > 1) model = model[1]
  
  gamma = F
  if(length(grep("gamma", model)) > 0) gamma = T
  
  ## for testing: testData = data_kmers %>% filter(Pup.ID %in% c(3152, 2986, 2927, 2841))

  data %>% select(one_of("Pup.ID", "seq.Gene")) %>%
    distinct() %>%
    left_join(data.frame(Pup.ID = unique(data$Pup.ID), 
                         indP = seq(1:length(unique(data$Pup.ID)))), by="Pup.ID") %>%
    ungroup() %>%
    select("indP") %>%
    t() %>% as.vector() -> indP
  
  data %>% 
    #left_join(merge_indP, by = c("Pup.ID")) %>% 
    select(one_of("seq.Gene","Pup.ID")) %>%
    distinct() %>%
    ungroup() %>%
    mutate(indG=seq(n())) -> merge_indG
  
  data %>% group_by(Pup.ID, seq.Gene) %>%
    mutate(indK=seq(n())) %>%
    group_by(Pup.ID) %>%
    #left_join(merge_indP, by = c("Pup.ID")) %>% 
    left_join(merge_indG, by = c("seq.Gene", "Pup.ID")) %>%
    select(-one_of("seq.ProbeSeq", "seq.TranscriptID", "seq.end5", "seq.end3", "seq.consensus", 
                   "seq.Chromosome", "seq.Position", "seq.refseq", "seq.altseq", "seq.alt", "seq.ref",
                   "CC_1","CC_2")) -> data
  
  y_gk <- data$fin1  
  N_gk <- data$sum   
  nP <- length(unique(data$Pup.ID))
  indG <- data$indG
  nG <- length(unique(indG))
  
  
  terms <- c("dir","RRIX", "Diet")
  encodeKmers <- getEncoding(data, terms)  
  
  
  data %>% select(one_of("Pup.ID","RRIX", "dir", "Diet")) %>%
    distinct() -> index
  
  indArray <- as.data.frame(sapply(terms, function(x){
    as.numeric(factor(as.vector(t(index[,paste(x)])), levels=encodeKmers$Level[which(encodeKmers$Variable == x)]))
  }))
 
  PO <- indArray$dir
  PO <- ifelse(PO == 1, 0.5, -0.5)
  indPORIX = indRIX = indArray$RRIX
  indDiet = indArray$Diet
  nDiet = length(unique(indDiet))
  nRIX = length(unique(indRIX))
  
  if(!identical(unique(data$Pup.ID), unique(index$Pup.ID))) warning("Check Pup.ID order")
  
  #jags.data <- list("y_g", "N_g", "nG")
  jags.params <- c("mu", "alpha", "mu_g", "betaRIX", "betaPORIX", "tauRIX", "tauPORIX")
  if(gamma) jags.params <- c("P", "Q", "mu_g")
  jags.init <- list("mu_g"=rep(mu_g0, nG),
                    #"mu"=matrix(rep(mu0, nP), nrow=1, ncol=nP), 
                    #"eta"=matrix(rep(mu0, nP), nrow=1, ncol=nP), 
                    "alpha"=rep(alpha0,nP), 
                    "betaRIX"=rep(0, nRIX), 
                    "betaPORIX"=rep(0, nRIX), 
                    "tauRIX"=alpha0, 
                    "tauPORIX"=alpha0)
  if(gamma) jags.init <- list("P"=rep(P0,nP), "Q"=rep(Q0,nP), "mu_g"=rep(mu_g0, nG))
  
  
  jags <- jags.model(model,
                     data = list("nG" = nG, "y_gk" = y_gk, "N_gk" = N_gk, "PO" = PO, 
                                 "indG" = indG, "nP"=nP, "indP"=indP, 
                                 "indRIX" = indRIX, "indPORIX" = indPORIX, "nRIX" = nRIX),
                     inits = jags.init, 
                     n.chains = nchains,
                     n.adapt = niter*0.1)
  
  update(jags, niter)
  
  res<- jags.samples(jags,
                     jags.params,
                     niter)
    
  if(gamma){
    mu.est <- as.mcmc.list(res$P/(res$P+res$Q)) #lapply(res, function(x) as.mcmc.list(x$P/(x$P + x$Q)))
  } else {
    mu.est <- res$mu #lapply(res, function(x) x$mu)
  }
  if(gamma){
    alpha.est <- as.mcmc.list(res$P+res$Q)      #lapply(res, function(x) as.mcmc.list(x$P + x$Q))
  } else {
    alpha.est <- res$alpha #lapply(res, function(x) x$alpha)
  }
  return(list(res=res, mu.est=mu.est, alpha.est=alpha.est))
}

testPups <- sample(pups, 45)

testData = data_kmers %>% filter(Pup.ID %in% testPups)

beta_reg_test <- jags.genes.run(data=testData, mu_g0=0.5, niter=50000, nchains=4,
                             P0=1, Q0=1, model="beta_beta_reg.jags")

beta_test <- jags.run(data=data_kmers, pups=pups[1:5], mu_g0=0.5, niter=2000, nchains=2,
                             P0=1, Q0=1, model="beta_beta.jags")

beta_kmers_test <- jags.genes.run(data=data_kmers, pups=pups[1:5], mu_g0=0.5, niter=2000, nchains=2, 
                            mu0=0.5, alpha0=100, model="beta_beta_kmers.jags")





#saveRDS(gamma_test, file.path(dataSource, "/all185gamma_2000_2chainz_29oct2018.rds") )

############################
plot(as.mcmc.list(beta_test$mu.est[[1]]))
plot(as.mcmc.list(beta_reg_test$mu.est))
plot(as.mcmc.list(gamma_kmers_test$mu.est[[1]]))


#########################
jags.summarize <- function(mcmc_output, pups, param=c("mu", "alpha"), nchains=4){
  est <- mcmc_output$mu.est
  if(length(param) > 1) param=param[1]
  if(param == "alpha") est <- mcmc_output$alpha.est
#browser()
  est_comb <- sapply(1:length(est), function(x) as.vector(unlist(est[[x]])), simplify=F)
  
  sep_int <- lapply(est, function(x) data.frame(HPDinterval(as.mcmc.list(x))))
  sep_int <- do.call("rbind", sep_int)
  colnames(sep_int)[1:2] <- c(paste0("lower.", nchains), paste0("upper.", nchains))
  
  comb_int <- lapply(est_comb, function(x) HPDinterval(as.mcmc(x)))
  comb_int <- do.call("rbind", comb_int)
  
  means <- lapply(est_comb, function(x) c(Mode(as.numeric(x)), mean(as.numeric(x)), median(as.numeric(x))))
  means <- do.call("rbind", means)
  colnames(means) <- c("Mode","Mean","Median")
  
  summary <- cbind(means, comb_int, sep_int)
  rownames(summary) <- pups
  return(summary)
}

beta_sum <- data.frame(data.matrix(jags.summarize(beta_test, pups)))
betaKmers_sum <- data.frame(data.matrix(jags.summarize(beta_kmers_test, pups)))
#beta_sum <- jags.summarize(beta_test, pups[1:5])

gamTemp <- gamma_sum[,c("Mean","lower","upper")]
colnames(gamTemp) <- c("X.Mean","X.Low","X.Up")
#betTemp <- beta_sum[,c("Mean","lower","upper")]
#colnames(betTemp) <- c("bet.mean","bet.low","bet.up")
#compare <- cbind(gamTemp, betTemp)
#compare$diff <- abs(compare$gam.mean - compare$bet.mean)
gamTemp$Pup.ID <- as.numeric(paste(rownames(gamTemp)))

plotTest <- data_genes %>% filter(Pup.ID %in% pups) %>%
  left_join(gamTemp, by="Pup.ID")

plot_gene_ratios <- function(plotGenes, save_path="./plots.pdf", plotOn=T){
  p <- list()
  dat <- list()
  for(i in 1:length(unique(plotGenes$CCs))){
    cc <- unique(plotGenes$CCs)[i]
    title <- paste(cc)
    
    rixdat <- plotGenes %>% filter(CCs == cc) %>%
      group_by() %>%
      mutate(minPos = as.numeric(minPos), maxPos = as.numeric(maxPos),
             meanPos = (minPos+maxPos)/2) %>%
      arrange(dir, Pup, minPos) 
    #rixdat$ratio <- as.numeric(rixdat$ratio)
    rixdat$Pup <- factor(rixdat$Pup, levels=unique(rixdat$Pup))

    p[[i]] = ggplot(data=rixdat, aes(x=meanPos)) + #, size=logSumTot
      #geom_ribbon(aes(x=minPos, ymin = gLB, ymax=gLB)) + #
      geom_point(aes(y=gRat, alpha=logSumTot, color=Xist)) + 
      geom_segment(aes(xend=meanPos, y = gLB, yend=gUB, alpha=logSumTot, color=Xist), size=1) + 
      ylab("Ratio") + xlab("Position") +
      theme_minimal() + 
      ggtitle(title) + 
      geom_ribbon(aes(ymin=X.Low, ymax=X.Up, alpha=0.1)) + 
      theme(axis.text.x = element_text(angle = 75, hjust = 1)) + 
      facet_wrap( ~ Pup, scales = "fixed") + ylim(c(0,1)) + 
      scale_colour_manual(values=c("royalblue", "sandybrown"))
    
    if("X.Mean" %in% colnames(plotGenes)){
      #rixdat[which(is.na(rixdat$X.mean)), c("X.mean", "X.Low", "X.Up")] <- 0
      
      p[[i]] = p[[i]] + geom_ribbon(aes(ymin=X.Low, ymax=X.Up), alpha=0.15) +
        geom_line(aes(y=X.Mean))
    }
    dat[[i]] <- rixdat
  }
  
  if(plotOn){
    pdf(file.path(save_path), width=8.5, height=11)
    print(p)
    dev.off()
  }
  
  return(list(p=p, dat=dat))
}



test <- plot_gene_ratios(plotTest, plotOn=T, 
  save_path=file.path(dataSource, "gamma_Xinact_plots_29oct2018.pdf"))

test$p[[1]]



#############################
y_g <- as.vector(t(data_genes %>% filter(Pup.ID == pups[1]) %>% group_by(sum1) %>% dplyr::select("sum1")))
N_g <- as.vector(t(data_genes %>% filter(Pup.ID == pups[1]) %>% group_by(sumTot) %>% dplyr::select("sumTot")))
nG <- length(y_g)

jags.data <- list("y_g", "N_g", "nG")
jags.params <- c("P", "Q", "mu_g")
jags.init <- list("P"=1, "Q"=1, "mu_g"=rep(0.5, nG))


jags <- jags.model('beta_gamma.jags',
                   data = list("y_g"=y_g, "N_g"=N_g, "nG"=nG),
                   inits = jags.init, 
                   n.chains = 4,
                   n.adapt = 100)

update(jags, 5000)

res <- jags.samples(jags,
                    jags.params,
                    5000)



beta_gamma <- "model {\n 
  P ~ dgamma(1, 0.001)\n
  Q ~ dgamma(1, 0.001)\n

  for(i in 1:nG){\n
    mu_g[i] ~ dbeta(P, Q)\n
  }\n

  for (i in 1:nG){\n
    y_g[i] ~ dbin(N_g[i], mu_g[i])\n
  }\n
}"

beta_beta <- "model {\n 
  mu ~ dbeta(1, 1)\n
  alpha ~ dgamma(1, 0.001)\n

  for(i in 1:nG){\n
    mu_g[i] ~ dbeta(mu*(1/alpha), (1-mu)*(1/alpha))\n
  }\n
  
  for (i in 1:nG){\n
    y_g[i] ~ dbeta(N_g[i], mu_g[i])\n
  }\n
}"
