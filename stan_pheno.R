library(tidyverse)
library(rstan)

setwd("C:/Users/Kathie/models_matnut/src")

source("lmer_functions_rna.R")
source("jags_functions.R")
source("prediction_functions.R")
source("summary_functions.R")


###### Read in data
dir <- "C:/Users/Kathie/Dropbox\ (ValdarLab)"
matnut <- readRDS(file.path(dir,'phenotype_analysis/matnut_data.rds'))

stanlist <- lapply(matnut$ptypes[5:20], function(x) stanSum(df=matnut$df, encoded=matnut$encoded, phenotype=x, 
                   randvar=c("DamID", "RIX", "DietRIX"), fixvar="Diet", POvar=c("RIX", "Diet:RIX"),
                   tryLam=c(-1, 0, .25, .33, .5, 1, 2, 3), normd=T, 
                   chains=1, iter=2000))

stanSum <- function(df, phenotype, encoded=NULL, 
                    tryLam=1, normd=T,
                    chains=2, iter=10000,
                    plot=T, contrasts=F,
                    randvar=NA, fixvar=NA, POvar=NA){
  
  formulas = getFormulas(fixvar, randvar, POvar)
  
  if(length(phenotype) > 1) {
    use <- which(colSums(df[,colnames(df) %in% phenotype], na.rm = T) != 0)
    use <- unique(c(use, which(apply(df[,colnames(df) %in% phenotype], 2, var) != 0) ))
    phenotype <- phenotype[use]
  }
  y.mat <- data.frame(df[, phenotype])
  df <- data.frame(df)
  colnames(y.mat) <- phenotype
  
  bcObject <- BC.model(y.mat = y.mat, data=df, indvariable=formulas$lmerform, 
                       transformParams = getMatnutTransformParams(tryLam = tryLam, normd = normd))
  transf <- data.frame(lambda = bcObject$lambda,  
                       pval = bcObject$p.val)
  
  form <- formulaWrapper$parseCovariateString(paste(formulas$jagsform))
  
  if(any(duplicated(form$fixef))){
    form$fixef[[which(duplicated(form$fixef))]] <- NULL
  }
  
  y.mat = bcObject$y.transform[[1]]
  y <- as.vector(y.mat[!is.na(y.mat)])
  N <- length(y)
  df = df[-which(is.na(y.mat)),]
  x_fx = data.frame(df[,fixvar])
  colnames(x_fx) = fixvar
  x_rd = data.frame(df[,randvar])
  colnames(x_rd) = randvar
  x_fx = do.call("cbind", 
                 sapply(1:ncol(x_fx), function(i) encoded$Index[match(x_fx[,i], encoded$Level)], 
                        simplify=F))
  colnames(x_fx) = fixvar
  
  x_rd = do.call("cbind", 
                 sapply(1:ncol(x_rd), function(i) encoded$Index[match(x_rd[,i], encoded$Level)], 
                        simplify=F))
  colnames(x_rd) = randvar
  modelMat = model.matrix(~ Diet + RIX + DietRIX, df)
  x_d  = modelMat[,grep("[^0-9]$", colnames(modelMat))]
  x_s  = modelMat[,grep("^RIX", colnames(modelMat))]
  x_sd = modelMat[,grep("^DietRIX", colnames(modelMat))]

  standat <-  list(N       = length(y),
                   y       = y, 
                   K_d     = ncol(x_d),
                   K_s     = ncol(x_s),
                   K_sd    = ncol(x_sd),
                   x_d     = x_d,
                   x_s     = x_s,
                   x_sd    = x_sd,
                   SPO     = as.vector(df$PO))
  
  fileName <- "../stan_SPO.stan"
  stan_code <- readChar(fileName, file.info(fileName)$size)
  #cat(stan_code)
  

  warmup=round((floor((iter*0.25)/(iter/10))*(iter/10)))
  thin=iter/1000
  
  stan_code <- readChar(fileName, file.info(fileName)$size)
  try(resStan <- stan(model_code = stan_code, data = standat,
                      chains = chains, iter = iter, warmup = warmup, thin = thin))
  
  return(resStan)
}


  stanmcmc<- As.mcmc.list(resStan)
  summcmc <- summary(stanmcmc)
  
  #armform <- y ~ Diet + PO*RIX + Diet*RIX + PO*Diet*RIX
  #test_stanarm = stan_glm(armform, data=df, prior = lasso(autoscale = T))
  
  par(mfrow=c(3,3))
  coda::traceplot(stanmcmc[[1]][,grep("SPO[[]", colnames(stanmcmc[[1]])),drop=F])
  coda::traceplot(stanmcmc[[1]][,grep("sigma|lambda", colnames(stanmcmc[[1]])),drop=F])
  tset <- sapply(1:(ncol(stanmcmc[[1]])-1), function(x) HPDinterval(stanmcmc[,x,drop=T]), simplify=F)

  
  
  
  
  
  
reg.jags <- jags.model(textConnection(madeModel$modelFin), data=madeModel$data, n.chains = chains, n.adapt = n.adapt)
update(reg.jags, n.iter=n.adapt)
fitJags[[phenotype[i]]]$fit <- coda.samples(reg.jags, variable.names = madeModel$paramUse, thin=thin, n.iter=n.iter)
fitJags[[phenotype[i]]]$madeModel <- madeModel
allnames <- c()
paramUse <- gsub("beta","", unique(unlist(
  lapply(strsplit(varnames(fitJags[[phenotype[i]]]$fit),'[[]'), function(x) x[[1]]))))
wantParam <- intersect(paramUse, c(paste(unique(encoded$Variable)), "grand_mu"))
  
  
  
  
  
  
  jagsFit <- runJagsModels_cov(df=df, chains=chains, testLMER = lmerFit, 
                               n.adapt=n.adapt, n.iter=n.iter, thin=thin, 
                               encoded=encoded, phenotype=phenotype, sq=sq, addS=addS)
  

  
  keep <- c("Diet", "DietRIX","PODietRIX","PORIX", "RIX")
  lmerSum <- lmer.getDecodedSummary(lmerFit, lmerFit$LMformula, phenotypes=phenotype)
  test_cm <- contrast_matrix(encoded)
  
  jagsSum <- list()
  allSummary <- list()
  null.test <- list()
  jagsLmer_compare <- list()
  plotCat <- list()
  madeMod <- list()
  all.mcmc <- list()
  
  contrastMat <- contrast_matrix(encoded, levels=c("RIX","Diet"))
  
  diet.test <- list()
  rix.test<- list()
  
  for(i in 1:length(lmerFit$lmerobj$y.transform)){
    all.mcmc[[phenotype[i]]] <- mcmc.stack(jagsFit[[i]]$fit)
    jagsSum[[i]] <- jags.getDecodedSummary(all.mcmc[[phenotype[i]]], encoded)
    jagsSum[[i]]$Level <- factor(jagsSum[[i]]$Level, levels=jagsSum[[i]]$Level)
    if(!any(is.na(lmerSum[[i]]))) {
      JLSum <- compareSummary(lmerSum[[i]], jagsSum[[i]], keep)
      jagsLmer_compare[[phenotype[i]]] <- data.frame(phenotype=rep(phenotype[i],nrow(JLSum)), JLSum) 
      keepLev <- unique(JLSum$Level)
    } else {
      jagsLmer_compare[[phenotype[i]]] <- NA
      plot=F
      keepLev <- jagsSum[[i]]$Level[which(jagsSum[[i]]$Variable %in% keep)]
    }
    
    null.test[[phenotype[i]]]  <- contrast.null(dfObject=all.mcmc[[phenotype[i]]], levels=keepLev)
    
    null_p <- data.frame(Level = rownames(null.test[[phenotype[i]]]),
                         pval = null.test[[phenotype[i]]]$pval, 
                         signif = null.test[[phenotype[i]]]$signif)
    null_p$Level <- factor(null_p$Level, levels=jagsSum[[i]]$Level)
    null_p <- null_p[order(null_p$Level),]
    mergeTable <- merge(jagsSum[[i]], null_p, by="Level")
    mergeTable$Level <- factor(mergeTable$Level, levels=encoded$Level[which(encoded$Variable %in% keep)])
    mergeTable$Variable <- factor(mergeTable$Variable, levels=unique(encoded$Variable)[which(encoded$Variable %in% keep)])
    mergeTable <- mergeTable[order(mergeTable$Level),]
    allSummary[[phenotype[i]]]  <- cbind(phenotype=rep(phenotype[i],nrow(mergeTable)), mergeTable[order(mergeTable$Variable, mergeTable$Level),])
    
    if(plot) plotCat[[phenotype[i]]] <- plotSummary(lmerSum[[i]], jagsSum[[i]], JLSum, phenotype)
    madeMod[[phenotype[i]]]  <- jagsFit[[i]]$madeModel
    
    if(contrasts){
      diet.test[[phenotype[i]]] <- contrast.test(data=all.mcmc[[phenotype[i]]], encoded, variable="Diet", byVar="RIX", contrast=contrastMat$Diet.contr)
      rix.test[[phenotype[i]]] <- contrast.test(data=all.mcmc[[phenotype[i]]], encoded, variable="RIX", byVar="Diet", contrast=contrastMat$RIX.contr)
    }
  }
  
  return(list(transformation=transf, mcmcObject=all.mcmc, plot=plotCat, allSummary=allSummary,
              null.test=null.test, diet.test=diet.test, rix.test=rix.test, S=S, modelDef=madeMod, 
              JL_compare=jagsLmer_compare, lmer_obj = lmerFit))
}


getFormulas = function(fixvar, randvar, POvar){
  effectlist <- list()
  ind <- 1
  
  if(!any(is.na(fixvar))){
    effectlist[[ind]] <- vector()
    for(i in 1:length(fixvar)){
      effectlist[[ind]] <- c(effectlist[[ind]], fixvar[i])
    }
    ind=ind+1
  }
  if(!any(is.na(randvar))){
    effectlist[[ind]] <- vector()
    for(i in 1:length(randvar)){
      effectlist[[ind]] <- c(effectlist[[ind]], paste0("(1 | ", randvar[i], ")"))
    }
    ind=ind+1
  }
  
  if(!any(is.na(POvar))){   
    effectlist[[ind]] <- vector()
    for(i in 1:length(POvar)){
      effectlist[[ind]] <- c(effectlist[[ind]], paste0("(PO + 0 | ", POvar[i], ")"))
    }
  }
  fitLMER <- list()
  lmerform <- paste("~ -1 +", paste(unlist(effectlist), collapse="+"))
  jagsform <- lmer.to.jags.form(lmerform)
  
  return(list(lmerform=lmerform, jagsform=jagsform))
}

