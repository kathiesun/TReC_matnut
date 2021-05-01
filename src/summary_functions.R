source("jags_functions.R")
source("lmer_functions_rna.R")
source("contrast_functions.R")
source("plot.hpd_ks.R")

getEncoding <- function(df, terms){
  encoded <- data.frame()
  for(i in 1:length(terms)){
    var <- which(colnames(df) == paste(terms[i]))
    if(is.factor(unlist(df[,var]))){
      temp = unlist(df[,var])
    } else {
      temp = as.factor(df[,var])
    }
    len <- length(levels(temp))
    
    tempdf <- data.frame(Level = as.character(levels(temp)), 
                         Index = 1:len, 
                         Variable = rep(paste(terms[i]),len))
    encoded <- rbind(encoded,tempdf)
  }
  return(encoded)
}

mcmc.stack <- function(coda.object){
  chains <- length(coda.object)
  all.mcmc <- list()
  if(chains == 1){
    all.mcmc[[1]] <- rbind(coda.object)
  } else {
    all.mcmc <- rbind(coda.object[[1]], coda.object[[2]])
    if (chains>2){
      for (j in 3:chains){
        all.mcmc <- rbind(all.mcmc, coda.object[[j]])
      }
    }
  }
  all.mcmc <- as.mcmc(all.mcmc)
  return(all.mcmc)
}


#####################
# Summaries of fits #
#####################


makeSummary <- function(datalist, phenotype=NA, tryLam=1, Match=F, sq=F, normd=T, addS=c("off", "force", "fit"), 
                        chains=2, n.adapt=2500, n.iter=10000, thin=10, plot=T,  
                        randvar=NA, fixvar=NA, POvar=NA, contrasts=F,
                        encoded=NA){
  if(class(datalist) == "list"){
    dataf <- datalist$df
    encoded <- datalist$encoded
    if(any(is.na(phenotype))){
      phenotype <- datalist$ptypes
    } 
  } else {
    dataf <- datalist
  }

  addS <- addS[1]
  S <- c()
  madeMod <- list()
  lmerFit <- runLMERModels_cov(phenotype=phenotype, dataf=dataf, tryLam=tryLam, normd=normd, Match=Match, 
                               randvar=randvar, fixvar=fixvar, POvar=POvar)
  if(length(grep("DamID", colnames(dataf))) > 1){
    lmerFit$JAGSformula =  "~ 1 + (1 | DamID.2) + (-1* | DamID.1) + (1 | RIX) + (1 | DietRIX)"
  } 
  jagsFit <- runJagsModels_cov(datalist=dataf, chains=chains, testLMER = lmerFit, 
                               n.adapt=n.adapt, n.iter=n.iter, thin=thin, 
                               encoded=encoded, phenotype=phenotype, sq=sq, addS=addS)
  
  transf <- data.frame(lambda = lmerFit$lmerobj$lambda,  
                       pval = lmerFit$lmerobj$p.val)
  
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

makeMatchedSummary <- function(datalist, phenotype=NA, tryLam=1, sq=F, normd=T, 
                               chains=2, n.adapt=2500, n.iter=10000, thin=10, plot=T,  
                               randvar=NA, fixvar=NA, POvar=NA, contrasts=F,
                               encoded=NA){
  if(class(datalist) == "list"){
    dataf <- datalist$df
    encoded <- datalist$encoded
    if(any(is.na(phenotype))){
      phenotype <- datalist$ptypes
    } 
  } else {
    dataf <- datalist
  }

  lmerFit <- runLMERModels_cov(phenotype=phenotype, dataf=dataf, tryLam=tryLam, normd=normd, Match=Match, 
                               randvar=randvar, fixvar=fixvar, POvar=POvar)
  if(length(grep("DamID", colnames(dataf))) > 1){
    lmerFit$JAGSformula =  "~ 1 + (1 | DamID.2) + (-1* | DamID.1) + (1 | RIX) + (1 | DietRIX)"
  } 
  jagsFit <- runJagsModels_cov(datalist=dataf, chains=chains, testLMER = lmerFit, 
                               n.adapt=n.adapt, n.iter=n.iter, thin=thin, 
                               encoded=encoded, phenotype=phenotype, sq=sq, addS=addS)
  
  transf <- data.frame(lambda = lmerFit$lmerobj$lambda,  
                       pval = lmerFit$lmerobj$p.val)
  
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

compareSummary <- function(lmerSum, jagsSum, keep=NA){
  if(any(is.na(keep))) keep <- unique(lmerSum$Variable)
  
  common <- intersect(as.character(lmerSum$Level), as.character(jagsSum$Level))
  lmerCom <- lmerSum[which(lmerSum$Level %in% common),]
  jagsCom <- jagsSum[which(jagsSum$Level %in% common),]
  merged <- merge(lmerCom, jagsCom, by="Level")
  Diff <- merged$Intercept - merged$mu
  Ratio <- Diff/merged$Intercept
  diffSummary <- cbind(merged, Diff, Ratio)
  diffSummary$Level <- factor(diffSummary$Level, levels=jagsCom$Level)
  diffSummary <- diffSummary[order(diffSummary$Level), c("Level", "Intercept", "mu", "Diff", "Ratio", "Variable.x")]
  colnames(diffSummary) <- c("Level", "LMER_est", "JAGS_est", "Diff", "Ratio", "Variable")

  getSummary <- diffSummary[which(diffSummary$Variable %in% keep),]
  
  compareSum <- getSummary
  
  
  return(compareSum)
}

plotSummary <- function(lmerSum, jagsSum, getSummary, phenotype){
  jagsSumTemp <- jagsSum[which(jagsSum$Level %in% getSummary$Level),]
  lmerSumTemp <- lmerSum[which(lmerSum$Level %in% getSummary$Level),]
  jagsPlot <- plot.inter.ci(med=jagsSumTemp$med, mu=jagsSumTemp$mu, 
                            hpd.narrow = cbind(jagsSumTemp$lower, jagsSumTemp$upper), 
                            hpd.wide = cbind(jagsSumTemp$lower.1, jagsSumTemp$upper.1), 
                            names=jagsSumTemp$Level, order=3,
                            col.midvals="white", pch.midvals="|", addline = F, wide=T, 
                            grouped=jagsSumTemp$Variable, ordered=F)
  comPlot <- jagsPlot$plot + geom_point(data=lmerSumTemp, aes(x=Level, y=Intercept), col="red", size=3) +
    ggtitle(paste("LMER and JAGS comparison for", phenotype))
  return(comPlot)
}





