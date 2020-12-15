source("./plot.hpd.R")


getEncoding <- function(df, terms){
  encoded <- data.frame()
  for(i in 1:length(terms)){
    var <- paste(terms[i])
    len <- length(levels(df[,var]))
    
    tempdf <- data.frame(Level = as.character(levels(df[,var])), 
                         Index = 1:len, 
                         Variable = rep(var,len))
    encoded <- rbind(encoded,tempdf)
  }
  return(encoded)
}

jags.getDecodedSummary <- function(mcmcObject, encoded, narrow=0.5, wide=0.95){
  listNames <- as.character(encoded$Level)
  mu    <- c(colMeans(mcmcObject))
  med   <- apply(coda::HPDinterval(mcmcObject, prob=0.01), 1, median)
  hpd.wide    <- coda::HPDinterval(mcmcObject, prob=wide)
  hpd.narrow  <- coda::HPDinterval(mcmcObject, prob=narrow)
  effectdata <- data.frame(cbind(mu,med,hpd.narrow,hpd.wide))
  effectdata$Level <- rownames(effectdata)
  jagsSummary <- effectdata[intersect(effectdata$Level, listNames),]
  jagsSummary$Variable <- as.character(encoded$Variable[match(jagsSummary$Level, encoded$Level)])
  return(jagsSummary)
}

contrasts.getDecodedSummary <- function(summaryObject, compare, narrow=0.5, wide=0.95){
  testThis <- ifelse(compare=="Diet", "diet.test", "rix.test")
  mcmcObject <- as.mcmc(do.call(cbind,summaryObject[[testThis]]$contrastHPD))
  tempPval <- summaryObject[[testThis]]$pval
  colnames(tempPval) <- paste0("a",colnames(summaryObject[[testThis]]$pval))
  tempsignif <- summaryObject[[testThis]]$signif
  colnames(tempsignif) <- paste0("a",colnames(summaryObject[[testThis]]$signif))
  
  pval <- cbind(melt(tempPval), signif=melt(tempsignif)[,"value"])
  pval$Var2 <- as.character(gsub("a","",pval$Var2))
  
  mu    <- c(colMeans(mcmcObject))
  med   <- apply(coda::HPDinterval(mcmcObject, prob=0.01), 1, median)
  hpd.wide    <- coda::HPDinterval(mcmcObject, prob=wide)
  hpd.narrow  <- coda::HPDinterval(mcmcObject, prob=narrow)
  effectdata <- data.frame(cbind(mu,med,hpd.narrow,hpd.wide))
  effectdata$Level1 <- unlist(strsplit(rownames(effectdata),"[.]"))[c(T, F)]
  effectdata$Level2 <- unlist(strsplit(rownames(effectdata),"[.]"))[c(F, T)]
  if(compare=="Diet"){
    pval$PO[grep("PO",pval$Var1)] <- "PO"
    pval$RIX <- gsub('\\D+','', pval$Var1)
    pval$Level1 <- unlist(strsplit(as.character(pval$Var2),"[.]"))[c(T, F)]
    pval$Level2 <- unlist(strsplit(as.character(pval$Var2),"[.]"))[c(F, T)]
    effectdata$RIX <- gsub('\\D+','', effectdata$Level1)
    effectdata$Level1 <- gsub('[0-9]','', effectdata$Level1)
    effectdata$Level2 <- gsub('[0-9]','', effectdata$Level2)
    effectdata$PO[grep("PO",effectdata$Level1)] <- "PO"
    effectdata$Level1 <- gsub('PO','', effectdata$Level1)
    effectdata$Level2 <- gsub('PO','', effectdata$Level2)
    
    contrastSum <- merge(effectdata, pval[,c("PO","RIX","Level1","Level2","value","signif")], by=c("PO","RIX","Level1","Level2"))
    contrastSum$RIX <- factor(contrastSum$RIX, levels=unique(effectdata$RIX))
    
    contrastSum$Level1 <- factor(contrastSum$Level1, levels=unique(pval$Level1))
    contrastSum$Level2 <- factor(contrastSum$Level2, levels=unique(pval$Level2))
    contrastSum <- contrastSum[order(contrastSum$PO, contrastSum$RIX, contrastSum$Level1, contrastSum$Level2),]
    
  } else {
    pval$PO[grep("PO",pval$Var1)] <- "PO"
    pval$Diet <- gsub('PO','', pval$Var1)
    pval$Level1 <- unlist(strsplit(as.character(pval$Var2),"[.]"))[c(T, F)]
    pval$Level2 <- unlist(strsplit(as.character(pval$Var2),"[.]"))[c(F, T)]
    pval[order(pval$PO, pval$Level1, pval$Level2),]
    effectdata$Diet <- gsub('\\d+','', effectdata$Level1)
    effectdata$PO[grep("PO",effectdata$Diet)] <- "PO"
    effectdata$Diet <- gsub('PO','', effectdata$Diet)
    effectdata$Level1 <- gsub('\\D+','', effectdata$Level1)
    effectdata$Level2 <- gsub('\\D+','', effectdata$Level2)
    
    contrastSum <- merge(effectdata, pval[,c("PO","Diet","Level1","Level2","value","signif")], 
                         by=c("PO","Diet","Level1","Level2"))
    contrastSum$Diet <- factor(contrastSum$Diet, levels=unique(effectdata$Diet))
    contrastSum$Level1 <- factor(contrastSum$Level1, levels=unique(pval$Level1))
    contrastSum$Level2 <- factor(contrastSum$Level2, levels=unique(pval$Level2))
    contrastSum <- contrastSum[order(contrastSum$PO, contrastSum$Diet, contrastSum$Level1, contrastSum$Level2),]
  }
  return(contrastSum)
}



lmer.getDecodedSummary <- function(lmerObject, fixed, rand=T, fix=T){
  lmerSummary <- data.frame(Intercept = numeric(), Level=factor(), Variable=factor())
  if(fix){
    for(i in 1: length(fixed)){
      vec <- fixef(lmerObject)
      levels <- names(vec)[grep(fixed[i], names(vec))]
      levels2 <- gsub(fixed[i],"",levels)
      for (j in 1: length(levels)){
        this_rand <- vec[j]
        tempSum <- data.frame(Intercept = as.numeric(vec[j]), 
                              Level = levels2[j], 
                              Variable = fixed[i])
        lmerSummary <- rbind(lmerSummary, tempSum)
      }
    }
  }
  if (rand==T){
    randEfs <- list()
    for(j in 1:length(names(ranef(lmerObject)))){
      tempName <- names(ranef(lmerObject))[j]
      randEfs[[gsub("`","",tempName)]] <- ranef(lmerObject)[[tempName]]
      colnames(randEfs[[tempName]])[which(colnames(randEfs[[tempName]]) == "(Intercept)")] <- ""
      for (k in 1:ncol(randEfs[[tempName]])){
        addToSum <- data.frame(Intercept = as.numeric(randEfs[[j]][,k]), 
                               Level = paste0(colnames(randEfs[[j]])[k], gsub(":","",rownames(randEfs[[j]]))), 
                               Variable = paste0(colnames(randEfs[[j]])[k], gsub(":","",tempName)))
        lmerSummary <- rbind(lmerSummary, addToSum)
      }
    }
  }
  return(lmerSummary)
}

###################################
# Model with different parameters #
###################################

runJagsModels_cov <- function(phenotype, dataf, allparam, encoded, testLMER=NULL,
                          OF=F, LD=F, SIH=F, FST=F, Stress=F, Match=F){
  
  check.reg <- list()
  for (i in 1:length(phenotype)){
    if(is.null(testLMER)){
      pheno = as.vector(scale(dataf[[phenotype]]))
    } else {
      #transPhenoTemp <- paste0(phenotype[i],"_trans")
      use <- which(is.na(dataf[,phenotype[i]]) == F)
      matnut_use <- dataf[use,]
      
      pheno <- testLMER[[phenotype[i]]]$phen_1$y.transformed  #testLMER[[phenotype[i]]]$y.transform
      
      #pheno = matnut_use[[transPhenoTemp]]
    }
    
    N <- length(pheno)
    y <- as.vector(pheno)
    
    covariates <- colnames(matched_df)[which(colnames(matched_df) %in% unique(encoded$Variable))]
    
    for(i in 1:length(random)){
      random_ef[[i]] <- as.numeric(dataf[[random[i]]])
      nrandom_ef[i] <- length(which(encoded$Variable == random[i])) 
    }
    for(i in 1:length(fixed)){
      fixed_ef[[i]] <- as.numeric(dataf[[fixed[i]]])
      nfixed_ef[i] <- length(which(encoded$Variable == fixed[i])) 
    }
    diet <- as.numeric(matnut_use$Diet)
    rix <- as.numeric(matnut_use$RIX)
    dietrix <- as.numeric(matnut_use$DietRIX)
    poe <- as.numeric(matnut_use$PO)
    be <- as.numeric(matnut_use$BehaviorBatch)
    dam <- as.numeric(matnut_use$DamID)
    sire <- as.numeric(matnut_use$SireID)
    #cage <- as.numeric(matnut_use$Cage)
  
    ndiet <- length(which(encoded$Variable == "Diet")) 
    nrix <- length(which(encoded$Variable == "RIX"))
    ndietrix <- length(which(encoded$Variable == "DietRIX"))
    nbatch <- length(which(encoded$Variable == "BehaviorBatch")) 
    ndam <- length(which(encoded$Variable == "DamID")) 
    nsir <- length(which(encoded$Variable == "SireID")) 
    #ncage <- length(which(encoded$Variable == "Cage"))
    
    if (Match){
      modelUse = modelMatch
      data = list("N"=N, "y"=y, "rix"=rix, "nrix"=nrix, "dietrix"=dietrix, "ndietrix"=ndietrix)
      paramUse = c("DietRIX", "RIX", "sigInv2", "tauDRInv2", "tauInv2")
    } else if (OF){
      modelUse = modelOF
      nbox <- length(which(encoded$Variable == "OFBox"))
      of <- as.numeric(matnut_use$OFBox)
      
      data = list("N"=N, "y"=y, "rix"=rix, "nrix"=nrix, "diet"=diet, "ndiet"=ndiet, "dietrix"=dietrix, "ndietrix"=ndietrix,
                  "poe"=poe, "nsir"=nsir, "sire"=sire,"ndam"=ndam, "dam"=dam, "be"=be, "nbatch"=nbatch, "of"=of, "nbox"=nbox)
      paramUse = allparam[-match(c("FSTChamber", "LDChamber","RestraintExperimenter","RestraintOrder", "SIHOrder"), allparam)]
    } else if (FST){
      modelUse = modelFST
      nbox <- length(which(encoded$Variable == "FSTChamber"))
      fst <- as.numeric(matnut_use$FSTChamber)
      
      data = list("N"=N, "y"=y, "rix"=rix, "nrix"=nrix, "diet"=diet, "ndiet"=ndiet, "dietrix"=dietrix, "ndietrix"=ndietrix,
                  "poe"=poe, "nsir"=nsir, "sire"=sire,"ndam"=ndam, "dam"=dam, "be"=be, "nbatch"=nbatch, "fst"=fst, "nbox"=nbox)
      paramUse = allparam[-match(c("OFBox", "LDChamber","RestraintExperimenter","RestraintOrder", "SIHOrder"), allparam)]
    } else if (LD){
      modelUse = modelLD
      nbox <- length(which(encoded$Variable == "LDChamber"))
      ld <- as.numeric(matnut_use$LDChamber)
      
      data = list("N"=N, "y"=y, "rix"=rix, "nrix"=nrix, "diet"=diet, "ndiet"=ndiet, "dietrix"=dietrix, "ndietrix"=ndietrix,
                  "poe"=poe, "nsir"=nsir, "sire"=sire,"ndam"=ndam, "dam"=dam, "be"=be, "nbatch"=nbatch, "ld"=ld, "nbox"=nbox)
      paramUse = allparam[-match(c("FSTChamber", "OFBox", "RestraintExperimenter","RestraintOrder", "SIHOrder"), allparam)]
    } else if (SIH){
      modelUse = modelSIH
      norder <- length(which(encoded$Variable == "SIHOrder"))
      sih <- as.numeric(matnut_use$SIHOrder)
      
      data = list("N"=N, "y"=y, "rix"=rix, "nrix"=nrix, "diet"=diet, "ndiet"=ndiet, "dietrix"=dietrix, "ndietrix"=ndietrix,
                  "poe"=poe, "nsir"=nsir, "sire"=sire,"ndam"=ndam, "dam"=dam, "be"=be, "nbatch"=nbatch, "sih"=sih, "norder"=norder)
      paramUse = allparam[-match(c("FSTChamber", "LDChamber","RestraintExperimenter","RestraintOrder", "OFBox"), allparam)]
    } else if (Stress){
      modelUse = modelStress
      norder <- length(which(encoded$Variable == "RestraintOrder"))
      nex <- length(which(encoded$Variable == "RestraintExperimenter"))
      restex <- as.numeric(matnut_use$RestraintExperimenter)
      restord <- as.numeric(matnut_use$RestraintOrder)

      data = list("N"=N, "y"=y, "rix"=rix, "nrix"=nrix, "diet"=diet, "ndiet"=ndiet, "dietrix"=dietrix, "ndietrix"=ndietrix,
                  "poe"=poe, "nsir"=nsir, "sire"=sire,"ndam"=ndam, "dam"=dam, 
                  "be"=be, "nbatch"=nbatch, "restex"=restex, "nex"=nex, "restord"=restord, "norder"=norder)
      paramUse = allparam[-match(c("FSTChamber", "LDChamber","OFBox","SIHOrder"), allparam)] #, "RestraintOrder"
    } else {
      modelUse = modelWeight
      data = list("N"=N, "y"=y, "rix"=rix, "nrix"=nrix, "diet"=diet, "ndiet"=ndiet, "dietrix"=dietrix, "ndietrix"=ndietrix,
                  "poe"=poe, "nsir"=nsir, "sire"=sire,"ndam"=ndam, "dam"=dam, "be"=be,"nbatch"=nbatch)
      paramUse = allparam[-match(c("FSTChamber", "LDChamber","RestraintExperimenter","RestraintOrder", 
                                   "SIHOrder", "OFBox"), allparam)]
      
    }
    
    #############
    # Run model #
    #############
    
    chains <- 2 
    reg.jags <- jags.model(textConnection(modelUse), data=data, n.chains = chains, n.adapt = 2500)
    update(reg.jags, n.iter=2500)
    check.reg[[phenotype[i]]] <- coda.samples(reg.jags, variable.names = paramUse, thin=10, n.iter=10000)
  }
  allnames <- c()
  paramUse <- unique(unlist(strsplit(varnames(check.reg[[i]]),'[[]'))[c(T, F)])
  wantParam <- intersect(paramUse, unique(encoded$Variable)) 
  for(i in 1:length(wantParam)){
    allnames <- c(allnames, as.character(encoded$Level[which(encoded$Variable == wantParam[i])])) 
  }
  for(i in 1:length(check.reg)){
    varnames(check.reg[[i]])[1:length(allnames)] <- allnames
  } 
  return(check.reg)
}

##############
# LMER model #
##############
runLMERModels_cov <- function(phenotype, dataf, tryLam=1, normd=T,
                          Match=F, OF=F, LD=F, SIH=F, FST=F, Stress=F, checkAnova = T){
  check.mod_results <- list()
  init <- "~ -1"
  dt <- "Diet"
  dietRix <- "(1 | RIX) + (1 | DietRIX)"
  randEf <- "(1 | DamID) + (1 | SireID) + (1 | BehaviorBatch)"
  podietrix <- "(1 | RIX) + (PO + 0 | RIX) + (1 | Diet:RIX) + (PO + 0 | Diet:RIX)"
  
  
  if(OF){
    indvariable <- as.formula(paste(init, dt, "OFBox", randEf, podietrix, sep="+"))
  } else if(LD){
    indvariable <- as.formula(paste(init, dt, "LDChamber", randEf, podietrix, sep="+")) 
  } else if(SIH){
    indvariable <- as.formula(paste(init, dt, "SIHOrder", randEf, podietrix, sep="+")) 
  } else if(FST){
    indvariable <- as.formula(paste(init, dt, "FSTChamber", randEf, podietrix, sep="+"))  
  } else if(Stress){
    indvariable <- as.formula(paste(init, dt, "RestraintExperimenter + RestraintOrder", randEf, podietrix, sep="+"))  
  } else {
    indvariable <- as.formula(paste(init, dt, randEf, podietrix, sep="+"))  
  }
  for (i in 1:length(phenotype)){  
    pheno <- paste(phenotype[i])
    use <- which(is.na(dataf[,pheno]) == F)
    matnut_use <- dataf[use,]
    #check.mod_results[[phenotype[i]]] <- BC.model(y.mat = matnut_use[,pheno], matnut_use, indvariable, 
    #                                                      transformParams = getMatnutTransformParams(normd=normd, tryLam = tryLam))
    if(Match){
      matnut_use$phen_use <- scale(matnut_use[,pheno])
      form <- formula(paste("phen_use", indvariable))
      check.mod_results[[phenotype[i]]]$phen_1 <- list()
      check.mod_results[[phenotype[i]]]$phen_1$fit <- lmer(form, data=matnut_use)
      check.mod_results[[phenotype[i]]]$phen_1$y.transformed <- as.vector(matnut_use$phen_use)
      
    } else {
      check.mod_results[[phenotype[i]]] <- fit.model.bc$fit(y.mat = matnut_use[,pheno], matnut_use, indvariable, 
                                                            transformParams = getMatnutTransformParams(normd=normd, tryLam = tryLam))
    }
    
    print(paste(pheno, "fit by LMER"))
    matnut_use <- cbind(matnut_use, y.trans=check.mod_results[[phenotype[i]]]$phen_1$y.transformed)
    form <- as.formula(paste0("y.trans", indvariable))
    lobj <- lmer(form, data=matnut_use)
    #check.mod_results[[phenotype[i]]]$rand_signif <- rand(lobj)
    #check.mod_results[[phenotype[i]]]$fix_signif <- fixef(lobj)
    check.mod_results[[phenotype[i]]]$phen_1$fit <- lobj
  }
  return(check.mod_results)
}

BC.model <- function(y.mat, data, invariable, transformParams=getMatnutTransformParams()){
  p_val <- c()
  y_temp <- list()
  store_fit <- list()
  if(transformParams$normalizeBeforeTransform == T){
    y.mat <- scale(y.mat)
  }
  for(i in 1:length(transformParams$lambdasToTry)){
    lambda <- transformParams$lambdasToTry[i]
    y_temp[[i]] <- c()
    if (lambda == 0){
      y.min <- min(y.mat, 0)
      lamb2 <- ifelse(y.min <= 0, abs(y.min)+1, 0)
      y_temp[[i]] <- log(y.mat+lamb2)
    } else {
      y.min <- min(y.mat) 
      lamb2 <- ifelse(lambda<1 && y.min<=0, abs(y.min)+1, 0)
      y_temp[[i]] <- ((y.mat+lamb2)^lambda - 1) / lambda
    }
    data_use <- cbind(data, transf = y_temp[[i]])
    form <- as.formula(paste0("transf", indvariable))
    fit <- lmer(form, data=data_use)
    store_fit[[i]] <- fit
    res <- residuals(fit)
    normed <- shapiro.test(res)
    warn <- fit@optinfo$conv$lme4$messages
    #warn2 <- is.na(fit@optinfo$warnings)
    p_val[i] <- ifelse(is.null(warn), normed$p.value, 0)
  }
  ind <- which(p_val == max(p_val))
  y.transform <- y_temp[[ind]]

  if(transformParams$normalizeAfterTransform == T){
    y.transform <- scale(y.transform)
  }
  phen_1 <- store_fit[[ind]]
  
  return(list(p.val=p_val[ind], lambda=transformParams$lambdasToTry[ind], y.transform=y.transform, phen_1=phen_1))
}


###############
# Get summary #
###############

allSummaries_cov <- function(phenotype, dataf, allparam, encoded, tryLam=1, Match=F, checkAnova=F)
{
  if(Match == F){
    OF <- ifelse(length(grep("OF",phenotype))==0,F,T)
    LD <- ifelse(length(grep("LD",phenotype))==0,F,T)
    SIH <- ifelse(length(grep("SIH",phenotype))==0,F,T)
    FST <- ifelse(length(grep("FST",phenotype))==0,F,T)
    Stress <- ifelse(length(grep("CORT",phenotype))==0,F,T)
    fix=T
  } else {
    OF=F
    LD=F
    SIH=F
    FST=F
    Stress=F
    fix=F
  }
  testLMER <- runLMERModels_cov(phenotype=phenotype, dataf, tryLam = tryLam, Match=Match,
                                OF=OF, LD=LD, SIH=SIH, FST=FST, Stress=Stress, checkAnova=checkAnova)
  fixed = ifelse(Match, "","Diet")
  lmerSum <- lmer.getDecodedSummary(lmerObject=testLMER[[phenotype]]$phen_1$fit, fixed=fixed, rand=T, fix=fix)
  transf <- c("lambda" = testLMER[[phenotype]]$phen_1$lambda,  #testLMER[[phenotype]]$lambda
              "pval" = testLMER[[phenotype]]$phen_1$best.pval)
  testJags <- runJagsModels_cov(phenotype=phenotype, dataf, allparam, encoded, testLMER, Match=Match,
                                OF=OF, LD=LD, SIH=SIH, FST=FST, Stress=Stress)
  all.reg <- testJags[[phenotype]]
  all.mcmc <- mcmc.stack(all.reg)
  jagsSum <- jags.getDecodedSummary(all.mcmc, encoded)
  allSummaries <- compareSummaries(lmerSum, jagsSum)
  
  keep <- c("Diet", "DietRIX","PODietRIX","PORIX", "RIX")
  
  getSummaries <- allSummaries[which(allSummaries$Variable %in% keep),]
  
  null.test <- contrast.null(dfObject=as.data.frame(all.mcmc), levels=unique(allSummaries$Level))
  if (!Match){
    diet.test <- contrast.test(dfObject=as.data.frame(all.mcmc), encoded, variable="Diet", byVar="RIX", contrast=diet.contr)
    rix.test <- contrast.test(dfObject=as.data.frame(all.mcmc), encoded, variable="RIX", byVar="Diet", contrast=rix.contr)
  } else {
    diet.test <- NULL
    rix.test <- NULL
  }
  
  jagsSum$Level <- factor(jagsSum$Level, levels=jagsSum$Level)
  
  null_p <- data.frame(Level = names(null.test$pval),
                       pval = null.test$pval, 
                       signif = null.test$signif)
  null_p$Level <- factor(null_p$Level, levels=jagsSum$Level)
  null_p <- null_p[order(null_p$Level),]
  mergeTable <- merge(jagsSum, null_p, by="Level")
  mergeTable$Level <- factor(mergeTable$Level, levels=encoded$Level)
  mergeTable$Variable <- factor(mergeTable$Variable, levels=unique(encoded$Variable))
  
  summary_table <- cbind(phenotype=rep(phenotype,nrow(mergeTable)), mergeTable[order(mergeTable$Variable, mergeTable$Level),])
  jagsLmer_compare <- cbind(phenotype=rep(phenotype,nrow(getSummaries)), getSummaries) 
  
  ####
  
  plotCat <- combinedCat(lmerSum, jagsSum, getSummaries, phenotype)
  return(list(transformation=transf, mcmcObject=all.mcmc, plot=plotCat, allSummaries=allSummaries,
              null.test=null.test, diet.test=diet.test, rix.test=rix.test, 
              summary_table=summary_table, JL_compare=jagsLmer_compare, lmer_obj = testLMER[[phenotype]]))
}



compareSummaries <- function(lmerSum, jagsSum){
  common <- intersect(lmerSum$Level, jagsSum$Level)
  lmerCom <- lmerSum[which(lmerSum$Level %in% common),]
  jagsCom <- jagsSum[which(jagsSum$Level %in% common),]
  merged <- merge(lmerCom, jagsCom, by="Level")
  Diff <- merged$Intercept - merged$mu
  Ratio <- Diff/merged$Intercept
  diffSummary <- cbind(merged, Diff, Ratio)
  diffSummary$Level <- factor(diffSummary$Level, levels=jagsCom$Level)
  diffSummary <- diffSummary[order(diffSummary$Level), c("Level", "Intercept", "mu", "Diff", "Ratio", "Variable.x")]
  colnames(diffSummary) <- c("Level", "LMER_est", "JAGS_est", "Diff", "Ratio", "Variable")
  
  return(diffSummary)
}

combinedCat <- function(lmerSum, jagsSum, getSummaries, phenotype){
  jagsSumTemp <- jagsSum[which(jagsSum$Level %in% getSummaries$Level),]
  lmerSumTemp <- lmerSum[which(lmerSum$Level %in% getSummaries$Level),]
  jagsPlot <- plot.inter.ci(med=jagsSumTemp$med, mu=jagsSumTemp$mu, 
                            hpd.narrow = cbind(jagsSumTemp$lower, jagsSumTemp$upper), 
                            hpd.wide = cbind(jagsSumTemp$lower.1, jagsSumTemp$upper.1), 
                            names=jagsSumTemp$Level, order=3,
                            col.midvals="white", pch.midvals="|", addline = F, wide=T, 
                            grouped=jagsSumTemp$Variable)
  
  comPlot <- jagsPlot + geom_point(data=lmerSumTemp, aes(x=Level, y=Intercept), col="red", size=3) +
    ggtitle(paste("LMER and JAGS comparison for", phenotype))
  return(comPlot)
}

contrast.test <- function(dfObject, encoded, variable, byVar, contrast, 
                          sig.val=0.05, sig=T)
{
  ### set up contrasts ###
  
  comparisons <- as.character(encoded$Level[which(encoded$Variable == variable)])
  byLevel <- encoded$Level[which(encoded$Variable == byVar)]
  tempNames <- encoded$Level[grep(variable, encoded$Variable)]
  tempLevels <- unique(encoded$Variable[grep(variable, encoded$Variable)])
  useNames <- matrix(NA, nrow=length(comparisons), 
                     ncol=(length(tempNames)/length(comparisons)) + 1)
  count=1
  for(i in 1:length(tempLevels)){
    if (tempLevels[i] == variable){
      useNames[,count] <- comparisons
      count=count+1
    } else {
      tempCategory <- encoded$Level[which(encoded$Variable == tempLevels[i])]
      if(length(grep(byLevel[1],tempCategory))>0){
        for(j in 1:length(byLevel)){
          if(variable == "Diet"){
            pasteThis = paste0(byLevel[j],"$")
          } else {pasteThis = byLevel[j]}
          string <- as.character(tempCategory[grep(pasteThis,tempCategory)])
          useNames[,count] <- rep(c(string, NA), length.out=nrow(useNames))
          
          count=count+1
        }
      } else { useNames[,count] <- as.character(tempCategory)
      count = count+1
      }
    }
  }
  
  ### generate things to return ###
  
  signif <- matrix(NA, nrow=ncol(useNames), ncol=nrow(contrast))
  pval <- matrix(NA, nrow=ncol(useNames), ncol=nrow(contrast))
  dir <- matrix(NA, nrow=ncol(useNames), ncol=nrow(contrast))
  title <- c()
  contrastHPD <- list()
  for (k in 1:ncol(useNames)){
    testname <- useNames[,k]
    #pval[[k]] <- vector()
    for (i in 1:nrow(contrast)){
      compareTemp <- comparisons[contrast[i,]]
      if(variable == "RIX"){
        compareThis = paste0(compareTemp,"$")
      } else {compareThis = compareTemp
      }
      compare <- unique(testname[c(grep(compareThis[1], testname), 
                                   grep(compareThis[2], testname))])
      lab <- paste0(compare[1], ".", compare[2])
      numSamps <- nrow(dfObject)
      if (k==1){
        title[i] <- lab }
      if (length(compare)==2){
        for (j in 1:numSamps){
          contrastHPD[[lab]][j] <- dfObject[,compare[1]][j] - dfObject[,compare[2]][j]
        }
        
        if (sig==T){
          over0 <- length(which(contrastHPD[[lab]]>0))
          pval[k,i] <- min(over0, numSamps - over0) / numSamps
          dir[k,i] <- ifelse(over0 > numSamps - over0, "+","-")
          if (pval[k,i] < sig.val/10) {
            signif[k,i] <- "**"
          } else if (pval[k,i] < sig.val){
            signif[k,i] <- "*"
          } else {signif[k,i] <- "-"}
        }
      }
    }
  }
  
  colnames(signif) <- title
  rownames(signif) <- gsub(comparisons[1],"",useNames[1,])
  colnames(pval) <- title
  rownames(pval) <- gsub(comparisons[1],"",useNames[1,])
  
  if (sig){
    return(list(signif=signif, pval=pval,contrastHPD=contrastHPD, dir=dir))
  } else {
    return(list(contrastHPD=contrastHPD, dir=dir))
  }
}

contrast.null <- function(dfObject, levels, sig.val=0.05, sig=T)
{
  ### set up contrasts ###
  useDF <- dfObject[,as.character(levels)]
  numSamps <- nrow(useDF)
  
  ### generate things to return ###
  
  signif <- c()
  pval <- c()
  dir <- c()
  for (i in 1:ncol(useDF)){
    over0 <- length(which(useDF[,i]>0))
    pval[i] <- min(over0, numSamps - over0) / numSamps
    dir[i] <- ifelse(over0 > numSamps - over0, "+","-")
    if (pval[i] < sig.val/10) {
      signif[i] <- "**"
    } else if (pval[i] < sig.val){
      signif[i] <- "*"
    } else {signif[i] <- "-"}
  }
  
  names(signif) <- colnames(useDF)
  names(pval) <- colnames(useDF)
  
  return(list(signif=signif, pval=pval, dir=dir))
}


getMatnutTransformParams <- function(normd=T, tryLam=1)
{
  out = list()
  out$lambdasToTry = tryLam
  out$lambdaPerY      = NULL
  if (normd){
    out$normalizeBeforeTransform = T
    out$normalizeAfterTransform  = T
  }
  else {
    out$normalizeBeforeTransform = F
    out$normalizeAfterTransform  = F
  }
  out$extremelb      = -3
  out$extremeub      = 3
  return(out)
}


getPlot.compare.helper <- function(dfplot,ptypes,byVar)
{
  
  limz = c(-3,3)
  dfplot$pval[-log10(dfplot$pval)>(max(limz)-.01)] = 10^(-(max(limz)-.01))
  
  dfplot$dir.string = ""
  
  dfplot$p.with.direction = -log10(dfplot$pval) * dfplot$direction
  dfplot$stringg = 32
  dfplot$stringg[dfplot$pval<.05] = 8 ##symbol for sig.                                                                                                                                            
  
  aplot = ggplot(dfplot, aes(x=effectName1, y=effectName2, 
                             fill = p.with.direction, label=stringg, shape = stringg))
  aplot = aplot + geom_tile()
  aplot = aplot + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))
  aplot = aplot + scale_shape_identity()
  aplot = aplot + xlab(paste(byVar,"1")) + ylab(paste(byVar,"2"))
  aplot = aplot + scale_fill_gradient2(low = "red", high = "blue",
                                       name = expression('-Log'[10]*'(p)'%.%'EffectDirection'),
                                       limits =limz)
  aplot = aplot + theme_classic() + geom_point(size = 2.5, color="white") + 
    facet_grid(Level~LevelPO) + ggtitle(paste(byVar,"contrast plot for",ptypes[i]))
  aplot = aplot 
  return(aplot)
}


