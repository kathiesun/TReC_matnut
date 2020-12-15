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



lmer.getDecodedSummary <- function(lmerObject, fixed, rand=T){
  lmerSummary <- data.frame(Intercept = numeric(), Level=factor(), Variable=factor())
  for(i in 1: length(fixed)){
    vec <- fixef(lmerObject$fit)
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
  if (rand==T){
    randEfs <- list()
    for(j in 1:length(names(ranef(lmerObject$fit)))){
      tempName <- names(ranef(lmerObject$fit))[j]
      randEfs[[gsub("`","",tempName)]] <- ranef(lmerObject$fit)[[tempName]]
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

runJagsModels <- function(phenotype, orderdat, allparam, encoded, testLMER=NULL,
                          DietRIX=T, PO=T, PODietRIX=T, Batch=F, Dam=F, Cage=F){
  
  check.reg <- list()
  for (i in 1:length(phenotype)){
    PhenoTemp = paste(phenotype[i])
    transPhenoTemp <- paste0(PhenoTemp,"_trans")
    use <- which(is.na(orderdat[,PhenoTemp]) == F)
    matnut_use <- orderdat[use,]
    
    matnut_use[[transPhenoTemp]] <- testLMER[[PhenoTemp]]$phen_1$y.transformed
    
    pheno = matnut_use[[transPhenoTemp]]
    if(is.null(testLMER)){
      pheno = as.vector(scale(pheno))
    }
    
    N <- length(pheno)
    y <- pheno
    
    diet <- as.numeric(matnut_use$Diet)
    rix <- as.numeric(matnut_use$RIX)
    batch <- as.numeric(matnut_use$BreedingBatch)
    dam <- as.numeric(matnut_use$DamID)
    cage <- as.numeric(matnut_use$Cage)
    
    
    ndiet <- length(which(encoded$Variable == "Diet")) 
    nrix <- length(which(encoded$Variable == "RIX"))
    nbatch <- length(which(encoded$Variable == "BreedingBatch")) 
    ndam <- length(which(encoded$Variable == "DamID")) 
    ncage <- length(which(encoded$Variable == "Cage"))
    
    if (DietRIX){
      dietrix <- as.numeric(matnut_use$DietRIX)
      ndietrix <- length(which(encoded$Variable == "DietRIX"))
      if (PO){
        poe <- matnut_use$PO
        data = list("N"=N, "rix"=rix, "y"=y, "nrix"=nrix, "diet"=diet, "ndiet"=ndiet, "poe"=poe, 
                    "dietrix"=dietrix, "ndietrix"=ndietrix)
        if (PODietRIX){
          if(Dam){
            if(Cage){
              if(Batch){
                modelUse = modelBreedBatch
                paramUse = allparam
                data = list("N"=N, "rix"=rix, "y"=y, "nrix"=nrix, "diet"=diet, "ndiet"=ndiet, "poe"=poe, 
                            "dietrix"=dietrix, "ndietrix"=ndietrix, 
                            "nbatch"=nbatch, "batch"=batch,"ncage"=ncage, "cage"=cage,"ndam"=ndam, "dam"=dam)
                
              } else {
                modelUse = modelCage
                paramUse = allparam[-match(c("BreedingBatch", "tauBBInv2"), allparam)]
                data = list("N"=N, "rix"=rix, "y"=y, "nrix"=nrix, "diet"=diet, "ndiet"=ndiet, "poe"=poe, 
                            "dietrix"=dietrix, "ndietrix"=ndietrix, 
                            "ncage"=ncage, "cage"=cage,"ndam"=ndam, "dam"=dam)
              }
            } else {
              modelUse = modelDam
              paramUse = allparam[-match(c("Cage", "tauCageInv2","BreedingBatch","tauBBInv2"), allparam)]
              data = list("N"=N, "rix"=rix, "y"=y, "nrix"=nrix, "diet"=diet, "ndiet"=ndiet, "poe"=poe, 
                          "dietrix"=dietrix, "ndietrix"=ndietrix, "ndam"=ndam, "dam"=dam)
            }
          } else { 
          modelUse = modelFull
          paramUse = allparam[-match(c("BreedingBatch", "tauBBInv2","DamID","tauDamInv2","Cage", 
                                       "tauCageInv2","DamID","tauDamInv2"), allparam)]
          }
        } else {
          modelUse = modelNOPODietRIX
          paramUse = allparam[-match(c("PODietRIX", "tauPDRInv2", "BreedingBatch", "tauBBInv2",
                                       "Cage", "tauCageInv2","DamID","tauDamInv2"), allparam)]
          data = list("N"=N, "rix"=rix, "y"=y, "nrix"=nrix, "diet"=diet, "ndiet"=ndiet, "poe"=poe, 
                      "dietrix"=dietrix, "ndietrix"=ndietrix)
        }
      } else {
        data = list("N"=N, "rix"=rix, "y"=y, "nrix"=nrix, "diet"=diet, "ndiet"=ndiet,
                    "dietrix"=dietrix, "ndietrix"=ndietrix)
        modelUse = modelNOPO
        paramUse <- allparam[-match(c("PODietRIX", "tauPDRInv2", "PORIX", "tauPEInv2", 
                                      "BreedingBatch", "tauBBInv2","DamID","tauDamInv2",
                                      "Cage", "tauCageInv2","DamID","tauDamInv2"), allparam)]
        
      }
    } else {
      data = list("N"=N, "rix"=rix, "y"=y, "nrix"=nrix, "diet"=diet, "ndiet"=ndiet)
      modelUse = modelNOrand
      paramUse <- allparam[-match(c("PODietRIX", "tauPDRInv2", "PORIX", "tauPEInv2", "DietRIX", 
                                    "tauDRInv2", "BreedingBatch", "tauBBInv2","DamID","tauDamInv2",
                                    "Cage", "tauCageInv2","DamID","tauDamInv2"), allparam)]
    }
    #############
    # Run model #
    #############
    chains <- 2 
    reg.jags <- jags.model(textConnection(modelUse), data=data, n.chains = chains, n.adapt = 2500)
    update(reg.jags, n.iter=2500)
    check.reg[[PhenoTemp]] <- coda.samples(reg.jags, variable.names = paramUse, thin=10, n.iter=10000)
    
  }
  allnames <- c()
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
runLMERModels <- function(phenotype, orderdat, tryLam, normd=T,
                          DietRIX=T, PO=T, PODietRIX=T, Batch=F, Cage=F, Dam=F){
  check.mod_results <- list()
  if(DietRIX){
    if(PO){
      if(PODietRIX){
        if(Dam){
          if(Cage){
            if(Batch){
              indvariable <- "~ -1 + Diet + (1 | Cage) + (1 | DamID) + (1 | BreedingBatch) + (1 | RIX) + (PO + 0 | RIX) + (1 | Diet:RIX) + (PO + 0 | Diet:RIX)" 
              } else { 
              indvariable <- "~ -1 + Diet + (1 | DamID) + (1 | Cage) + (1 | RIX) + (PO + 0 | RIX) + (1 | Diet:RIX) + (PO + 0 | Diet:RIX)" 
              } 
            } else {
            indvariable <- "~ -1 + Diet + (1 | DamID) + (1 | RIX) + (PO + 0 | RIX) + (1 | Diet:RIX) + (PO + 0 | Diet:RIX)"
          }
        } else {   
          indvariable <- "~ -1 + Diet + (1 | RIX) + (PO + 0 | RIX) + (1 | Diet:RIX) + (PO + 0 | Diet:RIX)" 
        }
      } else {
        indvariable <- "~ -1 + Diet + (1 | RIX) + (PO + 0 | RIX) + (1 | Diet:RIX)" 
      }
    } else {
      indvariable <- "~ -1 + Diet + (1 | RIX) + (1 | Diet:RIX)" 
    } 
  } else {
    indvariable <- "~ -1 + Diet + (1 | RIX)"
  }
  
  for (i in 1:length(phenotype)){  
    pheno <- paste(phenotype[i])
    use <- which(is.na(orderdat[,pheno]) == F)
    matnut_use <- orderdat[use,]
    check.mod_results[[phenotype[i]]] <- fit.model.bc$fit(y.mat = matnut_use[,pheno], matnut_use, indvariable,
                                         transformParams = getMatnutTransformParams(normd=normd, tryLam = tryLam))
  }
  return(check.mod_results)
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

dontuse.lmer.getDecodedSummary <- function(lmerObject, fixed, rand1, rand2){
  lmerSummary <- data.frame(est = numeric(), Level=character(), Variable=character())
  for(i in 1:length(rand1)){
    this_rand <- paste0(rand1[i], rownames(lmerObject$ran))
    tempSum <- data.frame(cbind(lmerObject$ran[,rand1[i]], this_rand, rep(paste0(rand1[i],"RIX"))))
    colnames(tempSum) <- c("est","Level","Variable")
    lmerSummary <- rbind(lmerSummary, tempSum)
    
  }
  for(i in 1:length(rand2)){
    levels <- colnames(lmerObject$ran)[grep(rand2[i], colnames(lmerObject$ran))]
    levels2 <- gsub(rand2[i],"",levels)
    for (j in 1: length(levels)){
      this_rand <- paste0(levels2[j], rownames(lmerObject$ran))
      tempSum <- data.frame(cbind(lmerObject$ran[,levels[j]], this_rand, rep(paste0(rand2[i],"RIX"))))
      colnames(tempSum) <- c("est","Level","Variable")
      lmerSummary <- rbind(lmerSummary, tempSum)
    }
  }
  for(i in 1: length(fixed)){
    vec <- fixef(lmerObject$fit)
    levels <- names(vec)[grep(rand2[i], names(vec))]
    levels2 <- gsub(rand2[i],"",levels)
    for (j in 1: length(levels)){
      this_rand <- vec[j]
      tempSum <- data.frame(cbind(vec[j], levels2[j], rand2[i]))
      colnames(tempSum) <- c("est","Level","Variable")
      lmerSummary <- rbind(lmerSummary, tempSum)
    }
  }
  lmerSummary$est <- as.numeric(levels(lmerSummary$est))[lmerSummary$est]
  return(lmerSummary)
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


