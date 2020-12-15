makeJagsModel_trynaFitS <- function(lmerForm, dataf, encoded, y, a=0.001,b=0.001,c=1/10000^2, mu=NA,
sq=F, addS=F, prior="gamma"){
    form <- formulaWrapper$parseCovariateString(lmerForm)
    #form$fixef[[which(duplicated(form$fixef))]] <- NULL
    MatchS=FALSE
    model <- paste0("model {\n sig2inv ~ d", prior, "(",a,",",b,")")
    strg <- ifelse(addS, "\nsign[i] <- ifelse(betaPORIX[indPORIX[i]]+betaPODietRIX[indPODietRIX[i]] < 0, -1, 1)", "")
    lik <- paste0("\n for (i in 1:N){", strg, "\nmu[i,1] <- ")
    n=length(y)
    data=list("N"=n, "y"=y)
    paramUse=c("sig2inv")
    randEf <- list()
    k=1
    loopS <- c()
    z <- 1
    for(i in 1:length(form$ranef)){
        for(j in 1:length(form$ranef[[i]]$components)){
            num <- k+j-1
            randEf[[num]] <- list()
            addon <- ifelse(form$ranef[[i]]$components[[j]] %in% colnames(dataf), form$ranef[[i]]$components[[j]], "")
            addon <- gsub("^s", "", addon)
            randEf[[num]]$name <- paste0(addon, paste(form$ranef[[i]]$group, collapse = ""))
            randEf[[num]]$index <- as.numeric(dataf[[randEf[[num]]$name]])
            randEf[[num]]$length <- length(which(encoded$Variable == randEf[[num]]$name))
            if (!is.null(form$ranef[[i]]$components[[j]]) && form$ranef[[i]]$components[[j]] %in% colnames(dataf)){
                randEf[[num]]$mult <- form$ranef[[i]]$components[[j]]
                #if(!addS){
                data[[form$ranef[[i]]$components[[j]]]] <- dataf[[form$ranef[[i]]$components[[j]]]]
                #}
            } else {
                randEf[[num]]$mult <- ""
                #rep(1,length(randEf[[num]]$index))
            }
            loop1 <- paste0("tau", randEf[[num]]$name," ~ d", prior, "(",a,",", b,")")
            if (randEf[[num]]$mult == "PO" && sq){
                loop2 <- paste0("for (i in 1:n",randEf[[num]]$name,"){ \n",
                "beta", randEf[[num]]$name, "[i] <- pow(unsc_", randEf[[num]]$name, "[i]*tau", randEf[[num]]$name,",0.5)",
                #"\n unsc_", <- pow(beta", randEf[[num]]$name, "[i], 2)/tau", randEf[[num]]$name,
                "\n unsc_", randEf[[num]]$name, "[i] ~ dchisq(1) }")
                #if(addS){ S <- "betaS[indRIX[i]]"}
            } else{
                loop2 <- paste0("for (i in 1:n",randEf[[num]]$name,"){ \n beta",
                randEf[[num]]$name,"[i] ~ dnorm(0, tau", randEf[[num]]$name, ") }")
            }
            if(randEf[[num]]$mult == ""){
                loopLik <- paste0("beta", randEf[[num]]$name, "[ind", randEf[[num]]$name, "[i]]")
            } else if(addS){
                loopS[z] <- paste0("beta", randEf[[num]]$name, "[ind", randEf[[num]]$name, "[i]]")
                loopLik <- NA
                z=z+1
            } else{
                loopLik <- paste0(randEf[[num]]$mult, "[i]*beta", randEf[[num]]$name, "[ind", randEf[[num]]$name, "[i]]")
            }
            
            model <- paste(model, loop1, loop2, sep="\n")
            if(!is.na(loopLik)){
                if(num == 1){
                    lik <- paste(lik, loopLik, sep="")
                } else {
                    lik <- paste(lik, loopLik, sep="+")
                }
            }
            data[[paste0("n", randEf[[num]]$name)]] = randEf[[num]]$length
            data[[paste0("ind", randEf[[num]]$name)]] = randEf[[num]]$index
            paramUse <- c(paramUse, paste0("beta", randEf[[num]]$name), paste0("tau",randEf[[num]]$name))
        }
        k=k+j
    }
    if(length(grep("grand_mu", lmerForm))>0){
        loop <- paste("grand_mu ~ dunif(0, ",mu, ")")
        model <- paste(model, loop, sep="\n")
        lik <- paste(lik, "grand_mu", sep="+")
        paramUse <- c(paramUse, "grand_mu")
    }
    fixEf <- list()
    k=1
    for(i in 1:length(form$fixef)){
        if(form$fixef[[i]] != "1"){
            if(gsub("[*]","",form$fixef[[i]]) %in% colnames(dataf)){
                fixEf[[k]] <- list()
                fixEf[[k]]$name <- paste(gsub("[*]","",form$fixef[[i]]), collapse = "")
                fixEf[[k]]$index <- as.numeric(dataf[[fixEf[[k]]$name]])
                fixEf[[k]]$length <- length(which(encoded$Variable == fixEf[[k]]$name))
                
                loop <- paste0("for (i in 1:n",fixEf[[k]]$name,"){ beta", fixEf[[k]]$name,"[i] ~ dnorm(0, ", c, ") }")
                loopLik <- paste0("beta",fixEf[[k]]$name, "[ind", fixEf[[k]]$name, "[i]]")
                model <- paste(model, loop, sep="\n")
                lik <- paste(lik, loopLik, sep="+")
                data[[paste0("n", fixEf[[k]]$name)]] = fixEf[[k]]$length
                data[[paste0("ind", fixEf[[k]]$name)]] = fixEf[[k]]$index
                paramUse <- c(paramUse, paste0("beta",fixEf[[k]]$name))
                k=k+1
            }
        }
    }
    if(addS){
        # lik <- paste(lik, "+ betaS[indRIX[i]]*abs(", paste0(loopS, collapse = "+"), ")")
        lik <- paste(lik, "+ PO[i]*sign[i]*(", paste0(loopS, collapse = "+"), ")")
        paramUse <- c(paramUse)
    }
    likFin <- paste(lik, "\n y[i]  ~ dnorm(mu[i,1], sig2inv) \n }\n}")
    #loop3 <- ifelse(addS, paste0("\nthetaS ~ dbeta(1,1)",
    #                             "\nfor (i in 1:nRIX){ \n",
    #                             "unscS[i] ~ dbern(thetaS) \n",
    #                             "betaS[i] <- ifelse(unscS[i] < 0.5, -1, 1) }"), "")
    #loop3 <- ifelse(addS, paste0("\nfor(j in 1:nRIX){",
    #                             "\nfor (i in 1:N){", 
    #                             "\nsumS[i,j] <- ifelse(PO[i] < 0, 0, ", paste0(loopS, collapse = "+"), ") }",
    #                             "\nbetaS[j] <- ifelse(sum(sumS[,j]) < 0, -1, 1) }"), "")
    #"betaS[j] <- ifelse(betaPORIX[j] < 0, -1, 1) }"), "")
    
    modelFin <- paste(model, likFin)
    return(list(modelFin=modelFin, data=data, paramUse=paramUse))
}



prediction_old <- function(mcmcOb, ptypes, encoded, dietOnly=F, Match=F){
  for (j in 1:length(ptypes)){
    ribbon9_md <- list()
    ribbon9 <- list()
    ribbon4 <- list()
    browser()
    for (i in 1:length(mcmcOb)){
      all.dist <- cbind(mcmcOb[[i]], rep(0,nrow(mcmcOb[[i]])), rep(0,nrow(mcmcOb[[i]]))) 
      colnames(all.dist)[c(length(colnames(all.dist))-1, length(colnames(all.dist)))] <- c("emp", "POemp")
      
      dietrix <- matrix(c(as.character(encoded$Level[which(encoded$Variable == "DietRIX")]),"emp", "emp"), 
                        ncol=9, nrow=4)
      dietrix[,9] <- c("emp","ME10", "PD9", "emp")
      nrix = length(which(encoded$Variable == "RIX"))
      ndiet = length(which(encoded$Variable == "Diet"))
      
      predPO <- list()
      predPOneg <- list()
      predPOpos <- list()
      predDiet <- list()
      for(i in 1:nrow(all.dist)){
        for (rixn in 1:nrix){
          for (dietn in 1:ndiet){
            
            rixlabs <- as.character(encoded$Level[which(encoded$Variable == "RIX")])
            dietlabs <- as.character(encoded$Level[which(encoded$Variable == "Diet")])
            
            temprix <- rixlabs[rixn]
            tempdiet <- dietlabs[dietn]
            tempinter <- dietrix[dietn,rixn]
            prLab <- paste0(tempdiet,temprix)
            
            if (dietOnly){
              predDiet[[prLab]][i] <- all.dist[i,tempdiet] + all.dist[i,dietrix[dietn,rixn]]
            } else if (Match){
              predPO[[prLab]][i] <- all.dist[i,temprix] + all.dist[i,dietrix[dietn,rixn]]
            } else{
              predPOneg[[paste0(prLab, ".neg")]][i] <- all.dist[i,temprix] + all.dist[i,tempdiet] + 
                all.dist[i,dietrix[dietn,rixn]] + all.dist[i,paste0("PO", temprix)]*(-0.5) + 
                all.dist[i,paste0("PO", tempinter)]*(-0.5)
              
              predPOpos[[paste0(prLab, ".pos")]][i] <- all.dist[i,temprix] + all.dist[i,tempdiet] + 
                all.dist[i,dietrix[dietn,rixn]] + all.dist[i,paste0("PO", temprix)]*(0.5) + 
                all.dist[i,paste0("PO", tempinter)]*(0.5)
            }
          }
        }
      }
    }
    
    if (dietOnly){
      predDiet <- do.call(rbind, predDiet)
      flydata <- as.data.frame(predDiet[-which(rownames(predDiet) %in% c("STD10","VDD10")),])
      flydata <- as.data.frame(allpredDiet)
    } else if (Match){
      predPO_temp <- do.call(rbind, predPO)
      flydata <- as.data.frame(predPO_temp[-which(rownames(predPO_temp) %in% c("STD10","VDD10")),])
    } else { 
      predPOneg <- do.call(rbind, predPOneg)
      predPOpos <- do.call(rbind, predPOpos)
      allpredPOneg <- as.data.frame(t(predPOneg[-which(rownames(predPOneg) %in% c("STD10.neg","VDD10.neg")),]))
      allpredPOpos <- as.data.frame(t(predPOpos[-which(rownames(predPOpos) %in% c("STD10.pos","VDD10.pos")),]))
      flydata <- cbind(allpredPOneg, allpredPOpos) 
    }
  } 
  mcmcOb <- as.mcmc(flydata)
  return(mcmcOb=mcmcOb)
}


OLD_runJagsModels_cov <- function(phenotype, dataf, allparam, encoded, testLMER=NULL,
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
    
    
    rix <- as.numeric(matnut_use$RIX)
    nrix <- length(which(encoded$Variable == "RIX"))
    dietrix <- as.numeric(matnut_use$DietRIX)
    ndietrix <- length(which(encoded$Variable == "DietRIX"))
    diet <- as.numeric(matnut_use$Diet)
    
    if(Match==F){
      ndiet <- length(which(encoded$Variable == "Diet")) 
      poe <- as.numeric(matnut_use$PO)
      be <- as.numeric(matnut_use$BehaviorBatch)
      dam <- as.numeric(matnut_use$DamID)
      sire <- as.numeric(matnut_use$SireID)
      #cage <- as.numeric(matnut_use$Cage)
      nbatch <- length(which(encoded$Variable == "BehaviorBatch")) 
      ndam <- length(which(encoded$Variable == "DamID")) 
      nsir <- length(which(encoded$Variable == "SireID")) 
      #ncage <- length(which(encoded$Variable == "Cage"))
    }
    
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
  
  test <- runLMERModels_cov(phenotype = ptypes[i], dataf=orderdat, tryLam=c(-1, 0, .25, .33, .5, 1, 2, 3))
  testJags <- runJagsModels_cov(ptypes[i], orderdat, encoded, testLMER = test)
  
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





###############
# unused code #
###############
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
    vec <- fixef(lmerObject)
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



# below code sorts first - DON'T USE
###########
sorted <- list()
for (i in 1:ncol(all.mcmc)){
  name <- varnames(all.mcmc)[i]
  sorted[[name]] <- sort(as.matrix(all.mcmc[,i]))
}

for (i in 1:nrow(rix.contr)){
  compare <- rixlabs[which(rix.contr[i,]!=0)]
  lab <- paste0(compare[1], compare[2])
  numSamps <- length(sorted[[compare[1]]])
  for (j in 1:numSamps){
    diff[[lab]][j] <- sorted[[compare[1]]][j] - sorted[[compare[2]]][j]
  }
  #plot(diff[[lab]])
  over0 <- length(which(diff[[lab]]>0))
  p_val[i] <- min(over0, numSamps - over0) / numSamps
  names(p_val)[i] <- lab
}
for (i in 1:nrow(diet.contr)){
  compare <- dietlabs[which(diet.contr[i,]!=0)]
  lab <- paste0(compare[1], compare[2])
  numSamps <- length(test.mcmc[,compare[1]])
  for (j in 1:numSamps){
    dietdiff[[lab]][j] <- test.mcmc[,compare[1]][j] - test.mcmc[,compare[2]][j]
  }
  plot(density(dietdiff[[lab]]), main=paste0("Plot of delta",lab))
  over0 <- length(which(dietdiff[[lab]]>0))
  diet_p_val[i] <- min(over0, numSamps - over0) / numSamps
  names(diet_p_val)[i] <- lab
}

dietcol <- match(dietlabs,varnames(all.reg))
rixcol <- match(rixlabs,varnames(all.reg))

rixef <- apply(all.mcmc[,rixcol], 2, median)
names(rixef) <- rixlabs
dietef <- apply(all.mcmc[,dietcol], 2, median)
names(dietef) <- dietlabs

tryplots$data$dietef <- dietef[match(tryplots$data$inter1, names(dietef))]
tryplots$data$rixef <- rixef[match(tryplots$data$inter2, names(rixef))]

add2plots <- tryplots + 
  geom_line(data=tryplots$data, aes(x=names, y=rixef, group=inter2), col="black") + 
  geom_point(data=tryplots$data, aes(x=names, y=dietef, group=inter2), col="darkgray", shape="|", size=3) +
  labs(title=paste("PO-by-RIX-by-diet effects of",phenotype[i])) 
print(add2plots)
}  
### POE effects ###
tryplots2 <- plot.hpd(all.reg, wanted=poe, order=2)
#print(tryplots)
add2plots2 <- tryplots2 + 
  geom_point(data=tryplots2$data, aes(x=names, y=rixef), col="darkgray", shape="|", size=3) +
  labs(title=paste("RIX-by-POE effects of",phenotype[i])) 
print(add2plots2)


###########################
# Didn't use this portion #
###########################

seed=245234
jags.init <- list(.RNG.seed=seed, .RNG.name="base::Mersenne-Twister",
                  rixef=rep(0,nrix), perixef=rep(0,nrix), beta=rep(0,ndiet), dietrixef=rep(0,nrix*ndiet))
#,tauInv2=0.01, sigInv2=0.01, tauDRInv2=0.01, tauPEInv2=0.01)
#,tau=100, sig=100) # if using unif prior on tau and sig
varparam <- c("sigInv2", "tauInv2","tauDRInv2", "tauPEInv2") 
effectparam <- c("rixef", "beta","dietrixef", "perixef") 

####################
# Use HPDplot code #
####################

num.plots <- length(phenotype)
my.plots <- vector(num.plots, mode='list')

for (i in 1:num.plots){
  
  all.reg <- list.reg[[i]]
  all.mcmc <- mcmc.stack(all.reg)
  all.mcmc <- cbind(all.mcmc, rep(0, nrow(all.mcmc)))
  
  
  #layout(matrix(c(1,1,2), ncol=1, nrow=3))  
  plot.hpd(all.mcmc, wanted=c(dietlabs,"",rixlabs,"",PO.R,"",rixFirst, "", poRIXFirst), cex=0.55)
  title(paste("HPD distributions for", phenotype[i]), outer = TRUE, line = -2)
  
  my.plots[[i]] <- recordPlot()
}

graphics.off()

pdf('hpd_all_singlepg.pdf', onefile=TRUE, width = 8.5, height = 11)
for (my.plot in my.plots) {
  replayPlot(my.plot)
}
graphics.off() 


#############
# Plot loop #
#############
num.plots <- length(phenotype)
my.plots <- vector(num.plots, mode='list')

for (i in 1:length(phenotype)){
  
  all.reg <- list.reg[[i]]
  all.mcmc <- mcmc.stack(all.reg)
  #HPDinterval(all.mcmc)
  
  ### Diet effects  ###
  tryplots <- plot.mine.hpd(all.reg, wanted=dietrix, order=2, addline=T)
  tryplots2 <- plot.mine.hpd(all.reg, wanted=podietrix, order=2, remove = "PO.", addline=T)
  tryplots3 <- plot.mine.hpd(all.reg, wanted=c(dietlabs, rixlabs, PO.R), order=3, remove = "PO.")
  
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(2, 3, heights = unit(c(0.5, 6),"null"))))  
  grid.text(paste("Estimates for",phenotype[i]), vp = viewport(layout.pos.row = 1, layout.pos.col = 1:3))
  print(tryplots, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))         
  print(tryplots2, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
  print(tryplots3, vp = viewport(layout.pos.row = 2, layout.pos.col = 3))
  
  my.plots[[i]] <- recordPlot()
  
} 

graphics.off()

pdf('hpd_all.pdf', onefile=TRUE, width = 11, height = 8.5)
for (my.plot in my.plots) {
  replayPlot(my.plot)
}
graphics.off() 




##################
# Some functions #
##################
rixtest <- function(all.mcmc, allnames, rix.contr, rixlabs, sig.val=0.05, sig=F)
{
  p_val <- list()
  signif <- matrix(NA, nrow=nrow(allnames), ncol=nrow(rix.contr))
  title <- c()
  rixdiff <- list()
  for (k in 1:nrow(allnames)){
    rixnames <- allnames[k,]
    test.mcmc <- all.mcmc[,rixnames]
    p_val[[k]] <- vector()
    for (i in 1:nrow(rix.contr)){
      compare <- rixnames[which(rix.contr[i,]!=0)]
      if (k==1){
        titlePieces <- rixlabs[which(rix.contr[i,]!=0)]
        title[i] <- paste0(titlePieces[1], titlePieces[2])  }
      lab <- paste0(compare[1], compare[2])
      numSamps <- length(test.mcmc[,compare[1]])
      for (j in 1:numSamps){
        rixdiff[[lab]][j] <- test.mcmc[,compare[1]][j] - test.mcmc[,compare[2]][j]
      }
      if (sig==T){
        over0 <- length(which(rixdiff[[lab]]>0))
        p_val[[k]][i] <- min(over0, numSamps - over0) / numSamps
        names(p_val[[k]])[i]<- lab
        if (p_val[[k]][i] < sig.val/(nrow(rix.contr)*10)) {
          signif[k,i] <- "**"
        } else if (p_val[[k]][i] < sig.val/nrow(rix.contr)){
          signif[k,i] <- "*"
        } else {signif[k,i] <- "-"}
      }
    }
  }
  colnames(signif) <- title
  rownames(signif) <- rownames(allnames)
  if (sig){
    return(list(signif=signif, p_val=p_val, title=title, rixdiff=rixdiff))
  } else {
    return(list(rixdiff=rixdiff, title=title))
  }
}

factorize <- function(data, variables){
  for (i in 1:length(variables)){
    var <- variables[i]
    newvar <- paste0(var, "_f")
    data[,newvar] <- factor(data[,var])
  }
  return(data=data)
}

diction <- function(data, variables){
  dictionary <- list()
  for (i in 1:length(variables)){
    var <- variables[i]
    newvar <- paste0(var, "_f")
    dictionary[[var]]$name <- levels(data[,newvar])
    dictionary[[var]]$index <- as.numeric(unique(data[,newvar]))
    dictionary[[var]]$length <- length(levels(data[,newvar]))
  }
  return(dictionary=dictionary) 
}

labelnames <- function(dict) {
  dietnames <- c()
  rixnames <- c()
  dietrixnames <- c()
  poenames <- c()
  podietrixnames <- c()
  allnames <- c()
  for (i in 1:dict$Diet$length){
    dietnames[i] <- dict$Diet$name[i]
    for (j in 1:dict$RIX$length){
      dietrixnames[((i-1)*dict$RIX$length)+j] <- paste0(dict$Diet$name[i],".RIX",
                                                        dict$RIX$name[dict$RIX$index][j])
      podietrixnames[((i-1)*dict$RIX$length)+j] <- paste0("PO.",dict$Diet$name[i],".RIX",
                                                          dict$RIX$name[dict$RIX$index][j])
    }
  }
  for (i in 1:dict$RIX$length){
    rixnames[i] <- paste0("RIX",dict$RIX$name[dict$RIX$index][i])
    poenames[i] <- paste0("PO.RIX",dict$RIX$name[dict$RIX$index][i])
  }
  allnames <- c(dietnames, dietrixnames, poenames, podietrixnames, rixnames)
  return(allnames)
} 


###################
# Trace/den plots #
###################

all.reg <- list.reg[[i]]
for (i in 1:length(allparam)){
  if(ncol(all.reg) > 1){
    indiv.reg <- coda.samples(reg.jags, variable.names = allparam[i], thin=10, n.iter=10000)
    for (j in 1:ncol(all.reg[[1]])){
      plot(indiv.reg[,j], main=paste("Plot of",allparam[i],j))
    }
  } else {
    plot(indiv.reg)
  }
}



###############
# Get summary #
###############
allSummaries <- function(phenotype, orderdat, allparam, encoded, tryLam=1){
  #phenotype=ptypes[i]
  testLMER <- runLMERModels(phenotype=phenotype, orderdat, tryLam = tryLam, Batch = Batch)
  lmerSum <- lmer.getDecodedSummary(lmerObject=testLMER[[phenotype]]$phen_1, fixed="Diet", rand=T)
  transf <- c("lambda" = testLMER[[phenotype]]$phen_1$lambda, 
              "pval" = testLMER[[phenotype]]$phen_1$best.pval)
  
  testJags <- runJagsModels(phenotype=phenotype, orderdat, allparam, encoded, testLMER, Batch = Batch)
  all.reg <- testJags
  all.mcmc <- mcmc.stack(all.reg)
  jagsSum <- jags.getDecodedSummary(all.mcmc, encoded)
  
  allSummaries <- compareSummaries(lmerSum, jagsSum)
  if(Batch){
    getSummaries <- allSummaries[-which(allSummaries$Variable %in% c("BreedingBatch","DamID")),]
  } else {
    getSummaries <- allSummaries
  }
  
  #diet.test <- contrast.test(dfObject=as.data.frame(all.mcmc), encoded, variable="Diet", byVar="RIX", contrast=diet.contr)
  null.test <- contrast.null(dfObject=as.data.frame(all.mcmc), levels=encoded$Level)
  #rix.test <- contrast.test(dfObject=all.mcmc, encoded, variable="RIX", byVar="Diet", contrast=rix.contr)
  
  ##########
  # Tables #
  ##########
  
  jagsSum$Level <- factor(jagsSum$Level, levels=jagsSum$Level)
  null_p <- data.frame(Level = names(null.test$pval),
                       pval = null.test$pval, 
                       signif = null.test$signif)
  null_p$Level <- factor(null_p$Level, levels=jagsSum$Level)
  null_p <- null_p[order(null_p$Level),]
  mergeTable <- merge(jagsSum, null_p, by="Level")
  mergeTable$Level <- factor(mergeTable$Level, levels=encoded$Level)
  mergeTable$Variable <- factor(mergeTable$Variable, levels=unique(encoded$Variable))
  
  summary_table[[phenotype]] <- cbind(phenotype=rep(phenotype,nrow(mergeTable)), mergeTable[order(mergeTable$Variable, mergeTable$Level),])
  jagsLmer_compare <- cbind(phenotype=rep(phenotype,nrow(getSummaries)), getSummaries) 
  
  ####
  
  #dietTab <- contrasts.getDecodedSummary("Diet")
  #dietDiet_table <- cbind(phenotype=rep(phenotype,nrow(dietTab)), dietTab)
  #rixTab <- contrasts.getDecodedSummary("RIX")
  #rixRix_table <- cbind(phenotype=rep(phenotype,nrow(rixTab)), rixTab)
  
  
  #########
  # Plots #
  #########
  
  plotCat <- combinedCat(lmerSum, jagsSum, getSummaries, phenotype)
  return(list(transformation=transf, mcmcObject=all.mcmc, plot=plotCat, allSummaries=allSummaries,
              null.test=null.test, #diet.test=diet.test, rix.test=rix.test, 
              getSummaries=getSummaries, jagsSummary=jagsSum, lmerSummary=lmerSum))
}


testPhen <- list()
for(i in 1:length(ptypes)){
  testPhen[[ptypes[i]]] <- allSummaries(phenotype = ptypes[i], orderdat, allparam, encoded, 
                                        tryLam=c(-1, 0, .25, .33, .5, .1, 2, 3))
}

psychPhen <- readRDS("./matnut_outputs/psychPhenResults.rds")


####################
# Compare rand efs #
####################

allSummaries_rand <- function(phenotype, orderdat, allparam, encoded, tryLam=1, 
                              Batch=F, Dam=F, Cage=F){
  phenotype=ptypes[i]
  testLMER <- runLMERModels(phenotype=phenotype, orderdat, tryLam = tryLam, 
                            Batch = Batch, Dam = Dam, Cage = Cage)
  lmerSum <- lmer.getDecodedSummary(lmerObject=testLMER[[phenotype]]$phen_1, fixed="Diet", rand=T)
  transf <- c("lambda" = testLMER[[phenotype]]$phen_1$lambda, 
              "pval" = testLMER[[phenotype]]$phen_1$best.pval)
  
  testJags <- runJagsModels(phenotype=phenotype, orderdat, allparam, encoded, testLMER,
                            Batch = Batch, Dam = Dam, Cage = Cage)
  all.reg <- testJags[[phenotype]]
  all.mcmc <- mcmc.stack(all.reg)
  jagsSum <- jags.getDecodedSummary(all.mcmc, encoded)
  
  allSummaries <- compareSummaries(lmerSum, jagsSum)
  if(Dam){
    if(Cage){
      if(Batch){
        getSummaries <- allSummaries[-which(allSummaries$Variable %in% c("BreedingBatch", "DamID","Cage")),]
      } else {
        getSummaries <- allSummaries[-which(allSummaries$Variable %in% c("Cage", "DamID")),]
      } 
    } else {
      getSummaries <- allSummaries[-which(allSummaries$Variable %in% "DamID"),]
    }
  } else {
    getSummaries <- allSummaries
  }
  
  null.test <- contrast.null(dfObject=as.data.frame(all.mcmc), levels=unique(allSummaries$Level))
  diet.test <- contrast.test(dfObject=as.data.frame(all.mcmc), encoded, variable="Diet", byVar="RIX", contrast=diet.contr)
  rix.test <- contrast.test(dfObject=as.data.frame(all.mcmc), encoded, variable="RIX", byVar="Diet", contrast=rix.contr)
  
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
              summary_table=summary_table, JL_compare=jagsLmer_compare))
}


### Run function ###
testPhen <- list()
for(i in 1:length(ptypes)){
  testPhen[[ptypes[i]]] <- allSummaries_rand(phenotype = ptypes[i], orderdat, allparam, encoded, Batch=TRUE, Dam=TRUE, Cage=FALSE, 
                                             tryLam=c(-1, 0, .25, .33, .5, .1, 2, 3))
}

### Make summary table ###
summaryAll <- data.frame()
jagsLmer_compare <- data.frame()

for(i in 1:length(ptypes)){
  summaryAll <- rbind(summaryAll, testPhen[[ptypes[i]]]$summary_table
                      [-which(testPhen[[ptypes[i]]]$summary_table$Variable %in% c("BreedingBatch","DamID")),])
  jagsLmer_compare <- rbind(jagsLmer_compare, testPhen[[ptypes[i]]]$JL_compare)
  
}
write.csv(summaryAll, file.path("./","matnut_outputs/",'allPheno_summaryTable.csv'), row.names = F)
write.csv(jagsLmer_compare, file.path("./","matnut_outputs/",'allPheno_JL_compare.csv'), row.names = F)




###################################
# Model with different parameters #
###################################

runJagsModels <- function(phenotype, orderdat, allparam, encoded, testLMER=NULL,
                          DietRIX=T, PO=T, PODietRIX=T, Batch=T, Dam=T, Cage=T){
  
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


###############
# Get summary #
###############

allSummaries_rand <- function(phenotype, orderdat, allparam, encoded, tryLam=1)
{
  phenotype=ptypes[i]
  testLMER <- runLMERModels(phenotype=phenotype, orderdat, tryLam = tryLam)
  lmerSum <- lmer.getDecodedSummary(lmerObject=testLMER[[phenotype]]$phen_1, fixed="Diet", rand=T)
  transf <- c("lambda" = testLMER[[phenotype]]$phen_1$lambda, 
              "pval" = testLMER[[phenotype]]$phen_1$best.pval)
  
  testJags <- runJagsModels(phenotype=phenotype, orderdat, allparam, encoded, testLMER,
                            Batch = Batch, Dam = Dam, Cage = Cage)
  all.reg <- testJags[[phenotype]]
  all.mcmc <- mcmc.stack(all.reg)
  jagsSum <- jags.getDecodedSummary(all.mcmc, encoded)
  
  allSummaries <- compareSummaries(lmerSum, jagsSum)
  if(Dam){
    if(Cage){
      if(Batch){
        getSummaries <- allSummaries[-which(allSummaries$Variable %in% c("BreedingBatch", "DamID","Cage")),]
      } else {
        getSummaries <- allSummaries[-which(allSummaries$Variable %in% c("Cage", "DamID")),]
      } 
    } else {
      getSummaries <- allSummaries[-which(allSummaries$Variable %in% "DamID"),]
    }
  } else {
    getSummaries <- allSummaries
  }
  
  null.test <- contrast.null(dfObject=as.data.frame(all.mcmc), levels=unique(allSummaries$Level))
  diet.test <- contrast.test(dfObject=as.data.frame(all.mcmc), encoded, variable="Diet", byVar="RIX", contrast=diet.contr)
  rix.test <- contrast.test(dfObject=as.data.frame(all.mcmc), encoded, variable="RIX", byVar="Diet", contrast=rix.contr)
  
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
              summary_table=summary_table, JL_compare=jagsLmer_compare))
}

