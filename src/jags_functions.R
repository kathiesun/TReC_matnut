source("formulaWrapper.R")
library(rjags)

lmer.to.jags.form <- function(lmerform){
  jagsform <- gsub("-1","1",lmerform)
  jagsform <- gsub("[+] 0", "", jagsform)
  return(jagsform)
}


jags.getDecodedSummary <- function(mcmcObject, encoded, narrow=0.5, wide=0.95){
  listNames <- as.character(encoded$Level)
  #mu    <- c(colMeans(mcmcObject))
  mu <- c()
  for(i in 1:ncol(mcmcObject)){
    tempdens <- density(mcmcObject[,i])
    mu[i] <- tempdens$x[which.max(tempdens$y)]
  }
  med   <- apply(coda::HPDinterval(mcmcObject, prob=0.01), 1, median)
  hpd.wide    <- coda::HPDinterval(mcmcObject, prob=wide)
  hpd.narrow  <- coda::HPDinterval(mcmcObject, prob=narrow)
  effectdata <- data.frame(cbind(mu,med,hpd.narrow,hpd.wide))
  effectdata$Level <- rownames(effectdata)
  jagsSummary <- effectdata[intersect(effectdata$Level, listNames),]
  jagsSummary$Variable <- as.character(encoded$Variable[match(jagsSummary$Level, encoded$Level)])
  return(jagsSummary)
}

makeJagsModel <- function(lmerForm, dataf, encoded, y, qVar=NA,
                          a=0.001,b=0.001,c=1/10000^2, mu=0, 
                          sq=F, absVal=F, addS=c("off", "force", "fit"), prior="gamma"){
  addS <- addS[1]
  form <- formulaWrapper$parseCovariateString(paste(lmerForm, collapse=""))
  #MatchS=FALSE
  model <- paste0("model {\n sig2inv ~ d", prior, "(",a,",",b,")\n",
                  "b0 ~ dnorm(",mu,",", a,")\n")
  lik <- paste0("\nmu[i,1] <- b0 + ")
  n=length(y)
  data=list("N"=n, "y"=y)
  paramUse=c("sig2inv","b0")
  randEf <- list()
  k=1
  loopS <- c()
  z <- 1
  hasRandEf <- ifelse(length(form$ranef)>0, T, F)
  hasFixEf = F
  isQuant = T
  fixEf <- list()
  k=1
  if(length(form$fixef) > 0){
    for(i in 1:length(form$fixef)){
      if(form$fixef[[i]] != "1"){
        tempFix <- ifelse(length(grep("[*]", form$fixef[[i]])) >0, 
                          unlist(strsplit(form$fixef[[i]], "[*]"))[which(unlist(strsplit(form$fixef[[i]], "[*]")) %in% colnames(dataf))],
                          form$fixef[[i]])
        if(tempFix %in% colnames(dataf)){        
          hasFixEf = T
          isQuant = ifelse(tempFix %in% qVar, T, F)
          addTemp = ifelse(isQuant, "","temp")
          fixEf[[k]] <- list()
          fixEf[[k]]$name <- tempFix
          fixEf[[k]]$index <- as.numeric(dataf[[fixEf[[k]]$name]])
          fixEf[[k]]$length <- length(which(encoded$Variable == gsub("[0-9]+", "", fixEf[[k]]$name)))
          
          loop <- paste0("for (i in 1:n",fixEf[[k]]$name,"){ beta", addTemp, 
                         fixEf[[k]]$name,"[i] ~ dnorm(0, ", c, ") }")
          loopLik <- paste0("beta",fixEf[[k]]$name, "[ind", fixEf[[k]]$name, "[i]]") 
          tempLoop = ifelse(isQuant, "", 
                            paste0("beta",fixEf[[k]]$name," <- betatemp",fixEf[[k]]$name, "- mean(betatemp",fixEf[[k]]$name,")"))
          model <- paste(model, loop, tempLoop, sep="\n")
          
          if(length(grep("1$", fixEf[[k]]$name)) > 0){
            lik <- paste(lik, loopLik, sep="-")
          } else if (k==1) {
            lik <- paste(lik, loopLik, sep="")
          } else {
            lik <- paste(lik, loopLik, sep="+")
          }
          data[[paste0("n", fixEf[[k]]$name)]] = fixEf[[k]]$length
          data[[paste0("ind", fixEf[[k]]$name)]] = fixEf[[k]]$index
          paramUse <- c(paramUse, paste0("beta",fixEf[[k]]$name))
          k=k+1
        }
      }
    }
  }
  
  k=1
  if(hasRandEf){
    for(i in 1:length(form$ranef)){
      for(j in 1:length(form$ranef[[i]]$components)){
        num <- k+j-1
        randEf[[num]] <- list()
        addon <- ifelse(form$ranef[[i]]$components[[j]] %in% colnames(dataf), form$ranef[[i]]$components[[j]], "")
        addon <- gsub("^s", "", addon)
        randEf[[num]]$name <- paste0(addon, paste(form$ranef[[i]]$group, collapse = ""))
        randEf[[num]]$index <- as.numeric(dataf[[randEf[[num]]$name]])
        randEf[[num]]$length <- length(which(encoded$Variable == gsub("[0-9.]+", "", randEf[[num]]$name)))
        
        isQuant = ifelse(randEf[[num]]$name %in% qVar, T, F)
        addTemp = ifelse(isQuant, "","temp")
        
        if (!is.null(form$ranef[[i]]$components[[j]]) && form$ranef[[i]]$components[[j]] != "1"){
          randEf[[num]]$mult <- form$ranef[[i]]$components[[j]]
          if(form$ranef[[i]]$components[[j]] %in% colnames(dataf)){
            data[[form$ranef[[i]]$components[[j]]]] <- dataf[[form$ranef[[i]]$components[[j]]]] 
          }
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
          loop2 <- paste0("for (i in 1:n",randEf[[num]]$name,"){ \n beta",addTemp, 
                          randEf[[num]]$name,"[i] ~ dnorm(0, tau", randEf[[num]]$name, ") }\n")
          tempLoop = ifelse(isQuant, "", 
                            paste0("beta",randEf[[num]]$name," <- betatemp",randEf[[num]]$name, "- mean(betatemp",randEf[[num]]$name,")"))
        }
        if(form$ranef[[i]]$components[[j]] %in% colnames(dataf)){
          loopLik <- paste0(randEf[[num]]$mult, "[i]*beta", randEf[[num]]$name, "[ind", randEf[[num]]$name, "[i]]") 
        } else if(addS == "fit"){
          loopS[z] <- paste0("beta", randEf[[num]]$name, "[ind", randEf[[num]]$name, "[i]]")
          loopLik <- NA
          z=z+1
        } else{
          loopLik <- paste0(randEf[[num]]$mult, "beta", randEf[[num]]$name, "[ind", randEf[[num]]$name, "[i]]") 
        }
        
        model <- paste(model, loop1, loop2, tempLoop, sep="\n")
        if(!is.na(loopLik)){
          if(num==1 & !hasFixEf){
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
  }
  
  ######
  #if(length(grep("grand_mu", lmerForm))>0){
  #  loop <- paste("grand_mu ~ dunif(0, ",mu, ")")
  #  model <- paste(model, loop, sep="\n")
  #  lik <- paste(lik, "grand_mu", sep="+")
  #  paramUse <- c(paramUse, "grand_mu")
  #}
  ######

  
  lik0 <- "\nfor (i in 1:N){" 
  if(addS == "fit"){
    data[["PO"]] <- dataf$PO
    addStr <- c("PORIX_net[i] <- S[i] *betaPORIX[indPORIX[i]]",
                "PODietRIX_net[i] <- S[i] *betaPODietRIX[indPODietRIX[i]]")
    addThisStr <- ifelse(rep(length(grep("PO", lmerForm)) > 0, length(addStr)), addStr, gsub("PO", "", addStr))
    addThis <- c("S[i] <- ifelse(PO[i] > 0, betaS[indRIX[i]], betaS[indRIX[i]]*-1)", 
                 addThisStr)
    lik0 <- paste(lik0, paste0(addThis, collapse = "\n"), sep = "\n")
    lik <- paste(lik, "+ PORIX_net[i] + PODietRIX_net[i]")
    lik <- ifelse(length(grep("PO", lmerForm)) > 0, lik, "\nmu[i,1] <- RIX_net[i] + DietRIX_net[i]")
    #  lik <- paste(lik, "+ PO[i]*sign[i]*(", paste0(loopS, collapse = "+"), ")")
    paramUse <- c(paramUse, "betaS")
  }
  
  likFin <- paste(lik0, lik, "\n y[i]  ~ dnorm(mu[i,1], sig2inv) \n }\n}")
  loop3 <- ifelse(addS == "fit", paste0("\nthetaS ~ dbeta(1,1)",
                                        "\nfor (i in 1:nRIX){ \n", 
                                        "unscS[i] ~ dbern(thetaS) \n",
                                        "betaS[i] <- ifelse(unscS[i] < 0.5, -0.5, 0.5) }"), "")
  #loop3 <- ifelse(addS, paste0("\nfor(j in 1:nRIX){",
  #                             "\nfor (i in 1:N){", 
  #                             "\nsumS[i,j] <- ifelse(PO[i] < 0, 0, ", paste0(loopS, collapse = "+"), ") }",
  #                             "\nbetaS[j] <- ifelse(sum(sumS[,j]) < 0, -1, 1) }"), "")
  #                             "betaS[j] <- ifelse(betaPORIX[j] < 0, -1, 1) }"), "")
  
  modelFin <- paste(model, loop3, likFin)
  return(list(modelFin=modelFin, data=data, paramUse=paramUse))
}

runJagsModels_cov <- function(datalist, chains=2, formula=NULL, testLMER=NULL,
                              n.adapt=2500, n.iter=20000, thin=10, qVar = NA, 
                              encoded=NA, phenotype=NA, mu=0, sq=F, addS=c("off", "force", "fit")){
  addS <- addS[1]
  fitJags <- list()
  if(class(datalist) == "list"){
    dataf <- datalist$df
    encoded <- datalist$encoded
    if(any(is.na(phenotype))){
      phenotype <- datalist$ptypes
    } 
  } else {
    dataf <- datalist
  }

  for (i in 1:length(phenotype)){
    if(is.null(testLMER)){
      datause <- dataf[!is.na(dataf[,phenotype[i]]),]
      pheno = as.vector(dataf[,phenotype[i]])
      lmerForm=formula
    } else if (class(testLMER$lmerobj$fits[[i]]) == "character"){
      datause <- dataf[!is.na(dataf[,phenotype[i]]),]
      pheno = as.vector(dataf[,phenotype[i]])
      lmerForm <- testLMER$JAGSformula
    } else {
      datause <- dataf[!is.na(dataf[,phenotype[i]]),]
      pheno <- testLMER$lmerobj$y.transform[[i]]
      lmerForm <- testLMER$JAGSformula
    }
    N <- length(pheno)
    y <- as.vector(pheno[!is.na(pheno)])

    madeModel <- makeJagsModel(lmerForm,dataf=datause,encoded,y,
                               mu=mu, sq=sq, addS=addS, qVar = NA)
      
    reg.jags <- jags.model(textConnection(madeModel$modelFin), data=madeModel$data, n.chains = chains, n.adapt = n.adapt)
    update(reg.jags, n.iter=n.adapt)
    fitJags[[phenotype[i]]]$fit <- coda.samples(reg.jags, variable.names = madeModel$paramUse, thin=thin, n.iter=n.iter)
    fitJags[[phenotype[i]]]$madeModel <- madeModel
    allnames <- c()
    paramUse <- gsub("beta","", unique(unlist(
      lapply(strsplit(varnames(fitJags[[phenotype[i]]]$fit),'[[]'), function(x) x[[1]]))))
    wantParam <- intersect(paramUse, c(paste(unique(encoded$Variable)), "grand_mu"))
    #if(madeModel$MatchS){
    #  wantParam <- c(wantParam, "s")
    #  allnames <- c(allnames, "s[1]", "s[2]")
    #} 
    for(j in 1:length(wantParam)){
      allnames <- c(allnames, as.character(encoded$Level[which(encoded$Variable == wantParam[j])])) 
    }
    varnames(fitJags[[phenotype[i]]]$fit)[grep("[[]", varnames(fitJags[[phenotype[i]]]$fit))] <- allnames
  }
  return(fitJags)
}

