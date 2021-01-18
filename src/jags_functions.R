source("formulaWrapper.R")
library(rjags)
library(coda)

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


CC_pairs <- c("CC001/011", "CC041/051", "CC017/004", "CC023/047", "CC026/006", "CC014/003", "CC035/062", "CC032/042", "CC005/040")

###################################

run_jags_regress <- function(data_kmers, 
                             niter=1000,
                             n.thin = 1, 
                             seg_regions=NULL,
                             reg=NULL,
                             save_dir=NULL,
                             STZ=T, use_gene=F,#no_gene=F, 
                             no_theta=F, alpha=NULL,
                             cond=NULL,
                             stan=F, stanMod=NULL){
  reg <- list()
  for(i in sort(unique(data_kmers$RRIX))){
    testPups <- as.vector(t(data_kmers %>% ungroup() %>%
                              filter(RRIX == i) %>%            #, !Pup.ID %in% problemPups
                              dplyr::select("Pup.ID") %>% distinct()))
    testData <- data_kmers %>% filter(Pup.ID %in% testPups)
    chr = unique(testData$seq.Chromosome)
    terms <- c("dir", "Diet", "DietRIX","RRIX", "PODIETRIX")
    if(!is.null(cond)){
      testData$PORIX = testData$RIX
      testData$PODIETRIX = paste(testData$DietRIX, testData$dir, sep="_")
      testData$PODIETRIX = factor(testData$PODIETRIX, 
                                  levels=apply(data.frame(expand.grid(sort(unique(testData$Diet)), sort(unique(testData$RRIX)), sort(unique(testData$dir)))),
                                               1, function(x) paste(x,collapse="_")))
      getCond = table((testData %>% select(PODIETRIX, Pup.ID) %>% distinct())$PODIETRIX)
      
      terms <- c("dir", "Diet", "DietRIX","RRIX", "PODIETRIX")
      cond = "PODIETRIX"
      if(any(getCond == 1)){
        terms = cond = "RRIX"
      }
    }
    
    #colnames(testData) <- toupper(colnames(testData))
    if(!is.null(seg_regions) & chr!= "X"){
      use_region = seg_regions[[paste(i)]] %>% filter(Chr == chr)
      rem = c()
      for(j in 1:nrow(use_region)){
        tmp = use_region[j,]
        rem = c(rem, which(testData$seq.Position > (tmp$start) & testData$seq.Position < (tmp$end)))
      }
      if(length(rem)>0) testData = testData[-rem,]
    }
    #terms <- c("dir", "Diet", "DietRIX","RRIX", "PODIETRIX")
    if(!stan){
      reg[[paste0("rix_",i)]] <- jags.genes.run(data=testData, mu_g0=0.5, niter=niter, n.thin=n.thin, 
                                                nchains=2, terms=terms, STZ=STZ, use_gene=use_gene, 
                                                no_theta=no_theta, alpha=alpha, quant=NULL, tau=NULL, 
                                                cond=cond, #no_gene=no_gene, 
                                                TIMBR=F, C_diet_4=C_diet_4, C_diet_2=C_diet_2,C_diet_3=C_diet_3)
    } else {
      if(length(which(!levels(testData$DietRIX) %in% testData$DietRIX)) > 0){
        newLev = levels(testData$DietRIX)[-which(!levels(testData$DietRIX) %in% testData$DietRIX)]
        testData$DietRIX = factor(testData$DietRIX, levels = newLev)
      }
      
      if(length(which(!levels(testData$seq.Gene) %in% testData$seq.Gene)) > 0){
        newLev = levels(testData$seq.Gene)[-which(!levels(testData$seq.Gene) %in% testData$seq.Gene)]
        testData$seq.Gene = factor(testData$seq.Gene, levels = newLev)
        testData$pup_gene = factor(paste(testData$Pup.ID, testData$seq.Gene, sep="_"))
      }

      encoded <- unique(getEncoding(testData, terms = unique(c(terms,"Pup.ID","seq.Gene","pup_gene"))))

      reg[[paste0("rix_",i)]] = list()
      reg[[paste0("rix_",i)]]$stan = stan_ase(df=testData, encoded=encoded, 
                                         mu_g0=0.5, t20=0.001, 
                                         nchains=2, iter=niter, 
                                         C_diet_4=C_diet_4, C_diet_2=C_diet_2,C_diet_3=C_diet_3,
                                         terms = terms,
                                         stanMod=stanMod)
      summ = data.frame(summary(reg[[paste0("rix_",i)]]$stan)$summary)
      if(length(grep("eta|ind", rownames(summ))) > 0) summ = summ[-grep("eta|ind", rownames(summ)),]
      summ$param = rownames(summ)
      summ$param[grep("mu_g[[]", rownames(summ))] = encoded$Level[which(encoded$Variable == "SEQ.GENE")]
      summ$param[grep("mu_p[[]", rownames(summ))] = encoded$Level[which(encoded$Variable == "PUP.ID")]
      
      reg[[paste0("rix_",i)]]$sum = summ
      
      freq_test = testData %>% group_by(seq.Gene) %>%
        summarize(y = sum(CC_1), N = sum(sum))
      freq_res = apply(freq_test, 1, function(x) 
          binom.test(as.numeric(x["y"]), as.numeric(x["N"]), p = 0.5))
      names(freq_res) = freq_test$seq.Gene
      reg[[paste0("rix_",i)]]$freq = freq_res
    }
  }
  
  return(reg) 
}


############# encoding stuff #################
getEncoding <- function(df, terms){
  encoded <- data.frame()
  #if(length(grep("gene", terms)) > 0) terms = terms[-grep("gene", terms)]
  for(i in 1:length(terms)){
    
    var <- toupper(paste(terms[i]))
    ind <- match(var, toupper(colnames(df)))
    if(!is.na(ind) > 0){
      useLev = unique(df[,ind])
      if(!is.null(levels(df[,ind]))) useLev = levels(df[,ind])[which(levels(df[,ind]) %in% useLev)]
      
      #vec <- factor(t(df[,ind]), levels=useLev)
      len <- length(unlist(useLev))       #length(levels(as.factor(vec)))
      tempdf <- data.frame(Level = as.character(paste(unlist(useLev))), 
                           Index = 1:len, 
                           Variable = rep(var,len))
      encoded <- rbind(encoded,tempdf)
    }
  }
  return(encoded)
}

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


array.summarize <- function(array, labels=NULL){
  if(class(array) == "list" & class(array[[1]]) %in% c("list", "mcarray")){   
    nme <- names(array)
    listLength <- length(array)
  } else {
    listed <- F
    nme <- NULL
    listLength <- 1
  }
  sum_list <- list()
  for(i in 1:listLength){
    comb_int = means = NULL
    len3 <- length(dim(as.array(array[[i]])))
    nchains <- dim(as.array(array[[i]]))[len3]
    if(len3 == 3){
      est_comb <- sapply(1:dim(array[[i]])[1], function(x){
        sapply(1:nchains, function(y){
          tmp = as.vector(unlist(array[[i]][x,,y]))
          tmp[!is.na(tmp)]}, simplify=F)
      }, simplify=F)
      sep_int <- lapply(est_comb, function(x){
        data.frame(lapply(x, function(y) HPDinterval(as.mcmc(y))))  })
    } else if (len3 == 4){
      est_comb <- sapply(1:dim(array[[i]])[1], function(x){
        sapply(1:dim(array[[i]])[2], function(y){
          sapply(1:nchains, function(z){
            tmp = as.vector(unlist(array[[i]][x,y,,z]))
            tmp[!is.na(tmp)]}, simplify=F)
        }, simplify=F) }, simplify=F)
      sep_int <- lapply(est_comb, function(x){
        do.call("rbind", lapply(x, function(y) 
          data.frame(lapply(y, function(z) HPDinterval(as.mcmc(z)))) ))  })
      comb_int <- lapply(est_comb, function(x) 
        do.call("rbind", lapply(x, function(y) HPDinterval(as.mcmc(unlist(y)))) ) )
      comb_int <- do.call("rbind", comb_int)
      means <- lapply(est_comb, function(x) 
        do.call("rbind", lapply(x, function(y) c(Mode(unlist(y)), mean(unlist(y)), median(unlist(y))))) )
      means <- do.call("rbind", means)
      colnames(means) <- c("Mode","Mean","Median")  
      
    } else {
      est_comb <- lapply(array, unlist)
      sep_int <- lapply(est_comb, function(x){
        data.frame(HPDinterval(as.mcmc(x)))  })
    }
    
    sep_int <- do.call("rbind", sep_int)
    colnames(sep_int)[1:2] <- c(paste0("lower.", nchains), paste0("upper.", nchains))
    
    #sep_ints <- sapply(1:dim(array)[1], function(x) HPDinterval(as.mcmc.list(array)), simplify=F)
    #names(sep_ints) <- paste0("ints",seq(1:nchains))
    #sep_int <- do.call("cbind", sep_ints)
    if(is.null(comb_int)){
      comb_int <- lapply(est_comb, function(x) HPDinterval(as.mcmc(unlist(x))))
      comb_int <- do.call("rbind", comb_int)
    }
    if(is.null(means)){
      means <- lapply(est_comb, function(x) c(Mode(unlist(x)), mean(unlist(x)), median(unlist(x))))
      means <- do.call("rbind", means)
      colnames(means) <- c("Mode","Mean","Median")
    }
    
    summary <- cbind(means, comb_int, sep_int)
    if(is.null(labels) | is.null(nme)){
      if(is.null(rownames(summary))) rownames(summary) <- seq(1:dim(array[[i]])[1])
    } else {
      v <- gsub("B_|PO", "", toupper(nme[i]))
      if(toupper(nme[i]) == "MU_P") v <- "PUP.ID"
      if(toupper(nme[i]) == "MU_G") v <- "PUP_GENE"
      if(toupper(nme[i]) == "MU_R") v <- "RRIX"
      if(toupper(nme[i]) == "B_DIETPORIX") v <- c("DIET","RRIX")
      v = ifelse(v == "RIX" & length(v) == 1, "RRIX", v)
      ind <- which(toupper(labels$Variable) %in% v)
      if(any(ind) == 0){
        if(is.null(rownames(summary))){
          rownames(summary) <- seq(1:dim(array[[i]])[1])
        } else {
          rownames(summary) <- gsub("var","", rownames(summary))
        }
      } else if(length(v) > 1){
        names = sapply(1:length(v), function(i)
          labels$Level[which(toupper(labels$Variable) %in% v[i])], simplify=F)
        names_ex = expand.grid(names[[2]], names[[1]])
        rownames(summary) = paste0(names_ex$Var1, "_", names_ex$Var2)
      } else {
        rownames(summary) <- unique(gsub(" ", "_", labels$Level[ind]))
      }
    }
    sum_list[[i]] <- summary
  }
  if(listLength == 1){
    sum_list <- sum_list[[1]]
  } else {
    names(sum_list) <- nme
  }
  
  return(summary = sum_list)
}

