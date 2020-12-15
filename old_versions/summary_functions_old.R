source("./plot.hpd.R")
source(file.path(".", "matnut", "jags_functions.R"))
source(file.path(".", "matnut", "lmer_functions.R"))
source(file.path(".", "matnut", "contrast_functions.R"))

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


makeSummary <- function(datalist, phenotype, tryLam=1, Match=F, sq=F, normd=T, addS=c("off", "force", "fit"), 
                        chains=2, n.adapt=2500, n.iter=10000, thin=10,  
                        encoded=NA){
  if(class(datalist) == "list"){
    dataf <- datalist$df
    encoded <- datalist$encoded
    if(is.na(phenotype)){
      phenotype <- datalist$ptypes
    } 
  } else {
    dataf <- datalist
  }
  addS <- addS[1]
  S <- c()
  
  lmerFit <- runLMERModels_cov(phenotype=phenotype, dataf=dataf, tryLam=tryLam, normd=normd, Match=Match)
  jagsFit <- runJagsModels_cov(datalist=dataf, chains=chains, testLMER = lmerFit, 
                               n.adapt=n.adapt, n.iter=n.iter, thin=thin, 
                               encoded=encoded, phenotype=phenotype, sq=sq, addS=addS)
browser()
  all.reg <- jagsFit[[phenotype]]$fit
  all.mcmc <- mcmc.stack(all.reg)
  
  rixes <- unique(dataf$RIX)
  for(i in 1:length(rixes)){  
    sumPO <- sum(colMeans(all.mcmc)[intersect(grep("^PO", colnames(all.mcmc)), grep(paste0(rixes[i], "$"), colnames(all.mcmc)))])
    S[i] <- ifelse(sumPO < 0, -1, 1)
  }
  
  if(addS == "estimate"){
    keepMCMC <- all.mcmc        #[,-grep("^d", colnames(all.mcmc))]
    
    PO_net <- data.frame(matrix(nrow=length(colnames(keepMCMC)[grep("^PO", colnames(keepMCMC))]), ncol=nrow(keepMCMC)), 
                         row.names=colnames(keepMCMC)[grep("^PO", colnames(keepMCMC))])
    
    rixes <- encoded$Level[which(encoded$Variable == "RIX")]
    diets <- encoded$Level[which(encoded$Variable == "Diet")]
    
    for(i in 1:length(rixes)){
      PO_net[paste0("PO",rixes[i]),] <- keepMCMC[,paste0("PO",rixes[i])] * keepMCMC[,paste0("betaS[",i,"]")]
      for(j in 1:length(diets)){
        if(length(grep(paste0("PO",diets[j], rixes[i]), colnames(keepMCMC))) > 0){
          temp <- keepMCMC[,paste0("PO", diets[j], rixes[i])] * keepMCMC[,paste0("betaS[",i,"]")]
          PO_net[paste0("PO", diets[j], rixes[i]),] <- temp
        }
      }
    }   
    
    PO_net <- as.mcmc(t(PO_net))
    colnames(all.mcmc)[grep("^PO",colnames(all.mcmc))] <- paste(colnames(all.mcmc)[grep("^PO",colnames(all.mcmc))],"drop", sep="_")
    
    temp1 <- cbind(all.mcmc, PO_net)
    temp2 <- as.mcmc(temp1[,-grep("drop", colnames(temp1))])
    
    all.mcmc <- temp2
    
  } else if(addS == "force"){
    rixes <- unique(dataf$RIX)
    for(i in 1:length(rixes)){  
      #sumPO <- sum(colMeans(all.mcmc)[intersect(grep("^PO", colnames(all.mcmc)), grep(paste0(rixes[i], "$"), colnames(all.mcmc)))])
      #S[i] <- ifelse(sumPO < 0, -1, 1)
      dataf[which(dataf$RIX == rixes[i]),"PO"] <- S[i]*dataf[which(dataf$RIX == rixes[i]),"PO"]
    }
    jagsFit <- runJagsModels_cov(datalist=dataf, chains=chains, testLMER = lmerFit, 
                                 n.adapt=n.adapt, n.iter=n.iter, thin=thin, 
                                 encoded=encoded, phenotype=phenotype, sq=sq)
    all.reg <- jagsFit[[phenotype]]$fit
    all.mcmc <- mcmc.stack(all.reg)
  }
  #browser()
  jagsSum <- jags.getDecodedSummary(all.mcmc, encoded)
   
  fixed = "Diet"    #ifelse(Match, "",
  lmerSum <- lmer.getDecodedSummary(lmerOb=lmerFit, phenotype)
  #if (MatchS){
  #  lmerSum <- lmerSum[-which(gsub("\\s","",lmerSum$Variable) == "s"), ]
  #}
  transf <- c("lambda" = lmerFit[[phenotype]]$phen_1$lambda,  #testLMER[[phenotype]]$lambda
              "pval" = lmerFit[[phenotype]]$phen_1$best.pval)

  allSummary <- compareSummary(lmerSum, jagsSum)
  
  keep <- c("Diet", "DietRIX","PODietRIX","PORIX", "RIX")
  
  getSummary <- allSummary[which(allSummary$Variable %in% keep),]
  
  test_cm <- contrast_matrix(encoded)
  null.test <- contrast.null(dfObject=as.data.frame(all.mcmc), levels=unique(allSummary$Level))
  #if (!Match){
  #  diet.test <- contrast.test(data=as.data.frame(all.mcmc), encoded, variable="Diet", byVar="RIX", contrast=test_cm$Diet.contr)
  #  rix.test <- contrast.test(data=as.data.frame(all.mcmc), encoded, variable="RIX", byVar="Diet", contrast=test_cm$RIX.contr)
  #} else {
    diet.test <- NULL
    rix.test <- NULL
  #}
  
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
  jagsLmer_compare <- cbind(phenotype=rep(phenotype,nrow(getSummary)), getSummary) 
  
  ####
  plotCat <- plotSummary(lmerSum, jagsSum, getSummary, phenotype)
  return(list(transformation=transf, mcmcObject=all.mcmc, plot=plotCat, allSummary=allSummary,
              null.test=null.test, diet.test=diet.test, rix.test=rix.test, S=S, 
              summary_table=summary_table, JL_compare=jagsLmer_compare, lmer_obj = lmerFit[[phenotype]]))
}

#colMeans(fitJags[[phenotype[i]]][[1]])[grep("betaRIX", colnames(fitJags[[phenotype[i]]][[1]]))]

compareSummary <- function(lmerSum, jagsSum){
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
  return(diffSummary)
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




