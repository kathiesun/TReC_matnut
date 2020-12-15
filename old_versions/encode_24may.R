source("./plot.hpd.R")
source("./jags_functions.R")
source("./lmer_functions.R")
source("./contrast_functions.R")

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

#####################
# Summaries of fits #
#####################


makeSummary <- function(data, phenotype, tryLam=1, Match=F, checkAnova=F, 
                             chains=2, n.adapt=2500, n.iter=10000, thin=10)
{
  dataf <- data$df
  encoded <- data$encoded

  lmerFit <- runLMERModels_cov(phenotype, dataf, tryLam, normd, Match, checkAnova)
  jagsFit <- runJagsModels_cov(phenotype, dataf, encoded, chains=chains, testLMER = lmerFit, 
                               n.adapt=n.adapt, n.iter=n.iter, thin=thin)
  fixed = ifelse(Match, "","Diet") 
  lmerSum <- lmer.getDecodedSummary(lmerObject=lmerFit, phenotype)
  transf <- c("lambda" = lmerFit[[phenotype]]$phen_1$lambda,  #testLMER[[phenotype]]$lambda
              "pval" = lmerFit[[phenotype]]$phen_1$best.pval)

  all.reg <- jagsFit[[phenotype]]
  all.mcmc <- mcmc.stack(all.reg)
  jagsSum <- jags.getDecodedSummary(all.mcmc, encoded)
  allSummary <- compareSummary(lmerSum, jagsSum)
  
  keep <- c("Diet", "DietRIX","PODietRIX","PORIX", "RIX")
  
  getSummary <- allSummary[which(allSummary$Variable %in% keep),]
  
  test_cm <- contrast_matrix(encoded)
  null.test <- contrast.null(dfObject=as.data.frame(all.mcmc), levels=unique(allSummary$Level))
  if (!Match){
    diet.test <- contrast.test(dfObject=as.data.frame(all.mcmc), encoded, variable="Diet", byVar="RIX", contrast=test_cm$Diet.contr)
    rix.test <- contrast.test(dfObject=as.data.frame(all.mcmc), encoded, variable="RIX", byVar="Diet", contrast=test_cm$RIX.contr)
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
  jagsLmer_compare <- cbind(phenotype=rep(phenotype,nrow(getSummary)), getSummary) 
  
  ####
  
  plotCat <- combinedCat(lmerSum, jagsSum, getSummary, phenotype)
  return(list(transformation=transf, mcmcObject=all.mcmc, plot=plotCat, allSummary=allSummary,
              null.test=null.test, diet.test=diet.test, rix.test=rix.test, 
              summary_table=summary_table, JL_compare=jagsLmer_compare, lmer_obj = lmerFit[[phenotype]]))
}


compareSummary <- function(lmerSum, jagsSum){
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




