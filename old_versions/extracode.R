####from fitting rna stuff

# fit models
fitLMER <- BC.model(y.mat=datause[,-noPhen], data=datause,
                    indvariable=indvariable, 
                    transformParams=getMatnutTransformParams(tryLam = tryLam, normd = T))

lmerSum <- lmer.getDecodedSummary(lmerOb=fitLMER, formula=indvariable, phenotypes=fitLMER$phenotypes)
fitLMER$JAGSformula <- lmer.to.jags.form(indvariable)
jagsFit <- runJagsModels_cov(datalist=datause, 
                             testLMER = fitLMER, encoded=encoded, phenotype=fitLMER$phenotypes, n.iter=100000)

jagsSum <- list()
all.y <- matrix(NA, nrow=nrow(datause), ncol=length(fitLMER$y.transform))
for(i in 1:length(fitLMER$y.transform)){
  all.y[,i] <- fitLMER$y.transform[[i]]
  all.mcmc <- mcmc.stack(jagsFit[[i]]$fit)
  jagsSum[[i]] <- jags.getDecodedSummary(all.mcmc, encoded)
}

allSummary <- list(lmer=lmerSum, jags=jagsSum)

allSummary$compare <- sapply(1:length(lmerSum), function(x) data.frame(compareSummary(allSummary$lmer[[x]],allSummary$jags[[x]])), 
                             simplify = F)


### from matching in summary

if(addS == "fit"){
  keepMCMC <- all.mcmc      
  
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
  
  temp1 <- as.mcmc(cbind(all.mcmc, PO_net))
  
  all.mcmc <- temp1
  mult <- all.mcmc[,-c(grep("^d", colnames(all.mcmc)), grep("drop", colnames(all.mcmc)))]
  nomult <- all.mcmc[,-c(grep("^d", colnames(all.mcmc)), 
                         setdiff(grep("^PO", colnames(all.mcmc)), grep("_", colnames(all.mcmc)))) ]
  colnames(nomult) <- gsub("_drop", "", colnames(nomult))
  
  all.mcmc <- nomult
} else if(addS == "force"){
  rixes <- unique(dataf$RIX)
  for(i in 1:length(rixes)){  
    sumPO <- sum(colMeans(all.mcmc)[intersect(grep("^PO", colnames(all.mcmc)), grep(paste0(rixes[i], "$"), colnames(all.mcmc)))])
    S[i] <- ifelse(sumPO < 0, -1, 1)
    dataf[which(dataf$RIX == rixes[i]),"PO"] <- S[i]*dataf[which(dataf$RIX == rixes[i]),"PO"]
  }
  
  jagsFit <- runJagsModels_cov(datalist=dataf, chains=chains, testLMER = lmerFit, 
                               n.adapt=n.adapt, n.iter=n.iter, thin=thin, 
                               encoded=encoded, phenotype=phenotype, sq=sq)
  all.reg <- jagsFit[[phenotype]]$fit
  all.mcmc <- mcmc.stack(all.reg)
  madeMod[[2]] <- jagsFit[[phenotype]]$madeModel
}

if(is.null(lmerFit)){
  allSummary <- NULL
  null.test=NULL
  diet.test=NULL 
  rix.test=NULL
  S=NULL
  transf=NULL
  JL_compare=NULL
  lmer_obj = NULL
  modelDef= NULL               #madeMod
  summary_table= NULL          #summary_table
  allSummary= NULL             #allSummary
  plot= NULL                   #plotCat
} else {
  
  #if (MatchS){
  #  lmerSum <- lmerSum[-which(gsub("\\s","",lmerSum$Variable) == "s"), ]
  #}
  
  #if (!Match){
  #  diet.test <- contrast.test(data=as.data.frame(all.mcmc), encoded, variable="Diet", byVar="RIX", contrast=test_cm$Diet.contr)
  #  rix.test <- contrast.test(data=as.data.frame(all.mcmc), encoded, variable="RIX", byVar="Diet", contrast=test_cm$RIX.contr)
  #} else {