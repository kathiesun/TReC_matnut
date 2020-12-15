library(rjags)
library(coda)
library(ggplot2)
library(data.table)


source("./matnut/jags_model.R")
source("./matnut/encode.R")
source("./lm/fitBoxCoxModels.R")
source("./mnp/simdata.R")

# setwd("C:/Users/kathie/Dropbox/source/nutridiallel2014/src")

testdata <- getTestCase()
################
# Prepare data #
################

matnut <- read.csv("../data/AllPhenotypes_Matnut5.csv")
#matnut <- testdata$cov.data
matnut[grep("a", matnut$Reciprocal),"PO"] <- 0.5
matnut[grep("b", matnut$Reciprocal),"PO"] <- -0.5
matnut[which(matnut$PO == 0.5), "pof"] <- "+"
matnut[which(matnut$PO == -0.5), "pof"] <- "-"
dietlabs <- c("STD", "ME", "PD", "VDD")
matnut$RIX <- factor(matnut$RIX, levels=c(1:4,6:10))
matnut$Diet <- factor(matnut$Diet, levels=dietlabs)
matnut$DamID <- factor(matnut$DamID)
#matnut$Diet <- factor(paste0("Diet",matnut$Diet))
#matnut$RIX <- factor(paste0("RIX", matnut$RIX))
matnut$BBoriginal <- matnut$BreedingBatch
matnut$BreedingBatch <- factor(paste0("bb",matnut$BBoriginal), levels=unique(paste0("bb",matnut$BBoriginal)))
matnut$Cage <- factor(matnut$Cage)

orderdat <- matnut[order(matnut$RIX, matnut$Diet, matnut$PO),]
orderdat$PORIX <- paste0("PO",orderdat$RIX)
orderdat$PORIX <- factor(orderdat$PORIX, levels=unique(orderdat$PORIX))
orderdat$DietRIX <- paste0(orderdat$Diet, orderdat$RIX)
orderdat$DietRIX <- factor(orderdat$DietRIX, levels=unique(orderdat$DietRIX))
orderdat$PODietRIX <- paste0("PO",orderdat$DietRIX)
orderdat$PODietRIX <- factor(orderdat$PODietRIX, levels=unique(orderdat$PODietRIX))

variables <- c("Diet", "RIX", "PORIX", "DietRIX", "PODietRIX","BreedingBatch", "DamID","Cage")
dispersions <- c("sigInv2", "tauInv2","tauDRInv2", "tauPEInv2", 
                 "tauPDRInv2","tauBBInv2", "tauDamInv2","tauCageInv2")

encoded <- getEncoding(orderdat, variables)

ptypes <- colnames(matnut)[44:63]
#ptypes <- c("OFTotalDistance","BasalCORT","StressCORT")

allparam <- c(variables, dispersions)
allparam <- allparam[order(allparam)]

#############
# Contrasts #
#############
ndiet=length(dietlabs)
nrix=length(encoded$Level[which(encoded$Variable == "RIX")])
diet.contr <- matrix(NA, nrow=choose(ndiet,2), ncol=2)
rix.contr <- matrix(NA, nrow=choose(nrix, 2), ncol=2)
count=1
ind=1
for(i in 1:ndiet){
  ind = i+1
  while(ind <= ndiet && count <= nrow(diet.contr)){
    diet.contr[count,] <- c(i, ind)   
    count=count+1
    ind=ind+1
  }
}

count=1
ind=1
for(i in 1:nrix){
  ind = i+1
  while(ind <= nrix && count <= nrow(rix.contr)){
    rix.contr[count,] <- c(i, ind)   
    count=count+1
    ind=ind+1
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
  
source("./matnut/encode.R")
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
  testPhen[[ptypes[i]]] <- allSummaries_rand(phenotype = ptypes[i], orderdat, allparam, encoded, Batch=F, Dam=F, Cage=F, 
                                             tryLam=c(-1, 0, .25, .33, .5, 1, 2, 3))
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

#############
# Contrasts #
#############

### Make summary tables ###
dietDiet_table <- list()
rixRix_table <- list()


for (i in 1:length(ptypes)){
  ### contrasts ###
  dietTab <- contrasts.getDecodedSummary(testPhen[[ptypes[i]]], "Diet")
  dietDiet_table[[ptypes[i]]] <- cbind(phenotype=rep(ptypes[i],nrow(dietTab)), dietTab)
  rixTab <- contrasts.getDecodedSummary(testPhen[[ptypes[i]]], "RIX")
  rixRix_table[[ptypes[i]]] <- cbind(phenotype=rep(ptypes[i],nrow(rixTab)), rixTab)
}

write.csv(do.call(rbind, dietDiet_table), file.path("./","matnut_outputs/",'allPheno_dietDiet_table.csv'), row.names = F)
write.csv(do.call(rbind, rixRix_table), file.path("./","matnut_outputs/",'allPheno_rixRix_table.csv'), row.names = F)

### Plots ###

byVar = "Diet"
contrastPlots <- list()
returnPval <- list()
for(i in 1:length(ptypes)){  
  useTest <- ifelse(byVar=="Diet","diet.test","rix.test")
  bindTemp <- testPhen[[ptypes[i]]][[useTest]]$dir
  rownames(bindTemp) <- rownames(testPhen[[ptypes[i]]][[useTest]]$pval)
  colnames(bindTemp) <- paste0("a",colnames(testPhen[[ptypes[i]]][[useTest]]$pval))
  meltTemp <- melt(bindTemp)
  colnames(meltTemp)[3] <- "Dir"
  tempPval <- testPhen[[ptypes[i]]][[useTest]]$pval
  colnames(tempPval) <- paste0("a",colnames(testPhen[[ptypes[i]]][[useTest]]$pval))
  
  meltPval <- merge(melt(tempPval), meltTemp)
  meltPval$Var2 <- factor(gsub("a","",meltPval$Var2), levels=unique(gsub("a","",meltPval$Var2)))
  
  
  meltPval$effectName1 <- unlist(strsplit(as.character(paste(meltPval$Var2)),"[.]"))[c(T, F)]
  meltPval$effectName2 <- unlist(strsplit(as.character(paste(meltPval$Var2)),"[.]"))[c(F, T)]
  meltPval$effectName1 <- factor(meltPval$effectName1, 
                                 levels=encoded$Level[which(encoded$Variable==byVar)])
  meltPval$effectName2 <- factor(meltPval$effectName2, 
                                 levels=encoded$Level[which(encoded$Variable==byVar)])
  meltPval$Level <- meltPval$Var1
  if(byVar=="Diet"){
    meltPval$Level <- gsub('\\D+','', meltPval$Var1) 
    meltPval$Level <- factor(meltPval$Level, levels=unique(meltPval$Level))
    meltPval$LevelPO <- gsub("[0-9ME]",'', meltPval$Var1)
  } else {
    meltPval$Level <- gsub('PO','', meltPval$Var1) 
    meltPval$Level <- factor(meltPval$Level, levels=unique(meltPval$Level))
    meltPval$LevelPO <- rep("", length(meltPval$Level))
    meltPval$LevelPO[grep("PO",meltPval$Var1)] <- "PO"
  }
  
  meltPval$dirNum <- ifelse(meltPval$Dir == "+", 1, -1)
  returnPval[[i]] <- meltPval
  
  dfplot <- data.frame(pval = returnPval[[i]]$value, effectName1 = returnPval[[i]]$effectName1, 
                       effectName2 = returnPval[[i]]$effectName2, Level=returnPval[[i]]$Level, 
                       LevelPO=returnPval[[i]]$LevelPO, direction=returnPval[[i]]$dirNum)
  contrastPlots[[i]] <- getPlot.compare.helper(dfplot,ptypes=ptypes,byVar=byVar)
  
}

##############
# Prediction #
##############

dietOnly = FALSE
ribbon9_Dietonly <- list()
ribbon9 <- list()
ribbon4 <- list()

for (j in 1:length(ptypes)){
  all.mcmc <- testPhen[[ptypes[j]]]$mcmcObject
  #all.reg <- list.reg[[j]]
  #all.mcmc <- mcmc.stack(all.reg)
  all.dist <- cbind(all.mcmc, rep(0,nrow(all.mcmc)), rep(0,nrow(all.mcmc))) #, diffdiet_post, diffRIX_post)
  colnames(all.dist)[c(length(colnames(all.dist))-1, length(colnames(all.dist)))] <- c("emp", "POemp")
  
  dietrix <- matrix(c(as.character(encoded$Level[which(encoded$Variable == "DietRIX")]),"emp", "emp"), 
                    ncol=9, nrow=4)
  dietrix[,9] <- c("emp","ME10", "PD9", "emp")
  
  predPOneg <- list()
  predPOpos <- list()
  predDiet <- list()
  
  for(i in 1:nrow(all.dist)){
    for (rixn in 1:nrix){
      for (dietn in 1:length(dietlabs)){
        
        rixlabs <- as.character(encoded$Level[which(encoded$Variable == "RIX")])
        dietlabs <- as.character(encoded$Level[which(encoded$Variable == "Diet")])
        
        temprix <- rixlabs[rixn]
        tempdiet <- dietlabs[dietn]
        tempinter <- dietrix[dietn,rixn]
        prLab <- paste0(tempdiet,temprix)
        
        if (!dietOnly){
          predPOneg[[paste0(prLab, ".neg")]][i] <- all.dist[i,temprix] + all.dist[i,tempdiet] + 
            all.dist[i,dietrix[dietn,rixn]] + all.dist[i,paste0("PO", temprix)]*(-0.5) + 
            all.dist[i,paste0("PO", tempinter)]*(-0.5)
          
          predPOpos[[paste0(prLab, ".pos")]][i] <- all.dist[i,temprix] + all.dist[i,tempdiet] + 
            all.dist[i,dietrix[dietn,rixn]] + all.dist[i,paste0("PO", temprix)]*(0.5) + 
            all.dist[i,paste0("PO", tempinter)]*(0.5)
        } else {
          predDiet[[prLab]][i] <- all.dist[i,tempdiet] + all.dist[i,dietrix[dietn,rixn]]
        }
      }
    }
  }
  if (!dietOnly){
    allpredPOneg <- do.call(cbind, predPOneg)
    allpredPOpos <- do.call(cbind, predPOpos)
    allpredPO <- list(allpredPOneg, allpredPOpos) #, all.dist[,rixlabs]
    flydata <- as.data.frame(c(allpredPO))
  } else {
    allpredDiet <- do.call(cbind,predDiet)
    flydata <- as.data.frame(allpredDiet)
  }
  
  flydata <- as.mcmc(flydata)
  mu    <- c(colMeans(flydata))
  med   <- apply(coda::HPDinterval(flydata, prob=0.01), 1, median)
  hpd.wide    <- coda::HPDinterval(flydata, prob=0.95)
  hpd.narrow  <- coda::HPDinterval(flydata, prob=0.5)
  printdata <- data.frame(cbind(mu,med,hpd.wide,hpd.narrow))
  printdata$dietrix <- do.call(rbind,strsplit(rownames(printdata),"[.]"))[,1] 
  printdata$diet <- gsub('[0-9]+', '', printdata$dietrix)
  printdata$rix <- gsub('[A-Z]+', '', printdata$dietrix)
 
  if (dietOnly){
    printdata <- printdata[-which(rownames(printdata) %in% c("STD10","VDD10","STD10","VDD10")),]
  } else {
    printdata$poe <- printdata$dietrix <- do.call(rbind,strsplit(rownames(printdata),"[.]"))[,2]
    printdata$poen <- ifelse(printdata$poe == "pos", 0.5, -0.5)
    printdata$rixpo <- paste0(printdata$rix, printdata$poe)
    printdata <- printdata[-which(rownames(printdata) %in% c("STD10.neg","VDD10.neg","STD10.pos","VDD10.pos")),]
  }
  
  dietmeans <- aggregate(printdata[,1:2], list(printdata$diet), mean)
  dietmeans <- dietmeans[order(dietmeans$mu),]
  
  printdata$diet <- factor(printdata$diet)    #, levels = dietmeans$Group.1)
  printdata$rix <- factor(printdata$rix, levels = as.character(encoded$Level[which(encoded$Variable == "RIX")]))
  
  if (dietOnly){
    ribbon9_Dietonly[[ptypes[j]]] <- ggplot(printdata, aes(color=rix)) + #geom_point(aes(y=mu, x=diet, alpha=poe)) + 
      theme_minimal() + 
      geom_line(aes(y=mu, x=diet, group=rix)) +  
      geom_ribbon(color=NA, alpha=0.7, aes(x=diet, ymin=lower.1, ymax=upper.1,group=rix,
                                           fill=rix, color=rix)) + 
      geom_hline(yintercept = 0, colour="lightsteelblue", size=0.5, linetype="longdash") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      scale_fill_discrete(guide=FALSE) + scale_color_discrete(guide=FALSE) + 
      facet_wrap(~rix) + 
      ylab("Predicted effect") + xlab("") + ggtitle(ptypes[j])
    
  } else {
    ribbon9[[ptypes[j]]] <- ggplot(printdata, aes(color=rix)) + #geom_point(aes(y=mu, x=diet, alpha=poe)) + 
      theme_minimal() + 
      geom_line(aes(y=mu, x=diet, group=poe, alpha=poe)) +  
      geom_ribbon(color=NA, aes(x=diet, ymin=lower.1, ymax=upper.1, group=poe, alpha=poe, 
                                fill=rix, color=rix)) + 
      geom_hline(yintercept = 0, colour="lightsteelblue", size=0.5, linetype="longdash") +
      scale_alpha_discrete(range = rev(c(0.45, 0.75)), name="POE") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      scale_fill_discrete(guide=FALSE) + scale_color_discrete(guide=FALSE) + 
      facet_wrap(~rix) + 
      ylab("Predicted effect") + xlab("") + ggtitle(ptypes[j])
    
    ribbon4[[ptypes[j]]] <- ggplot(printdata, aes(color=diet)) + 
      geom_line(aes(y=mu, x=rix, group=poe, alpha=poe)) + 
      geom_ribbon(color=NA, aes(x=rix, ymin=lower.1, ymax=upper.1, group=poe, alpha=poe, 
                                fill=diet, color=diet)) + 
      geom_hline(yintercept = 0, colour="lightsteelblue", size=0.5, linetype="longdash") +
      scale_alpha_discrete(range = rev(c(0.45, 0.75)), name="POE") +
      scale_fill_discrete(guide=FALSE) + scale_color_discrete(guide=FALSE) + 
      facet_wrap(~diet) + theme_minimal() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      ylab("Predicted effect") + xlab("") + ggtitle(ptypes[j])
  }
}

#########
# Plots #
#########

graphics.off()
pdf(file.path("./","matnut_outputs/","allPheno_DietContrastPlots.pdf"), 
    onefile=TRUE, width = 8.5, height = 11)
for(i in 1:length(ptypes)){  
  print(contrastPlots[[i]])
}
graphics.off()

graphics.off()
pdf(file.path("./","matnut_outputs/","allPheno_JLComparison.pdf"), 
    onefile=TRUE, width = 8.5, height = 11)
for(i in 1:length(ptypes)){  
  print(testPhen[[ptypes[i]]]$plot)
}
graphics.off() 

graphics.off()
pdf(file.path("./","matnut_outputs/","allPheno_ribbon.pdf"), onefile=TRUE, width = 14, height = 11)
for(i in 1:length(ptypes)){  
  print(ribbon9[[ptypes[i]]])
  print(ribbon4[[ptypes[i]]])
}
graphics.off()

graphics.off()
pdf(file.path("./","matnut_outputs/","allPheno_DietOnlyribbon.pdf"), onefile=TRUE, width = 14, height = 11)
for(i in 1:length(ptypes)){  
  print(ribbon9_Dietonly[[ptypes[i]]])
}
graphics.off()

###############
# Make tables #
###############

dietDiet_table <- list()
jagsLmer_compare <- list()
rixRix_table <- list()
summary_table <- list()

psychPhen <- testPhen

for(i in 1:length(ptypes)){
  pheno <- ptypes[i]
  
  ####
  tempTable <- jags.getDecodedSummary(psychPhen[[pheno]]$mcmcObject, encoded, narrow=0.5, wide=0.95)
  tempTable$Level <- factor(tempTable$Level, levels=tempTable$Level)
  null_p <- data.frame(Level = names(psychPhen[[pheno]]$null.test$pval), 
                       pval = psychPhen[[pheno]]$null.test$pval, 
                       signif = psychPhen[[pheno]]$null.test$signif)
  null_p$Level <- factor(null_p$Level, levels=tempTable$Level)
  null_p <- null_p[order(null_p$Level),]
  mergeTable <- merge(tempTable, null_p, by="Level")
  mergeTable$Level <- factor(mergeTable$Level, levels=encoded$Level)
  mergeTable$Variable <- factor(mergeTable$Variable, levels=unique(encoded$Variable))
  
  summary_table[[pheno]] <- cbind(phenotype=rep(pheno,nrow(mergeTable)), mergeTable[order(mergeTable$Variable, mergeTable$Level),])
  
  
  ####
  #dietTab <- contrasts.getDecodedSummary(psychPhen$OFTotalDistance, "Diet")
  #dietDiet_table[[pheno]] <- cbind(phenotype=rep(pheno,nrow(dietTab)), dietTab)
  #rixTab <- contrasts.getDecodedSummary(psychPhen$OFTotalDistance, "RIX")
  #rixRix_table[[pheno]] <- cbind(phenotype=rep(pheno,nrow(rixTab)), rixTab)
  
  
  
}











