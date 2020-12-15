library(rjags)
library(coda)
library(ggplot2)
library(data.table)
library(lmerTest)

source(file.path(".", "matnut", "summary_functions.R"))
source("./matnut/prediction_functions.R")
source("./matnut/matching_functions3.R")
source("./matnut/matching_functions.R")
source("./mnp/simdata.R")

# setwd("C:/Users/kathie/Dropbox/source/nutridiallel2014/src")

###### USE THIS!! ######
matnut <- readRDS(file.path(dataSource,'matnut_data.rds'))
matnut <- readRDS(file.path("./","matnut_outputs/",'matnut_data.rds'))

####################
#   Run function   #
####################

myPhen <- list()
for(i in 1:length(matnut$ptypes)){
  myPhen[[matnut$ptypes[i]]] <- makeSummary(datalist=matnut, phenotype=matnut$ptypes[i],
                                            tryLam=c(-1, 0, .25, .33, .5, 1, 2, 3), normd=T, 
                                            chains=2, n.adapt=20000, n.iter=50000, thin=10, sq=F, addS=)
}

#saveRDS(myPhen, file.path("./","matnut_outputs/",'allPheno_modelfits.rds'))
testPhen <- readRDS("./matnut_outputs/allPheno_modelfits.rds")
matnut <- readRDS(file.path(dataSource,'allPheno_modelfits.rds'))


multimp2 <- list()
### matching ###
for(i in (length(matnut$ptypes)-2):length(matnut$ptypes)){    #length(matnut$ptypes)){
  multimp2[[matnut$ptypes[i]]] <- match.multimp(data=matnut, ptypes=matnut$ptypes[i], N=2, addS="fit",
                                                matchon="Cage", matchoff="PO", idcol="ID")
}


#saveRDS(multimp, file.path("./","matnut_outputs/","mult_imp_match.rds"))


############################
# Generate contrast tables #
############################

dietDiet_table <- list()
rixRix_table <- list()
ptypes <- matnut$ptypes
encoded <- matnut$encoded

for (i in 1:length(ptypes)){
  ### contrasts ###
  dietTab <- contrasts.getDecodedSummary(testPhen[[matnut$ptypes[i]]], "Diet")
  dietDiet_table[[matnut$ptypes[i]]] <- cbind(phenotype=rep(matnut$ptypes[i],nrow(dietTab)), dietTab)
  #rixTab <- contrasts.getDecodedSummary(testPhen[[ptypes[i]]], "RIX")
  #rixRix_table[[ptypes[i]]] <- cbind(phenotype=rep(ptypes[i],nrow(rixTab)), rixTab)
}

combineContrasts <- do.call(rbind, dietDiet_table)
contrasts_PO <- combineContrasts[which(combineContrasts$PO == "PO"),]
contrasts_noPO <- combineContrasts[-which(combineContrasts$PO == "PO"),]

rix <- encoded$Level[which(encoded$Variable == "RIX")]
diet <- encoded$Level[which(encoded$Variable == "Diet")]

summary_PO <- matrix(0, nrow=length(ptypes), ncol=length(rix))
colnames(summary_PO) <- paste0("RIX",rix)
rownames(summary_PO) <- ptypes
summary_noPO <- matrix(0, nrow=length(ptypes), ncol=length(rix))
colnames(summary_noPO) <- paste0("RIX",rix)
rownames(summary_noPO) <- ptypes


for (i in 1:length(ptypes)){
  temp <- contrasts_PO[which(contrasts_PO$phenotype == ptypes[i] & as.character(contrasts_PO$RIX) != ""),]
  tempNO <- contrasts_noPO[which(contrasts_noPO$phenotype == ptypes[i] & as.character(contrasts_noPO$RIX) != ""),]
  whichRIX <- paste0("RIX",temp$RIX)
  for(j in 1:nrow(temp)){
    if(temp$value[j] < 0.05){
      summary_PO[i,whichRIX[j]] <- summary_PO[i,whichRIX[j]]+ 1 
    }
  }
  whichRIX_no <- paste0("RIX",tempNO$RIX)
  for(j in 1:nrow(tempNO)){
    if(tempNO$value[j] < 0.05){
      summary_noPO[i,whichRIX[j]] <- summary_noPO[i,whichRIX_no[j]]+ 1 
    }
  }
}

summary_noPO <- rbind(summary_noPO, Total=apply(summary_noPO, 2, sum))
summary_PO <- rbind(summary_PO, Total=apply(summary_PO, 2, sum))
colnames(summary_PO) <- paste0("PO_", colnames(summary_PO))


#write.csv(cbind(summary_noPO, summary_PO), file.path("./","matnut_outputs/",'signif_counts.csv'), row.names = T)
#write.csv(do.call(rbind, dietDiet_table), file.path("./","matnut_outputs/",'allPheno_dietDiet_table.csv'), row.names = F)
#write.csv(do.call(rbind, rixRix_table), file.path("./","matnut_outputs/",'allPheno_rixRix_table.csv'), row.names = F)

###########
#  Plots  #
###########


### Flag Plots ###

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


### Mega-plot ###

keep <- c("Diet", "DietRIX","PODietRIX","PORIX", "RIX")
myMCMC <- testPhen$WeightPND60$mcmcObject[, which(colnames(testPhen$WeightPND60$mcmcObject) %in% 
                                                    encoded$Level[which(encoded$Variable %in% keep)])]
ribbon_plot(myMCMC, ptypes, encoded)


#ptypes <- names(psychPhenResults)
temp_rib<- list()
for(i in 1:5){#length(ptypes)){
  temp_rib[[ptypes[i]]]$predict <- prediction(mcmcOb=testPhen[[ptypes[i]]]$mcmcObject, ptypes=ptypes[i], encoded, Match=F)
  temp_rib[[ptypes[i]]]$rib <- ribbon_plot(mcmcOb=temp_rib[[ptypes[i]]]$predict, ptypes=ptypes[i], encoded)
}

all_phenDF <- data.frame()

for(i in 1:length(ptypes)){
  tempDF <- data.frame(phen = rep(names(temp_rib)[i], nrow(temp_rib[[names(temp_rib)[i]]]$savedata[[names(temp_rib)[i]]])), 
                       temp_rib[[names(temp_rib)[i]]]$savedata[[names(temp_rib)[i]]])
  all_phenDF <- rbind(all_phenDF, tempDF)
}

allPlot <- ggplot(all_phenDF, aes(color=rix)) + 
  geom_line(aes(y=mu, x=diet, group=poe, alpha=poe)) +  
  geom_ribbon(color=NA, aes(x=diet, ymin=lower.1, ymax=upper.1, group=poe, alpha=poe, 
                            fill=rix, color=rix)) + 
  geom_hline(yintercept = 0, colour="lightsteelblue", size=0.5, linetype="longdash") +
  scale_alpha_discrete(range = rev(c(0.45, 0.75)), name="POE") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_fill_discrete(guide=FALSE) + scale_color_discrete(guide=FALSE) + 
  facet_grid(rix~phen) + theme_minimal() + 
  ylab("Predicted effects") + xlab("Diet") #+ ggtitle(ptypes[j])


### More plots ###

graphics.off()
pdf(file.path("./","matnut_outputs/","allPheno_DietContrastPlots.pdf"), 
    onefile=TRUE, width = 8.5, height = 11)
for(i in 1:length(ptypes)){  
  print(contrastPlots[[i]])
}
graphics.off()

graphics.off()
pdf(file.path("./","matnut_outputs/","allPheno_JLComparison_1jun.pdf"), 
    onefile=TRUE, width = 8.5, height = 11)
for(i in 1:length(ptypes)){  
  print(myPhen[[ptypes[i]]]$plot)
}
graphics.off() 

graphics.off()
pdf(file.path("./","matnut_outputs/","matchedS2_12jul2017.pdf"), onefile=TRUE, width = 14, height = 11)
for(i in 1:length(matnut$ptypes)){  
  print(multimp2[[matnut$ptypes[i]]]$ribPlot$ribbon)
  #print(ribbon4[[ptypes[i]]])
}
graphics.off()

graphics.off()
pdf(file.path("./","matnut_outputs/","allPheno_DietOnlyribbon.pdf"), onefile=TRUE, width = 14, height = 11)
for(i in 1:length(ptypes)){  
  print(ribbon9_Dietonly[[ptypes[i]]]$ribbon9)
}
graphics.off()

### Matched ribbon ###
graphics.off()
pdf(file.path("./","matnut_outputs/","matchedribbon_MIMP.pdf"), onefile=TRUE, width = 14, height = 11)
for(i in 1:length(ptypes)){
  print(multimp[[ptypes[i]]]$plot$ribbon)
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
  shortTable <- mergeTable[which(mergeTable$Variable %in% c("Diet","RIX","DietRIX","PORIX","PODietRIX")),]
  
  summary_table[[pheno]] <- cbind(phenotype=rep(pheno,nrow(shortTable)), 
                                  shortTable[order(shortTable$Variable, shortTable$Level),])
  ####
  #dietTab <- contrasts.getDecodedSummary(psychPhen$OFTotalDistance, "Diet")
  #dietDiet_table[[pheno]] <- cbind(phenotype=rep(pheno,nrow(dietTab)), dietTab)
  #rixTab <- contrasts.getDecodedSummary(psychPhen$OFTotalDistance, "RIX")
  #rixRix_table[[pheno]] <- cbind(phenotype=rep(pheno,nrow(rixTab)), rixTab)
}


write.csv(do.call(rbind, summary_table), file.path("./","matnut_outputs/",'allPheno_summarytable.csv'), row.names = F)


#########
# Fixed #
#########

storeAn <- list()
#p_vals <- matrix(NA, nrow=length(ptypes), ncol=5, dimnames=list(ptypes))
#colnames(p_vals) <- c("RIX","Diet","Diet:RIX", "PO:RIX", "PO:DietRIX")

rix <- encoded$Level[which(encoded$Variable == "RIX")]
porix <- paste0("RIX",rix,":PO")
podietrix <- paste0("PO:DietRIX",unlist(strsplit(as.character(encoded$Level[which(encoded$Variable == "PODietRIX")]),'PO'))[c(F,T)])
p_vals <- matrix(NA, nrow=length(ptypes), ncol=(length(podietrix) + length(porix)), dimnames=list(ptypes))
colnames(p_vals) <- c(porix, podietrix)

for(i in 1:length(ptypes)){
  pheno <- paste(ptypes[i])
  use <- which(is.na(orderdat[,pheno]) == F)
  matnut_use <- orderdat[use,]
  
  OF <- ifelse(length(grep("OF",pheno))==0,F,T)
  LD <- ifelse(length(grep("LD",pheno))==0,F,T)
  SIH <- ifelse(length(grep("SIH",pheno))==0,F,T)
  FST <- ifelse(length(grep("FST",pheno))==0,F,T)
  Stress <- ifelse(length(grep("CORT",pheno))==0,F,T)
  
  if(OF){
    indvariable <- "~ -1 + (1 | DamID) + (1 | SireID) + (1 | BehaviorBatch) + OFBox + RIX + Diet + Diet:RIX + (0+PO|RIX) + (0+PO|Diet:RIX)" 
  } else if(LD){
    indvariable <- "~ -1 + (1 | DamID) + (1 | SireID) + (1 | BehaviorBatch) + LDChamber + RIX + Diet + Diet:RIX + (0+PO|RIX) + (0+PO|Diet:RIX)" 
  } else if(SIH){
    indvariable <- "~ -1 + (1 | DamID) + (1 | SireID) + (1 | BehaviorBatch) + SIHOrder + RIX + Diet + Diet:RIX + (0+PO|RIX) + (0+PO|Diet:RIX)" 
  } else if(FST){
    indvariable <- "~ -1 + (1 | DamID) + (1 | SireID) + (1 | BehaviorBatch) + FSTChamber + RIX + Diet + Diet:RIX + (0+PO|RIX) + (0+PO|Diet:RIX)" 
  } else if(Stress){
    indvariable <- "~ -1 + (1 | DamID) + (1 | SireID) + (1 | BehaviorBatch) + RestraintExperimenter + RestraintOrder + 
    RIX + Diet + Diet:RIX + (0+PO|RIX) + (0+PO|Diet:RIX)" 
  } else {
    indvariable <- "~ -1 + (1 | DamID) + (1 | SireID) + (1 | BehaviorBatch) + RIX + Diet + Diet:RIX + (0+PO|RIX) + (0+PO|Diet:RIX)" 
  }
  
  form <- formula(paste0(pheno, indvariable))
  lmObj <- lmer(form, data=matnut_use)
  ano <- anova(lmObj, type=1)
  poef <- fixef(lmObj)
  for(j in 1:ncol(p_vals)){
    val <- ifelse(is.element(colnames(p_vals)[j], names(poef)), poef[which(names(poef) %in% colnames(p_vals)[j])], NA)
    p_vals[i,j] <- val
  }
  #p_vals[i,] <- ano$'Pr(>F)'[c((length(ano$'Pr(>F)')-4):length(ano$'Pr(>F)'))]
  
  storeAn[[i]] <- lmObj
}

################
# Prepare data #
################
matnut <- read.csv("../data/AllPhenotypes_Matnut5.csv")

matnut[grep("a", matnut$Reciprocal),"PO"] <- 0.5
matnut[grep("b", matnut$Reciprocal),"PO"] <- -0.5
matnut[which(matnut$PO == 0.5), "pof"] <- "+"
matnut[which(matnut$PO == -0.5), "pof"] <- "-"
dietlabs <- c("STD", "ME", "PD", "VDD")
matnut$RIX <- factor(matnut$RIX, levels=c(1:4,6:10))
matnut$Diet <- factor(matnut$Diet, levels=dietlabs)


### covariates ###

matnut$BBoriginal <- matnut$BreedingBatch
matnut$BreedingBatch <- factor(paste0("br",matnut$BBoriginal), levels=unique(paste0("br",matnut$BBoriginal)))
matnut$BEoriginal <- matnut$BehaviorBatch
matnut$BehaviorBatch <- factor(paste0("be",matnut$BEoriginal), levels=unique(paste0("be",matnut$BEoriginal)))
matnut$Cage <- factor(matnut$Cage)
matnut$DamID <- factor(paste0("d",matnut$DamID))
matnut$SireID <- factor(paste0("s",matnut$SireID))

matnut$OFBox <- factor(paste0("of",matnut$OFBox))
matnut$LDChamber <- factor(paste0("ld",matnut$LDChamber))
matnut$SIHOrder <- factor(paste0("sih",matnut$SIHOrder))
matnut$FSTChamber <- factor(paste0("fst",matnut$FSTChamber))
matnut$RestraintOrder <- factor(paste0("ro",matnut$RestraintOrder))
matnut$RestraintExperimenter <- factor(gsub(" ","",matnut$RestraintExperimenter))

orderdat <- matnut[order(matnut$RIX, matnut$Diet, matnut$PO),]
orderdat$PORIX <- paste0("PO",orderdat$RIX)
orderdat$PORIX <- factor(orderdat$PORIX, levels=unique(orderdat$PORIX))
orderdat$DietRIX <- paste0(orderdat$Diet, orderdat$RIX)
orderdat$DietRIX <- factor(orderdat$DietRIX, levels=unique(orderdat$DietRIX))
orderdat$PODietRIX <- paste0("PO",orderdat$DietRIX)
orderdat$PODietRIX <- factor(orderdat$PODietRIX, levels=unique(orderdat$PODietRIX))

variables <- c("Diet", "RIX", "PORIX", "DietRIX", "PODietRIX",
               "BehaviorBatch","DamID","SireID",
               "OFBox","LDChamber","SIHOrder","FSTChamber","RestraintOrder","RestraintExperimenter")
dispersions <- c("tauInv2","sigInv2", "tauPEInv2","tauDRInv2", "tauPDRInv2", 
                 "tauBEInv2","tauDamInv2","tauSirInv2")
#"tauOFInv2","tauLDInv2","tauSIHInv2","tauFSTInv2","tauROInv2","tauREInv2")

encoded <- getEncoding(orderdat, variables)
#write.csv(encoded, file.path("./","matnut_outputs/",'encoding.csv'), row.names = F)


ptypes <- colnames(matnut)[44:63]
#ptypes <- colnames(matnut)[c(44,48,53,54,55,61,62,63)]
#ptypes <- c("OFTotalDistance","BasalCORT","StressCORT")

allparam <- c(variables, dispersions)
allparam <- allparam[order(allparam)]
matnut <- list(df=orderdat, parameters=allparam, ptypes=ptypes, encoded=encoded)

saveRDS(matnut, file.path("./","matnut_outputs/",'matnut_data.rds'))
# many of the following functions expect list of data with 1) data, 2) parameters, 3) phenotypes, 4) encoding


##############
#  Sim data  #
##############

testdata <- getTestCase()
resultCompare <- testdata$effects[-c(which(testdata$effects$coef == 0), grep("tau", testdata$effects$variable)),]

firstDiet <- c()
for(i in 1:length(which(resultCompare$variable == "Diet"))){
  firstDiet[i] <- switch(as.numeric(resultCompare$level[i]), "AD", "BD", "CD", "DD")
}
firstDiet <- c(firstDiet, rep(1:9, 2))
use <- grep("[.]", resultCompare$level)
diet <- as.numeric(unlist(strsplit(resultCompare$level[use], "[.]"))[c(F, T)])
rix <- unlist(strsplit(resultCompare$level[use], "[.]"))[c(T, F)]
dietLet <- c()
for(i in 1:length(use)){
  dietLet[i] <- c(switch(diet[i], "AD", "BD", "CD", "DD"))
}

dietLet <- c(firstDiet, paste0(dietLet, rix))
labelNew <- ifelse(1:94 %in% grep("PO", resultCompare$variable), paste0("PO", dietLet), dietLet)

resultCompare$Level <- labelNew

mergeTab <- merge(resultCompare, testPhen$JL_compare, by="Level")
## compares real answers with lmer and jags estimates


testdataDF <- data.frame(testdata$cov.data)
testdataDF$DietOld <- testdataDF$Diet
for(i in 1:nrow(testdataDF)){
  testdataDF$Diet_use[i] <- switch(as.numeric(testdataDF$Diet)[i], "AD", "BD", "CD", "DD")
}
testdataDF$DamID <- factor(paste0("dd", testdataDF$DamID),levels=unique(paste0("dd", testdataDF$DamID)))
testdataDF$Diet <- as.factor(testdataDF$Diet_use)
testdataDF$DietRIX <- as.factor(paste0(testdataDF$Diet, testdataDF$RIX))
testdataDF$PORIX <- as.factor(paste0("PO", testdataDF$RIX))
testdataDF$PODietRIX <- as.factor(paste0("PO", testdataDF$DietRIX))
encodeTest <- getEncoding(testdataDF, c("DamID", "Diet", "RIX", "DietRIX", "PORIX", "PODietRIX")) #colnames(testdataDF)[c(1:3, 11:13)])
testPhen <- makeSummary(datalist=testdataDF, phenotype = "outcome", tryLam=c(1), normd=F, 
                        chains=2, n.adapt=20000, n.iter=50000, thin=10, encoded=encodeTest, sq=F)

predTest <- testdataDF[,c(7, 2, 3, 6, 8, 9)]

for(i in 1:nrow(predTest)){
  predTest$jagsDiet[i] <- testPhen$JL_compare$JAGS_est[which(testPhen$JL_compare$Level == paste(predTest$Diet[i]))]
  predTest$jagsRIX[i] <- testPhen$JL_compare$JAGS_est[which(testPhen$JL_compare$Level == paste(predTest$RIX[i]))]
  predTest$lmerDiet[i] <- testPhen$JL_compare$LMER_est[which(testPhen$JL_compare$Level == paste(predTest$Diet[i]))]
  predTest$lmerRIX[i] <- testPhen$JL_compare$LMER_est[which(testPhen$JL_compare$Level == paste(predTest$RIX[i]))]
}

predTest$jagsEps <- predTest$outcome - (predTest$jagsDiet + predTest$jagsRIX)
predTest$lmerEps <- predTest$outcome - (predTest$lmerDiet + predTest$lmerRIX)
## sum of effects and calculation of epsilon from estimates

jagsMu <- c()
for(i in 1:ncol(testPhen$mcmcObject)){
  jagsMu[i] <- density(testPhen$mcmcObject[,i])$x[which.max(density(testPhen$mcmcObject[,i])$y)]
}
source(file.path(".", "matnut", "summary_functions.R"))


