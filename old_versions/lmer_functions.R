source("./lm/formulaWrapper.R")
source("./lm/fitBoxCoxModels.R")

src <- ifelse(length(grep("C:", getwd())) >0, file.path("C:","Users","kathie","Documents","source_ks"), 
              ifelse(length(grep("kys6", getwd())) >0, file.path("~","source_ks/matnut_src"), file.path("~","Documents","matnut_src")))

source(file.path(src, "boxcox_functions.R"))


lmer.getDecodedSummary <- function(lmerOb, phenotype, formula=NA){
  lmerSummary <- data.frame(Intercept = numeric(), Level=factor(), Variable=factor())
  if(is.na(formula)){
    formula <- attributes(lmerOb[[phenotype]]$phen_1$fit@frame)$formula[3]
  }
  form <- formulaWrapper$parseCovariateString(paste(formula))
    #lmerOb[[phenotype]]$formula)
  if(any(duplicated(form$fixef))){
    form$fixef[[which(duplicated(form$fixef))]] <- NULL
  }
  
  lmerObject <- lmerOb[[phenotype]]$phen_1$fit
  if(class(lmerObject) == "lm"){
    vec <- summary(lmerObject)$coefficients[,"Estimate"]
  } else {
    vec <- fixef(lmerObject)
  }

  no <- c("0", "1", "-1")
  for(i in 1: length(form$fixef)){
    if(length(intersect(form$fixef[[i]], no)) == 0){
      levels <- names(vec)[grep(paste0(gsub("^form$fixef[[i]][*]\\s","",form$fixef[[i]]), collapse=""), names(vec))]
      levels2 <- gsub(paste0(gsub("[*]\\s","",form$fixef[[i]]), collapse=""),"",levels)
      if (length(levels)>0){
        for (j in 1: length(levels)){
          this_rand <- vec[j]
          tempSum <- data.frame(Intercept = as.numeric(vec[j]), 
                                Level = levels2[j], 
                                Variable = paste0(gsub("[*]","",form$fixef[[i]]), collapse=""))
          lmerSummary <- rbind(lmerSummary, tempSum)
        }
      }
    }
  }
  if (length(form$ranef) > 0){
    randEfs <- list()
    for(j in 1:length(names(ranef(lmerObject)))){
      tempName <- names(ranef(lmerObject))[j]
      randEfs[[gsub("`","",tempName)]] <- ranef(lmerObject)[[tempName]]
      colnames(randEfs[[tempName]])[which(colnames(randEfs[[tempName]]) == "(Intercept)")] <- ""
      for (k in 1:ncol(randEfs[[tempName]])){
        addToSum <- data.frame(Intercept = as.numeric(randEfs[[j]][,k]), 
                               Level = as.character(paste0(gsub("^s", "", colnames(randEfs[[j]])[k]), 
                                                           gsub(":","",rownames(randEfs[[j]])))), 
                               Variable = as.character(paste0(gsub("^s", "", colnames(randEfs[[j]])[k]), 
                                                              gsub(":","",tempName))))
        lmerSummary <- rbind(lmerSummary, addToSum)
      }
    }
  }
  return(lmerSummary)
}



runLMERModels_cov <- function(phenotype, dataf, tryLam=1, randvar=NA, fixvar=NA, POvar=NA, 
                              normd=T, randDam=F, Match=F){
  # put in another separate file
  effectlist <- list()
  ind <- 1
  
  if(!is.na(fixvar)){
    effectlist[[ind]] <- vector()
    for(i in 1:length(fixvar)){
      effectlist[[ind]] <- c(effectlist[[ind]], fixvar[i])
    }
    ind=ind+1
  }
  if(!is.na(randvar)){
    effectlist[[ind]] <- vector()
    for(i in 1:length(randvar)){
      effectlist[[ind]] <- c(effectlist[[ind]], paste0("(1 | ", randvar[i], ")"))
    }
    ind=ind+1
  }
  
  if(!is.na(POvar)){   
    for(i in 1:length(POvar)){
    effectlist[[ind]] <- vector()
      effectlist[[ind]] <- c(effectlist[[ind]], paste0("(PO + 0 | ", POvar[i], ")"))
    }
  }
  
  fitLMER <- list()
#browser()
  #dt <- "Diet"
  #rixef <- "(1 | RIX)"
  #dietRix <- "(1 | RIX) + (1 | DietRIX)"
  #dietRixS <- "(s | RIX) + (s | DietRIX)"
  #randEf <- ifelse(randDam, "(1 | DamID)", "")      # + (1 | SireID) + (1 | BehaviorBatch)"
  #podietrix <- "(1 | RIX) + (PO + 0 | RIX) + (1 | Diet:RIX) + (PO + 0 | Diet:RIX)"
  #JAGSpodietrix <- "(PO + 1 | RIX) +  (PO + 1 | Diet:RIX)"
  #if(length(grep("OF",phenotype))>=1){
  #  fixEf="OFBox"
  #} else if(length(grep("LD",phenotype))>=1){
  #  fixEf="LDChamber"
  #} else if(length(grep("SIH",phenotype))>=1){
  #  fixEf="SIHOrder"
  #} else if(length(grep("FST",phenotype))>=1){
  #  fixEf="FSTChamber"
  #} else if(length(grep("CORT",phenotype))>=1){
  #  fixEf=paste("RestraintExperimenter", "RestraintOrder", sep="+")  
  #} else {
  #  fixEf=NA  
  #}
  #if(Match){
  #  indvariable <- paste(init, dietRix, sep="+")
  #  lmerForm <- paste("~ 1 + ", dietRix)      #ifelse(MatchS, paste("~ 1 + ", dietRixS), 
  #} else if (is.na(fixEf)){
  #  indvariable <- paste(init, dt, podietrix, sep="+")            #randEf, 
  #  init2 <- paste("~", dt)
  #  lmerForm <- paste(init2, JAGSpodietrix, sep="+")              #randEf,
  #} else {
  #  indvariable <- paste(init, dt, fixEf, podietrix, sep="+")     #randEf,
  #  init2 <- paste("~", dt)
  #  lmerForm <- paste(init2, fixEf, JAGSpodietrix, sep="+")       #randEf,
  #}

  lmerform <- paste("~ -1 +", paste(unlist(effectlist), collapse="+"))
  jagsform <- gsub("-1","1",lmerform)
  jagsform <- gsub("[+] 0", "", jagsform)
  
  #for (i in 1:length(phenotype)){ 
    #pheno <- paste(phenotype[i])
    #use <- which(is.na(dataf[,phenotype]) == F)
    use <- which(colSums(dataf[,phenotype]) != 0)
    phenotype <- phenotype[use]
    y.mat <- dataf[complete.cases(dataf), phenotype]
    
    #matnut_use <- matnut_use[,use2]
    #fitLMER[[phenotype[i]]]$bc <- BC.model(y.mat = matnut_use[,pheno], matnut_use, indvariable, 
    #                                       transformParams = getMatnutTransformParams(tryLam = tryLam))
    
    #if(Match || (tryLam == 1 && normd == F)){
    #  matnut_use$y.trans <- ifelse(rep(normd, nrow(matnut_use)), scale(matnut_use[,pheno]), matnut_use[,pheno])
      #form <- formula(paste("phen_use", indvariable))
    #  fitLMER[[phenotype[i]]]$phen_1 <- list()
    #  fitLMER[[phenotype[i]]]$phen_1$y.transformed <- as.vector(matnut_use$y.trans)
    #} else {
  browser()
      fitLMER <- fit.model.bc$fit(y.mat = y.mat, dataf, lmerform, 
                                                  transformParams = getMatnutTransformParams(tryLam = tryLam, normd = normd))
    #}
    #print(paste(pheno, "fit by LMER"))
    #matnut_use <- cbind(matnut_use, y.trans=fitLMER[[phenotype[i]]]$phen_1$y.transformed)
    #form <- as.formula(paste0("y.trans", lmerform))
    #if (is.na(randvar)){
    #  lobj <- lm(form, data=matnut_use)
    #} else {
    #  lobj <- lmer(form, data=matnut_use)
    #}
    
    #fitLMER[[phenotype[i]]]$phen_1$fit <- lobj
    fitLMER$formula <- lmerform
  #}
  return(fitLMER)
}

