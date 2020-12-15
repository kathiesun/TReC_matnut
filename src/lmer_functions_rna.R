source("formulaWrapper.R")
source("boxcox_functions.R")
#source("./loadParams.R")
#source("./parallel/accumulator.R")

require(car)

lmer.getDecodedSummary <- function(lmerOb, formula, phenotypes=NA){

  form <- formulaWrapper$parseCovariateString(paste(formula))
  
  if(any(duplicated(form$fixef))){
    form$fixef[[which(duplicated(form$fixef))]] <- NULL
  }
  
  if (all(is.na(phenotypes))){
    phenotypes <- lmerOb$lmerobj$phenotypes
  }
  
  lmerSummary <- list()
  for(p in 1:length(phenotypes)){
    
    lmerObject <- lmerOb$lmerobj$fits[[p]]
 
    if (class(lmerObject) == "character"){
      lmerSummary[[phenotypes[p]]] <- NA
    } else {
      
      lmerSummary[[phenotypes[p]]] <- data.frame(Intercept = numeric(), Level=factor(), Variable=factor())
      
      if(class(lmerObject) == "lm"){
        vec <- as.data.frame(summary(lmerObject)$coefficients)$Estimate
        names(vec) <- rownames(summary(lmerObject)$coefficients)
      } else {  
        vec <- fixef(lmerObject)
      }
      no <- c("0", "1", "-1")
      for(i in 1: length(form$fixef)){
        if(length(intersect(form$fixef[[i]], no)) == 0){
          levels <- names(vec)[grep(paste0(gsub("^form$fixef[[i]][*]\\s","",form$fixef[[i]]), collapse=""), names(vec))]
          if(length(levels) == 1){
            levels2 <- levels
          } else {
            levels2 <- gsub(paste0(gsub("[*]\\s","",form$fixef[[i]]), collapse=""),"",levels)
          }
          if (length(levels2)>0){
            for (j in 1: length(levels2)){
              this_rand <- vec[j]
              tempSum <- data.frame(Intercept = as.numeric(vec[grep(levels2[j], names(vec))]), 
                                    Level = levels2[j], 
                                    Variable = paste0(gsub("[*]","",form$fixef[[i]]), collapse=""))
              lmerSummary[[phenotypes[p]]] <- rbind(lmerSummary[[phenotypes[p]]], tempSum)
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
            lmerSummary[[phenotypes[p]]] <- rbind(lmerSummary[[phenotypes[p]]], addToSum)
          }
        }
      }
    }
  }
  
  return(as.list(lmerSummary))
}



runLMERModels_cov <- function(phenotype, dataf, tryLam=1, randvar=NA, fixvar=NA, POvar=NA, 
                              normd=T, randDam=F, Match=F){
  effectlist <- list()
  ind <- 1
  
  if(!any(is.na(fixvar))){
    effectlist[[ind]] <- vector()
    for(i in 1:length(fixvar)){
      effectlist[[ind]] <- c(effectlist[[ind]], fixvar[i])
    }
    ind=ind+1
  }
  if(!any(is.na(randvar))){
    effectlist[[ind]] <- vector()
    for(i in 1:length(randvar)){
      effectlist[[ind]] <- c(effectlist[[ind]], paste0("(1 | ", randvar[i], ")"))
    }
    ind=ind+1
  }
  
  if(!any(is.na(POvar))){   
    effectlist[[ind]] <- vector()
    for(i in 1:length(POvar)){
      effectlist[[ind]] <- c(effectlist[[ind]], paste0("(PO + 0 | ", POvar[i], ")"))
    }
  }
  fitLMER <- list()
  lmerform <- paste("~ -1 +", paste(unlist(effectlist), collapse="+"))
  jagsform <- lmer.to.jags.form(lmerform)
  
  if(length(phenotype) >1) {
    use <- which(colSums(dataf[,colnames(dataf) %in% phenotype], na.rm = T) != 0)
    use <- unique(c(use, which(apply(dataf[,colnames(dataf) %in% phenotype], 2, var) != 0) ))
    phenotype <- phenotype[use]
  }
  y.mat <- data.frame(dataf[, phenotype])
  dataf <- data.frame(dataf)
  colnames(y.mat) <- phenotype
  
  fitLMER$lmerobj <- BC.model(y.mat = y.mat, data=dataf, indvariable=lmerform, 
                              transformParams = getMatnutTransformParams(tryLam = tryLam, normd = normd))
  
  fitLMER$LMformula <- lmerform
  fitLMER$JAGSformula <- jagsform
  return(fitLMER)
  
}

