source("boxcox_functions.R")

stanSum <- function(df, phenotype, encoded=NULL, 
                    tryLam=1, normd=T,
                    chains=2, iter=10000,
                    plot=T, contrasts=F,
                    randvar=NA, fixvar=NA, POvar=NA){
  
  formulas = getFormulas(fixvar, randvar, POvar)
  
  if(length(phenotype) > 1) {
    use <- which(colSums(df[,colnames(df) %in% phenotype], na.rm = T) != 0)
    use <- unique(c(use, which(apply(df[,colnames(df) %in% phenotype], 2, var) != 0) ))
    phenotype <- phenotype[use]
  }
  y.mat <- data.frame(df[, phenotype])
  df <- data.frame(df)
  phenotype = ifelse(length(grep("^[0-9]",phenotype)>0), paste0("X",phenotype), phenotype)
  colnames(y.mat) <- phenotype
  
  bcObject <- BC.model(y.mat = y.mat, data=df, indvariable=formulas$lmerform, 
                       transformParams = getMatnutTransformParams(tryLam = tryLam, normd = normd))
  transf <- data.frame(lambda = bcObject$lambda,  
                       pval = bcObject$p.val)
  
  form <- formulaWrapper$parseCovariateString(paste(formulas$jagsform))
  
  if(any(duplicated(form$fixef))){
    form$fixef[[which(duplicated(form$fixef))]] <- NULL
  }
  
  y.mat = bcObject$y.transform[[1]]
  y <- as.vector(y.mat[!is.na(y.mat)])
  N <- length(y)
  if(length(which(is.na(y.mat))) > 0) df = df[-which(is.na(y.mat)),]
  x_fx = data.frame(df[,fixvar])
  colnames(x_fx) = fixvar
  x_rd = data.frame(df[,randvar])
  colnames(x_rd) = randvar
  x_fx = do.call("cbind", 
                 sapply(1:ncol(x_fx), function(i) encoded$Index[match(x_fx[,i], encoded$Level)], 
                        simplify=F))
  colnames(x_fx) = fixvar
  
  x_rd = do.call("cbind", 
                 sapply(1:ncol(x_rd), function(i) encoded$Index[match(x_rd[,i], encoded$Level)], 
                        simplify=F))
  colnames(x_rd) = randvar
  
  RIX = df$RIX
  Diet = df$Diet
  DietRIX = df$DietRIX
  
  
  contrasts(RIX) = contr.sum(length(unique(df$RIX)))
  contrasts(Diet) = contr.sum(length(unique(df$Diet)))
  contrasts(DietRIX) = contr.sum(length(unique(df$DietRIX)))
  
  modelMat = model.matrix(~ 0 + Diet + RIX, df)     # + DietRIX
  modelMat = modelMat[,-1]
  x_d  = modelMat[,grep("[^0-9]$", colnames(modelMat))]
  x_s  = modelMat[,grep("^RIX", colnames(modelMat))]
  #x_sd = modelMat[,grep("^DietRIX", colnames(modelMat))]
  x_d = model.matrix(~ 0 + Diet, df)
  x_s = model.matrix(~ 0 + RIX, df)
  
  standat <-  list(N       = length(y),
                   y       = y, 
                   K_d     = ncol(x_d),
                   K_s     = ncol(x_s),
                   #K_sd    = ncol(x_sd),
                   x_d     = x_d,
                   x_s     = x_s,
                   #x_sd    = x_sd,
                   SPO     = as.vector(df$PO))
  
  fileName <- "../stan_SPO.stan"  
  fileName <- "../stan_SPO_fullMat.stan"

  stan_code <- readChar(fileName, file.info(fileName)$size)
  #cat(stan_code)
  
  
  warmup=round((floor((iter*0.25)/(iter/10))*(iter/10)))
  thin=iter/1000
  
  stan_code <- readChar(fileName, file.info(fileName)$size)
  try(resStan <- stan(model_code = stan_code, data = standat,
                      chains = chains, iter = iter, warmup = warmup, thin = thin))
  
  return(resStan)
}



getFormulas = function(fixvar, randvar, POvar){
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
  
  return(list(lmerform=lmerform, jagsform=jagsform))
}



