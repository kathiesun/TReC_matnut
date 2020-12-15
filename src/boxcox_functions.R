source("formulaWrapper.R")
library(data.table)
library(lmerTest)
library(car)

BC.model <- function(y.mat, data, indvariable, transformParams=getMatnutTransformParams(), 
                     print=T, anova=T, offset=0.1, sameLambda=F){
  p_val <- c()
  y_temp <- list()
  store_fit <- list()
  store_phen <- list()
  anovaOut <- list()
  
  indx <- which(apply(as.matrix(y.mat), 2, sum) == 0)
  indx <- c(indx, which(apply(as.matrix(y.mat), 2, var) == 0))
  #indx <- c(indx, which(apply(y.mat, 2, function(x) any(is.na(x) | is.nan(x)))))
  y.orig <- y.mat
  
  if(length(indx)>0){
    y.mat <- y.mat[,-indx]
    data <- data[,-which(colnames(data) %in% names(indx))]
  }
  if(dim(y.mat)[2] > 1){
    temp.mat <- as.list(y.mat)
    y.mat <- lapply(temp.mat, function(x) x[!is.na(x)])
  } else {
    temp.mat <- y.mat
    y.mat <- list()
    y.mat[[1]] <- temp.mat
  }
  pvals <- matrix(NA, nrow=length(transformParams$lambdasToTry), ncol=length(y.mat))
  
  if(transformParams$normalizeBeforeTransform == T){
    y.mat <- lapply(y.mat, scale)
  }
  
  for(i in 1:length(transformParams$lambdasToTry)){
    lambda <- transformParams$lambdasToTry[i]
    y_temp[[i]] <- list()
    if (lambda == 0){
      add.vals <- ifelse(lapply(y.mat, function(x) min(x[!is.na(x)])) <= 0, as.list(abs(unlist(lapply(y.mat, min, na.rm=T))) + offset), 
                         as.list(rep(0, length(y.mat))))
      for(m in 1:length(y.mat)){
        y_temp[[i]][[m]] <- log(y.mat[[m]] + add.vals[[m]])
      }
      #temp <- sapply(1:length(y.mat), function(x) y.mat[[x]] + add.vals[[x]])
      #y_temp[[i]] <- lapply(temp, function(x) log(x[!is.na(x)]))
    } else {
      #lamb2 <- ifelse(lambda<1 && min(x) <=0, abs(min(x))+1, 0)
      if(lambda<1){
        add.vals <- ifelse(lapply(y.mat, function(x) min(x[!is.na(x)])) <= 0, 
                           as.list(abs(unlist(lapply(y.mat, function(x) min(x, na.rm = T)))) + offset), 
                           as.list(rep(0, length(y.mat))))
      } else if (lambda == 1){
        add.vals <- as.list(rep(1, length(y.mat)))
      } else {
        add.vals <- as.list(rep(0, length(y.mat)))
      }
      for(m in 1:length(y.mat)){
        y_temp[[i]][[m]] <- ((y.mat[[m]] + add.vals[[m]])^lambda - 1) / lambda
      }
      
      #temp <- sapply(1:length(y.mat), function(x) y.mat[[x]] + add.vals[[x]])
      #<- lapply(y.mat, function(x)
      #                     (x^lambda - 1) / lambda
      #                     )
    }
    if(length(indx > 0)){
      phen <- colnames(y.orig[,-indx])
    } else {
      phen <- colnames(y.orig)
    }
    
    names(y_temp[[i]]) <- phen
    if(transformParams$normalizeAfterTransform == T){
      y_temp[[i]] <- lapply(y_temp[[i]], scale)
    }

    #if(is.vector(y.orig) || dim(y.orig)[2] == 1) phen <- "y"
    

    #data_use <- cbind(data, y_temp[[i]])
    #y <- as.matrix(data_use[,which(colnames(data_use) %in% phen)])
    #x_lab <- unlist(formulaWrapper$parseCovariateString(paste(indvariable))$fixef)
    #x_lab <- c(x_lab, unlist(formulaWrapper$parseCovariateString(paste(indvariable))$ranef))
    #x <- data.matrix(data_use[,which(colnames(datause) %in% x_lab)])
    rand=T
    fit <- list()
    anovaOut[[i]] <- list()
    if(length(formulaWrapper$parseCovariateString(paste(indvariable))$ranef) == 0){
      rand = F
      data$PO <- as.factor(data$PO)
    }
    #grepl(pattern="\\([a-zA-Z0-9\\.]+\\|[a-zA-Z0-9\\.]+\\)", x=indvariable, perl=TRUE)
    
    for(j in 1:length(phen)){
      didItWork <- "YES"
      form <- as.formula(paste0(phen[j], indvariable))
      temp_data <- data[which(is.na(data[,phen[j]]) == F),-which(colnames(data) %in% phen)]
      use_data <- cbind(temp_data, y_temp[[i]][[j]][!is.na(y_temp[[i]][[j]])])
      
      colnames(use_data)[length(colnames(use_data))] <- phen[j]
      #complete.cases(data_use[,phen[j]])

      if (!rand){
        fit[[j]] <- tryCatch(lm(form, use_data), error=function(e) didItWork="NO")
      } else {
        fit[[j]] <- tryCatch(lmer(form, use_data), error=function(e) didItWork="NO")
      }
      
      didItWork <- tryCatch(Anova(fit[[j]], type="III"), error=function(e) didItWork="NO",
                            warning=function(w) didItWork="NO")
      
      if(anova) anovaOut[[i]][[j]] <- didItWork
        
      if("anova" %in% class(didItWork)){
        res <- residuals(fit[[j]])
        pvals[i,j] <- shapiro.test(res)$p.value
        printThis <- paste("everything worked:", i, "phen", j)
      } else {
        pvals[i,j] <- 0
        printThis <- "nope"
      }
      if(print) print(printThis)
    }

    store_fit[[i]] <- fit
    store_phen[[i]] <- phen
  }
  max_p <- apply(pvals, 2, which.max)

  if(sameLambda) {
    ind_max <- as.numeric(names(sort(-table(max_p)))[1])
    p_max <- pvals[ind_max,]
    lamb_max <- transformParams$lambdasToTry[ind_max]
  } else {
    ind_max <- apply(pvals, 2, which.max)    
    p_max <- sapply(1:ncol(pvals), function(x) pvals[ind_max[x],x])
    lamb_max <- sapply(1:ncol(pvals), function(x) transformParams$lambdasToTry[ind_max[x]])
  }
  
  y.transform <- list()
  fits <- list()
  phenotypes <- list()
  anovas <- list()
  for(m in 1:ncol(pvals)){
    y.transform[[m]] <- y_temp[[ind_max[m]]][[m]]
    fits[[m]] <- store_fit[[ind_max[m]]][[m]]
    anovas[[m]] <- anovaOut[[ind_max[m]]][[m]]
    #phenotypes <- store_phen[[ind_max[m]]][[m]]
  }
  
  #pvals <- pvals[ind,]
  #temp <- sapply(1:length(y.mat), function(x) y.mat[[x]] + add.vals[[x]])
  
  return(list(p.val=p_max, lambda=lamb_max, y.transform=y.transform, 
              fits=fits, anovaOut = anovas, phenotypes=phen))
}


getMatnutTransformParams <- function(normd=T, tryLam=1)
{
  out = list()
  out$lambdasToTry = tryLam
  out$lambdaPerY      = NULL
  if (normd){
    out$normalizeBeforeTransform = T
    out$normalizeAfterTransform  = T
  }
  else {
    out$normalizeBeforeTransform = F
    out$normalizeAfterTransform  = F
  }
  out$extremelb      = -3
  out$extremeub      = 3
  return(out)
}

