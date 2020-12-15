
contrasts.getDecodedSummary <- function(summaryObject, compare=NULL, narrow=0.5, wide=0.95){ 
  # remove redundancy
  testThis <- ifelse(compare=="Diet", "diet.test", "rix.test")
  mcmcObject <- as.mcmc(do.call(cbind,summaryObject$contrastHPD))
  tempPval <- summaryObject$pval
  colnames(tempPval) <- paste0("a",colnames(summaryObject$pval))
  tempsignif <- summaryObject$signif
  colnames(tempsignif) <- paste0("a",colnames(summaryObject$signif))
  
  pval <- cbind(melt(tempPval), signif=melt(tempsignif)[,"value"])
  pval$Var2 <- as.character(gsub("a","",pval$Var2))
  
  mu    <- c(colMeans(mcmcObject))
  med   <- apply(coda::HPDinterval(mcmcObject, prob=0.01), 1, median)
  hpd.wide    <- coda::HPDinterval(mcmcObject, prob=wide)
  hpd.narrow  <- coda::HPDinterval(mcmcObject, prob=narrow)
  effectdata <- data.frame(cbind(mu,med,hpd.narrow,hpd.wide))
  effectdata$Level1 <- unlist(strsplit(rownames(effectdata),"[.]"))[c(T, F)]
  effectdata$Level2 <- unlist(strsplit(rownames(effectdata),"[.]"))[c(F, T)]
  if(compare=="Diet"){
    pval$PO[grep("PO",pval$Var1)] <- "PO"
    pval$RIX <- gsub('\\D+','', pval$Var1)
    pval$Level1 <- unlist(strsplit(as.character(pval$Var2),"[.]"))[c(T, F)]
    pval$Level2 <- unlist(strsplit(as.character(pval$Var2),"[.]"))[c(F, T)]
    effectdata$RIX <- gsub('\\D+','', effectdata$Level1)
    effectdata$Level1 <- gsub('[0-9]','', effectdata$Level1)
    effectdata$Level2 <- gsub('[0-9]','', effectdata$Level2)
    effectdata$PO[grep("PO",effectdata$Level1)] <- "PO"
    effectdata$Level1 <- gsub('PO','', effectdata$Level1)
    effectdata$Level2 <- gsub('PO','', effectdata$Level2)
    
    contrastSum <- merge(effectdata, pval[,c("PO","RIX","Level1","Level2","value","signif")], by=c("PO","RIX","Level1","Level2"))
    contrastSum$RIX <- factor(contrastSum$RIX, levels=unique(effectdata$RIX))
    
    contrastSum$Level1 <- factor(contrastSum$Level1, levels=unique(pval$Level1))
    contrastSum$Level2 <- factor(contrastSum$Level2, levels=unique(pval$Level2))
    contrastSum <- contrastSum[order(contrastSum$PO, contrastSum$RIX, contrastSum$Level1, contrastSum$Level2),]
    
  } else {
    pval$PO[grep("PO",pval$Var1)] <- "PO"
    pval$Diet <- gsub('PO','', pval$Var1)
    pval$Level1 <- unlist(strsplit(as.character(pval$Var2),"[.]"))[c(T, F)]
    pval$Level2 <- unlist(strsplit(as.character(pval$Var2),"[.]"))[c(F, T)]
    pval[order(pval$PO, pval$Level1, pval$Level2),]
    effectdata$Diet <- gsub('\\d+','', effectdata$Level1)
    effectdata$PO[grep("PO",effectdata$Diet)] <- "PO"
    effectdata$Diet <- gsub('PO','', effectdata$Diet)
    effectdata$Level1 <- gsub('\\D+','', effectdata$Level1)
    effectdata$Level2 <- gsub('\\D+','', effectdata$Level2)
    
    contrastSum <- merge(effectdata, pval[,c("PO","Diet","Level1","Level2","value","signif")], 
                         by=c("PO","Diet","Level1","Level2"))
    contrastSum$Diet <- factor(contrastSum$Diet, levels=unique(effectdata$Diet))
    contrastSum$Level1 <- factor(contrastSum$Level1, levels=unique(pval$Level1))
    contrastSum$Level2 <- factor(contrastSum$Level2, levels=unique(pval$Level2))
    contrastSum <- contrastSum[order(contrastSum$PO, contrastSum$Diet, contrastSum$Level1, contrastSum$Level2),]
  }
  return(contrastSum)
}



contrast.null <- function(dfObject, levels, sig.val=0.05, sig=T)
{
  ### set up contrasts ###
  useDF <- dfObject[,as.character(levels)]
  numSamps <- nrow(useDF)
  
  ### generate things to return ###
  signif <- c()
  pval <- c()
  dir <- c()
  for (i in 1:ncol(useDF)){
    over0 <- length(which(useDF[,i]>0))
    pval[i] <- min(over0, numSamps - over0) / numSamps
    dir[i] <- ifelse(over0 > numSamps - over0, "+","-")
    if (pval[i] < sig.val/10) {
      signif[i] <- "**"
    } else if (pval[i] < sig.val){
      signif[i] <- "*"
    } else {signif[i] <- "-"}
  }
  
  #names(signif) <- colnames(useDF)
  #names(pval) <- colnames(useDF)
  contrasts <- data.frame(signif = signif, pval = pval, dir = dir)
  rownames(contrasts) <- colnames(useDF)
  
  return(contrasts)
}

contrast.test <- function(data, encoded, variable, byVar, contrast, sig.val=0.05, sig=T)
{
  ### set up contrasts ###
  comparisons <- as.character(encoded$Level[which(encoded$Variable == variable)])
  byLevel <- encoded$Level[which(encoded$Variable == byVar)]
  tempNames <- encoded$Level[grep(variable, encoded$Variable)]
  tempLevels <- unique(encoded$Variable[grep(variable, encoded$Variable)])
  useNames <- matrix(NA, nrow=length(comparisons), 
                     ncol=(length(tempNames)/length(comparisons)) + 1)
  count=1
  for(i in 1:length(tempLevels)){
    if (tempLevels[i] == variable){
      useNames[,count] <- comparisons
      count=count+1
    } else {
      tempCategory <- encoded$Level[which(encoded$Variable == tempLevels[i])]
      if(length(grep(byLevel[1],tempCategory))>0){
        for(j in 1:length(byLevel)){
          if(variable == "Diet"){
            pasteThis = paste0(byLevel[j],"$")
          } else {pasteThis = byLevel[j]}
          string <- as.character(tempCategory[grep(pasteThis,tempCategory)])
          useNames[,count] <- rep(c(string, NA), length.out=nrow(useNames))
          
          count=count+1
        }
      } else { useNames[,count] <- as.character(tempCategory)
      count = count+1
      }
    }
  }
  ### generate things to return ###
  
  signif <- matrix(NA, nrow=ncol(useNames), ncol=nrow(contrast))
  pval <- matrix(NA, nrow=ncol(useNames), ncol=nrow(contrast))
  dir <- matrix(NA, nrow=ncol(useNames), ncol=nrow(contrast))
  title <- c()
  contrastHPD <- list()
  for (k in 1:ncol(useNames)){
    testname <- useNames[,k]
    #pval[[k]] <- vector()
    for (i in 1:nrow(contrast)){
      compareTemp <- comparisons[contrast[i,]]
      if(variable == "RIX"){
        compareThis = paste0(compareTemp,"$")
      } else {compareThis = compareTemp
      }
      compare <- unique(testname[c(grep(compareThis[1], testname), 
                                   grep(compareThis[2], testname))])
      lab <- paste0(compare[1], ".", compare[2])
      numSamps <- nrow(data)
      if (k==1){
        title[i] <- lab }
      if (length(compare)==2){
        for (j in 1:numSamps){
          contrastHPD[[lab]][j] <- data[,compare[1]][j] - data[,compare[2]][j]
        }
        
        if (sig==T){
          over0 <- length(which(contrastHPD[[lab]]>0))
          pval[k,i] <- min(over0, numSamps - over0) / numSamps
          dir[k,i] <- ifelse(over0 > numSamps - over0, "+","-")
          if (pval[k,i] < sig.val/10) {
            signif[k,i] <- "**"
          } else if (pval[k,i] < sig.val){
            signif[k,i] <- "*"
          } else {signif[k,i] <- "-"}
        }
      }
    }
  }
  
  colnames(signif) <- title
  rownames(signif) <- gsub(comparisons[1],"",useNames[1,])
  colnames(pval) <- title
  rownames(pval) <- gsub(comparisons[1],"",useNames[1,])
  
  if (sig){
    return(list(signif=signif, pval=pval,contrastHPD=contrastHPD, dir=dir))
  } else {
    return(list(contrastHPD=contrastHPD, dir=dir))
  }
}


###################
# Contrast matrix #
###################

contrast_matrix <- function(encoded, levels=c("RIX","Diet")){
  nlevels=length(levels)
  output <- list()
  for(j in 1:nlevels){
    nobs=length(encoded$Level[which(encoded$Variable == levels[j])])
    tempmat <- matrix(NA, nrow=choose(nobs,2), ncol=2)
    count=1
    ind=1
    for(i in 1:nobs){
      ind = i+1
      while(ind <= nobs && count <= nrow(tempmat)){
        tempmat[count,] <- c(i, ind) 
        count=count+1
        ind=ind+1
      }
    }
    output[[paste(levels[j], "contr", sep=".")]] <- tempmat
  }
  return(output)
}
