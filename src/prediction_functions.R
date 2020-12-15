##############
# Prediction #
##############

prediction <- function(mcmcList, ptypes, encoded, Match=F){
  pred <- list()
  for (a in 1:length(ptypes)){
    if(length(mcmcList) == length(ptypes)){
      mcmcOb <- mcmcList[[a]]
    } else if(class(mcmcList) == "data.frame") {
      mcmcOb <- mcmcList#[grep(ptypes[a], names(mcmcList))]
      #colnames() = gsub(paste0(ptypes[a],"."), "", colnames(mcmcOb))
    } else{
      mcmcOb <- mcmcList[[which(names(mcmcList) == ptypes[a])]]
    }
    ribbon <- list()
    predPOneg <- list()
    predPOpos <- list()

    its <- 1
    #for (i in 1:its){
    #  if(class(mcmcOb) == "list"){
    #    mcmcOb <- mcmcOb[[i]]
    #    its <- length(mcmcOb)
    #  } 
      encoded$Variable = factor(encoded$Variable, levels=unique(encoded$Variable))
      #encoded$Variable = factor(encoded$Variable, levels=unique(encoded$Variable))
      levelz <- unique(encoded$Variable[match(colnames(mcmcOb), encoded$Level)])
      hi_level <- ifelse(max(as.numeric(levelz), na.rm=T) < 5, levels(encoded$Variable)[max(as.numeric(levelz), na.rm=T)], 
                         levels(encoded$Variable)[5])
      hi_variables <- gsub("PO", "", encoded$Level[which(encoded$Variable == hi_level)])
      
      if (length(grep("PO", encoded$Level[which(encoded$Variable == hi_level)])) > 0){
        dimnames <- c(paste(hi_variables,"neg", sep="_"), paste(hi_variables,"pos", sep="_"))
      } else {
        dimnames <- hi_variables
      }
      
      pred[[ptypes[a]]] <- data.frame()     #row.names=dimnames
      
      for(j in 1:length(hi_variables)){
        temp_var <- hi_variables[j]
        for(k in 1:nrow(mcmcOb)){
          labs <- c(gsub("PO", "", temp_var), gsub("\\d|PO", "", temp_var), gsub("\\D", "", temp_var), "grand_mu")
          POlabs <- paste0("PO", labs)
          ind <- match(labs, colnames(mcmcOb))[!is.na(match(labs, colnames(mcmcOb)))]
          POind <- match(POlabs, colnames(mcmcOb))[!is.na(match(POlabs, colnames(mcmcOb)))]
          temp_pred <- 0
          for(m in 1:length(ind)){
            temp_pred <- temp_pred + mcmcOb[k, ind[m]]
          }
          if(length(POind) > 0){
            temppos = tempneg = temp_pred
            for(l in 1:length(POind)){
              temppos <- temppos + mcmcOb[k, POind[l]]*0.5
              tempneg <- tempneg + mcmcOb[k, POind[l]]*-0.5
            }
            pred[[ptypes[a]]][k, as.character(paste(temp_var,"neg",sep="_"))] <- temppos
            pred[[ptypes[a]]][k, as.character(paste(temp_var,"pos",sep="_"))] <- tempneg
          } else {
            pred[[ptypes[a]]][k, as.character(temp_var)] <- temp_pred
          }
        }
      }
      pred[[ptypes[a]]] <- as.mcmc(pred[[ptypes[a]]])
    #}
    print(paste("pred finished:",a))
  }
  return(pred)
}
