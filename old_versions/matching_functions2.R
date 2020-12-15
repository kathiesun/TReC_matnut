source("./matching_ks.R")
source("./lm/formulaWrapper.R")
source("./matnut/summary_functions.R")
source("./matnut/prediction_functions.R")

############
# Matching #
############

match.multimp <- function(data, matchon, matchoff, idcol="ID", N=1, MatchS=F, 
                          encoded=NA, allparam=NA, ptypes=NA)
{
  multimp <- list()
  plots <- list()
  matched_mcmc <- list()
  
  if(class(data) == "list"){
    df <- data$df
    allparam <- data$parameters
    encoded <- data$encoded
    if(is.na(ptypes)){
      ptypes <- data$ptypes
    } 
  } else {
    df <- data
  }
  
  for(i in 1:length(ptypes)){
    df <- df[-which(is.na(df[,ptypes[i]])), ]
    multimp[[ptypes[i]]] <- list()
    tempimp <- list()
    matched_df <- list()
    for(j in 1:N){
      multimp[[ptypes[i]]]$keepSum <- list()
      match_init <- generateMatching(df=df, matchon=matchon, matchoff=matchoff, idcol=idcol)
      covariates <- allparam[-grep(paste("tau","sig",matchoff, sep="|"), allparam)]
      matched <- match_init$out
      for(k in 1:length(covariates)){
        cname <- covariates[k]
        for(m in 1:2){
          matched[[paste0(cname,".",m)]] = df[[cname]][match(matched[,m], df[[idcol]])]
        }
        if(all.equal(matched[[paste0(cname,".",1)]], matched[[paste0(cname,".",2)]])==T){
          matched <- matched[,-which(colnames(matched) == paste0(cname,".",2))]
        } else {
          matched <- matched[,-grep(cname, colnames(matched))]
        }
      }
      colnames(matched)[3:length(colnames(matched))] <- unlist(strsplit(colnames(matched)[3:length(colnames(matched))],"[.]"))[c(T, F)]
      matched <- matched[,-which(duplicated(colnames(matched)))]
      temp_match <- getDeltaForPhen(original=df, matched, idcol="ID", phen=ptypes[i])
      temp_match[,ptypes[i]] <- temp_match$y
      matched_df[[j]] <- temp_match
      temp_jags <- makeSummary(datalist=temp_match, encoded=encoded, phenotype = ptypes[i], MatchS=MatchS, Match=T, 
                               tryLam=c(-1, 0, .25, .33, .5, 1, 2, 3),
                               chains=2, n.adapt=2500, n.iter=10000, thin=10)
      RIX_names <- gsub("PO", "", encoded$Level[encoded$Variable == "PORIX"])
      tempDF <- as.data.frame(temp_jags$mcmcObject)
      keepNames <- colnames(tempDF)[which(colnames(tempDF) %in% RIX_names)]
      tempimp[[j]] <- as.mcmc(tempDF[,keepNames])
    }
    RIXef_tot <- do.call(rbind, tempimp)
    matched_tot <- do.call(rbind, matched_df)
    RIX_means <- data.frame("RIX_means" = colMeans(RIXef_tot))
    RIX_means[,"RIX"] <- RIX_names
    for(j in 1:N){
      for(k in 1:length(RIX_names)){   
        matched_df[[j]][which(matched_df[[j]]$RIX == RIX_names[k]), "RIX_mean"] <- RIX_means[which(RIX_means$RIX == RIX_names[k]), "RIX_means"]
      }
      matched_df[[j]][,"s"] <- ifelse(matched_df[[j]]$RIX_mean < 0, -1, 1)
      matched_df[[j]][,ptypes[i]] <- matched_df[[j]]$y * matched_df[[j]][,"s"]
      meanPhen <- mean(matched_df[[j]]$y * matched_df[[j]][,"s"])
      #add intercept for second run??
      
      formUse <- paste("~ 1", gsub("~ 1", "grand_mu", temp_jags$lmer_obj$formula), sep="+")
      matched_jags <- runJagsModels_cov(datalist=matched_df[[j]], chains=2, formula=formUse, testLMER=NULL,
                                    n.adapt=2500, n.iter=10000, thin=10, mu=(2*meanPhen),
                                    encoded=encoded, phenotype=ptypes[i])
      #matched_jags <- makeSummary(datalist=matched_df[[j]], encoded=encoded, phenotype = ptypes[i], MatchS=MatchS, Match=T, 
      #                         tryLam=c(-1, 0, .25, .33, .5, 1, 2, 3),
      #                         chains=2, n.adapt=2500, n.iter=10000, thin=10)
      matched_mcmc[[j]] <- as.data.frame(mcmc.stack(matched_jags[[ptypes[i]]]))
    }
    ## assigning PORIX mean effect to each pair with that RIX
    boundOb <- do.call(rbind, matched_mcmc)
    multimp[[ptypes[i]]]$matched_df <- matched_df
    multimp[[ptypes[i]]]$imp  <- boundOb
    multimp[[ptypes[i]]]$grand_mu <- ifelse(length(boundOb[,"grand_mu"]) > 0, mean(boundOb[,"grand_mu"]), NA)
    multimp[[ptypes[i]]]$pred <- as.mcmc(prediction(mcmcOb=boundOb, ptypes=ptypes[i], encoded, Match=T))
    multimp[[ptypes[i]]]$ribPlot <- ribbon_plot(mcmcOb=multimp[[ptypes[i]]]$pred, N=N, ptypes=ptypes[i], encoded, Match=T)
  }
  if (i==1){
    multimp <- multimp[[ptypes[i]]]
  }
  return(multimp=multimp)
}

