source("matching.R")
#source("matching/generate.R")
source("./matnut/summary_functions.R")
source("./lm/formulaWrapper.R")
source("./matnut/prediction_functions.R")
source("./matnut/tables_n_plots.R")
source("./matnut/plot.hpd_ks.R")



#src <- ifelse(length(grep("C:", getwd())) >0, file.path("C:","Users","kathie","Documents","source_ks"), 
#              ifelse(length(grep("kys6", getwd())) >0, file.path("~","source_ks/matnut_src"), file.path("~","Documents","matnut_src")))


############
# Matching #
############

match.multimp <- function(data, matchon, matchoff, idcol="ID", N=3, addS=c("off", "force", "fit", "estimate"), 
                          randvar = NA, fixvar = NA, 
                          chains=2, n.adapt=2500, n.iter=10000, thin=10,
                          encoded=NA, allparam=NA, ptypes=NA, sq=F){
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
    if(length(which(is.na(df[,ptypes[i]]))) > 0){
      df <- df[-which(is.na(df[,ptypes[i]])), ]
    }
    multimp[[ptypes[i]]] <- list()
    tempimp <- list()
    matched_df <- list()
    
    if(addS == "force"){
      temp_jags <- makeSummary(datalist=df, encoded=encoded, phenotype = ptypes[i], 
                               randvar = randvar, fixvar = fixvar,
                               tryLam=c(-1, 0, .25, .33, .5, 1, 2, 3), addS = addS, 
                               chains=chains, n.adapt=n.adapt, n.iter=n.iter, thin=thin)
      rixes <- encoded$Level[which(encoded$Variable == "RIX")]
      df[, "PO_old"] <- df[, "PO"]
      for(j in 1:length(rixes)){
        df[which(df$RIX == paste(rixes[j])), "PO"] <- df[which(df$RIX == paste(rixes[j])), "PO_old"] * temp_jags$S[j]
      }
      addS = "off"
    }
    
    for(j in 1:N){
      multimp[[ptypes[i]]]$keepSum <- list()
      match_init <- matching$generateMatching(df=df, matchon=matchon, matchoff=matchoff, idcol=idcol)
      covariates <- allparam[-grep(paste("tau","sig",matchoff, sep="|"), allparam)]
      covariates = setdiff(c(randvar, fixvar), 
                           unlist(strsplit(colnames(match_init$out), "[.]"))[c(T, F)])
      matched <- match_init$out
      if(length(covariates) > 0){
        for(k in 1:length(covariates)){
          cname <- gsub(":","",covariates[k])
          for(m in 1:2){
            matched[[paste0(cname,".",m)]] = df[[cname]][match(matched[,m], df[[idcol]])]
          }
          if(all.equal(matched[[paste0(cname,".",1)]], matched[[paste0(cname,".",2)]])==T){
            matched <- matched[,-which(colnames(matched) == paste0(cname,".",2))]
          } else if (cname != "DamID"){
            matched <- matched[,-grep(cname, colnames(matched))]
          }
        }
      }
      
      matched[[paste0("PO",".",1)]] = df[["PO"]][match(matched[,1], df[[idcol]])]
      keepInd <- c("ID.1","ID.2","DamID.1","DamID.2")         
      colnames(matched)[-which(colnames(matched) %in% keepInd)] <- unlist(strsplit(colnames(matched)[-which(colnames(matched) %in% keepInd)],"[.]"))[c(T, F)]
      
      if(length(duplicated(colnames(matched)))>0){
        matched <- matched[,-which(duplicated(colnames(matched)))]
      } 
#browser()
      temp_match <- matching$getDeltaForPhen(original=df, matched, idcol="ID", phen=ptypes[i])
      temp_match[,ptypes[i]] <- temp_match$y
      matched_df[[j]] <- temp_match

      #colnames(temp_match) <- gsub("[.]","",colnames(temp_match))
      matched_jags <- makeSummary(datalist=temp_match, encoded=encoded, phenotype = ptypes[i], Match=T,
                               tryLam=c(-1, 0, .25, .33, .5, 1, 2, 3), addS=addS,
                               randvar = randvar, fixvar = fixvar,
                               chains=chains, n.adapt=n.adapt, n.iter=n.iter, thin=thin)
      #RIX_names <- gsub("PO", "", encoded$Level[encoded$Variable == "PORIX"])
      #tempDF <- as.data.frame(matched_jags$mcmcObject)
      #keepNames <- colnames(tempDF)[which(colnames(tempDF) %in% RIX_names)]
      #tempimp[[j]] <- as.mcmc(tempDF[,keepNames])
      matched_mcmc[[j]] <- as.data.frame(matched_jags$mcmcObject)
      
      if(addS == "fit"){
        betaS <- colnames(matched_mcmc[[j]])[grep("betaS", colnames(matched_mcmc[[j]]))]
        keepNames <- colnames(matched_mcmc[[j]])[which(colnames(matched_mcmc[[j]]) %in% encoded$Level)]
        for(k in 1:length(keepNames)){
          whichRIX <- match(gsub("[A-z]", "", keepNames[k]), encoded$Level[which(encoded$Variable == "RIX")])
          matched_mcmc[[j]][,keepNames[k]] <- matched_mcmc[[j]][,keepNames[k]] * matched_mcmc[[j]][,betaS[whichRIX]]
        }
      }
    }
    #browser()
    #if(addS == "estimate"){
      #keepNames <- colnames(matched_mcmc[[1]])[which(colnames(matched_mcmc[[1]]) %in% encoded$Level)]
      #matched_mcmc[[j]][,keepNames] <- matched_mcmc$keepnames
    #  for(k in 1:length(RIX_names)){   
    #    matched_df[[j]][which(matched_df[[j]]$RIX == RIX_names[k]), "RIX_mean"] <- RIX_means[which(RIX_means$RIX == RIX_names[k]), "RIX_means"]
    #  }
    #  matched_df[[j]][,"s"] <- ifelse(matched_df[[j]]$RIX_mean < 0, -1, 1)
    #  matched_df[[j]][,ptypes[i]] <- matched_df[[j]]$y * matched_df[[j]][,"s"]
    #  meanPhen <- mean(matched_df[[j]]$y * matched_df[[j]][,"s"])
    #}
    if(addS == "estimate"){
      RIXef_tot <- do.call(rbind, tempimp)
      #matched_tot <- do.call(rbind, matched_df)
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
        
        formUse <- paste("~ 1", temp_jags$lmer_obj$formula, sep="+") #gsub("~ 1", "grand_mu"...
        matched_jags <- runJagsModels_cov(datalist=matched_df[[j]], chains=2, formula=formUse, testLMER=NULL,
                                          chains=chains, n.adapt=n.adapt, n.iter=n.iter, thin=thin, 
                                          mu=(2*meanPhen), addS="off",
                                          sq=sq, encoded=encoded, phenotype=ptypes[i])
        #matched_jags <- makeSummary(datalist=matched_df[[j]], encoded=encoded, phenotype = ptypes[i], MatchS=MatchS, Match=T, 
        #                         tryLam=c(-1, 0, .25, .33, .5, 1, 2, 3),
        #                         chains=2, n.adapt=2500, n.iter=10000, thin=10)
        
        matched_mcmc[[j]] <- as.data.frame(mcmc.stack(matched_jags[[ptypes[i]]]))
      }
    }
    
    # assigning PORIX mean effect to ecah pair with that RIX
    boundOb <- do.call(rbind, matched_mcmc)
    colnames(boundOb) = gsub(paste0(ptypes[i],"."),"",colnames(boundOb)) 
    boundOb = as.mcmc(boundOb)
    
    summary = jags.getDecodedSummary(boundOb, encoded)
    multimp[[ptypes[i]]]$summary <- summary
    multimp[[ptypes[i]]]$levels <- factor(summary$Level, levels=summary$Level)
    
    multimp[[ptypes[i]]]$matched_df <- matched_df
    multimp[[ptypes[i]]]$mcmc_ob  <- boundOb
    #multimp[[ptypes[i]]]$grand_mu <- ifelse(length(boundOb[,"grand_mu"]) > 0, mean(boundOb[,"grand_mu"]), NA)
    multimp[[ptypes[i]]]$pred <- as.mcmc(prediction(mcmcList=as.data.frame(boundOb), ptypes=ptypes[i], encoded, Match=T))
    multimp[[ptypes[i]]]$plot <- plot.inter.ci(med=summary$med, mu=summary$mu, 
                                               hpd.narrow = cbind(summary$lower, summary$upper), 
                                               hpd.wide = cbind(summary$lower.1, summary$upper.1), 
                                               names=multimp[[ptypes[i]]]$levels, order=3,
                                               col.midvals="white", pch.midvals="|", addline = F, wide=T, 
                                               grouped=summary$Variable, ordered=F)
    multimp[[ptypes[i]]]$ribPlot <- ribbon_plot(mcmcList=multimp[[ptypes[i]]]$pred, 
                                                N=N, ptypes=ptypes[i], encoded, Match=T)
      #ribbon_plot(mcmcOb=multimp[[ptypes[i]]]$pred, N=N, ptypes=ptypes[i], encoded, Match=T)
      #grand_mean=multimp[[ptypes[i]]]$grand_mu
    
  }
  #browser()
  if (i==1){
    multimp <- multimp[[ptypes[i]]]
  }
  return(multimp=multimp)
}

