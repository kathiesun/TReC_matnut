source("./matching_ks.R")
source("./lm/formulaWrapper.R")

src <- ifelse(length(grep("C:", getwd())) >0, file.path("C:","Users","kathie","Documents","source_ks"), 
              file.path("~","Documents","matnut_src"))
source(file.path(src, "summary_functions.R"))
source(file.path(src, "prediction_functions.R"))

############
# Matching #
############

match.multimp_direct <- function(data, matchon, matchoff, idcol="ID", N=1, MatchS=F, 
                          encoded=NA, allparam=NA, ptypes=NA)
{
  multimp <- list()
  #comb_multimp <- list()
  plots <- list()
  
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
  #browser()
  for(i in 1:length(ptypes)){
    if(length(which(is.na(df[,ptypes[i]]))) > 0){
      df <- df[-which(is.na(df[,ptypes[i]])), ]
    }
    multimp[[ptypes[i]]] <- list()
    #comb_multimp[[ptypes[i]]] <- list()
    tempimp <- list()
    for(j in 1:N){
      multimp[[ptypes[i]]]$keepSum <- list()
      match_init <- generateMatching(df=df, matchon=matchon, matchoff=matchoff, idcol=idcol)
      #covariates <- c("DamID", "SireID","BreedingBatch","BehaviorBatch","Diet","RIX",
      #                "FSTChamber", "LDChamber","RestraintExperimenter","RestraintOrder", 
      #                "SIHOrder", "OFBox","DietRIX")
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

      matched_df <- getDeltaForPhen(original=df, matched, idcol="ID", phen=ptypes[i])
      #if (MatchS){
        rix_ave <- aggregate(y ~ RIX, data=matched_df, FUN = function(x) M=mean(x))
        rix_ave[,"s"] <- as.numeric(ifelse(rix_ave$y<0, -1, 1))
        matched_df[,"s"] <- rix_ave$s[match(matched_df$RIX, rix_ave$RIX)] 
        matched_df[,ptypes[i]] <- matched_df$y * matched_df$s
        #keepNames <- c("s[1]", "s[2]")
        #formTemp <- formula(paste(ptypes[i], "~RIX"))
        #cross_ave <- aggregate(formTemp, data=matched_df, FUN = function(x) M=mean(x))
        #invisible(plot(density(rix_ave$y)))
        #plots$diff <- recordPlot()
        #invisible(plot(density(cross_ave[[ptypes[i]]])))
        #plots$std <- recordPlot()
      #} else {
      #  matched_df[,ptypes[i]] <- matched_df$y
      #}
      keepNames <- c()
      browser()
      temp_jags <- makeSummary(datalist=matched_df, encoded=encoded, phenotype = ptypes[i], MatchS=MatchS, Match=T, 
                               tryLam=c(-1, 0, .25, .33, .5, 1, 2, 3),
                               chains=2, n.adapt=2500, n.iter=10000, thin=10)
      tempDF <- as.data.frame(temp_jags$mcmcObject)
      keepNames <- c(colnames(tempDF)[which(colnames(tempDF) %in% encoded$Level)], keepNames)
      tempimp[[j]] <- as.mcmc(tempDF[,keepNames])
      #multimp[[ptypes[i]]]$keepSum[[j]] <- temp_jags$plot
      #multimp[[ptypes[i]]]$combined <- Reduce("+", multimp[[ptypes[i]]]) / length(multimp[[ptypes[i]]])
    }
    #comb_multimp[[ptypes[i]]] <- as.mcmc(multimp[[ptypes[i]]][,which(colnames(multimp[[ptypes[i]]]) %in% encoded$Level)])
    boundOb <- do.call(rbind, tempimp)
    boundObSum <- jags.getDecodedSummary(as.mcmc(boundOb), encoded)
    multimp[[ptypes[i]]]$catPlot <- plot.inter.ci(med=boundObSum$med, mu=boundObSum$mu, 
                                    hpd.narrow = cbind(boundObSum$lower, boundObSum$upper), 
                                    hpd.wide = cbind(boundObSum$lower.1, boundObSum$upper.1), 
                                    names=boundObSum$Level, order=3,
                                    col.midvals="white", pch.midvals="|", addline = F, wide=T, 
                                    grouped=boundObSum$Variable, ordered=F)
    multimp[[ptypes[i]]]$imp  <- boundOb
    multimp[[ptypes[i]]]$pred <- as.mcmc(prediction(mcmcOb=boundOb, ptypes=ptypes[i], encoded, Match=T))
    multimp[[ptypes[i]]]$ribPlot <- ribbon_plot(mcmcOb=multimp[[ptypes[i]]]$pred, N=N, ptypes=ptypes[i], encoded, Match=T)
  }
  if (i==1){
    multimp <- multimp[[ptypes[i]]]
  }
  return(multimp=multimp)
}

