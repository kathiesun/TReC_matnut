library(data.table)
source("plot.hpd_ks.R")
source("prediction_functions.R")

#######################
# Generate all tables #
#######################

writeTables <- function(results, df, ptypes=NULL, encoded=NULL, print=F, output="."){
  dietDiet_table <- list()
  rixRix_table <- list()
  
  jagsLmer_compare <- list()
  summary_table <- list()
  
  if (class(df) == "list"){
    if(is.null(ptypes)) ptypes <- df$ptypes
    if(is.null(encoded)) encoded <- df$encoded
  }
  
  for (i in 1:length(ptypes)){
    ### summary ###
    tempTable <- jags.getDecodedSummary(results$mcmcObject[[i]], encoded, narrow=0.5, wide=0.95)
    tempTable$Level <- factor(tempTable$Level, levels=tempTable$Level)
    null_p <- data.frame(Level = rownames(results$null.test[[i]]), 
                         pval = results$null.test[[i]]$pval, 
                         signif = results$null.test[[i]]$signif)
    null_p$Level <- factor(null_p$Level, levels=tempTable$Level)
    null_p <- null_p[order(null_p$Level),]
    mergeTable <- merge(tempTable, null_p, by="Level")
    mergeTable$Level <- factor(mergeTable$Level, levels=encoded$Level)
    mergeTable$Variable <- factor(mergeTable$Variable, levels=unique(encoded$Variable))
    shortTable <- mergeTable[which(mergeTable$Variable %in% c("Diet","RIX","DietRIX","PORIX","PODietRIX")),]
    
    summary_table[[ptypes[i]]] <- cbind(phenotype=rep(ptypes[i],nrow(shortTable)), 
                                        shortTable[order(shortTable$Variable, shortTable$Level),])
    ### contrasts ###
    if(!is.null(results$diet.test)){
    dietTab <- contrasts.getDecodedSummary(results$diet.test[[ptypes[i]]], compare = "Diet")
      dietDiet_table[[ptypes[i]]] <- cbind(phenotype=rep(matnut$ptypes[i],nrow(dietTab)), dietTab)
      rixTab <- contrasts.getDecodedSummary(results$rix.test[[i]], compare = "RIX")
      rixRix_table[[ptypes[i]]] <- cbind(phenotype=rep(ptypes[i],nrow(rixTab)), rixTab)
    }
  }
  summary_PO <- NULL
  
  if(!is.null(results$diet.test)){
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
  }
  
  diet_contr <- rbindlist(dietDiet_table)
  colnames(diet_contr) <- c("Phenotype","PO","RIX", "Level1", "Level2", "Mean","Median", 
                            "Lower_50", "Upper_50", "Lower_95","Upper_95","Value","Signif")
  rix_contr=rbindlist(rixRix_table)
  colnames(rix_contr) <- c("Phenotype","PO","Diet", "Level1", "Level2", "Mean","Median", 
                           "Lower_50", "Upper_50", "Lower_95","Upper_95","Value","Signif")
  summary_table = rbindlist(summary_table)
  colnames(summary_table) <- c("Phenotype","Level", "Mean","Median", 
                               "Lower_50", "Upper_50", "Lower_95","Upper_95","Variable","Pval","Signif")
  
  if(print){
    if(!is.null(summary_PO)){
      write.csv(cbind(summary_noPO, summary_PO), file.path(paste0(output,'/signif_counts.csv')), row.names = T)
      write.csv(diet_contr, file.path(paste0(output,'/allPheno_dietDiet_table.csv')), row.names = F)
      write.csv(rix_contr, file.path(paste0(output,'/allPheno_rixRix_table.csv')), row.names = F)
    }
    
    write.csv(summary_table, file.path(paste0(output,'/allPheno_summarytable.csv')), row.names = F)
  }
  return(list(diet_contr=diet_contr, rix_contr=rix_contr, 
              signif_counts=cbind(summary_noPO, summary_PO), summary_table = summary_table)) 
}


###########
#  Plots  #
###########


### Flag Plots ###

flagPlots <- function(results, df, ptypes=NULL, encoded=NULL, print=F, byVar = c("Diet","RIX"), output="."){
  if(length(byVar) > 1) byVar = byVar[1]
  contrastPlots <- list()
  returnPval <- list()
  if (class(df) == "list"){
    if(is.null(ptypes)) ptypes <- df$ptypes
    if(is.null(encoded))encoded <- df$encoded
  } 

  for(i in 1:length(ptypes)){  
    useTest <- ifelse(byVar=="Diet","diet.test","rix.test")
    bindTemp <- results[[useTest]][[ptypes[i]]]$dir
    rownames(bindTemp) <- rownames(results[[useTest]][[ptypes[i]]]$pval)
    colnames(bindTemp) <- paste0("a",colnames(results[[useTest]][[ptypes[i]]]$pval))
    meltTemp <- melt(bindTemp)
    colnames(meltTemp)[3] <- "Dir"
    tempPval <- results[[useTest]][[ptypes[i]]]$pval
    colnames(tempPval) <- paste0("a",colnames(results[[useTest]][[ptypes[i]]]$pval))
    
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
    contrastPlots[[i]] <- getPlot.compare.helper(dfplot,ptypes=ptypes[i],byVar=byVar)
  }
  
  if(print) {
    graphics.off()
    pdf(file.path(paste0(output,"/allPheno_",byVar,"_ContrastPlots.pdf")), 
        onefile=TRUE, width = 8.5, height = 11)
    for(i in 1:length(ptypes)){  
      print(contrastPlots[[i]])
    }
    graphics.off()
  }
  return(contrastPlots)
}

generateRibbonPlots <- function(results, df, ptypes=NULL, encoded=NULL, print=F, output="."){
  if (class(df) == "list"){
    if(is.null(ptypes)) ptypes <- df$ptypes
    if(is.null(encoded)) encoded <- df$encoded
  }
  ribbon_plots <- list()
  keep <- c("Diet", "DietRIX","PODietRIX","PORIX", "RIX")
  ribbon_plots <- ribbon_plot(mcmcList=results$mcmcObject, ptypes, encoded)
  browser()
  temp_pred <- prediction(mcmcList=results$mcmcObject, ptypes=ptypes, encoded, Match=F)
  temp_pred_M <- prediction(mcmcList=results$mcmcObject, ptypes=ptypes, encoded, Match=T)
  
  temp_rib <- ribbon_plot(mcmcList=temp_pred, ptypes=ptypes, encoded)
  all_phenDF <- data.frame()
  
  for(i in 1:length(ptypes)){
    phen <- ptypes[i]
    len <- length(temp_rib$savedata[[phen]])
    tempDF <- data.frame(phen = rep(phen, nrow(temp_rib$savedata[[phen]][[len]])), 
                         temp_rib$savedata[[phen]][[len]])
    all_phenDF <- rbind(all_phenDF, tempDF)
    
    #myMCMC <- results$mcmcObject[[ptypes[i]]][, which(colnames(results$mcmcObject[[ptypes[i]]]) %in% 
    #                                                    encoded$Level[which(encoded$Variable %in% keep)])]
    
  }
  allPlot <- ggplot(all_phenDF, aes(color=rix)) + 
    geom_line(aes(y=mu_tot, x=diet, group=poe, alpha=poe)) +  
    geom_ribbon(color=NA, aes(x=diet, ymin=lower.1, ymax=upper.1, group=poe, alpha=poe, 
                              fill=rix, color=rix)) + 
    geom_hline(yintercept = 0, colour="lightsteelblue", size=0.5, linetype="longdash") +
    scale_alpha_discrete(range = rev(c(0.45, 0.75)), name="POE") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_fill_discrete(guide=FALSE) + scale_color_discrete(guide=FALSE) + 
    facet_grid(rix~phen) + theme_minimal() + 
    ylab("Predicted effects") + xlab("Diet") #+ ggtitle(ptypes[j])
  
  if(print){
    graphics.off()
    pdf(file.path(paste0(output,"/ribbonPlots.pdf")), onefile=TRUE, width = 14, height = 11)
    for(i in 1:length(ptypes)){  
      print(ribbon_plots$ribbon[[ptypes[i]]])
    }
    
    pdf(file.path(paste0(output,"/predictionPlots.pdf")), onefile=TRUE, width = 14, height = 11)
    for(i in 1:length(ptypes)){  
      print(temp_rib$ribbon[[ptypes[i]]])
    }

    pdf(file.path(paste0(output,"/allPlot.pdf")), onefile=TRUE, width = 14, height = 11)
    print(allPlot)
    graphics.off()
  }
  return(list(ribbon_plots=ribbon_plots, predictive_plots = temp_rib, allPlot = allPlot))
}


