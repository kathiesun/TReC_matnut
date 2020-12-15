<<<<<<< HEAD
library(data.table)
library(lmerTest)

setwd("/nas02/home/k/y/kys6/source_dan/nutridiallel2014/mnp2017/src")
source(file.path("..","..","..","..","source_ks","matnut_src",'parArgs.R'))
lmerFunc <- file.path("..","..","..","..","source_ks","matnut_src",'lmer_functions_rna.R')
source(lmerFunc) 
cov_short <- readRDS(file.path("..","..","..","..","Data","matnut_outputs",'cov_short_data.rds'))



rixes <- cov_short$encoded$Level[cov_short$encoded$Variable == "RIX"]
phenotypes <- colnames(cov_short$df)[grep("ENS", colnames(cov_short$df))]
cov_short$df[,phenotypes] <- sqrt(cov_short$df[,phenotypes])

rnaseq <- list()
tryLam=c(0, 0.5)

#for(j in 1:length(rixes)){
j=1
  datause <- cov_short$df[which(paste(cov_short$df$RIX) == rixes[j]), ]
    #for(i in 1:100){  #length(phenotypes)){  
      #lmerFit <- lm(get(phenotypes[i]) ~ PO + Diet, data=datause)
      
      lmerFit <- runLMERModels_cov(dataf=datause, tryLam=tryLam, phenotype=cov_short$ptypes,
                                   fixvar=c("PO", "Diet"), 
                                   parallelArgs = getBestParArgs(filesToSource = lmerFunc))
      #lmerSum <- lmer.getDecodedSummary(lmerOb=lmerFit, 
      #                                 lmerFit$formula, phenotypes=names(lmerFit$lmerobj))
      
      #anovaFit <- anova(lmerFit)
      #summaryFit <- summary(lmerFit)

      rnaseq <- list(summary = lmerSum, fit = lmerFit)
      #print(paste("RIX",j,"gene",i))
    #}
  saveRDS(lmerFit, paste0("rnaseq.par",j,".RDS"))
	
  rm(rnaseq)
  gc()     
#}



=======
library(data.table)
library(lmerTest)

setwd("/nas02/home/k/y/kys6/source_dan/nutridiallel2014/mnp2017/src")
source(file.path("..","..","..","..","source_ks","matnut_src",'parArgs.R'))
lmerFunc <- file.path("..","..","..","..","source_ks","matnut_src",'lmer_functions_rna.R')
source(lmerFunc) 
cov_short <- readRDS(file.path("..","..","..","..","Data","matnut_outputs",'cov_short_data.rds'))



rixes <- cov_short$encoded$Level[cov_short$encoded$Variable == "RIX"]
phenotypes <- colnames(cov_short$df)[grep("ENS", colnames(cov_short$df))]
#cov_short$df[,phenotypes] <- sqrt(cov_short$df[,phenotypes])

rnaseq <- list()
tryLam=c(0, 0.5)

library(foreach)
library(snow)
library(doSNOW)



cl <- makeCluster(4, type="SOCK")
registerDoSNOW(cl)
            
output.lines <- foreach(i = (1:128)) %dopar% {
  hello.world(i)
}
            
            
output.lines <- foreach(i = 1:length(rixes)) %dopar% {
  datause <- cov_short$df[which(paste(cov_short$df$RIX) == rixes[i]), ]
  #lmFit <- lm(as.matrix(datause[,grep("ENS",colnames(datause))]) ~ PO + Diet, data=datause)
  lmerFit <- runLMERModels_cov(dataf=datause[,1:100], tryLam=tryLam, phenotype=cov_short$ptypes,
                               fixvar=c("PO", "Diet"))
}

datause <- cov_short$df[which(paste(cov_short$df$RIX) == rixes[j]), ]
phenOnly <- datause[,grep("ENS", colnames(datause))]
lmFit <- lm(as.matrix(datause[,grep("ENS",colnames(datause))]) ~ PO + Diet, data=datause)
PO <- datause$PO
Diet <- datause$Diet
#for(j in 1:){
j=1
  datause <- cov_short$df[which(paste(cov_short$df$RIX) == rixes[j]), ]
    #for(i in 1:100){  #length(phenotypes)){  
      #lmerFit <- lm(get(phenotypes[i]) ~ PO + Diet, data=datause)
      
      lmerFit <- runLMERModels_cov(dataf=datause, tryLam=tryLam, phenotype=cov_short$ptypes,
                                   fixvar=c("PO", "Diet"))
      #lmerSum <- lmer.getDecodedSummary(lmerOb=lmerFit, 
      #                                 lmerFit$formula, phenotypes=names(lmerFit$lmerobj))
      
      #anovaFit <- anova(lmerFit)
      #summaryFit <- summary(lmerFit)

      #rnaseq <- list(summary = lmerSum, fit = lmerFit)
      #print(paste("RIX",j,"gene",i))
    #}
  saveRDS(lmerFit, paste0("rnaseq.par",j,".RDS"))
	
  rm(rnaseq)
  gc()     
#}


>>>>>>> eb5ca7c96557ca30eb215d35a723a8d65ed94959
