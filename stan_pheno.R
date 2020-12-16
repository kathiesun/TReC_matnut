library(tidyverse)
library(rstan)

setwd("C:/Users/Kathie/TReC_matnut/src")

source("lmer_functions_rna.R")
source("stan_pheno_functions.R")
source("prediction_functions.R")
source("summary_functions.R")


###### Read in data
#dir <- "C:/Users/Kathie/Dropbox\ (ValdarLab)"
dir <- "/nas/depts/006/valdar-lab/users/sunk/"
matnut = read.csv(file.path(dir,'matnut_main/AllMice_GeneExpression_SSupdated_11.27.19.csv'))
#matnut <- readRDS(file.path(dir,'phenotype_analysis/matnut_data.rds'))
genes = read.csv("deseq2/priorityTryGenes_16dec2020.csv")



matnut = matnut %>% select(-contains("X"))
gene_count <- read.csv(file.path(dir,'/hisat2_stringtie_fasta/stringtie/round2/gene_count_matrix.csv'))
#gene_count <- read.csv(file.path(dir,'/trec/gene_count_matrix.csv'))

colnames(gene_count)[1] = "gene_id"
gene_count$Gene.ID = do.call("rbind",(strsplit(as.character(gene_count$gene_id), "[|]")))[,1]
gene_count$Gene.Name = do.call("rbind", (strsplit(as.character(gene_count$gene_id), "[|]")))[,2]

samples = colnames(gene_count)[grep("Pup", colnames(gene_count))]
matnut$ID = paste0("Pup.ID_",matnut$Pup.ID)
matnut$Diet = gsub(" $", "", matnut$Diet)
matnut$RIX = gsub("a|b","",matnut$Reciprocal)
matnut$PO = ifelse(gsub("[0-9]","",matnut$Reciprocal) == "a", 0.5, -0.5)
matnut$DietRIX = gsub(" ", ".",paste0(matnut$Diet, matnut$RIX))
matnut$DietRIXPOq = paste0(matnut$DietRIX,"_",matnut$PO)
matnut$RIX = factor(matnut$RIX, levels = c(1:4,6:10))
map = read.csv(file.path(dir,"mini/combined_map_cs_gbrs_allmarkers_17mar2020.csv"))

annot = read.table(file.path(dir, "mm10/ref_sequence/Mus_musculus.GRCm38.96.gtf"), skip=5, fill=T)
#annot = read.table(file.path(dir, "matnut_main/Mus_musculus.GRCm38.96.gtf"), skip=5, fill=T)

annot_genes = annot %>% filter(V15 == "gene_name") %>%
  dplyr::select(one_of(paste0("V",c(1,4,5,7,10,16)))) %>%
  distinct()

colnames(annot_genes) = c("Chr","Start","End","Strand","Gene.ID","Gene.Name")
annot_genes = annot_genes %>% mutate(Start = as.numeric(paste(Start)), End = as.numeric(paste(End)))
# 
if(length(which(duplicated(annot_genes$Gene.Name))) > 0){
  annot_genes = annot_genes[-which(duplicated(annot_genes$Gene.Name)),]
}

annot_genes = annot_genes[which(annot_genes$Gene.ID %in% gene_count$Gene.ID |annot_genes$Gene.Name %in% gene_count$Gene.Name),]


regions_list = readRDS(file.path(dir,"mini/seg_regions_by_genotyping_perRIX_kmerBased_23nov2019.rds"))
haplofiles = list.files(file.path(dir,"mini/pup_haplo_blocks_by_CC_parent_dec2019"), 
                        pattern="haploBlocks.rds", full.names = T)
phased_par_CC_haplotypes = lapply(haplofiles, readRDS)
tmp = do.call("rbind", strsplit(haplofiles, "_"))
names(phased_par_CC_haplotypes) = tmp[,ncol(tmp)-1]

maxn = floor(nrow(annot_genes)/500)
its_1 = ((it-1)*500)+1
its_2 = ifelse(it==maxn, nrow(annot_genes), it*500)

print(paste(its_1,its_2))
stanlist <- lapply(matnut$ptypes[7], function(x) stanSum(df=matnut$df, encoded=matnut$encoded, phenotype=x, 
                   randvar=c("DamID", "RIX", "DietRIX"), fixvar="Diet", POvar=c("RIX", "Diet:RIX"),
                   tryLam=c(-1, 0, .25, .33, .5, 1, 2, 3), normd=T, 
                   chains=1, iter=2000))


  stanmcmc<- As.mcmc.list(resStan)
  summcmc <- summary(stanmcmc)
  
  #armform <- y ~ Diet + PO*RIX + Diet*RIX + PO*Diet*RIX
  #test_stanarm = stan_glm(armform, data=df, prior = lasso(autoscale = T))
  
  par(mfrow=c(3,3))
  coda::traceplot(stanmcmc[[1]][,grep("SPO[[]", colnames(stanmcmc[[1]])),drop=F])
  coda::traceplot(stanmcmc[[1]][,grep("sigma|lambda", colnames(stanmcmc[[1]])),drop=F])
  tset <- sapply(1:(ncol(stanmcmc[[1]])-1), function(x) HPDinterval(stanmcmc[,x,drop=T]), simplify=F)

  
  
  
  
  
  
reg.jags <- jags.model(textConnection(madeModel$modelFin), data=madeModel$data, n.chains = chains, n.adapt = n.adapt)
update(reg.jags, n.iter=n.adapt)
fitJags[[phenotype[i]]]$fit <- coda.samples(reg.jags, variable.names = madeModel$paramUse, thin=thin, n.iter=n.iter)
fitJags[[phenotype[i]]]$madeModel <- madeModel
allnames <- c()
paramUse <- gsub("beta","", unique(unlist(
  lapply(strsplit(varnames(fitJags[[phenotype[i]]]$fit),'[[]'), function(x) x[[1]]))))
wantParam <- intersect(paramUse, c(paste(unique(encoded$Variable)), "grand_mu"))
  
  
  
  
  
  
  jagsFit <- runJagsModels_cov(df=df, chains=chains, testLMER = lmerFit, 
                               n.adapt=n.adapt, n.iter=n.iter, thin=thin, 
                               encoded=encoded, phenotype=phenotype, sq=sq, addS=addS)
  

  
  keep <- c("Diet", "DietRIX","PODietRIX","PORIX", "RIX")
  lmerSum <- lmer.getDecodedSummary(lmerFit, lmerFit$LMformula, phenotypes=phenotype)
  test_cm <- contrast_matrix(encoded)
  
  jagsSum <- list()
  allSummary <- list()
  null.test <- list()
  jagsLmer_compare <- list()
  plotCat <- list()
  madeMod <- list()
  all.mcmc <- list()
  
  contrastMat <- contrast_matrix(encoded, levels=c("RIX","Diet"))
  
  diet.test <- list()
  rix.test<- list()
  
  for(i in 1:length(lmerFit$lmerobj$y.transform)){
    all.mcmc[[phenotype[i]]] <- mcmc.stack(jagsFit[[i]]$fit)
    jagsSum[[i]] <- jags.getDecodedSummary(all.mcmc[[phenotype[i]]], encoded)
    jagsSum[[i]]$Level <- factor(jagsSum[[i]]$Level, levels=jagsSum[[i]]$Level)
    if(!any(is.na(lmerSum[[i]]))) {
      JLSum <- compareSummary(lmerSum[[i]], jagsSum[[i]], keep)
      jagsLmer_compare[[phenotype[i]]] <- data.frame(phenotype=rep(phenotype[i],nrow(JLSum)), JLSum) 
      keepLev <- unique(JLSum$Level)
    } else {
      jagsLmer_compare[[phenotype[i]]] <- NA
      plot=F
      keepLev <- jagsSum[[i]]$Level[which(jagsSum[[i]]$Variable %in% keep)]
    }
    
    null.test[[phenotype[i]]]  <- contrast.null(dfObject=all.mcmc[[phenotype[i]]], levels=keepLev)
    
    null_p <- data.frame(Level = rownames(null.test[[phenotype[i]]]),
                         pval = null.test[[phenotype[i]]]$pval, 
                         signif = null.test[[phenotype[i]]]$signif)
    null_p$Level <- factor(null_p$Level, levels=jagsSum[[i]]$Level)
    null_p <- null_p[order(null_p$Level),]
    mergeTable <- merge(jagsSum[[i]], null_p, by="Level")
    mergeTable$Level <- factor(mergeTable$Level, levels=encoded$Level[which(encoded$Variable %in% keep)])
    mergeTable$Variable <- factor(mergeTable$Variable, levels=unique(encoded$Variable)[which(encoded$Variable %in% keep)])
    mergeTable <- mergeTable[order(mergeTable$Level),]
    allSummary[[phenotype[i]]]  <- cbind(phenotype=rep(phenotype[i],nrow(mergeTable)), mergeTable[order(mergeTable$Variable, mergeTable$Level),])
    
    if(plot) plotCat[[phenotype[i]]] <- plotSummary(lmerSum[[i]], jagsSum[[i]], JLSum, phenotype)
    madeMod[[phenotype[i]]]  <- jagsFit[[i]]$madeModel
    
    if(contrasts){
      diet.test[[phenotype[i]]] <- contrast.test(data=all.mcmc[[phenotype[i]]], encoded, variable="Diet", byVar="RIX", contrast=contrastMat$Diet.contr)
      rix.test[[phenotype[i]]] <- contrast.test(data=all.mcmc[[phenotype[i]]], encoded, variable="RIX", byVar="Diet", contrast=contrastMat$RIX.contr)
    }
  }
  
  return(list(transformation=transf, mcmcObject=all.mcmc, plot=plotCat, allSummary=allSummary,
              null.test=null.test, diet.test=diet.test, rix.test=rix.test, S=S, modelDef=madeMod, 
              JL_compare=jagsLmer_compare, lmer_obj = lmerFit))
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

