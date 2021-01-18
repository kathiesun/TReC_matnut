source("boxcox_functions.R")
library(rstan)

stanSum <- function(df, phenotype, encoded=NULL, 
                    tryLam=1, normd=T,
                    nchains=2, iter=10000,
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
                      nchains = nchains, iter = iter, warmup = warmup, thin = thin))
  
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



stan_ase <- function(df, encoded=NULL, 
                     mu_g0=0.5, t20=0.001, 
                     terms=NULL, STZ=T, use_gene=F, 
                     nchains=2, iter=10000, 
                     plot=T, contrasts=F,
                     C_diet_4=NULL, C_diet_2=NULL,C_diet_3=NULL,
                     stanMod = NULL){
  terms = toupper(terms)
  colnames(df) <- toupper(colnames(df))
  
  if(is.null(encoded)){
    encoded <- unique(getEncoding(df, terms = unique(c(terms,"PUP.ID","SEQ.GENE","PUP_GENE"))))
  }
  
  df$PUP.ID <- factor(df$PUP.ID, encoded$Level[which(encoded$Variable == "PUP.ID")])
  #data$RIX_DIR_DIET = paste(data$RRIX, data$DIR, data$DIET, sep="_")
  df %>% ungroup() %>% arrange(PUP.ID) -> df
  
  df %>% dplyr::select(one_of("PUP.ID")) %>%    
    left_join(data.frame(PUP.ID = unique(df$PUP.ID), 
                         INDP = seq(1:length(unique(df$PUP.ID)))), by="PUP.ID") %>%
    distinct() -> merge_indP
  
  df %>% 
    dplyr::select(one_of("SEQ.GENE","PUP.ID")) %>%
    distinct() %>%
    ungroup() %>%
    mutate(INDG=seq(n())) -> merge_indG
  
  merge_indG$indP = as.numeric(merge_indG$PUP.ID)
  
  df %>% group_by(PUP.ID, SEQ.GENE) %>%
    mutate(indK=seq(n())) %>%
    group_by(PUP.ID) %>%
    left_join(merge_indP, by = c("PUP.ID")) %>% 
    left_join(merge_indG, by = c("SEQ.GENE", "PUP.ID")) %>%
    dplyr::select(-one_of(toupper(c("seq.ProbeSeq", "seq.TranscriptID", "seq.end5", "seq.end3", "seq.consensus", 
                                    "seq.Chromosome", "seq.Position", "seq.refseq", "seq.altseq", "seq.alt", "seq.ref")))) %>%
    arrange(PUP.ID, SEQ.GENE) %>%
    ungroup() -> df
  
  df %>% group_by(SEQ.GENE, PUP.ID) %>%
    summarize(sum_CC_1 = sum(CC_1), sum_CC_2 = sum(CC_2), sumTot = sum(SUM)) %>%
    arrange(PUP.ID, SEQ.GENE) %>%
    left_join(merge_indP, by="PUP.ID") %>%
    left_join(merge_indG, by = c("SEQ.GENE", "PUP.ID")) -> gene_tot

  if(use_gene){
    y_gk = gene_tot$sum_CC_1
    N_gk = gene_tot$sumTot
    indG = gene_tot$INDG
    indP = gene_tot$INDP
  } else {
    y_gk <- df$CC_1
    N_gk <- df$SUM
    indG <- df$INDG
    indP =  df$INDP
  }
  
  nGP = length(unique(df$INDG))
  nP <- length(unique(df$PUP.ID))
  nG <- length(unique(df$SEQ.GENE))
  df$PO <- ifelse(df$DIR == "a", 0.5, -0.5)
  
  mu_0 <- mean(y_gk/N_gk)
  tau_0 <- 1/var(y_gk/N_gk)
  
  indP    = model.matrix(~ 0 + PUP.ID, df)     # + DietRIX
  indG    = model.matrix(~ 0 + SEQ.GENE, df)     # + DietRIX
  #indRIX  = model.matrix(~ 0 + RIX, df)     # + DietRIX
  indGP   = model.matrix(~ 0 + PUP_GENE, df)     # + DietRIX
  
  
  indDiet = model.matrix(~ 0 + DIET, df)     # + DietRIX
  nDiet   = ncol(indDiet)
  C_diet  = as.matrix(get(paste0("C_diet_",nDiet)))
  XC      = indDiet %*% C_diet
  
  indDR   = model.matrix(~ 0 + DIETRIX, df)     # + DietRIX
  nDR     = ncol(indDR)
  indDR   = indDR[,-ncol(indDR)]

  map_gp_p = model.matrix(~ 0 + PUP.ID, unique(df %>% select(SEQ.GENE, PUP.ID, PUP_GENE)))
  map_g_gp = t(model.matrix(~ 0 + SEQ.GENE, unique(df %>% select(SEQ.GENE, PUP.ID, PUP_GENE))))
  map_g_p = map_g_gp %*% map_gp_p
  
  standat <-  list(N           = length(y_gk),
                   y_gk        = y_gk, 
                   N_gk        = N_gk,
                   nP          = nP,
                   indP        = indP,
                   nG          = nG,
                   indG        = indG,
                   nGP         = nGP,
                   indGP       = indGP,
                   map_gp_p    = map_gp_p,
                   map_g_gp    = map_g_gp,
                   map_g_p     = map_g_p,
                   weight_g    = 1/rowSums(map_g_gp), 
                   weight_g_p  = 1/rowSums(map_g_p), 
                   SPO         = as.vector(df$PO))
    
  warmup=round((floor((iter*0.25)/(iter/10))*(iter/10)))
  thin=max(iter/1000,1)
  
  fileName <- stanMod

  stan_code <- readChar(fileName, file.info(fileName)$size)
  #cat(stan_code)
  resStan = NULL

  try(resStan <- stan(model_code = stan_code, data = standat,
                      chains = nchains, iter = iter, warmup = warmup, thin = thin))

  return(resStan)
}


run_stan_regress <- function(data_kmers, 
                             niter=1000,
                             n.thin = 1, 
                             seg_regions=NULL,
                             reg=NULL,
                             save_dir=NULL,
                             STZ=T, use_gene=F,#no_gene=F, 
                             no_theta=F, alpha=NULL,
                             cond=NULL,
                             stan=T, stanMod=NULL){
  reg <- list()
  for(i in sort(unique(data_kmers$RRIX))){
    testPups <- as.vector(t(data_kmers %>% ungroup() %>%
                              filter(RRIX == i) %>%            #, !Pup.ID %in% problemPups
                              dplyr::select("Pup.ID") %>% distinct()))
    testData <- data_kmers %>% filter(Pup.ID %in% testPups)
    chr = unique(testData$seq.Chromosome)
    terms <- c("dir", "Diet", "DietRIX","RRIX", "PODIETRIX")
    if(!is.null(cond)){
      testData$PORIX = testData$RIX
      testData$PODIETRIX = paste(testData$DietRIX, testData$dir, sep="_")
      testData$PODIETRIX = factor(testData$PODIETRIX, 
                                  levels=apply(data.frame(expand.grid(sort(unique(testData$Diet)), sort(unique(testData$RRIX)), sort(unique(testData$dir)))),
                                               1, function(x) paste(x,collapse="_")))
      getCond = table((testData %>% select(PODIETRIX, Pup.ID) %>% distinct())$PODIETRIX)
      
      terms <- c("dir", "Diet", "DietRIX","RRIX", "PODIETRIX")
      cond = "PODIETRIX"
      if(any(getCond == 1)){
        terms = cond = "RRIX"
      }
    }
    
    #colnames(testData) <- toupper(colnames(testData))
    if(!is.null(seg_regions) & chr!= "X"){
      use_region = seg_regions[[paste(i)]] %>% filter(Chr == chr)
      rem = c()
      for(j in 1:nrow(use_region)){
        tmp = use_region[j,]
        rem = c(rem, which(testData$seq.Position > (tmp$start) & testData$seq.Position < (tmp$end)))
      }
      if(length(rem)>0) testData = testData[-rem,]
    }
    #terms <- c("dir", "Diet", "DietRIX","RRIX", "PODIETRIX")
    if(length(which(!levels(testData$DietRIX) %in% testData$DietRIX)) > 0){
      newLev = levels(testData$DietRIX)[-which(!levels(testData$DietRIX) %in% testData$DietRIX)]
      testData$DietRIX = factor(testData$DietRIX, levels = newLev)
    }
    
    if(length(which(!levels(testData$seq.Gene) %in% testData$seq.Gene)) > 0){
      newLev = levels(testData$seq.Gene)[-which(!levels(testData$seq.Gene) %in% testData$seq.Gene)]
      testData$seq.Gene = factor(testData$seq.Gene, levels = newLev)
      testData$pup_gene = factor(paste(testData$Pup.ID, testData$seq.Gene, sep="_"))
    }
    
    encoded <- unique(getEncoding(testData, terms = unique(c(terms,"Pup.ID","seq.Gene","pup_gene"))))
    
    reg[[paste0("rix_",i)]] = list()
    if(stan){
      reg[[paste0("rix_",i)]]$stan = stan_ase(df=testData, encoded=encoded, 
                                              mu_g0=0.5, t20=0.001, 
                                              nchains=2, iter=niter, 
                                              C_diet_4=C_diet_4, C_diet_2=C_diet_2,C_diet_3=C_diet_3,
                                              terms = terms,
                                              stanMod=stanMod)
      summ = data.frame(summary(reg[[paste0("rix_",i)]]$stan)$summary)
      if(length(grep("eta|ind", rownames(summ))) > 0) summ = summ[-grep("eta|ind", rownames(summ)),]
      summ$param = rownames(summ)
      summ$param[grep("mu_g[[]", rownames(summ))] = encoded$Level[which(encoded$Variable == "SEQ.GENE")]
      summ$param[grep("mu_p[[]", rownames(summ))] = encoded$Level[which(encoded$Variable == "PUP.ID")]
      
      reg[[paste0("rix_",i)]]$sum = summ
    }
    
    freq_test = testData %>% group_by(seq.Gene) %>%
      summarize(y = sum(CC_1), N = sum(sum))
    freq_res = apply(freq_test, 1, function(x) 
      binom.test(as.numeric(x["y"]), as.numeric(x["N"]), p = 0.5))
    names(freq_res) = freq_test$seq.Gene
    reg[[paste0("rix_",i)]]$freq = freq_res
  }
  
  return(reg) 
}


############# encoding stuff #################
getEncoding <- function(df, terms){
  encoded <- data.frame()
  #if(length(grep("gene", terms)) > 0) terms = terms[-grep("gene", terms)]
  for(i in 1:length(terms)){
    
    var <- toupper(paste(terms[i]))
    ind <- match(var, toupper(colnames(df)))
    if(!is.na(ind) > 0){
      useLev = unique(df[,ind])
      if(!is.null(levels(df[,ind]))) useLev = levels(df[,ind])[which(levels(df[,ind]) %in% useLev)]
      
      #vec <- factor(t(df[,ind]), levels=useLev)
      len <- length(unlist(useLev))       #length(levels(as.factor(vec)))
      tempdf <- data.frame(Level = as.character(paste(unlist(useLev))), 
                           Index = 1:len, 
                           Variable = rep(var,len))
      encoded <- rbind(encoded,tempdf)
    }
  }
  return(encoded)
}

    