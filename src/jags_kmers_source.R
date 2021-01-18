library(tidyverse)
library(rjags)
library(coda)

CC_pairs <- c("CC001/011", "CC041/051", "CC017/004", "CC023/047", "CC026/006", "CC014/003", "CC035/062", "CC032/042", "CC005/040")

###################################

run_jags_regress <- function(data_kmers, 
                             niter=1000,
                             n.thin = 1, 
                             seg_regions=NULL,
                             reg=NULL,
                             save_dir=NULL,
                             STZ=T, use_gene=F,#no_gene=F, 
                             no_theta=F, alpha=NULL,
                             cond=NULL,
                             stan=F, stanMod=NULL){
  reg <- list()
  for(i in sort(unique(data_kmers$RRIX))){
    testPups <- as.vector(t(data_kmers %>% ungroup() %>%
                              filter(RRIX == i) %>%            #, !Pup.ID %in% problemPups
                              dplyr::select("Pup.ID") %>% distinct()))
    testData <- data_kmers %>% filter(Pup.ID %in% testPups)
    chr = unique(testData$seq.Chromosome)
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
    if(!stan){
      reg[[paste0("rix_",i)]] <- jags.genes.run(data=testData, mu_g0=0.5, niter=niter, n.thin=n.thin, 
                                                nchains=2, terms=terms, STZ=STZ, use_gene=use_gene, 
                                                no_theta=no_theta, alpha=alpha, quant=NULL, tau=NULL, 
                                                cond=cond, #no_gene=no_gene, 
                                                TIMBR=F, C_diet_4=C_diet_4, C_diet_2=C_diet_2,C_diet_3=C_diet_3)
    } else {
      if(length(which(!levels(testData$DietRIX) %in% testData$DietRIX)) > 0){
        newLev = levels(testData$DietRIX)[-which(!levels(testData$DietRIX) %in% testData$DietRIX)]
        testData$DietRIX = factor(testData$DietRIX, levels = newLev)
      }
      
      if(length(which(!levels(testData$seq.Gene) %in% testData$seq.Gene)) > 0){
        newLev = levels(testData$seq.Gene)[-which(!levels(testData$seq.Gene) %in% testData$seq.Gene)]
        testData$seq.Gene = factor(testData$seq.Gene, levels = newLev)
        testData$pup_gene = factor(paste(testData$Pup.ID, testData$seq.Gene, sep="_"))
      }
      reg[[paste0("rix_",i)]] = stan_ase(df=testData, encoded=NULL, 
                                         mu_g0=0.5, t20=0.001, 
                                         chains=2, iter=niter, 
                                         C_diet_4=C_diet_4, C_diet_2=C_diet_2,C_diet_3=C_diet_3,
                                         terms = NULL,
                                         alpha = NULL)
      
    }
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

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


array.summarize <- function(array, labels=NULL){
  if(class(array) == "list" & class(array[[1]]) %in% c("list", "mcarray")){   
    nme <- names(array)
    listLength <- length(array)
  } else {
    listed <- F
    nme <- NULL
    listLength <- 1
  }
  sum_list <- list()
  for(i in 1:listLength){
    comb_int = means = NULL
    len3 <- length(dim(as.array(array[[i]])))
    nchains <- dim(as.array(array[[i]]))[len3]
    if(len3 == 3){
      est_comb <- sapply(1:dim(array[[i]])[1], function(x){
        sapply(1:nchains, function(y){
          tmp = as.vector(unlist(array[[i]][x,,y]))
          tmp[!is.na(tmp)]}, simplify=F)
      }, simplify=F)
      sep_int <- lapply(est_comb, function(x){
        data.frame(lapply(x, function(y) HPDinterval(as.mcmc(y))))  })
    } else if (len3 == 4){
      est_comb <- sapply(1:dim(array[[i]])[1], function(x){
        sapply(1:dim(array[[i]])[2], function(y){
          sapply(1:nchains, function(z){
            tmp = as.vector(unlist(array[[i]][x,y,,z]))
            tmp[!is.na(tmp)]}, simplify=F)
        }, simplify=F) }, simplify=F)
      sep_int <- lapply(est_comb, function(x){
        do.call("rbind", lapply(x, function(y) 
          data.frame(lapply(y, function(z) HPDinterval(as.mcmc(z)))) ))  })
      comb_int <- lapply(est_comb, function(x) 
        do.call("rbind", lapply(x, function(y) HPDinterval(as.mcmc(unlist(y)))) ) )
      comb_int <- do.call("rbind", comb_int)
      means <- lapply(est_comb, function(x) 
        do.call("rbind", lapply(x, function(y) c(Mode(unlist(y)), mean(unlist(y)), median(unlist(y))))) )
      means <- do.call("rbind", means)
      colnames(means) <- c("Mode","Mean","Median")  
      
    } else {
      est_comb <- lapply(array, unlist)
      sep_int <- lapply(est_comb, function(x){
        data.frame(HPDinterval(as.mcmc(x)))  })
    }
    
    sep_int <- do.call("rbind", sep_int)
    colnames(sep_int)[1:2] <- c(paste0("lower.", nchains), paste0("upper.", nchains))
    
    #sep_ints <- sapply(1:dim(array)[1], function(x) HPDinterval(as.mcmc.list(array)), simplify=F)
    #names(sep_ints) <- paste0("ints",seq(1:nchains))
    #sep_int <- do.call("cbind", sep_ints)
    if(is.null(comb_int)){
      comb_int <- lapply(est_comb, function(x) HPDinterval(as.mcmc(unlist(x))))
      comb_int <- do.call("rbind", comb_int)
    }
    if(is.null(means)){
      means <- lapply(est_comb, function(x) c(Mode(unlist(x)), mean(unlist(x)), median(unlist(x))))
      means <- do.call("rbind", means)
      colnames(means) <- c("Mode","Mean","Median")
    }
    
    summary <- cbind(means, comb_int, sep_int)
    if(is.null(labels) | is.null(nme)){
      if(is.null(rownames(summary))) rownames(summary) <- seq(1:dim(array[[i]])[1])
    } else {
      v <- gsub("B_|PO", "", toupper(nme[i]))
      if(toupper(nme[i]) == "MU_P") v <- "PUP.ID"
      if(toupper(nme[i]) == "MU_G") v <- "PUP_GENE"
      if(toupper(nme[i]) == "MU_R") v <- "RRIX"
      if(toupper(nme[i]) == "B_DIETPORIX") v <- c("DIET","RRIX")
      v = ifelse(v == "RIX" & length(v) == 1, "RRIX", v)
      ind <- which(toupper(labels$Variable) %in% v)
      if(any(ind) == 0){
        if(is.null(rownames(summary))){
          rownames(summary) <- seq(1:dim(array[[i]])[1])
        } else {
          rownames(summary) <- gsub("var","", rownames(summary))
        }
      } else if(length(v) > 1){
        names = sapply(1:length(v), function(i)
          labels$Level[which(toupper(labels$Variable) %in% v[i])], simplify=F)
        names_ex = expand.grid(names[[2]], names[[1]])
        rownames(summary) = paste0(names_ex$Var1, "_", names_ex$Var2)
      } else {
        rownames(summary) <- unique(gsub(" ", "_", labels$Level[ind]))
      }
    }
    sum_list[[i]] <- summary
  }
  if(listLength == 1){
    sum_list <- sum_list[[1]]
  } else {
    names(sum_list) <- nme
  }
  
  return(summary = sum_list)
}

########################################################

jags.genes.run <- function(data, mu_g0=0.5, t20=0.001, n.thin=1, 
                           niter=1000, nchains=2, terms=NULL, quant=NULL, 
                           P0=1, Q0=1, STZ=TRUE, use_gene=F, gamma = F, 
                           no_theta=FALSE, alpha=NULL, tau=NULL, TIMBR=F, 
                           C_diet_4=NULL, C_diet_2=NULL,C_diet_3=NULL,
                           add_impute=NULL, cond=NULL){
  
  colnames(data) <- toupper(colnames(data))
  
  encodeKmers <- unique(getEncoding(df = data, term = unique(c(terms,"PUP.ID","SEQ.GENE","PUP_GENE"))))
  
  data$PUP.ID <- factor(data$PUP.ID, encodeKmers$Level[which(encodeKmers$Variable == "PUP.ID")])
  #data$RIX_DIR_DIET = paste(data$RRIX, data$DIR, data$DIET, sep="_")
  data %>% ungroup() %>% arrange(PUP.ID) -> data
  
  data %>% dplyr::select(one_of("PUP.ID")) %>%    
    left_join(data.frame(PUP.ID = unique(data$PUP.ID), 
                         INDP = seq(1:length(unique(data$PUP.ID)))), by="PUP.ID") %>%
    distinct() -> merge_indP
  
  data %>% 
    dplyr::select(one_of("SEQ.GENE","PUP.ID")) %>%
    distinct() %>%
    ungroup() %>%
    mutate(INDG=seq(n())) -> merge_indG
  
  if(!is.null(cond)){
    if(!is.null(quant)){
      if(quant%in% cond) cond = cond[-match(quant, cond)]
    }
    data %>%
      dplyr::select(one_of(paste(cond))) %>%
      arrange(get(cond)) %>%
      #mutate(RIX_DIR_DIET = paste(RRIX, DIR, DIET, sep="_")) %>%
      distinct() %>%
      ungroup() %>%
      mutate(INDC=seq(n())) %>%
      right_join(data[,c("PUP.ID", paste(cond))]) %>%
      distinct() -> merge_indC
    
    data = data %>% left_join(merge_indC, by = c("PUP.ID", paste(cond)))
  }
  
  merge_indG$indP = as.numeric(merge_indG$PUP.ID)
  
  data %>% group_by(PUP.ID, SEQ.GENE) %>%
    mutate(indK=seq(n())) %>%
    group_by(PUP.ID) %>%
    left_join(merge_indP, by = c("PUP.ID")) %>% 
    left_join(merge_indG, by = c("SEQ.GENE", "PUP.ID")) %>%
    dplyr::select(-one_of(toupper(c("seq.ProbeSeq", "seq.TranscriptID", "seq.end5", "seq.end3", "seq.consensus", 
                                    "seq.Chromosome", "seq.Position", "seq.refseq", "seq.altseq", "seq.alt", "seq.ref",
                                    "CC_1","CC_2")))) %>%
    arrange(PUP.ID, SEQ.GENE) %>%
    ungroup() -> data
  
  
  use_PORIX=use_diet=use_RIX=two_diets=F
  singleRIX=singlePO=T
  if(length(unique(data$RRIX)) > 1){
    use_RIX = T
    singleRIX = F
    terms = unique(c(terms, "RRIX"))
  }
  
  data %>% group_by(SEQ.GENE, PUP.ID) %>%
    summarize(sum_fin1 = sum(FIN1), sum_fin2 = sum(FIN2), sumTot = sum(SUM)) %>%
    arrange(PUP.ID, SEQ.GENE) %>%
    left_join(merge_indP, by="PUP.ID") %>%
    left_join(merge_indG, by = c("SEQ.GENE", "PUP.ID")) -> gene_tot
  
  if(!is.null(cond)){
    gene_tot = gene_tot %>% left_join(merge_indC, by = c("PUP.ID")) 
  }
  
  if(use_gene){
    y_gk = gene_tot$sum_fin1
    N_gk = gene_tot$sumTot
    indG = gene_tot$INDG
    indP = gene_tot$INDP
  } else {
    y_gk <- data$FIN1
    N_gk <- data$SUM
    indG <- data$INDG
    indP = data$INDP
    indP = gene_tot$INDP
  }
  
  if(!is.null(cond)){
    indC = merge_indC$INDC[unique(indP)]
  } else {
    indC = vector(rep(1, length(unique(indP))))
  }
  
  nG <- length(unique(indG))
  nP <- length(unique(data$PUP.ID))
  nC <- length(unique(indC))
  
  
  mu_0 <- mean(y_gk/N_gk)
  tau_0 <- 1/var(y_gk/N_gk)
  
  betaStuff <- c()
  beta_0 <- list(mu_r = rep(0.5, length(unique(data$RRIX))))  
  
  dataStuff = list("y_gk" = y_gk, "N_gk" = N_gk, 
                   "nP" = nP, "indP" = indP)
  if(!is.null(cond)){
    dataStuff$nC = nC
    dataStuff$indC = indC
  }
  if(!use_gene){
    dataStuff$nG = nG
    dataStuff$indG = indG
  }
  
  indArray = NULL
  if(!is.null(terms)){
    data %>% dplyr::select(one_of("PUP.ID",toupper(terms),quant)) %>%
      dplyr::select(-contains("PUP_GENE")) %>%
      distinct() -> index
    indexPup = index
    indArray <- as.data.frame(sapply(colnames(index), function(x){
      as.numeric(factor(as.vector(t(index[,paste(x)])), 
                        levels=encodeKmers$Level[which(encodeKmers$Variable == x)]))
    }))
    indArrayPup=indArray
    if(!is.null(cond)){
      indArray = indArray %>% select(-"PUP.ID") %>% distinct()
      index = index %>% select(-"PUP.ID") %>% distinct()
    }
    
    if(length(quant) > 0) indArrayPup[,quant] = indexPup[,quant]
    
    if(length(grep("PO|DIR", toupper(terms))) > 0) use_PORIX=T
    if(length(grep("DIET", toupper(terms))) > 0) use_diet=T
    if(length(grep("RIX", toupper(terms))) > 0) use_RIX=T
    
    if(use_PORIX){
      PO <- indArray$DIR
      if(STZ){
        PO <- ifelse(PO == 1, 0.5, -0.5)
      } else {
        PO <- ifelse(PO == 1, 0, 1)
      }
      
      if(length(unique(PO)) == 1) use_PORIX = F
      
      indPORIX = indRIX = indArrayPup$RRIX
      if(is.null(indPORIX)) indPORIX = rep(1, length(PO))
      nRIX = length(unique(indPORIX))
      
      singleRIX = ifelse(nRIX == 1, T, F)
      use_RIX = F           
      
      betaStuff <- c(betaStuff, "b_PORIX")
      dataStuff <- c(dataStuff, list("PO" = PO, "nRIX" = nRIX, "indRIX" = indRIX))
    }
    
    
    
    if(!use_RIX & singleRIX) {
      indArrayPup$RRIX <- rep(1, nrow(indArrayPup))
      dataStuff$indRIX = indArrayPup$RRIX
      dataStuff$nRIX = 1
    } else if(use_RIX){
      dataStuff$indRIX = indArrayPup$RRIX
      dataStuff$nRIX = length(unique(dataStuff$indRIX))
    }
    
    if(use_diet){
      indDiet = indArray$DIET
      nDiet = length(unique(indDiet))
      
      if (any(is.na(indDiet))) nDiet = length(unique(indDiet[-which(is.na(indDiet))]))
      
      if(nDiet < 3) two_diets = T
      
      X_diet = matrix(0, nrow=length(indDiet), ncol=nDiet)
      for(i in 1:length(indDiet)){
        X_diet[i,indDiet[i]] <- 1
      }
      if(TIMBR){
        require(TIMBR)
        C_diet = TIMBR:::sumtozero.contrast(nDiet)
      } else {
        C_diet = as.matrix(get(paste0("C_diet_",nDiet)))
      }
      XC <- X_diet %*% C_diet
      betaStuff <- c(betaStuff, "b_diet")  
      
      if(STZ) {
        beta_0 <- c(beta_0, list("a_diet"=matrix(0, ncol=1, nrow=(nDiet-1)))) 
        dataStuff <- c(dataStuff, list("XC" = XC, "C_diet" = C_diet, "nDiet"=nDiet))    
        
      } else {
        beta_0 <- c(beta_0, list("b_diet"=c(NA, rep(0, nDiet-1))))    
        dataStuff <- c(dataStuff, list("indDiet"=indDiet, "nDiet"=nDiet))
      }
    }
    
    if(all(use_diet, use_PORIX)){
      
      if(STZ){
        indDietRIX = cbind(indArray$DIET, indArray$RRIX)    
        nDietRIX = prod(apply(indDietRIX, 2, max))
        beta_0$a_dietPORIX=matrix(0,nDiet-1, nRIX)
      } else {
        indDietRIX = indArray$DIETRIX
        nDietRIX = length(unique(indDietRIX))
        beta_0$b_dietPORIX=c(NA, rep(0, nDietRIX-1))
      }
      
      betaStuff <- c(betaStuff, "b_dietPORIX")    
    }
    #if(!identical(unique(data$PUP.ID), unique(index$PUP.ID))) warning("Check Pup.ID order")
  }
  use_impute = NULL
  
  if(!is.null(quant)){
    D = as.vector(unlist(indexPup[,quant]))
    
    beta_0$b_D = 0
    betaStuff <- c(betaStuff, "b_D")
    
    if(any(is.na(D))){
      dataStuff$D_temp = D + 0.5
      use_impute = add_impute
      betaStuff <- c(betaStuff, "D")
    } else {
      dataStuff$D = D
    }
    
  }
  useModel <- makeSTZModel(diet=use_diet, RIX=use_RIX, PORIX=use_PORIX, 
                           singleRIX, two_diets, use_gene, no_theta, alpha, tau, cond,
                           STZ, quant,add_impute = use_impute)
  alphastr <- ifelse(!no_theta, "theta", 
                     ifelse(length(grep("dgamma", alpha))>0, "dgamma", 
                            ifelse(length(grep("dunif", alpha))>0, "alpha", 
                                   ifelse(is.numeric(alpha), as.character(alpha), 
                                          paste(unlist(strsplit(as.character(alpha)," ","")), collapse="")))))
  write.table(useModel, paste0("makeModelTemp_",alphastr,niter,".txt"), quote = F, row.names = F, col.names = F)
  
  # check: cat(makeSTZModel(use_diet, use_RIX, use_PORIX, singleRIX, tau="sfsdf"))
  
  jags.params <- c("mu_p","mu_r",betaStuff)   #"tau2",
  
  if(gamma) jags.params <- c("P", "Q", "mu_g")
  jags.init <- list()     #"tau2"=tau_0
  
  if(!use_gene){
    jags.params = c(jags.params, "mu_g")
    jags.init$mu_g = rep(mu_g0, nG)
  }
  if(!is.null(beta_0)) jags.init<- c(jags.init, beta_0)
  
  
  if(gamma) jags.init <- list("P"=rep(P0,nP), "Q"=rep(Q0,nP), "mu_g"=rep(mu_g0, nG))
  if(length(grep("alpha~", gsub(" " ,"", useModel))) > 0){
    jags.params <- c(jags.params, "alpha","alpha0","mu_c") 
    jags.init[["alpha"]] = tau_0        #tau_0     #t20
    jags.init[["alpha0"]] = 10        #rep(tau_0,nP)
  }
  
  jags <- jags.model(paste0("makeModelTemp_",alphastr,niter,".txt"),
                     data = dataStuff, 
                     inits = jags.init, 
                     n.chains = nchains,
                     n.adapt = niter*0.1)
  
  update(jags, niter)
  
  res<- jags.samples(jags,
                     thin = n.thin, 
                     jags.params,
                     niter)
  
  
  if(gamma){
    mu.est <- as.mcmc.list(res$P/(res$P+res$Q)) 
  } else {
    mu.est <- res$mu_p 
  }
  if(gamma){
    alpha.est <- as.mcmc.list(res$P+res$Q)      
  } else {
    alpha.est <- res$alpha 
  }
  covarSummary <- array.summarize(res, encodeKmers)
  
  return(list(res=res, data = list(data = dataStuff, inits = jags.init),
              codes = encodeKmers, indices = indArray, 
              model=useModel, summary = covarSummary))
}


################## STZ model ########################

makeSTZModel <- function(diet=T, RIX=T, PORIX=T, singleRIX=F, 
                         two_diets=F, use_gene=F, no_theta=F, 
                         alpha=NULL,tau=NULL,cond=NULL, 
                         STZ=T, quant=NULL,
                         add_impute = NULL){
  if(is.numeric(alpha)){
    alpha_str  <- paste0("alpha=", alpha)
    alpha0_str <- paste0("alpha0=", alpha)
  } else if(is.character(alpha)){
    alpha_str  = alpha
    alpha0_str = alpha
  } else {
    #alpha_str = "alpha ~ dunif(1E-3,1E3)"
    alpha_str  = "alpha ~ dgamma(0.01, 0.01)"
    alpha0_str = "alpha0 ~ dunif(1,1000)"
  }
  if(is.numeric(tau)){
    tau_str <- paste0("tau2=", tau)
  } else if(is.character(tau)){
    tau_str = tau
  } else {
    #tau_str = "tau2 ~ dunif(1E-3,1E3)"
    tau_str = "tau2 ~ dgamma(0.01, 0.01)"
    
  }
  baseMod_kmer <- "model {
  for(i in 1:length(y_gk)){
    y_gk[i] ~ dbin(mu_g[indG[i]], N_gk[i])
  }\n"
  baseMod_kmer_ng <- "model {
  for(i in 1:length(y_gk)){
    y_gk[i] ~ dbin(mu_p[indP[i]], N_gk[i])
  }\n"
  baseMod_gene <- paste("for(g in 1:nG){
    mu_g[g] ~ dbeta((mu_p[indP[g]]*alpha), ((1-mu_p[indP[g]])*alpha)) T(0.001,0.999)                      
  }", alpha_str, sep="\n")
  
  baseMod_gene_ng = ""
  
  baseMod_rix = "\nfor(r in 1:nRIX){
  mu_r[r] ~ dbeta(1,1) T(0.001,0.999)
  b0[r] = log(mu_r[r]/(1-mu_r[r]))
  }\n"
  
  baseMod_cond = "\nfor(c in 1:nC){
  mu_c[c] = exp(eta_c[c]) / (1+exp(eta_c[c]))
  eta_c[c] <- b0[indRIX[c]]"
  #  eta_c[c] = log(mu_c[c] / (1-mu_c[c]))
  
  
  modRIX <- "+ b_RIX[indRIX[p]]" 
  modPORIX <- "+ PO[p]*b_PORIX[indRIX[p]]"
  modDietPORIX <- "+ PO[p]*(XC[p,] %*% a_dietPORIX[,indRIX[p]])"
  modD <- ifelse(is.null(quant), "", "+ D[p]*b_D")
  baseMod_addtau = paste0("\n}\n",tau_str)
  
  if(STZ){
    modDiet <- "+ XC[p,] %*% a_diet"
    modDiet2 <- "\nprec_diet ~ dgamma(0.01, 0.01)   
  for(i in 1:(nDiet-1)){
    a_diet[i,1] ~ dnorm(0, prec_diet)
  }\n"
    modDiet2 <- "\nfor(i in 1:(nDiet-1)){
    a_diet[i,1] ~ dnorm(0, 0.0001)
  }\n"
  } else {
    modDiet <- "+ b_diet[indDiet[p]]"
    modDiet2 <- "\nb_diet[1] <- 0
  prec_diet ~ dgamma(0.01, 0.01)   
  for(i in 2:nDiet){
    b_diet[i] ~ dnorm(0, prec_diet)
  }\n"
    modDiet2 <- "\nb_diet[1] <- 0
  for(i in 2:nDiet){
    b_diet[i] ~ dnorm(0, 0.0001)
  }\n"
  }
  
  
  baseMod_pup = "\nfor(p in 1:nP){
  mu_p[p] <- exp(eta[p])/(1+exp(eta[p]))
  eta[p] ~ dnorm(theta[p], tau2)
  theta[p] <- b0[indRIX[p]]"
  
  
  modDiet3 = modDiet3_TD = modDietPORIX3 = modDietPORIX3_TD = ""
  
  if(diet & STZ){
    modDiet3 <- "b_diet <- C_diet %*% a_diet \n"
    modDiet3_TD <- "b_diet <- C_diet * a_diet \n"
    if(PORIX){
      modDietPORIX3 <- "b_dietPORIX <- C_diet %*% a_dietPORIX \n"
      modDietPORIX3_TD <- "b_dietPORIX <- C_diet * a_dietPORIX \n"
    }
  }
  
  #baseMod1_ng <- paste0("model {\n",
  #tau_str, 
  #"\nb0 ~ dnorm(0, 0.001)
  #for(p in 1:nP){
  #  eta[p] ~ dnorm(theta[p], 1/tau_2)
  #  theta[p] <- b0")
  
  baseMod1_nt <- paste0("model {
  for(g in 1:nG){
    mu_g[g] ~ dbeta((mu[indP[g]]*alpha), ((1-mu[indP[g]])*alpha)) T(0.001,0.999)                      
  }\n", alpha_str,
                        "\nb0 ~ dnorm(0,0.0001)
  for(p in 1:nP){
    eta[p] <- b0")
  
  baseMod1_nt <- paste0("model {
  for(g in 1:nG){
    mu_g[g] ~ dbeta((mu[indP[g]]*alpha)+1, ((1-mu[indP[g]])*alpha)+1) T(0.001,0.999)                      
  }\n", alpha_str,#b0 ~ dnorm(0,1)
                        "\n", tau_str,
                        "\nmu_r ~ dbeta(1,1 T(0.001,0.999))",
                        "\nb0 = log(mu_r/(1-mu_r))")
  
  
  modRIX2 <- "prec_RIX ~ dgamma(0.01, 0.01)
  for(i in 1:nRIX){
    b_RIX[i] ~ dnorm(0, prec_RIX)
  }\n"
  
  modRIX2 <- "for(i in 1:nRIX){
    b_RIX[i] ~ dnorm(0, 0.0001)
  }\n"
  
  modPORIX2 <- "prec_PORIX ~ dgamma(0.01, 0.01)
  for(i in 1:nRIX){
    b_PORIX[i] ~ dnorm(0, prec_PORIX)
  }\n"
  
  modPORIX2 <- "for(i in 1:nRIX){
    b_PORIX[i] ~ dnorm(0, 0.0001)
  }\n"
  
  modDietPORIX2 <- "prec_dietPORIX ~ dgamma(0.01, 0.01)
  for(r in 1:nRIX){
    for(i in 1:(nDiet-1)){
      a_dietPORIX[i,r] ~ dnorm(0, prec_dietPORIX)
    }
  }\n"
  modDietPORIX2 <- "for(r in 1:nRIX){
    for(i in 1:(nDiet-1)){
      a_dietPORIX[i,r] ~ dnorm(0, 0.0001)
    }
  }\n"
  
  
  single_modPORIX2 <- "b_PORIX ~ dnorm(0, 0.0001)\n"
  
  if(singleRIX){
    modRIX = modRIX2 = ""
    modPORIX2 = single_modPORIX2
  }
  
  if(two_diets){
    modDiet3 = modDiet3_TD
    modDietPORIX3 = modDietPORIX3_TD
  }
  if(!diet) modDiet = modDiet2 = ""
  if(!RIX) modRIX = modRIX2 = ""
  if(!PORIX) modPORIX = modPORIX2 = ""
  if(any(c(!diet, !PORIX))){        #, !(RIX|singleRIX), , singleRIX
    modDietPORIX = modDietPORIX3 = modDietPORIX2 = ""
  } 
  if(use_gene){
    baseMod_kmer = baseMod1_kmer_ng
    baseMod_gene = baseMod_gene_ng
  }
  if(no_theta){
    baseMod1 = baseMod1_nt
  }
  baseMod_D = ""
  if(!is.null(quant)){
    baseMod_D = "\nprec_D ~ dgamma(0.01, 0.01)   
      b_D ~ dnorm(0, prec_D)"
    baseMod_D = "\nb_D ~ dnorm(0, 0.0001)"
  }
  
  
  linReg = paste0(modDiet, modRIX, modPORIX, modDietPORIX, modD)
  
  if(!is.null(cond)){
    baseMod_pup = paste("\nfor(p in 1:nP){
  mu_p[p] ~ dbeta((mu_c[indC[p]]*alpha0), ((1-mu_c[indC[p]])*alpha0)) T(0.001,0.999) 
  }", alpha0_str, sep="\n")
    
    linReg = gsub("[p]","c",linReg)
    #linReg = ""
    baseMod_cond = paste(baseMod_cond, linReg,"\n}")
    baseMod_addtau = linReg = ""
    #modDiet2 = modDiet3 = modPORIX2 = modDietPORIX2 = modDietPORIX3 = baseMod_D = ""
  } else {
    baseMod_cond = ""
    
  }
  
  return(paste0(baseMod_kmer, add_impute, baseMod_gene, 
                baseMod_pup, linReg, baseMod_cond, 
                baseMod_addtau,     #modRIX,
                modDiet2, modDiet3, modPORIX2, modDietPORIX2, modDietPORIX3, baseMod_D, baseMod_rix, #modRIX2, 
                "\n}"))
}

