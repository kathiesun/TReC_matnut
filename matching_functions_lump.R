source("matching.R")
source("./matnut/summary_functions.R")
source("./lm/formulaWrapper.R")
source("./matnut/prediction_functions.R")
source("./matnut/tables_n_plots.R")
options(mc.cores = parallel::detectCores())
library(bayesplot)
library(rstan)
############
# Matching #
############
match.multimp <- function(data, matchon, matchoff, idcol="ID", N=10, 
                          randvar = NA, fixvar = NA, tryLam = 1, Nclust = 2, 
                          chains=1, n.iter=30000, thin=10, clust = NULL, 
                          encoded=NA, allparam=NA, ptypes=NA, 
                          sq=F, fixPO = T, fixClust = T, p="p", beta=T){
  matched_df <- list()
  multimp <- list()
  plots = list()
  
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
  
  match_all <- matching$generateEveryMatch(df=df, matchon=matchon, matchoff=matchoff, idcol=idcol)
  match_all$out$PO = df[["PO"]][match(match_all$out$ID.1, df[[idcol]])]
  
  for(i in 1:length(ptypes)){
    if(length(which(is.na(df[,ptypes[i]]))) > 0){
      df <- df[-which(is.na(df[,ptypes[i]])), ]
    }
  
    
    tot_pairs = df %>% group_by(get(matchon)) %>%
      summarize(pos = length(which(PO == 0.5)), neg = length(which(PO == -0.5))) %>%
      mutate(prod = pos*neg) %>%
      summarize(sum=sum(prod))
    which_larger = df %>% mutate(PO = factor(PO)) %>%
      group_by(RIX, PO) %>%
      summarize(mean = mean(get(ptypes[i]), na.rm=T),
                med = median(get(ptypes[i]), na.rm=T))
    larger_PO = unlist(lapply(unique(which_larger$RIX), function(x){
      tmp = filter(which_larger, RIX == x)
      ifelse(tmp$med[which(tmp$PO == 0.5)] <= tmp$med[which(tmp$PO == -0.5)],
             ifelse(tmp$mean[which(tmp$PO == -0.5)] > tmp$mean[which(tmp$PO == 0.5)], -0.5, 0.5), 0.5)
    }))
    names(larger_PO) = unique(which_larger$RIX)
    
    firstMatch_fix = lapply(names(larger_PO), function(x){
      tmp = filter(match_all$out, RIX.1 == x)
      out = tmp
      if(larger_PO[[paste(x)]] == -0.5){
        out$PO = "-0.5"
      } else {
        out$ID.2 = tmp$ID.1
        out$ID.1 = tmp$ID.2
      }
      return(out)
    })
    firstMatch_fix = do.call("rbind",firstMatch_fix)
    
    if(!fixPO){
      firstMatch_fix = match_all$out
    }
    
    firstPass <- matching$getDeltaForPhen(original=df, firstMatch_fix, idcol="ID", phen=ptypes[i])
    
    firstPass = firstPass %>% group_by(RIX.1) %>%
      mutate(norm = scale(y)) %>%
      summarize(mean_y = mean(y, na.rm=T),mean_norm = mean(norm, na.rm=T),
                sq_mean_y = mean_y^2, sq_mean_norm = mean_norm^2)
 
    hiRIX = firstPass$RIX.1[which.max(firstPass$sq_mean_norm)]
    hiRIXind = encoded$Index[match(hiRIX, encoded$Level[which(encoded$Variable == "RIX")])]
    findLo = firstPass$sq_mean_norm
    findLo[which(firstPass$RIX.1 == hiRIX)] = NA
    loRIX = firstPass$RIX.1[which.min(findLo)]
    loRIXind = encoded$Index[match(loRIX, encoded$Level[which(encoded$Variable == "RIX")])]
    
    if(is.null(clust)){
      clust = rep(NA,length(unique(df$RIX)))
      clust[hiRIXind]=1 # highest value assigned to alt cluster, mu_a
      clust[loRIXind]=0 # smallest value assigned to null cluster, mu_0
      if(!fixClust){
        clust = rep(NA,length(unique(df$RIX)))
      }
    }
    
    matched_df[[ptypes[i]]] <- list()
    multimp[[ptypes[i]]] <- list()
    multimp[[ptypes[i]]]$trans = list()
    matched_jags = fit_stan = list()
    
    #covariates <- allparam[-grep(paste("tau","sig",matchoff, sep="|"), allparam)]
    
    for(j in 1:N){
      match_init <- matching$generateMatching(df=df, matchon=matchon, matchoff=matchoff, idcol=idcol)
      vars = ifelse(any(is.na(c(randvar, fixvar))), c(randvar, fixvar)[which(!is.na(c(randvar, fixvar)))],
                    c(randvar, fixvar))
      covariates = setdiff(vars, unlist(strsplit(colnames(match_init$out), "[.]"))[c(T, F)])
      
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
      
      matched[[paste0("PO",".",1)]] = df[["PO"]][match(matched$ID.1, df[[idcol]])]
      keepInd <- c("ID.1","ID.2","DamID.1","DamID.2")         
      colnames(matched)[-which(colnames(matched) %in% keepInd)] <- 
        unlist(strsplit(colnames(matched)[-which(colnames(matched) %in% keepInd)],"[.]"))[c(T, F)]
      
      if(length(duplicated(colnames(matched)))>0){
        matched <- matched[,-which(duplicated(colnames(matched)))]
      } 

      matched_fix = lapply(names(larger_PO), function(x){
        tmp = filter(matched, RIX == x)
        out = tmp
          if(larger_PO[[paste(x)]] == -0.5){
            out$PO = "-0.5"
          } else {
            out$ID.2 = tmp$ID.1
            out$ID.1 = tmp$ID.2
          }
        return(out)
      })
      matched_fix = do.call("rbind",matched_fix)
      
      if(!fixPO){
        matched_fix = matched
      }
      
      temp_match <- matching$getDeltaForPhen(original=df, matched_fix, idcol="ID", phen=ptypes[i])
      #temp_match[,ptypes[i]] <- temp_match$y
      temp_match$pair = paste(temp_match$ID.1, temp_match$ID.2, sep="_")
      temp_match = temp_match %>% ungroup() %>% mutate(norm_y = (y-mean(y))/sd(y)) 
      
      check_sums = temp_match %>% group_by(RIX) %>% summarize(sum_y = sum(y), sum_norm = sum(norm_y))
      apply(check_sums[,2:3], 2, sum)
      
      use_match = temp_match  
      y.mat = data.frame(y = as.matrix(use_match$y))
      #colnames(y.mat) = ptypes[i]
      
      lmerobj <- BC.model(y.mat = y.mat, data=use_match, indvariable= "~ 1", 
                          transformParams = getMatnutTransformParams(tryLam = tryLam, normd = T))
      use_match = use_match %>% rename("orig_y"="y")
      use_match = cbind(use_match, data.frame(lmerobj$y.transform[[1]]))
      
      matched_df[[ptypes[i]]][[j]] = use_match
      
      
      y = use_match$y         #orig_y
      n = length(y)
      indRIX = encoded$Index[match(use_match$RIX, encoded$Level[which(encoded$Variable == "RIX")])]
      nRIX = length(unique(use_match$RIX))
      secondPass = use_match %>% group_by(RIX) %>% summarize(means = mean(y))
      #sign = ifelse(firstPass$mean_norm  > 0, 1, -1)
      
      dataList = list(
        y = y ,
        n = n ,
        nRIX = nRIX , 
        indRIX = indRIX ,
        clust = clust,
        sign = rep(NA,length(unique(df$RIX)))
      )
      if(!beta){
        dataList$Nclust = Nclust 
        dataList$onesRepNclust = rep(1,Nclust)
        paramUse = c("clust","muOfClust","tau","pClust")
        dataList$clust = dataList$clust+1
      } else {
        paramUse = c("clust","muOfRIX","muOfRIXsq","mu_a","mu_a_sq","tau", "mu_a_abs")
        if(!is.numeric(p)){
          paramUse = c(paramUse, "p")
        }
      }
      
      if(fixPO){
        sq=F
      } else {
        sq=T
      }
      if(fixClust){
        tau_mu0 ="1e-2"
        tau_mua ="1e-2"
      } else {
        tau_mu0 ="1e0"
        tau_mua ="1e0"
      }
      modelFin = make_model(tau_mu0, tau_mua, beta=beta, sq=sq, p=p)
      reg.jags <- jags.model(textConnection(modelFin), data=dataList, n.chains = 1, n.adapt = round(n.iter/5))
      update(reg.jags, n.iter=round(n.iter/5))
      fit_jags <- coda.samples(reg.jags, variable.names = paramUse, thin=thin, n.iter=n.iter)
      
      nK = 3
      #######################
      initf <- function(chain_id = 1) {
        list(sigma = 1, 
             mu_a = c(0.01),
             p = matrix(1/nK, nRIX, nK),
             mu_a_sq = 0.01^2,
             mu_a_abs = 0.01)
             #p = c(1/2, 1/2),
             #muOfRIX = rep(0, dataList$nRIX))
      } 
      
      n_chains <- 1
      init_ll <- lapply(1:n_chains, function(id) initf(chain_id = id))
      fit_stan[[j]] = stan(file="matnut/phenModel_discrete.stan", data=c(dataList, nK=nK), init = init_ll,
                           iter=n.iter, chains=n_chains, thin=thin, cores = 2)
      #post_stan <- as.array(fit_stan)
      #relab = relabel(dataList=dataList, fit=fit, tau_mu0, tau_mua)
      #relabs = relabel_stephens(dataList, fit)
      #fit[[1]] = as.mcmc(data.frame(apply(relab$mu_out, 2, as.mcmc)))
      matched_jags[[j]] = as.mcmc.list(fit_jags)
      
      multimp[[ptypes[i]]]$trans[[j]] = lmerobj
      #multimp[[ptypes[i]]]$clusters[[j]] = relab$clusters
    }
    #browser()
    
    full_list = lapply(unlist(dimnames(matched_jags[[1]][[1]])), function(x) 
      as.mcmc.list(lapply(matched_jags, function(y) as.mcmc(y[[1]][,x]))))
    names(full_list) = unlist(dimnames(matched_jags[[1]][[1]]))
  
    #lapply(full_list,summary)
    multimp[[ptypes[i]]]$fit = full_list
    multimp[[ptypes[i]]]$model = modelFin
    
  }
 
  if (i==1){
    multimp <- multimp[[ptypes[i]]]
    matched_df <- matched_df[[ptypes[i]]]
  }
  return(list(multimp=multimp, matched_df=matched_df, fit_stan=fit_stan))
}   

make_model = function(tau_mu0, tau_mua, beta=F, sq=F,
                      p="p"){
  m1 = "model {
  ## Likelihood
  for (i in 1:n){
  mu[i] <- muOfRIX[ indRIX[i] ]
  y[i]  ~ dnorm( mu[i], tau ) 
  }"
  
  if(!beta){
    m2 = "for (r in 1:nRIX) {
    muOfRIX[r] <- muOfClust[ clust[r] ]
    clust[r] ~ dcat( pClust[1:Nclust] )
    }"
    m3 = "## Prior
    for (c in 1:Nclust) {
    muOfClust[c] ~ dnorm( 0, 1e-3 )
    }
    
    pClust[1:Nclust] ~ ddirich( onesRepNclust )"
  } else {
    p_pri = ifelse(is.numeric(p), "","p ~ dbeta( 1,1 )")
    
    if(sq){
      m2 = paste("for (r in 1:nRIX) {
    muOfRIXsq[r] = clust[r]*mu_a_sq
    muOfRIX[r] = pow( muOfRIXsq[r], 0.5) * sign[r]
    clust[r] ~ dbin(",p,", 1)      #
    sign_tmp[r] ~ dbin(0.5, 1)
    sign[r] = (sign_tmp[r] - 0.5) * 2
    }")
      m2 = paste("for (r in 1:nRIX) {
    muOfRIXsq[r] = clust[r]*mu_a_sq
    #muOfRIXsq[r] = pow( muOfRIX[r], 2)
    #muOfRIX[r] = clust[r]*mu_a
    muOfRIX[r] = pow( muOfRIXsq[r], 0.5)
    clust[r] ~ dbin(",p,", 1)      #
    }")
      m3 = paste0("## Prior\n",
                  p_pri,
                  "\nmu_0 = 0 
    mu_a ~ dnorm( 0, ", tau_mua,")
    mu_a_sq = pow( mu_a, 2)
    mu_a_abs = pow( mu_a_sq, 0.5)")      #
    } else {
      m2 = paste("for (r in 1:nRIX) {
      muOfRIX[r] = mu_0*(1-clust[r]) + clust[r]*mu_a
      clust[r] ~ dbin(",p,", 1)       #0.95
    }")
      m3 = paste0("## Prior\n",
                  p_pri,
    "\nmu_0 ~ dnorm( 0, ", tau_mu0," )
    mu_a ~ dnorm( 0, ", tau_mua," ) T(mu_0, )")
    }
  }
  m4 = "tau ~ dgamma( 0.01,0.01 )
}"
  
  m = paste(m1, m2, m3, m4, sep="\n")
  return(m)
}

#fit = list()
#fit[[1]] = list(mu_0 = matched_results[[1]]$multimp$fit$mu_0[[1]], 
#                mu_a = matched_results[[1]]$multimp$fit$mu_a[[1]],
#                tau = matched_results[[1]]$multimp$fit$tau[[1]])
relabel = function(fit, dataList, tau_mu0, tau_mua){
  require(mclust)
  #thetas = kmeans(c(fit[[1]][,"mu_0"], fit[[1]][,"mu_a"]), 3)
  thetas = Mclust(c(fit[[1]][,"mu_0"], fit[[1]][,"mu_a"]), G=3, model="V")
  means = thetas$parameters$mean
  clust = thetas$classification
  km = ifelse(any(thetas$parameters$variance$sigmasq > 0.5), T, F)
  
  if(km){
    thetas = kmeans(c(fit[[1]][,"mu_0"], fit[[1]][,"mu_a"]), 3)
    means = thetas$centers
    clust = thetas$cluster
  }
  
  mu_0 = means[which.min(abs(means))]
  mus_a  = sort(setdiff(means, mu_0))
  all_est = data.frame(est = c(fit[[1]][,"mu_0"], fit[[1]][,"mu_a"]), clust = clust)

  #ggplot(all_est, aes(x=est, group=clust, color=as.factor(clust))) + geom_line(stat="density")
  L_mat = sapply(1:dataList$nRIX, function(r){
    y = dataList$y[which(dataList$indRIX == r)]
    sd_use = mean(fit[[1]][,"tau"])^-0.5
    unlist(lapply(c(mus_a[1], mu_0, mus_a[2]), function(mu) {
      -log(prod(unlist(lapply(y, function(n)
        dnorm(n, mean=mu, sd=sd_use)))))
    }))
  })
  
  rix_clus = apply(L_mat, 2, which.min)
  rix_clus[which(rix_clus == 3)] = 1
  rix_clus[which(rix_clus == 2)] = 0
  #if(length(unique(rix_clus)) == 1){
  #  mu_out = cbind(unlist(fit[,"mu_0"]), unlist(fit[,"mu_a"]),
  #                 unlist(fit[,"p"]), unlist(fit[,"tau"]))
  #  colnames(mu_out) = c("mu_0", "mu_a","p","tau")
  #} else {
    sd0 = 1/(as.numeric(tau_mu0)^0.5)
    sda = 1/(as.numeric(tau_mua)^0.5)
    sd0 = sda = unlist(fit[[1]][,"tau"])^-0.5
    clus_probs_1 = apply(cbind(dnorm(unlist(fit[[1]][,"mu_0"]), mu_0, sd0), 
                               apply(cbind(dnorm(unlist(fit[[1]][,"mu_a"]), mus_a[1], sda),
                                           dnorm(unlist(fit[[1]][,"mu_a"]), mus_a[2], sda)),1,max)
    ), 1, prod)
    clus_probs_2 = apply(cbind(dnorm(unlist(fit[[1]][,"mu_a"]), mu_0, sd0), 
                               apply(cbind(dnorm(unlist(fit[[1]][,"mu_0"]), mus_a[1], sda),
                                           dnorm(unlist(fit[[1]][,"mu_0"]), mus_a[2], sda)),1,max)
    ), 1, prod)
    which_clus = apply(cbind(-log(clus_probs_1), -log(clus_probs_2)), 1, which.min)
    
    mu_out = data.frame(t(sapply(1:length(which_clus), function(t) {  
      if(which_clus[t] == 1){
        mus = cbind(unlist(fit[,"mu_0"])[t], unlist(fit[,"mu_a"])[t])
        names(mus) = c("mu_0", "mu_a")
        #tmp = unlist(fit[t,grep("clust", unlist(dimnames(fit[[1]])))])
        p = unlist(fit[[1]][,"p"])[t]
        names(p) = "p"
        tau = unlist(fit[[1]][,"tau"])[t]
        names(tau) = "tau"
        c(mus, p, tau)      #tmp, 
      } else {
        mus = cbind(unlist(fit[,"mu_a"])[t], unlist(fit[,"mu_0"])[t])
        names(mus) = c("mu_0", "mu_a")
        #tmp = unlist(fit[t,grep("clust", unlist(dimnames(fit[[1]])))])
        #tmp[which(tmp == 0)] = 1
        #tmp[which(unlist(fit[t,grep("clust", unlist(dimnames(fit[[1]])))]) == 1)] = 0
        p = 1-unlist(fit[,"p"])[t]
        names(p) = "p"
        tau = unlist(fit[,"tau"])[t]
        names(tau) = "tau"
        c(mus, p, tau)      #tmp, 
      }
    })))
  #}
    thetas = Mclust(mu_out[,"mu_a"], G=2, model="E")
    mus_a  = thetas$parameters$mean
    clust = thetas$classification
    
    if(km){
      thetas = kmeans(mu_out[,"mu_a"], 2)
      mus_a = thetas$centers
      clust = thetas$cluster
    }
    all_est = data.frame(est = mu_out[,"mu_a"], clust = clust)
    
    if(order(mus_a)[1] == 1){
      all_est$est[which(all_est$clust == 1)] = all_est$est[which(all_est$clust == 1)]*-1
    } else {
      all_est$est[which(all_est$clust == 2)] = all_est$est[which(all_est$clust == 2)]*-1
    }
    
    mu_out[,"mu_a_comb"] = all_est$est
  return(list(clusters = L_mat, mu_out = mu_out))
}


relabel_stephens = function(dataList, fit){
  require(label.switching)
  Qmat = array(rep(NaN), c(nrow(fit[[1]]), length(dataList$y), 
                           length(grep("mu", unlist(dimnames(fit[[1]]))))    ))
  
  for(t in 1: nrow(fit[[1]])){
    for(n in 1:length(dataList$y)){
      mu_0_use = ifelse(dataList$y[n] > 0, abs(unlist(fit[t,"mu_0"])), -1*abs(unlist(fit[t,"mu_0"])))
      #mu_0_use = unlist(fit[t,"mu_0"])
      pr_nul1 = dnorm(dataList$y[n], mean=mu_0_use, sd=(unlist(fit[t,"tau"])^-0.5)) * 
        (1-unlist(fit[t,"p"]))
      
      mu_a_use = ifelse(dataList$y[n] > 0, abs(unlist(fit[t,"mu_a"])), -1*abs(unlist(fit[t,"mu_a"])))
      #mu_a_use = unlist(fit[t,"mu_a"])
      pr_alt1 = dnorm(dataList$y[n], mean=mu_a_use, sd=(unlist(fit[t,"tau"])^-0.5)) * 
        unlist(fit[t,"p"])
      Qmat[t,n,1] = pr_nul1/(pr_nul1+pr_alt1) 
      Qmat[t,n,2] = pr_alt1/(pr_nul1+pr_alt1) 
      
    }
  }
  
  perm = stephens(Qmat)
  
  mus_orig = cbind(unlist(fit[,"mu_0"]), unlist(fit[,"mu_a"]))
  mus_fixd = t(sapply(1:nrow(perm$permutations), function(i)
    c(mus_orig[i, perm$permutations[i,1]], mus_orig[i, perm$permutations[i,2]])))
  colnames(mus_fixd) = c("mu_0", "mu_a")
  ps_fixd = sapply(1:nrow(perm$permutations), function(i)
    ifelse(perm$permutations[i,1] == 1, unlist(fit[,"p"])[i], 1-unlist(fit[,"p"])[i]))
  clust_fixed = t(sapply(1:nrow(perm$permutations), function(i){
    if(perm$permutations[i,1] == 1){
      tmp = fit[[1]][i,grep("clust", dimnames(fit[[1]])[[2]])]
    } else{
      tmp = fit[[1]][i,grep("clust", dimnames(fit[[1]])[[2]])]
      tmp[which(tmp == 1)] = 0
      tmp[which(fit[[1]][i,grep("clust", dimnames(fit[[1]])[[2]])] == 0)] = 1
    }
    tmp
  }))
  
  new_mcmc = cbind(clust_fixed, mus_fixd, ps_fixd, tau = fit[[1]][,"tau"])
  return(list(new = new_mcmc, perm = perm))
}

