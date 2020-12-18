library(rjags)
library(coda)
library(ggplot2)
library(RColorBrewer)
library(ggmap)
library(data.table)
library(lmerTest)
library(tidyverse)
library(gridExtra)
library(grid)
#setwd("~/matnut/src")
setwd("C:/Users/Kathie/matnut/src")

source("./matnut/lmer_functions_rna.R")
source("./matnut/jags_functions.R")
source("./matnut/boxcox_functions.R")
#source("./matnut/matching_functions3.R")
source("./matnut/matching_functions_lump.R")

source("./matnut/prediction_functions.R")
source("./matnut/summary_functions.R")
###### Read in data
dir <- "C:/Users/Kathie/Dropbox\ (ValdarLab)"
matnut <- readRDS(file.path(dir,'phenotype_analysis/matnut_data.rds'))
matnut_jags = readRDS(file.path(dir, "phenotype_analysis/out/complete_phen_model.rds"))
matnut_stan = readRDS(file.path(dir, "phenotype_analysis/out/stan20PhenOut_15dec2020.rds"))

####################
#   Run function   #
####################

myPhen <- makeSummary(datalist=matnut, phenotype=matnut$ptypes[5:7], 
                        randvar=c("DamID", "RIX", "Diet:RIX"), fixvar="Diet", POvar=c("RIX", "Diet:RIX"),
                        tryLam=c(-1, 0, .25, .33, .5, 1, 2, 3), 
                        chains=2, n.adapt=20000, n.iter=100000)
#n.adapt=20000, n.iter=50000

#saveRDS(myPhen, file.path("./","matnut_outputs/",'allPheno_modelfits.rds'))
testPhen <- readRDS(file.path(dataSource, 'allPheno_modelfits.rds'))

matched_results <- list()
### matching ###
for(phen in c(1,2,4,13,18)) {      # 1:length(matnut$ptypes)
  matnut_use = matnut
  matnut_use$df$DietRIXPOq = paste(matnut_use$df$DietRIX, matnut_use$df$PO,sep="_")
  tab = sort(table((matnut_use$df %>% filter(!is.na(get(matnut$ptypes[phen]))))$DietRIXPOq))
  rem = unique(unlist(strsplit(names(tab)[which(tab < 5)],"_"))[c(T,F)])
  matnut_use$df = matnut_use$df %>% filter(!DietRIX %in% rem)
  matched_results[[matnut$ptypes[phen]]] <- match.multimp(data=matnut_use, ptypes=matnut$ptypes[phen], 
                                                N=5,chains=1, n.iter=10000, 
                                                tryLam=c(-1, 0, .25, .33, .5, 1, 2, 3), thin=10,
                                                randvar=NA, fixvar="Diet",  
                                                matchon=c("Cage","RIX","Diet"), matchoff="PO", idcol="ID",
                                                fixPO = F, fixClust = F, p="p", beta=T, clust=NULL)
  
  #plot_rix = matnut_use$df[-which(is.na( matnut_use$df$ptypes[i]))]
}


#do.call("cbind",lapply(1:3, function(i) 
#  c(sapply(1:9, function(j) 
#   mean(fit_stan[[1]]@sim$samples[[1]][[paste0("p[",j,",",i,"]")]]))) ))

summcmc <- lapply(matnut_stan, function(x) if(!is.na(nrow(x))) summary(As.mcmc.list(x)))

sum_stan = sapply(1:length(summcmc),function(i){
  tmp = NULL
  if(!is.null(names(summcmc[[i]]))){
    tmp = data.frame(do.call("cbind",list(summcmc[[i]][[1]],summcmc[[i]][[2]])))
    tmp$phen = names(summcmc)[[i]]
    tmp$param = rownames(tmp)
    tmp = tmp[-grep("raw",tmp$param),]
    tmp$Level = tmp$Variable = ""
    tmp$Level[grep("SPO", tmp$param)] = encoded$Level[which(encoded$Variable == "PORIX")][c(2:length(grep("SPO", tmp$param)),1)]
    tmp$Level[grep("_s[[1-9]", tmp$param)] = encoded$Level[which(encoded$Variable == "RIX")][c(2:length(grep("_s[[1-9]", tmp$param)),1)]
    tmp$Level[grep("_d[[1-9]", tmp$param)] = encoded$Level[which(encoded$Variable == "Diet")][c(2:length(grep("_d[[1-9]", tmp$param)),1)]
    tmp$Level[grep("_sd[[1-9]", tmp$param)] = encoded$Level[which(encoded$Variable == "DietRIX")][c(2:length(grep("_sd[[1-9]", tmp$param)),1)]
    tmp$Level[grep("_sdp[[1-9]", tmp$param)] = encoded$Level[which(encoded$Variable == "PODietRIX")][c(2:length(grep("_sdp[[1-9]", tmp$param)),1)]
    tmp$Variable[grep("SPO", tmp$param)] = "PORIX"
    tmp$Variable[grep("_s[[1-9]", tmp$param)] = "RIX"
    tmp$Variable[grep("_d[[1-9]", tmp$param)] = "Diet"
    tmp$Variable[grep("_sd[[1-9]", tmp$param)] = "DietRIX"
    tmp$Variable[grep("_sdp[[1-9]", tmp$param)] = "PODietRIX"
  }
  tmp
})
names(sum_stan) = names(summcmc)
plot_tmp = sum_stan$WeightPND21
plot_tmp = plot_tmp %>% filter(Level!= "")

 <- plot.inter.ci(med=plot_tmp$X50., mu=plot_tmp$Mean, 
                          hpd.narrow = cbind(plot_tmp$X25., plot_tmp$X75.), 
                          hpd.wide = cbind(plot_tmp$X2.5., plot_tmp$X97.5.), 
                          names=plot_tmp$Level, order=3,
                          col.midvals="white", pch.midvals="|", addline = F, wide=T, 
                          grouped=plot_tmp$Variable, ordered=F)
comPlot <- jagsPlot$plot + geom_point(data=lmerSumTemp, aes(x=Level, y=Intercept), col="red", size=3) +
  ggtitle(paste("LMER and JAGS comparison for", phenotype))




param_df = do.call("rbind",param_df)
rownames(param_df) = NULL
allKeep[[i]] = list(summary = summcmc, sig_params = param_df)
#mins = unlist(lapply(dat, function(x) summary(as.mcmc.list(x$multi

mcmc_compare_plot_list = sapply(1:length(matched_results), function(i) 
  compare_muA_stan_jags(results_object = matched_results[[i]], phen=names(matched_results)[[i]]),
  simplify=F)
names(mcmc_compare_plot_list) = names(matched_results)

p_plot_list = lapply(matched_results, function(x) plot_p(x$fit_stan))

matched_df_ob_list = lapply(matched_results, function(x) {
  tmp = x$matched_df
  for(i in 1:length(tmp)){
    tmp[[i]]$chain = i
  }
  do.call("rbind", tmp)
})

plot_raw_delts_list = sapply(1:length(matched_results), function(i)
  plot_deltas(matched_df_ob_list[[i]], p_plot_list[[i]]$cluster))
names(plot_raw_delts_list) = names(matched_results)


i=1
grid.arrange(mcmc_compare_plot_list[[i]]$plot)
p_plot_list[[i]]$plot
grid.arrange(plot_raw_delts_list[[i]])


compare_muA_stan_jags = function(results_object, phen = NULL){
  stan_lst = lapply(results_object$fit_stan, function(x) as.array(x)[,,"mu_abs"])

  stan_it = length(stan_lst[[1]])
  jags_it = length(results_object$multimp$fit$mu_a[[1]])
  stan_n = length(stan_lst)
  
  jags_n = length(results_object$multimp$fit$mu_a_abs)
  
  stan_array = array(NA, c(stan_it,stan_n,1))
  jags_array = array(NA, c(jags_it,jags_n,1))
  stan_df = matrix(NA, nrow=stan_it*length(stan_lst), ncol=2)
  jags_df = matrix(NA, nrow=jags_it*length(results_object$multimp$fit$mu_a), ncol=2)
  colnames(stan_df) = colnames(jags_df) = c("est","chain")
  stan_df = data.frame(stan_df)
  jags_df = data.frame(jags_df)
  
  np_cp = c()
  for(i in 1:length(results_object$fit_stan)){
    if(i == 1){
      np_cp = nuts_params(results_object$fit_stan[[1]])
    } else {
      np_tmp = nuts_params(results_object$fit_stan[[i]])
      np_tmp$Chain = i
      np_cp = rbind(np_cp, np_tmp)
    }
    tmp = stan_lst[[i]]
    if(mean(tmp) < 0) tmp = -tmp
    stan_array[,i,1] = tmp
    stan_df$est[(((i-1)*stan_it)+1):(i*stan_it)] = tmp
    stan_df$chain[(((i-1)*stan_it)+1):(i*stan_it)] = i
    
    tmp = results_object$multimp$fit$mu_a[[i]]
    if(mean(tmp) < 0) tmp = -tmp
    jags_array[,i,1] = tmp
    jags_df$est[(((i-1)*jags_it)+1):(i*jags_it)] = tmp
    jags_df$chain[(((i-1)*jags_it)+1):(i*jags_it)] = i
  }
  stan_df$chain = as.factor(stan_df$chain)
  jags_df$chain = as.factor(jags_df$chain)
  
  dimnames(stan_array) = dimnames(jags_array) = list(NULL, c(1:5),"mu_a")
  ms = mcmc_trace(stan_array, pars="mu_a", np=np_cp)
  ms = ms + theme_classic() + scale_color_brewer(palette="Spectral") + 
    ggtitle(paste("Traceplot stan",phen))
  mj = mcmc_trace(jags_array, pars="mu_a")
  mj = mj + theme_classic() + scale_color_brewer(palette="Spectral") + 
    ggtitle(paste("Traceplot jags",phen))
  ds = ggplot(stan_df, aes(x=est, group=chain, col=chain)) + 
    geom_density() + 
    theme_classic() + 
    scale_colour_brewer(palette="Spectral") + 
    ggtitle(paste("Density curves stan", phen))
  dj = ggplot(jags_df, aes(x=est, group=chain, col=chain)) + 
    geom_density() + 
    theme_classic() + 
    scale_colour_brewer(palette="Spectral") + 
    ggtitle(paste("Density curves jags", phen))
  plot_sum = arrangeGrob(ms,mj,ds,dj,nrow=2,ncol=2)
  
  jags_sum = jags_df %>% group_by(chain) %>%
    summarize(mean=mean(est),
              med = median(est),
              lo = quantile(est, probs=0.025), 
              hi = quantile(est, probs=0.975),
              spread = hi-lo,
              est = "jags")
  stan_sum = stan_df %>% group_by(chain) %>%
    summarize(mean=mean(est),
              med = median(est),
              lo = quantile(est, probs=0.025), 
              hi = quantile(est, probs=0.975),
              spread = hi-lo,
              est = "stan")
  sum = rbind(jags_sum, stan_sum)
  return(list(plot=plot_sum, df=sum))
}     


plot_deltas = function(matched_df_ob, clusters){
  n_chains = length(unique(matched_df_ob$chain))
  matched_df_ob = matched_df_ob %>%
    mutate(RIX = factor(RIX, levels=c(1:4,6:10)),
           chain_RIX = factor(paste(chain,RIX,sep="_"), 
                       levels=apply(expand.grid(unique(matched_df_ob$chain), c(1:4,6:10)), 1, function(x)
                         paste(x, collapse="_")))) %>%
    arrange(chain, RIX, Cage)
  
  p1 = ggplot(data = matched_df_ob, aes(group=chain_RIX)) + 
    theme_classic() + 
    geom_density(aes(x=y, color = RIX)) +
    scale_color_manual(values = hue_pal()(9))
  
  p_cols = rep(ifelse(is.na(clusters$clus), "#808080", 
                      ifelse(clusters$clus == 0, "#26b9ef",  "#8d00b0")) ,each=n_chains)
  matched_df_ob$clus = clusters$clus[match(matched_df_ob$RIX,clusters$rix)]
  matched_df_ob$clus[which(is.na(matched_df_ob$clus))] = 2
  matched_df_ob$clus = factor(matched_df_ob$clus, levels=c(0,1,2))
  
  p2 = ggplot(data = matched_df_ob, aes(group=chain_RIX)) + 
    theme_classic() + 
    geom_density(aes(x=y, color = clus)) +
    scale_color_manual(values=c(hue_pal()(2),"#808080"),
                       labels=c("mu_0","mu_a","unclear"))
  
  p <- arrangeGrob(p1, p2, name = paste("Raw deltas over",n_chains,"chains"))
  return(p)
}
matched_results$WeightPND21$matched_df[[1]]

plot_p = function(stan_list){
  n_chains = length(stan_list)
  p_results = p_sum_list = list()
  for(c in 1:n_chains){
    summary = data.frame(summary(stan_list[[c]])[[1]])
    p_sum = summary[grep("p[[]", rownames(summary)),]
    p_results[[c]] = data.frame(matrix(NA, nrow=(nrow(p_sum)*6), ncol=4))
    rix_mat = data.frame(rix = c(1:4,6:10), num = 1:9)
    nEst = length(c("mean","X2.5.","X25.","X50.","X75.","X97.5."))
    for(i in 1:nrow(p_sum)){
      ind = (i-1)*nEst + 1
      estCol = which(colnames(p_sum) %in% c("mean","X2.5.","X25.","X50.","X75.","X97.5."))
      colnames(p_results[[c]]) = c("rix","mu","est","type")
      rix_tmp = unlist(strsplit(rownames(p_sum)[i],"[[]|,|[]]"))[2] 
      
      p_results[[c]]$rix[ind:(ind+nEst-1)] = rix_mat$rix[match(rix_tmp, rix_mat$num)]
      p_results[[c]]$mu[ind:(ind+nEst-1)] = unlist(strsplit(rownames(p_sum)[i],"[[]|,|[]]"))[3] 
      p_results[[c]]$est[ind:(ind+nEst-1)] = unlist(p_sum[i,estCol])
      p_results[[c]]$type[ind:(ind+nEst-1)] = c("mean","X2.5.","X25.","X50.","X75.","X97.5.")
    }
    p_results[[c]]$chain = c
    rix_tmp = unlist(do.call("rbind",strsplit(rownames(p_sum),"[[]|,|[]]"))[,2])
    p_sum$rix = rix_mat$rix[match(rix_tmp, rix_mat$num)]
    p_sum$mu = as.numeric(unlist(do.call("rbind",strsplit(rownames(p_sum),"[[]|,|[]]"))[,3]))
    colnames(p_sum)[grep("X", colnames(p_sum))] = c("lolo", "lo", "med", "hi", "hihi")
    p_sum_list[[c]] = p_sum
    p_sum_list[[c]]$chain = c
    
  }
  p_results = do.call("rbind", p_results)
  p_sum_list = do.call("rbind", p_sum_list)
  p_sum_list$chain = as.factor(p_sum_list$chain)
  p_sum_list %>% group_by(rix, chain) %>% 
    arrange(rix, chain, desc(mean)) %>%
    slice(1) %>%
    ungroup() %>% group_by(rix, mu) %>%
    tally() -> clusters_simp
  
  certain_cluster = (clusters_simp %>% tally() %>%
    filter(n==1))$rix
  clusters = data.frame(matrix(NA, nrow=length(unique(p_sum_list$rix)), ncol=3))
  colnames(clusters) = c("rix", "mu", "clus")
  clusters$rix = sort(unique(p_sum_list$rix))
  clusters$mu[which(clusters$rix %in% certain_cluster)] = 
    clusters_simp$mu[match(certain_cluster,clusters_simp$rix)]
  clusters$clus = ifelse(is.na(clusters$mu), NA, ifelse(clusters$mu == 2, 0, 1))
  clusters$rix = factor(clusters$rix, levels=rix_mat$rix)
  #p_results %>% group_by(rix, mu, type) %>% 
  #  summarize(est = mean(est)) %>% 
  #  filter(type == "mean")
    
  p = ggplot(data=p_sum_list) + 
    facet_grid(~ rix) + 
    theme_classic() + 
    scale_color_viridis_d() + 
    scale_alpha(range=c(0.1,0.2)) + 
    geom_ribbon(aes(x=mu, ymin=lolo, ymax=hihi, alpha=0.1, group=chain)) + 
    geom_ribbon(aes(x=mu, ymin=lo, ymax=hi, alpha=0.2,group=chain)) + 
    geom_line(aes(x=mu, y=med, color=chain)) + 
    geom_point(aes(x=mu, y=mean, color=chain)) + 
    scale_x_continuous(breaks=c(1,2,3)) 
  return(list(plot = p, cluster = clusters))
}


saveRDS(matched_results, file.path(dir,"phenotype_analysis/out/jags_stan_phen_10k_20oct2020.rds"))
#multimp = readRDS(file.path(dir,"phenotype_analysis/out/mult_imp_match_7jul2020.rds"))


plot_clust(matched_results$WeightPND21$multimp$fit$`mu_0`, 
           matched_results$WeightPND21$multimp$fit$`mu_a`)

get_clusters(matched_results$WeightPND60$multimp$fit)
get_clusters(matched_results$OFTotalDistance$multimp$fit)

x=matched_results$WeightPND21
plots_phen = lapply(matched_results, function(x){  
  clusters = get_clusters(x$multimp$fit)
  
  #if(length(unique(clusters)) == 1){
  #  clus = plot_clust(x$multimp$fit[["mu_a"]],
  #                    x$multimp$fit[["mu_a_sq"]])
    clus = plot_clust_single(x$multimp$fit[["mu_a_sq"]])
  #} else {
  #  clus = plot_clust(x$multimp$fit[["mu_0"]],
  #                    lapply(x$multimp$fit[["mu_a"]], abs))
  #}
  #plot(mcmc.list(lapply(x$multimp$fit$mu_a, function(x) (x^2)^0.5)))
  #plot(mcmc.list(lapply(x$multimp$fit$mu_a, function(x) abs(x))))
  
  #clusters = get_clusters(x$multimp$fit)$medians  
      ### only picks up one side
  rix = plot_rix(x$matched_df, clusters$medians)
  return(list(clus, rix))
})


require(gridExtra)
windows()
grid.arrange(plots_phen[[2]][[2]][[1]],plots_phen[[2]][[2]][[2]], top="Raw data", ncol = 2, nrow = 1)

grid.arrange(plots_phen[[1]][[2]][[1]],plots_phen[[1]][[2]][[2]], top="Raw data", ncol = 2, nrow = 1)
grid.arrange(plots[[2]][[2]][[1]],plots[[2]][[2]][[2]], top="Raw data", ncol = 2, nrow = 1)


pdf("C:/Users/Kathie/Dropbox (ValdarLab)/phenotype_analysis/plots_8sep2020.pdf")
for(n in 1:length(plots)){
  grid.arrange(plots[[n]][[2]][[1]],plots[[n]][[2]][[2]], top=names(plots)[n], ncol = 2, nrow = 1)
  print(plots[[n]][[1]])
  plot(matched_results[[n]]$multimp$fit$mu_a_sq[[1]], main=names(plots)[n])
  plot(matched_results[[n]]$multimp$fit$mu_a, main = names(plots)[n])
}
dev.off()

### c(1,2,4,5,13,17,18)
muClust1 = data.frame(do.call("cbind", matched_results[[n]]$multimp$fit$mu_0))
colnames(muClust1) = paste0("it",seq(1:ncol(muClust1)))
muClust1_long = gather(muClust1, "it","est")
ggplot(muClust1_long, aes(x=est, col=it)) + geom_density() + xlim(-5,5)

muClust2 = data.frame(do.call("cbind", matched_results[[n]]$multimp$fit$mu_a))
colnames(muClust2) = paste0("it",seq(1:ncol(muClust2)))
muClust2_long = gather(muClust2, "it","est")
ggplot(muClust2_long, aes(x=est, col=it)) + geom_density()
#plot(full_list[[10]])


plots$WeightPND60[[1]]
plots$OFTotalDistance[[1]]

##############################################################
plot_clust = function(muList1, muList2){
  muClust1 = data.frame(do.call("cbind", muList1))
  muClust2 = data.frame(do.call("cbind", muList2))
  meanClust1 = colMeans(muClust1)
  meanVarClust1 = mean(apply(muClust1, 2, var))
  Bm1 = sum(unlist(lapply(meanClust1, function(x) 
    ((x - mean(meanClust1))^2)/length(meanClust1))))
  Tm1 = meanVarClust1 + (1+1/length(meanClust1))*Bm1
  
  
  meanClust2 = colMeans(muClust2)
  meanVarClust2 = mean(apply(muClust2, 2, var))
  Bm2 = sum(unlist(lapply(meanClust2, function(x) 
    ((x - mean(meanClust2))^2)/length(meanClust2))))
  Tm2 = meanVarClust2 + (1+1/length(meanClust2))*Bm2
  
  
  min = max(-5, min(min(HPDinterval(as.mcmc(muClust1), prob = 0.999)[,1]),
                    min(HPDinterval(as.mcmc(muClust2), prob = 0.999)[,1])))
  max = min(8, max(max(HPDinterval(as.mcmc(muClust1), prob = 0.999)[,2]),
                   max(HPDinterval(as.mcmc(muClust2), prob = 0.999)[,2])))
  ribs_min1 = HPDinterval(as.mcmc(rowMeans(muClust1)), prob = 0.95)
  mean_est1 = mean(rowMeans(muClust1))
  ribs_min2 = HPDinterval(as.mcmc(rowMeans(muClust2)), prob = 0.95)
  mean_est2 = mean(rowMeans(muClust2))
  

  colnames(muClust1) = paste0("it",seq(1:ncol(muClust1)))
  muClust_long1 = gather(muClust1, "it","est")
  muClust_long1$it = factor(muClust_long1$it, levels=colnames(muClust1))
  muClust_long1$mean = F
  muClust_long1$clust = 1
  
  muClust_long1 = rbind(muClust_long1, 
                       data.frame(it = "mean", mean=T, clust = 1, 
                                  est = rowMeans(muClust1)))
  
  colnames(muClust2) = paste0("it",seq(1:ncol(muClust2)))
  muClust_long2 = gather(muClust2, "it","est")
  muClust_long2$it = factor(muClust_long2$it, levels=colnames(muClust2))
  muClust_long2$mean = F
  muClust_long2$clust = 2
  
  muClust_long2 = rbind(muClust_long2, 
                        data.frame(it = "mean", mean=T, clust = 2, 
                                   est = rowMeans(muClust2)))
  mcClust = rbind(muClust_long1, muClust_long2)
  mcClust$clustcol = paste0(mcClust$mean, mcClust$clust)
  mcClust$clustcol = factor(mcClust$clustcol, levels=c("FALSE1", "TRUE1", "FALSE2", "TRUE2"))
  coluse = c("#07575B","#C4DFE6","#4B7447","#A2C523")
  p = ggplot(mcClust %>% filter(clustcol == "TRUE1"), aes(x=est, group=it)) +
    geom_line(stat='density', col=coluse[1], size=1) + 
    geom_vline(xintercept=mean_est1, linetype="dashed", col=coluse[1]) + 
    annotate("rect", xmin=ribs_min1[1], xmax=ribs_min1[2], 
             ymin=-Inf, ymax=Inf, alpha=0.1, fill=coluse[1])
  p = p + 
    geom_line(stat='density', data=mcClust %>% filter(clustcol == "FALSE1"),
                 aes(group=it), col=coluse[2])
  p = p + 
    geom_line(stat='density', data=mcClust %>% filter(clustcol == "TRUE2"),
                aes(group=it), col=coluse[3], size=1) + 
    geom_vline(xintercept=mean_est2, linetype="dashed", col=coluse[3]) + 
    annotate("rect", xmin=ribs_min2[1], xmax=ribs_min2[2], 
             ymin=-Inf, ymax=Inf, alpha=0.1, fill=coluse[3])
  p = p + 
    geom_line(stat='density', data=mcClust %>% filter(clustcol == "FALSE2"),
                 aes(group=it), col=coluse[4]) +
    theme_bw() + 
    xlim(min, max) + aes(ymin=0)
  return(p)
}

plot_clust_single = function(muList1){
  muClust1 = data.frame(do.call("cbind", muList1))
  meanClust1 = colMeans(muClust1)
  meanVarClust1 = mean(apply(muClust1, 2, var))
  Bm1 = sum(unlist(lapply(meanClust1, function(x) 
    ((x - mean(meanClust1))^2)/length(meanClust1))))
  Tm1 = meanVarClust1 + (1+1/length(meanClust1))*Bm1
  
  min = 0  #max(-5, HPDinterval(as.mcmc(muClust1), prob = 0.999)[,1]*0.9)
  max = min(5, max(HPDinterval(as.mcmc(muClust1), prob = 0.975)[,2])*1.5)
  ribs_min1 = HPDinterval(as.mcmc(rowMeans(muClust1)), prob = 0.95)
  mean_est1 = mean(rowMeans(muClust1))
  
  colnames(muClust1) = paste0("it",seq(1:ncol(muClust1)))
  muClust_long1 = gather(muClust1, "it","est")
  muClust_long1$it = factor(muClust_long1$it, levels=colnames(muClust1))
  muClust_long1$mean = F

  muClust_long1 = rbind(muClust_long1, 
                        data.frame(it = "mean", mean=T,est = rowMeans(muClust1)))
  mcClust = muClust_long1
  coluse = c("#07575B","#C4DFE6","#4B7447","#A2C523")
  p = ggplot(mcClust %>% filter(mean), aes(x=est, group=it)) +
    geom_line(stat='density', col=coluse[1], size=1) + 
    geom_vline(xintercept=mean_est1, linetype="dashed", col=coluse[1]) + 
    annotate("rect", xmin=ribs_min1[1], xmax=ribs_min1[2], 
             ymin=-Inf, ymax=Inf, alpha=0.1, fill=coluse[1])
  p = p + 
    geom_line(stat='density', data=mcClust %>% filter(!mean),
              aes(group=it), col=coluse[2]) + 
    theme_bw() + 
    xlim(min, max) + aes(ymin=0)
  return(p)
}


plot_rix <- function(dfList, clusters){
  dfList = sapply(1:length(dfList), function(x) {
    dfList[[x]]$n = x
    tmp = dfList[[x]] %>% group_by(RIX) %>% tally()
    tmp$prop = tmp$n/max(tmp$n)
    dfList[[x]]$alph = tmp$prop[match(dfList[[x]]$RIX,tmp$RIX)]
    as.data.frame(dfList[[x]])
  }, simplify=F)
  
  
  fitDf = do.call("rbind",dfList)
  fitDf$RIX = factor(fitDf$RIX, ordered=T)
  fitDf = fitDf %>% arrange(n, RIX)
  
  max = max(fitDf$norm_y)*1.1
  min = min(fitDf$norm_y)*1.1
  coluse = c("#999999", "#E69F00", "#56B4E9", "#009E73",
             "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#b4eeb4")
  coluse = c('#e41a1c','#377eb8','#4daf4a','#984ea3',
             '#ff7f00','#ffff33','#a65628','#f781bf','#999999')
  names(coluse) = levels(fitDf$RIX)
  fitDf$clust = as.factor(clusters[match(as.numeric(fitDf$RIX), 1:length(clusters))])
  alphaMap = fitDf %>% select(RIX, alph) %>% distinct()

  p = ggplot(fitDf, aes(x=norm_y,group=n))
  for(i in levels(fitDf$RIX)){
    p = p + geom_line(stat='density', data=fitDf %>% filter(RIX == i),
                      alpha=alphaMap$alph[which(alphaMap$RIX == i)], aes(color=clust))
  }
  p = p + xlim(min, max) + aes(ymin=0) + theme_classic() +
    labs(x = "Delta",
         y = "Density",
         color = "Cluster #") 
  
  p2 = ggplot(fitDf, aes(x=norm_y,group=n))  
  for(i in levels(fitDf$RIX)){
    p2 = p2 + geom_line(stat='density', data=fitDf %>% filter(RIX == i),
                        alpha=alphaMap$alph[which(alphaMap$RIX == i)], aes(color=RIX))
  }
  p2 = p2 + xlim(min, max) + aes(ymin=0) + theme_classic() +
    labs(x = "Delta",
         y = "Density",
         color = "RIX") +
    scale_color_manual(values = coluse)
    #scale_color_manual(values = rainbow(9), labels=levels(fitDf$RIX))

  return(list(p_clus = p, p_rix = p2))
}

get_clusters = function(mcList){
  clustLists = grep("clust", names(mcList))
  lists = lapply(clustLists, function(x) 
    data.frame(do.call("cbind", mcList[[x]]))  )
  names(lists) = names(mcList)[clustLists]
  means = unlist(lapply(lists, function(x) 
    summary(as.mcmc(colMeans(x)))[["statistics"]][["Mean"]]))
  medians = unlist(lapply(lists, function(x) 
    median(summary(as.mcmc(x))[["quantiles"]][,"50%"])))
  return(list(means=means, medians=medians))
}

get_clusters_2 = function(clus_mat){
  get_min_likelihood = lapply(clus_mat, function(x){
    t(apply(x,2,function(y){
      c(alt = min(y[1], y[3]), nul = y[2])
    }))
  })
  
  sum_cluster = t(sapply(1:nrow(get_min_likelihood[[1]]), function(r){
    alt = mean(unlist(lapply(get_min_likelihood, function(x) x[r,1])))
    nul = mean(unlist(lapply(get_min_likelihood, function(x) x[r,2])))
    c(alt=alt, nul=nul)
  }))
  
  clusters = apply(sum_cluster, 1, which.min)
  clusters[which(clusters == 2)] = 0
  return(clusters)
}

###########################################################################

graphics.off()
pdf(file.path(dir,"phenotype_analysis/out/mult_imp_match_catplots.pdf"), 
    onefile=TRUE, width = 8.5, height = 11)
for(i in 1:length(multimp)){  
  print(multimp[[i]]$plot$plot)
}

dev.off()
pdf(file.path(dir,"phenotype_analysis/out/mult_imp_match_predribplots.pdf"), 
    onefile=TRUE, width = 8.5, height = 11)
for(i in 1:length(multimp)){  
  print(multimp[[i]]$ribPlot$ribbon[[1]])
}
dev.off()

multimp = readRDS(file.path(dir, "phenotype_analysis/mult_imp_match.rds"))

tables <- writeTables(results = testPhen, df=matnut)



###############
# Make tables #
###############

dietDiet_table <- list()
jagsLmer_compare <- list()
rixRix_table <- list()
summary_table <- list()

psychPhen <- testPhen

for(i in 1:length(ptypes)){
  pheno <- ptypes[i]
  
  ####
  tempTable <- jags.getDecodedSummary(psychPhen[[pheno]]$mcmcObject, encoded, narrow=0.5, wide=0.95)
  tempTable$Level <- factor(tempTable$Level, levels=tempTable$Level)
  null_p <- data.frame(Level = names(psychPhen[[pheno]]$null.test$pval), 
                       pval = psychPhen[[pheno]]$null.test$pval, 
                       signif = psychPhen[[pheno]]$null.test$signif)
  null_p$Level <- factor(null_p$Level, levels=tempTable$Level)
  null_p <- null_p[order(null_p$Level),]
  mergeTable <- merge(tempTable, null_p, by="Level")
  mergeTable$Level <- factor(mergeTable$Level, levels=encoded$Level)
  mergeTable$Variable <- factor(mergeTable$Variable, levels=unique(encoded$Variable))
  shortTable <- mergeTable[which(mergeTable$Variable %in% c("Diet","RIX","DietRIX","PORIX","PODietRIX")),]
  
  summary_table[[pheno]] <- cbind(phenotype=rep(pheno,nrow(shortTable)), 
                                  shortTable[order(shortTable$Variable, shortTable$Level),])
  ####
  #dietTab <- contrasts.getDecodedSummary(psychPhen$OFTotalDistance, "Diet")
  #dietDiet_table[[pheno]] <- cbind(phenotype=rep(pheno,nrow(dietTab)), dietTab)
  #rixTab <- contrasts.getDecodedSummary(psychPhen$OFTotalDistance, "RIX")
  #rixRix_table[[pheno]] <- cbind(phenotype=rep(pheno,nrow(rixTab)), rixTab)
}


write.csv(do.call(rbind, summary_table), file.path("./","matnut_outputs/",'allPheno_summarytable.csv'), row.names = F)




###########
#  Plots  #
###########


### Flag Plots ###

flagPlots <- function(results, df, ptypes=NULL, encoded=NULL, write=F, byVar = c("Diet","RIX")){
  if(is.null(byVar)) byVar = "Diet"
  contrastPlots <- list()
  returnPval <- list()
  if (class(df) == "list"){
    ptypes <- df$ptypes
    encoded <- df$encoded
  }
  
  for(i in 1:length(ptypes)){  
    useTest <- ifelse(byVar=="Diet","diet.test","rix.test")
    bindTemp <- results[[ptypes[i]]][[useTest]]$dir
    rownames(bindTemp) <- rownames(results[[ptypes[i]]][[useTest]]$pval)
    colnames(bindTemp) <- paste0("a",colnames(results[[ptypes[i]]][[useTest]]$pval))
    meltTemp <- melt(bindTemp)
    colnames(meltTemp)[3] <- "Dir"
    tempPval <- results[[ptypes[i]]][[useTest]]$pval
    colnames(tempPval) <- paste0("a",colnames(results[[ptypes[i]]][[useTest]]$pval))
    
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
    contrastPlots[[i]] <- getPlot.compare.helper(dfplot,ptypes=ptypes,byVar=byVar)
  }
  return(contrastPlots)
}



### Mega-plot ###

keep <- c("Diet", "DietRIX","PODietRIX","PORIX", "RIX")
myMCMC <- testPhen$WeightPND60$mcmcObject[, which(colnames(testPhen$WeightPND60$mcmcObject) %in% 
                                                    encoded$Level[which(encoded$Variable %in% keep)])]
ribbon_plot(myMCMC, ptypes, encoded)


#ptypes <- names(psychPhenResults)
temp_rib<- list()
for(i in 1:5){#length(ptypes)){
  temp_rib[[ptypes[i]]]$predict <- prediction(mcmcOb=testPhen[[ptypes[i]]]$mcmcObject, ptypes=ptypes[i], encoded, Match=F)
  temp_rib[[ptypes[i]]]$rib <- ribbon_plot(mcmcOb=temp_rib[[ptypes[i]]]$predict, ptypes=ptypes[i], encoded)
}

all_phenDF <- data.frame()

for(i in 1:length(ptypes)){
  tempDF <- data.frame(phen = rep(names(temp_rib)[i], nrow(temp_rib[[names(temp_rib)[i]]]$savedata[[names(temp_rib)[i]]])), 
                       temp_rib[[names(temp_rib)[i]]]$savedata[[names(temp_rib)[i]]])
  all_phenDF <- rbind(all_phenDF, tempDF)
}

allPlot <- ggplot(all_phenDF, aes(color=rix)) + 
  geom_line(aes(y=mu, x=diet, group=poe, alpha=poe)) +  
  geom_ribbon(color=NA, aes(x=diet, ymin=lower.1, ymax=upper.1, group=poe, alpha=poe, 
                            fill=rix, color=rix)) + 
  geom_hline(yintercept = 0, colour="lightsteelblue", size=0.5, linetype="longdash") +
  scale_alpha_discrete(range = rev(c(0.45, 0.75)), name="POE") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_fill_discrete(guide=FALSE) + scale_color_discrete(guide=FALSE) + 
  facet_grid(rix~phen) + theme_minimal() + 
  ylab("Predicted effects") + xlab("Diet") #+ ggtitle(ptypes[j])


### More plots ###

graphics.off()
pdf(file.path("./","matnut_outputs/","allPheno_DietContrastPlots.pdf"), 
    onefile=TRUE, width = 8.5, height = 11)
for(i in 1:length(ptypes)){  
  print(contrastPlots[[i]])
}
graphics.off()

graphics.off()
pdf(file.path("./","matnut_outputs/","allPheno_JLComparison_1jun.pdf"), 
    onefile=TRUE, width = 8.5, height = 11)
for(i in 1:length(ptypes)){  
  print(myPhen[[ptypes[i]]]$plot)
}
graphics.off() 

graphics.off()
pdf(file.path("./","matnut_outputs/","matchedS2_12jul2017.pdf"), onefile=TRUE, width = 14, height = 11)
for(i in 1:length(matnut$ptypes)){  
  print(multimp2[[matnut$ptypes[i]]]$ribPlot$ribbon)
  #print(ribbon4[[ptypes[i]]])
}
graphics.off()

graphics.off()
pdf(file.path("./","matnut_outputs/","allPheno_DietOnlyribbon.pdf"), onefile=TRUE, width = 14, height = 11)
for(i in 1:length(ptypes)){  
  print(ribbon9_Dietonly[[ptypes[i]]]$ribbon9)
}
graphics.off()

### Matched ribbon ###
graphics.off()
pdf(file.path("./","matnut_outputs/","matchedribbon_MIMP.pdf"), onefile=TRUE, width = 14, height = 11)
for(i in 1:length(ptypes)){
  print(multimp[[ptypes[i]]]$plot$ribbon)
}
graphics.off()


###############
# Make tables #
###############

dietDiet_table <- list()
jagsLmer_compare <- list()
rixRix_table <- list()
summary_table <- list()

psychPhen <- testPhen

for(i in 1:length(ptypes)){
  pheno <- ptypes[i]
  
  ####
  tempTable <- jags.getDecodedSummary(psychPhen[[pheno]]$mcmcObject, encoded, narrow=0.5, wide=0.95)
  tempTable$Level <- factor(tempTable$Level, levels=tempTable$Level)
  null_p <- data.frame(Level = names(psychPhen[[pheno]]$null.test$pval), 
                       pval = psychPhen[[pheno]]$null.test$pval, 
                       signif = psychPhen[[pheno]]$null.test$signif)
  null_p$Level <- factor(null_p$Level, levels=tempTable$Level)
  null_p <- null_p[order(null_p$Level),]
  mergeTable <- merge(tempTable, null_p, by="Level")
  mergeTable$Level <- factor(mergeTable$Level, levels=encoded$Level)
  mergeTable$Variable <- factor(mergeTable$Variable, levels=unique(encoded$Variable))
  shortTable <- mergeTable[which(mergeTable$Variable %in% c("Diet","RIX","DietRIX","PORIX","PODietRIX")),]
  
  summary_table[[pheno]] <- cbind(phenotype=rep(pheno,nrow(shortTable)), 
                                  shortTable[order(shortTable$Variable, shortTable$Level),])
  ####
  #dietTab <- contrasts.getDecodedSummary(psychPhen$OFTotalDistance, "Diet")
  #dietDiet_table[[pheno]] <- cbind(phenotype=rep(pheno,nrow(dietTab)), dietTab)
  #rixTab <- contrasts.getDecodedSummary(psychPhen$OFTotalDistance, "RIX")
  #rixRix_table[[pheno]] <- cbind(phenotype=rep(pheno,nrow(rixTab)), rixTab)
}


write.csv(do.call(rbind, summary_table), file.path("./","matnut_outputs/",'allPheno_summarytable.csv'), row.names = F)


#########
# Fixed #
#########

storeAn <- list()
#p_vals <- matrix(NA, nrow=length(ptypes), ncol=5, dimnames=list(ptypes))
#colnames(p_vals) <- c("RIX","Diet","Diet:RIX", "PO:RIX", "PO:DietRIX")

rix <- encoded$Level[which(encoded$Variable == "RIX")]
porix <- paste0("RIX",rix,":PO")
podietrix <- paste0("PO:DietRIX",unlist(strsplit(as.character(encoded$Level[which(encoded$Variable == "PODietRIX")]),'PO'))[c(F,T)])
p_vals <- matrix(NA, nrow=length(ptypes), ncol=(length(podietrix) + length(porix)), dimnames=list(ptypes))
colnames(p_vals) <- c(porix, podietrix)

for(i in 1:length(ptypes)){
  pheno <- paste(ptypes[i])
  use <- which(is.na(orderdat[,pheno]) == F)
  matnut_use <- orderdat[use,]
  
  OF <- ifelse(length(grep("OF",pheno))==0,F,T)
  LD <- ifelse(length(grep("LD",pheno))==0,F,T)
  SIH <- ifelse(length(grep("SIH",pheno))==0,F,T)
  FST <- ifelse(length(grep("FST",pheno))==0,F,T)
  Stress <- ifelse(length(grep("CORT",pheno))==0,F,T)
  
  if(OF){
    indvariable <- "~ -1 + (1 | DamID) + (1 | SireID) + (1 | BehaviorBatch) + OFBox + RIX + Diet + Diet:RIX + (0+PO|RIX) + (0+PO|Diet:RIX)" 
  } else if(LD){
    indvariable <- "~ -1 + (1 | DamID) + (1 | SireID) + (1 | BehaviorBatch) + LDChamber + RIX + Diet + Diet:RIX + (0+PO|RIX) + (0+PO|Diet:RIX)" 
  } else if(SIH){
    indvariable <- "~ -1 + (1 | DamID) + (1 | SireID) + (1 | BehaviorBatch) + SIHOrder + RIX + Diet + Diet:RIX + (0+PO|RIX) + (0+PO|Diet:RIX)" 
  } else if(FST){
    indvariable <- "~ -1 + (1 | DamID) + (1 | SireID) + (1 | BehaviorBatch) + FSTChamber + RIX + Diet + Diet:RIX + (0+PO|RIX) + (0+PO|Diet:RIX)" 
  } else if(Stress){
    indvariable <- "~ -1 + (1 | DamID) + (1 | SireID) + (1 | BehaviorBatch) + RestraintExperimenter + RestraintOrder + 
    RIX + Diet + Diet:RIX + (0+PO|RIX) + (0+PO|Diet:RIX)" 
  } else {
    indvariable <- "~ -1 + (1 | DamID) + (1 | SireID) + (1 | BehaviorBatch) + RIX + Diet + Diet:RIX + (0+PO|RIX) + (0+PO|Diet:RIX)" 
  }
  
  form <- formula(paste0(pheno, indvariable))
  lmObj <- lmer(form, data=matnut_use)
  ano <- anova(lmObj, type=1)
  poef <- fixef(lmObj)
  for(j in 1:ncol(p_vals)){
    val <- ifelse(is.element(colnames(p_vals)[j], names(poef)), poef[which(names(poef) %in% colnames(p_vals)[j])], NA)
    p_vals[i,j] <- val
  }
  #p_vals[i,] <- ano$'Pr(>F)'[c((length(ano$'Pr(>F)')-4):length(ano$'Pr(>F)'))]
  
  storeAn[[i]] <- lmObj
}

################
# Prepare data #
################
matnut <- read.csv("../data/AllPhenotypes_Matnut5.csv")

matnut[grep("a", matnut$Reciprocal),"PO"] <- 0.5
matnut[grep("b", matnut$Reciprocal),"PO"] <- -0.5
matnut[which(matnut$PO == 0.5), "pof"] <- "+"
matnut[which(matnut$PO == -0.5), "pof"] <- "-"
dietlabs <- c("STD", "ME", "PD", "VDD")
matnut$RIX <- factor(matnut$RIX, levels=c(1:4,6:10))
matnut$Diet <- factor(matnut$Diet, levels=dietlabs)


### covariates ###

matnut$BBoriginal <- matnut$BreedingBatch
matnut$BreedingBatch <- factor(paste0("br",matnut$BBoriginal), levels=unique(paste0("br",matnut$BBoriginal)))
matnut$BEoriginal <- matnut$BehaviorBatch
matnut$BehaviorBatch <- factor(paste0("be",matnut$BEoriginal), levels=unique(paste0("be",matnut$BEoriginal)))
matnut$Cage <- factor(matnut$Cage)
matnut$DamID <- factor(paste0("d",matnut$DamID))
matnut$SireID <- factor(paste0("s",matnut$SireID))

matnut$OFBox <- factor(paste0("of",matnut$OFBox))
matnut$LDChamber <- factor(paste0("ld",matnut$LDChamber))
matnut$SIHOrder <- factor(paste0("sih",matnut$SIHOrder))
matnut$FSTChamber <- factor(paste0("fst",matnut$FSTChamber))
matnut$RestraintOrder <- factor(paste0("ro",matnut$RestraintOrder))
matnut$RestraintExperimenter <- factor(gsub(" ","",matnut$RestraintExperimenter))

orderdat <- matnut[order(matnut$RIX, matnut$Diet, matnut$PO),]
orderdat$PORIX <- paste0("PO",orderdat$RIX)
orderdat$PORIX <- factor(orderdat$PORIX, levels=unique(orderdat$PORIX))
orderdat$DietRIX <- paste0(orderdat$Diet, orderdat$RIX)
orderdat$DietRIX <- factor(orderdat$DietRIX, levels=unique(orderdat$DietRIX))
orderdat$PODietRIX <- paste0("PO",orderdat$DietRIX)
orderdat$PODietRIX <- factor(orderdat$PODietRIX, levels=unique(orderdat$PODietRIX))

variables <- c("Diet", "RIX", "PORIX", "DietRIX", "PODietRIX",
               "BehaviorBatch","DamID","SireID", "PO",
               "OFBox","LDChamber","SIHOrder","FSTChamber","RestraintOrder","RestraintExperimenter")
dispersions <- c("tauInv2","sigInv2", "tauPEInv2","tauDRInv2", "tauPDRInv2", 
                 "tauBEInv2","tauDamInv2","tauSirInv2")
#"tauOFInv2","tauLDInv2","tauSIHInv2","tauFSTInv2","tauROInv2","tauREInv2")

encoded <- getEncoding(orderdat, variables)
#write.csv(encoded, file.path("./","matnut_outputs/",'encoding.csv'), row.names = F)


ptypes <- colnames(matnut)[44:63]
#ptypes <- colnames(matnut)[c(44,48,53,54,55,61,62,63)]
#ptypes <- c("OFTotalDistance","BasalCORT","StressCORT")

allparam <- c(variables, dispersions)
allparam <- allparam[order(allparam)]
matnut <- list(df=orderdat, parameters=allparam, ptypes=ptypes, encoded=encoded)

saveRDS(matnut, file.path("./","matnut_outputs/",'matnut_data.rds'))
# many of the following functions expect list of data with 1) data, 2) parameters, 3) phenotypes, 4) encoding


##############
#  Sim data  #
##############

testdata <- getTestCase()
resultCompare <- testdata$effects[-c(which(testdata$effects$coef == 0), grep("tau", testdata$effects$variable)),]

firstDiet <- c()
for(i in 1:length(which(resultCompare$variable == "Diet"))){
  firstDiet[i] <- switch(as.numeric(resultCompare$level[i]), "AD", "BD", "CD", "DD")
}
firstDiet <- c(firstDiet, rep(1:9, 2))
use <- grep("[.]", resultCompare$level)
diet <- as.numeric(unlist(strsplit(resultCompare$level[use], "[.]"))[c(F, T)])
rix <- unlist(strsplit(resultCompare$level[use], "[.]"))[c(T, F)]
dietLet <- c()
for(i in 1:length(use)){
  dietLet[i] <- c(switch(diet[i], "AD", "BD", "CD", "DD"))
}

dietLet <- c(firstDiet, paste0(dietLet, rix))
labelNew <- ifelse(1:94 %in% grep("PO", resultCompare$variable), paste0("PO", dietLet), dietLet)

resultCompare$Level <- labelNew

mergeTab <- merge(resultCompare, testPhen$JL_compare, by="Level")
## compares real answers with lmer and jags estimates


testdataDF <- data.frame(testdata$cov.data)
testdataDF$DietOld <- testdataDF$Diet
for(i in 1:nrow(testdataDF)){
  testdataDF$Diet_use[i] <- switch(as.numeric(testdataDF$Diet)[i], "AD", "BD", "CD", "DD")
}
testdataDF$DamID <- factor(paste0("dd", testdataDF$DamID),levels=unique(paste0("dd", testdataDF$DamID)))
testdataDF$Diet <- as.factor(testdataDF$Diet_use)
testdataDF$DietRIX <- as.factor(paste0(testdataDF$Diet, testdataDF$RIX))
testdataDF$PORIX <- as.factor(paste0("PO", testdataDF$RIX))
testdataDF$PODietRIX <- as.factor(paste0("PO", testdataDF$DietRIX))
encodeTest <- getEncoding(testdataDF, c("DamID", "Diet", "RIX", "DietRIX", "PORIX", "PODietRIX")) #colnames(testdataDF)[c(1:3, 11:13)])
testPhen <- makeSummary(datalist=testdataDF, phenotype = "outcome", tryLam=c(1), normd=F, 
                        chains=2, n.adapt=20000, n.iter=50000, thin=10, encoded=encodeTest, sq=F)

predTest <- testdataDF[,c(7, 2, 3, 6, 8, 9)]

for(i in 1:nrow(predTest)){
  predTest$jagsDiet[i] <- testPhen$JL_compare$JAGS_est[which(testPhen$JL_compare$Level == paste(predTest$Diet[i]))]
  predTest$jagsRIX[i] <- testPhen$JL_compare$JAGS_est[which(testPhen$JL_compare$Level == paste(predTest$RIX[i]))]
  predTest$lmerDiet[i] <- testPhen$JL_compare$LMER_est[which(testPhen$JL_compare$Level == paste(predTest$Diet[i]))]
  predTest$lmerRIX[i] <- testPhen$JL_compare$LMER_est[which(testPhen$JL_compare$Level == paste(predTest$RIX[i]))]
}

predTest$jagsEps <- predTest$outcome - (predTest$jagsDiet + predTest$jagsRIX)
predTest$lmerEps <- predTest$outcome - (predTest$lmerDiet + predTest$lmerRIX)
## sum of effects and calculation of epsilon from estimates

jagsMu <- c()
for(i in 1:ncol(testPhen$mcmcObject)){
  jagsMu[i] <- density(testPhen$mcmcObject[,i])$x[which.max(density(testPhen$mcmcObject[,i])$y)]
}
source(file.path(".", "matnut", "summary_functions.R"))



