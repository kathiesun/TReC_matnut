muClust1 = data.frame(do.call("cbind", full_list[["mu_0_new"]]))
colnames(muClust1) = paste0("it",seq(1:ncol(muClust1)))
muClust1_long = gather(muClust1, "it","est")
ggplot(muClust1_long, aes(x=est, col=it)) + geom_density() + xlim(-5,5)

muClust2 = data.frame(do.call("cbind", full_list[["mu_a_new"]]))
colnames(muClust2) = paste0("it",seq(1:ncol(muClust2)))
muClust2_long = gather(muClust2, "it","est")
ggplot(muClust2_long, aes(x=est, col=it)) + geom_density()
plot(full_list[["mu_0_new"]])

#multimp[[ptypes[i]]]$madeModel <- madeModel


###############################################




for(t in 1: nrow(fit[[1]])){
  for(n in 1:length(dataList$y)){
    nul1 = dnorm(dataList$y[n], mean=unlist(fit[t,"mu_0"]), sd=(unlist(fit[t,"tau"])^-0.5)) * 
      (1-unlist(fit[t,"p"]))
    nul2 = dnorm(dataList$y[n], mean=(-1*unlist(fit[t,"mu_0"])), sd=(unlist(fit[t,"tau"])^-0.5)) * 
      (1-unlist(fit[t,"p"]))
    nul = max(nul1, nul2)
    alt1 = dnorm(dataList$y[n], mean=unlist(fit[t,"mu_a"]), sd=(unlist(fit[t,"tau"])^-0.5)) * 
      unlist(fit[t,"p"])
    alt2 = dnorm(dataList$y[n], mean=(-1*unlist(fit[t,"mu_a"])), sd=(unlist(fit[t,"tau"])^-0.5)) * 
      unlist(fit[t,"p"])
    alt = max(alt1, alt2)
    Q1mat[n,1,t] = nul/(nul+alt) 
    Q1mat[n,2,t] = alt/(nul+alt) 
    
    
    nul1 = dnorm(dataList$y[n], mean=unlist(fit[t,"mu_a"]), sd=(unlist(fit[t,"tau"])^-0.5)) * 
      (1-unlist(fit[t,"p"]))
    nul2 = dnorm(dataList$y[n], mean=(-1*unlist(fit[t,"mu_a"])), sd=(unlist(fit[t,"tau"])^-0.5)) * 
      (1-unlist(fit[t,"p"]))
    nul = max(nul1, nul2)
    alt1 = dnorm(dataList$y[n], mean=unlist(fit[t,"mu_0"]), sd=(unlist(fit[t,"tau"])^-0.5)) * 
      unlist(fit[t,"p"])
    alt2 = dnorm(dataList$y[n], mean=(-1*unlist(fit[t,"mu_0"])), sd=(unlist(fit[t,"tau"])^-0.5)) * 
      unlist(fit[t,"p"])
    alt = max(alt1, alt2)
    Q2mat[n,1,t] = nul/(nul+alt) 
    Q2mat[n,2,t] = alt/(nul+alt) 
    
  }
}


v1 = sapply(1:dim(Q1mat)[3], function(t) {
  q1_0 = sum(Q1mat[,1,t])/dim(Q1mat)[1]
  q1_a = sum(Q1mat[,2,t])/dim(Q1mat)[1]
  
  sum((log(Q1mat[,1,t]/q1_0) * Q1mat[,1,t]) + (log(Q1mat[,2,t]/q1_a) * Q1mat[, 2,t]))
})

v2 = sapply(1:dim(Q2mat)[3], function(t) {
  q2_0 = sum(Q2mat[,1,t])/dim(Q2mat)[1]
  q2_a = sum(Q2mat[,2,t])/dim(Q2mat)[1]
  
  sum((log(Q2mat[,1,t]/q2_0) * Q2mat[,1,t]) + (log(Q2mat[,2,t]/q2_a) * Q2mat[, 2,t]))
})





v1 = sapply(1:dim(Q1mat)[3], function(t) {
  q1_0 = sum(Q1mat[,1,t])/dim(Q1mat)[1]
  q1_a = sum(Q1mat[,2,t])/dim(Q1mat)[1]
  
  sum((log(Q1mat[,1,t]/q1_0) * Q1mat[,1,t]) + (log(Q1mat[,2,t]/q1_a) * Q1mat[, 2,t]))
})

v2 = sapply(1:dim(Q2mat)[3], function(t) {
  q2_0 = sum(Q2mat[,1,t])/dim(Q2mat)[1]
  q2_a = sum(Q2mat[,2,t])/dim(Q2mat)[1]
  
  sum((log(Q2mat[,1,t]/q2_0) * Q2mat[,1,t]) + (log(Q2mat[,2,t]/q2_a) * Q2mat[, 2,t]))
})









q1_min = which.min(c(-sum(log(Q1mat[,1,t]) + log(Q1mat[,2,t]))))
q2_min = which.min(c(-sum(log(Q2mat[,1,t])), -sum(log(Q2mat[,2,t]))))
which_q = which.min(c(c(-sum(log(Q1mat[,1,t])), -sum(log(Q1mat[,2,t])))[q1_min], 
                      c(-sum(log(Q2mat[,1,t])), -sum(log(Q2mat[,2,t])))[q2_min]))
mu0_mua = ifelse(which_q == 1, ifelse(q1_min == 1, 0, "a"), 
                 ifelse(q2_min == 1, 0, "a"))
c(which_q, mu0_mua)


v_comb = cbind(v1, v2)
which_min = apply(v_comb, 1, which.min)

mu_out = t(sapply(1:length(v1), function(t) {  
  if(which_min[t] == 1){
    cbind(unlist(fit[,"mu_0"])[t], unlist(fit[,"mu_a"])[t])
  } else {
    cbind(unlist(fit[,"mu_a"])[t], unlist(fit[,"mu_0"])[t])
  }
}))


-log(sum(Q2mat[,1,1]))
-log(sum(Q2mat[,2,1]))


####################################################################



data("mcmc_output")




L_mat = apply(fit[[1]], 1, function(t){
  do.call("rbind", lapply(dataList$y, function(n){
    mu_0_use = ifelse(n > 0, abs(t["mu_0"]), -1*abs(t["mu_0"]))
    pr_nul1 = dnorm(n, mean=mu_0_use, sd=t["tau"]^-0.5) * (1-t["p"])
    
    mu_a_use = ifelse(n > 0, abs(t["mu_a"]), -1*abs(t["mu_a"]))
    pr_alt1 = dnorm(n, mean=mu_a_use, sd=t["tau"]^-0.5) * (t["p"])
    
    mu_0_use = ifelse(n > 0, abs(t["mu_a"]), -1*abs(t["mu_a"]))
    pr_nul2 = dnorm(n, mean=mu_0_use, sd=t["tau"]^-0.5) * (1-t["p"])
    
    mu_a_use = ifelse(n > 0, abs(t["mu_0"]), -1*abs(t["mu_0"]))
    pr_alt2 = dnorm(n, mean=mu_a_use, sd=t["tau"]^-0.5) * (t["p"])
    
    
    c(pr_nul1, pr_alt1, pr_nul2, pr_alt2)
  }))
})




apply(Lmat, 2 function(t){
  for(i in 1:207){
    
  }
})





q_ij_0 = do.call("rbind", lapply(dataList$y, function(n){
  apply(t(apply(fit[[1]], 1, function(t){
    pr_nul1 = dnorm(n, mean=t["mu_0"], sd=t["tau"]^-0.5) * (1-t["p"])
    pr_alt1 = dnorm(n, mean=t["mu_a"], sd=t["tau"]^-0.5) * (t["p"])
    c(pr_nul1/(pr_nul1 + pr_alt1), pr_alt1/(pr_nul1 + pr_alt1))
  })),2,sum)/nrow(fit[[1]])
}))


dim(p_mat_0)


n=length(dataList$y)
to_flip = apply(L_mat, 2, function(t){
  v0 = -log(prod(apply(cbind(t[1:n], t[(n+1):(2*n)]),1, max)))
  v0_grp = apply(cbind(t[1:n], t[(n+1):(2*n)]),1, which.max)
  v1 = -log(prod(apply(cbind(t[(2*n+1):(3*n)], t[(3*n+1):(4*n)]),1, max)))
  v1_grp = apply(cbind(t[(2*n+1):(3*n)], t[(3*n+1):(4*n)]),1, which.max)
  grp_fin = v1_grp
  if(v0 < v1) grp_fin = v0_grp
  list(flip = ifelse(v0 < v1, "v0", "v1"), grp_fin)
})

flips = unlist(lapply(to_flip, function(x) x[[1]]))
mu_out = t(sapply(1:length(flips), function(t) {  
  if(flips[t] == "v0"){
    cbind(unlist(fit[,"mu_0"])[t], unlist(fit[,"mu_a"])[t])
  } else {
    cbind(unlist(fit[,"mu_a"])[t], unlist(fit[,"mu_0"])[t])
  }
}))


Q1mat = array(rep(NaN), c(length(dataList$y), Nclust, nrow(fit[[1]])))
Q2mat = array(rep(NaN), c(length(dataList$y), Nclust, nrow(fit[[1]])))

for(t in 1: nrow(fit[[1]])){
  for(n in 1:length(dataList$y)){
    mu_0_use = ifelse(dataList$y[n] > 0, abs(unlist(fit[t,"mu_0"])), -1*abs(unlist(fit[t,"mu_0"])))
    pr_nul1 = dnorm(dataList$y[n], mean=mu_0_use, sd=(unlist(fit[t,"tau"])^-0.5)) * 
      (1-unlist(fit[t,"p"]))
    
    mu_a_use = ifelse(dataList$y[n] > 0, abs(unlist(fit[t,"mu_a"])), -1*abs(unlist(fit[t,"mu_a"])))
    pr_alt1 = dnorm(dataList$y[n], mean=mu_a_use, sd=(unlist(fit[t,"tau"])^-0.5)) * 
      unlist(fit[t,"p"])
    Qmat[n,1,t] = -log(pr_nul1*pr_alt1)
    
    
    mu_0_use = ifelse(dataList$y[n] > 0, abs(unlist(fit[t,"mu_a"])), -1*abs(unlist(fit[t,"mu_a"])))
    pr_nul2 = dnorm(dataList$y[n], mean=mu_0_use, sd=(unlist(fit[t,"tau"])^-0.5)) * 
      (1-unlist(fit[t,"p"]))
    
    mu_a_use = ifelse(dataList$y[n] > 0, abs(unlist(fit[t,"mu_a"])), -1*abs(unlist(fit[t,"mu_a"])))
    pr_alt2 = dnorm(dataList$y[n], mean=mu_a_use, sd=(unlist(fit[t,"tau"])^-0.5)) * 
      unlist(fit[t,"p"])
    Qmat[n,2,t] = -log(pr_nul2*pr_alt2)
    
  }
}




which_min = apply(Qmat, 1, which.min)

mu_out = t(sapply(1:nrow(which_min), function(t) {  
  if(which_min[t] == 1){
    cbind(unlist(fit[,"mu_0"])[t], unlist(fit[,"mu_a"])[t])
  } else {
    cbind(unlist(fit[,"mu_a"])[t], unlist(fit[,"mu_0"])[t])
  }
}))

###########################################
alt / (nul+alt)

nul = dnorm(dataList$y[n], mean=unlist(fit[1,"mu_0"]), sd=(unlist(fit[1,"tau"])^-0.5)) * 
  unlist(fit[1,paste0("clust[",dataList$indRIX[n],"]")])
alt = dnorm(dataList$y[n], mean=unlist(fit[1,"mu_a"]), sd=(unlist(fit[1,"tau"])^-0.5)) * 
  (1-unlist(fit[1,paste0("clust[",dataList$indRIX[n],"]")]))

nul = dnorm(dataList$y[n], mean=unlist(fit[1,"mu_a"]), sd=(unlist(fit[1,"tau"])^-0.5)) * 
  unlist(fit[1,paste0("clust[",dataList$indRIX[n],"]")])
alt = dnorm(dataList$y[n], mean=unlist(fit[1,"mu_0"]), sd=(unlist(fit[1,"tau"])^-0.5)) * 
  (1-unlist(fit[1,paste0("clust[",dataList$indRIX[n],"]")]))
