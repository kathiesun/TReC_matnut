dir = "/nas/depts/006/valdar-lab/users/sunk/trec"
library(rjags)

files = list.files(dir, full.names=T)

allKeep = list()
for(i in 1:length(files)){
	dat = readRDS(files[i])
	mins = unlist(lapply(dat, function(x) summary(as.mcmc.list(x$multimp$fit$mu_a_abs))[[2]][1]))
	keep = as.vector(which(mins > 0.01))
	kdat = lapply(keep, function(x) dat[[x]])
	names(kdat) = names(dat)[keep]
	allKeep = c(allKeep, kdat)
}

######################################


reg_trec = readRDS(file.path(dir, "trec/keep_genes_14sep2020.rds"))

meds = unlist(lapply(reg_trec, function(x) 
  HPDinterval(as.mcmc(unlist(x$multimp$fit$mu_a_abs)), prob=0.0001)))
meds = sort(unlist(meds), decreasing = T)

mins = unlist(lapply(reg_trec, function(x) 
  HPDinterval(as.mcmc(unlist(x$multimp$fit$mu_a_abs)))[1]))
mins = sort(unlist(mins), decreasing = T)
maxs = unlist(lapply(reg_trec, function(x) 
  HPDinterval(as.mcmc(unlist(x$multimp$fit$mu_a_abs)))[2]))
maxs = sort(unlist(maxs), decreasing = T)
modes = unlist(lapply(reg_trec, function(x){
  dd <- density(unlist(x$multimp$fit$mu_a_abs))
  dd$x[which.max(dd$y)]
  }))
modes = sort(modes, decreasing = T)

mids = intersect(names(meds)[1:500], names(modes)[1:500])
ends = intersect(names(mins)[1:500], names(maxs)[1:500])
intersect(mids, ends)

lapply(intersect(mids, ends), function(x) 
  plot(as.mcmc.list(reg_trec[[x]]$multimp$fit$mu_a_abs)))

plots = lapply(intersect(mids, ends), function(x){  
  clusters = get_clusters(reg_trec[[x]]$multimp$fit)
  
  #if(length(unique(clusters)) == 1){
  #  clus = plot_clust(x$multimp$fit[["mu_a"]],
  #                    x$multimp$fit[["mu_a_sq"]])
  clus = plot_clust_single(reg_trec[[x]]$multimp$fit[["mu_a_abs"]])
  #} else {
  #  clus = plot_clust(x$multimp$fit[["mu_0"]],
  #                    lapply(x$multimp$fit[["mu_a"]], abs))
  #}
  #plot(mcmc.list(lapply(x$multimp$fit$mu_a, function(x) (x^2)^0.5)))
  #plot(mcmc.list(lapply(x$multimp$fit$mu_a, function(x) abs(x))))
  
  #clusters = get_clusters(x$multimp$fit)$medians  
  ### only picks up one side
  rix = plot_rix(reg_trec[[x]]$matched_df, clusters$medians)
  return(list(clus, rix))
})


pdf("C:/Users/Kathie/Dropbox (ValdarLab)/phenotype_analysis/plots_trec_18sep2020.pdf")
for(n in 1:length(plots)){
  grid.arrange(plots[[n]][[2]][[1]],plots[[n]][[2]][[2]], top=names(plots)[n], ncol = 2, nrow = 1)
  print(plots[[n]][[1]])
  plot(reg_trec[[n]]$multimp$fit$mu_a_abs[[1]], main=names(plots)[n])
  plot(reg_trec[[n]]$multimp$fit$mu_a, main = names(plots)[n])
}
dev.off()
