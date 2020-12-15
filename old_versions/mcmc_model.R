library(lme4)
library(ggplot2)
source("./lm/fitBoxCoxModels.R")
library(effects)
library(MCMCglmm)

matnut <- read.csv("../data/AllPhenotypes_Matnut4.csv")
matnut[grep("a", matnut$Reciprocal),"poe"] <- 0.5
matnut[grep("b", matnut$Reciprocal),"poe"] <- -0.5
matnut[, "diet_f"] <- factor(matnut$Diet, levels=c("Standard", "LowProtein ", "MethylEnriched", "VitaminDDeficient"))
matnut[, "rix_f"] <- factor(matnut$RIX, levels=c(1,2,3,4,6,7,8,9,10))
matnut$breedbat_f = factor(matnut$BreedingBatch)
matnut$dam_f = factor(matnut$DamID)

orderdat <- matnut[order(matnut$rix_f, matnut$diet_f),]

ptypes <- colnames(matnut)[44:63]
vec <- c(1,2,18:20)
phenotype <- ptypes[vec]

##################

prior <- list(R = list(V = 1e-16, nu = -2), 
              G = list(G1 = list(V = 1e-16, nu = -2), 
                       G2 = list(V=diag(4), nu=4, alpha.mu=rep(0,4), alpha.V=diag(4)*1000),
                       G3 = list(V = 1e-16, nu = -2)))
              #,G4 = list(V=diag(8), nu=8, alpha.mu=rep(0,8), alpha.V=diag(8)*1000)))
mod_all <- list()
all_plot_dat <-list()

##################
for (i in 1:length(phenotype)){
  use <- which(is.na(orderdat[,phenotype[i]]) == F)
  matnut_use <- orderdat[use,]
  pheno = matnut_use[,phenotype[i]]
  pheno = scale(pheno)
  
  mod1<- MCMCglmm(fixed = pheno ~ diet_f - 1,
                  family = "gaussian", prior=prior,
                  random = ~ rix_f + us(diet_f):rix_f + us(poe):rix_f,
                  nitt = 6000, burnin = 1000, thin=5, data = matnut_use, pr=T)
  
  dietrix <- intersect(grep("diet",colnames(mod1$Sol)), grep("rix",colnames(mod1$Sol)))
  onlyrix <- setdiff(grep("rix",colnames(mod1$Sol)), grep(".rix",colnames(mod1$Sol)))
  onlydiet <- setdiff(grep("diet",colnames(mod1$Sol)), grep(".rix",colnames(mod1$Sol)))
  #diet <- c(rep("LowProtein", 9), rep("MethylEnriched",9),rep("Standard", 9), rep("VitaminDDeficient",9))
  #rix <- ordered(rep(c(1:4,6:10),4), levels=c(1:4,6:10))
  rixef <- rep(apply(mod1$Sol[,onlyrix], 2, median),4)
  dietef <- rep(apply(mod1$Sol[,onlydiet], 2, median), each=9)
  means <- apply(mod1$Sol[,dietrix], 2, median)
  #xax <- ordered(colnames(mod1$Sol)[dietrix], levels=colnames(mod1$Sol)[dietrix])
  #plotdat <- as.data.frame(cbind(phenotype[i], xax, as.factor(diet), rix))
  #colnames(plotdat) <- c("pheno", "xax", "diet", "rix")
  interval <- data.frame(HPDinterval(mod1$Sol)[dietrix,])
  interval_thin <- data.frame(HPDinterval(mod1$Sol, prob=0.5)[dietrix,])
  names(interval_thin) <- c("lowThin", "upThin")
  all_plot_dat[[i]] <- cbind(means, interval, interval_thin, rixef, dietef)
  all_plot_dat[[i]]$names <- colnames(mod1$Sol)[dietrix]
}

## copy rand ef plot from lmer
#tot_copy_rand <- rbindlist(all_plot_dat)

plot_rand_mcmc <- list()
i=5
for(i in length(all_plot_dat)){

  tot_copy_rand <- all_plot_dat[[i]]
  #tot_copy_rand$rix <- factor(rep(c(1:4,6:10),4), order(c(1:10)))
  rix <- rep(ordered(rep(c(1:4,6:10),4), levels=c(1:4,6:10)))
  tot_copy_rand$rix <- rix
  diets <- ordered(rep(c("Std","lowPr","Me","VitD"), each=9), levels=c("Std","lowPr","Me","VitD"))
  tot_copy_rand$diets <- diets
  tot_copy_rand <- tot_copy_rand[order(tot_copy_rand$rix, tot_copy_rand$diets),]
  tot_copy_rand$names <- factor(tot_copy_rand$names, levels=tot_copy_rand$names)
  
  #########################
  # This plot below works #
  #########################
  ## individual cat plots
  plot_rand_mcmc[[i]] <- ggplot(data=tot_copy_rand) + 
    geom_point(aes(x=names, y=means, group=diets, colours=diets)) +
    geom_line(aes(x=names, y=rixef, group=rix), col="black") + 
    geom_point(aes(x=names, y=dietef, group=rix), col="lightgray", shape="|", size=3) + 
    labs(title=paste("Random effects by RIX in",phenotype[[i]])) + 
    geom_segment(aes(x=names, xend=names, y=lower, yend=upper, group=diets, col=diets), size=1) +
    geom_segment(aes(x=names, xend=names, y=lowThin, yend=upThin, group=diets, col=diets), size=1.5) +
    theme_bw() + coord_flip()
  print(plot_rand_mcmc[[i]])

}

tot_copy_rand <- tot_copy_rand[order(tot_copy_rand$diets, tot_copy_rand$rix),]
tot_copy_rand$names <- factor(tot_copy_rand$names, levels=tot_copy_rand$names)

plot_rand_mcmc <- ggplot() + 
  geom_point(data=tot_copy_rand, aes(x=names, y=means, group=rix, col=rix)) +
  geom_line(data=tot_copy_rand, aes(x=names, y=dietef, group=diets), col="black") + 
  geom_point(data=tot_copy_rand, aes(x=names, y=rixef, group=diets), col="lightgray", size=2, alpha=0.8) + 
  labs(title=paste("Random effects by RIX in",phenotype[[1]])) + 
  geom_segment(data=tot_copy_rand, aes(x=names, xend=names, y=lower, yend=upper, group=rix, col=rix)) +
  theme_bw() + coord_flip()
print(plot_rand_mcmc)

####################################
# Not sure if anything below works #
####################################


#tot_copy_rand$diet <- factor(rep(diets, each=9), diets)
#tot_copy_rand$pheno <- rep(phenotype, each=9)

copy_rand_mcmc <- ggplot() + 
  geom_point(data=tot_copy_rand, aes(x=rix, y=means, group=diet, col=diet)) +
  geom_line(data=tot_copy_rand, aes(x=rix, y=means, group=diet, col=diet)) + 
  labs(title=paste("Random effects by RIX from mcmc output")) + 
  theme_bw() + facet_wrap(~pheno, scales = "free_y")
print(copy_rand_mcmc) 



plot_randci_mcmc <- ggplot() + 
  geom_point(data=tot_copy_rand, aes(x=xax, y=means, group=diets, col=diets)) +
  #geom_point(data=tot_copy_rand, aes(x=xax, y=rixef), col="black") +
  geom_line(data=tot_copy_rand, aes(x=xax, y=rixef, group=diets), col="black") + 
  labs(title=paste("Standardized random effects by RIX from mcmc output"), 
       xlab=xax_2) + 
  geom_segment(data=tot_copy_rand, aes(x=xax, xend=xax, y=lower, yend=upper, group=diets, col=diets)) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + facet_wrap(~pheno) 
  #+ coord_flip()
print(plot_randci_mcmc) 

library(mcmcplots)
caterplot(all_plot_dat[[1]], reorder = F)
#traceplot(mod1$Sol)
#HPDinterval(mod1$Sol)
