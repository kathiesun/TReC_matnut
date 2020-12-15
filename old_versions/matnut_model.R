library(lme4)
library(ggplot2)
source("./lm/fitBoxCoxModels.R")
library(effects)
library(MCMCglmm)

getMatnutTransformParams <- function(normd=T, tryLam)
{
  out = list()
  out$lambdasToTry = tryLam
  out$lambdaPerY      = NULL
  if (normd){
    out$normalizeBeforeTransform = T
    out$normalizeAfterTransform  = T
  }
  else {
    out$normalizeBeforeTransform = F
    out$normalizeAfterTransform  = F
  }
  out$extremelb      = -3
  out$extremeub      = 3
  return(out)
}

#######################
matnut <- read.csv("../data/AllPhenotypes_Matnut4.csv")
matnut[grep("a", matnut$Reciprocal),"poe"] <- 0.5
matnut[grep("b", matnut$Reciprocal),"poe"] <- -0.5
matnut[, "diet_f"] <- factor(matnut$Diet) 
matnut[, "rix_f"] <- factor(matnut$RIX)
matnut$breedbat_f = factor(matnut$BreedingBatch)
matnut$dam_f = factor(matnut$DamID)



ptypes <- colnames(matnut)[44:63]
vec <- c(1,2,13:15,18:20)
phenotype <- ptypes[vec]

dat_rand <- list()
ranef_plots <- list()
fixef_plots <- list()
dat_fix <- list()

########################
plotRanEfs <- function(normd=T, plots, phenotype, fix=T, indvariables, tryLam){
  mod_results <- list()
  dat_rand <- list()
  ranef_plots <- list()
  fixef_plots <- list()
  dat_fix <- list()
  
  for (i in 1:length(phenotype)){
    pheno <- phenotype[i]
    use <- which(is.na(matnut[,pheno]) == F)
    matnut_use <- matnut[use,]
    #mod_results <- lmer(matnut_use[,pheno] ~ -1 + diet_f + 
    #                        (1 | rix_f) + (poe - 1 | rix_f) + (diet_f - 1 | rix_f), matnut_use)
    #randoms[[i]] <- ranef(mod_results[[i]])
    mod_results[[i]] <- fit.model.bc$fit(y.mat=matnut_use[,pheno], matnut_use, indvariables, 
                                         transformParams = getMatnutTransformParams(normd=normd, tryLam = tryLam))
    if (!is.null(mod_results[[i]]$phen_1$fit)){
      ran <- as.data.frame(ranef(mod_results[[i]]$phen_1$fit)$rix_f)
      ran2 <- cbind(factor(rownames(ran), order(c(1:10))), ran)
      colnames(ran2)[1] <- "rix"
    
      melt_dat <- melt(ran2, id="rix")
      diet <- grep("diet",melt_dat$"variable")
      poe <- grep("poe",melt_dat$"variable")
      
      if(plots == "diet"){
        dat_rand[[i]] <- cbind(pheno, melt_dat[diet,])
      } else {
        dat_rand[[i]] <- cbind(pheno, melt_dat[poe,])
      }
      
      if (normd){
        transtr <- "Standardized"
      } else {
        transtr <- ""
      }
      
      p <- ggplot() + 
        geom_point(data=dat_rand[[i]], aes(x=rix, y=value, group=variable, col=variable)) +
        labs(x="RIX", y="", title=paste(transtr, pheno)) + 
        theme_bw()
      ranef_plots[[i]] <- p

      if (fix){
        conf <- as.data.frame(confint(mod_results[[i]]$phen_1$fit, "beta_", method="Wald"))
        point <- as.data.frame(fixef(mod_results[[i]]$phen_1$fit))
        dat_fix[[i]] <- cbind(pheno, levels(matnut$diet_f), conf, point)
        names(dat_fix[[i]]) <- c("pheno","Diet", "low", "high", "Mean")
        
        pfix <- ggplot() +
          geom_point(data=dat_fix[[i]], aes(x=Diet, y=Mean), size=3) +
          geom_segment(data=dat_fix[[i]], aes(x=Diet, xend=Diet, y=low, yend=high)) +
          labs(x="Diet", y="", title=paste(transtr, pheno)) + 
          theme_bw()
        fixef_plots[[i]] <- pfix
      }
    }
  }
  tot_dat_fix <- rbindlist(dat_fix)
  tot_dat_rand <- rbindlist(dat_rand)
  return(list(tot_dat_rand=tot_dat_rand, tot_dat_fix=tot_dat_fix, mod_results=mod_results))
}

###############
# start plots #
###############

normd=F
diet_normd <- plotRanEfs(normd=normd, plots="diet", phenotype=phenotype, fix=F, tryLam = 1, #c(-3, -2, -1, -0.5, 0, 0.5, 1, 2, 3),
                         "~ -1 + diet_f + (1 | rix_f) + (poe - 1 | rix_f) + (diet_f - 1 | rix_f)")

if(normd==T){transtr="Standardized"} else{transtr=""}

plot_fix_tot <- ggplot() +
  geom_point(data=diet_normd$tot_dat_fix, aes(x=Diet, y=Mean), size=3) +
  geom_segment(data=tot_dat_fix, aes(x=Diet, xend=Diet, y=low, yend=high)) +
  labs(title=paste(transtr, "fixed effects of diet")) + 
  theme_bw() + facet_wrap(~pheno)
print(plot_fix_tot)

plot_rand_tot <- ggplot() + 
  geom_point(data=diet_normd$tot_dat_rand, aes(x=rix, y=value, group=variable, col=variable)) +
  geom_line(data=diet_normd$tot_dat_rand, aes(x=rix, y=value, group=variable, col=variable)) + 
  labs(title=paste(transtr, "Random Effects by RIX from lmer output")) + 
  theme_bw() + facet_wrap(~pheno, scales = "free_y")
print(plot_rand_tot)

for (i in 1:20){
  ifelse(!is.null(mod_results[[i]]$phen_1$fit), print(mod_results[[i]]$phen_1$lambda), print(paste(i,"didn't work")))
}

##################
prior <- list(R = list(V = 1e-16, nu = -2), 
              G = list(G1 = list(V = 1e-16, nu = -2), 
                       G2 = list(V=diag(4), nu=4, alpha.mu=rep(0,4), alpha.V=diag(4)*1000),
                       G3 = list(V = 1e-16, nu = -2),
                       G4 = list(V=diag(8), nu=8, alpha.mu=rep(0,8), alpha.V=diag(8)*1000)))
mod_all <- list()
all_plot_dat <-list()
for (i in 1:length(phenotype)){
  use <- which(is.na(matnut[,phenotype[i]]) == F)
  matnut_use <- matnut[use,]
  pheno = matnut_use[,phenotype[i]]
  pheno = scale(pheno)
  
  mod1<- MCMCglmm(fixed = pheno ~ diet_f - 1,
                   family = "gaussian", prior=prior,
                   random = ~ rix_f + us(diet_f):rix_f + us(poe):rix_f + us(diet_f*poe):rix_f,
                   nitt = 6000, burnin = 1000, thin=5, data = matnut_use, pr=T)
  
  dietrix <- intersect(grep("diet",colnames(mod1$Sol)), grep("rix",colnames(mod1$Sol)))
  onlyrix <- setdiff(grep("rix",colnames(mod1$Sol)), grep(".rix",colnames(mod1$Sol)))
  diet <- c(rep("LowProtein", 9), rep("MethylEnriched",9),rep("Standard", 9), rep("VitaminDDeficient",9))
  rix <- ordered(rep(c(1:4,6:10),4), levels=c(1:4,6:10))
  rixef <- rep(apply(mod1$Sol[,onlyrix], 2, median),4)
  means <- apply(mod1$Sol[,dietrix], 2, median)
  xax <- ordered(colnames(mod1$Sol)[dietrix], levels=colnames(mod1$Sol)[dietrix])
  plotdat <- as.data.frame(cbind(phenotype[i], xax, as.factor(diet), rix))
  colnames(plotdat) <- c("pheno", "xax", "diet", "rix")
  interval <- HPDinterval(mod1$Sol)[dietrix,]
  all_plot_dat[[i]] <- cbind(plotdat, means, interval, rixef)
}

## copy rand ef plot from lmer
tot_copy_rand <- rbindlist(all_plot_dat)

#tot_copy_rand$rix <- factor(rep(c(1:4,6:10),4), order(c(1:10)))
rix <- rep(ordered(rep(c(1:4,6:10),4), levels=c(1:4,6:10)),length(phenotype))
tot_copy_rand$rix <- rix
xax <- rep(ordered(c("L1","L2","L3","L4","L6","L7","L8","L9","L10",
                     "M1","M2","M3","M4","M6","M7","M8","M9","M10",
                     "S1","S2","S3","S4","S6","S7","S8","S9","S10",
                     "V1","V2","V3","V4","V6","V7","V8","V9","V10"), 
                   levels=c("L1","L2","L3","L4","L6","L7","L8","L9","L10",
                            "M1","M2","M3","M4","M6","M7","M8","M9","M10",
                            "S1","S2","S3","S4","S6","S7","S8","S9","S10",
                            "V1","V2","V3","V4","V6","V7","V8","V9","V10")),
                   length(phenotype))
xax_2 <- rep(ordered(c("L1","L2","L3","L4","L6","L7","L8","L9","L10",
                            "M1","M2","M3","M4","M6","M7","M8","M9","M10",
                            "S1","S2","S3","S4","S6","S7","S8","S9","S10",
                            "V1","V2","V3","V4","V6","V7","V8","V9","V10"), 
                     levels=c("L1","M1", "S1","V1", "L2","M2", "S2","V2", "L3","M3", "S3","V3",
                              "L4","M4", "S4","V4", "L6","M6", "S6","V6", "L7","M7", "S7","V7",
                              "L8","M8", "S8","V8", "L9","M9", "S9","V9", "L10","M10", "S10","V10")),
                    length(phenotype))
tot_copy_rand$xax <- xax_2
tot_copy_rand$diet <- as.factor(c(rep("LowProtein", 9), rep("MethylEnriched",9),
                                  rep("Standard", 9), rep("VitaminDDeficient",9)))

copy_rand_mcmc <- ggplot() + 
  geom_point(data=tot_copy_rand, aes(x=rix, y=means, group=diet, col=diet)) +
  geom_line(data=tot_copy_rand, aes(x=rix, y=means, group=diet, col=diet)) + 
  labs(title=paste("Random effects by RIX from mcmc output")) + 
  theme_bw() + facet_wrap(~pheno, scales = "free_y")
print(copy_rand_mcmc) 



plot_randci_mcmc <- ggplot() + 
  geom_point(data=tot_copy_rand, aes(x=xax, y=means, group=rix, col=diet)) +
  #geom_point(data=tot_copy_rand, aes(x=xax, y=rixef), col="black") +
  geom_line(data=tot_copy_rand, aes(x=xax, y=rixef, group=rix), col="black") + 
  labs(title=paste("Random effects by RIX from mcmc output"), 
       xlab=c("L1","M1", "S1","V1", "L2","M2", "S2","V2", "L3","M3", "S3","V3",
                       "L4","M4", "S4","V4", "L6","M6", "S6","V6", "L7","M7", "S7","V7",
                       "L8","M8", "S8","V8", "L9","M9", "S9","V9", "L10","M10", "S10","V10")) + 
  geom_segment(data=tot_copy_rand, aes(x=xax, xend=xax, y=lower, yend=upper, group=diet, col=diet)) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + facet_wrap(~pheno) + 
  coord_flip()
print(plot_randci_mcmc) 

## individual cat plots
plot_rand_mcmc <- ggplot() + 
  geom_point(data=plotdat, aes(x=xax, y=means, group=diet, col=diet)) +
  #geom_line(data=plotdat, aes(x=xax, y=means, group=diet, col=diet)) + 
  labs(title=paste("Random effects by RIX in",phenotype[[1]])) + 
  #geom_segment(data=plotdat, aes(x=xax, xend=xax, y=lower, yend=upper, group=diet, col=diet)) +
  theme_bw()
print(plot_rand_mcmc)


#traceplot(mod1$Sol)
#HPDinterval(mod1$Sol)
