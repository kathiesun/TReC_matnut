library(lme4)
library(ggplot2)
source("./lm/fitBoxCoxModels.R")
library(multcomp)
library(MCMCglmm)

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

normd=T
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



summary(glht(diet_normd$mod_results[[1]]$phen_1$fit, mcp(diet_f="Tukey")))

