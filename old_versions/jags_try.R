library(rjags)
library(coda)
library(data.table)
library(xtable)
library(grid)
library(beanplot)

source("./plot.hpd.R")
source("./jags_model.R")
source("./matnut/encode.R")
# setwd("C:/Users/kathie/Dropbox/source/nutridiallel2014/src")

################
# Prepare data #
################

matnut <- read.csv("../data/AllPhenotypes_Matnut5.csv")
matnut[grep("a", matnut$Reciprocal),"poe"] <- 0.5
matnut[grep("b", matnut$Reciprocal),"poe"] <- -0.5
matnut[which(matnut$poe == 0.5), "pof"] <- "+"
matnut[which(matnut$poe == -0.5), "pof"] <- "-"

dietlabs <- c("STD", "ME", "PD", "VDD")
#rixlabs <- c("RIX1", "RIX2", "RIX3", "RIX4", "RIX6", "RIX7", "RIX8", "RIX9", "RIX10")
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
matnut$RIX <- factor(matnut$RIX, levels=c(1:4,6:10))
matnut$Diet <- factor(matnut$Diet, levels=dietlabs)
matnut$DamID <- factor(matnut$DamID)
matnut$BreedingBatch <- factor(matnut$BreedingBatch)

orderdat <- matnut[order(matnut$RIX, matnut$Diet, matnut$poe),]
orderdat$RIXPO <- paste0("PO.",orderdat$RIX)
orderdat$RIXPO <- factor(orderdat$RIXPO, levels=unique(orderdat$RIXPO))
orderdat$DietRIX <- paste0(orderdat$Diet, orderdat$RIX)
orderdat$DietRIX <- factor(orderdat$DietRIX, levels=unique(orderdat$DietRIX))
orderdat$DietRIXPO <- paste0("PO.",orderdat$DietRIX)
orderdat$DietRIXPO <- factor(orderdat$DietRIXPO, levels=unique(orderdat$DietRIXPO))

encoded <- getEncoding(orderdat, c("Diet", "RIX", "RIXPO", "DietRIX", "DietRIXPO"))

ptypes <- colnames(matnut)[44:63]
vec <- c(1,2,15,18:20)
phenotype <- "BasalCORT" #ptypes[vec]

#for(i in 1:length(phenotype)){
#  hist(matnut[,phenotype[i]], main=paste("Histogram of",phenotype[i]))
#}
# check: c("RIX", "Reciprocal", "Diet", "Diet_or", "Diet_f", "RIX_f", "poe", "names")
######################
# Loop through model #
######################
list.reg <- list()
for (i in 1:length(phenotype)){
  orderdat <- matnut
  i=1
  use <- which(is.na(orderdat[,phenotype[i]]) == F)
  matnut_use <- orderdat[use,]
  pheno = matnut_use[,phenotype[i]]
  pheno = scale(pheno)
  
  N <- length(pheno)
  y <- pheno
  
  diet <- as.numeric(matnut_use$Diet)
  rix <- as.numeric(matnut_use$RIX)
  nrix <- length(unique(rix))
  ndiet <- length(unique(diet))
  poe <- matnut_use$poe
  dietrix <- as.numeric(matnut_use$DietRIX)
  ndietrix <- length(unique(dietrix))
  
  data = list("N"=N, "rix"=rix, "y"=y, "nrix"=nrix, "diet"=diet, "ndiet"=ndiet, 
              "poe"=poe, "dietrix"=dietrix, "ndietrix"=ndietrix)
  
#############
# Run model #
#############
  
  allparam <- c("sigInv2", "tauInv2","tauDRInv2", "tauPEInv2", "tauPDRInv2", "RIX", "Diet", 
                "DietRIX", "RIXPO", "DietRIXPO")
  allparam <- allparam[order(allparam)] 
  chains <- 2 
  reg.jags <- jags.model(textConnection(model), data=data, n.chains = chains, n.adapt = 2500)
  update(reg.jags, n.iter=2500)
  list.reg[[i]] <- coda.samples(reg.jags, variable.names = allparam, thin=10, n.iter=10000)
 
}
##########
# Labels #
##########

allnames <- c()
wantParam <- c("RIX","RIXPO","Diet","DietRIX","DietRIXPO")
for(i in 1:length(wantParam)){
  wantParam <- allparam[which(wantParam %in% allparam)]
  allnames <- c(allnames, as.character(encoded$Level[which(encoded$Variable == wantParam[i])])) 
}


#for(i in 1:length(list.reg))
#{
#  varnames(list.reg[[1]])[1:length(allnames)] <- allnames
#}  


#pick.params <- allnames[-grep("Inv", allnames)]
hasPO <- allnames[grep("PO", allnames)]
noPO <- allnames[-grep("PO", allnames)]
PO.R <- allnames[grep("PO.R", allnames)]

useset <- noPO
lowp <- useset[grep("PD.", useset)]
methyl <- useset[grep("ME.", useset)]
std <- useset[grep("STD.", useset)]
vitd <- useset[grep("VDD.", useset)]

dietrix <- c(std, methyl, lowp, vitd)


RIX10 <- dietrix[intersect(grep("1",dietrix),grep("10",dietrix))]
RIX1 <- dietrix[setdiff(grep("1",dietrix),grep("10",dietrix))]
RIX2 <- dietrix[grep("2",dietrix)]
RIX3 <- dietrix[grep("3",dietrix)]
RIX4 <- dietrix[grep("4",dietrix)]
RIX6 <- dietrix[grep("6",dietrix)]
RIX7 <- dietrix[grep("7",dietrix)]
RIX8 <- dietrix[grep("8",dietrix)]
RIX9 <- dietrix[grep("9",dietrix)]

######
podietrix <- paste0("PO.", dietrix) 
po.lowp <- podietrix[grep("PD.", podietrix)]
po.methyl <- podietrix[grep("ME.", podietrix)]
po.std <- podietrix[grep("STD.", podietrix)]
po.vitd <- podietrix[grep("VDD.", podietrix)]

po.RIX1 <- paste0("PO.", RIX1) 
po.RIX2 <- paste0("PO.", RIX2) 
po.RIX3 <- paste0("PO.", RIX3) 
po.RIX4 <- paste0("PO.", RIX4) 
po.RIX6 <- paste0("PO.", RIX6) 
po.RIX7 <- paste0("PO.", RIX7)
po.RIX8 <- paste0("PO.", RIX8) 
po.RIX9 <- paste0("PO.", RIX9) 
po.RIX10 <- paste0("PO.", RIX10)

rixFirst <- c(RIX1,RIX2,RIX3,RIX4,RIX6,RIX7,RIX8,RIX9,RIX10)
poRIXFirst <- c(po.RIX1,po.RIX2,po.RIX3,po.RIX4,po.RIX6,po.RIX7,po.RIX8,po.RIX9,po.RIX10)


##################
# Compare groups #
##################
diet.contr <- matrix(rbind(c(1,1,0,0), c(1,0,1,0), c(1,0,0,1), c(0,1,1,0), c(0,1,0,1),c(0,0,1,1)), nrow=6)
rix.contr <- matrix(rbind(c(1,1,0,0,0,0,0,0,0), c(1,0,1,0,0,0,0,0,0), c(1,0,0,1,0,0,0,0,0),
                          c(1,0,0,0,1,0,0,0,0), c(1,0,0,0,0,1,0,0,0), c(1,0,0,0,0,0,1,0,0),
                          c(1,0,0,0,0,0,0,1,0), c(1,0,0,0,0,0,0,0,1), 
                          c(0,1,1,0,0,0,0,0,0), c(0,1,0,1,0,0,0,0,0), c(0,1,0,0,1,0,0,0,0), 
                          c(0,1,0,0,0,1,0,0,0), c(0,1,0,0,0,0,1,0,0), c(0,1,0,0,0,0,0,1,0),
                          c(0,1,0,0,0,0,0,0,1), c(0,0,1,1,0,0,0,0,0), c(0,0,1,0,1,0,0,0,0), 
                          c(0,0,1,0,0,1,0,0,0), c(0,0,1,0,0,0,1,0,0), c(0,0,1,0,0,0,0,1,0),
                          c(0,0,1,0,0,0,0,0,1), c(0,0,0,1,1,0,0,0,0), c(0,0,0,1,0,1,0,0,0),
                          c(0,0,0,1,0,0,1,0,0), c(0,0,0,1,0,0,0,1,0), c(0,0,0,1,0,0,0,0,1),
                          c(0,0,0,0,1,1,0,0,0), c(0,0,0,0,1,0,1,0,0), c(0,0,0,0,1,0,0,1,0),
                          c(0,0,0,0,1,0,0,0,1), c(0,0,0,0,0,1,1,0,0), c(0,0,0,0,0,1,0,1,0),
                          c(0,0,0,0,0,1,0,0,1), c(0,0,0,0,0,0,1,1,0), c(0,0,0,0,0,0,1,0,1),
                          c(0,0,0,0,0,0,0,1,1)), nrow=36, ncol=9)  
for(i in 1:length(diet))
diet_diff <- list()
RIX_diff <- list()
for(i in 1:1){   #length(phenotype)
  all.reg <- list.reg[[i]]
  all.mcmc <- mcmc.stack(all.reg)
  
  compDiets <- rbind(dietlabs, RIX1,RIX2,RIX3,RIX4,RIX6,RIX7,RIX8,RIX9,RIX10,
                     po.RIX1,po.RIX2,po.RIX3,po.RIX4,po.RIX6,po.RIX7,po.RIX8,po.RIX9,po.RIX10)
  diet_diff[[phenotype[i]]] <- rixtest(all.mcmc, compDiets, diet.contr, dietlabs)
  
  compRIX <- rbind(rixlabs, PO.R) #, std,lowp,methyl,vitd,po.std,po.lowp,po.methyl,po.vitd)
  RIX_diff[[phenotype[i]]] <- rixtest(all.mcmc, compRIX, rix.contr, rixlabs)
}

signif_tables <- list()
for(i in 1:length(phenotype)){
  #signif_tables[[phenotype[i]]]$diet_diff <- xtable(diet_diff[[phenotype[i]]]$signif[1:10,], 
  #                                                  label=paste("Significant diet contrasts for", phenotype[i]))
  #signif_tables[[phenotype[i]]]$podiet_diff <- xtable(diet_diff[[phenotype[i]]]$signif[11:19,], 
  #                                                    label=paste("Significant PO-diet contrasts for", phenotype[i]))
  signif_tables[[phenotype[i]]]$RIX_diff1 <- xtable(RIX_diff[[phenotype[i]]]$signif[,1:6], 
                                                   label=paste("Significant RIX contrasts for", phenotype[i]))
  signif_tables[[phenotype[i]]]$RIX_diff2 <- xtable(RIX_diff[[phenotype[i]]]$signif[,7:12], 
                                                    label=paste("Significant RIX contrasts for", phenotype[i]))
  signif_tables[[phenotype[i]]]$RIX_diff3 <- xtable(RIX_diff[[phenotype[i]]]$signif[,13:18], 
                                                    label=paste("Significant RIX contrasts for", phenotype[i]))
  signif_tables[[phenotype[i]]]$RIX_diff4 <- xtable(RIX_diff[[phenotype[i]]]$signif[,19:24], 
                                                    label=paste("Significant RIX contrasts for", phenotype[i]))
  signif_tables[[phenotype[i]]]$RIX_diff5 <- xtable(RIX_diff[[phenotype[i]]]$signif[,25:30], 
                                                    label=paste("Significant RIX contrasts for", phenotype[i]))
  signif_tables[[phenotype[i]]]$RIX_diff6 <- xtable(RIX_diff[[phenotype[i]]]$signif[,31:36], 
                                                    label=paste("Significant RIX contrasts for", phenotype[i]))
}


tempdietrix <- c()
tempodietrix <- c()
tempdiets <- c()
tempdiet1 <- c()
tempdiet2 <- c()
tempdietrixcont <- c()
tempodietrixcont <- c()

#finNames <- matrix(NA, nrow=length(rixlabs), ncol=2+length(dietlabs)*2 + nrow(diet.contr)*2)
#for (i in 1:length(rixlabs)){  
#  temprix <- rixlabs[i]
#  tempo <- paste0("PO.",rixlabs[i])
#  for (j in 1:length(dietlabs)){
#    tempdietrix[j] <- paste0(dietlabs[j],".",rixlabs[i])
#    tempodietrix[j] <- paste0("PO.", dietlabs[j],".",rixlabs[i])
#  }
#  for (k in 1:nrow(diet.contr)){
#    tempdiets <- dietlabs[which(diet.contr[k,] == 1)]
#    tempdiet1 <- tempdietrix[grep(tempdiets[1],tempdietrix)]
#    tempdiet2 <- tempdietrix[grep(tempdiets[2],tempdietrix)]
#    tempdietrixcont[k] <- paste0(tempdiet1,tempdiet2)
#    tempodietrixcont[k] <- paste0("PO.",tempdiet1,"PO.",tempdiet2)
#  }
#  finNames[i,] <- c(temprix, tempo, tempdietrix[1:j], tempodietrix[1:j], 
#                    tempdietrixcont[1:k], tempodietrixcont[1:k])
#}             


group = c(rep(1, length(dietlabs)), rep(2,9), rep(3,9), rep(4,6))
group2 = c(1,1,rep(2,length(dietlabs)), rep(3,6), rep(4,4), rep(5,6))

for (k in 1:1){#length(phenotype)
  tryplots <- list()
  my.plots <- list()
  all.reg <- list.reg[[k]]
  all.mcmc <- mcmc.stack(all.reg)
  diffdiet_post <- do.call(cbind, diet_diff[[phenotype[k]]]$rixdiff)
  diffRIX_post <- do.call(cbind, RIX_diff[[phenotype[k]]]$rixdiff)
  all.dist <- cbind(all.mcmc, diffdiet_post, diffRIX_post)
  
  for (i in 1:length(rixlabs)){
    tryplots[[i]] <- plot.mine.hpd(all.dist, wanted=finNames[i,], addline=T, group=group2)
    if (i %% 3 == 0){
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(2, 3, heights = unit(c(0.5, 8),"null"))))  
      grid.text(paste("Estimates for",phenotype[k]), vp = viewport(layout.pos.row = 1, layout.pos.col = 1:3))
      print(tryplots[[i-2]], vp = viewport(layout.pos.row = 2, layout.pos.col = 1))         
      print(tryplots[[i-1]], vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
      print(tryplots[[i]], vp = viewport(layout.pos.row = 2, layout.pos.col = 3))
      my.plots[[i/3]] <- recordPlot()
    }
  }
  tryplots[[length(rixlabs)+1]] <- plot.mine.hpd(all.dist, addline=T,group=group,  
                                                 wanted=c(dietlabs,rixlabs, PO.R, 
                                                          names(diet_diff[["WeightPND21"]]$rixdiff[1:6])))
  tryplots[[length(rixlabs)+2]] <- plot.mine.hpd(all.dist, addline=F, 
                                                 wanted=c(names(RIX_diff[["WeightPND21"]]$rixdiff[1:36])))
      
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(2, 2, heights = unit(c(0.5, 8),"null"))))  
  grid.text(paste("Estimates for",phenotype[k]), vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2))
  print(tryplots[[length(rixlabs)+1]], vp = viewport(layout.pos.row = 2, layout.pos.col = 1))  
  print(tryplots[[length(rixlabs)+2]], vp = viewport(layout.pos.row = 2, layout.pos.col = 2))    
  
  my.plots[[(i/3)+1]] <- recordPlot()

  graphics.off()
  
  pdf(paste0(phenotype[k],'hpd_all.pdf'), onefile=TRUE, width = 14, height = 11)
  for (my.plot in my.plots) {
    replayPlot(my.plot)
  }
  graphics.off() 
}


##############
# Prediction #
##############

all.reg <- list.reg[[1]]
all.mcmc <- mcmc.stack(all.reg)
#diffdiet_post <- do.call(cbind, diet_diff[[phenotype[1]]]$rixdiff)
#diffRIX_post <- do.call(cbind, RIX_diff[[phenotype[1]]]$rixdiff)
all.dist <- cbind(all.mcmc, rep(0,nrow(all.mcmc)), rep(0,nrow(all.mcmc))) #, diffdiet_post, diffRIX_post)
colnames(all.dist)[c(length(colnames(all.dist))-1, length(colnames(all.dist)))] <- c("emp", "POemp")

dietrix <- matrix(c(as.character(encoded$Level[which(encoded$Variable == "DietRIX")]),"emp", "emp"), ncol=9, nrow=4)
dietrix[,9] <- c("emp","ME10", "PD9", "emp")

predPOneg <- list()
predPOpos <- list()

for(i in 1:nrow(all.dist)){
  for (rixn in 1:nrix){
    for (dietn in 1:length(dietlabs)){
      
      rixlabs <- as.character(encoded$Level[which(encoded$Variable == "RIX")])
      dietlabs <- as.character(encoded$Level[which(encoded$Variable == "Diet")])
      
      temprix <- rixlabs[rixn]
      tempdiet <- dietlabs[dietn]
      tempinter <- dietrix[dietn,rixn]
      prLab <- paste0(tempdiet,temprix)
      
      predPOneg[[paste0(prLab, ".pos")]][i] <- all.dist[i,temprix] + all.dist[i,tempdiet] + 
        all.dist[i,dietrix[dietn,rixn]] + all.dist[i,paste0("PO", temprix)]*(-0.5) + 
        all.dist[i,paste0("PO", tempinter)]*(-0.5)
      predPOpos[[paste0(prLab, ".neg")]][i] <- all.dist[i,temprix] + all.dist[i,tempdiet] + 
        all.dist[i,dietrix[dietn,rixn]] + all.dist[i,paste0("PO", temprix)]*(0.5) + 
        all.dist[i,paste0("PO", tempinter)]*(0.5)
    }
  }
}


allpredPOneg <- do.call(cbind, predPOneg)
allpredPOpos <- do.call(cbind, predPOpos)
allpredPO <- list(allpredPOneg, allpredPOpos) #, all.dist[,rixlabs]
#groupPred <- rep(1:length(rixlabs), each=length(dietlabs))
#predplots <- plot.mine.hpd(allpredPO, addline=T, group=groupPred, wanted= colnames(allpredPOneg), wide=F)


flydata <- as.data.frame(c(allpredPO))
flydata <- as.mcmc(flydata)
mu    <- c(colMeans(flydata))
med   <- apply(coda::HPDinterval(flydata, prob=0.01), 1, median)
hpd.wide    <- coda::HPDinterval(flydata, prob=0.95)
hpd.narrow  <- coda::HPDinterval(flydata, prob=0.5)
printdata <- data.frame(cbind(mu,med,hpd.wide,hpd.narrow))
printdata$dietrix <- do.call(rbind,strsplit(rownames(printdata),"[.]"))[,1] 
printdata$diet <- gsub('[0-9]+', '', printdata$dietrix)
printdata$rix <- gsub('[A-Z]+', '', printdata$dietrix)
printdata$poe <- printdata$dietrix <- do.call(rbind,strsplit(rownames(printdata),"[.]"))[,2]
printdata$poen <- ifelse(printdata$poe == "pos", 0.5, -0.5)
printdata$rixpo <- paste0(printdata$rix, printdata$poe)


#factor(rep(dietlabs,length.out=length(printdata$mu)), levels=dietmeans$Group.1)
#printdata$diet <- rep("",length(printdata$mu))
#printdata$diet[grep("STD",colnames(flydata))] <- "STD"
#printdata$diet[grep("ME",colnames(flydata))] <- "ME"
#printdata$diet[grep("PD",colnames(flydata))] <- "PD"
#printdata$diet[grep("VDD",colnames(flydata))] <- "VDD"
#printdata$rixes <- factor(rep(rixlabs, each=length(dietlabs)), levels=rixlabs)
#printdata$names <- factor(rep(colnames(allpredPOneg),2), levels = colnames(allpredPOneg))
#printdata$poe <- rep(c("+","-"), each=ncol(allpredPOneg))
#printdata$poen <- rep(c(0.5,-0.5), each=ncol(allpredPOneg))
#printdata$rixpo <- paste0(printdata$rixes, printdata$poe)
#printdata$rixpo <- factor(printdata$rixpo, 
#                          levels = rev(paste0(rep(rixlabs, each=2), 
#                                              rep(c("+","-"),length.out=length(rixlabs)*2))))
# rev(paste0(rep(rixlabs, each=2), rep(c("+","-"),length.out=length(rixlabs)*2)))
#printdata <- data.frame(printdata)[order(printdata$rixpo),]
printdata <- printdata[-which(rownames(printdata) %in% c("STD10.neg","VDD10.neg","STD10.pos","VDD10.pos")),]

dietmeans <- aggregate(printdata[,1:2], list(printdata$diet), mean)
dietmeans <- dietmeans[order(dietmeans$mu),]
#rixmeans <- aggregate(printdata[,1:2], list(printdata$rixes), mean)
#rixmeans <- rixmeans[order(rixmeans$mu),]

printdata$diet <- factor(printdata$diet, levels = dietmeans$Group.1)
printdata$rix <- factor(printdata$rix, levels = as.character(encoded$Level[which(encoded$Variable == "RIX")]))

melt.pred <- melt(as.matrix(flydata))
melt.pred$RIX <- rep(rixlabs, each=(4*2000))
melt.pred$diet <- rep(dietlabs, each=2000, length.out=length(melt.pred$RIX))
melt.pred$PO <- rep(c("+","-"), each=(4*2000*9))

###############
# Raw effects #
###############

all.reg <- list.reg[[1]]
all.mcmc <- mcmc.stack(all.mcmc)


##### start plots #####
#### for raw effects 
all.dist <- as.matrix(all.mcmc)

melt.rixes <- melt(all.dist[,which(colnames(all.dist) %in% listNames)])

melt.diet <- melt(all.dist[,1:4])
melt.diet$Var2 <- factor(melt.diet$Var2, levels=dietlabs[rank(effectdata$mu[1:4])])
ggplot(data.frame(melt.rixes), aes(x = Var2, y=value, fill=Var2)) +
  geom_violin() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_fill_discrete(guide=FALSE) + theme_minimal() + xlab("Variables") + ylab("Estimated effect size")
for(i in 1:nrow(melt.rixes)){
  if(i %in% grep("PO", melt.rixes$Var2)){
    melt.rixes[i,"RIX"] <- unlist(strsplit(as.character(melt.rixes$Var2[i]),"[.]"))[2]
    melt.rixes[i,"PO"] <- "PO"
  } else {
    melt.rixes[i,"RIX"] <- as.character(melt.rixes$Var2[i])
    melt.rixes[i,"PO"] <- "NA"
  }
}

melt.rixes$Var2 <- factor(melt.rixes$Var2, 
                          levels=rev(paste0(rep(c("PO.",""),length.out=length(rixlabs)*2),
                                            rep(rixlabs[rank(effectdata$mu[86:94])], each=2))))
melt.rixes$RIX <- factor(melt.rixes$RIX, levels=rev(rixlabs[rank(effectdata$mu[86:94])]))
ggplot(data.frame(melt.rixes), aes(x = Var2, y=value, fill=RIX, alpha=PO)) +
  geom_violin() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_fill_discrete(guide=FALSE) + theme_minimal() + xlab("RIX") + ylab("Estimated effect") +
  geom_hline(yintercept = 0, colour="lightsteelblue", size=0.5, linetype="longdash") +
  scale_alpha_discrete(range = rev(c(0.75,1)), guide=F)

beanplot(value ~ PO*RIX, data=melt.rixes, ll = 0.04,
         main = "Estimated effects of RIX and PO", side = "both", xlab="RIX",
         col = list("lightseagreen", c("paleturquoise", "gray30")),
         axes=F)
axis(1, at=c(1:9),labels=levels(melt.rixes$RIX))
axis(2)
legend("bottomright", fill = c("lightseagreen", "paleturquoise"),
       legend = c("No PO", "PO"), box.lty=0)

##### start plots #####
### for predicted effects 

melt.pred$RIX <- factor(melt.pred$RIX, levels=rixmeans$Group.1)
melt.pred$diet <- factor(melt.pred$diet, levels=dietmeans$Group.1)

ggplot(data.frame(melt.pred), aes(x = RIX, y=value, alpha=PO, fill=RIX)) +
  geom_violin() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_fill_discrete() + theme_minimal() + xlab("Diet") + ylab("Estimated effect") +
  geom_hline(yintercept = 0, colour="lightsteelblue", size=0.5, linetype="longdash") +
  scale_alpha_discrete(range = rev(c(0.6,1)), name="POE") +
  scale_color_discrete() + 
  facet_wrap(~diet) +
  ylab("Predicted effect") + xlab("") 


##### start multi bean plot

par(mfrow = c(2, 2))
for (i in 1:4){
  thisdiet = dietlabs[i]
  beanplot(value ~ PO*RIX, data=melt.pred[which(melt.pred$diet==thisdiet),], ll = 0.04,
           main = thisdiet, side = "both", xlab="",
           col = list("lightseagreen", c("paleturquoise", "gray30")),
           axes=F)
  axis(2)
  if (i>=3){
    axis(1, at=c(1:9),labels=levels(melt.pred$RIX))
    legend("bottomright", fill = c("lightseagreen", "paleturquoise"),
           legend = c("+", "-"), box.lty=0)
  }
}


beanplot(value ~ PO*RIX, data=melt.pred[which(melt.pred$diet=="ME"),], ll = 0.04,
         main = "Estimated effects of RIX and PO", side = "both", xlab="RIX",
         col = list("lightseagreen", c("paleturquoise", "gray30")),
         axes=F)
axis(1, at=c(1:9),labels=levels(melt.pred$RIX))
axis(2)
legend("bottomright", fill = c("lightseagreen", "paleturquoise"),
       legend = c("No PO", "PO"), box.lty=0)

beanplot(value ~ PO*RIX, data=melt.pred[which(melt.pred$diet=="PD"),], ll = 0.04,
         main = "Estimated effects of RIX and PO", side = "both", xlab="RIX",
         col = list("lightseagreen", c("paleturquoise", "gray30")),
         axes=F)
axis(1, at=c(1:9),labels=levels(melt.pred$RIX))
axis(2)
legend("bottomright", fill = c("lightseagreen", "paleturquoise"),
       legend = c("No PO", "PO"), box.lty=0)

beanplot(value ~ PO*RIX, data=melt.pred[which(melt.pred$diet=="VDD"),], ll = 0.04,
         main = "Estimated effects of RIX and PO", side = "both", xlab="RIX",
         col = list("lightseagreen", c("paleturquoise", "gray30")),
         axes=F)
axis(1, at=c(1:9),labels=levels(melt.pred$RIX))
axis(2)
legend("bottomright", fill = c("lightseagreen", "paleturquoise"),
       legend = c("No PO", "PO"), box.lty=0)




heatma <- ggplot(printdata, aes(rixpo, diet)) + 
  geom_tile(aes(fill = mu), colour = "honeydew") + 
  scale_fill_gradient(low = "honeydew", high = "aquamarine4")

ribbon9 <- ggplot(printdata, aes(color=rix)) + #geom_point(aes(y=mu, x=diet, alpha=poe)) + 
  theme_minimal() + 
  geom_line(aes(y=mu, x=diet, group=poe, alpha=poe)) +  
  geom_ribbon(color=NA, aes(x=diet, ymin=lower.1, ymax=upper.1, group=poe, alpha=poe, 
                  fill=rix, color=rix)) + 
  geom_hline(yintercept = 0, colour="lightsteelblue", size=0.5, linetype="longdash") +
  scale_alpha_discrete(range = rev(c(0.45, 0.75)), name="POE") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_fill_discrete(guide=FALSE) + scale_color_discrete(guide=FALSE) + 
  facet_wrap(~rix) +
  ylab("Predicted effect") + xlab("") 

ribbon4 <- ggplot(printdata, aes(color=diet)) + 
  geom_line(aes(y=mu, x=rix, group=poe, alpha=poe)) + 
  geom_ribbon(color=NA, aes(x=rix, ymin=lower.1, ymax=upper.1, group=poe, alpha=poe, 
                            fill=diet, color=diet)) + 
  geom_hline(yintercept = 0, colour="lightsteelblue", size=0.5, linetype="longdash") +
  scale_alpha_discrete(range = rev(c(0.45, 0.75)), name="POE") +
  scale_fill_discrete(guide=FALSE) + scale_color_discrete(guide=FALSE) + 
  facet_wrap(~diet) + theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylab("Predicted effect") + xlab("") 

boxes4col <- ggplot(printdata, aes(color=poe)) + ylim(-2.5, 2.5) +
  geom_segment(aes(y=poen, yend=poen, x=lower.1, xend=upper.1, group=poe), size=1.5) + 
  geom_segment(aes(y=poen, yend=poen, x=lower, xend=upper, group=poe), size=0.5) + 
  geom_point(aes(y=poen, x=med, group=poen), size=3) +
  geom_point(aes(y=poen, x=mu, group=poen), size=3, shape="|", color="white") +
  facet_grid(rixes~diet) +
  geom_vline(yintercept = 0, colour="lightsteelblue", size=0.5, linetype="longdash") +
  scale_alpha_discrete(range = rev(c(0.45, 0.75)), guide=FALSE) + 
  theme_minimal() + scale_color_brewer(palette="Paired",guide=FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  xlab("Predicted effect") + ylab("") 

boxes9col <- ggplot(printdata, aes(color=poe)) + ylim(-3, 3) +
  geom_segment(aes(y=poen, yend=poen, x=lower.1, xend=upper.1, group=poe), size=1.5) + 
  geom_segment(aes(y=poen, yend=poen, x=lower, xend=upper, group=poe), size=0.5) + 
  geom_point(aes(y=poen, x=med, group=poen), size=3) +
  geom_point(aes(y=poen, x=mu, group=poen), size=3, shape="|", color="white") +
  facet_grid(diet~rixes) +
  geom_vline(yintercept = 0, colour="lightsteelblue", size=0.5, linetype="longdash") +
  scale_alpha_discrete(range = rev(c(0.45, 0.75)), guide=FALSE) + 
  theme_minimal() + scale_color_brewer(palette="Paired",guide=FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  xlab("Predicted effect") + ylab("") 
  
  
graphics.off()

pdf(paste0(phenotype[1],'ribbon.pdf'), onefile=TRUE, width = 14, height = 11)
print(ribbon9)
print(ribbon4)
#print(boxes4col)
#print(boxes9col)
graphics.off() 



#+ (1 | RIX) + (PO - 1 | RIX) + (Diet - 1 | RIX) + (Diet:PO - 1 | RIX)")
##########
# Labels #
##########

allnames <- c()
for(i in 1:length(variables)){
  wantParam <- allparam[which(variables %in% allparam)]
  allnames <- c(allnames, as.character(encoded$Level[which(encoded$Variable == wantParam[i])])) 
}

for(i in 1:length(check.list.reg)){
  varnames(check.list.reg[[i]])[1:length(allnames)] <- allnames
}  


######################
# Loop through model #
######################
list.reg <- list()
for (i in 1:length(phenotype)){  
  use <- which(is.na(orderdat[,phenotype[i]]) == F)
  matnut_use <- orderdat[use,]
  pheno = matnut_use[,phenotype[i]]
  pheno = scale(pheno)
  
  N <- length(pheno)
  y <- pheno
  
  diet <- as.numeric(matnut_use$Diet)
  rix <- as.numeric(matnut_use$RIX)
  nrix <- length(which(encoded$Variable == "RIX")) 
  ndiet <- length(which(encoded$Variable == "Diet")) 
  poe <- matnut_use$PO
  dietrix <- as.numeric(matnut_use$DietRIX)
  ndietrix <- length(which(encoded$Variable == "DietRIX"))  
  
  data = list("N"=N, "rix"=rix, "y"=y, "nrix"=nrix, "diet"=diet, "ndiet"=ndiet, 
              "poe"=poe, "dietrix"=dietrix, "ndietrix"=ndietrix)
  
  #############
  # Run model #
  #############
  
  chains <- 2 
  reg.jags <- jags.model(textConnection(model), data=data, n.chains = chains, n.adapt = 2500)
  update(reg.jags, n.iter=2500)
  list.reg[[i]] <- coda.samples(reg.jags, variable.names = allparam, thin=10, n.iter=10000)
  
}

