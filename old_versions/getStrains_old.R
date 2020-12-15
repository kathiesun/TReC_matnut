

=======
  library(rjags)
library(coda)
library(ggplot2)
library(data.table)
library(lmerTest)
library("DESeq2")

# setwd("C:/Users/kathie/Dropbox/source/nutridiallel2014/src")
# setwd("C:/Users/kathie/Dropbox/source/nutridiallel2014/mnp2017/src")
# setwd("~/Dropbox/source/nutridiallel2014/src")

src <- ifelse(length(grep("C:", getwd())) >0, file.path("C:","Users","kathie","Documents","source_ks"), 
              ifelse(length(grep("kys6", getwd())) >0, file.path("~","source_ks/matnut_src"), file.path("~","Documents","matnut_src")))

source(file.path(src,"summary_functions.R"))
source(file.path(src,"parArgs.R"))
source(file.path(src,"lmer_functions_rna.R"))
#source(file.path(src,"prediction_functions.R"))
#source(file.path(src,"matching_functions3.R"))
#source("./mnp/simdata.R")
########## MIN CODE #############

cov_short <- readRDS(file.path("..","..","..","data","matnut_outputs",'cov_short_data.rds'))

# many of the following functions expect list of data with 1) data, 2) parameters, 3) phenotypes, 4) encoding

#summary(glm(vsn(ENSMUSG00000000001) ~ RIX + Diet + PORIX, data=orderdat, family = quasipoisson))
#length(which(colSums2(y.mat) == 0))

source(file.path("..","..","..","..","source_ks","matnut_src",'lmer_functions_rna.R')) 
# for killdevil


rnaseq <- list()
tryLam=c(-1, 0, .25, .33, .5, 1, 2, 3)

rixes <- cov_short$encoded$Level[cov_short$encoded$Variable == "RIX"]

rnaseq <- list()
for(j in 1:length(rixes)){
  datause <- cov_short$df[which(paste(cov_short$df$RIX) == rixes[j]), 1:100]
  lmerFit <- runLMERModels_cov(dataf=datause, tryLam=tryLam, phenotype=cov_short$ptypes[1:100],
                               fixvar=c("PO", "Diet"))
  lmerSum <- lmer.getDecodedSummary(lmerOb=lmerFit, lmerFit$formula, phenotypes=names(lmerFit$lmerobj))
  
  
  rnaseq[[j]] <- list(summary=lmerSum, lmerObj=lmerFit)     
}



#####################
matnut <- readRDS(file.path("..","..","..","data","matnut_outputs",'matnut_data.rds'))
behaviors <- read.csv("../../../data/matnut_outputs/2015-10_behavior_pups.csv")
stmap <-read.csv("../../../data/matnut_outputs/cc_strainNameMap.csv")
stmap[,"StrainNum"] <- gsub("[a-zA-Z]", "", stmap$Alias)
sep <- unlist(strsplit(as.character(stmap$StrainName),"/"))
stmap[,"simStrain"] <- sep[grep("[0-9]",sep)]
behaviors[,"DamStrain"] <- stmap[,"simStrain"][match(behaviors$Dam.Line, stmap$StrainNum)]
behaviors[,"SireStrain"] <- stmap[,"simStrain"][match(behaviors$Sire.Line, stmap$StrainNum)]

mapping <- behaviors[,c(5:7,10:11,13:14,30:31)]

colnames(mapping)[which(colnames(mapping) == "Pup.ID")] <- "ID"
cov_mapped <- merge(matnut$df, mapping[which(as.numeric(paste(mapping$ID)) %in% matnut$df$ID),], by="ID")
identical(as.character(cov_mapped$Cross.x), as.character(cov_mapped$Cross.y))

all.cc <- unique(cov_short_mapped$SireStrain.y)

cov_mapped[1:10,c(1,4,5,6,8,68:70,77:78)]
var1 <- read.csv(paste(src,"var_1286800.csv", sep="/"))
var1[,"strain_short"] <- unlist(strsplit(as.character(var1$strain), "_"))[grep("^CC",unlist(strsplit(as.character(var1$strain), "_")))]
var1_short <- var1[which(var1$variant_id == 1286800),]

setdiff(all.cc, var1_short$strain_short)

var1_complete <- matrix(NA, ncol=ncol(var1_short), nrow=length(all.cc))
colnames(var1_complete) <- colnames(var1_short)
var1_complete[1:nrow(var1_short), ] <- as.matrix(var1_short)
var1_complete[(nrow(var1_short)+1):nrow(var1_complete), "strain_short"]  <- setdiff(all.cc, var1_short$strain_short)


which_al <- ifelse(length(unique(var1_short$allele_1[which(var1_short$consequence_1 == "reference")])) > 0, "1", "2")
not_al <- ifelse(which_al == 1, 2, 1)

ref <- unique(var1_short[,paste0("allele_",which_al)][which(var1_short[,paste0("consequence_",which_al)] == "reference")])
mut <- unique(var1_short$allele_1[which(var1_short$consequence_1 != "reference")])

var1_complete[(nrow(var1_short)+1):nrow(var1_complete),paste0("allele_",which_al)] <- as.character(ref)
var1_complete[(nrow(var1_short)+1):nrow(var1_complete),paste0("allele_",not_al)] <- as.character(mut)



behaviors[,c("dam_allele_1","dam_allele_2",
             "sire_allele_1","sire_allele_2")] <- c(var1_complete[,"allele_1"][match(cov_short_mapped$DamStrain, var1_complete[,"strain_short"])],
                                                    var1_complete[,"allele_2"][match(cov_short_mapped$DamStrain, var1_complete[,"strain_short"])],
                                                    var1_complete[,"allele_1"][match(cov_short_mapped$SireStrain.y, var1_complete[,"strain_short"])],
                                                    var1_complete[,"allele_2"][match(cov_short_mapped$SireStrain.y, var1_complete[,"strain_short"])])


#################

counts <- fread(paste0(src,"/collated.kb.matrix.out"))
cov <- fread(paste0(src,"/covariates.txt"))
counts <- counts[,-1]
counts[1:10,1:10]
class(counts)
dim(counts)

#which(as.numeric(gsub("[A-Z]","",colnames(counts))) %in% cov$recordID)
y.mat <- data.frame(row.names = colnames(counts)[-1])
y.mat <- as.data.frame(t(counts[,-1]))
colnames(y.mat) <- counts$locus

cov_short <- data.frame(row.names=rownames(y.mat))
cov_short <- as.data.frame(cov[match(gsub("[A-Z]","",rownames(cov_short)), cov$recordID),])
rownames(cov_short) <- rownames(y.mat)
cov_short <- cov_short[-which(cov_short$Diet == ""),]

### covariates ###

#colnames(cov_short)

cov_short[grep("a", cov_short$dir),"PO"] <- 0.5
cov_short[grep("b", cov_short$dir),"PO"] <- -0.5
cov_short[which(cov_short$PO == 0.5), "pof"] <- "+"
cov_short[which(cov_short$PO == -0.5), "pof"] <- "-"

cov_short$Diet <- ifelse(cov_short$Diet == "Standard", "STD",
                         ifelse(cov_short$Diet == "Methyl Enriched", "ME",
                                ifelse(cov_short$Diet == "Low Protein", "PD",
                                       "VDD")))

dietlabs <- c("STD", "ME", "PD", "VDD")
cov_short$RIXdir <- cov_short$RIX
cov_short$RIX <- factor(cov_short$RRIX, levels=c(1:4,6:10))
cov_short$Diet <- factor(cov_short$Diet, levels=dietlabs)


cov_short$BBoriginal <- cov_short$Breeding.Batch
cov_short$BreedingBatch <- factor(paste0("br",cov_short$BBoriginal), levels=unique(paste0("br",cov_short$BBoriginal)))
cov_short$Cage <- factor(cov_short$Pup.Cage.num)
cov_short$DamID <- factor(paste0("d",cov_short$Dam.ID))
#cov_short$SireID <- factor(paste0("s",cov_short$Sire.ID))

orderdat <- cov_short[order(cov_short$RIX, cov_short$Diet, cov_short$PO),]
orderdat$PORIX <- paste0("PO",orderdat$RIX)
orderdat$PORIX <- factor(orderdat$PORIX, levels=unique(orderdat$PORIX))
orderdat$DietRIX <- paste0(orderdat$Diet, orderdat$RIX)
orderdat$DietRIX <- factor(orderdat$DietRIX, levels=unique(orderdat$DietRIX))
orderdat$PODietRIX <- paste0("PO",orderdat$DietRIX)
orderdat$PODietRIX <- factor(orderdat$PODietRIX, levels=unique(orderdat$PODietRIX))

variables <- c("Diet", "RIX", "PORIX", "DietRIX", "PODietRIX","DamID")
dispersions <- c("tauInv2","sigInv2", "tauPEInv2","tauDRInv2", "tauPDRInv2",
                 "tauBEInv2","tauDamInv2")
#"tauOFInv2","tauLDInv2","tauSIHInv2","tauFSTInv2","tauROInv2","tauREInv2")

encoded <- getEncoding(orderdat, variables)
#write.csv(encoded, file.path("./","cov_short_outputs/",'encoding.csv'), row.names = F)


#ptypes <- colnames(cov_short)[44:63]
#ptypes <- colnames(cov_short)[c(44,48,53,54,55,61,62,63)]
#ptypes <- c("OFTotalDistance","BasalCORT","StressCORT")

#vsted <- varianceStabilizingTransformation(as.matrix(orderdat[,grep("^ENSMUS", colnames(orderdat))]+1)[,1:50])
#ptypes <- colnames(vsted)[grep("ENSMUS",colnames(vsted))]
allparam <- c(variables, dispersions)
allparam <- allparam[order(allparam)]
y.mat_0 <- y.mat[,-which(colSums(as.matrix(y.mat)) == 0)]
merged <- cbind(orderdat, y.mat_0[match(rownames(orderdat), rownames(y.mat_0)),1:50])
#merged <- cbind(orderdat, vsted)
ptypes <- colnames(merged)[grep("ENSMUS",colnames(merged))]

cov_short <- list(df=merged, parameters=allparam, ptypes=ptypes, encoded=encoded)

#saveRDS(cov_short, file.path(".","matnut_outputs",'cov_short_data.rds'))
# many of the following functions expect list of data with 1) data, 2) parameters, 3) phenotypes, 4) encoding

#summary(glm(vsn(ENSMUSG00000000001) ~ RIX + Diet + PORIX, data=orderdat, family = quasipoisson))
#length(which(colSums2(y.mat) == 0))


rnaseq <- list()
tryLam=c(-1, 0, .25, .33, .5, 1, 2, 3)


#for(i in 1:length(cov_short$ptypes)){
#phenotype <- cov_short$ptypes[i]

#rnaseq[[cov_short$ptypes[i]]] <- makeSummary(datalist=cov_short, phenotype = cov_short$ptypes[i],
#                                            tryLam=trylam, normd=T,
#                                            chains=2, n.adapt=20000, n.iter=100000, thin=10, sq=F, addS="off")


rixes <- encoded$Level[encoded$Variable == "RIX"]


rnaseq <- list()
for(j in 1:length(rixes)){
  datause <- cov_short$df[which(paste(cov_short$df$RIX) == rixes[j]), ]
  #if(sum(datause[,phenotype]) > 0){
  lmerFit <- runLMERModels_cov(dataf=datause, tryLam=tryLam, phenotype=cov_short$ptypes,
                               fixvar=c("PO", "Diet"))
  lmerSum <- lmer.getDecodedSummary(lmerOb=lmerFit, lmerFit$formula, phenotypes=names(lmerFit$lmerobj))
  #lmerPlot <- lmerSum[which(lmerSum$Variable != "DamID"),]
  #lmerPlot$Level <- as.character(lmerPlot$Level)
  #lmerPlot$Level[which(lmerPlot$Variable == "PO")] <- "POE"
  
  #transf <- c("lambda" = lmerFit[[phenotype]]$phen_1$lambda,  
  #            "pval" = lmerFit[[phenotype]]$phen_1$best.pval)
  #lmerPlot <- transform(lmerPlot,Level=factor(Level,levels=unique(Level)), stat="identity")
  #catplot <- ggplot(data=lmerPlot, aes()) + 
  #  geom_segment(aes(x=Level, xend=Level, y = Intercept, yend=0), col="#000000", size=1) +
  #  geom_point(aes(x=Level, y=Intercept), col="red", size=3) +
  #  coord_flip() + theme_minimal() + 
  #  geom_hline(yintercept = 0, colour="gray50", size=0.5, linetype="longdash") +
  #  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  #  scale_alpha_continuous(range=c(0.5,1), guide=FALSE) +
  #  ggtitle(paste("LMER estimate for ", phenotype))
  rnaseq[[j]] <- list(summary=lmerSum, lmerObj=lmerFit)     #transform=transf, plot=catplot, 
}
print(paste(ptypes[i], "finished"))


#-1, 0, .25, .33, .5, 1, 2, 3

varianceStabilizingTransformation(as.matrix(y.mat+1))


##########
library(data.table)

df1 = data.table(PupID = LETTERS[1:10],
                 x = 1:10)

df2 = data.table(samplename = sample(replace = F, size=10, LETTERS[1:10]),
                 y = rnorm(n=10))

merged = df1[df2, on=c(PupID="samplename")]
print(merged)



df1 = data.table(PupID = LETTERS[1:10],
                 Diet = sample(size = 10, c("BEANS", "RICE"), replace = T),
                 x = 1:10)

df2 = data.table(samplename = sample(replace = F, size=10, LETTERS[1:10]),
                 Diet       = sample(size = 10, c("BEANS", "RICE"), replace = T),
                 y = rnorm(n=10))

merged = df1[df2, on=c(PupID="samplename", Diet="Diet")]

print(merged)


df2 = data.table(samplename = sample(replace = T, size=100000, LETTERS[1:10]),
                 y = rnorm(n=100000))

accumulated = df2[,list(mysum= sum(y)),by="samplename"]


################

library(foreach)
library(doSNOW)
data <- read.csv('dataset.csv')
cl <- makeCluster( mpi.universe.size(), type='MPI' )
clusterExport(cl,c('data'))
registerDoSNOW(cl)
results <- foreach( i = c(25,25,25,25) ) %dopar% {
  kmeans( x=data, centers=4, nstart=i )
}

temp.vector <- sapply( results, function(result) 
{ result$tot.withinss } )
result <- results[[which.min(temp.vector)]]

print(result)
stopCluster(cl)
mpi.exit()
=======
  library(rjags)
library(coda)
library(ggplot2)
library(data.table)
library(lmerTest)
library(car)
library(contrast)

setwd("~/matnut/src")
source("./matnut/summary_functions.R")
source("./matnut/parArgs.R")
source("./matnut/lmer_functions_rna.R")
source("./matnut/jags_functions.R")
source("./matnut/boxcox_functions.R")
dataSource <- file.path("..", "..","Dropbox","data","matnut_outputs")
cov_short <- readRDS(paste0(dataSource,'/cov_short_data.rds'))

variables <- c("Diet", "RIX", "PORIX", "DietRIX", "PODietRIX",
               "DamID","PO")

cov_short$df$PO <- as.factor(cov_short$df$PO)
cov_short$encoded <- getEncoding(cov_short$df, variables)

########## MIN CODE #############

#tryLam=c(-1, 0, .25, .33, .5, 1, 2, 3)
tryLam=c(0, 0.5)

rixes <- cov_short$encoded$Level[cov_short$encoded$Variable == "RIX"]
phenInd <- grep("ENS",colnames(cov_short$df))


system.time({
  rnaseq <- list()
  for(j in 1:3){    #length(rixes)){
    
    datause <- cov_short$df[which(paste(cov_short$df$RIX) == rixes[j]), ]
    fitLMER <- BC.model(y.mat=datause[,phenInd[1:50]], data=datause[,-phenInd], 
                        indvariable="~ 1+ PO + Diet", 
                        transformParams=getMatnutTransformParams(tryLam = c(0, 0.5), normd = T))
    fitLMER$formula <- "~ PO + Diet"
    lmerSum <- lmer.getDecodedSummary(lmerOb=fitLMER, formula="~ PO + Diet", phenotypes=fitLMER$phenotypes)
    
    jagsFit <- runJagsModels_cov(datalist=datause[,1:tail(which(colnames(datause) %in% fitLMER$phenotypes), n=1)], testLMER = fitLMER, 
                                 encoded=cov_short$encoded, phenotype=fitLMER$phenotypes, n.iter=200000)
    jagsSum <- list()
    for(i in 1:length(jagsFit)){
      all.reg <- jagsFit[[i]]$fit
      all.mcmc <- mcmc.stack(all.reg)
      
      jagsSum[[i]] <- jags.getDecodedSummary(all.mcmc, cov_short$encoded)
    }
    rnaseq[[j]] <- list(lmerSummary=lmerSum, lmerObj=fitLMER, jagsObj=jagsFit, jagsSummary=jagsSum)     
  }
})




for(i in 1:9){
  temp <- readRDS(paste0(dataSource,"/rna/rnaseq_out2_", rixes[i], ".RDS"))
  assign(paste0("rnaseq.",rixes[i]), temp)
}

mat <- list()
pdf("../../Data/matnut_outputs/lm_VDD_STD_compare.pdf")
for(i in 1:9){
  mat[[i]] <- matrix(NA, nrow=length(get(paste0("rnaseq.",rixes[i]))$lmerObj$fits), ncol=4)
  
  for(j in 1:nrow(mat[[i]])){
    mat[[i]][j,] <- summary(get(paste0("rnaseq.",rixes[i]))$lmerObj$fits[[j]])$coefficients["DietVDD",]
  }
  
  plot(mat[[i]][,1], -log10(mat[[i]][,4]), xlab="Effect size", ylab="-log10(p-value)",
       main=paste("Effect of VDD in RIX", rixes[i]))
}
dev.off()
#### reading in data from killdevil

for(i in 1:9){
  assign(paste0("sig.p",i), numeric(length(rnaseq.1)))
  name <- paste0("rnaseq.",i)
  for(j in 1:length(get(name))){
    temp <- ifelse(get(name)[[j]]$anova$`Pr(>F)`[1] < 0.05/length(rnaseq.1), 1, 0)
    name2 <- paste0("sig.p",i)
    assign(name2, c(get(name2),temp))
  }
}
# many of the following functions expect list of data with 1) data, 2) parameters, 3) phenotypes, 4) encoding

#summary(glm(vsn(ENSMUSG00000000001) ~ RIX + Diet + PORIX, data=orderdat, family = quasipoisson))
#length(which(colSums2(y.mat) == 0))

# for killdevil source(file.path("..","..","..","..","source_ks","matnut_src",'lmer_functions_rna.R')) 


cov_short <- readRDS(file.path("..","..","..","..","data","matnut_outputs",'cov_short_data.rds'))
for(i in 1:9){
  temp <- readRDS(paste0(src,"/data/","rnaseq.", i, ".RDS"))
  assign(paste0("rnaseq.",i), temp)
}

rnaseq.1$fit$lmerobj$ENSMUSG00000000001$anovaWrapper$an$`Pr(>F)`[[1]]

for(i in 1:9){
  assign(paste0("sig.p",i), numeric(length(rnaseq.1)))
  name <- paste0("rnaseq.",i)
  for(j in 1:length(get(name))){
    temp <- ifelse(get(name)[[j]]$anova$`Pr(>F)`[1] < 0.05/length(rnaseq.1), 1, 0)
    name2 <- paste0("sig.p",i)
    assign(name2, c(get(name2),temp))
  }
}
# many of the following functions expect list of data with 1) data, 2) parameters, 3) phenotypes, 4) encoding

#summary(glm(vsn(ENSMUSG00000000001) ~ RIX + Diet + PORIX, data=orderdat, family = quasipoisson))
#length(which(colSums2(y.mat) == 0))

# for killdevil source(file.path("..","..","..","..","source_ks","matnut_src",'lmer_functions_rna.R')) 



tryLam=c(-1, 0, .25, .33, .5, 1, 2, 3)

rixes <- cov_short$encoded$Level[cov_short$encoded$Variable == "RIX"]

rnaseq <- list()
for(j in 1:length(rixes)){
  datause <- cov_short$df[which(paste(cov_short$df$RIX) == rixes[j]), ]
  lmerFit <- runLMERModels_cov(dataf=datause, tryLam=0, phenotype=cov_short$ptypes,
                               fixvar=c("PO", "Diet"))
  lmerSum <- lmer.getDecodedSummary(lmerOb=lmerFit, lmerFit$formula, phenotypes=names(lmerFit$lmerobj))
  
  rnaseq[[j]] <- list(summary=lmerSum, lmerObj=lmerFit)     
}


saveRDS(rnaseq, "~/Dropbox/data/matnut_outputs/rnaseq_jags.RDS")
#####################

behaviors <- read.csv("../data/2015-10_behavior_pups.csv")
stmap <-read.csv("../data/cc_strainNameMap.csv")
stmap[,"StrainNum"] <- gsub("[a-zA-Z]", "", stmap$Alias)
sep <- unlist(strsplit(as.character(stmap$StrainName),"/"))
stmap[,"simStrain"] <- sep[grep("[0-9]",sep)]
behaviors[,"DamStrain"] <- stmap[,"simStrain"][match(behaviors$Dam.Line, stmap$StrainNum)]
behaviors[,"SireStrain"] <- stmap[,"simStrain"][match(behaviors$Sire.Line, stmap$StrainNum)]

mapping <- behaviors[,c(5:7,10:11,13:14,30:31)]

colnames(mapping)[which(colnames(mapping) == "Pup.ID")] <- "ID"
cov_short_mapped <- merge(cov_short$df, mapping[which(mapping$ID %in% cov_short$df$ID),], by="ID")
identical(as.character(test$Cross.x), as.character(test$Cross.y))

all.cc <- unique(cov_short_mapped$SireStrain.y)

cov_short_mapped[1:10,c(1,4,5,6,8,68:70,77:78)]
var1 <- read.csv(paste(src,"var_1286800.csv", sep="/"))
var1[,"strain_short"] <- unlist(strsplit(as.character(var1$strain), "_"))[grep("^CC",unlist(strsplit(as.character(var1$strain), "_")))]
var1_short <- var1[which(var1$variant_id == 1286800),]

setdiff(all.cc, var1_short$strain_short)

var1_complete <- matrix(NA, ncol=ncol(var1_short), nrow=length(all.cc))
colnames(var1_complete) <- colnames(var1_short)
var1_complete[1:nrow(var1_short), ] <- as.matrix(var1_short)
var1_complete[(nrow(var1_short)+1):nrow(var1_complete), "strain_short"]  <- setdiff(all.cc, var1_short$strain_short)


which_al <- ifelse(length(unique(var1_short$allele_1[which(var1_short$consequence_1 == "reference")])) > 0, "1", "2")
not_al <- ifelse(which_al == 1, 2, 1)

ref <- unique(var1_short[,paste0("allele_",which_al)][which(var1_short[,paste0("consequence_",which_al)] == "reference")])
mut <- unique(var1_short$allele_1[which(var1_short$consequence_1 != "reference")])

var1_complete[(nrow(var1_short)+1):nrow(var1_complete),paste0("allele_",which_al)] <- as.character(ref)
var1_complete[(nrow(var1_short)+1):nrow(var1_complete),paste0("allele_",not_al)] <- as.character(mut)



behaviors[,c("dam_allele_1","dam_allele_2",
             "sire_allele_1","sire_allele_2")] <- c(var1_complete[,"allele_1"][match(cov_short_mapped$DamStrain, var1_complete[,"strain_short"])],
                                                    var1_complete[,"allele_2"][match(cov_short_mapped$DamStrain, var1_complete[,"strain_short"])],
                                                    var1_complete[,"allele_1"][match(cov_short_mapped$SireStrain.y, var1_complete[,"strain_short"])],
                                                    var1_complete[,"allele_2"][match(cov_short_mapped$SireStrain.y, var1_complete[,"strain_short"])])


#################

counts <- fread(paste0(src,"/collated.kb.matrix.out"))
cov <- fread(paste0(src,"/covariates.txt"))
counts <- counts[,-1]
counts[1:10,1:10]
class(counts)
dim(counts)

#which(as.numeric(gsub("[A-Z]","",colnames(counts))) %in% cov$recordID)
y.mat <- data.frame(row.names = colnames(counts)[-1])
y.mat <- as.data.frame(t(counts[,-1]))
colnames(y.mat) <- counts$locus

cov_short <- data.frame(row.names=rownames(y.mat))
cov_short <- as.data.frame(cov[match(gsub("[A-Z]","",rownames(cov_short)), cov$recordID),])
rownames(cov_short) <- rownames(y.mat)
cov_short <- cov_short[-which(cov_short$Diet == ""),]

### covariates ###

#colnames(cov_short)

cov_short[grep("a", cov_short$dir),"PO"] <- 0.5
cov_short[grep("b", cov_short$dir),"PO"] <- -0.5
cov_short[which(cov_short$PO == 0.5), "pof"] <- "+"
cov_short[which(cov_short$PO == -0.5), "pof"] <- "-"

cov_short$Diet <- ifelse(cov_short$Diet == "Standard", "STD",
                         ifelse(cov_short$Diet == "Methyl Enriched", "ME",
                                ifelse(cov_short$Diet == "Low Protein", "PD",
                                       "VDD")))

dietlabs <- c("STD", "ME", "PD", "VDD")
cov_short$RIXdir <- cov_short$RIX
cov_short$RIX <- factor(cov_short$RRIX, levels=c(1:4,6:10))
cov_short$Diet <- factor(cov_short$Diet, levels=dietlabs)


cov_short$BBoriginal <- cov_short$Breeding.Batch
cov_short$BreedingBatch <- factor(paste0("br",cov_short$BBoriginal), levels=unique(paste0("br",cov_short$BBoriginal)))
cov_short$Cage <- factor(cov_short$Pup.Cage.num)
cov_short$DamID <- factor(paste0("d",cov_short$Dam.ID))
#cov_short$SireID <- factor(paste0("s",cov_short$Sire.ID))

orderdat <- cov_short[order(cov_short$RIX, cov_short$Diet, cov_short$PO),]
orderdat$PORIX <- paste0("PO",orderdat$RIX)
orderdat$PORIX <- factor(orderdat$PORIX, levels=unique(orderdat$PORIX))
orderdat$DietRIX <- paste0(orderdat$Diet, orderdat$RIX)
orderdat$DietRIX <- factor(orderdat$DietRIX, levels=unique(orderdat$DietRIX))
orderdat$PODietRIX <- paste0("PO",orderdat$DietRIX)
orderdat$PODietRIX <- factor(orderdat$PODietRIX, levels=unique(orderdat$PODietRIX))

variables <- c("Diet", "RIX", "PORIX", "DietRIX", "PODietRIX","DamID")
dispersions <- c("tauInv2","sigInv2", "tauPEInv2","tauDRInv2", "tauPDRInv2",
                 "tauBEInv2","tauDamInv2")
#"tauOFInv2","tauLDInv2","tauSIHInv2","tauFSTInv2","tauROInv2","tauREInv2")

encoded <- getEncoding(orderdat, variables)
#write.csv(encoded, file.path("./","cov_short_outputs/",'encoding.csv'), row.names = F)


#ptypes <- colnames(cov_short)[44:63]
#ptypes <- colnames(cov_short)[c(44,48,53,54,55,61,62,63)]
#ptypes <- c("OFTotalDistance","BasalCORT","StressCORT")

#vsted <- varianceStabilizingTransformation(as.matrix(orderdat[,grep("^ENSMUS", colnames(orderdat))]+1)[,1:50])
#ptypes <- colnames(vsted)[grep("ENSMUS",colnames(vsted))]
allparam <- c(variables, dispersions)
allparam <- allparam[order(allparam)]
y.mat_0 <- y.mat[,-which(colSums(as.matrix(y.mat)) == 0)]
merged <- cbind(orderdat, y.mat_0[match(rownames(orderdat), rownames(y.mat_0)),1:50])
#merged <- cbind(orderdat, vsted)
ptypes <- colnames(merged)[grep("ENSMUS",colnames(merged))]

cov_short <- list(df=merged, parameters=allparam, ptypes=ptypes, encoded=encoded)

#saveRDS(cov_short, file.path(".","matnut_outputs",'cov_short_data.rds'))
# many of the following functions expect list of data with 1) data, 2) parameters, 3) phenotypes, 4) encoding

#summary(glm(vsn(ENSMUSG00000000001) ~ RIX + Diet + PORIX, data=orderdat, family = quasipoisson))
#length(which(colSums2(y.mat) == 0))


rnaseq <- list()
tryLam=c(-1, 0, .25, .33, .5, 1, 2, 3)


#for(i in 1:length(cov_short$ptypes)){
#phenotype <- cov_short$ptypes[i]

#rnaseq[[cov_short$ptypes[i]]] <- makeSummary(datalist=cov_short, phenotype = cov_short$ptypes[i],
#                                            tryLam=trylam, normd=T,
#                                            chains=2, n.adapt=20000, n.iter=100000, thin=10, sq=F, addS="off")


rixes <- encoded$Level[encoded$Variable == "RIX"]


rnaseq <- list()
for(j in 1:length(rixes)){
  datause <- cov_short$df[which(paste(cov_short$df$RIX) == rixes[j]), ]
  #if(sum(datause[,phenotype]) > 0){
  lmerFit <- runLMERModels_cov(dataf=datause, tryLam=tryLam, phenotype=cov_short$ptypes,
                               fixvar=c("PO", "Diet"))
  lmerSum <- lmer.getDecodedSummary(lmerOb=lmerFit, lmerFit$formula, phenotypes=names(lmerFit$lmerobj))
  #lmerPlot <- lmerSum[which(lmerSum$Variable != "DamID"),]
  #lmerPlot$Level <- as.character(lmerPlot$Level)
  #lmerPlot$Level[which(lmerPlot$Variable == "PO")] <- "POE"
  
  #transf <- c("lambda" = lmerFit[[phenotype]]$phen_1$lambda,  
  #            "pval" = lmerFit[[phenotype]]$phen_1$best.pval)
  #lmerPlot <- transform(lmerPlot,Level=factor(Level,levels=unique(Level)), stat="identity")
  #catplot <- ggplot(data=lmerPlot, aes()) + 
  #  geom_segment(aes(x=Level, xend=Level, y = Intercept, yend=0), col="#000000", size=1) +
  #  geom_point(aes(x=Level, y=Intercept), col="red", size=3) +
  #  coord_flip() + theme_minimal() + 
  #  geom_hline(yintercept = 0, colour="gray50", size=0.5, linetype="longdash") +
  #  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  #  scale_alpha_continuous(range=c(0.5,1), guide=FALSE) +
  #  ggtitle(paste("LMER estimate for ", phenotype))
  rnaseq[[j]] <- list(summary=lmerSum, lmerObj=lmerFit)     #transform=transf, plot=catplot, 
}
print(paste(ptypes[i], "finished"))


#-1, 0, .25, .33, .5, 1, 2, 3

varianceStabilizingTransformation(as.matrix(y.mat+1))


##########
library(data.table)

df1 = data.table(PupID = LETTERS[1:10],
                 x = 1:10)

df2 = data.table(samplename = sample(replace = F, size=10, LETTERS[1:10]),
                 y = rnorm(n=10))

merged = df1[df2, on=c(PupID="samplename")]
print(merged)



df1 = data.table(PupID = LETTERS[1:10],
                 Diet = sample(size = 10, c("BEANS", "RICE"), replace = T),
                 x = 1:10)

df2 = data.table(samplename = sample(replace = F, size=10, LETTERS[1:10]),
                 Diet       = sample(size = 10, c("BEANS", "RICE"), replace = T),
                 y = rnorm(n=10))

merged = df1[df2, on=c(PupID="samplename", Diet="Diet")]

print(merged)


df2 = data.table(samplename = sample(replace = T, size=100000, LETTERS[1:10]),
                 y = rnorm(n=100000))

accumulated = df2[,list(mysum= sum(y)),by="samplename"]


=======
  library(rjags)
library(coda)
library(ggplot2)
library(data.table)
library(lmerTest)
library("DESeq2")

# setwd("C:/Users/kathie/Dropbox/source/nutridiallel2014/src")
# setwd("C:/Users/kathie/Dropbox/source/nutridiallel2014/mnp2017/src")
# setwd("~/Dropbox/source/nutridiallel2014/src")

src <- ifelse(length(grep("C:", getwd())) >0, file.path("C:","Users","kathie","Documents","source_ks"), 
              ifelse(length(grep("kys6", getwd())) >0, file.path("~","source_ks/matnut_src"), file.path("~","Documents","matnut_src")))

source(file.path(src,"summary_functions.R"))
source(file.path(src,"parArgs.R"))
source(file.path(src,"lmer_functions_rna.R"))
#source(file.path(src,"prediction_functions.R"))
#source(file.path(src,"matching_functions3.R"))
#source("./mnp/simdata.R")
########## MIN CODE #############

cov_short <- readRDS(file.path("..","..","..","data","matnut_outputs",'cov_short_data.rds'))

# many of the following functions expect list of data with 1) data, 2) parameters, 3) phenotypes, 4) encoding

#summary(glm(vsn(ENSMUSG00000000001) ~ RIX + Diet + PORIX, data=orderdat, family = quasipoisson))
#length(which(colSums2(y.mat) == 0))

source(file.path("..","..","..","..","source_ks","matnut_src",'lmer_functions_rna.R')) 
# for killdevil


rnaseq <- list()
tryLam=c(-1, 0, .25, .33, .5, 1, 2, 3)

rixes <- cov_short$encoded$Level[cov_short$encoded$Variable == "RIX"]

rnaseq <- list()
for(j in 1:length(rixes)){
  datause <- cov_short$df[which(paste(cov_short$df$RIX) == rixes[j]), 1:100]
  lmerFit <- runLMERModels_cov(dataf=datause, tryLam=tryLam, phenotype=cov_short$ptypes[1:100],
                               fixvar=c("PO", "Diet"))
  lmerSum <- lmer.getDecodedSummary(lmerOb=lmerFit, lmerFit$formula, phenotypes=names(lmerFit$lmerobj))
  
  
  rnaseq[[j]] <- list(summary=lmerSum, lmerObj=lmerFit)     
}



#####################
matnut <- readRDS(file.path("..","..","..","data","matnut_outputs",'matnut_data.rds'))
behaviors <- read.csv("../../../data/matnut_outputs/2015-10_behavior_pups.csv")
stmap <-read.csv("../../../data/matnut_outputs/cc_strainNameMap.csv")
stmap[,"StrainNum"] <- gsub("[a-zA-Z]", "", stmap$Alias)
sep <- unlist(strsplit(as.character(stmap$StrainName),"/"))
stmap[,"simStrain"] <- sep[grep("[0-9]",sep)]
behaviors[,"DamStrain"] <- stmap[,"simStrain"][match(behaviors$Dam.Line, stmap$StrainNum)]
behaviors[,"SireStrain"] <- stmap[,"simStrain"][match(behaviors$Sire.Line, stmap$StrainNum)]

mapping <- behaviors[,c(5:7,10:11,13:14,30:31)]

colnames(mapping)[which(colnames(mapping) == "Pup.ID")] <- "ID"
cov_mapped <- merge(matnut$df, mapping[which(as.numeric(paste(mapping$ID)) %in% matnut$df$ID),], by="ID")
identical(as.character(cov_mapped$Cross.x), as.character(cov_mapped$Cross.y))

all.cc <- unique(cov_short_mapped$SireStrain.y)

cov_mapped[1:10,c(1,4,5,6,8,68:70,77:78)]
var1 <- read.csv(paste(src,"var_1286800.csv", sep="/"))
var1[,"strain_short"] <- unlist(strsplit(as.character(var1$strain), "_"))[grep("^CC",unlist(strsplit(as.character(var1$strain), "_")))]
var1_short <- var1[which(var1$variant_id == 1286800),]

setdiff(all.cc, var1_short$strain_short)

var1_complete <- matrix(NA, ncol=ncol(var1_short), nrow=length(all.cc))
colnames(var1_complete) <- colnames(var1_short)
var1_complete[1:nrow(var1_short), ] <- as.matrix(var1_short)
var1_complete[(nrow(var1_short)+1):nrow(var1_complete), "strain_short"]  <- setdiff(all.cc, var1_short$strain_short)


which_al <- ifelse(length(unique(var1_short$allele_1[which(var1_short$consequence_1 == "reference")])) > 0, "1", "2")
not_al <- ifelse(which_al == 1, 2, 1)

ref <- unique(var1_short[,paste0("allele_",which_al)][which(var1_short[,paste0("consequence_",which_al)] == "reference")])
mut <- unique(var1_short$allele_1[which(var1_short$consequence_1 != "reference")])

var1_complete[(nrow(var1_short)+1):nrow(var1_complete),paste0("allele_",which_al)] <- as.character(ref)
var1_complete[(nrow(var1_short)+1):nrow(var1_complete),paste0("allele_",not_al)] <- as.character(mut)



behaviors[,c("dam_allele_1","dam_allele_2",
             "sire_allele_1","sire_allele_2")] <- c(var1_complete[,"allele_1"][match(cov_short_mapped$DamStrain, var1_complete[,"strain_short"])],
                                                    var1_complete[,"allele_2"][match(cov_short_mapped$DamStrain, var1_complete[,"strain_short"])],
                                                    var1_complete[,"allele_1"][match(cov_short_mapped$SireStrain.y, var1_complete[,"strain_short"])],
                                                    var1_complete[,"allele_2"][match(cov_short_mapped$SireStrain.y, var1_complete[,"strain_short"])])


#################

counts <- fread(paste0(src,"/collated.kb.matrix.out"))
cov <- fread(paste0(src,"/covariates.txt"))
counts <- counts[,-1]
counts[1:10,1:10]
class(counts)
dim(counts)

#which(as.numeric(gsub("[A-Z]","",colnames(counts))) %in% cov$recordID)
y.mat <- data.frame(row.names = colnames(counts)[-1])
y.mat <- as.data.frame(t(counts[,-1]))
colnames(y.mat) <- counts$locus

cov_short <- data.frame(row.names=rownames(y.mat))
cov_short <- as.data.frame(cov[match(gsub("[A-Z]","",rownames(cov_short)), cov$recordID),])
rownames(cov_short) <- rownames(y.mat)
cov_short <- cov_short[-which(cov_short$Diet == ""),]

### covariates ###

#colnames(cov_short)

cov_short[grep("a", cov_short$dir),"PO"] <- 0.5
cov_short[grep("b", cov_short$dir),"PO"] <- -0.5
cov_short[which(cov_short$PO == 0.5), "pof"] <- "+"
cov_short[which(cov_short$PO == -0.5), "pof"] <- "-"

cov_short$Diet <- ifelse(cov_short$Diet == "Standard", "STD",
                         ifelse(cov_short$Diet == "Methyl Enriched", "ME",
                                ifelse(cov_short$Diet == "Low Protein", "PD",
                                       "VDD")))

dietlabs <- c("STD", "ME", "PD", "VDD")
cov_short$RIXdir <- cov_short$RIX
cov_short$RIX <- factor(cov_short$RRIX, levels=c(1:4,6:10))
cov_short$Diet <- factor(cov_short$Diet, levels=dietlabs)


cov_short$BBoriginal <- cov_short$Breeding.Batch
cov_short$BreedingBatch <- factor(paste0("br",cov_short$BBoriginal), levels=unique(paste0("br",cov_short$BBoriginal)))
cov_short$Cage <- factor(cov_short$Pup.Cage.num)
cov_short$DamID <- factor(paste0("d",cov_short$Dam.ID))
#cov_short$SireID <- factor(paste0("s",cov_short$Sire.ID))

orderdat <- cov_short[order(cov_short$RIX, cov_short$Diet, cov_short$PO),]
orderdat$PORIX <- paste0("PO",orderdat$RIX)
orderdat$PORIX <- factor(orderdat$PORIX, levels=unique(orderdat$PORIX))
orderdat$DietRIX <- paste0(orderdat$Diet, orderdat$RIX)
orderdat$DietRIX <- factor(orderdat$DietRIX, levels=unique(orderdat$DietRIX))
orderdat$PODietRIX <- paste0("PO",orderdat$DietRIX)
orderdat$PODietRIX <- factor(orderdat$PODietRIX, levels=unique(orderdat$PODietRIX))

variables <- c("Diet", "RIX", "PORIX", "DietRIX", "PODietRIX","DamID")
dispersions <- c("tauInv2","sigInv2", "tauPEInv2","tauDRInv2", "tauPDRInv2",
                 "tauBEInv2","tauDamInv2")
#"tauOFInv2","tauLDInv2","tauSIHInv2","tauFSTInv2","tauROInv2","tauREInv2")

encoded <- getEncoding(orderdat, variables)
#write.csv(encoded, file.path("./","cov_short_outputs/",'encoding.csv'), row.names = F)


#ptypes <- colnames(cov_short)[44:63]
#ptypes <- colnames(cov_short)[c(44,48,53,54,55,61,62,63)]
#ptypes <- c("OFTotalDistance","BasalCORT","StressCORT")

#vsted <- varianceStabilizingTransformation(as.matrix(orderdat[,grep("^ENSMUS", colnames(orderdat))]+1)[,1:50])
#ptypes <- colnames(vsted)[grep("ENSMUS",colnames(vsted))]
allparam <- c(variables, dispersions)
allparam <- allparam[order(allparam)]
y.mat_0 <- y.mat[,-which(colSums(as.matrix(y.mat)) == 0)]
merged <- cbind(orderdat, y.mat_0[match(rownames(orderdat), rownames(y.mat_0)),1:50])
#merged <- cbind(orderdat, vsted)
ptypes <- colnames(merged)[grep("ENSMUS",colnames(merged))]

cov_short <- list(df=merged, parameters=allparam, ptypes=ptypes, encoded=encoded)

#saveRDS(cov_short, file.path(".","matnut_outputs",'cov_short_data.rds'))
# many of the following functions expect list of data with 1) data, 2) parameters, 3) phenotypes, 4) encoding

#summary(glm(vsn(ENSMUSG00000000001) ~ RIX + Diet + PORIX, data=orderdat, family = quasipoisson))
#length(which(colSums2(y.mat) == 0))


rnaseq <- list()
tryLam=c(-1, 0, .25, .33, .5, 1, 2, 3)


#for(i in 1:length(cov_short$ptypes)){
#phenotype <- cov_short$ptypes[i]

#rnaseq[[cov_short$ptypes[i]]] <- makeSummary(datalist=cov_short, phenotype = cov_short$ptypes[i],
#                                            tryLam=trylam, normd=T,
#                                            chains=2, n.adapt=20000, n.iter=100000, thin=10, sq=F, addS="off")


rixes <- encoded$Level[encoded$Variable == "RIX"]


rnaseq <- list()
for(j in 1:length(rixes)){
  datause <- cov_short$df[which(paste(cov_short$df$RIX) == rixes[j]), ]
  #if(sum(datause[,phenotype]) > 0){
  lmerFit <- runLMERModels_cov(dataf=datause, tryLam=tryLam, phenotype=cov_short$ptypes,
                               fixvar=c("PO", "Diet"))
  lmerSum <- lmer.getDecodedSummary(lmerOb=lmerFit, lmerFit$formula, phenotypes=names(lmerFit$lmerobj))
  #lmerPlot <- lmerSum[which(lmerSum$Variable != "DamID"),]
  #lmerPlot$Level <- as.character(lmerPlot$Level)
  #lmerPlot$Level[which(lmerPlot$Variable == "PO")] <- "POE"
  
  #transf <- c("lambda" = lmerFit[[phenotype]]$phen_1$lambda,  
  #            "pval" = lmerFit[[phenotype]]$phen_1$best.pval)
  #lmerPlot <- transform(lmerPlot,Level=factor(Level,levels=unique(Level)), stat="identity")
  #catplot <- ggplot(data=lmerPlot, aes()) + 
  #  geom_segment(aes(x=Level, xend=Level, y = Intercept, yend=0), col="#000000", size=1) +
  #  geom_point(aes(x=Level, y=Intercept), col="red", size=3) +
  #  coord_flip() + theme_minimal() + 
  #  geom_hline(yintercept = 0, colour="gray50", size=0.5, linetype="longdash") +
  #  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  #  scale_alpha_continuous(range=c(0.5,1), guide=FALSE) +
  #  ggtitle(paste("LMER estimate for ", phenotype))
  rnaseq[[j]] <- list(summary=lmerSum, lmerObj=lmerFit)     #transform=transf, plot=catplot, 
}
print(paste(ptypes[i], "finished"))


#-1, 0, .25, .33, .5, 1, 2, 3

varianceStabilizingTransformation(as.matrix(y.mat+1))


##########
library(data.table)

df1 = data.table(PupID = LETTERS[1:10],
                 x = 1:10)

df2 = data.table(samplename = sample(replace = F, size=10, LETTERS[1:10]),
                 y = rnorm(n=10))

merged = df1[df2, on=c(PupID="samplename")]
print(merged)



df1 = data.table(PupID = LETTERS[1:10],
                 Diet = sample(size = 10, c("BEANS", "RICE"), replace = T),
                 x = 1:10)

df2 = data.table(samplename = sample(replace = F, size=10, LETTERS[1:10]),
                 Diet       = sample(size = 10, c("BEANS", "RICE"), replace = T),
                 y = rnorm(n=10))

merged = df1[df2, on=c(PupID="samplename", Diet="Diet")]

print(merged)


df2 = data.table(samplename = sample(replace = T, size=100000, LETTERS[1:10]),
                 y = rnorm(n=100000))

accumulated = df2[,list(mysum= sum(y)),by="samplename"]


################

library(foreach)
library(doSNOW)
data <- read.csv('dataset.csv')
cl <- makeCluster( mpi.universe.size(), type='MPI' )
clusterExport(cl,c('data'))
registerDoSNOW(cl)
results <- foreach( i = c(25,25,25,25) ) %dopar% {
  kmeans( x=data, centers=4, nstart=i )
}

temp.vector <- sapply( results, function(result) 
{ result$tot.withinss } )
result <- results[[which.min(temp.vector)]]

print(result)
stopCluster(cl)
mpi.exit()