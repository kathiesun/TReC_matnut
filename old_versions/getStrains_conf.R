library(rjags)
library(coda)
library(ggplot2)
library(data.table)
library(lmerTest)
library(car)
library(contrast)
library(limma)
library(DESeq2)
library(metaRNASeq)
library(metafor)
library(metap)
library(refGenome)


setwd("~/matnut/src")
source("./matnut/summary_functions.R")
source("./matnut/parArgs.R")
source("./matnut/lmer_functions_rna.R")
source("./matnut/jags_functions.R")
source("./matnut/boxcox_functions.R")
dataSource <- file.path("..","..","..", "Dropbox\ (Valdarlab)","outputs","matnut_outputs")
#dataSource <- file.path("~/Data/matnut_outputs")
cov_short <- readRDS(paste0(dataSource,'/cov_short_data.rds'))
variables <- c("Diet", "RIX", "PORIX", "DietRIX", "PODietRIX","PO")

cov_short$encoded <- getEncoding(cov_short$df, variables)
rixes <- cov_short$encoded$Level[cov_short$encoded$Variable == "RIX"]
diets <- cov_short$encoded$Level[cov_short$encoded$Variable == "Diet"]
dietCombos <- matrix(NA, nrow=6, ncol=2)
n <- length(diets)
add=0
counter=1
for(i in 1:n){
  add = 1
  while(i+add <= n){
    dietCombos[counter,] <- c(paste(diets[i]), paste(diets[i+add]))
    add=add+1
    counter=counter+1
  }
}
phenInd <- grep("ENS",colnames(cov_short$df))

### DESeq ###
#############
assay <- data.matrix(t(cov_short$df[,phenInd]))
assay <- assay[which(rowMeans(assay) > 0),]
cov_short$df$diet_dir <- as.factor(paste0(cov_short$df$Diet, cov_short$df$dir))
cov_short$df$dir <- as.factor(cov_short$df$dir)
colData <- data.frame(cov_short$df[,-phenInd], stringsAsFactors = T)
colData$BreedingBatch
de <- DESeqDataSetFromMatrix(assay, colData, design= ~ dir + Diet + BreedingBatch)

require(sva)
mod  <- model.matrix(~ dir + Diet, colData)
mod0 <- model.matrix(~ 1, colData)
svseq <- svaseq(assay, mod, mod0, n.sv = 2)

desva <- de
desva$SV1 <- svseq$sv[,1]
desva$SV2 <- svseq$sv[,2]
design(desva) <- ~ SV1 + SV2 + dir + Diet
de <- desva

#vsd <- vst(de, blind = FALSE)  
#meanSdPlot(assay(vsd), ranks = FALSE)

for(i in 1:length(rixes)){
  de_temp <- de[, de$RIX == paste(rixes[i])]
  de_temp$Diet <- droplevels(colData(de_temp)$Diet)  
  #de_temp$BreedingBatch <- droplevels(colData(de_temp)$BreedingBatch)  
  dds <- DESeq(de_temp)
  res <- results( dds, contrast = c("dir", "b", "a") )
  res$gene <- rownames(res)
  res$variable <- rep("dir")
  res$contrast <- rep("b-a")
  res$rix <- rep(as.numeric(paste(rixes[i])))
  
  keep <- res[which(res$padj < 0.05),]
  
  if (i==1) {
    sigResults <- data.frame(keep)
    allResults <- data.frame(res)
  } else {
    if(nrow(keep) > 0) sigResults <- rbind(sigResults, keep)
    allResults <- rbind(allResults, res)
  }
  
  for(j in 1:nrow(dietCombos)){
    if(all(c(dietCombos[j,2], dietCombos[j,1]) %in% unique(de_temp$Diet)) ){
      res <- results( dds, contrast = c("Diet", dietCombos[j,2], dietCombos[j,1]) )
      res$gene <- rownames(res)
      res$variable <- rep("diet")
      res$contrast <- rep(paste(dietCombos[j,2], dietCombos[j,1], sep="-"))
      res$rix <- rep(as.numeric(paste(rixes[i])))
      allResults <- rbind(allResults, res)
      
      keep <- res[which(res$padj < 0.05),]
      
      if (nrow(keep) > 0) sigResults <- rbind(sigResults, keep)
    }
  }
  if(i == length(rixes)){
    rownames(sigResults) <- c()
    rownames(allResults) <- c()
  }
}

saveRDS(sigResults, paste0(dataSource,'/rna/deseq_sig_results_sva.rds'))
saveRDS(allResults, paste0(dataSource,'/rna/deseq_all_results_sva.rds'))


#####################
##  Meta-analyses  ##
#####################

allResults <- readRDS(paste0(dataSource,'/rna/deseq_all_results_sva.rds'))

get.pvals <- function(allResults, contrast=c("dir", "diet")){
  if(contrast == "dir") cont = "b-a"
  else {
    cont = c("ME-STD","PD-STD","VDD-STD","PD-ME","VDD-ME","VDD-PD")
  }
  rixes <- unique(allResults$rix)
  pval <- list()
  efsz <- list()
  mat.rawpval <- list()
  mat.efsz <- list()
  mat.efse <- list()
  
  for(k in 1:length(cont)){
    pval[[cont[k]]] <- list()
    efsz[[cont[k]]] <- list()
    mat.rawpval[[cont[k]]] <- matrix(NA, nrow = length(unique(allResults$gene)), ncol = length(rixes))
    rownames(mat.rawpval[[cont[k]]]) = unique(allResults$gene)
    mat.efsz[[cont[k]]] <- matrix(NA, nrow = length(unique(allResults$gene)), ncol = length(rixes))
    rownames(mat.efsz[[cont[k]]]) = unique(allResults$gene)
    mat.efse[[cont[k]]] <- matrix(NA, nrow = length(unique(allResults$gene)), ncol = length(rixes))
    rownames(mat.efse[[cont[k]]]) = unique(allResults$gene)
    
    for(j in 1:length(rixes)){
      subset <- allResults[intersect(which(allResults$rix == rixes[j]), which(allResults$contrast == cont[k])),]
      efsz[[cont[k]]][[j]] <- cbind(subset$log2FoldChange, subset$lfcSE)
      names(efsz[[cont[k]]][[j]]) <- subset$gene
      pval[[cont[k]]][[j]] <- cbind(subset$pvalue, subset$padj)
      names(pval[[cont[k]]][[j]]) <- subset$gene
      mat.rawpval[[cont[k]]][match(subset$gene, rownames(mat.rawpval[[cont[k]]])), j] <- subset$pvalue
      mat.efsz[[cont[k]]][match(subset$gene, rownames(mat.efsz[[cont[k]]])), j] <- subset$log2FoldChange
      mat.efse[[cont[k]]][match(subset$gene, rownames(mat.efse[[cont[k]]])), j] <- subset$lfcSE
    }
  }
  return(list(efsz = efsz, pval = pval, mat.rawpval = mat.rawpval, mat.efsz = mat.efsz, mat.efse = mat.efse))
  #efsz.vec <- t(sapply(1:length(efsz), function(x) efsz[[x]][,1]))
  #mat.rawpval <- do.call("cbind", rawpval)
}

dir <- get.pvals(allResults, contrast="dir")
diets <- get.pvals(allResults, contrast="diet")
  

  
# Using metaRNASeq
## dir ##

tmp_pvals <- lapply(apply(dir$mat.rawpval[[1]],2,as.list), unlist)
#adj_pvals <- lapply(tmp_pvals, function(x) lapply(x, function(y) y/Finv))
fishcomb <- fishercomb(tmp_pvals, BHth = 0.05)
par(mfrow=c(1,1))
hist(fishcomb$rawpval, breaks=100, col="grey", main="Fisher method", xlab = "Raw p-values (meta-analysis)")
dir$mat.rawpval[[1]][fishcomb$DEindices,]
  ## all 55 fdr-sig results
dir$mat.rawpval[[1]][order(fishcomb$adjpval)[1:50],]
  ## top 55 results ranked by combined adj p-value
  ## first 17 p-vals all "0"

par(mfrow=c(3,3))
sig <- c()
for(i in 1:9){
  hist(dir$pval$`b-a`[[i]][,1], main = paste("Unadj p-val of RIX ", rixes[i]))
  #sig <- c(sig, length(which(dir$pval$`b-a`[[i]][,2] < 0.1)))
}


# Using metap
## dir ##
schweder(dir$mat.rawpval[[1]][sort.list(fishcomb$rawpval)[2],], drawline = c("bh", "ls", "ab"),
         ls.control = list(frac = 0.5), ab.control = list(a = 0, b = 0.01))
    ## 6, 12, 14 look the best?

fishMetap <- sapply(1:nrow(dir$mat.rawpval[[1]]), function(x) {tmp <- dir$mat.rawpval[[1]][x,!is.na(dir$mat.rawpval[[1]][x,])]
                                                      if(length(tmp) > 2) sumlog(tmp)$p 
                                                      else 1}  )

varvec <- apply(dir$mat.rawpval[[1]],2,function(x) var(x, na.rm = T))
weights <- sapply(1:length(varvec), function(x) (sum(varvec^-1)^-1)/varvec[x])

keepNA <- apply(dir$mat.rawpval[[1]], 1, function(x) sum(is.na(x)) <= 3)
full.rawpval <- dir$mat.rawpval[[1]][keepNA,]

stoufcomb <- apply(full.rawpval, 1, function(x) sumz(x[!is.na(x)], weights=weights[!is.na(x)])$p)
adj_stouf <- p.adjust(stoufcomb, method = "BH")
sigGenes <- length(which(adj_stouf < 0.05))
length(which(rownames(full.rawpval)[sort.list(adj_stouf)[1:sigGenes]] %in% rownames(full.rawpval)[fishcomb$DEindices]) )
#length(which(rownames(full.rawpval)[sort.list(fishMetap)[1:55]] %in% rownames(full.rawpval)[fishcomb$DEindices]) )
full.rawpval[sort.list(stoufcomb)[1:10],]

row <- full.rawpval[sort.list(adj_stouf)[3],]
schweder(row[!is.na(row)], drawline = c("bh", "ls", "ab"),
         ls.control = list(frac = 0.5), ab.control = list(a = 0, b = 0.01))

## same as 14 above!!


#simulated data
N=10000
a=b=1
psim <- list()
psim <- sapply(1:9, function(x) rbeta(N, 1,1), simplify = F)
psim_2 <- sapply(1:9, function(x) rbeta(N, 1.2,1), simplify = F)
psim_3 <- sapply(1:9, function(x) rbeta(N, (abs(rnorm(1,0, 0.36)) + 1),1), simplify = F)

par(mfrow=c(3,3))
lapply(psim_3, function(x) hist(x))


fishcomb <- fishercomb(psim_3, BHth = 0.05)
par(mfrow=c(1,1))
hist(fishcomb$rawpval, breaks=100, col="grey", main="Fisher method", 
     xlab = "Raw p-values for simulated data")


## applying correction
vec <- dir$pval$`b-a`[[i]][,1]
vec <- fishcomb$rawpval



Finv <- unlist(lapply(psim_3[[1]], function(y) qchisq(y, 1)))
med <- median.default(Finv, na.rm=T)/qchisq(0.5, 1)
#adj_pvals <- unlist(lapply(unlist(lapply(vec, function(x) qnorm(x, 1)/med)), function(y) pnorm(y)) )
adj_pvals <- unlist(lapply(Finv, function(y) pchisq(y/med, 1)))
par(mfrow=c(3,3))
lapply(adj_pvals, function(x) hist(x))
plot(dchisq(seq(0,1,by=0.1), 1), type="l")

Finv <- lapply(psim_3, function(x) unlist(lapply(x, function(y) qchisq(y, 1))))
med <- unlist(lapply(Finv, function(x) median.default(x, na.rm=T)/qchisq(0.5, 1)))
#adj_pvals <- unlist(lapply(unlist(lapply(vec, function(x) qnorm(x, 1)/med)), function(y) pnorm(y)) )
adj_pvals <- sapply(1:9, function(x) unlist(lapply(Finv[[x]], function(y) pchisq(y/med[[x]], 1))), simplify=F)
par(mfrow=c(3,3))
lapply(adj_pvals, function(x) hist(x))

dir$mat.rawpval[[1]][fishcomb$DEindices,]
dir$mat.rawpval[[1]][order(fishcomb$adjpval)[1:55],]


########## Annotate

sig_genes <- rownames(full.rawpval)[sort.list(adj_stouf)[1:sigGenes]]

sig_genes_mat <- as.data.frame(full.rawpval[sort.list(stoufcomb)[1:sigGenes], ])
sig_genes_mat$gene_id = sig_genes

# create ensemblGenome object for storing Ensembl genomic annotation data
ens <- ensemblGenome()

# read GTF file into ensemblGenome object
read.gtf(ens, paste0(dataSource, "/Mus_musculus.GRCm38.91.gtf"))
all_genes <- extractSeqids(ens, ensPrimAssembly())
my_gene <- getGenePositions(all_genes)
sig_genes_temp <- my_gene[which(my_gene$gene_id %in% sig_genes),] 
cols_keep <-  c("id","seqid","start","end","gene_id","gene_name","gene_source","gene_biotype")


all_sig_anno <- data.table(inner_join(sig_genes_mat, sig_genes_temp[, cols_keep], by = "gene_id"))


## using metafor

cont = c("ME-STD","PD-STD","VDD-STD","PD-ME","VDD-ME","VDD-PD")
res.fe <- list()
pvals <- list()

for(contrast in cont){
  efsz <- diets$mat.efsz[[contrast]]
  efse <- diets$mat.efse[[contrast]]
  
  res.fe[[contrast]] <- sapply(1:nrow(efsz), function(x)
    if(length(efsz[x,][!is.na(efsz[x,])]) > 3 ){
      metafor::rma(yi=efsz[x,][!is.na(efsz[x,])], 
                   sei=efse[x,][!is.na(efse[x,])], method="DL")
    } else {
      NA
    }, simplify=F)
  
  pvals[[contrast]] <- c()
  for(i in 1:length(res.fe[[contrast]])){
    if(any(!is.na(res.fe[[contrast]][[i]]))){
      pvals[[contrast]][i] <- getElement(res.fe[[contrast]][[i]], "pval")
    }
  }
}

ind <- 1
forest.rma(res.fe[[contrast]][[order(pvals[[contrast]])[ind]]], slab=rixes[1:res.fe[[contrast]][[order(pvals[[contrast]])[ind]]]$k], atransf=exp)


subset[[1]]$gene[order(pvals)][1:25]
mat.rawpval[order(pvals)[1:25],]
efsz.vec[order(pvals)[1:25],]

########### all rixes together ########

rnaout <- list()
for(i in 1:58){
  rnaout[[i]] <- readRDS(paste0("rnaseq_allRIX_",i,".rds"))
}

add <- list()
counter <- 1
for(i in 1:58){
  for (j in 1:length(rnaout[[i]]$allSummary)){
    keep <- intersect(grep("6",rnaout[[i]]$allSummary[[j]]$Level), 
                      which(rnaout[[i]]$allSummary[[j]]$pval < 0.05))
    if(length(keep > 0)){
      add[[counter]] <- rnaout[[i]]$allSummary[[j]][keep, ]
      counter <- counter+1
    }
  }
}
signif6 <- do.call("rbind",add)
length(which(signif6$pval<(0.05/29000*9)))
keepPhen <- signif6[order(abs(signif6$mu), decreasing = T)
                    [1:length(which(signif6$pval<(0.05/29000*9)))],]$phenotype
keepPhen <- signif6[which(signif5$pval == 0),]
#ENSMUSG00000004542

JLs <- list()
for(i in 1:length(rnaout)){
  JLs[[i]] <- do.call("rbind",rnaout[[i]]$JL_compare)
}
JLs_fix <- lapply(JLs, function(x) ifelse(dim(x)[2] == 1, NULL, x))
JLs_fix <- sapply(1:58, function(x) ifelse(dim(JLs[[x]])[2] == 1, NULL, JLs[[x]]))
for(i in 1:58){ 
  if(dim(JLs[[i]])[2] == 1) JLs_fix[[i]] <- NULL
  else JLs_fix[[i]] <- JLs[[i]] 
}

allJL <- do.call("rbind",JLs_fix)
#keepJL[intersect(which(abs(keepJL$Ratio) != Inf), grep("^6$",keepJL$Level)),]
#nonInf[order(abs(nonInf$Diff), decreasing=T)[1:20],]
allJL[which(allJL$phenotype %in% unique(keepPhen)[1:5]),]



########## MIN CODE #############

#tryLam=c(-1, 0, .25, .33, .5, 1, 2, 3)
tryLam=c(1)

rixes <- cov_short$encoded$Level[cov_short$encoded$Variable == "RIX"]
phenInd <- grep("ENS",colnames(cov_short$df))
rankSet <- matrix(NA, nrow=nrow(cov_short$df), ncol=length(phenInd))

# start with poe when adding poe

flag <- c()
mouse <- list()
for(j in 1:length(rixes)){
  phenInd <- grep("ENS",colnames(cov_short$df))
  datau <- cov_short$df[which(paste(cov_short$df$RIX) == rixes[j]), ]
  colsums <- colSums(datau[,phenInd])
  #flag <- c(flag, which(colsums == 0))
  datause <- datau#[,-flag]
  phenInd <- grep("ENS",colnames(datause))
  
  mouse_i <- matrix(NA, nrow=length(phenInd), ncol=nrow(datause))
  for(m in 1:nrow(datause)){v
    mouse_i[,m] <- unlist(as.vector(datause[m,phenInd]))
  }
  mouse[[j]] <- mouse_i
  #mouse_j <- apply(datause, 1, function(x) unlist(as.vector(datause[x,phenInd])))
  #hist(mouse_i[,1])
  
  #rankSet[which(paste(cov_short$df$RIX) == rixes[j]), ] <- apply(datause[,phenInd], 2, 
  #                                                               function(x) rank(x, ties.method = "min")/length(!is.na(x)))
}

maxes <- lapply(mouse, function(x) apply(x, 2, max))
unlist(lapply(maxes, max))
plot(mouse[[1]][,14], mouse[[5]][,5])
lines(x=seq(1,2e10, by=1000), y=seq(1,2e10, by=1000), col="gray")

# alignment to transcriptome or rob's (kmer based method) selective alignment of reads that are trustworthy (which doesn't exist)
# plot expression from one mouse in highly expressed rix with one in poorly expressed rix, find differential expression ones, 
## investigate if it's due to mapping/snp density
# rix specific snps in transcriptome and realign

rankSet <- cbind(cov_short$df[,-phenInd], rankSet)
colnames(rankSet) <- colnames(cov_short$df)
#temp <- matrix(rankSet[,phenInd[23]])
#colnames(temp) <- colnames(rankSet[phenInd[23]])
indvariable="~ -1 + Diet + (0 + PO | RIX) + (1| Diet:RIX)"
indvariable="~ -1 + PO + Diet + RIX + (0+PO|Diet:RIX)"
indvariable="~ -1 + Diet + (PO + 0 | RIX) + (PO + 0|Diet:RIX)"


############ RANK WHOLE DATASET

datause <- cov_short$df[, -which(colSums(cov_short$df[,phenInd]) == 0)]
phenInd <- grep("ENS",colnames(datause))

rankData <- cbind( datause[,-phenInd],
                           t(apply(datause[,phenInd], 1, function(x) rank(x, ties.method = "min")/length(!is.na(x)))) )

fitLMER2 <- BC.model(y.mat=rankData[,phenInd[1:100]], data=rankData[,1:phenInd[100]], 
                     indvariable=indvariable, 
                     transformParams=getMatnutTransformParams(tryLam = tryLam, normd = F))

lmerSum <- lmer.getDecodedSummary(lmerOb=fitLMER2, formula=indvariable, phenotypes=fitLMER2$phenotypes)

fitLMER2$JAGSformula <- lmer.to.jags.form(indvariable)

jagsFit2 <- runJagsModels_cov(datalist=rankData[,1:tail(which(colnames(datause) %in% fitLMER2$phenotypes), n=1)], 
                              testLMER = fitLMER2, encoded=cov_short$encoded, phenotype=fitLMER2$phenotypes, n.iter=200000)

jagsSum <- list()
all.y <- matrix(NA, nrow=nrow(rankSet), ncol=length(fitLMER2$y.transform))
for(i in 1:length(fitLMER2$y.transform)){
  all.y[,i] <- fitLMER2$y.transform[[i]]
  all.mcmc <- mcmc.stack(jagsFit2[[i]]$fit)
  jagsSum[[i]] <- jags.getDecodedSummary(all.mcmc, cov_short$encoded)
}

allSummary <- list(lmer=lmerSum, jags=jagsSum)

allSummary$compare <- sapply(1:length(lmerSum), function(x) data.frame(compareSummary(allSummary$lmer[[x]],allSummary$jags[[x]])), 
                             simplify = F)

PO6 <- lapply(allSummary$compare, function(x) x[which(x$Level == c("PO6", "POVDD6")),"LMER_est"])

saveRDS(fitLMER, "xRIX_rna_27Nov2017.RDS")


system.time({
  rnaseq <- list()
  for(j in 1:length(rixes)){
    datause <- cov_short$df[which(paste(cov_short$df$RIX) == rixes[j]), ]
    rankSet[which(paste(cov_short$df$RIX) == rixes[j]), ] <- apply(datause[,phenInd], 2, 
                                                                   function(x) rank(x, ties.method = "min")/length(!is.na(x)))
    
    fitLMER <- BC.model(y.mat=datause[,phenInd[1:50]], data=datause[,1:phenInd[50]], 
                        indvariable="~ -1 + PO + Diet", 
                        transformParams=getMatnutTransformParams(tryLam = c(0, 0.5), normd = T))
    lmerSum <- lmer.getDecodedSummary(lmerOb=fitLMER, formula="~ -1 + PO + Diet", phenotypes=fitLMER$phenotypes)
    fitLMER$formula <- "~ PO + Diet"
    
    
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

allFits <- list()

for(i in 2:9){
  temp <- readRDS(paste0(dataSource,"/rna/rnaseq_out2_", rixes[i], ".RDS"))
  allFits[[i]] <- list()
  allFits[[i]]$summary <- lapply(temp$lmerObj$fits, function(x) summary(x)$coefficients)
  allFits[[i]]$anova <- temp$lmerObj$anovaOut
}

allFits[[1]]$summary[[1]]

getPvals <- function(rix){
  pval <- unlist(lapply(allFits[[rix]]$summary, function(x) x['DietVDD','Pr(>|t|)']))
  prop <- length(which(pval < 0.05)) / length(pval)
  hist(pval)
  return(list(pvals=pval, propo=prop))
}

pvals5 <- getPvals(5)
sorted5 <- sort(pvals5$pvals, index.return=T)


pvals_other <- getPvals(8)



mat <- list()
for(i in 1:9){
  mat[[i]] <- matrix(NA, nrow=length(get(paste0("rnaseq.",rixes[i]))$lmerObj$fits), ncol=4)
  
  for(j in 1:nrow(mat[[i]])){
    mat[[i]][j,] <- summary(get(paste0("rnaseq.",rixes[i]))$lmerObj$fits[[j]])$coefficients["DietVDD",]
  }
  
  plot(mat[[i]][,1], -log10(mat[[i]][,4]), xlab="Effect size", ylab="-log10(p-value)",
       main=paste("Effect of VDD in RIX", rixes[i]))
}

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


