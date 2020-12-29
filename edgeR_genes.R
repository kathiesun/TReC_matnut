library(dplyr)
library(ggplot2)
library(rstan)
library(edgeR)
library(DESeq2)
library(sva)

setwd("C:/Users/Kathie/TReC_matnut/src")

source("lmer_functions_rna.R")
source("stan_pheno_functions.R")
source("prediction_functions.R")
source("summary_functions.R")


###### Read in data
dir <- "C:/Users/Kathie/Dropbox\ (ValdarLab)"

#matnut <- readRDS(file.path(dir,'phenotype_analysis/matnut_data.rds'))
matnut = read.csv(file.path(dir,'matnut_main/AllMice_GeneExpression_SSupdated_11.27.19.csv'))
gene_count <- read.csv(file.path(dir,'/trec/gene_count_matrix.csv'))
genes = read.csv("../deseq2/priorityTryGenes_16dec2020.csv", header=F)

colnames(gene_count)[1] = "gene_id"
gene_count$Gene.ID = do.call("rbind",(strsplit(as.character(gene_count$gene_id), "[|]")))[,1]
gene_count$Gene.Name = do.call("rbind", (strsplit(as.character(gene_count$gene_id), "[|]")))[,2]

samples = colnames(gene_count)[grep("Pup", colnames(gene_count))]
matnut$ID = paste0("Pup.ID_",matnut$Pup.ID)
matnut$Diet = gsub(" $", "", matnut$Diet)
matnut$Diet = factor(matnut$Diet, levels=c("Standard","Low Protein","Methyl Enriched","Vitamin D Deficient"))
matnut$RIX = gsub("a|b","",matnut$Reciprocal)
matnut$RIX = factor(matnut$RIX, levels=c(1:4,6:10))
matnut$PO = ifelse(gsub("[0-9]","",matnut$Reciprocal) == "a", 0.5, -0.5)
matnut$DietRIX = gsub(" ", ".",paste0(matnut$Diet, matnut$RIX))
levs = apply(expand.grid(gsub(" ",".",levels(matnut$Diet)), levels(matnut$RIX)),1,function(x) paste(x, collapse=""))
matnut$DietRIX = factor(matnut$DietRIX, levels = levs[which(levs%in%matnut$DietRIX)])
matnut$DietRIXPOq = paste0(matnut$DietRIX,"_",matnut$PO)
matnut$RIX = factor(matnut$RIX, levels = c(1:4,6:10))
map = read.csv(file.path(dir,"mini/combined_map_cs_gbrs_allmarkers_17mar2020.csv"))
encoded <- getEncoding(matnut, terms = c("RIX","Diet","DietRIX"))

gene_count = data.frame(gene_count[-which(duplicated(gene_count$Gene.Name)),])
rownames(gene_count) = gene_count$Gene.Name
counts = gene_count[, grep("Pup.ID", colnames(gene_count))]
#counts = gene_count[which(gene_count$Gene.Name %in% genes$V1), grep("Pup.ID", colnames(gene_count))]

t_counts = data.frame(t(gene_count %>% filter(Gene.Name %in% genes$V1) %>%
                          select(contains("Pup.ID")))) 
tmp = gene_count$Gene.Name[which(gene_count$Gene.Name %in% genes$V1)] 
tmp[46] = "Snhg14_v2"
colnames(t_counts) = tmp 
t_counts$ID = unlist(as.character(rownames(t_counts)))
matnut_use = right_join(matnut, t_counts, "ID")
colData = matnut[match(colnames(counts),matnut$ID),c("Breeding.Batch","Behavior.Batch","RIX","Reciprocal","Diet",
                 "Dam.ID","ID","PO","DietRIX","DietRIXPOq")]

#######################################################

d0 <- DGEList(counts)
d0$samples = cbind(d0$samples, colData)
samplenames = d0$samples$ID

cpm0 <- cpm(d0)
lcpm0 <- cpm(d0, log=TRUE)

L <- mean(d0$samples$lib.size) * 1e-6
M <- median(d0$samples$lib.size) * 1e-6
c(L, M)
summary(lcpm0)

table(rowSums(d0$counts==0)==(ncol(d0$counts)))

keep.exprs <- filterByExpr(d0, group="Reciprocal")
d1 <- d0[keep.exprs,, keep.lib.sizes=FALSE]
dim(d0)
dim(d1)
d2 = d1[-grep("MSTRG", rownames(d1)),, keep.lib.sizes=FALSE]
dim(d2)

cutoff <- 1
drop <- which(apply(cpm(d2), 1, max) < cutoff)
d3 <- d2
if(length(drop) > 0) d3 <- d2[-drop,] 
dim(d3)


lcpm.cutoff <- log2(10/M + 2/L)
library(RColorBrewer)
nsamples <- ncol(d0)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(lcpm0[,1]), col=col[1], lwd=2, ylim=c(0,1.15), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm0[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
#legend("topright", samplenames, text.col=col, bty="n")
lcpm <- cpm(d3, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.2), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)

for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
#legend("topright", samplenames, text.col=col, bty="n")

snames <- colnames(d3$counts) # Sample names
d3 <- calcNormFactors(d3)


############## looking at PO only #################


d3$samples$RRIX = d3$samples$RIX

mm = model.matrix(~ 0 + Reciprocal, d3$samples) 
mm[,grep("RRIX", colnames(mm))] = apply(mm[,grep("RRIX", colnames(mm))], 2, function(x) x * d3$samples$PO)
cont.wt = makeContrasts(
  "Reciprocal10a-Reciprocal10b",
  "Reciprocal1a-Reciprocal1b",
  "Reciprocal2a-Reciprocal2b",
  "Reciprocal3a-Reciprocal3b",
  "Reciprocal4a-Reciprocal4b",
  "Reciprocal6a-Reciprocal6b",
  "Reciprocal7a-Reciprocal7b",
  "Reciprocal8a-Reciprocal8b",
  "Reciprocal9a-Reciprocal9b",
  levels = mm
)
#cont.wt = cont.wt * 0.5
y <- voom(d3, mm, plot = T)
vfit = lmFit(y, mm)
vfit <- contrasts.fit(vfit, contrasts=cont.wt)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")


summary(decideTests(efit))


############# RIX, diet, PO, and interactions model ############# 
RIX = colData$RIX
RRIX = d3$samples$RRIX
PO = d3$samples$PO 
Diet = d3$samples$Diet
contrasts(RIX) = contr.sum(9)
contrasts(Diet) = contr.sum(4)
mm = model.matrix(~ 0 + RIX + Diet + PO:RIX)    #+ RIX:Diet
y <- voom(d3, mm, plot = T)
vfit = lmFit(y, mm)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")

summary(decideTests(efit,p.value=0.1))
voom_sig = which(unlist(apply(decideTests(efit,p.value=0.1), 1, function(x) any(x[grep("PO",colnames(efit))] != 0))))
##################  DESeq2  ##########################
colData = colData %>% arrange(RIX, Diet, PO)
RIX = colData$RIX
RRIX = colData$RIX
Diet = colData$Diet
PO = colData$PO
DietRIX = colData$DietRIX
contrasts(RIX) = contr.sum(9)
contrasts(RRIX) = contr.sum(9)
contrasts(Diet) = contr.sum(4)
contrasts(DietRIX) = contr.sum(length(unique(colData$DietRIX)))
mm = model.matrix(~ 0 + RIX + Diet + PO:RIX)    

all.zero <- apply(mm, 2, function(x) all(x==0))
idx <- which(all.zero)
if(length(idx) > 0) mm <- mm[,-idx]

rownames(colData) = colData$ID
cts = counts[,match(colData$ID, colnames(counts))]
all(rownames(colData) == colnames(cts))
dds <- DESeqDataSetFromMatrix(countData = cts, 
                              colData = colData,
                              design = mm)

## pre-filtering
keep <- rowSums(counts(dds) >= 10) >= 10
dds <- dds[keep,]
remove = grep("MSTRG", rownames(counts(dds)))
if(length(remove) > 0) dds = dds[-remove,]

dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)

base = "RIX + Diet + PO:RIX"
dat = counts(dds, normalized=T)
dat = dat[rowMeans(dat) > 1,]
mod  <- model.matrix(as.formula(paste("~ 0 +", base)), colData(dds))
mod0 <- model.matrix(~   1                           , colData(dds))
svseq <- svaseq(dat, mod, mod0, n.sv = 1)  
ddssva = dds
ddssva$SV1 <- svseq$sv[,1]
  
design(ddssva) <- as.formula(paste("~ 0 + SV1 +", base))
dds = ddssva

dds <- DESeq(dds)
dds = readRDS(file.path(dir,"de_results/deseq_28dec2020.rds"))
res <- results(dds)

#saveRDS(dds, file.path(dir,"de_results/deseq_28dec2020.rds"))


PO_res = resultsNames(dds)[grep("PO", resultsNames(dds))]
POres_list = lapply(PO_res, function(x) results(dds, name = x))
POres_Ordered <- lapply(POres_list, function(x) x[order(x$pvalue),])
sig_genes_list = lapply(POres_Ordered, function(x) rownames(x)[which(x$padj < 0.1)])
sig_genes_short = lapply(sig_genes_list, function(x) 
  if(length(grep("Gm|Rik|[.]", x)) > 0) x[-grep("Gm|Rik|[.]", x)])
deseq_sig = sort(table(unlist(lapply(POres_Ordered, function(x) rownames(x[which(x$padj < 0.1),])))), decreasing = T)

genes = unique(unlist(lapply(sig_genes_short, function(x) x[1:min(5,length(x))])))

ddsPlot = dds
ddsPlot$PO = factor(ddsPlot$PO)
plotCts = plotCounts(ddsPlot, gene=sig_genes_short[[1]][3], intgroup=c("PO","RIX","Diet"), returnData = T)
ggplot(plotCts, aes(x=PO, y=count, col=Diet)) + 
  geom_point(position = position_jitter(width = 0.15)) + 
  theme_classic() + 
  facet_wrap(~RIX)


#######################################################
annot = read.table(file.path(dir, "mm10/ref_sequence/Mus_musculus.GRCm38.96.gtf"), skip=5, fill=T)
#annot = read.table(file.path(dir, "matnut_main/Mus_musculus.GRCm38.96.gtf"), skip=5, fill=T)

annot_genes = annot %>% filter(V15 == "gene_name") %>%
  dplyr::select(one_of(paste0("V",c(1,4,5,7,10,16)))) %>%
  distinct()

colnames(annot_genes) = c("Chr","Start","End","Strand","Gene.ID","Gene.Name")
annot_genes = annot_genes %>% mutate(Start = as.numeric(paste(Start)), End = as.numeric(paste(End)))
# 
if(length(which(duplicated(annot_genes$Gene.Name))) > 0){
  annot_genes = annot_genes[-which(duplicated(annot_genes$Gene.Name)),]
}

annot_genes = annot_genes[which(annot_genes$Gene.ID %in% gene_count$Gene.ID |annot_genes$Gene.Name %in% gene_count$Gene.Name),]

#annot_genes = annot_genes %>% filter(Gene.Name %in% genes$V1)

##############################


#regions_list = readRDS(file.path(dir,"mini/seg_regions_by_genotyping_perRIX_kmerBased_23nov2019.rds"))
#haplofiles = list.files(file.path(dir,"mini/pup_haplo_blocks_by_CC_parent_dec2019"), 
#                        pattern="haploBlocks.rds", full.names = T)
#phased_par_CC_haplotypes = lapply(haplofiles, readRDS)
#tmp = do.call("rbind", strsplit(haplofiles, "_"))
#names(phased_par_CC_haplotypes) = tmp[,ncol(tmp)-1]

maxn = floor(nrow(annot_genes)/10)
its_1 = ((it-1)*10)+1
its_2 = ifelse(it==maxn, nrow(annot_genes), it*10)

print(paste(its_1,its_2))

stanlist <- lapply(as.character(genes$V1)[its_1:its_2], function(x) 
  
  
  tmp = stanSum(df=matnut_use, encoded=encoded, phenotype=x, 
                randvar=c("RIX", "DietRIX"), fixvar="Diet", POvar=c("RIX", "DietRIX"),
                tryLam=c(-1, 0, .25, .33, .5, 1, 2, 3), normd=T, 
                chains=3, iter=6000)
)
names(stanlist) = as.character(genes$V1)[its_1:its_2]
saveRDS(stanlist, paste0(dir,"/trec/tot_expr_priorityGenes_",it,"_16dec2020.rds"))
#  stanmcmc<- As.mcmc.list(resStan)
#  summcmc <- summary(stanmcmc)

#armform <- y ~ Diet + PO*RIX + Diet*RIX + PO*Diet*RIX
#test_stanarm = stan_glm(armform, data=df, prior = lasso(autoscale = T))

#  par(mfrow=c(3,3))
#  coda::traceplot(stanmcmc[[1]][,grep("SPO[[]", colnames(stanmcmc[[1]])),drop=F])
#  coda::traceplot(stanmcmc[[1]][,grep("sigma|lambda", colnames(stanmcmc[[1]])),drop=F])
#  tset <- sapply(1:(ncol(stanmcmc[[1]])-1), function(x) HPDinterval(stanmcmc[,x,drop=T]), simplify=F)







#reg.jags <- jags.model(textConnection(madeModel$modelFin), data=madeModel$data, n.chains = chains, n.adapt = n.adapt)
#update(reg.jags, n.iter=n.adapt)
#fitJags[[phenotype[i]]]$fit <- coda.samples(reg.jags, variable.names = madeModel$paramUse, thin=thin, n.iter=n.iter)
#fitJags[[phenotype[i]]]$madeModel <- madeModel
#allnames <- c()
#paramUse <- gsub("beta","", unique(unlist(
#  lapply(strsplit(varnames(fitJags[[phenotype[i]]]$fit),'[[]'), function(x) x[[1]]))))
#wantParam <- intersect(paramUse, c(paste(unique(encoded$Variable)), "grand_mu"))

colData = matnut[match(snames,matnut$ID),c("Breeding.Batch","Behavior.Batch","RIX","Reciprocal","Diet",
                                           "Dam.ID","ID","PO","DietRIX","DietRIXPOq")]
rownames(colData) = colData$ID
dds = DESeqDataSetFromMatrix(countData = counts, colData = colData, 
                             design = ~ RIX + PO*RIX + Diet)    # + Diet:PO)
counts = counts[rowSums(counts(dds)) >= 10,]
dds <- estimateSizeFactors(dds)
dat = counts(dds, normalized=T)
dat = dat[rowMeans(dat) > 1,]
mod  <- model.matrix(~   RIX + PO*RIX + Diet, colData(dds))
mod0 <- model.matrix(~   1                  , colData(dds))
svseq <- svaseq(dat, mod, mod0, n.sv = 1)  
ddssva = dds
ddssva$SV1 <- svseq$sv[,1]
design(ddssva) <- as.formula(~ SV1 + RIX + PO*RIX + Diet)
dds = ddssva
dds = DESeq(dds)
POresnames = resultsNames(dds)[grep("PO", resultsNames(dds))]
POres_05 = lapply(POresnames, function(x) results(dds, name=x, alpha=0.05))

res[which(rownames(res) == "Gm23935"),]
#dds = nbinomWaldTest(dds)

########### EDA ###########

hist(log10(counts + 1))
abline(v=1, col="blue", lwd=3)
length(which(rs < 10))

par(mfrow=c(1,2))
plot(countData[,1:2])
logcts <- log10(countData + 1)
plot(logcts[,1:2])

hist(logcts[,1])

rv <- rowVars(logcts)
o <- order(rv, decreasing=TRUE)[1:500]
pc <- prcomp(t(logcts[o,]))
plot(pc$x[,1:2])
dat <- data.frame(pc$x[,1:2], RRIX=factor(colData$RRIX), Diet=factor(colData$Die))
ggplot(dat, aes(PC1, PC2, col=RRIX)) + geom_point() + theme_classic()
ggplot(dat, aes(PC1, PC2, col=Diet)) + geom_point() + theme_classic()


plot(pc$sdev[1:10]^2 / sum(pc$sdev^2), type="b", ylab="% Var")


######## RUN DESEQ2 #########

## new way
#dds_string_sva = run_DESeq2(colData, countData, sva=T, interact=F)
rem_seg_dds  = run_DESeq2(colData, countData, sva=T, interact=F, 
                          contrast="PO", annot=annot_genes, regions_list=regions_list)
keepSeg_dds_diet  = run_DESeq2(colData, countData, sva=T, interact=F, 
                               contrast="Diet", annot=annot_genes, regions_list=NULL)
keep_seg_dds = run_DESeq2(colData, countData, sva=T, interact=F, 
                          contrast="Diet", annot=annot_genes, regions_list=NULL)
rem_seg_res  = compare_diets_PO(dds_lst=rem_seg_dds, fdrtool_adjust=T, 
                                adjust_p=F, ref_diet = "Vitamin.D.Deficient")
remSeg_res_diet = compare_diets_PO(dds_lst=remSeg_dds_diet, fdrtool_adjust=T, 
                                   adjust_p=F, ref_diet = "Vitamin.D.Deficient")
keepSeg_res_diet = compare_diets_PO(dds_lst=keepSeg_dds_diet, fdrtool_adjust=T, 
                                    adjust_p=F, ref_diet = "Vitamin.D.Deficient")


remSeg_dds_interDiet  = run_DESeq2(colData, countData, sva=T, interact=T, 
                                   contrast="Diet", annot=annot_genes, regions_list=regions_list)
remSeg_res_diet = compare_diets_PO(dds_lst=remSeg_dds_interDiet, fdrtool_adjust=T, 
                                   adjust_p=F, ref_diet = "Vitamin.D.Deficient")
interResLst = lapply(remSeg_dds_interDiet, function(x){
  tmp = results(x)
  dtmp = data.frame(tmp)
  dtmp$gene = rownames(tmp)
  fdr <- fdrtool(dtmp$stat, statistic="normal", plot=F)
  dtmp$pval_raw = dtmp$pvalue
  dtmp$padj_deseq2 = dtmp$padj
  dtmp$pvalue = fdr$pval
  dtmp$padj  <- p.adjust(dtmp$pvalue, method = "BH")
  dtmp %>% arrange(pvalue)
})




