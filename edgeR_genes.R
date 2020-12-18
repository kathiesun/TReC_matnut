library(tidyverse)
library(rstan)
library(edgeR)
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

d0 <- DGEList(counts)
d0 <- calcNormFactors(d0)
d0


cutoff <- 2
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0
if(length(drop) > 0) d <- d0[-drop,] 
dim(d)

snames <- colnames(counts) # Sample names

mm = model.matrix(~ 0 + Diet + RIX + PO*RIX + DietRIX + PO*DietRIX, matnut[match(snames,matnut$ID),])

y <- voom(d, mm, plot = T)
############################################
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

annot_genes = annot_genes %>% filter(Gene.Name %in% genes$V1)

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


dds = DESeqDataSetFromMatrix(countData = counts, colData = matnut[match(snames,matnut$ID),], 
                             design = ~ RRIX + Diet + PO_cat)    # + Diet:PO)
countData = countData[rowSums(counts(dds)) >= 10,]
dds <- estimateSizeFactors(dds)
dat = counts(dds, normalized=T)
dat = dat[rowMeans(dat) > 1,]
mod  <- model.matrix(~   RRIX + Diet + PO_cat, colData(dds))
mod0 <- model.matrix(~   1            , colData(dds))
svseq <- svaseq(dat, mod, mod0, n.sv = 1)  
ddssva = dds
ddssva$SV1 <- svseq$sv[,1]
design(ddssva) <- as.formula(~ SV1 + RRIX + Diet + PO_cat)
dds = ddssva
dds = DESeq(dds)
dds = nbinomWaldTest(dds)

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




