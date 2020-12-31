library(tidyverse)
library(rstan)
options(mc.cores = parallel::detectCores())

#setwd("C:/Users/Kathie/TReC_matnut/src")

setwd("~/TReC_matnut/src")

source("lmer_functions_rna.R")
source("stan_pheno_functions.R")
#source("prediction_functions.R")
source("summary_functions.R")

args <- commandArgs(trailingOnly = TRUE) 
it = as.numeric(args[1])

###### Read in data
#dir <- "C:/Users/Kathie/Dropbox\ (ValdarLab)"
dir <- "/nas/depts/006/valdar-lab/users/sunk/"
matnut = read.csv(file.path(dir,'matnut_main/AllMice_GeneExpression_SSupdated_11.27.19.csv'))
#matnut <- readRDS(file.path(dir,'phenotype_analysis/matnut_data.rds'))
#genes = read.csv("../deseq2/priorityTryGenes_16dec2020.csv", header=F)
#genes = genes$V1[-which(genes$V1 == "Gm23935")]


matnut = matnut %>% select(-contains("X"))
gene_count <- read.csv(file.path(dir,'/hisat2_stringtie_fasta/stringtie/round2/gene_count_matrix.csv'))
#gene_count <- read.csv(file.path(dir,'/trec/gene_count_matrix.csv'))

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

annot_genes = annot_genes %>% filter(Gene.Name %in% genes)

t_counts = data.frame(t(gene_count %>% filter(Gene.Name %in% annot_genes$Gene.Name) %>%
			select(contains("Pup.ID")))) 
tmp = gene_count$Gene.Name[which(gene_count$Gene.Name %in% annot_genes$Gene.Name)] 
#tmp[46] = "Snhg14_v2"
colnames(t_counts) = tmp 
t_counts$ID = unlist(as.character(rownames(t_counts)))
matnut_use = right_join(matnut, t_counts, "ID")
colnames(matnut_use) = gsub("[-]",".",colnames(matnut_use))
genes = gsub("[-]",".",genes)
#regions_list = readRDS(file.path(dir,"mini/seg_regions_by_genotyping_perRIX_kmerBased_23nov2019.rds"))
#haplofiles = list.files(file.path(dir,"mini/pup_haplo_blocks_by_CC_parent_dec2019"), 
#                        pattern="haploBlocks.rds", full.names = T)
#phased_par_CC_haplotypes = lapply(haplofiles, readRDS)
#tmp = do.call("rbind", strsplit(haplofiles, "_"))
#names(phased_par_CC_haplotypes) = tmp[,ncol(tmp)-1]

maxn = floor(nrow(annot_genes)/10)
its_1 = ((it-1)*10)+1
its_2 = ifelse(it==maxn, nrow(annot_genes), it*10)

#its_1 = 1
#its_2 = length(genes)
print(paste(its_1,its_2))

if(length(which(!genes %in% colnames(matnut_use))) > 0){
  genes = genes[-which(!genes %in% colnames(matnut_use))]
}
stanlist <- lapply(as.character(genes)[its_1:its_2], function(x) 
  stanSum(df=matnut_use, encoded=encoded, phenotype=x, 
                   randvar=c("RIX", "DietRIX"), fixvar="Diet", POvar=c("RIX", "DietRIX"),
                   tryLam=c(-1, 0, .25, .33, .5, 1, 2, 3), normd=T, 
                   chains=1, iter=2000)
	)

names(stanlist) = as.character(genes)[its_1:its_2]


saveRDS(stanlist, paste0(dir,"/trec/tot_expr_priorityGenes_16dec2020.rds"))
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
  
  
  
  
  
