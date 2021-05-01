options(digits=4)

setwd("C:/Users/Kathie/TReC_matnut/src")
library(tidyverse)
library(fdrtool)
library(DESeq2)
library(rjags)
library(apeglm)

source("prediction_functions.R")
source("summary_functions.R")
source("ase_summary_source.R")
#source("stan_pheno_functions.R")


dir <- "C:/Users/Kathie/Dropbox\ (ValdarLab)"
#gecco <- read.csv(file.path(dir, "de_results/fullGeccoRnaDump.csv"))

C_diet_4 = read.csv(file.path(dir, "variant_data/C_matrix_4diets.csv"), header=F)
C_diet_2 = read.csv(file.path(dir, "variant_data/C_matrix_2diets.csv"), header=F)
C_diet_3 = matrix(c(0.9659258,-0.2588190,-0.7071068,-0.2588190,0.9659258,-0.7071068),nrow=3,ncol=2)

#gene_count <- read.csv(file.path(dir,'/trec/gene_count_matrix.csv'))
#colnames(gene_count)[1] = "gene_id"
#gene_count$Gene.ID = do.call("rbind",(strsplit(as.character(gene_count$gene_id), "[|]")))[,1]
#gene_count$Gene.Name = do.call("rbind", (strsplit(as.character(gene_count$gene_id), "[|]")))[,2]

gene_count <- read.csv(file.path(dir,'/trec/gene_count_matrix_hetsOnly.csv'))
rownames(gene_count) = gene_count$Gene.Name
samples = colnames(gene_count)[grep("Pup", colnames(gene_count))]

matnut = read.csv(file.path(dir,'matnut_main/AllMice_GeneExpression_SSupdated_11.27.19.csv'))
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

#gene_count = data.frame(gene_count[-which(duplicated(gene_count$Gene.Name)),])
#rownames(gene_count) = gene_count$Gene.Name
bwt_sizes = read.csv(file.path(dir, "matnut_main/meta_bwtsize_jun2019.csv"))
bwt_sizes$ID = paste0("Pup.ID_",bwt_sizes$Pup.ID)


counts = gene_count[, grep("Pup.ID", colnames(gene_count))]

colData = matnut[match(colnames(counts),matnut$ID),
                 c("Pup.ID","Breeding.Batch","Behavior.Batch","RIX","Reciprocal","Diet", 
                   "Dam.ID","ID","PO","DietRIX","DietRIXPOq")]
ie_genes = read.csv("C:/Users/Kathie/Dropbox (ValdarLab)/imprinted_genes/2014_05_14_allImprintedGenes.csv")

##########################

regFiles = list.files(pattern = "summary_mnt_50k_29dec2020.rds", file.path(dir, "variant_data/regression_outputs/"),
                      full.names = T)

regs = lapply(regFiles, readRDS)
mu_g = do.call("rbind", lapply(regs, function(x) x[[1]]$summary$mu_g))
mu_g = mu_g[,-intersect(grep("[.]",colnames(mu_g)), grep("^[a-z]", colnames(mu_g)))]
mu_g$Gene.Name = unlist(strsplit(rownames(mu_g),"_"))[c(F,T)]
mu_g$Pup.ID  = as.integer(unlist(strsplit(rownames(mu_g),"_"))[c(T,F)])
mu_g = mu_g %>% left_join(colData, by="Pup.ID")

mu_g$lower_opp = mu_g$lower
mu_g$upper_opp = mu_g$upper
mu_g$lower_opp[which(mu_g$PO < 0)] = 1-mu_g$upper[which(mu_g$PO < 0)]
mu_g$upper_opp[which(mu_g$PO < 0)] = 1-mu_g$lower[which(mu_g$PO < 0)]

imp_perc = mu_g %>% #filter(lower > 0.5 | upper < 0.5) %>%
  group_by(Gene.Name, RIX) %>% 
  summarize(over = sum(lower_opp > 0.5), under = sum(upper_opp < 0.5), total = n(),
            perc = max(over/total, under/total), 
            ratio = max(over/total, under/total) / min(over/total, under/total)) %>%
  filter(perc > 0.25, ratio > 3, total > 3) %>%
  arrange(desc(perc), desc(ratio),desc(total))
  #left_join(mu_g %>% group_by(Gene.Name) %>% tally(), by="Gene.Name") %>%
  #mutate(imp_perc = n.x / n.y) %>%
  #arrange(desc(imp_perc)) %>% 
  #filter(n.x > 3, imp_perc > 0.1)
imp_perc %>% filter(Gene.Name %in% ie_genes$mgi_symbol)

#write.csv(imp_perc, file.path(dir, "trec/priority_ase_genes_11jan2021.csv"))

## look specifically at known imprinted genes
tail((mu_g %>% 
  filter(Gene.Name %in% ie_genes$mgi_symbol) %>%
  group_by(Gene.Name, RIX) %>% 
  summarize(over = sum(lower_opp > 0.5), under = sum(upper_opp < 0.5), total = n(),
            perc = max(over/total, under/total), 
            ratio = max(over/total, under/total) / min(over/total, under/total)) %>%
  filter(perc > 0, !is.nan(ratio), total > 3) %>%
  arrange(desc(perc), desc(ratio),desc(total))), n=50)

mu_g_ie = mu_g %>% filter(Gene.Name %in% imp_perc$Gene.Name)


#imp_genes = imp_perc$Gene.Name[(which(imp_perc$Gene.Name %in% unlist(sig_genes_short)))]
keep_genes = imp_perc$Gene.Name[(which(imp_perc$Gene.Name %in% unlist(sig_genes_list)))]


print(mu_g %>% filter(lower > 0.55 | upper < 0.45) %>%
    group_by(Gene.Name) %>% tally() %>%
    left_join(mu_g %>% group_by(Gene.Name) %>% tally(), by="Gene.Name") %>%
    mutate(imp_perc = n.x / n.y) %>%
    arrange(desc(imp_perc)) %>% 
    filter(n.x > 3, imp_perc > 0.1), n=30)


gregg_genes = read.csv("../gregg_POgenes.csv")
colnames(gregg_genes)[1] = "Gene.Name"

#gregg_genes = gregg_genes %>% filter(sample != "E15", ie_status != "known") %>%
#  select(-ucsc_ID_SNP) %>% distinct()

gregg_ie  = unique(gregg_genes$Gene.Name)[which(unique(gregg_genes$Gene.Name) %in% ie_genes$mgi_symbol)]
gregg_ase = unique(gregg_genes$Gene.Name)[which(unique(gregg_genes$Gene.Name) %in% imp_perc$Gene.Name)]
gregg_des = unique(gregg_genes$Gene.Name)[which(unique(gregg_genes$Gene.Name) %in% unlist(sig_genes_list))]

#de_genes = read.table(file.path(dir, "trec/priority_deseq_genes_26jan2021.txt"))
de_genes = read.csv(file.path(dir, "trec/priority_deseq_genes_30jan2021.csv"))
ase_genes = read.csv(file.path(dir, "trec/priority_ase_genes_11jan2021.csv"))

all_genes = unique(c(gregg_genes$Gene.Name, de_genes$gene, ase_genes$Gene.Name, ie_genes$mgi_symbol))

################################################################


indiv_pups = list.files(file.path(dir, "mini/pup_haplo_blocks_by_CC_parent_jan2021"), pattern="haploBlocks", full.names = T)

phased_CC_haplotype = lapply(indiv_pups, readRDS)
tmp = do.call("rbind", lapply(indiv_pups, function(x) unlist(strsplit(x, "_"))))
names(phased_CC_haplotype) = paste0("Pup.ID_", tmp[,ncol(tmp)-1])

seg_regions = readRDS(file.path(dir,"mini/seg_regions_by_genotyping_perRIX_kmerBased_23nov2019.rds"))
masterSnps_files = list.files(file.path(dir, "variant_data/kmers_to_run"), pattern="masterSnps_chr", full.names = T)
count_files <- list.files(file.path(dir,"variant_data/kmer_counts"), ".csv", full.names = T)

## CC labels
#lab = read.csv("C:/Users/Kathie/Dropbox (ValdarLab)/variant_data/matched_v2_4jun2018.csv")
lab = read.csv(file.path(dir, "matnut_main/matched_v2_4jun2018.csv"))
lab = read.csv(file.path(dir, "matnut_main/CC_labels_paternity.csv"))
#lab %>% dplyr::select(one_of("Pup.ID","RRIX", "CC.1","CC.2")) %>%
#  filter(!is.na(Pup.ID)) -> lab


rix_6_haps = do.call("rbind", lapply(lab$Pup.ID[which(lab$RIX == 6)], function(x) phased_CC_haplotype[[paste0("Pup.ID_", x)]]$full_founders %>% filter(chr==7, start>41730195, end<67030195)))
write.csv(rix_6_haps, file.path(dir, "rix_6_haps.csv"))


## Pup demographic info
pupInfo = matnut[,c(1,2,4,5,8,10,12,13,16,17)]
pupInfo %>% mutate(dir = gsub("[0-9]","", Reciprocal)) %>%
  mutate(Diet = gsub(" ","", Diet)) -> pupInfo

data_kmers = ratios_lst = list()

#if(j == 1){
  ## snp info

for(c in 1:19){

  masterSnps = readRDS(masterSnps_files[grep(paste0("chr",c,"_"), masterSnps_files)])
  masterSnps = masterSnps[[1]]
  masterSnps <- do.call("cbind", masterSnps)
  masterSnps$seq.consensus <- paste0(masterSnps$seq.end5, masterSnps$seq.end3)
  #useSnps = masterSnps %>% filter(seq.Gene %in% all_genes)
  
  ## count data 
  chr_files <- count_files[grep(paste0("chr", c,"_"), count_files)]
  ref_file <- chr_files[grep("ref", chr_files)]
  alt_file <- chr_files[grep("alt", chr_files)]
  reforig <- unique(read.csv(ref_file, header = T))
  altorig <- unique(read.csv(alt_file, header = T))
  reforig$X = altorig$X <- NULL
  #reforig = reforig %>% filter(k.mer %in% useSnps$seq.refseq)
  #altorig = altorig %>% filter(k.mer %in% useSnps$seq.altseq)
  reforig$pup = reforig$pup.id 
  altorig$pup = altorig$pup.id 
  
  data_kmers[[c]] = process_and_plot(chr=c, 
                                snp_info=masterSnps, 
                                sample_info=pupInfo, 
                                RIX_info=lab, 
                                ref_counts=reforig, alt_counts=altorig, 
                                phased_CC_haplotype=phased_CC_haplotype, 
                                use_gene=F,
                                problemPups = c(1404, 1716, 1371, 569, 1911, 1951, 1015),
                                seg_regions=seg_regions)
 
  #data_kmers$podietrix = paste(data_kmers$dir, data_kmers$DietRIX, sep="_")
  print(c)

  print(paste(c, "done"))
  
}



pvals = do.call("rbind", lapply(ratios_lst, function(z) 
  do.call("rbind", sapply(1:length(z), function(y) {
    tmpp = data.frame(do.call("rbind", lapply(z[[y]]$freq, function(x) {   
      tmp = data.frame(x$p.value)
      tmp$rix = names(z)[y]
      tmp
    })))
    tmpp$gene = rownames(tmpp)
    tmpp
  }, simplify=F))
))

pval_list = lapply(unique(pvals$gene), function(x) pvals[which(pvals$gene == x),])
#library(metaRNASeq)
#pval_comb = fishercomb(pval_list)

files = list.files(file.path(dir, "trec/data_kmers_from_process_and_plot/"), pattern=".txt",
                   full.names = T)


#data_kmers_list = saveRDS(file.path(dir, "trec/data_kmers_from_process_and_plot/data_kmers_list_out_29jan2021.rds"))
data_kmers_list = readRDS(file.path(dir, "trec/data_kmers_from_process_and_plot/data_kmers_list_out_29jan2021.rds"))


#lab = unique(do.call("rbind", lapply(data_kmers_list, function(x) 
#  x[,c("Pup.ID","CCs","RIX","Reciprocal","DamLine_NewCC_ID","SireLine_NewCC_ID","CC_lab","RRIX","CC.1","CC.2" )]
#  )))



#data_kmers_list = lapply(files, read.table, sep="\n")
#data_kmers_list = lapply(data_kmers_list, function(x) 
#  do.call("rbind", data.frame(apply(x, 1, function(y) strsplit(y," ")))))

#for(i in 1:length(data_kmers_list)){
#  colnames(data_kmers_list[[i]]) = data_kmers_list[[i]][1,]
#  data_kmers_list[[i]] = data_kmers_list[[i]][-1,]
#  rownames(data_kmers_list[[i]]) = NULL
#  data_kmers_list[[i]] = data.frame(data_kmers_list[[i]])
#  data_kmers_list[[i]][,c("CC_1","CC_2","sum","kRat","kLB","kUB")] = 
#    apply(data_kmers_list[[i]][,c("CC_1","CC_2","sum","kRat","kLB","kUB")], 2, as.numeric)
#  data_kmers_list[[i]]$Pup.ID = as.numeric(data_kmers_list[[i]]$Pup.ID)
#  data_kmers_list[[i]] = data_kmers_list[[i]] %>% left_join(lab[, -which(colnames(lab) == "RRIX")], by="Pup.ID")
#  data_kmers_list[[i]]$DamLine_NewCC_ID = unlist(strsplit(data_kmers_list[[i]]$DamLine_NewCC_ID, "/"))[c(T,F)]
#  data_kmers_list[[i]]$SireLine_NewCC_ID = unlist(strsplit(data_kmers_list[[i]]$SireLine_NewCC_ID, "/"))[c(T,F)]
#  data_kmers_list[[i]]$l_count = data_kmers_list[[i]]$CC_1
#  data_kmers_list[[i]]$r_count = data_kmers_list[[i]]$CC_2
#  data_kmers_list[[i]]$mat_count = data_kmers_list[[i]]$CC_1
#  data_kmers_list[[i]]$pat_count = data_kmers_list[[i]]$CC_2
  
#  data_kmers_list[[i]]$mat_count[which(data_kmers_list[[i]]$DamLine_NewCC_ID == data_kmers_list[[i]]$CC.2)] = 
#    data_kmers_list[[i]]$CC_2[which(data_kmers_list[[i]]$DamLine_NewCC_ID == data_kmers_list[[i]]$CC.2)]
#  data_kmers_list[[i]]$pat_count[which(data_kmers_list[[i]]$SireLine_NewCC_ID == data_kmers_list[[i]]$CC.1)] = 
#    data_kmers_list[[i]]$CC_1[which(data_kmers_list[[i]]$SireLine_NewCC_ID == data_kmers_list[[i]]$CC.1)]
#  data_kmers_list[[i]]$CC_1 = data_kmers_list[[i]]$mat_count
#  data_kmers_list[[i]]$CC_2 = data_kmers_list[[i]]$pat_count
#  data_kmers_list[[i]]$mat_count = data_kmers_list[[i]]$pat_count = NULL
#  data_kmers_list[[i]]$CC_lab = paste(data_kmers_list[[i]]$DamLine_NewCC_ID,data_kmers_list[[i]]$SireLine_NewCC_ID, sep="/")
#  data_kmers_list[[i]] = data_kmers_list[[i]] %>% select(one_of("seq.Gene","seq.Chromosome","seq.rsId","seq.Position",
#                            "CC_1","CC_2","sum","CC1_hap","CC2_hap",
#                            "Pup.ID","CCs","logSum","Breeding.Batch","Behavior.Batch","RIX","Reciprocal","Diet",
#                            "DamLine_NewCC_ID","Dam.ID","SireLine_NewCC_ID", "Sire.ID",
#                            "CC_lab","RRIX","pup_gene","DietRIX","CC.1","CC.2","l_count","r_count"))
#}

#######################################################
runAPEASE = function(r, data_kmers_list){
  test_dat = do.call("rbind", lapply(data_kmers_list, function(x) x %>% filter(RRIX == r)))
  test_dat = test_dat %>% group_by(Pup.ID, seq.Gene) %>% 
    summarize(mat_tot = sum(CC_1), pat_tot = sum(CC_2))
  genes = unique(test_dat$seq.Gene)
  pups = unique(test_dat$Pup.ID)
  mat_tots = test_dat %>% select(c(Pup.ID, seq.Gene, mat_tot)) %>%
    spread(seq.Gene, mat_tot)
  mat_tots = t(mat_tots)
  colnames(mat_tots) = paste0("mat_",mat_tots[1,])
  pat_tots = test_dat %>% select(c(Pup.ID, seq.Gene, pat_tot)) %>%
    spread(seq.Gene, pat_tot)
  pat_tots = t(pat_tots)
  colnames(pat_tots) = paste0("pat_",pat_tots[1,])
  pat_tots = pat_tots[-1,]
  mat_tots = mat_tots[-1,]
  ase_sums = mat_tots + pat_tots
  
  ase_tots = cbind(mat_tots,pat_tots)
  ase_tots_comp = ase_tots[complete.cases(ase_tots),]
  ase_sums_comp = ase_sums[complete.cases(ase_tots),]
  mat_comp = mat_tots[complete.cases(ase_tots),]
  
  colData = pupInfo %>% filter(RIX == r, Pup.ID %in% pups) %>%
    mutate(ID = paste0("Pup.ID_",Pup.ID), 
           RIX = factor(RIX, levels=c(1:4,6:10)),
           Diet = factor(Diet, levels=c("Standard","LowProtein","MethylEnriched","VitaminDDeficient")),
           dir = factor(dir)) 
  rownames(colData) = colData$ID
  
  colnames(ase_sums_comp) = gsub("mat_","Pup.ID_", colnames(ase_sums_comp))
  colnames(mat_comp) = gsub("mat_","Pup.ID_", colnames(mat_comp))
  reorder_cols = match(colData$ID, colnames(ase_sums_comp))[which(!is.na(match(colData$ID, colnames(ase_sums_comp))))]
  ase_sums_comp = ase_sums_comp[,reorder_cols]
  mat_comp = mat_comp[,reorder_cols]
  
  
  all(rownames(colData) == colnames(ase_sums_comp))
  
  dds <- DESeqDataSetFromMatrix(countData = ase_sums_comp,
                                colData = colData,
                                design = ~ 1)
  dds <- DESeq(dds)
  res <- results(dds, name="Intercept")
  
  #resLFC <- lfcShrink(dds, coef="Intercept", type="apeglm")
  #resLFC
  
  n <- length(pups)
  f <- factor(rep(1:2,each=n))
  
  theta.hat <- 1000 # rough initial estimate of dispersion
  #x <- model.matrix(~ RIX, data=colData)
  x = matrix(rep(1, n))
  niter=5
  for (i in 1:niter) {
    param <- cbind(theta.hat, ase_sums_comp)
    #param <- cbind(theta.hat, ase_tots_comp[,1:n])
    
    fit.mle <- apeglm(Y=mat_comp, x=x, log.lik=NULL, param=param,
                      no.shrink=TRUE, log.link=FALSE, method="betabinCR")
    theta.hat <- bbEstDisp(success=mat_comp, size=ase_sums_comp,
                           x=x, beta=fit.mle$map,
                           minDisp=.01, maxDisp=5000)
  }
  
  coef <- 1
  #xlab <- "mean of total counts"
  #plot(rowMeans(ase_sums_comp), fit.mle$map[,1], log="x", xlab=xlab, ylab="log odds")
  
  
  mle <- cbind(fit.mle$map[,coef], fit.mle$sd[,coef])
  param <- cbind(theta.hat, ase_sums_comp)
  fit2 <- apeglm(Y=mat_comp, x=x, log.lik=NULL, param=param,
                 coef=coef, mle=mle, threshold=0.7,
                 log.link=FALSE, method="betabinCR")
  
  #ylim <- c(-1,1.5)
  s.val <- svalue(fit2$thresh) # small-or-false-sign value
  #cols <- ifelse(s.val < .01, "red", "black")
  #plot(rowMeans(ase_sums_comp), fit2$map[,coef], main="apeglm",
  #     log="x", xlab=xlab, ylab="log odds", col=cols, ylim=ylim)
  #abline(h=0,col=rgb(1,0,0,.5))
  fit2$map = cbind(fit2$map, abs(fit2$map[,coef]))
  fit2$map = fit2$map[order(fit2$map[,coef+1], decreasing = T),]
  
  return(fit2)
}


fit2_list = list()
for(x in c(1:4,6:10)){
  fit2_list[[x]] = runAPEASE(x, data_kmers_list)
}

lapply(fit2_list, function(x) x$map[which(x$map[,2] > 0.619),])
## log-odds 0.619 = ratio 0.35


sig_genes = unique(unlist(lapply(lapply(fit2_list, function(x) x$map[which(x$map[,2] > logit(0.65)),]), rownames)))

###################################

for (r in c(2:4,6:10)){
  test_dat = do.call("rbind", lapply(data_kmers_list, function(x) x %>% filter(RRIX == r)))


  test_dat = test_dat %>% group_by(Pup.ID, seq.Gene, Diet, Reciprocal, pup_gene) %>% 
    summarize(mat_tot = sum(CC_1), pat_tot = sum(CC_2)) %>%
    mutate(Pup.ID = factor(Pup.ID), 
           seq.Gene = factor(seq.Gene), 
           Diet = factor(Diet), 
           Reciprocal = factor(Reciprocal))
  test_dat$pup_gene = factor(test_dat$pup_gene)
  
  test_dat$tot = test_dat$mat_tot + test_dat$pat_tot
  encoded = getEncoding(test_dat, terms=c("Pup.ID","seq.Gene", "Diet", "Reciprocal", "pup_gene") )

  fit_list = list()
  for (g in encoded$Level[which(encoded$Variable == "seq.Gene")]){
    temp_dat = test_dat %>% filter(seq.Gene == g)
    try(fit_list[[g]] <- vglm(cbind(mat_tot, pat_tot) ~ 1, betabinomial (zero = 2),
                              data = temp_dat, trace = TRUE),
        silent=T)
    
    #coef(fit2)
  }
  saveRDS(fit_list, paste0("C:/Users/Kathie/Dropbox (ValdarLab)/ase_vgam/ase_intOnly_rix",r,"_1may2021.rds"))
  
  
  
}


#test_dat[,1:5] = apply(test_dat[,1:5], 2, as.factor)
#test_dat$Pup.ID = gsub(" ", "",test_dat$Pup.ID)
genes = unique(test_dat$seq.Gene)
pups = unique(test_dat$Pup.ID)


mat_tots = test_dat %>% dplyr::select(c(Pup.ID, seq.Gene, mat_tot)) %>%
  spread(seq.Gene, mat_tot)
mat_tots = t(mat_tots)
colnames(mat_tots) = paste0("mat_",mat_tots[1,])
pat_tots = test_dat %>% dplyr::select(c(Pup.ID, seq.Gene, pat_tot)) %>%
  spread(seq.Gene, pat_tot)
pat_tots = t(pat_tots)
colnames(pat_tots) = paste0("pat_",pat_tots[1,])
pat_tots = pat_tots[-1,]
mat_tots = mat_tots[-1,]
ase_sums = mat_tots + pat_tots

ase_tots = cbind(mat_tots,pat_tots)


temp_p = encoded[which(encoded$Variable == "Pup.ID"),]
test_dat$indP = temp_p$Index[match(test_dat$Pup.ID,temp_p$Level)]
temp_g = encoded[which(encoded$Variable == "seq.Gene"),]
test_dat$indG = temp_g$Index[match(test_dat$seq.Gene,temp_g$Level)]
temp_pg = encoded[which(encoded$Variable == "pup_gene"),]
test_dat$indPG = temp_pg$Index[match(test_dat$pup_gene,temp_pg$Level)]
temp_PO = encoded[which(encoded$Variable == "Reciprocal"),]
test_dat$indPO = temp_g$Index[match(test_dat$Reciprocal,temp_PO$Level)]

jagsdat <-  list(N       = nrow(test_dat),
                 y       = test_dat$mat_tot, 
                 n       = test_dat$mat_tot + test_dat$pat_tot,
                 nP      = length(unique(test_dat$Pup.ID)),
                 X           = matrix(rep(1, nrow(test_dat)), ncol=1),
                 indP        = test_dat$indP,
                 nG          = length(unique(test_dat$seq.Gene)),
                 indG        = test_dat$indG,
                 nGP         = length(unique(test_dat$pup_gene)), 
                 indSPO      = test_dat$indPO,
                 a_s = 10, b_s = 10)

jagsModel <-"
model {
## sampling
for (i in 1:N){
   y[i] ~ dbetabin(p_g[i]*theta, (1-p_g[i])*theta, n[i])
}
## priors
for(i in 1:N){
  for (g in 1:nG){
     p_g[indG[indP[i]]] ~ ilogit(X[i,] %*% betas)
  }
}

theta ~ dgamma(a_s, b_s)
betas ~ dt(mu, tau, 10)

}
"
write.table(jagsModel, "temp_ASE_genomewide.txt", quote = F, row.names = F, col.names = F)
line_inits <- list(list("p_g" = 0.5, "betas"= 1, "theta" = 1))
model <- jags.model("temp_ASE_genomewide.txt", data = jagsdat, inits=inits, n.chains=2)
update(jagsModel, n.iter=5000)
samples <- coda.samples(model, variable.names=c("p_g", "theta", "betas"),
                        n.iter=5000)


jags <- jags.model(model,
                   data = list("y_g"=y_g, "N_g"=N_g, "nG"=nG),
                   inits = jags.init, 
                   n.chains = nchains,
                   n.adapt = niter*0.1)

update(jags, niter)

res[[i]] <- jags.samples(jags,
                         jags.params,
                         niter)

library(VGAM)


########### -----------------------
ratios_lst = list()
data_kmers_list_df = do.call("rbind",data_kmers_list)
ratios_lst = run_stan_regress(data_kmers=data_kmers_list_df, 
                              niter=1000, n.thin=5,  
                              seg_regions=seg_regions,
                              save_dir=NULL, 
                              STZ=T, use_gene=F,
                              no_theta=F, alpha=NULL,
                              stan=F, stanMod = "ase_beta_binom_disperse.stan")   #ase_mu_g_regr.stan, ase_mu_g_simple.stan

data_kmers_list = data_kmers_list %>% arrange(RIX, Reciprocal, Pup.ID)
encoded = getEncoding(data_kmers_list, terms=c("RIX","Reciprocal", "Diet", "DietRIX", "seq.Gene", "pup_gene") )
stan_ase(data_kmers_list, encoded=data_kmers_list, 
                     mu_g0=0.5, t20=0.001, 
                     terms=c(Reciprocal), STZ=T, use_gene=F, 
                     nchains=2, iter=1000, 
                     plot=T, contrasts=F,
                     C_diet_4=NULL, C_diet_2=NULL,C_diet_3=NULL,
                     stanMod = NULL)


data_kmers_list_df = do.call("rbind", data_kmers_list)

lapply(names(ratios_lst[[1]]), function(y) 
  print(hist(unlist(lapply(ratios_lst, function(x) unlist(lapply(x[[y]]$freq, function(z) z$p.value)))))))


binom_test_out = binom_test_0.3 = list()
for(r in unique(data_kmers_list_df$RRIX)){
  keep_genes = data_kmers_list_df %>% filter(RRIX == r) %>% group_by(seq.Gene, CC_lab) %>%
    tally() %>% group_by(seq.Gene) %>% tally() %>% filter(n == 2)
  tmp = data_kmers_list_df %>% filter(RRIX == r) %>% group_by(seq.Gene, CC_lab) %>%
    filter(seq.Gene %in% keep_genes$seq.Gene) %>%
    summarize(mat_counts = as.numeric(sum(CC_1)),
              pat_counts = as.numeric(sum(CC_2)), 
              tot_counts = as.numeric(sum(sum)),
              ratio      = as.numeric(mat_counts / tot_counts))
  #binom_test_out[[r]] = apply(tmp, 1, function(x) 
  binom_test_out[[r]] = list()
  #binom_test_0.3[[r]] = list()
    for(i in 1:nrow(tmp)){
      x=tmp[i,]
      binom_test_out[[r]][[tmp$seq.Gene[i]]] = binom.test(x=unlist(x["mat_counts"]), n=unlist(x["tot_counts"]))
      #if(x["mat_counts"] <= x["pat_counts"]){
      #  binom_test_0.3[[r]][[tmp$seq.Gene[i]]] = binom.test(x=unlist(x["mat_counts"]), n=unlist(x["tot_counts"]), 
      #                                                      p=0.3, alternative="less")
      #} else {
      #  binom_test_0.3[[r]][[tmp$seq.Gene[i]]] = binom.test(x=unlist(x["mat_counts"]), n=unlist(x["tot_counts"]), 
      #                                                      p=0.7, alternative="greater")
      #}
    }
}

lapply(binom_test_out, function(x) print(hist(unlist(lapply(x, function(y) y$p.value)))))

manhattan.plot<-function(chr, pos, pvalue,
                         sig.level=NA, annotate=NULL, ann.default=list(), 
                         should.thin=T, thin.pos.places=2, thin.logp.places=2, 
                         main=NULL, xlab="Chromosome", ylab=expression(-log[10](p-value)),
                         region_list=NULL, 
                         col=c("gray","darkgray"), panel.extra=NULL, pch=20, cex=0.8,...) {
  data_kmers_list_plot = list()
  for(i in 1:length(data_kmers_list)){
    tmp=data.frame(data_kmers_list[[i]] %>% group_by(seq.Chromosome, seq.Gene, Pup.ID, 
                                                 CCs, CC.1, CC.2, CC1_hap, CC2_hap, 
                                                 Breeding.Batch,Behavior.Batch,         
                                                 RIX,Reciprocal,Diet) %>%
      summarize(pos = mean(as.numeric(seq.Position)), 
                mat_count = sum(CC_1), pat_count = sum(CC_2),
                l_count = sum(l_count), r_count = sum(r_count),
                total_count=sum(mat_count,pat_count)))
    
    tmp_list = list()
    tmp$fac=NA
    for(j in 1:nrow(tmp)){
      x = tmp[j,]
      fac = bwt_sizes$bwt_size[which(bwt_sizes$Pup.ID == as.numeric(x["Pup.ID"]))]/1e6
      tmp[j, grep("count", colnames(tmp))] = tmp[j, grep("count", colnames(tmp))]/fac
      tmp$fac[j] = fac
    }
    data_kmers_list_plot[[i]] = tmp
  }
}

#data_kmers_plot_df[,grep("count", colnames(data_kmers_plot_df))] = data_kmers_plot_df[,grep("count", colnames(data_kmers_plot_df))]*10

data_kmers_plot_df = do.call("rbind", data_kmers_list_plot)
check = data_kmers_plot_df %>% filter(seq.Gene %in% c("Ndn", "Ube3a","Snhg14","Snrpn",
                              "RP24-274B19.1","RP23-54G8.4","RP23-54G8.1","Magel")) %>%
  gather(key="count_type", value="count", "mat_count":"total_count")

head(check)
#write.csv(data_kmers_plot_df, file.path(dir, "data_kmers_normed_plot_df.csv"))
data_kmers_plot_df = read.csv(file.path(dir, "data_kmers_normed_plot_df.csv"))

#saveRDS(data_kmers_list, file.path(dir, "trec/data_kmers_from_process_and_plot/data_kmers_list_out_29jan2021.rds"))

ratios_lst = readRDS(file.path(dir, "trec/data_kmers_from_process_and_plot/ratios_lst_out_29jan2021.rds"))

total_counts = lapply(data_kmers_list, function(x) 
  x %>% group_by(seq.Gene, seq.Chromosome, Pup.ID, Reciprocal, Diet) %>%
    summarize(sum_1 = sum(CC_1), sum_2 = sum(CC_2), sum_tot = sum(sum), 
              pos = mean(as.numeric(seq.Position)),
              ratio = sum_1 / sum_tot))

total_counts = do.call("rbind", total_counts)
total_counts$stringtie_cts = NA
total_counts$RRIX = gsub("[a-z]","",total_counts$Reciprocal)
total_counts$PO   = gsub("[0-9]","",total_counts$Reciprocal)



ase_total_counts = matrix(NA, ncol=length(colData$ID), nrow=nrow(gene_count))
colnames(ase_total_counts) = colnames(gene_count)
rownames(ase_total_counts) = rownames(gene_count)

for(p in colnames(ase_total_counts)){
  tmp = total_counts %>% filter(Pup.ID == gsub("Pup.ID_","",p))
  matches = match(tmp$seq.Gene,rownames(ase_total_counts))
  ase_total_counts[matches[!is.na(matches)],which(colnames(ase_total_counts) == p)] = tmp$sum_tot[!is.na(matches)]
}

ase_total_counts = read.csv(file.path(dir,"trec/ase_total_count_matrix.csv"))
rownames(ase_total_counts) = ase_total_counts$X
ase_total_counts$X = NULL

for(i in grep("Pup.ID", colnames(gene_count))){
  coln = grep("Pup.ID", colnames(gene_count))[i]
  pup = unlist(strsplit(colnames(gene_count)[coln], "_"))[c(F,T)]
  tmp = total_counts[which(total_counts$Pup.ID == pup),]
  #gene_count[match(tmp$seq.Gene, rownames(gene_count)), coln]
  total_counts$stringtie_cts[which(total_counts$Pup.ID == pup)] = gene_count[match(tmp$seq.Gene, rownames(gene_count)), coln]
}

plot_tot_counts = total_counts %>% filter(!is.na(stringtie_cts)) %>%
  group_by(Pup.ID) %>% mutate(norm_tot = sum_tot / sd(sum_tot), 
                              norm_str = stringtie_cts / sd(stringtie_cts),
                              log_norm_tot = log(norm_tot+0.5), 
                              log_norm_str = log(norm_str+0.5),
                              ratio = norm_str / (norm_tot + norm_str))
p_list = list()
for(i in unique(plot_tot_counts$seq.Chromosome)){
  pdf = plot_tot_counts %>% filter(seq.Chromosome == i)
  p_list[[i]] = ggplot(pdf, aes(x=pos)) + 
    geom_line(aes(y=log_norm_tot), col="blue") + 
    geom_line(aes(y=log_norm_str), col="red") + 
    #geom_line(aes(y=ratio), col="blue") + 
    #ylim(0, 5) + 
    theme_classic() + 
    facet_wrap( ~ RRIX)
}


counts_per_pup = lapply(unique(total_counts$Pup.ID), function(x)
  total_counts %>% filter(Pup.ID == x, !is.na(stringtie_cts)))

##############  fisher  ##################


binom_test_pvals = lapply(ratios_lst, function(x){   
  adj = do.call("rbind", sapply(1:length(x), function(y) {
    tmp = data.frame(do.call("rbind", lapply(x[[y]]$freq, function(z) 
      c(z$statistic, z$parameter, z$p.value, as.vector(z$conf.int), z$estimate))))
    colnames(tmp) = c("n.success", "n.trial", "p.value", "lower", "upper", "est")
    tmp$Gene.Name = rownames(tmp)
    tmp$RIX = names(x)[[y]]
    tmp
  }, simplify=F)) 
  adj$padj = as.numeric(p.adjust(adj$p.value, method = "BH"))
  adj
})


ratios_lst = binom_test_out

binom_test_pvals = lapply(ratios_lst, function(x){   
  tmp = data.frame(do.call("rbind", lapply(x, function(z) 
    c(z$statistic, z$parameter, z$p.value, as.vector(z$conf.int), z$estimate))))
  colnames(tmp) = c("n.success", "n.trial", "p.value", "lower", "upper", "est")
  tmp$Gene.Name = rownames(tmp)
  tmp$padj = as.numeric(p.adjust(tmp$p.value, method = "BH"))
  zscore_pos = 1-qnorm(tmp$p.value)
  zscore_neg =  ( tmp$est-0.5 )  / sqrt((0.5*0.5) / tmp$n.trial)
  
  tmp$zscore = ifelse(tmp$est < 0.5, zscore_pos, -zscore_pos) #qnorm(tmp$p.value)#ifelse(tmp$est < 0.5, zscore_neg, zscore_pos)
  tmp$zscore = zscore_neg
  tmp$var = (tmp$est*(1-tmp$est)) / tmp$n.trial  #tmp$est * (1-tmp$est)
  tmp
}) 
names(binom_test_pvals) = names(ratios_lst)

binom_test_pvals = do.call("rbind", binom_test_pvals)
binom_test_pvals$RIX = unlist(lapply(strsplit(rownames(binom_test_pvals),"[.]"), function(x) x[1]))
binom_test_pvals$offset = abs(binom_test_pvals$est - 0.5)
binom_test_pvals_sig = binom_test_pvals %>% filter(padj < 0.1)
fdrz = fdrtool(binom_test_pvals$zscore, statistic="normal")
binom_test_pvals$fdrp = fdrz$pval
pval_list = lapply(unique(binom_test_pvals$Gene.Name), function(x) binom_test_pvals %>% filter(Gene.Name == x))
names(pval_list) = unique(binom_test_pvals$Gene.Name)
pmat = data.frame(apply(do.call("rbind", lapply(pval_list, function(g){
  weights = sqrt(g$n.trial) / sum(sqrt(g$n.trial))
  data.frame(souffer_stat = sum(g$zscore * weights) / sum(weights^2), n=length(weights))
})), 2, as.numeric))
fdrp = fdrtool(pmat$souffer_stat)
pmat$p_fromz = pnorm(abs(pmat$souffer_stat), lower.tail = F)
pmat$p = fdrp$pval
pmat$qval = fdrp$qval
pmat$lfdr = fdrp$lfdr
fisher = sumlog(pmat$fdrp)
genes = rownames(pmat %>% filter(qval<0.1))# %>% arrange(desc(weighted_z)))


pmat_fish = do.call("rbind", lapply(pval_list, function(g){
  data.frame(stat = -2*sum(log(g$fdrp)), df = 2*nrow(g), p = pchisq(-2*sum(log(g$fdrp)), df=2*nrow(g), lower.tail = F))
}))
fdr_fish = fdrtool(pmat_fish$p, statistic ="pvalue")
#colnames(pmat) = c("fisher_stat", "n")
#rownames(pmat) = unique(binom_test_pvals$Gene.Name)

#pmat = data.frame(do.call("rbind", lapply(pval_list, function(g){
#  data.frame(weighted_z = sum(g$zscore*(1/g$var)) / sqrt(sum((1/g$var)^2)), n_rix = nrow(g))
#})))
#pmat=do.call("rbind",pmat)
#colnames(pmat) = c("weight_stouffer_stat", "n")
rownames(pmat) = unique(binom_test_pvals$Gene.Name)
pmat$stouffer_p = unlist(lapply(1:nrow(pmat), function(x)
  pnorm(pmat$weighted_z[x])))
pmat %>% arrange(desc(abs(weighted_z)))
pmat %>% arrange(stouffer_p)

hist(pmat$souffer_stat, breaks=2000, xlim=c(-10,10))
curve(dnorm(x)*nrow(pmat), from=-5, to=5, col="blue", add=T)

curve(dchisq(x, df = length(unique(binom_test_pvals$RIX))*2)*nrow(pmat), from=0, to=150, col="blue", add=T)

pmat$fisher_p = unlist(lapply(1:nrow(pmat), function(x)
  pchisq(pmat$fisher_stat[x], df=(as.numeric(paste(pmat$n[x]))*2),lower.tail = F)))

fdr = fdrtool(x = pmat$fisher_p, statistic = "pvalue")
pmat$padj = as.numeric(p.adjust(pmat$fisher_p, method = "BH"))

pmat$padj_qval = as.numeric(fdr$qval)
pmat$padj_fdr  = as.numeric(fdr$lfdr)


pmat$ie    = ifelse(rownames(pmat) %in% ie_genes$mgi_symbol, T, F)
pmat$gregg = ifelse(rownames(pmat) %in% gregg_genes$Gene.Name, T, F)
pmat$de    = ifelse(rownames(pmat) %in% de_genes$gene, T, F)
pmat$de    = ifelse(rownames(pmat) %in% genes, T, F)

pmat %>% arrange(desc(abs(fisher_stat)))

dim(pmat %>% filter(ie, padj_fdr<0.1) %>% arrange(desc(fisher_stat)))
pmat %>% filter(!ie, padj_fdr<0.05, de) %>% arrange(desc(fisher_stat))

## Gnas ### Wars, Meg3 (not anymore)
## 100 ie but not de
## 114 de but not ie

genes = rownames(pmat %>% filter(padj_fdr < 1/length(unique(total_counts$seq.Gene)), !ie, n > 2) %>% arrange(padj_fdr))
genes = rownames(pmat %>% filter(padj_fdr < 1e-4, !ie, de) %>% arrange(desc(fisher_stat)))
genes = rownames(pmat %>% filter(!ie, gregg) %>% arrange(desc(fisher_stat)))
genes = rownames(pmat %>% filter(!ie, padj_qval<0.05) %>% arrange(desc(fisher_stat)))
genes = rownames(pmat %>% filter(padj_qval<0.1))# %>% arrange(desc(weighted_z)))
genes = rownames(pmat %>% filter(de) %>% arrange(padj_qval))

hi_rat = data_kmers_plot_df %>% ##ungroup() %>%
  group_by(seq.Gene) %>%
  #summarize(sum_1 = sum(mat_count), sum_tot = sum(total_count),n = n()) %>%
  mutate(ratio_m = mat_count / total_count) %>%
  summarize(mean_rat = mean(ratio_m), abs_rat = abs(mean_rat - 0.5),
            n=n(), mean_counts = mean(total_count), max_count = max(total_count)) %>%
  arrange(desc(abs_rat)) %>%
  filter(abs_rat > 0.150, mean_rat < 1, n > 2, max_count > 10)
genes = hi_rat$seq.Gene[which(hi_rat$seq.Gene %in% genes)]
genes = genes[which(genes %in% hi_rat$seq.Gene)]
a = length(which(genes %in% ie_genes$mgi_symbol)) + 1 
b = length(genes)-a
c = length(ie_genes$mgi_symbol) - a
d = dim(pmat %>% filter(padj_qval < 0.1))[1] - c - a - b

fisher.test(matrix(c(a,b,c,d), nrow=2, byrow = T))
geneUse = "Eef2"
geneUse = "Inpp5f"
geneUse = "Igf2"
geneUse = "Gatad1"
geneUse = "Inpp5f"

p_list = list()
use_genes = data.frame(old = c("Ube3a", "Snhg14", "Snrpn", "RP24-274B19.1",
                               "RP23-54G8.4", "RP23-54G8.1","Ndn","Magel2"),
                       new = c("Ube3a", "Ube3a-ATS", "Snrpn", "A230057D06Rik",
                               "Gm2357","A330076H08Rik","Ndn","Magel2"))
#use_genes$old
for(g in genes[1:100]){
  geneUse = g
  geneNew = use_genes$new[which(use_genes$old == geneUse)]
  plot_gene = data_kmers_plot_df %>% filter(seq.Gene == geneUse) %>%
    left_join(lab[,-which(colnames(lab) %in% c("RIX", "Reciprocal", "Diet","CCs","CC.1","CC.2"))], by="Pup.ID") %>%
    mutate(ratio = mat_count / total_count,
           PO = gsub("[0-9]", "", Reciprocal))
  check_recip = table(plot_gene %>% select(Reciprocal, CCs) %>% distinct())
  rem_CC = names(which(colSums(check_recip) < 2))
  if(length(rem_CC) > 0){
    plot_gene = plot_gene%>% filter(!CCs %in% rem_CC)
  }
  if(nrow(plot_gene)>0){
    plot_gene_long = gather(plot_gene, data_type, value, mat_count:total_count, factor_key=TRUE) %>% 
      filter(data_type %in% c("mat_count","pat_count","total_count")) %>%
      mutate(DamLine_NewCC_ID = gsub("CC0","",DamLine_NewCC_ID),
             CCs = gsub("CC0","",CCs))
    jitter <- position_jitter(width = 0.1, height = 0)
    jitter2 <- position_jitter(width = 0.2, height = 0)
    
    colors <- c("Maternal" = "red", "Paternal" = "blue", "Total"="purple")
    colors <- c("firebrick2","royalblue","purple")
    
    p_list[[g]] = ggplot(data=plot_gene_long, aes(x=DamLine_NewCC_ID)) +   
      geom_boxplot(data=plot_gene_long %>% filter(data_type != "total_count"), 
                   aes(y=value, fill=data_type)) + 
      geom_point(data=plot_gene_long %>% filter(data_type == "total_count"), 
                 aes(y=value, col=data_type), size = 3, alpha = 0.3, position=jitter) + 
      facet_grid(~ CCs, scales="free_x") + 
      scale_fill_manual(values=colors, name = "Allele-specific",
                        labels=c("Maternal","Paternal")) + 
      scale_color_manual(name="Total \nexpression", labels = "", values="purple") +
      labs(title = geneNew, y = "Counts", x="Maternal strain") +
      theme_classic()
  }
}



pdf(file.path(dir,'trec/ase_plots_3mar2021.pdf'))
p_list
dev.off()

pdf(file.path(dir,'trec/ase_DE_gene_plots.pdf'))
p_list
dev.off()

pdf(file.path(dir,'trec/ase_IE_gene_plots.pdf'))
p_list
dev.off()

#Fbn1, Lrrc48, Ndnf
#p_list[[g]] = ggplot(plot_gene, aes(x=DamLine_NewCC_ID)) +   
#  geom_point(aes(y=sum_tot,col="Total"), size = 3, alpha = 0.3, position=jitter2) + 
#  geom_point(aes(y=sum_1,col="Maternal"), position=jitter) + 
#  geom_point(aes(y=sum_2,col="Paternal"), position=jitter) + 
#geom_point(aes(y = ratio)) +
#  facet_grid(~ CCs, scales="free_x") + 
#  labs(title = geneUse, color = "Expression",
#       y = "Counts", x="Maternal strain") +
#  scale_color_manual(values = colors) + 
#  theme_classic()


fdr = fdrtool(x = pmat$fisher_p, statistic = "pvalue")
pmat$padj_qval = fdr$qval
pmat$padj_fdr = fdr$pval
fdr = fdrtool(x = pmat$praw, statistic = "pvalue")
pmat$praw_qval = fdr$qval
pmat$praw_fdr = fdr$pval
pmat$padj_bh = p.adjust(pmat$padj_fdr, method = "BH")
pmat$praw_bh = p.adjust(pmat$praw_fdr, method = "BH")


  j=1
  resStan_tmp = ratios_lst[[j]]
  sum_df = summary(resStan_tmp)$summary
  simp_res = do.call("rbind", sapply(1:length(ratios_lst_simp), function(i){
    tmp = data.frame(summary(ratios_lst_simp[[i]])$summary)
    tmp$rix = names(ratios_lst_simp)[i]
    tmp
    }, simplify=F)
  )
  sum_df[-grep("ind|eta|_a|_b|weight", rownames(sum_df)),]
  
  #saveRDS(data_kmers, 
  #        file.path(dir, paste0("/regression_outputs/chr",i,"_data_mnt_29dec2020.rds")))
#} else {
#  data_kmers = readRDS(file.path(dir, paste0("/regression_outputs/chr",i,"_data_mnt_29dec2020.rds")))
#}
  
  



