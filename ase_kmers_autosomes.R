setwd("C:/Users/Kathie/TReC_matnut/src")
library(tidyverse)
library(fdrtool)

source("prediction_functions.R")
source("summary_functions.R")
source("ase_summary_source.R")
source("stan_pheno_functions.R")


dir <- "C:/Users/Kathie/Dropbox\ (ValdarLab)"
gecco <- read.csv(file.path(dir, "de_results/fullGeccoRnaDump.csv"))

C_diet_4 = read.csv(file.path(dir, "variant_data/C_matrix_4diets.csv"), header=F)
C_diet_2 = read.csv(file.path(dir, "variant_data/C_matrix_2diets.csv"), header=F)
C_diet_3 = matrix(c(0.9659258,-0.2588190,-0.7071068,-0.2588190,0.9659258,-0.7071068),nrow=3,ncol=2)


matnut = read.csv(file.path(dir,'matnut_main/AllMice_GeneExpression_SSupdated_11.27.19.csv'))
gene_count <- read.csv(file.path(dir,'/trec/gene_count_matrix.csv'))

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

de_genes = read.table(file.path(dir, "trec/priority_deseq_genes_11jan2021.txt"))
ase_genes = read.csv(file.path(dir, "trec/priority_ase_genes_11jan2021.csv"))

all_genes = unique(c(gregg_genes$Gene.Name, de_genes$V1, ase_genes$Gene.Name, ie_genes$mgi_symbol))

################################################################


indiv_pups = list.files(file.path(dir, "mini/pup_haplo_blocks_by_CC_parent_dec2019"), pattern="haploBlocks", full.names = T)

phased_CC_haplotype = lapply(indiv_pups, readRDS)
tmp = do.call("rbind", lapply(indiv_pups, function(x) unlist(strsplit(x, "_"))))
names(phased_CC_haplotype) = paste0("Pup.ID_", tmp[,ncol(tmp)-1])

seg_regions = readRDS(file.path(dir,"mini/seg_regions_by_genotyping_perRIX_kmerBased_23nov2019.rds"))
masterSnps_files = list.files(file.path(dir, "variant_data/kmers_to_run"), pattern="masterSnps_chr", full.names = T)
count_files <- list.files(file.path(dir,"variant_data/kmer_counts"), ".csv", full.names = T)

## CC labels
#lab = read.csv("C:/Users/Kathie/Dropbox (ValdarLab)/variant_data/matched_v2_4jun2018.csv")
lab = read.csv(file.path(dir, "matnut_main/matched_v2_4jun2018.csv"))
lab %>% dplyr::select(one_of("Pup.ID","RRIX", "CC.1","CC.2")) %>%
  filter(!is.na(Pup.ID)) -> lab

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
  useSnps = masterSnps %>% filter(seq.Gene %in% all_genes)
  
  ## count data 
  chr_files <- count_files[grep(paste0("chr", c,"_"), count_files)]
  ref_file <- chr_files[grep("ref", chr_files)]
  alt_file <- chr_files[grep("alt", chr_files)]
  reforig <- unique(read.csv(ref_file, header = T))
  altorig <- unique(read.csv(alt_file, header = T))
  reforig$X = altorig$X <- NULL
  reforig = reforig %>% filter(k.mer %in% useSnps$seq.refseq)
  altorig = altorig %>% filter(k.mer %in% useSnps$seq.altseq)
  reforig$pup = reforig$pup.id 
  altorig$pup = altorig$pup.id 
  
  data_kmers[[c]] = process_and_plot(chr=c, 
                                snp_info=useSnps, 
                                sample_info=pupInfo, 
                                RIX_info=lab, 
                                ref_counts=reforig, alt_counts=altorig, 
                                phased_CC_haplotype=phased_CC_haplotype, 
                                use_gene=F,
                                problemPups = c(1404, 1716, 1371, 569, 1911, 1951, 1015),
                                seg_regions=seg_regions)
 
  #data_kmers$podietrix = paste(data_kmers$dir, data_kmers$DietRIX, sep="_")
  print(c)
  ratios_lst[[c]] = run_stan_regress(data_kmers=data_kmers[[c]], 
                                niter=10000, n.thin=5,  
                                seg_regions=seg_regions,
                                save_dir=NULL, 
                                STZ=T, use_gene=F,
                                no_theta=F, alpha=NULL,
                                stan=F, stanMod = "ase_mu_g_simple.stan")
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
pval_comb = fishercomb(pval_list)

##############  fisher  ##################
pmat = do.call("rbind", lapply(pval_list, function(g){
  data.frame(-2*sum(log(g$x.p.value)), nrow(g))
}))

pmat=do.call("rbind",pmat)
colnames(pmat) = c("fisher_stat", "n")
pmat$gene = unique(pvals$gene)

hist(pmat$fisher_stat, breaks=900)
curve(dchisq(x, df = length(unique(pvals$rix))*2)*nrow(pmat), from=0, to=150, col="blue", add=T)

pmat$fisher_p = unlist(lapply(1:nrow(pmat), function(x)
  pchisq(pmat$fisher_stat[x], df=(as.numeric(paste(pmat$n[x]))*2),lower.tail = F)))

fdr = fdrtool(x = pmat$fisher_p, statistic = "pvalue")
pmat$padj_qval = fdr$qval
pmat$padj_fdr  = fdr$lfdr
pmat$padj = p.adjust(pmat$fisher_p, method = "BH")

pmat$ie    = ifelse(pmat$gene %in% ie_genes$mgi_symbol, T, F)
pmat$gregg = ifelse(pmat$gene %in% gregg_genes$Gene.Name, T, F)
pmat$de    = ifelse(pmat$gene %in% de_genes$V1, T, F)
pmat %>% arrange(fisher_p)

pmat %>% filter(!ie, de) %>% arrange(fisher_p)
## Gnas, Wars, Meg3
## 100 ie but not de
## 114 de but not ie




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
