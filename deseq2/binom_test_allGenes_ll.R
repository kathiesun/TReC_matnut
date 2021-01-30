setwd("~/TReC_matnut/src")
library(tidyverse)
library(fdrtool)

source("prediction_functions.R")
source("summary_functions.R")
source("ase_summary_source.R")
source("stan_pheno_functions.R")

args <- commandArgs(trailingOnly = TRUE)
c = as.numeric(args[1])

dir <- "/nas/depts/006/valdar-lab/users/sunk/"

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

#ie_genes = read.csv("C:/Users/Kathie/Dropbox (ValdarLab)/imprinted_genes/2014_05_14_allImprintedGenes.csv")
#gregg_genes = read.csv("../gregg_POgenes.csv")
#colnames(gregg_genes)[1] = "Gene.Name"
#de_genes = read.table(file.path(dir, "trec/priority_deseq_genes_11jan2021.txt"))
#ase_genes = read.csv(file.path(dir, "trec/priority_ase_genes_11jan2021.csv"))
#all_genes = unique(c(gregg_genes$Gene.Name, de_genes$V1, ase_genes$Gene.Name, ie_genes$mgi_symbol))

################################################################

indiv_pups = list.files(file.path(dir, "mini/jan2021_pup_haplo_blocks_by_CC_parent"), pattern="haploBlocks", full.names = T)

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

#for(c in 1:19){
  
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
  write.table(data_kmers[[c]], paste0(dir, "//ase_autosomes/data_kmers_from_process_and_plot/chr_",c,"_29jan2021.txt"), quote=F, row.names=F)

 # ratios_lst[[c]] = run_stan_regress(data_kmers=data_kmers[[c]], 
 #                                    niter=10000, n.thin=5,  
 #                                    seg_regions=seg_regions,
 #                                    save_dir=NULL, 
 #                                    STZ=T, use_gene=F,
 #                                    no_theta=F, alpha=NULL,
 #                                    stan=F, stanMod = "ase_mu_g_simple.stan")
 # print(paste(c, "done"))
  
#}

