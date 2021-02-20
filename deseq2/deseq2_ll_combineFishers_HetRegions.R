setwd("~/TReC_matnut/src")
library(tidyverse)
library(fdrtool)
library(DESeq2)

source("prediction_functions.R")
source("summary_functions.R")
source("ase_summary_source.R")
source("stan_pheno_functions.R")

args <- commandArgs(trailingOnly = TRUE)
c = as.numeric(args[1])
n=1000

dir <- "/nas/depts/006/valdar-lab/users/sunk/"
annot = read.table(file.path(dir, "mm10/ref_sequence/Mus_musculus.GRCm38.96.gtf"), skip=5, fill=T)
annot_genes = annot %>% filter(V15 == "gene_name") %>%
	dplyr::select(one_of(paste0("V",c(1,4,5,7,10,16)))) %>%
	distinct()

colnames(annot_genes) = c("Chr","Start","End","Strand","Gene.ID","Gene.Name")
annot_genes = annot_genes %>% mutate(Start = as.numeric(paste(Start)), End = as.numeric(paste(End)))
    # 
if(length(which(duplicated(annot_genes$Gene.Name))) > 0){
	annot_genes = annot_genes[-which(duplicated(annot_genes$Gene.Name)),]
}


indiv_pups = list.files(file.path(dir, "mini/pup_haplo_blocks_by_CC_parent_jan2021"), pattern="haploBlocks", full.names = T)

phased_CC_haplotype = lapply(indiv_pups, readRDS)
tmp = do.call("rbind", lapply(indiv_pups, function(x) unlist(strsplit(x, "_"))))
names(phased_CC_haplotype) = paste0("Pup.ID_", tmp[,ncol(tmp)-1])

dds_list = readRDS(file.path(dir,"de_results/dds_list_allReg_sepRIX_30jan2021.rds"))

POres_list = lapply(dds_list, function(x) results(x, name="PO"))
POres_Ordered <- lapply(POres_list, function(x) x[order(x$pvalue),])
all_genes = unique(unlist(lapply(POres_list, function(x) rownames(x))))
use_genes = all_genes[((c-1)*n + 1):min(length(all_genes), (c-1)*n + n)]

sig_genes_list = lapply(POres_list, function(x) x[which(rownames(x) %in% use_genes),])
#sig_genes_short = lapply(sig_genes_list, function(x) 
#  if(length(grep("Gm|Rik|[.]", x)) > 0) x[-grep("Gm|Rik|[.]", x)])
#deseq_sig = sort(table(unlist(lapply(POres_Ordered, function(x) rownames(x[which(x$padj < 0.1),])))), decreasing = T)

for(i in 1:length(sig_genes_list)){
  sig_genes_list[[i]] = data.frame(sig_genes_list[[i]])
  if(nrow(sig_genes_list[[i]]) > 0){
	sig_genes_list[[i]]$gene = rownames(sig_genes_list[[i]])
  	sig_genes_list[[i]]$rix = names(sig_genes_list)[i]
  	rownames(sig_genes_list[[i]]) = NULL
  }
}

sig_genes_df = do.call("rbind", sig_genes_list)

lab = read.csv(file.path(dir, "matnut_main/CC_labels_paternity.csv"))

which_homs_list = pmat_list = list()

for(g in unique(sig_genes_df$gene)){
  p_tmptmp = lapply(POres_list, function(x) x[which(rownames(x) == g),])
  p_tmp = do.call("rbind", p_tmptmp)
  p_tmp$gene = g
  p_tmp$RIX = names(p_tmptmp)[which(unlist(lapply(p_tmptmp, nrow)) > 0)]
  #chec_het = data.frame(counts = t(gene_count[which(rownames(gene_count) == g),]), 
  #                      Pup = colnames(gene_count_red))
  #chec_het$Pup.ID = as.numeric(unlist(strsplit(chec_het$Pup, "_"))[c(F,T)])
  gene_det = annot_genes %>% filter(Gene.Name == g)
  if(gene_det$Chr %in% c(1:19)){
    pup_haps = do.call("rbind", lapply(phased_CC_haplotype, function(x) {  
      tmp = x$full_founders %>% filter(chr == gene_det$Chr, start < gene_det$Start, end > gene_det$End) 
      if(!names(x$founder_by_parent)[1] %in% tmp$cc){
        tmp1 = x$full_founders %>% filter(chr == gene_det$Chr, cc == names(x$founder_by_parent)[1])
        start = tmp1[which.min(abs(tmp1$start - gene_det$Start)),]
        end   = tmp1[which.min(abs(tmp1$end - gene_det$End)),]
        if(abs(start[,"start"] - gene_det$End) < abs(end[,"end"] - gene_det$Start)){
          tmp[1,] = start
        } else {
          tmp[1,] = end
        }
      }
      if(!names(x$founder_by_parent)[2] %in% tmp$cc){
        tmp2 = x$full_founders %>% filter(chr == gene_det$Chr, cc == names(x$founder_by_parent)[2])
        start = tmp2[which.min(abs(tmp2$start - gene_det$Start)),]
        end   = tmp2[which.min(abs(tmp2$end - gene_det$End)),]
        if(abs(start[,"start"] - gene_det$End) < abs(end[,"end"] - gene_det$Start)){
          tmp[2,] = start
        } else {
          tmp[2,] = end
        }
      }
      tmp$group = as.numeric(tmp$group)
      if("N" %in% tmp$found[1]){
        newfound1 = x$full_founders %>% filter(cc == tmp$cc[1], group %in% c(tmp$group[1]-1,tmp$group[1]+1))
        tmp$found[1] = paste(unique(newfound1$found), collapse="_")
      } 
      if("N" %in% tmp$found[2]){
        newfound2 = x$full_founders %>% filter(cc == tmp$cc[2], group %in% c(tmp$group[2]-1,tmp$group[2]+1))
        tmp$found[2] = paste(unique(newfound2$found), collapse="/")
      }
      c(founders = paste(tmp$found, collapse="/"), CCs = paste(tmp$cc, collapse="/"), Pup.ID = tmp$pup[1])
    }))
    lab_to_rix = unique(lab[,c("RIX","CCs")])
    which_homs = data.frame(unique(pup_haps[,c("founders","CCs")]))
    which_homs$hom = unlist(lapply(which_homs$founders, function(y) any(table(unlist(strsplit(y, "_|/"))) > 1)))
    which_homs$RIX = as.character(lab_to_rix$RIX[match(which_homs$CCs, lab_to_rix$CCs)])
    which_homs = which_homs %>% group_by(RIX, CCs) %>%
      summarize(hom_tot = ifelse(any(hom),T,F),
                founders_all = paste(founders,collapse=":"))
    p_tmp = data.frame(p_tmp) %>% left_join(which_homs, by="RIX")
    pmat = data.frame(het_stat = -2*sum(log(p_tmp$pvalue[which(!p_tmp$hom_tot)])), 
                      het_n    = length(which(!p_tmp$hom_tot)),
                      hom_stat = -2*sum(log(p_tmp$pvalue[which(p_tmp$hom_tot)])), 
                      hom_n    = ifelse(length(which(p_tmp$hom_tot)) > 0, length(which(p_tmp$hom_tot)), 0)
    )
    rownames(pmat) = g
    
    which_homs_list[[g]] = which_homs
    pmat_list[[g]] = pmat
  }
  
}


saveRDS(which_homs_list, paste0(dir,"/de_results/het_regions/which_homs_list_",c,".rds"))
saveRDS(pmat_list, paste0(dir,"/de_results/het_regions/pmat_list_",c,".rds"))

