library(rjags)
library(coda)
library(ggplot2)
library(RColorBrewer)
#library(ggmap)
library(data.table)
library(lmerTest)
library(tidyverse)
setwd("~/matnut/src")
#setwd("C:/Users/Kathie/matnut/src")
args <- commandArgs(trailingOnly = TRUE) 
it = as.numeric(args[1])

source("./matnut/matching_functions_lump.R")

#source("./matnut/lmer_functions_rna.R")
#source("./matnut/jags_functions.R")
#source("./matnut/boxcox_functions.R")
#source("./matnut/prediction_functions.R")
#source("./matnut/summary_functions.R")

###### Read in data
#dir <- "C:/Users/Kathie/Dropbox\ (ValdarLab)"
dir = "/nas/depts/006/valdar-lab/users/sunk"

matnut = read.csv(file.path(dir,'matnut_main/AllMice_GeneExpression_SSupdated_11.27.19.csv'))
matnut = matnut %>% select(-contains("X"))
gene_count <- read.csv(file.path(dir,'/hisat2_stringtie_fasta/stringtie/round2/gene_count_matrix.csv'))
#gene_count <- read.csv(file.path(dir,'/trec/gene_count_matrix.csv'))

colnames(gene_count)[1] = "gene_id"
gene_count$Gene.ID = do.call("rbind",(strsplit(as.character(gene_count$gene_id), "[|]")))[,1]
gene_count$Gene.Name = do.call("rbind", (strsplit(as.character(gene_count$gene_id), "[|]")))[,2]

samples = colnames(gene_count)[grep("Pup", colnames(gene_count))]
matnut$ID = paste0("Pup.ID_",matnut$Pup.ID)
matnut$Diet = gsub(" $", "", matnut$Diet)
matnut$RIX = gsub("a|b","",matnut$Reciprocal)
matnut$PO = ifelse(gsub("[0-9]","",matnut$Reciprocal) == "a", 0.5, -0.5)
matnut$DietRIX = gsub(" ", ".",paste0(matnut$Diet, matnut$RIX))
matnut$DietRIXPOq = paste0(matnut$DietRIX,"_",matnut$PO)
matnut$RIX = factor(matnut$RIX, levels = c(1:4,6:10))
map = read.csv(file.path(dir,"mini/combined_map_cs_gbrs_allmarkers_17mar2020.csv"))

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


regions_list = readRDS(file.path(dir,"mini/seg_regions_by_genotyping_perRIX_kmerBased_23nov2019.rds"))
haplofiles = list.files(file.path(dir,"mini/pup_haplo_blocks_by_CC_parent_dec2019"), 
                        pattern="haploBlocks.rds", full.names = T)
phased_par_CC_haplotypes = lapply(haplofiles, readRDS)
tmp = do.call("rbind", strsplit(haplofiles, "_"))
names(phased_par_CC_haplotypes) = tmp[,ncol(tmp)-1]

maxn = floor(nrow(annot_genes)/500)
its_1 = ((it-1)*500)+1
its_2 = ifelse(it==maxn, nrow(annot_genes), it*500)

print(paste(its_1,its_2))
####################
#   Run function   #
####################


matched_results <- list()
### matching ###c(1,2,4,5,13,17,18)
for(g in its_1:its_2) {    #c(1,2,4,5,10,11,13,16,18)
  matnut_use = matnut[match(samples, matnut$ID),]
  matnut_use$gene = unlist(gene_count[g,samples])
  if(length(unique(matnut_use$gene)) > 5){
    if(length(-which(is.na(matnut_use$gene))) > 0){
      tab = sort(table(matnut[-which(is.na(matnut_use$gene)),"DietRIX"]))
    } else {
      tab = sort(table(matnut[,"DietRIX"]))
    }
    rem = unique(unlist(strsplit(names(tab)[which(tab < 5)],"_"))[c(T,F)])
    matnut_use = matnut_use %>% filter(!DietRIX %in% rem)
    enc = getEncoding(matnut_use, terms = c("RIX","Diet","PO"))
    gene = c(paste(gene_count$Gene.ID[g]), paste(gene_count$Gene.Name[g]))
    gene_det = annot_genes[which(annot_genes$Gene.ID %in% gene | annot_genes$Gene.Name %in% gene) ,]
    if(nrow(gene_det) > 0){
      haplos = do.call("rbind", lapply(sort(unique(matnut_use$RIX)), function(x){
        pups = (matnut_use %>% filter(RIX == paste(x)))$Pup.ID
        haps = unique(do.call("rbind", lapply(pups, function(y){
          mat = phased_par_CC_haplotypes[[paste(y)]]$founder_by_parent[[1]] %>% filter(chr == gene_det$Chr)
          mat_ind = which(mat$start < gene_det$Start & mat$end > gene_det$End)
          pat = phased_par_CC_haplotypes[[paste(y)]]$founder_by_parent[[2]] %>% 
            filter(chr == gene_det$Chr)
          pat_ind = which(pat$start < gene_det$Start & pat$end > gene_det$End)
          if(length(mat_ind) == 1 & length(pat_ind) == 1){
            c(mat$found[mat_ind], pat$found[pat_ind])
          } else {
            c("-","-")
          }
        })))
        apply(haps, 2, function(z) paste(unique(z), collapse="/"))
      }))
      
      hets = apply(haplos, 1, function(x) ifelse(any(strsplit(x, "/")[[1]] %in% strsplit(x, "/")[[2]]), 0, 1))
      segs = apply(haplos, 1, function(x) ifelse(any(unlist(strsplit(x, "/")) == "-"), NA, 0))
      homog = hets + segs
      homog = ifelse(homog == 0, 0, 1)

      try(matched_results[[gene[2]]] <- match.multimp(data=matnut_use, ptypes="gene", clust=homog, 
                                                      N=5,chains=1, n.iter=20000, encoded = enc, 
                                                      tryLam=c(-1, 0, .25, .33, .5, 1, 2, 3), thin=10,
                                                      randvar=NA, fixvar="Diet",  
                                                      matchon=c("RIX","Diet"), matchoff="PO", idcol="ID",
                                                      fixPO = F, fixClust = F, p="p", beta=T))
    }
  }
  #plot_rix = matnut_use$df[-which(is.na( matnut_use$df$ptypes[i]))]
}

saveRDS(matched_results, paste0(dir,"/trec/tot_expr_POE_",it,"_23nov2020.rds"))

