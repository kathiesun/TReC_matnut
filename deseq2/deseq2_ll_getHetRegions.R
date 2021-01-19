options(stringsAsFactors = FALSE)
library(tidyverse)

#setwd("C:/Users/Kathie/rna_seq/kmerSearch")
#dir <- "C:/Users/Kathie/Dropbox\ (ValdarLab)/variant_data"

setwd("~/TReC_matnut/src")
dir="/nas/depts/006/valdar-lab/users/sunk/"
source("summary_functions.R")

args <- commandArgs(trailingOnly = TRUE)
p = as.numeric(args[1])
#chr = c(1:19, "X")
#c = chr[i]


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
length(which(!gene_count$Gene.Name %in% annot_genes$Gene.Name))
gene_count = gene_count %>% filter(Gene.Name %in% annot_genes$Gene.Name)

if(any(duplicated(gene_count$Gene.Name))) gene_count = data.frame(gene_count[-which(duplicated(gene_count$Gene.Name)),])
gene_annot = left_join(gene_count, annot_genes, by="Gene.Name")

indiv_pups = list.files(file.path(dir, "mini/pup_haplo_blocks_by_CC_parent_dec2019"), pattern="haploBlocks", full.names = T)
phased_CC_haplotype = lapply(indiv_pups, readRDS)
tmp = do.call("rbind", lapply(indiv_pups, function(x) unlist(strsplit(x, "_"))))
names(phased_CC_haplotype) = paste0("Pup.ID_", tmp[,ncol(tmp)-1])

#lapply(names(phased_CC_haplotype), function(x) {
#  mat = phased_CC_haplotype[[x]]$founder_by_parent[[1]]
#  pat = phased_CC_haplotype[[x]]$founder_by_parent[[2]]
#  apply(gene_annot, 1, function(y) {
#    mat_tmp = mat %>% filter(as.character(chr) == as.character(y["Chr"]), start < y["Start"], end > y["End"])
#    pat_tmp = pat %>% filter(as.character(chr) == as.character(y["Chr"]), start < y["Start"], end > y["End"])
#    ifelse(mat_tmp$found == pat_tmp$found, NA, 
#           ifelse((nchar(mat_tmp$found) != 1 | nchar(pat_tmp$found) != 1), NA, y[[x]]))
#  })
#})


#het_counts = list()
#for (x in names(phased_CC_haplotype)[p]) {
x = names(phased_CC_haplotype)[p]
  het_counts = rep(NA, length=nrow(gene_annot))
  mat = phased_CC_haplotype[[x]]$founder_by_parent[[1]]
  pat = phased_CC_haplotype[[x]]$founder_by_parent[[2]]
  for(i in 1:nrow(gene_annot)){
    y=gene_annot[i,]
    mat_tmp = mat %>% filter(as.character(chr) == as.character(y["Chr"]), start < y["Start"], end > y["End"])
    pat_tmp = pat %>% filter(as.character(chr) == as.character(y["Chr"]), start < y["Start"], end > y["End"])
    if(nrow(mat_tmp) > 0 & nrow(pat_tmp) > 0){
      if(mat_tmp$found != pat_tmp$found & (nchar(mat_tmp$found) == 1 | nchar(pat_tmp$found) == 1)){
        het_counts[i] = y[[x]]
      }
    }
    
    #if(i%%10 == 0) print(paste(x,i))
  }
#}

#het_counts = data.frame(do.call("cbind", het_counts))


write.table(het_counts, paste0(dir,"/trec/geneCounts_", x, "_hetRegions.txt"),
            row.names=F, col.names = F )
