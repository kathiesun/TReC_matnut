library(qtl2) 
#library(qtl2convert)
library(tidyverse)
setwd("C:/Users/Kathie/TReC_matnut/src")


dir <- "C:/Users/Kathie/Dropbox\ (ValdarLab)"
dir <- "/nas/depts/006/valdar-lab/users/sunk"
source("C:/Users/Kathie/rna_seq/compare_geno/genotyping_compare_func.R")

problemPups = c(1404, 1716, 1371, 569, 1911, 1951, 1015)

matnut = read.csv(file.path(dir,'matnut_main/AllMice_GeneExpression_SSupdated_11.27.19.csv'))
matnut = matnut %>% filter(!Pup.ID %in% problemPups)
matnut$ID = paste0("Pup.ID_",matnut$Pup.ID)
matnut$Pup.ID = as.character(matnut$Pup.ID)
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
map = map[,-1]
use_map = lapply(unique(map$Chr), function(c){
  tmp = map$pos[which(map$Chr == c)]
  names(tmp) = map$Gene[which(map$Chr == c)]
  tmp})
names(use_map) = unique(map$Chr)

encoded <- getEncoding(matnut, terms = c("RIX","Diet","DietRIX"))


genoprobs_minimuga <- readRDS(file.path(dir, "mini/interp_rqtl_allChr_genoprob_11mar2020.rds"))
#genoprobs_minimuga = lapply(genoprobs_minimuga, function(x) x[-which(rownames(genoprobs_minimuga[[1]]) %in% problemPups),,])

matnut = matnut[match(rownames(genoprobs_minimuga[[1]]),matnut$Pup.ID),]
matnut = matnut[!is.na(matnut$Pup.ID), ]
unlist(lapply(genoprobs_minimuga, function(x) identical(rownames(x), matnut$Pup.ID)))

bin_pheno = matrix(ifelse(matnut$RIX == 6, 1, 0), ncol=1)
colnames(bin_pheno) = "sigRIX"
rownames(bin_pheno) = matnut$Pup.ID

out_bin <- scan1(genoprobs_minimuga, bin_pheno, model="binary")
pks = find_peaks(out_bin, use_map)

par(mar=c(5.1, 4.1, 1.1, 1.1))
ymx <- maxlod(out_bin) # overall maximum LOD score
plot(out_bin, map, lodcolumn="sigRIX", col="slateblue", ylim=c(0, ymx*1.02))
