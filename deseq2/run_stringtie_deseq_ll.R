options(stringsAsFactors = FALSE)
library(DESeq2)
library(tidyverse)
library(sva)

setwd("C:/Users/Kathie/rna_seq/deseq")
setwd("~/rna_seq/deseq2")
dir="/nas/depts/006/valdar-lab/users/sunk/"
source("deseq2_functions.R")

args <- commandArgs(trailingOnly = TRUE)
r = as.numeric(args[1])
method = as.numeric(args[2])
## method = 0 for StringTie, method = 1 for kmers

#######  READ IN SAMPLE DATA  ##########
samples <- read.csv(file.path(dir, "variant_data", "2015-10_expression_pups.csv"))
rownames(samples) = samples$Pup.ID
samples$Diet <- gsub(" $", "", samples$Diet)
samples$Diet <- factor(samples$Diet, levels=c("Standard", "Low Protein","Methyl Enriched","Vitamin D Deficient"))
samples$PO <- ifelse(factor(gsub("[0-9]","", samples$RIX)) == "a", 0.5, -0.5)
samples$PO_cat <- factor(gsub("[0-9]","",samples$RIX))

problemPups <- c(1404, 1716, 1371, 569, 1911, 1951, 1015)

######## READ IN STRINGTIE DATA #########
if(method == 0){
  countData <- as.matrix(read.csv(file.path(dir,"string_pipe_out/gene_count_matrix.csv"), row.names="X"))
  annot = read.table("/nas/depts/006/valdar-lab/users/sunk/mm10/ref_sequence/Mus_musculus.GRCm38.96.gtf", skip=5, fill=T)
  
  annot_genes = annot %>% filter(V15 == "gene_name") %>%
    dplyr::select(one_of(paste0("V",c(1,4,5,7,10,16)))) %>%
    distinct()
  
  colnames(annot_genes) = c("Chr","Start","End","Strand","Gene.ID","Gene.Name")
  if(length(which(duplicated(annot_genes$Gene.Name))) > 0){
    annot_genes = annot_genes[-which(duplicated(annot_genes$Gene.Name)),]
  }
  
  rownames(countData) = annot_genes$Gene.Name[match(rownames(countData), annot_genes$Gene.ID)]
  countData = countData[which(rownames(countData) %in% annot_genes$Gene.Name[which(annot_genes$Chr %in% c(1:19))]),]
}

########  READ IN KMER DATA  ###########
if(method == 1) {
  invar_info <- read.csv(file.path(dir,"variant_data/raw_data/GeneProbesInvariant.csv"))
  invar_info %>% filter(Chromosome %in% c(1:19)) -> invar_info
  rix_data = readRDS(file.path(dir, "variant_data/invar_counts_by_rix_5aug2019.rds"))
  rixes = factor(tolower(names(rix_data)), levels=paste0("rix",c(1:4,6:10)))
  rixes = rixes[order(rixes)]
  for(rix in rixes){
    user = rix_data[[toupper(rix)]]
    if(all(unlist(lapply(user, function(p) identical(row.names(user[[1]]), row.names(p)))))){
      comb = do.call("cbind", lapply(user, function(p) p["sum"]))
    }
    countRIX = data.frame(apply(comb, c(1,2), as.numeric))
    colnames(countRIX) <- names(user)
    countRIX$gene = unlist(user[[1]]$Gene.Name)
    if(rix == rixes[1]){
      countData = countRIX
    } else {
      countData = full_join(countData, countRIX, by="gene")
    }
  } 
  rownames(countData) = countData$gene
  countData = countData[which(countData$gene %in% invar_info$Gene.Name),]
  countData = countData[,-grep("gene", colnames(countData))]
  countData = countData[,-grep(paste(paste0("_",problemPups), collapse="|"), colnames(countData))]
  countData = as.matrix(countData)
}

######## RUN DESEQ2 #########

use_pups <- as.numeric(gsub("Pup.ID_","",colnames(countData)))
if(length(which(use_pups %in% problemPups)) > 0) use_pups <- use_pups[-which(use_pups %in% problemPups)]

colData <- samples %>% 
  dplyr::filter(Pup.ID %in% use_pups) %>%
  dplyr::select(one_of("Pup.ID","RRIX","Diet","PO","PO_cat","Breeding.Batch","Behavior.Batch")) %>%
  mutate(Breeding.Batch  = as.factor(Breeding.Batch), Behavior.Batch = as.factor(Behavior.Batch))
rownames(colData) = paste0("Pup.ID_", colData$Pup.ID)

countData <- countData[, rownames(colData)]
all(rownames(colData) == colnames(countData))

## new way
sva_nointer = run_DESeq2(colData, countData, sva=T, interact=F, rix=r)

if(method == 1){
	nosva_nointer = run_DESeq2(colData, countData, sva=F, interact=F, rix=r)
	nosva_inter = run_DESeq2(colData, countData, sva=F, interact=T, rix=r)
	sva_inter = run_DESeq2(colData, countData, sva=T, interact=T, rix=r)

	dds_lst = list(sva_nointer=sva_nointer[[1]], 
        	       nosva_nointer=nosva_nointer[[1]], 
		       nosva_inter=nosva_inter[[1]], 
		       sva_inter=sva_inter[[1]])
} else {
	dds_lst = sva_nointer[[1]]
}
meth_title = ifelse(method == 1, "_kmer_", "_string_")
saveRDS(dds_lst, file.path(dir, paste0("variant_data/regression_outputs/rix",r,meth_title,"deseq_sva_ddsLst_19feb2020.rds")))
