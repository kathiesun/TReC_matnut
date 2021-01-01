library(tidyverse)

setwd("C:/Users/Kathie/TReC_matnut//src")

source("./matnut/prediction_functions.R")
source("./matnut/summary_functions.R")

dir <- "C:/Users/Kathie/Dropbox\ (ValdarLab)"
gecco <- read.csv(file.path(dir, "de_results/fullGeccoRnaDump.csv"))

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

##########################

regFiles = list.files(pattern = "summary_mnt_50k_29dec2020.rds", file.path(dir, "variant_data/regression_outputs/"),
                      full.names = T)

regs = lapply(regFiles[grep("chr17", regFiles)], readRDS)
mu_g = do.call("rbind", lapply(regs, function(x) x[[1]]$summary$mu_g))
mu_g$Gene.Name = unlist(strsplit(rownames(mu_g),"_"))[c(F,T)]
mu_g$Pup.ID  = as.integer(unlist(strsplit(rownames(mu_g),"_"))[c(T,F)])
mu_g = mu_g %>% left_join(colData, by="Pup.ID")

print(mu_g %>% filter(lower > 0.55 | upper < 0.45) %>%
    group_by(Gene.Name) %>% tally() %>%
    left_join(mu_g %>% group_by(Gene.Name) %>% tally(), by="Gene.Name") %>%
    mutate(imp_perc = n.x / n.y) %>%
    arrange(desc(imp_perc)) %>% 
    filter(n.x > 3, imp_perc > 0.1), n=30)
