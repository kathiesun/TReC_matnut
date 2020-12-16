library(tidyverse)
#library(tximport)
#library(tximportData)
library(DESeq2)
library(GenomicFeatures)
library(Biostrings)
#library(refGenome)

setwd("C:/Users/Kathie/rna_seq/kmerSearch")
setwd("~/rna_seq/kmerSearch)/kmerSearch/")
dir <- "C:/Users/Kathie/Dropbox\ (ValdarLab)"

meta <- read.csv(file.path(dir, "matnut_main/metaseq.csv"))
ie1 <- read.table(file.path(dir, "imprinted_genes/imprinted_genes_mott.txt"), header = T, stringsAsFactors = F)
ie2 <- read.table(file.path(dir, "imprinted_genes/imprinted_genes_jirtle.txt"), header = T, sep=",", stringsAsFactors = F)
ie3 <- read.table(file.path(dir, "imprinted_genes/GSE27016_imprinted_genes.readme"), header = F, sep=";",stringsAsFactors = F)

firstup <- function(x) {
  x <- tolower(x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
ie1 <- unlist(strsplit(ie1$Gene,"[(|)]"))
ie2_all_alias <- unlist(lapply(unlist(strsplit(ie2$Aliases, ";")), firstup))
ie2 <- unlist(lapply(ie2$Gene[which(ie2$Status == "Imprinted")], firstup))
ie3$V1 <- gsub(",","", ie3$V1)
ie3 <- ie3$V1[-grep("PMID", ie3$V1)]
ie3 <- unlist(strsplit(ie3, "[_|/]"))
ie3 <- ie3[-grep("^v[0-9]$", ie3)]
compare_ie <- unique(unlist(union(union(ie2, ie3), ie1))) 
compare_ie <- compare_ie[order(compare_ie)]

invar_info <- read.csv(file.path(dir,"variant_data/kmers_based_on_CC/kmers_to_run/GeneProbesInvariant.csv"))
invar_info = invar_info %>% dplyr::select(c("Gene.Name","Gene.ID")) %>% distinct()

dds_sea <- readRDS(file.path(dir, "de_results/dds_lst-5aug2019.rds"))
dds_cs <- readRDS(file.path(dir, "de_results/dds_lst_cs-5aug2019.rds"))
dds_mm10 <- readRDS(file.path(dir, "de_results/dds_lst_mm10transcript-15sep2019.rds"))
names(dds_cs) <- tolower(names(dds_cs))


dds_cs

#for(i in 1:length(dds_cs)){
#  rownames(dds_cs[[i]]) <- invar_info$Gene.ID[match(rownames(dds_cs[[i]]), invar_info$Gene.Name)]
#}

###### sample data ######
samples <- read.csv(file.path(dir, "variant_data", "2015-10_expression_pups.csv"))
rownames(samples) = samples$Pup.ID
samples$Diet <- gsub(" $", "", samples$Diet)
samples$Diet <- factor(samples$Diet, levels=c("Standard", "Low Protein","Methyl Enriched","Vitamin D Deficient"))
samples$PO <- ifelse(factor(gsub("[0-9]","", samples$RIX)) == "a", 0.5, -0.5)
samples$PO_cat <- factor(gsub("[0-9]","",samples$RIX))

pups = unlist(lapply(dds_cs, colnames))
pups = gsub("[A-Z|_|.]","",toupper(pups))

samples_use <- samples %>% 
  filter(Pup.ID %in% as.numeric(pups)) %>%
  dplyr::select(one_of("Pup.ID","RRIX","Diet","PO","PO_cat","Breeding.Batch","Behavior.Batch")) %>%
  mutate(Breeding.Batch  = as.factor(Breeding.Batch), Behavior.Batch = as.factor(Behavior.Batch))
#coldata <- data.frame(files=use_files, names=names(use_files), samples_use)

#samples_use %>% group_by(RRIX, PO, Diet) %>% tally() %>% filter(n < 3 )
#samples_use %>% group_by(RRIX, PO) %>% tally() 

####### ref fa and gecco ###########
ref_fa <- readDNAStringSet(file.path(dir, "matnut_main/Mus_musculus.GRCm38.cdna.all.fa"))

seq_name = names(ref_fa)                                                                                                             
seq_det = lapply(seq_name, function(x) {                                                                                            
  unlist(strsplit(x, " "))[c(1,3,4,7)]                                                                                                       
  #keep = grep("chr|gene:|gene_symbol|ENS", tmp)                                                                                        
  #tmp[keep]                                                                                                             
})                                                                                                                                
seq_det <- do.call("rbind", seq_det)   
tx2gene <- data.frame(TXNAME=seq_det[,1], GENEID=gsub("gene:","",seq_det[,3]))    
transcripts <- do.call("rbind", strsplit(seq_det[,1], "[.]"))                                                                                                 
chr <- do.call("rbind", strsplit(seq_det[,2], "[:]"))                                                                                            
gene <- do.call("rbind", strsplit(seq_det[,3], "[:|.]"))                                                                  
gene_sym <- do.call("rbind", strsplit(seq_det[,4], "[:]"))                                                                  
seq_det <- cbind(tx2gene, transcripts, chr, gene, gene_sym)                                                                                                      
colnames(seq_det) = c("TXNAME","GENEID", "transcript","tr_num", "del","build","chr","start","end",
                      "dir","del2","ensemblID","gene_num","del3","gene")                                                                                                                             
seq_det <- seq_det[,-grep("del", colnames(seq_det))]
seq_det <- data.frame(seq_det)        
seq_det %>% mutate(start = as.numeric(paste(start)),                                                                         
                   end = as.numeric(paste(end)),
                   sequence = paste(ref_fa)) -> seq_det_full

gecco <- read.csv(file.path(dir, "de_results/fullGeccoRnaDump.csv"))
gecco <- gecco %>% left_join(seq_det, by=c("ensembl_id"="ensemblID")) 
gecco_uniq <- gecco %>% dplyr::select(ensembl_id, sample_id, TReC, gene) %>%
  distinct()
gecco %>% 
  filter(!is.na(gene)) %>%
  group_by(gene) %>% 
  summarize(counts_gecco=mean(TReC)) %>%
  #rename("ensembl_id"="ensemblID") %>%
  dplyr::select(gene, counts_gecco) -> gecco_sum

#####################  Analyze any set of dds_lst result  ###########################

compare_diets_PO <- function(dds_lst, ref_diet="Standard"){
  res_lst_PO <- lapply(dds_lst, function(x) results(x, name="PO"))
  diets <- gsub(" ",".",as.character(unique(dds_lst[[1]]$Diet)))
  other_diets = diets[-which(diets == ref_diet)]
  res_lst_diet = list()
  res_lst_diet <- lapply(dds_lst, function(x){
    out = list()
    for(d in other_diets){
      if(length(grep(ref_diet, resultsNames(x)) ) > 0 & length(grep(d, resultsNames(x)) ) > 0){
        out[[d]] = results(x, contrast=c("Diet", ref_diet, d))
      }
    }
    return(out)
  })
  return(list(res_lst_PO = res_lst_PO, res_lst_diet = res_lst_diet))
}

dds_cs_compare = compare_diets_PO(dds_lst, ref_diet = "Vitamin.D.Deficient")
dds_mm10_compare = compare_diets_PO(dds_mm10, ref_diet = "Vitamin.D.Deficient")
dds_string_compare = compare_diets_PO(dds_string, ref_diet = "Vitamin.D.Deficient")

pval <- 0.05


cs_PO_0.05 <- lapply(dds_cs_compare$res_lst_PO, function(x){
  x$gene = rownames(x)
  x%>% as.data.frame() %>% filter(padj < pval) %>% arrange(padj)
})
  
cs_genes <- table(unlist(lapply(cs_PO_0.05, function(x) x$gene)))
cs_genes[which(cs_genes > 1)]

mm10_PO_0.05 <- lapply(dds_mm10_compare$res_lst_PO, function(x){
  x$gene = rownames(x)
  x%>% as.data.frame() %>% filter(padj < pval) %>% arrange(padj)
})

mm10_genes <- table(unlist(lapply(mm10_PO_0.05, function(x) x$gene)))
mm10_genes[which(mm10_genes > 1)]

all_cs_genes = table(unlist(lapply(cs_res_PO, function(x) rownames(x))))
cs_genes_consistent = names(all_cs_genes)[which(all_cs_genes == max(all_cs_genes))]

genes_consistent = intersect(cs_genes_consistent, rownames(res_PO$rix1))

str_cs_compare = list()
for(r in names(res_PO)){
  tmp_string = res_PO[[r]][match(genes_consistent, rownames(res_PO[[r]])),]
  tmp_cs = cs_res_PO[[r]][match(genes_consistent, rownames(cs_res_PO[[r]])),]
  colnames(tmp_string) = paste0(colnames(tmp_string),"_str")
  colnames(tmp_cs) = paste0(colnames(tmp_cs),"_cs")
  tmp_string$gene = tmp_cs$gene = genes_consistent

  str_cs_compare[[r]] = left_join(as.data.frame(tmp_string), as.data.frame(tmp_cs), by="gene")
}

sapply(1:9, function(r) 
  cor(x=str_cs_compare[[r]]$pvalue_str, y=str_cs_compare[[r]]$pvalue_cs, use="pairwise.complete.obs"))
for(r in names(str_cs_compare)){
  plot(x=str_cs_compare[[r]]$pvalue_str, y=str_cs_compare[[r]]$pvalue_cs,
       main = paste("P-value comparison in",r), xlab = "StringTie", ylab = "K-mers")
}


string_PO_0.05 <- lapply(res_PO, function(x){
  x$gene = rownames(x)
  x%>% as.data.frame() %>% filter(padj < pval) %>% arrange(padj)
})

string_genes <- table(unlist(lapply(string_PO_0.05, function(x) x$gene)))
plot_genes = string_genes[which(string_genes > 1)]
for(x in 1:length(string_PO_0.05)){
  if(nrow(string_PO_0.05[[x]]) > 0){
    string_PO_0.05[[x]]$rix = names(string_PO_0.05)[x]
  }
}
string_PO_0.05 = do.call("rbind", string_PO_0.05)
string_PO_0.05 = string_PO_0.05 %>% #filter(gene %in% names(plot_genes)) %>%
  group_by(gene) %>% arrange(-baseMean)
plot_genes = string_PO_0.05$gene[-grep("Gm|mt|rik", string_PO_0.05$gene)]
pdf("C:/Users/Kathie/Dropbox (ValdarLab)/figures_and_updates/images_meetings/top_string_PO_counts.pdf")
par(mfrow=c(2,1))
for(i in 1:length(plot_genes)){
  gene = plot_genes[i]
  gene="Fam205a2"
  rix= tolower(string_PO_0.05$rix[i])    #c(4,4,6,8,8,9)
  rix="rix6"
  pd <- plotCounts(dds_string[[rix]], gene=gene, intgroup = "PO",returnData = T)
  plot(x=as.factor(pd$PO), y=pd$count,
       main=paste("Normalized counts in",gene, rix))       
  points(x=jitter(as.numeric(as.factor(pd$PO))), y=pd$count)  
  #anova(lm(data=pd, count ~ PO))
}
dev.off()
### 99 sig genes

dds_cs_compare$res_lst_diet$rix6$Standard[order(dds_cs_compare$res_lst_diet$rix6$Standard$padj),]
dds_mm10_compare$res_lst_diet$rix6$Standard[order(dds_mm10_compare$res_lst_diet$rix6$Standard$padj),]
res_diet$rix6$Standard[order(res_diet$rix6$Standard$padj),]

dds_cs_compare$res_lst_diet$rix6$Standard$Gene = rownames(dds_cs_compare$res_lst_diet$rix6$Standard) 
sig_VDD_STD_rix6_cs = as.data.frame(dds_cs_compare$res_lst_diet$rix6$Standard) %>% filter(padj < pval)

dds_mm10_compare$res_lst_diet$rix6$Standard$Gene = rownames(dds_mm10_compare$res_lst_diet$rix6$Standard) 
sig_VDD_STD_rix6_mm10 = as.data.frame(dds_mm10_compare$res_lst_diet$rix6$Standard) %>% filter(padj < pval)

res_diet$rix6$Standard$Gene = rownames(res_diet$rix6$Standard) 
sig_VDD_STD_rix6_string = as.data.frame(res_diet$rix6$Standard) %>% filter(padj < pval) %>%
  arrange(padj)


length(which(sig_VDD_STD_rix6_mm10$Gene %in% sig_VDD_STD_rix6_string$Gene))
length(unique(sig_VDD_STD_rix6_mm10$Gene))

length(which(sig_VDD_STD_rix6_cs$Gene %in% sig_VDD_STD_rix6_string$Gene))
length(unique(sig_VDD_STD_rix6_string$Gene))
length(unique(sig_VDD_STD_rix6_cs$Gene))



#############   START ANALYZING SEA   ###############
pval <- 0.05
#name = "PO"

### DIET ###
res_sea_VDDSTD <- lapply(dds_sea, function(x){
  if(length(grep("Vitamin.D.Deficient", resultsNames(x))) > 0){
    results(x, contrast=c("Diet", "Vitamin.D.Deficient", "Standard"))
  }
})
res_sea_VDDSTD_pval<- lapply(res_sea_VDDSTD, function(x) x[which(x$padj < pval), ])
sea_rix6_VDDSTD_sig <- res_sea_VDDSTD_pval$rix6[order(res_sea_VDDSTD_pval$rix6$pvalue),]
unassgn <- lapply(res_sea_VDDSTD_pval, function(x) counts(x)[grep("gene", rownames(counts(x))),])
resLFC <- lapply(dds_sea, function(x) lfcShrink(x, coef="PO", type="apeglm"))


#### PO ####
res_sea_PO <- lapply(dds_sea, function(x) results(x, name="PO"))  
res_sea_PO_pval <- lapply(res_sea_PO, function(x){
  tmp <- x[which(x$padj < 0.1), ]
  tmp[order(tmp$pvalue),] })
res_sea_PO_pval_all <- sapply(1:length(res_sea_PO_pval), function(x){
  tmp <- data.frame(res_sea_PO_pval[[x]]) 
  tmp %>% mutate(gene=rownames(res_sea_PO_pval[[x]]), rix=names(res_sea_PO_pval)[[x]]) -> tmp
  tmp
}, simplify=F)
res_sea_PO_pval_all <- do.call("rbind", res_sea_PO_pval_all)

res_sea_PO_pval[[which(lapply(res_sea_PO_pval, nrow) == 0)]] <- NULL

res_sea_PO_match <- data.frame(gene = rownames(res_sea_PO_pval[[1]]),
                      rix = names(res_sea_PO_pval)[[1]], 
                      pval = res_sea_PO_pval[[1]]$padj, stringsAsFactors = F)
        
for(i in 2:length(res_sea_PO_pval)){
  res_sea_PO_match <- rbind(res_sea_PO_match, data.frame(gene = rownames(res_sea_PO_pval[[i]]),
                                       rix = names(res_sea_PO_pval)[[i]], 
                                       pval = res_sea_PO_pval[[i]]$padj, stringsAsFactors = F))
  #yes <- c(yes, ifelse(length(grep(paste(cs_sea_PO, collapse = "|"), rownames(resPO_set[[i]]))>1), 1, 0))
}

###8,1,4,9

allPO_sea <- data.frame(table(res_sea_PO_match), stringsAsFactors = F) %>% 
  mutate(gene = as.character(gene)) %>%
  filter(Freq > 0) %>% arrange(-Freq, pval)

concord <- unique(compare_ie[which(compare_ie %in% res_sea_PO_pval_all$gene)])
plot_concord <- data.frame(res_sea_PO_pval_all, stringsAsFactors = F) %>% 
  filter(gene %in% concord)

#unassgn_ref <- lapply(dds_ref, function(x) counts(x)[grep("gene", rownames(counts(x))),])
#resLFC <- lapply(dds_sea, function(x) lfcShrink(x, coef="PO", type="apeglm"))

#dds <- dds_sea$rix6
#res <- res_sea_VDDSTD$rix6

gene=rownames(res_sea_PO_pval[[paste0("rix",rix)]])[1]
dlist <- list()
for(i in 1:nrow(plot_concord)){
  dlist[[i]] <- plotCounts(dds_sea[[plot_concord$rix[i]]], gene=plot_concord$gene[i], 
                           intgroup="PO", returnData = T)
  dlist[[i]]$gene = plot_concord$gene[i]
  dlist[[i]]$rix = gsub("rix", "", plot_concord$rix[i])
}
dlist <- do.call("rbind", dlist)
dlist$PO <- as.factor(dlist$PO)
dlist$gene <- factor(dlist$gene, levels=unique(dlist$gene))
d <- ggplot(dlist, aes(x=PO, y=count, col=rix)) +    
  geom_point(position=position_jitter(w=0.05,h=0)) + 
  scale_y_log10() + #breaks=c(25,100,400)
  #ggtitle(paste(gene, "in RIX", gsub("rix","",rix))) + 
  theme_minimal() + 
  facet_wrap( ~ gene, scales = "free_y", ncol=2)

## which.min(res$padj)

###### plotting for sea ######

pdf(file.path(dir,"sea60_PO_0.1_28may2019.pdf"))
count=1
for(i in 1:10){
  dlist <- list()
  for(j in 1:6){
    if(count <= nrow(allPO_sea)){
      rix = as.character(allPO_sea$rix[count])
      gene = as.character(allPO_sea$gene[count])
      #rownames(resPO[[rix]][order(resPO[[rix]]$padj)[rank],])
      #rix=paste0("rix",i)
      
      dlist[[j]] <- plotCounts(dds_sea[[rix]], gene=gene, intgroup="PO", returnData = T)
      dlist[[j]]$gene = gene
      dlist[[j]]$rix = rix
      count=count+1
    }
  }
  
  
  dlist <- do.call("rbind", dlist)
  dlist$PO <- as.factor(dlist$PO)
  dlist$gene <- factor(dlist$gene, levels=unique(dlist$gene))
  d <- ggplot(dlist, aes(x=PO, y=count, col=rix)) +    
    geom_point(position=position_jitter(w=0.05,h=0)) + 
    scale_y_log10() + #breaks=c(25,100,400)
    #ggtitle(paste(gene, "in RIX", gsub("rix","",rix))) + 
    theme_minimal() + 
    facet_wrap( ~ gene, scales = "free_y", ncol=2)
  print(d)
}

dev.off()

#################   MORE SEA PLOTS   ######################
i=4
rank=1
rix=paste0("rix",i)
gene=rownames(resPO[[rix]][order(resPO[[rix]]$padj)[rank],])
gene="Igf2"
d <- plotCounts(dds_sea[[rix]], 
                gene=gene, intgroup="PO", returnData = T)
d$PO <- as.factor(d$PO)
title = rownames(resPO[[rix]][order(resPO[[rix]]$padj)[rank],])
ggplot(d, aes(x=PO, y=count)) +    
  geom_point(position=position_jitter(w=0.05,h=0)) + 
  scale_y_log10() + #breaks=c(25,100,400)
  #ggtitle(paste(title, "in RIX", i)) + 
  theme_minimal()

d <- plotCounts(dds_ref[[rix]], gene=rownames(res[[rix]][order(res[[rix]]$padj)[rank],]), intgroup="PO", returnData = T)
d$PO <- as.factor(d$PO)
title = rownames(res[[rix]][order(res[[rix]]$padj)[rank],])
ggplot(d, aes(x=PO, y=count)) +    
  geom_point(position=position_jitter(w=0.05,h=0)) + 
  scale_y_log10() + #breaks=c(25,100,400)
  ggtitle(paste(title, "in RIX", i)) + 
  theme_minimal()



#####################    START ANALYZING CS    ############################
pval = 0.05
#res <- results(dds, c("PO_cat","b","a"))
res_cs_PO <- lapply(dds_cs, function(x) results(x, name="PO"))#
lapply(dds_cs, resultsNames)
res_cs_VDDSTD <- lapply(dds_cs, function(x){
    if(length(grep("Standard", resultsNames(x)) ) > 0 & length(grep("Vitamin.D.Deficient", resultsNames(x)) ) > 0){
      results(x, contrast=c("Diet", "Vitamin.D.Deficient", "Standard"))
    }
  })
res_cs_VDDLP <- lapply(dds_cs, function(x){
  if(length(grep("Low.Protein", resultsNames(x)) ) > 0 & length(grep("Vitamin.D.Deficient", resultsNames(x)) ) > 0){
    results(x, contrast=c("Diet", "Vitamin.D.Deficient", "Low.Protein"))
  }
})
res_cs_VDDME_rix6 <- results(dds_cs$rix6, contrast=c("Diet", "Vitamin.D.Deficient", "Methyl.Enriched"))
res_cs_VDDLP_rix6 <- results(dds_cs$rix6, contrast=c("Diet", "Vitamin.D.Deficient", "Low.Protein"))
res_cs_VDDSTD_rix6 <- results(dds_cs$rix6, contrast=c("Diet", "Vitamin.D.Deficient", "Standard"))


me <- data.frame(res_cs_VDDME_rix6) %>% mutate(gene = rownames(res_cs_VDDME_rix6), compare = "VDDME") %>% filter(padj < 0.05) %>% arrange(padj)
std <- data.frame(res_cs_VDDSTD_rix6) %>% mutate(gene = rownames(res_cs_VDDSTD_rix6), compare = "VDDSTD") %>% filter(padj < 0.05) %>% arrange(padj)
lp <- data.frame(res_cs_VDDLP_rix6) %>% mutate(gene = rownames(res_cs_VDDLP_rix6), compare = "VDDLP") %>% filter(padj < 0.05) %>% arrange(padj)  

all_sig_vdd <- rbind(me, std, lp)


lapply(res_cs_PO, function(x) x[which(rownames(x) == "Cyp27b1"),])

  
  
results(dds_cs, contrast=c("Diet", "Vitamin.D.Deficient", "Standard"))
#res
res_cs_PO_pval <- lapply(res_cs_PO, function(x){
  tmp <- x[which(x$padj < pval), ]
  tmp[order(tmp$pvalue),] })
res_cs_PO_pval_all <- sapply(1:length(res_cs_PO_pval), function(x){
  tmp <- data.frame(res_cs_PO_pval[[x]]) 
  tmp %>% mutate(gene=rownames(res_cs_PO_pval[[x]]), rix=names(res_cs_PO_pval)[[x]]) -> tmp
}, simplify=F)
res_cs_PO_pval_all <- do.call("rbind", res_cs_PO_pval_all)

res_cs_PO_match <- data.frame(gene = rownames(res_cs_PO_pval[[1]]),
                      rix = names(res_cs_PO_pval)[[1]],
                      pval = res_cs_PO_pval[[1]]$padj, stringsAsFactors = F)

for(i in 2:length(res_cs_PO_pval)){
  if(length(rownames(res_cs_PO_pval[[i]])) > 0){
    res_cs_PO_match <- rbind(res_cs_PO_match, data.frame(gene = rownames(res_cs_PO_pval[[i]]),
                                               rix = names(res_cs_PO_pval)[[i]],
                                               pval = res_cs_PO_pval[[i]]$padj, stringsAsFactors = F))
    #yes <- c(yes, ifelse(length(grep(paste(cs_sea_PO, collapse = "|"), rownames(resPO_set[[i]]))>1), 1, 0))
  }
}

###8,1,4,9

allPO_cs <- data.frame(table(res_cs_PO_match))%>% filter(Freq > 0) %>% arrange(-Freq)
cs_sea_PO <- allPO_cs$gene[which(allPO_cs$gene %in% allPO_sea$gene)]
allPO_cs[which(allPO_cs$gene %in% cs_sea_PO),]
allPO_sea[which(allPO_sea$gene %in% cs_sea_PO),]

concord <- unique(compare_ie[which(compare_ie %in% res_cs_PO_pval_all$gene)])
plot_concord <- res_cs_PO_pval_all %>% filter(gene %in% concord) %>%
  mutate(gene = as.character(gene), rix = as.character(rix))
#length(which(igenes %in% rownames(res_cs_PO$rix1)))

dlist <- list()
for(i in 1:nrow(plot_concord)){
  dlist[[i]] <- plotCounts(dds_cs[[plot_concord$rix[i]]], gene=plot_concord$gene[i], 
                           intgroup="PO", returnData = T)
  dlist[[i]]$gene = plot_concord$gene[i]
  dlist[[i]]$rix = gsub("rix", "", plot_concord$rix[i])
}
dlist <- do.call("rbind", dlist)
dlist$PO <- as.factor(dlist$PO)
dlist$gene <- factor(dlist$gene, levels=unique(dlist$gene))
d <- ggplot(dlist, aes(x=PO, y=count, col=rix)) +    
  geom_point(position=position_jitter(w=0.05,h=0)) + 
  scale_y_log10() + #breaks=c(25,100,400)
  #ggtitle(paste(gene, "in RIX", gsub("rix","",rix))) + 
  theme_minimal() + 
  facet_wrap( ~ gene, scales = "free_y", ncol=2)



gene = "Cd81"  ## 8
gene = "Usp29" ## 8
gene = "Lrp1"    ## 4
gene = "Lrrtm1"  ## 4

rix = 1
d <- plotCounts(dds_lst[[paste0("rix",rix)]], gene=gene, 
                intgroup="PO", returnData = T)
d$PO <- as.factor(d$PO)
ggplot(d, aes(x=PO, y=count)) +    
  geom_point(position=position_jitter(w=0.05,h=0)) + 
  scale_y_log10() + #breaks=c(25,100,400)
  ggtitle(paste(gene, "in RIX", rix)) + 
  theme_minimal()


pval <- 0.1
table(unlist(lapply(res_cs_PO, function(x) rownames(x[which(x$padj < pval),]))))

##### plotting for cs ######
#allPO_cs
pdf(file.path(dir,"cs_PO_0.1_28may2019.pdf"))
count=1
for(i in 1:ceiling(nrow(allPO_cs)/6)){
  dlist <- list()
  for(j in 1:6){
    if(count <= nrow(allPO_cs)){
      rix = as.character(allPO_cs$rix[count])
      gene = as.character(allPO_cs$gene[count])
      
      #rownames(resPO[[rix]][order(resPO[[rix]]$padj)[rank],])
      #rix=paste0("rix",i)
      
      dlist[[j]] <- plotCounts(dds_lst[[rix]], gene=gene, intgroup="PO", returnData = T)
      dlist[[j]]$gene = gene
      dlist[[j]]$rix = rix
      count=count+1
    }
  }
  
  
  dlist <- do.call("rbind", dlist)
  dlist$PO <- as.factor(dlist$PO)
  dlist$gene <- factor(dlist$gene, levels=unique(dlist$gene))
  d <- ggplot(dlist, aes(x=PO, y=count, col=rix)) +    
    geom_point(position=position_jitter(w=0.05,h=0)) + 
    scale_y_log10() + #breaks=c(25,100,400)
    #ggtitle(paste(gene, "in RIX", gsub("rix","",rix))) + 
    theme_minimal() + 
    facet_wrap( ~ gene, scales = "free_y", ncol=2)
  print(d)
}

dev.off()


#### DIET ####
res_cs_VDDSTD <- lapply(dds_cs, function(x){
  if(length(grep("Vitamin.D.Deficient", resultsNames(x))) > 0){
    results(x, contrast=c("Diet", "Vitamin.D.Deficient", "Standard"))
  }
})
res_cs_VDDSTD_pval<- lapply(res_cs_VDDSTD, function(x) x[which(x$padj < pval), ])
cs_rix6_VDDSTD_sig <- res_cs_VDDSTD_pval$rix6[order(res_cs_VDDSTD_pval$rix6$pvalue),]


#############  COMPARE RIX 6 VDD vs STD  ##############


use_cs <- as_tibble(counts(dds_lst$rix6, normalized=T)) %>% mutate(gene = rownames(counts(dds_cs$rix6))) %>%
  gather(pup, counts, -gene) %>%
  mutate(pup=gsub("^X", "", pup), method = "cs")
use_mm10 <- as_tibble(counts(dds_mm10$rix6, normalized=T)) %>% mutate(gene = rownames(counts(dds_mm10$rix6))) %>%
  gather(pup, counts, -gene) %>%
  mutate(pup=paste0("Pup.ID_", pup), method = "mm10")
use_string <- as_tibble(counts(dds_string$rix6, normalized=T)) %>% mutate(gene = rownames(counts(dds_string$rix6)))  %>% gather(pup, counts, -gene) %>%
  mutate(method = "string")
#use_ref <- as_tibble(counts(dds_ref$rix6)) %>% mutate(gene = rownames(counts(dds_ref$rix6)))  %>% gather(pup, counts, -gene) %>%
#  mutate(method = "ref")
use_gecco <- gecco_uniq %>% rename("sample_id" = "pup", "TReC" = "counts") %>% dplyr::select(-ensembl_id) %>%
  mutate(method = "gecco") %>% dplyr::select(gene, pup, counts, method) %>%
  filter(!is.na(gene))
compare <- rbind(use_cs, use_mm10, use_string, use_gecco) #use_ref, 


factors <- data.frame(method=c("cs", "mm10", "string", "gecco"), #, "ref"), 
                 means=c(mean(use_cs$counts),mean(use_mm10$counts),mean(use_string$counts),mean(use_gecco$counts)), #mean(use_ref$counts)),
                 stringsAsFactors = F)
factors$fac <- factors$means/min(factors$means)

n=length(use_gecco$gene)
use1 <- intersect(use_gecco$gene[order(-use_gecco$counts)[1:n]], use_cs$gene[order(-use_cs$counts)])
use2 <- intersect(use_string$gene[order(-use_string$counts)], use1)
#use_ref$gene[order(-use_ref$counts)[1:n]]

compare %>% filter(gene %in% use2) -> plotMat

#d %>% group_by(method, gene) %>% summarise(mean=mean(counts)) %>% filter(method != "ref") -> divide
#for(g in unique(d$gene)){
#  div <- divide[which(divide$gene == g),]
#  d$counts[which(d$gene == g)] <- unlist(apply(d[which(d$gene == g),], 1, function(x) as.numeric(x["counts"])/(div$mean[which(div$method == x["method"])]/min(div$mean))))
#}
plotMat$fcounts <- rep(NA)
for(m in unique(plotMat$method)){
  plotMat$fcounts[which(plotMat$method == m)] <- plotMat$counts[which(plotMat$method == m)]/factors$fac[which(as.character(factors$method) == m)]
}


#plotMat %>% left_join(means, by=c("gene", "method")) -> plotMat
plotMat <- data.frame(plotMat)
d <- ggplot(plotMat, aes(x=gene, y=fcounts+0.5, col=method)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  #geom_line(data=means, aes(x=gene, y=mean, col=method)) +
  theme_minimal() +
  scale_y_log10() + 
  stat_summary(aes(y=fcounts+0.5,group=method), geom="line", fun.y = 'mean') +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))# + 

d

plotMat %>% group_by(gene, method) %>% summarize(mean=mean(fcounts)) -> means

  ggplot(means, aes(x=gene, y=mean, col=method)) + 
  geom_point() + 
  theme_minimal() +  
  stat_summary(aes(y=mean,group=method), geom="line") +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))# + 

#means %>% filter(gene != "Atp1a3") -> means
corr <- means %>% spread(method, mean)
cor.test(corr$cs, corr$string)
cor.test(corr$cs, corr$gecco)
cor.test(corr$gecco, corr$string)
ggplot(corr, aes(x=gecco, y=string)) + 
  geom_point()
cor(corr[,-1])

#############  COMPARE RIX 6 VDD vs STD  ##############
#############  @ GENE LEVEL  ##############

sea_rix6_VDDSTD_sig
cs_rix6_VDDSTD_sig
#ref_rix6_VDDSTD

write.csv(as.data.frame(sea_rix6_VDDSTD_sig), 
          file=file.path(dir, "sea_VDDSTD_0.001_27may2019.csv"))
write.csv(as.data.frame(cs_rix6_VDDSTD_sig), 
          file=file.path(dir, "cs_VDDSTD_0.001_27may2019.csv"))



if(method == "sea"){
  df <- sea_rix6_VDDSTD_sig
  deUse <- dds_sea$rix6
} else if (method == "cs"){
  df <- cs_rix6_VDDSTD_sig
  deUse <- dds_cs$rix6
} else {
  df <- ref_rix6_VDDSTD
  deUse <- dds
}
genes <- rownames(df)

genes <- intersect(genes, rownames(sea_rix6_VDDSTD_sig)[which(rownames(sea_rix6_VDDSTD_sig) %in% rownames(cs_rix6_VDDSTD_sig))])
genes <- genes[unlist(lapply(genes, function(x) gecco_sum$counts_gecco[gecco_sum$gene == x] > 1000))]

################## starting here 27 sept 2019 ####################
method = "sea"

#if(method == "sea"){
#  genes <- intersect(rownames(sea_rix6_VDDSTD_sig), rownames(cs_rix6_VDDSTD_sig))
#} else {
#  genes <- intersect(rownames(cs_rix6_VDDST_sig), rownames(sea_rix6_VDDSTD_sig))
#}

genes = sig_VDD_STD_rix6_string$Gene


### 'factors' made on 435
plotGenes <- function(ranks){
  plot_dat <- list()
  ind = 1
  for(i in ranks){
    plot_dat[[ind]] <- 
      rbind(data.frame(pup = names(counts(dds_string$rix6)[which(rownames(counts(dds_string$rix6)) == genes[i]),]),
                       counts = counts(dds_string$rix6, normalized=T)[which(rownames(counts(dds_string$rix6)) == genes[i]),], 
                       method = "string", stringsAsFactors = F),
            data.frame(pup = paste0("Pup.ID_", names(counts(dds_mm10$rix6)[which(rownames(counts(dds_mm10$rix6)) == genes[i]),])),
                       counts = counts(dds_mm10$rix6)[which(rownames(counts(dds_mm10$rix6)) == genes[i]),], 
                       method = "mm10", stringsAsFactors = F))
    if(genes[i] %in% gecco_uniq$gene){
      plot_dat[[ind]] <- rbind(plot_dat[[ind]], 
                               data.frame(pup = seq(1:length(gecco_uniq$TReC[which(gecco_uniq$gene == genes[i])])),
                               counts = gecco_uniq$TReC[which(gecco_uniq$gene == genes[i])],
                               method = "gecco", stringsAsFactors = F))
    }
            
    if(genes[i] %in% rownames(counts(dds_cs$rix6))){
      plot_dat[[ind]] <- rbind(plot_dat[[ind]], 
                               data.frame(pup = names(counts(dds_cs$rix6)[which(rownames(counts(dds_cs$rix6)) == genes[i]),]), 
                                          counts = counts(dds_cs$rix6, normalized=T)[which(rownames(counts(dds_cs$rix6)) == genes[i]),], 
                                          method = "cs", stringsAsFactors = F))
    }
    
    plot_dat[[ind]]$fcounts <- apply(plot_dat[[ind]] , 1, function(x) as.numeric(x["counts"])/as.numeric(factors$fac[which(factors$method == x["method"])]))
    plot_dat[[ind]]$PO <- apply(plot_dat[[ind]], 1, function(x) if(!x["method"] == "gecco"){
      colData(dds_cs$rix6)$PO[which(colData(dds_cs$rix6)$Pup.ID == as.character(x["pup"]))]
    } else {
      "gecco"
    })
    
    plot_dat[[ind]]$Diet <- apply(plot_dat[[ind]] , 1, function(x) if(!x["method"] == "gecco"){
      as.character(colData(dds_cs$rix6)$Diet[which(colData(dds_cs$rix6)$Pup.ID == x["pup"])])
    } else {
      "gecco"
    })
    plot_dat[[ind]]$gene = genes[i]
    ind=ind+1
  }
  
  plot_dat <- do.call("rbind", plot_dat)
  
  #d <- plotCounts(deUse, gene=genes[rank], intgroup="Diet", returnData = T)
  #gec <- gecco_sum$counts_gecco[which(gecco_sum$gene == genes[rank])]
  
  #$PO <- as.factor(d$PO)
  title = paste("Counts of genes in RIX 6")
  d <- ggplot(plot_dat, aes(x=method, y=fcounts+0.5, col=method)) +    
    geom_point(position=position_jitter(w=0.05,h=0)) + 
    scale_y_continuous(trans='log10') + #breaks=c(25,100,400)
    #ggtitle(title) + 
    facet_wrap( ~ gene, ncol=2, scales="free_y") + 
    theme_minimal() + 
    theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0.5))
  return(d)
  
}

#plotCounts(dds, gene=genes[rank], intgroup="Diet")
pdf(file.path(dir,"method_compare_VDDSTD_0.01_siggenes_28sep2019.pdf"))
for(j in 1:5){#ceiling(length(genes)/4)){
  ranks = ((j-1)*4 + 1): ((j-1)*4 + 4)
  print <- plotGenes(rank = ranks)
  print(print)
}
dev.off()

int1 <- intersect(rownames(counts(dds_string$rix6)),rownames(counts(dds_mm10$rix6)))
int2 <- intersect(int1, rownames(counts(dds_cs$rix6)))
int3 <- intersect(int2, gecco_uniq$gene)

corr_dat <- data.frame(string = counts(dds_string$rix6, normalized=T)[match(int3, rownames(counts(dds_string$rix6)))],
                       gecco = gecco_uniq$TReC[match(int3, gecco_uniq$gene)],
                       cs = counts(dds_cs$rix6, normalized=T)[match(int3, rownames(counts(dds_cs$rix6)))],
                       mm10 = counts(dds_mm10$rix6, normalized=T)[match(int3, rownames(counts(dds_mm10$rix6)))])
rownames(corr_dat) = int3
cor(corr_dat)
require(lattice)
levelplot(cor(corr_dat))


 ###########################

sea_rix6 <- as.character(rownames(counts(dds_sea$rix6)[which(unlist(apply(counts(dds_sea$rix6), 1, sum) > 10)),]))
    ## 10418 rows > 10 counts
    ### 13058 with new data
cs_rix6 <- as.character(rownames(counts(dds_cs$rix6)[which(unlist(apply(counts(dds_cs$rix6), 1, sum) > 10)),]))
    ## 14931 rows > 10 counts
ref_rix6 <- as.character(rownames(counts(dds_ref$rix6)[which(unlist(apply(counts(dds_ref$rix6), 1, sum) > 10)),]))
    ## 18047 rows
gecco_names <- as.character(unique(gecco_sum$gene)[which(gecco_sum$counts_gecco > 10)])
    ## 23075; 15361 > 10 counts

sea_cs <- setdiff(intersect(sea_rix6, cs_rix6), gecco_names)
    ## 351; 303 > 10 counts
sea_gecco <- setdiff(intersect(sea_rix6, gecco_names), cs_rix6)
    ## 2263; 1566 > 10 counts
cs_gecco <- setdiff(intersect(cs_rix6, gecco_names), sea_rix6)
    ## 5261; 5195 > 10 counts
sea_only <- setdiff(sea_rix6, union(gecco_names, cs_rix6))
    ## 1873; 1734 > 10 counts
gecco_only <- setdiff(gecco_names, union(sea_rix6, cs_rix6))
    ## 6833; 1785 > 10 counts
cs_only <- setdiff(cs_rix6, union(sea_rix6, gecco_names))
    ## 2560; 2618 > 10 counts

all_match <- intersect(intersect(sea_rix6, cs_rix6), gecco_names)
kmerize <- c(sea_gecco, sea_only)
kmerize <- allPO_sea$gene[which(!allPO_sea$gene %in% rownames(dds_cs$rix1))]

sea_only_cts <- data.frame(counts(dds_sea$rix6)[which(rownames(counts(dds_sea$rix6)) %in% sea_only),])
sea_only_cts$gene = as.character(rownames(sea_only_cts))
sea_only_cts<- gather(sea_only_cts, pup, counts, -gene)
#rownames(sea_only_cts[order(-rowMeans2(sea_only_cts))[1:20],])

cs_only_cts <- data.frame(counts(dds_cs$rix6)[which(rownames(counts(dds_cs$rix6)) %in% cs_only),])
cs_only_cts$gene = as.character(rownames(cs_only_cts))
cs_only_cts<- gather(cs_only_cts, pup, counts, -gene)
#rownames(cs_only_cts[order(-rowMeans2(cs_only_cts))[1:20],])

gecco_only_cts <- data.frame(gecco_sum[which(gecco_sum$gene %in% gecco_only),])
#gecco_only_cts$gene[order(-gecco_only_cts$counts_gecco)[1:20]]

all_match_sea <- data.frame(counts(dds_sea$rix6, normalized=T)[which(rownames(counts(dds_sea$rix6)) %in% all_match),])
all_match_sea$gene = rownames(all_match_sea)
all_match_sea<- gather(all_match_sea, pup, counts, -gene) %>%
  mutate(method = "sea")


all_match_cs <- data.frame(counts(dds_cs$rix6, normalized=T)[which(rownames(counts(dds_cs$rix6)) %in% all_match),])
all_match_cs$gene = rownames(all_match_cs)
all_match_cs <-gather(all_match_cs, pup, counts, -gene) %>%
  mutate(method = "cs")

all_match_gecco <- gecco_uniq[which(gecco_uniq$gene %in% all_match),c(4,2,3)]
all_match_gecco$method = "gecco"
colnames(all_match_gecco) <- colnames(all_match_cs)


all_match_counts <- rbind(all_match_gecco, all_match_sea, all_match_cs)
cor(all_match_cs$counts, all_match_sea$counts)

all_match_counts %>% group_by(gene, method) %>%
  mutate(counts = as.numeric(counts)) %>% 
  summarize(mean=mean(counts, na.rm=T)) %>%
  spread(method, mean) -> matched_means

cor(matched_means[,-1])
facs <- colMeans(matched_means[,-1])/min(colMeans(matched_means[,-1]))
matched_means[,-1] <- do.call("cbind", sapply(1:3, function(x) matched_means[,-1][,x] / facs[x]))

cor.test(matched_means$cs, matched_means$sea)
cor.test(matched_means$gecco, matched_means$sea)
cor.test(matched_means$gecco, matched_means$cs)

ggplot(matched_means, aes(x=cs, y=sea)) + # %>% filter(!gene %in% c("Fth1", "Lars2", "Sspn", "Cmss1", "Il31ra", "Jarid2"))
  geom_point() + 
  #scale_y_log10() +
  theme_minimal()
ggplot(matched_means, aes(x=gecco, y=sea)) + 
  geom_point() + 
  #scale_y_log10() +
  theme_minimal()
ggplot(matched_means, aes(x=gecco, y=cs)) + 
  geom_point() + 
  #scale_y_log10() +
  theme_minimal()

#, all_match_gecco$counts)

################

all_match <- intersect(rownames(sea_rix6_VDDSTD_sig), rownames(cs_rix6_VDDSTD_sig))
match_sea_rix6_VDDSTD <- sea_rix6_VDDSTD_sig[which(rownames(sea_rix6_VDDSTD_sig) %in% all_match),]
#ref_rix6_VDDSTD[which(rownames(ref_rix6_VDDSTD) %in% all_match),]
match_cs_rix6_VDDSTD <- cs_rix6_VDDSTD_sig[which(rownames(cs_rix6_VDDSTD_sig) %in% all_match),]
match_cs_rix6_VDDSTD[order(-match_cs_rix6_VDDSTD$log2FoldChange),]
match_sea_rix6_VDDSTD[order(-match_sea_rix6_VDDSTD$log2FoldChange),]
#compare_counts <- list()
#pup <- 2147

#counts_sea <- data.frame(counts_sea = counts(dds_sea$rix6)[which(rownames(dds_sea$rix6) %in% all_match), paste(pup)])
#counts_sea$Gene.Name <- rownames(counts_sea)

#counts_ref <- data.frame(counts_ref = counts(dds_ref$rix6)[which(rownames(dds_ref$rix6) %in% all_match), paste(pup)])
#counts_ref$Gene.Name <- rownames(counts_ref)

#counts_cs <- data.frame(counts_cs = counts(dds)[which(rownames(dds) %in% all_match), paste(pup)])
#counts_cs$Gene.Name <- rownames(counts_cs)

#gecco_uniq %>% filter(gene %in% all_match) -> gecco_mapped
#gecco_counts <- gecco_mapped[which(gecco_mapped$Gene.Name %in% all_match), -1]
#compare <- left_join(gecco_counts, counts_sea, by="Gene.Name") %>%
#  left_join(counts_ref, by="Gene.Name") %>%
#  left_join(counts_cs, by="Gene.Name")

(all_match_counts %>% filter(gene%in%all_match) %>% 
  group_by(method, gene) %>% 
  summarise(mean=mean(counts)) %>%
  filter(method != "gecco") %>%
  arrange(-mean))$gene -> top_genes
top_genes <- unique(top_genes)

fac <- all_match_counts %>% group_by(method) %>% summarise(mean=mean(counts)) %>%
  mutate(fac = mean / min(mean))

#all_match_counts$fcounts <- apply(all_match_counts, 1, function(x) x["counts"] / fac$fac[match(x["method"], fac$method)])
  
for(i in 1:3){#nrow(all_match_counts)){
  m <- fac$method[i]
  all_match_counts[which(all_match_counts$method == m),"fcounts"] = all_match_counts[which(all_match_counts$method == m),"counts"] / fac$fac[i]
}
  
  
n=30
compare <- all_match_counts %>% filter(gene %in% top_genes[1:n])
compare$gene <- factor(compare$gene, levels=top_genes[1:n])


ggplot(compare, aes(x=gene, y=fcounts, col=method)) + 
  geom_point(position=position_jitter(w=0.2,h=0)) + 
  #geom_line(data=means, aes(x=gene, y=mean, col=method)) +
  theme_minimal() +
  scale_y_log10() + 
  stat_summary(aes(y=fcounts,group=method), geom="line", fun.y = 'mean') +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))# + 

compare %>% group_by(gene, method) %>% summarize(mean=mean(fcounts)) -> means

ggplot(means, aes(x=gene, y=mean, col=method)) + 
  geom_point() + 
  theme_minimal() +  
  stat_summary(aes(y=mean,group=method), geom="line", fun.y = 'mean') +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))# + 

rank=which(top_genes == "Atp1a3")
#all_match_counts$gene <- factor(all_match_counts$gene, levels=top_genes)
title =paste("Counts of",top_genes[rank],"in RIX 6")
plotMat <- all_match_counts %>% filter(gene == top_genes[rank]) %>%
  mutate(pup = gsub("^X","",pup), fcounts = as.numeric(paste(fcounts)))
plotMat$Diet <- apply(plotMat, 1, function(x) if(!x["method"] == "gecco"){
  as.character(colData(dds_cs$rix6)$Diet[which(colData(dds_cs$rix6)$Pup.ID == x["pup"])])
} else {
  "gecco"
})
ggplot(plotMat, aes(x=Diet, y=fcounts, col=method)) +    
  geom_point(position=position_jitter(w=0.05,h=0)) + 
  #scale_y_log10() + #breaks=c(25,100,400)
  ggtitle(title) + 
  theme_minimal()


############################
Npas4
Dok3
Slc6a9
################  KMERIZE  #####################

to_kmer <- allPO_sea[which(!allPO_sea$gene %in% rownames(dds_cs$rix1)), ]
kmer_list <- list()

for(rix in unique(allPO_sea$rix)){
  gene=to_kmer$gene[which(to_kmer$rix == rix)]
  fa_file = file.path(dir, paste0("../pseudogenomes/scallop_transcripts/RIX",gsub("[a-z]", "", rix),"_unique.fa"))
  
  kmer_list[[rix]] <- kmerize_it(gene=gene, fa_file=fa_file, ref_seq=seq_det_full)
}

for(rix in names(kmer_list)){
  print <- unlist(kmer_list[[rix]]$gen_kmers)
  write.table(print, paste0(dir, "/../variant_data/csbio_rna_expression/novel_sea_kmers_",rix,".txt"),
              quote = F, row.names=F, col.names = F)
}


kmerize_it <- function(gene, klen = 25, gap = 3, fa_file="unique.fa",
                       ref_seq=NULL){
  if(is.null(ref_seq)) ref_seq = data.frame(gene=NA)
  rix_fa <- readDNAStringSet(fa_file)
  
  rix_seq <- lapply(names(rix_fa) , function(x) as.character(rix_fa[[x]]))
  names(rix_seq) = names(rix_fa)      
  rix_seq <- data.frame(gene =  names(rix_fa), 
                        sequence = unlist(rix_seq), stringsAsFactors = F)
  rix_seq$sequence <- as.character(rix_seq$sequence)
  kmer_genes <- gene[which(gene %in% union(rix_seq$gene, ref_seq$gene))]
  which_fa <- unlist(lapply(kmer_genes, function(x) ifelse(x %in% rix_seq$gene, 1, 0)))
  get_kmers <- sapply(1:length(which_fa), function(x) ifelse(which_fa[x] == 0, ref_seq$sequence[which(ref_seq$gene == kmer_genes[x])],
                                                             rix_seq$sequence[which(rix_seq$gene == kmer_genes[x])]))
  get_kmers <- data.frame(genes=kmer_genes,
                          seq = unlist(get_kmers), stringsAsFactors = F)
  
  #keep <- which(rownames(sea_rix6_VDDSTD_sig) %in% kmerize)
  #get_kmers <- get_kmers[keep,]
  
  gen_kmers <- lapply(get_kmers$seq, function(x){
    len=nchar(x)
    chars=unlist(strsplit(x,""))
    count=1
    p <- list()
    for(i in 1:(floor(len/gap)-floor(klen/gap))){
      max <- ifelse(count < len, count+klen-1, count)
      p[[i]] <- paste(chars[seq(count,max)], collapse="")
      count = count+gap
    }
    unlist(p)
  })
  names(gen_kmers) <- get_kmers$genes
  return(list(gen_kmers = gen_kmers, get_kmers=get_kmers))
}

##########   CHECK NOVEL COUNTS   #############

rix="rix1"


novel_sal_counts <- read.csv(paste0(dir, "/../variant_data/csbio_rna_expression/salmon_counts_",rix,"_5jun2019.csv"),
                             header=F)
names(novel_sal_counts) <- c("Pup.ID", "seq","dir","counts")
novel_sal_counts$Pup.ID <- do.call("rbind", strsplit(as.character(novel_sal_counts$Pup.ID), "/|_"))[,8]
first_kmers <- unlist(lapply(kmer_list[[rix]]$gen_kmers, function(x) x[[1]]))
novel_sal_counts$gene <- NA
first_pup <- novel_sal_counts %>% filter(Pup.ID == unique(novel_sal_counts$Pup.ID)[1])
for(i in 1:length(first_kmers)){
  first_pup$gene[which(first_pup$seq %in% first_kmers[i])] <- as.character(names(first_kmers)[i])
}

intervals <- which(!is.na(first_pup$gene))
for(i in 1:length(intervals)){
  if(i == length(intervals)){
    first_pup$gene[intervals[i]:nrow(first_pup)] <- first_pup$gene[intervals[i]]
  } else if(intervals[i+1]-intervals[i] > 1){
    first_pup$gene[intervals[i]:(intervals[i+1]-1)] <- first_pup$gene[intervals[i]]
  }
}

#sapply(1:27, function(check)
#all.equal(novel_sal_counts$seq[1:nrow(first_pup)], novel_sal_counts$seq[(check*nrow(first_pup)+1):(check*nrow(first_pup)+nrow(first_pup))])
#)

novel_sal_counts$gene <- rep(first_pup$gene, length(unique(novel_sal_counts$Pup.ID)))
novel_sal_counts$ordered <- rownames(novel_sal_counts)
novel_sal_counts %>% group_by(Pup.ID, gene, seq) %>%
  mutate(sum=sum(counts)) %>% 
  filter((as.numeric(ordered)-1) %% 9 == 0) %>%
  arrange(ordered) -> novel_sal_sums

print <- as.character((novel_sal_sums %>% filter(Pup.ID == "0379"))$seq)
write.table(print, paste0(dir, "/../variant_data/csbio_rna_expression/check_rix1_kmers_in_dna.txt"),
            quote = F, row.names=F, col.names = F)




## trim 50bp and try everything
## vst then do linear thing
########### Compare variances of genes for each method ############

cs_rowvars <- lapply(dds_cs, function(x) apply(counts(x, normalized=T), 1, function(y) ifelse(mean(y) > 10, sd(y)/mean(y), NA)))
sea_rowvars <- lapply(dds_sea, function(x) apply(counts(x, normalized=T), 1, function(y) ifelse(mean(y) > 10, sd(y)/mean(y), NA)))
unlist(lapply(cs_rowvars, function(x) mean(x, na.rm=T)*100))
unlist(lapply(sea_rowvars, function(x) mean(x, na.rm=T)*100))


#################   VDD vs other analysis   ####################
colData(dds_sea$rix6)
dds_sea$rix6$DietGroup <- as.factor(ifelse(dds_sea$rix6$Diet == "Vitamin D Deficient", "VDD", "Other"))
design(dds_sea$rix6) =  ~ DietGroup + PO + DietGroup:PO
dds_sea_VDDOth <- DESeq(dds_sea$rix6)
res_sea_VDDOth <- results(dds_sea_VDDOth, contrast=c("DietGroup", "VDD", "Other"))
res_sea_VDDOth <- res_sea_VDDOth[order(res_sea_VDDOth$padj),]
d <- plotCounts(dds_sea_VDDOth, gene=rownames(res_sea_VDDOth)[6], 
                intgroup="DietGroup")#, returnData = T)
overlap <- intersect(rownames(sea_rix6_VDDSTD_sig), rownames(res_sea_VDDOth[which(res_sea_VDDOth$padj < pval),]))
onlyOth <- setdiff(rownames(res_sea_VDDOth[which(res_sea_VDDOth$padj < pval),]), rownames(sea_rix6_VDDSTD_sig))

length(rownames(res_sea_VDDOth[which(res_sea_VDDOth$padj < pval),]))
d <- plotCounts(dds_sea_VDDOth, gene=onlyOth[5], 
                intgroup="DietGroup")

#######

dds_cs$rix6$DietGroup <- as.factor(ifelse(dds_cs$rix6$Diet == "Vitamin D Deficient", "VDD", "Other"))
design(dds_cs$rix6) =  ~ DietGroup + PO + DietGroup:PO
dds_cs_VDDOth <- DESeq(dds_cs$rix6)
res_cs_VDDOth <- results(dds_cs_VDDOth, contrast=c("DietGroup", "VDD", "Other"))
res_cs_VDDOth <- res_cs_VDDOth[order(res_cs_VDDOth$padj),]
d <- plotCounts(dds_cs_VDDOth, gene=rownames(res_cs_VDDOth)[6], 
                intgroup="DietGroup")#, returnData = T)
overlap <- intersect(rownames(cs_rix6_VDDSTD_sig), rownames(res_cs_VDDOth[which(res_cs_VDDOth$padj < pval),]))
onlyOth <- setdiff(rownames(res_cs_VDDOth[which(res_cs_VDDOth$padj < pval),]), rownames(cs_rix6_VDDSTD_sig))

length(rownames(res_cs_VDDOth[which(res_cs_VDDOth$padj < pval),]))
d <- plotCounts(dds_cs_VDDOth, gene=onlyOth[5], 
                intgroup="DietGroup")


nrow(res_cs_VDDOth[which(res_cs_VDDOth$padj < pval),])

overlap_cs_sea <- intersect(rownames(res_cs_VDDOth)[which(res_cs_VDDOth$padj < pval)], 
          rownames(res_sea_VDDOth)[which(res_sea_VDDOth$padj < pval)])
onlySeaOth <- setdiff(rownames(res_sea_VDDOth[which(res_sea_VDDOth$padj < pval),]), rownames(sea_rix6_VDDSTD_sig)) 
onlyCSOth <- setdiff(rownames(res_cs_VDDOth[which(res_cs_VDDOth$padj < pval),]), rownames(cs_rix6_VDDSTD_sig))

d <- plotCounts(dds_sea_VDDOth, gene=onlySeaOth[2], 
                intgroup="DietGroup")#, returnData = T)
d <- plotCounts(dds_cs_VDDOth, gene=onlySeaOth[2], 
                intgroup="DietGroup")#, returnData = T)
d <- plotCounts(dds_cs_VDDOth, gene=onlyCSOth[3], 
                intgroup="DietGroup")#, returnData = T)
d <- plotCounts(dds_sea_VDDOth, gene=onlyCSOth[3], 
                intgroup="DietGroup")#, returnData = T)

#################   checking/visualizing data   ####################

#norm_VDDOth <- assay(normTransform(dds_VDDOth))
vsd <- vst(dds_sea$rix6, blind=FALSE)
vsd_means <- data.frame(sd = unlist(apply(assay(vsd), 1, sd)),
                        means = unlist(apply(assay(vsd), 1, mean)),
                        genes = rownames(assay(vsd)))
vsd_means %>% arrange(-sd) %>% mutate(rank=rank(means)) -> vsd_means
vsd_means$density <- get_density(vsd_means$rank, vsd_means$sd)
ggplot(vsd_means) + geom_point(aes(rank, sd, color = density*nrow(vsd_means))) + theme_minimal()

plotPCA(vsd[which(res_sea_VDDSTD$rix6$padj < 0.1),], intgroup=c("Diet"))
pcaData <- plotPCA(vsd, intgroup=c("Diet", "PO"), returnData=T)
#[which(res_lst$rix6$padj < 0.1),]

percentVar <- round(100 * attr(pcaData, "percentVar"))
pcaData$PO <- factor(pcaData$PO)
ggplot(pcaData, aes(PC1, PC2, color=Diet, shape=PO)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  theme_minimal() +
  coord_fixed()

results(dds_sea$rix6, name=c("Diet_Vitamin.D.Deficient_vs_Standard"))


##########
require(MASS)
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

set.seed(1)
dat <- data.frame(
  x = c(
    rnorm(1e4, mean = 0, sd = 0.1),
    rnorm(1e3, mean = 0, sd = 0.1)
  ),
  y = c(
    rnorm(1e4, mean = 0, sd = 0.1),
    rnorm(1e3, mean = 0.1, sd = 0.2)
  )
)

##### DON'T REDO!!! #####
dds_cs <- list()
cs_map <- read.csv("C:/Users/Kathie/Dropbox (ValdarLab)/variant_data/csbio_rna_expression/map_invariant_probes.csv")

for(rix in unique(samples$RRIX)){
  cs_counts <- read.csv(paste0("C:/Users/Kathie/Dropbox (ValdarLab)/variant_data/csbio_rna_expression/invariantCounts_rix",rix,
                               "_22May2019.csv"), header=F)
  cs_counts$pups <- do.call("rbind", strsplit(as.character(cs_counts[,1]), "_|/"))[,8]
  cs_counts %>% left_join(cs_map, by=c("V2" = "Kmer")) -> joined_counts
  
  joined_counts %>% group_by(pups, Gene.Name) %>% tally() -> n_na
  joined_counts %>% group_by(pups, Gene.Name) %>% summarize(sum = sum(V4)) -> sums
  counts <- sums %>% left_join(n_na, by=c("pups","Gene.Name"))
  
  pup_list <- list()
  pups <- samples_use$Pup.ID[which(samples_use$RRIX == rix)]
  for(p in as.character(pups)){
    pup_list[[p]] <- counts$sum[which(as.numeric(paste(counts$pups)) == p)]
    names(pup_list[[p]]) <- counts$Gene.Name[which(as.numeric(paste(counts$pups)) == p)]
    #count_mat[,p] <- unlist(sapply(1:nrow(count_mat), function(x){
    #  counts[intersect(which(counts$pups == p), which(counts$Gene.Name == rownames(count_mat)[x])), "sum"]  }))
  }
  
  count_mat <- do.call("cbind", pup_list)
  colnames(count_mat) <- pups
  rownames(count_mat) <- names(pup_list[[1]])
  
  samples_use %>% filter(RRIX == rix) %>%
    mutate(Pup.ID = as.character(Pup.ID)) -> samples_de
  samples_de <- samples_de[match(samples_de$Pup.ID, colnames(count_mat)),]
  
  #rownames(count_mat) <- unique(counts$Gene.Name)
  colnames(count_mat) <- gsub("^0", "", colnames(count_mat))
  count_mat <- count_mat[,which(colnames(count_mat) %in% samples_de$Pup.ID)]
  
  dds <- DESeqDataSetFromMatrix(countData = count_mat,
                                colData = samples_de,
                                design = ~ Diet + PO + Diet:PO)
  dds
  dds <- DESeq(dds) 
  dds <- nbinomWaldTest(dds)
  dds_cs[[paste0("rix",rix)]] <- dds
}


################## OLD CODE ########################


compare_counts <- list()
pup <- 2272
rix <- samples_use$RRIX[which(samples_use$Pup.ID == pup)]

tmp  <- read.table(paste0(dir, "/../variant_data/csbio_rna_expression/rna_exp_Pup.ID_",pup,".csv"), sep="\n")
tmp <- apply(tmp, 1, function(x) {
  z <- unlist(strsplit(as.character(x), ","))
  c(z, rep(NA, 43-length(z)))
})

tmp <- t(tmp)
colnames(tmp) <- paste(tmp[1,])
tmp <- tmp[-1,]
tmp <- data.frame(tmp)
tmp[,grep("probe", colnames(tmp))] <- apply(tmp[,grep("probe", colnames(tmp))], 2, as.numeric)
tmp$sum_cs <- apply(tmp[,grep("probe", colnames(tmp))], 1, function(x) sum(x, na.rm=T))
tmp$na_cs <- apply(tmp[,grep("probe", colnames(tmp))], 1, function(x) length(which(is.na(x) == F)))
tmp$counts_cs <- tmp$sum_cs/(tmp$na_cs/2)
tmp <- tmp[order(-tmp$counts_cs),]
tmp <- dplyr::select(tmp, -contains("probe"))

## how many counts from CSBio not in sea?
length(unique(tmp$ensemblID[which(!tmp$ensemblID[which(tmp$sum_cs > 10)] %in% seq_det$ensemblID)])) 
tmp[which(!tmp$ensemblID[which(tmp$sum_cs > 10)] %in% seq_det$ensemblID),]
## 385
length(unique(seq_det$ensemblID[which(!seq_det$ensemblID %in% tmp$ensemblID)]))
## 20193

sea <- data.frame(gene=rownames(dds_sea[[paste0("rix",rix)]]),
                  counts_sea=as.numeric(counts(dds_sea[[paste0("rix",rix)]])[,paste(pup)], normalized=T))
sea <- data.frame(gene=rownames(dds_sea[[paste0("rix",rix)]]),
                  counts_sea=as.numeric(apply(counts(dds_sea[[paste0("rix",rix)]], normalized=T), 1, mean)))
sea <- sea[order(-sea$counts),]

length(which(!sea$gene %in% tmp$gene))
length(which(!tmp$gene %in% sea$gene))


#ensemblID = unlist(strsplit(rownames(dds_ref[[paste0("rix",rix)]]),"[.]"))[c(T,F)]
ref <- data.frame(GENEID=as.character(rownames(dds_ref[[paste0("rix",rix)]])),
                  counts_ref=as.numeric(counts(dds_ref[[paste0("rix",rix)]])[,paste(pup)]))
ref <- ref[order(-ref$counts),]
ref <- left_join(ref, seq_det[,c("GENEID","ensemblID","gene")], by="GENEID") %>%
  dplyr::select(-GENEID) %>%
  distinct()
length(which(!ref$gene %in% tmp$gene))
length(which(!tmp$gene %in% ref$gene))

compare <- full_join(tmp, sea, by="gene") %>%
  full_join(gecco_sum, by="ensemblID") %>%
  left_join(ref, by="ensemblID")
compare %>% mutate(rank_cs = rank(-counts_cs), rank_sea=rank(-counts_sea), 
                   rank_gecco=rank(-counts_gecco),  rank_ref=rank(-counts_ref)) -> compare
compare[,grep("counts", colnames(compare))] <- scale(compare[,grep("counts", colnames(compare))], center=F, 
                                                     scale = apply(compare[,grep("counts", colnames(compare))], 2, 
                                                                   function(x) mean(x, na.rm=T)/mean(compare$counts_ref, na.rm=T)))

#compare <- cbind(compare, compare_orig[,-grep("counts", colnames(compare_orig))])
compare_short <- compare %>% #filter(!is.na(counts_cs), !is.na(counts_sea), !is.na(counts_gecco), !is.na(counts_ref)) %>%
  arrange(rank_ref)
wilcox.test(compare_short$counts_cs, compare_short$counts_gecco)
wilcox.test(compare_short$counts_sea, compare_short$counts_ref)


n=20
plot_compare <- gather(compare_short, key = "method", value = "counts", counts_cs, counts_sea, counts_gecco, counts_ref) %>%
  #gather(key = "method", value = "rank", rank_cs, rank_sea) %>%
  filter(gene %in% compare_short$gene[1:n]) %>%
  mutate(method = gsub("counts_", "", method))

plot_compare$rank <- apply(plot_compare, 1, function(x) as.numeric(paste(
  ifelse(x["method"] == "cs", x["rank_cs"], 
         ifelse(x["method"] == "sea", x["rank_sea"], 
                ifelse(x["method"] == "gecco", x["rank_gecco"], x["rank_ref"]) ))))) 
plot_compare %>% dplyr::select(-rank_sea, -rank_cs, -rank_gecco) -> plot_compare

methods = c("CSBio", "Gecco", "Scallop-Salmon", "Salmon-mm10")
title = paste("Top", n, "genes by", methods[4])
ggplot(plot_compare, aes(x=gene, y=counts, fill=method)) +    
  geom_bar(stat="identity", position=position_dodge()) +
  #geom_segment(aes(x=gene, xend=gene, y=rep(1, 50), yend=counts_cs+1), col="red") + 
  #geom_segment(aes(x=gene, xend=gene, y=rep(1, 50), yend=counts_sea+1), col="blue") + 
  #scale_y_log10() + #breaks=c(25,100,400)
  ggtitle(title) + 
  theme_minimal() +  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))



#compare[order(compare$rank_sea)[1:50],]
wilcox.test(compare$ref, compare$sea)
compare_counts[[paste0("Pup.ID_",pup)]] <- full_join(tmp, sea, by="gene")
compare_counts[[paste0("Pup.ID_",pup)]][order(-compare_counts[[paste0("Pup.ID_",pup)]]$counts_sea)[1:50],]

length(which(csbio_counts[[paste0("Pup.ID_",pup)]]$gene %in% sea_counts[[paste0("Pup.ID_",pup)]]$gene))





#wilcox.test(compare$counts_cs, compare$counts_gecco)
#wilcox.test(compare$counts_sea, compare$counts_ref)

plot_compare <- gather(compare, key = "method", value = "counts", counts_cs, counts_sea, counts_gecco, counts_ref) %>%
  #gather(key = "method", value = "rank", rank_cs, rank_sea) %>%
  #filter(gene %in% compare_short$gene[1:n]) %>%
  mutate(method = gsub("counts_", "", method))

plot_compare$rank <- apply(plot_compare, 1, function(x) as.numeric(paste(
  ifelse(x["method"] == "cs", x["rank_cs"], 
         ifelse(x["method"] == "sea", x["rank_sea"], 
                ifelse(x["method"] == "gecco", x["rank_gecco"], x["rank_ref"]) ))))) 
plot_compare %>% dplyr::select(-rank_sea, -rank_cs, -rank_gecco) -> plot_compare

methods = c("CSBio", "Gecco", "Scallop-Salmon", "Salmon-mm10")
###### to get DESeq2 data from scratch ######
txmap <- read.table(file.path(dir, "mm10/scallop/gtf_merge/gffall.merged_scallop.gtf.tmap"), header=T)
tx2gene <- data.frame(TXNAME=txmap$qry_id, GENEID=txmap$ref_gene_id)   
head(tx2gene)
txi <- tximport(use_files, type="salmon", tx2gene=tx2gene)

dds_lst <- list()
res <- list()
keepCounts_list <- list()
res0.1 <- list()
for(i in unique(coldata$RRIX)){
  coldata_use <- coldata %>% filter(RRIX == i)
  use <- which(colnames(txi$counts) %in% coldata_use$names)
  txi_use <- list(abundance=txi$abundance[,use], 
                  counts=txi$counts[,use], 
                  length=txi$length[,use], 
                  countsFromAbundance=txi$countsFromAbundance)
  dedata <- DESeqDataSetFromTximport(txi_use,
                                     colData = coldata_use,
                                     design = ~ Diet + PO )#+ Diet:PO
  dds <- DESeq(dedata) 
  dds <- nbinomWaldTest(dds)
  #res <- results(dds, c("PO_cat","b","a"))
  res <- results(dds, name="PO")
  res[[paste0("rix",i)]] <- res
  res0.1[[paste0("rix",i)]] <- res[which(res$padj < 0.1),]
  keepCounts <- rowSums(counts(dds)) >= 10
  keepCounts <- counts(dds)[which(keepCounts==T),]
  keepCounts <- keepCounts[order(-rowSums(keepCounts)),]
  keepCounts_list[[paste0("rix",i)]] <- keepCounts
  dds_lst[[paste0("rix",i)]] <- dds
}

###########################################################

dir <- "C:/Users/Kathie/Dropbox (ValdarLab)/variant_data/csbio_rna_expression/"
cs_files <- list.files(dir, full.names = T)
cs <- lapply(cs_files, read.csv)
cs_counts<- lapply(cs, function(tmp){
  #tmp <- apply(tmp, 1, function(x) {
  #  z <- unlist(strsplit(as.character(x), ","))
  #  c(z, rep(NA, 43-length(z)))
  #)
  
  #tmp <- t(tmp)
  #colnames(tmp) <- paste(tmp[1,])
  #tmp <- tmp[-1,]
  #tmp <- data.frame(tmp)
  tmp[,grep("probe", colnames(tmp))] <- apply(tmp[,grep("probe", colnames(tmp))], 2, function(x) as.numeric(paste(x)))
  tmp$sum_cs <- apply(tmp[,grep("probe", colnames(tmp))], 1, function(x) sum(x, na.rm=T))
  tmp$na_cs <- apply(tmp[,grep("probe", colnames(tmp))], 1, function(x) length(which(is.na(x) == F)))
  tmp$counts_cs <- tmp$sum_cs/(tmp$na_cs/2)
  #tmp <- tmp[order(-tmp$counts_cs),]
  dplyr::select(tmp, -contains("probe"))
})
unlist(sapply(2:length(cs_counts), function(x) all.equal(cs_counts[[1]]$gene, cs_counts[[x]]$gene)))

cs_mat <- lapply(cs_counts, function(x) x[,"sum_cs"])
cs_mat <- do.call("cbind", cs_mat)
rownames(cs_mat) <- cs_counts[[1]]$gene
colnames(cs_mat) <- gsub("rna_exp_|.csv", "", do.call("rbind", strsplit(cs_files, "/"))[,9])

##### THIS IS NOT CURRENTLY WORKING #####
ref_gtf <- read.table(file.path(dir, "../pseudogenomes/ensembl_merge/Mus_musculus.GRCm38.96_nohead.gtf"), sep="\\")
ref_gtf <- apply(ref_gtf, 1, function(x) strsplit(x, "\t|;"))

#### tx2gene for ref mm10 ####                                                                                                                         

ref_fa <- readDNAStringSet(file.path(dir, "../pseudogenomes/ensembl_merge/Mus_musculus.GRCm38.cdna.all.fa"))
txdb <- makeTxDbFromGFF(file.path(dir, "../pseudogenomes/ensembl_merge/Mus_musculus.GRCm38.96.gtf"), format="gtf", 
                        dataSource="Ensembl", organism="Mus musculus")
# k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")
k <- keys(txdb, keytype = "CDSNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")
txdump <- as.list(txdb)

#k <- keys(txdb, keytype = "TXNAME")
#cols=columns(txdb)#[grep("EXON|GENE|TX", columns(txdb))]
#tx <- select(txdb, k, "TXNAME","GENEID")

unlist(strsplit(rownames(dds_ref[[1]]), "[.]"))[c(T,F)][(which(!unlist(strsplit(rownames(dds_ref[[1]]), "[.]"))[c(T,F)] %in% seq_det$ensemblID))]
which(!rownames(dds_ref[[1]]) %in% as.character(seq_det$GENEID))




################  discarded from deseq2_analysis ########################
#  dists = apply(minf, 1, function(x){
#    min(abs(find_gene$End - as.numeric(x["seq.Position"])), abs(as.numeric(x["seq.Position"]) - find_gene$Start))
#  })
#  minf = minf[order(dists),]
#} else {
#  min_1=which.min(abs(find_gene$Start - masterKmers$seq.Position))
#  min_2=which.min(abs(find_gene$End - masterKmers$seq.Position))
#  if(min_1 != min_2){
#    if(min(abs(find_gene$Start - masterKmers$seq.Position)) < min(abs(find_gene$End - masterKmers$seq.Position))){
#      minf = rbind(masterKmers[min_1,], masterKmers[min_2,])
#    } else {
#      minf = rbind(masterKmers[min_2,], masterKmers[min_1,])
#    }
#  } else {
#    minf = masterKmers[min_1,]
#  }



kmer_compare = data.frame(t(c(1,1,1,1,1,1,1,1)))
#pd %>% mutate(Pup.ID = gsub("Pup.ID_", "", pd$Pup.ID)) %>% 
#  left_join(samples_use)
