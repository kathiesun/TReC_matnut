options(stringsAsFactors = FALSE)
library(tidyverse)
library(DESeq2)
library(GenomicFeatures)
library(Biostrings)
library(fdrtool)
library(gridExtra)

setwd("~/rna_seq/")
dir <- "C:/Users/Kathie/Dropbox\ (ValdarLab)"
source("deseq2/deseq2_functions.R")

######## READ IN DESEQ2 RESULTS ######## 
invar_info <- read.csv(file.path(dir,"variant_data/kmers_based_on_CC/kmers_to_run/GeneProbesInvariant.csv"))
invar_info = invar_info %>% dplyr::select(c("Gene.Name","Gene.ID")) %>% distinct()

#dds_string <- readRDS(file.path(dir, "de_results/dds_lst_string_sva_11oct2019.rds"))
dds_cs_old <- readRDS(file.path(dir, "de_results/dds_lst_cs-5aug2019.rds"))
dds_cs <- readRDS(file.path(dir, "de_results/dds_lst_cs_nointeraction-16oct2019.rds"))

dds_lst_files = list.files("C:/Users/Kathie/Dropbox (ValdarLab)/de_results/stringtie_inter_sva_combos",full.names = T)
string_files = dds_lst_files[grep("string_deseq", dds_lst_files)]
cs_files = dds_lst_files[grep("kmer_deseq", dds_lst_files)]

dds_string <- lapply(string_files, function(x) readRDS(x))
rixes = do.call("rbind", strsplit(string_files,"[/|_]"))
rixes = rixes[,grep("rix", rixes[1,])]
names(dds_string) = rixes

dds_cs_lst <- lapply(cs_files, function(x) readRDS(x))
rixes = do.call("rbind", strsplit(cs_files,"[/|_]"))
rixes = rixes[,grep("rix", rixes[1,])]
names(dds_cs_lst) = rixes

dds_cs = lapply(dds_cs_lst, function(x) x$sva_nointer)
names(dds_cs) = rixes

result_sv_nointer = lapply(dds_cs_lst, function(x) results(x$sva_nointer, contrast = c("PO_cat","a","b")))
result_nsv_nointer = lapply(dds_cs_lst, function(x) results(x$nosva_nointer, contrast = c("PO_cat","a","b")))
result_sv_inter = lapply(dds_cs_lst, function(x) results(x$sva_inter, name = "PO"))
result_nsv_inter = lapply(dds_cs_lst, function(x) results(x$nosva_inter, name = "PO"))

pdf("C:/Users/Kathie/Dropbox (ValdarLab)/figures_and_updates/images_meetings/sva_inter_kmer_comparisons_pval.pdf")
par(mfrow=c(2,2))
for(i in 1:length(result_sv_nointer)){
  hist(result_sv_nointer[[i]]$pvalue, 
       main = paste("P-values in",names(result_sv_nointer)[i],"SV1, no interaction"),xlab="")
  hist(result_nsv_nointer[[i]]$pvalue, 
       main = paste("P-values in",names(result_sv_nointer)[i],"no SV, no interaction"),xlab="")
  hist(result_sv_inter[[i]]$pvalue, 
       main = paste("P-values in",names(result_sv_nointer)[i],"SV1, w/ interaction"),xlab="")
  hist(result_nsv_inter[[i]]$pvalue, 
       main = paste("P-values in",names(result_sv_nointer)[i],"no SV, w/ interaction"),xlab="")
}
dev.off()

#dds_mm10 <- readRDS(file.path(dir, "de_results/dds_lst_mm10transcript-15sep2019.rds"))
names(dds_cs_old) <- tolower(names(dds_cs_old))
identical(names(dds_cs), names(dds_string))

###### sample data ######
problemPups <- c(1404, 1716, 1371, 569)
samples <- read.csv(file.path(dir, "variant_data", "2015-10_expression_pups.csv"))
rownames(samples) = samples$Pup.ID
samples$Diet <- gsub(" $", "", samples$Diet)
samples$Diet <- factor(samples$Diet, levels=c("Standard", "Low Protein","Methyl Enriched","Vitamin D Deficient"))
samples$PO <- ifelse(factor(gsub("[0-9]","", samples$RIX)) == "a", 0.5, -0.5)
samples$PO_cat <- factor(gsub("[0-9]","",samples$RIX))

pups = unlist(lapply(rem_seg_dds, colnames))
pups = gsub("[A-Z|_|.]","",toupper(pups))
if(length(which(as.numeric(pups) %in% problemPups)) > 0) pups = as.numeric(pups)[-which(as.numeric(pups) %in% problemPups)]

samples_use <- samples %>% 
  filter(Pup.ID %in% as.numeric(pups)) %>%
  dplyr::select(one_of("Pup.ID","RRIX","Diet","PO","PO_cat","Breeding.Batch","Behavior.Batch","Dam.Line")) %>%
  mutate(Breeding.Batch  = as.factor(Breeding.Batch), Behavior.Batch = as.factor(Behavior.Batch))


##### GC content analysis ########
#txdb <- makeTxDbFromGFF("/home/kathie/Downloads/gencode.vM23.chr_patch_hapl_scaff.annotation.gtf.gz")
#saveDb(txdb, file=paste0(dir,"matnut_main/gencode.vM23.chr_patch_hapl_scaff.annotation.sqlite"))

txdb <- loadDb(paste0(dir,"matnut_main/gencode.vM23.chr_patch_hapl_scaff.annotation.sqlite"))
ebg <- exonsBy(txdb, by="gene")
head(names(ebg))

table(names(ebg) %in% seq_det_red$ensemblID)
table(seq_det_red$ensemblID[which(seq_det_red$gene %in% rownames(dds))] %in% names(ebg))
genes_in_set = seq_det_red$ensemblID[which(seq_det_red$gene %in% rownames(dds))]
keep = names(ebg[which(names(ebg) %in% genes_in_set)])
ebg = unlist(ebg)
ebg <- ebg[which(names(ebg) %in% keep)]
e <- ebg[which(names(ebg) == names(ebg)[1])]
#library(rafalib)
plotRanges <- function(e) {
  l <- length(e)
  r <- ranges(range(e))
  nullplot(start(r), end(r), 0, l+1)
  segments(start(e), 1:l, end(e), 1:l, lwd=5)
}
plotRanges(e)
plotRanges(GenomicRanges::reduce(e))
ebg.red <- GenomicRanges::reduce(ebg)
library(BSgenome.Mmusculus.UCSC.mm10)
dna <- extractTranscriptSeqs(Mmusculus, ebg)
all(sum(width(ebg.red)) == width(dna))

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
gene <- gsub("gene:", "", seq_det[,3])

gene_sym <- do.call("rbind", strsplit(seq_det[,4], "[:]"))                                                                  
seq_det <- cbind(tx2gene, transcripts, chr, gene, gene_sym)                                                                                                      
colnames(seq_det) = c("TXNAME","GENEID", "transcript","tr_num", "del","build",
                      "chr","start","end","dir",#"del2",
                      "ensemblID",#"gene_num",
                      "del3","gene")                                                                                                                             
seq_det <- seq_det[,-grep("del", colnames(seq_det))]
seq_det <- data.frame(seq_det)        
seq_det %>% mutate(start = as.numeric(paste(start)),                                                                         
                   end = as.numeric(paste(end)),
                   sequence = paste(ref_fa)) -> seq_det_full
seq_det_full %>% dplyr::select("chr", "start", "end", "ensemblID", "gene") %>%
  group_by(chr, gene) %>%
  arrange(start) %>%
  dplyr::slice(1, n()) -> tmp

seq_det_red = unique(seq_det_full[c("chr", "ensemblID", "gene")]) %>% 
  filter(chr %in% c(1:19,"X")) %>% 
  mutate(chr = factor(chr, levels=c(1:19,"X")))
seq_det_red$start = unlist(lapply(seq_det_red$gene, function(x)
  min(tmp$start[which(tmp$gene == x)])))
seq_det_red$end = unlist(lapply(seq_det_red$gene, function(x)
  max(tmp$end[which(tmp$gene == x)])))
seq_det_red = seq_det_red %>% arrange(chr, start)

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

################  KMER DATA  ###################

all_snps <- read.csv(file.path(dir, "variant_data/raw_data/SNPsInExonsTSL1.csv"), 
                     colClasses=c("FounderSDP"="character"))
allSDP = do.call("rbind", strsplit(all_snps$FounderSDP,""))
allSDP = data.frame(cbind(all_snps[, c("Chromosome","GeneID", "Gene", "rsId", "Position")], allSDP))
allSDP$chr_pos = paste0(allSDP$Chromosome, "_", allSDP$Position)
allSDP %>% group_by(Chromosome, Gene) %>% 
  dplyr::slice(1, n())

###### GENERATE RESULTS ######
dds_cs_compare = compare_diets_PO(dds_lst = dds_cs, ref_diet = "Vitamin.D.Deficient")
dds_cs_old_compare = compare_diets_PO(dds_lst = dds_cs_old, ref_diet = "Vitamin.D.Deficient")

dds_string_compare = compare_diets_PO(dds_lst = dds_string, ref_diet = "Vitamin.D.Deficient")

dds_mm10_compare = compare_diets_PO(dds_lst = dds_mm10, ref_diet = "Vitamin.D.Deficient")

#################   PO ANALYSIS   #########################

######## imprinted genes ######## 
meta <- read.csv(file.path(dir, "matnut_main/metaseq.csv"))
ie1 <- read.table(file.path(dir, "imprinted_genes/imprinted_genes_mott.txt"), header = T, stringsAsFactors = F)
ie2 <- read.table(file.path(dir, "imprinted_genes/imprinted_genes_jirtle.txt"), header = T, sep=",", stringsAsFactors = F)
ie3 <- read.table(file.path(dir, "imprinted_genes/GSE27016_imprinted_genes.readme"), header = F, sep=";",stringsAsFactors = F)
ie_genes = read.csv("C:/Users/Kathie/Dropbox (ValdarLab)/imprinted_genes/2014_05_14_allImprintedGenes.csv")

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
compare_ie = unique(c(compare_ie, ie_genes$mgi_symbol))
compare_ie = compare_ie[order(compare_ie)]
compare_ie = compare_ie[-which(compare_ie == "klf14")]

compare_ie[grep("MIR", toupper(compare_ie))] = toupper(compare_ie[grep("MIR", toupper(compare_ie))])
compare_ie = gsub("MIR[0-9]", "Mir-", compare_ie)
compare_ie = unique(gsub("MIR", "Mir", compare_ie))
compare_ie = compare_ie[-grep("@|[*]|SNORD", compare_ie)]
#write.csv(compare_ie, file.path(dir,"imprinted_genes/use_ie_gene_names.csv"), row.names = F)

## check imprinted genes
sig_imprint = lapply(rem_seg_res$res_PO, function(x) 
  do.call("rbind", lapply(compare_ie, function(y)
    if(length(grep(paste0("^",y,"$"), rownames(x))) > 0){
      tmp = data.frame(x[grep(paste0("^",y,"$"), rownames(x)),])
      tmp$gene = y
      tmp
    })
  ))
sig_imprint_string = lapply(dds_string_compare$res_PO, function(x) 
  do.call("rbind", lapply(compare_ie, function(y)
    if(length(grep(paste0("^",y,"$"), rownames(x))) > 0){
      tmp = data.frame(x[grep(paste0("^",y,"$"), rownames(x)),])
      tmp$gene = y
      tmp
    })
  ))

new_imprint = unlist(lapply(sig_imprint, function(x) rownames(x[which(x$padj < 0.05),])))
string_imprinted = do.call("rbind", lapply(sig_imprint_string, function(x) x[which(x$padj < 0.05),]))
rownames(string_imprinted) = NULL
plot_genes = new_imprint#(string_imprinted %>% arrange(padj))$gene

###### plot known imprinted genes that are significant for PO
pdf("C:/Users/Kathie/Dropbox (ValdarLab)/figures_and_updates/images_meetings/known_imprinted_PO_counts_18nov2019.pdf")
pdf("C:/Users/Kathie/Dropbox (ValdarLab)/figures_and_updates/images_meetings/fishcomb_top10_PO_counts_18nov2019.pdf")
pdf("C:/Users/Kathie/Dropbox (ValdarLab)/figures_and_updates/images_meetings/multRIX_PO_counts_18nov2019.pdf")

par(mfrow=c(3,3))
for(i in 1:length(plot_genes)){
  gene = plot_genes[i]
  #rix=paste0("rix",c(4,4,6,8,8,9))[i]    #c(4,4,6,8,8,9)
  for(rix in names(rem_seg_dds)){
    if(gene %in% rownames(rem_seg_dds[[rix]])){
      pd <- plotCounts(rem_seg_dds[[rix]], gene=gene, intgroup = "PO",returnData = T)
      boxplot(count ~ as.factor(PO), data=pd, outline=F, ylim=c(min(pd$count), max(pd$count)),
              main=paste(gene, rix), xlab = "PO")      
      points(x=jitter(as.numeric(as.factor(pd$PO))), y=pd$count)  
    } else {
      plot(0,0)
    }
    
  }
  
  #anova(lm(data=pd, count ~ PO))
}
dev.off()

###### plot raw or adj p-value distributions for results
pdf("C:/Users/Kathie/Dropbox (ValdarLab)/figures_and_updates/images_meetings/raw_pvalues_cs_PO.pdf")
#par(mfrow=c(2,1))
for(i in 1:length(dds_cs_compare$res_PO)){
  hist(dds_cs_compare$res_PO[[i]]$pvalue ,
       main = paste("Raw p-values in",names(dds_cs_compare$res_PO)[i]),xlab="")
}
dev.off()

########  CHECK SIGNIFICANT GENE SDP'S  ########
#phased_par_CC_haplotype = readRDS(file.path(dir, "mini/phased_par_CC_haplotype_joined_sep2019.rds"))
haplofiles = list.files("C:/Users/Kathie/Dropbox (ValdarLab)/mini/pup_haplo_blocks_by_CC_parent", 
                        pattern="haploBlocks.rds", full.names = T)
phased_par_CC_haplotypes = lapply(haplofiles, readRDS)
names(phased_par_CC_haplotypes) = do.call("rbind", strsplit(haplofiles, "_"))[,7]
CC_lab = read.csv(file.path(dir, "variant_data/matched_v2_4jun2018.csv"))
CC_lab = CC_lab %>% dplyr::select(c("Pup.ID","RIX","CC.1","CC.2"))


### result objects: dds_string_compare, dds_cs_compare
## 15 nov 2019 objects: rem_seg_res (filt rowSum > 10, remove seg), dds_string_compare (no filt, keep regions)
sig_PO_old = check_sig_PO_SDP(PO_results_lst=dds_string_compare$res_PO, sample_data = samples_use, 
                                 phased_CC_haplotypes=phased_par_CC_haplotypes,
                                 sequence_details = seq_det_red, CC_labels=CC_lab)
sig_PO_new = check_sig_PO_SDP(PO_results_lst=rem_seg_res$res_PO, sample_data = samples_use, 
                                 phased_CC_haplotypes=phased_par_CC_haplotypes,
                                 sequence_details = seq_det_red, CC_labels=CC_lab)
sig_PO_cs = check_sig_PO_SDP(PO_results_lst=dds_cs_compare$res_PO, sample_data = samples_use, 
                             phased_CC_haplotypes=phased_par_CC_haplotypes,
                             sequence_details = seq_det_red, CC_labels=CC_lab)
sig_PO_mm10 = check_sig_PO_SDP(PO_results_lst=dds_mm10_compare$res_lst_PO, sample_data = samples_use, 
                               phased_CC_haplotypes=phased_par_CC_haplotypes,
                               SDP_data=allSDP, CC_labels=CC_lab)
#sig_PO_cs %>% filter(chr %in% c(1:19)) -> sig_PO_cs
#sig_PO_string %>% filter(chr %in% c(1:19)) -> sig_PO_string


## check combined analysis genes 
sig_PO_string_all = check_sig_PO_SDP(PO_results_lst=dds_string_compare$res_PO, sample_data = samples_use, 
                             phased_CC_haplotypes=phased_par_CC_haplotypes,
                             sequence_details = seq_det_red, CC_labels=CC_lab, alpha=1)

haps = dim(sig_PO_string[intersect(which(sig_PO_string$CC.1 == sig_PO_string$CC.2), which(!is.na(sig_PO_string$chr))),])[1]
haps = dim(sig_PO_string[intersect(which(sig_PO_string$CC.1 == sig_PO_string$CC.2), which(!is.na(sig_PO_string$chr))),])[1]
haps = dim(sig_PO_string[which(sig_PO_string$CC.1 == sig_PO_string$CC.2),])[1]

haps / length(sig_PO_string$gene)

## homozygous
sig_PO_old %>% filter(!is.na(chr), CC.1 != "", CC.2 != "") %>%
  filter(CC.1 == CC.2) -> PO_but_homozyg
length(unique(PO_but_homozyg$gene))/length(unique(sig_PO_old$gene))
sig_PO_new %>% filter(!is.na(chr), CC.1 != "", CC.2 != "") %>%
  filter(CC.1 == CC.2) -> PO_but_homozyg
length(unique(PO_but_homozyg$gene))/length(unique(sig_PO_new$gene))

## NAs
length(unique(sig_PO_old$gene[which(is.na(sig_PO_old$chr))]))/length(unique(sig_PO_old$gene))
length(unique(sig_PO_new$gene[which(is.na(sig_PO_new$chr))]))/length(unique(sig_PO_new$gene))

##
sig_PO_old %>% filter(CC.1 != CC.2) -> PO_hetero_old
length(unique(PO_hetero_old$gene))/length(unique(sig_PO_old$gene))
sig_PO_new %>% filter(CC.1 != CC.2) -> PO_hetero_new
length(unique(PO_hetero_new$gene))/length(unique(sig_PO_new$gene))
## 7740 + 3049 + 1318 > 11644??

PO_but_homozyg[-intersect(grep(",",PO_but_homozyg$CC.1),grep(",",PO_but_homozyg$CC.2)),]

######  WHICH SIGNIFICANT GENES SHOW UP IN MULTIPLE RIX #########  

pval = 0.05

old_PO_0.05 <- lapply(dds_string_compare$res_PO, function(x){
  x$gene = rownames(x)
  x%>% as.data.frame() %>% filter(padj < pval) %>% arrange(padj)
})

old_genes <- table(unlist(lapply(old_PO_0.05, function(x) x$gene)))
old_genes[which(old_genes > 1)]

new_PO_0.05 <- lapply(names(rem_seg_res$res_PO), function(x){
  rem_seg_res$res_PO[[paste(x)]]$gene = rownames(rem_seg_res$res_PO[[paste(x)]])
  rem_seg_res$res_PO[[paste(x)]]$rix = x
  rem_seg_res$res_PO[[paste(x)]] %>% as.data.frame() %>% filter(padj < pval) %>% arrange(padj)
})
new_PO_df = do.call("rbind", new_PO_0.05)
new_PO_df %>% arrange(padj, pval_raw, gene, rix) -> new_PO_df

old_PO_0.05 <- lapply(names(dds_string_compare$res_PO), function(x){
  dds_string_compare$res_PO[[paste(x)]]$gene = rownames(dds_string_compare$res_PO[[paste(x)]])
  dds_string_compare$res_PO[[paste(x)]]$rix = x
  dds_string_compare$res_PO[[paste(x)]] %>% as.data.frame() %>% filter(padj < pval) %>% arrange(padj)
})
old_PO_df = do.call("rbind", old_PO_0.05)
old_PO_df %>% arrange(padj, pval_raw, gene, rix) -> old_PO_df

new_genes <- table(unlist(lapply(new_PO_0.05, function(x) x$gene)))
plot_genes = new_genes[-grep("Rik|[.]|Gm", names(new_genes))]
plot_genes = new_genes[which(new_genes > 2)]
new_PO_df %>% filter(gene %in% names(plot_genes)) -> which_rix



## start plotting
for(g in names(plot_genes)){
  pd = lapply(names(rem_seg_dds), function(x){
    if(g %in% rownames(rem_seg_dds[[paste(x)]])){
      pd = plotCounts(rem_seg_dds[[paste(x)]], gene=g, intgroup = "PO",returnData = T)
      pd$Pup.ID = rownames(pd)
      pd$RIX = x
      pd
    }
  })
  pd = do.call("rbind", pd)
  find_gene = annot_genes[which(annot_genes$Gene.Name == g),]
  which_rix = (new_PO_df %>% filter(gene == g))$rix
  
  
  masterKmers = readRDS(paste0(dir,"/variant_data/kmers_to_run/masterSnps_chr",find_gene$Chr,"_26aug2019.rds"))
  masterKmers = do.call("cbind", masterKmers[[1]])
  maxmarginal = readRDS(file.path(dir, "/mini/interp_rqtl_allChr_maxmarginal_19nov2019.rds"))
  maxmarginal = maxmarginal[[paste(find_gene$Chr)]]
  if(find_gene$Gene.Name %in% masterKmers$seq.Gene){
    minf = data.frame(masterKmers[which(masterKmers$seq.Gene == find_gene$Gene.Name),])
    dists = apply(minf, 1, function(x){
      min(abs(find_gene$End - as.numeric(x["seq.Position"])), abs(as.numeric(x["seq.Position"]) - find_gene$Start))
    })
    minf = minf[order(dists),]
  } else {
    min_1=which.min(abs(find_gene$Start - masterKmers$seq.Position))
    min_2=which.min(abs(find_gene$End - masterKmers$seq.Position))
    if(min_1 != min_2){
      if(min(abs(find_gene$Start - masterKmers$seq.Position)) < min(abs(find_gene$End - masterKmers$seq.Position))){
        minf = rbind(masterKmers[min_1,], masterKmers[min_2,])
      } else {
        minf = rbind(min_2, min_1)
      }
    } else {
      minf = masterKmers[min_1,]
    }
  }
  
  kmer_compare = minf[,grep("sdp", colnames(minf))]
  #annot_cc = apply(pd, 1, function(x){
  c_this = list()
  for(i in 1:nrow(pd)){
    x = pd[i,]
    pup=gsub("[^0-9]", "", x["Pup.ID"])
    if(pup %in% names(phased_par_CC_haplotypes)){
      if(samples_use$PO_cat[which(as.character(samples_use$Pup.ID) == pup)] == "a"){
        cc_1 = phased_par_CC_haplotypes[[pup]]$founder_by_parent[[1]]
        cc_2 = phased_par_CC_haplotypes[[pup]]$founder_by_parent[[2]]
        po = "a" 
      } else {
        cc_1 = phased_par_CC_haplotypes[[pup]]$founder_by_parent[[2]]
        cc_2 = phased_par_CC_haplotypes[[pup]]$founder_by_parent[[1]]
        po = "b" 
      }
      c_this_tmp = data.frame(cc_1_found = cc_1$found[which(cc_1$chr == find_gene$Chr & 
                                                              cc_1$start < find_gene$Start & cc_1$end > find_gene$End)],
                              cc_1 = unique(cc_1$cc),
                              cc_2_found = cc_2$found[which(cc_2$chr == find_gene$Chr & 
                                                              cc_2$start < find_gene$Start & cc_2$end > find_gene$End)],
                              cc_2 = unique(cc_2$cc), po=po)
      convert = data.frame(let = LETTERS[1:8], num = 1:8)
      found_num_1 = convert$num[which(convert$let == c_this_tmp$cc_1_found)]
      found_num_2 = convert$num[which(convert$let == c_this_tmp$cc_2_found)]
      c_this[[i]] = cbind(Pup.ID = paste0("Pup.ID_",pup),
                          c_this_tmp, 
                          kmer_compare[,c(paste0("sdp.",found_num_1),paste0("sdp.",found_num_2))])
      colnames(c_this[[i]])[7:8] = c("sdp_1", "sdp_2")
    }
  }
  c_this <- do.call("rbind", c_this)
  c_this$het = ifelse(c_this$sdp_1 == c_this$sdp_2, F, T)
  pd %>% left_join(c_this, by="Pup.ID") -> pd
  #pd$xlab = paste0(pd$cc_1, "(",pd$cc_1_found,"/",pd$sdp_1,") x ", pd$cc_2, "(",pd$cc_2_found,"/",pd$sdp_2,")")
  pd$xlab = paste0(pd$cc_1, "(",pd$cc_1_found,"/",pd$sdp_1,")")
  pd$which_rix = ifelse(pd$RIX %in% which_rix, T, F)
  
  for(fix in which(is.na(pd$cc_1_found))){
    tmp = pd[fix,]
    replace = which(pd$PO == tmp$PO & pd$RIX == tmp$RIX & !is.na(pd$cc_1_found))[1]
    pd[fix, 5:ncol(pd)] = pd[replace, 5:ncol(pd)]
  }
  
  col_use = add.alpha(c("lightblue"), 0.4)
  pdf(paste0(dir, "/figures_and_updates/images_meetings/multRIX_PO_counts_",g,"_25nov2019.pdf"),
      width=9, height=7)
  par(oma = c(3, 1, 1, 1), mfrow=c(3,3))
  for(r in c(1:4,6:10)){
    if(r %in% unique(pd$RIX)){
      pd_use = pd %>% filter(RIX == r) %>% arrange(PO)
      pd_use$xlab = factor(pd_use$xlab)
      xlab = paste(pd_use$cc_1[which(pd_use$po == "a")[1]], "x",  pd_use$cc_2[which(pd_use$po == "a")[1]])
      boxplot(pd_use$count ~ pd_use$xlab , 
              outline=F, ylim=c(min(pd_use$count), max(pd_use$count)), 
              main=paste(g, "RIX", r), xlab = xlab,
              col=ifelse(pd_use$which_rix, col_use, "white"))  
      points(x=jitter(as.numeric(pd_use$xlab)), y=pd_use$count, 
             col=ifelse(pd_use$sdp_1 == 1, "seagreen1", "slategray4"), pch=20)  
    } else {
      plot(0,0)
    }
  }
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend("bottom", c("Signif RIX", "Alt allele"), xpd = TRUE, horiz = TRUE, 
         inset = c(0,0), bty = "n", pch = c(15,20), 
         col = c(col_use, "seagreen1"),cex=1.5, text.col = "black")
  dev.off()
}
####### plot counts from some significant genes 

new_PO_df = new_PO_df %>% #filter(baseMean > 10) %>% #filter(gene %in% names(plot_genes)) %>%
  group_by(gene) %>% arrange(pval_raw, -baseMean)
unique((new_PO_df %>%    #gene %in% names(plot_genes), 
  arrange(pval_raw))$gene) -> plot_genes
#plot_genes = new_PO_0.05$gene[-grep("Gm|mt|rik", new_PO_0.05$gene)]
pdf("C:/Users/Kathie/Dropbox (ValdarLab)/figures_and_updates/images_meetings/top_new_PO_counts_remSegv2_19nov2019.pdf")
par(mfrow=c(2,2))
for(i in 1:nrow(new_PO_df[1:48,])){
  gene = new_PO_df$gene[i]
  rix= new_PO_df$rix[i]   #c(4,4,6,8,8,9)
  pd <- plotCounts(rem_seg_dds[[paste(rix)]], gene=gene, intgroup = "PO",returnData = T)
  boxplot(count ~ as.factor(PO), data=pd, outline=F, ylim=c(min(pd$count), max(pd$count)),
       main=paste("Normalized counts in",gene, rix), xlab = "PO")       
  points(x=jitter(as.numeric(as.factor(pd$PO))), y=pd$count)  
  #anova(lm(data=pd, count ~ PO))
}
dev.off()

######  CORRELATIONS BETWEEN STRINGTIE AND KMERS  ######

all_cs_genes = table(unlist(lapply(dds_cs_compare$res_PO, function(x) rownames(x))))
cs_genes_consistent = names(all_cs_genes)[which(all_cs_genes == max(all_cs_genes))]
all_string_genes = table(unlist(lapply(dds_string_compare$res_PO, function(x) rownames(x))))
string_genes_consistent = names(all_string_genes)[which(all_string_genes == max(all_string_genes))]

genes_consistent = intersect(cs_genes_consistent, string_genes_consistent)
## 18634 unique genes with kmers in every RIX
## 27782 unique genes in StringTie in every RIX
## 14369 unique overlaps

str_cs_compare = list()
for(r in names(dds_string_compare$res_PO)){
  tmp_string = dds_string_compare$res_PO[[r]][match(genes_consistent, rownames(dds_string_compare$res_PO[[r]])),]
  tmp_cs = dds_cs_compare$res_PO[[r]][match(genes_consistent, rownames(dds_cs_compare$res_PO[[r]])),]
  colnames(tmp_string) = paste0(colnames(tmp_string),"_str")
  colnames(tmp_cs) = paste0(colnames(tmp_cs),"_cs")
  tmp_string = data.frame(tmp_string)
  tmp_string$gene = tmp_cs$gene = genes_consistent
  
  str_cs_compare[[r]] = left_join(as.data.frame(tmp_string), as.data.frame(tmp_cs), by="gene")
}

sapply(1:9, function(r) 
  cor(x=str_cs_compare[[r]]$pval_raw_str, y=str_cs_compare[[r]]$pval_raw_cs, use="pairwise.complete.obs"))
for(r in names(str_cs_compare)){
  plot(x=str_cs_compare[[r]]$pvalue_str, y=str_cs_compare[[r]]$pvalue_cs,
       main = paste("P-value comparison in",r), xlab = "StringTie", ylab = "K-mers")
}

##############  DIET ANALYSES  ####################
pval = 0.01
#dds_cs_compare$res_diet$rix6$Standard[order(dds_cs_compare$res_diet$rix6$Standard$padj),]
#dds_string_compare$res_diet$rix6$Standard[order(dds_string_compare$res_diet$rix6$Standard$padj),]
lapply(keepSeg_res_diet$res_diet, function(x) if("Standard" %in% names(x)) print(hist(x$Standard$pval_raw)))
### VDD-STD in RIX6 only; includes X chr data
### dds_string
print=F
sig_VDD_STD_list = list()
for(r in names(keepSeg_res_diet$res_diet)){
  if(!is.null(keepSeg_res_diet$res_diet[[r]]$Standard)){
    VDD_STD_string = as.data.frame(keepSeg_res_diet$res_diet[[r]]$Standard) 
    VDD_STD_string$Gene = rownames(keepSeg_res_diet$res_diet[[r]]$Standard)
    sig_VDD_STD_string = VDD_STD_string %>% filter(padj_deseq2 < pval) %>% arrange(padj)
    sig_VDD_STD_list[[r]] = sig_VDD_STD_string
    
    if(print){
      write.csv(sig_VDD_STD_string, paste0(dir, "/de_results/signif_datasets/sig_VDD_STD_keepSeg_",r,"_string_22nov2019.csv"), quote = F)
    }
    plot_genes = sig_VDD_STD_string$Gene[which(sig_VDD_STD_string$padj_deseq2 == min(sig_VDD_STD_string$padj_deseq2, na.rm=T))]
    if(length(plot_genes) > 48) plot_genes = plot_genes[1:48]
    plot_genes = sig_VDD_STD_string$Gene[1:min(nrow(sig_VDD_STD_string), 48)]
    print(dim(sig_VDD_STD_string)[1])
  
    if(print){
      #par(mfrow=c(2,2))
      pdf(paste0("C:/Users/Kathie/Dropbox (ValdarLab)/figures_and_updates/images_meetings/sig_vdd_std_",r,"_keepSeg_countPlots_22nov2019.pdf"))
      #p = list()
      par(mfrow=c(2,2))
      for(i in 1:length(plot_genes)){
        gene = plot_genes[i]
        rix= r
        pd <- plotCounts(keepSeg_dds_diet[[paste(r)]], gene=gene, intgroup = "Diet",returnData = T)
        boxplot(count ~ Diet, data=pd, outline=F, ylim=c(min(pd$count), max(pd$count)),
                main=paste("Normalized counts in",gene, rix), xlab = "Diet")      
        points(x=jitter(as.numeric(pd$Diet)), y=pd$count)  
      
        #if(i%%4 == 0){
        #  quad = grid.arrange(p[[i-3]], p[[i-2]], p[[i-1]], p[[i]], ncol=2, nrow=2)
        #  print(quad)
        #}
      }
      
      dev.off()
    }
    
  }
}

all_sig_genes = lapply(names(sig_VDD_STD_list),function(x) data.frame(gene = sig_VDD_STD_list[[x]]$Gene,
                                                                      rix  = x ))
all_sig_genes = do.call("rbind", all_sig_genes)
tab_genes = table(unlist(lapply(sig_VDD_STD_list, function(x) x$Gene)))
tab_genes = tab_genes[which(tab_genes > 1)]
tab_genes = tab_genes[order(-tab_genes)]

all_sig_genes %>% filter(gene %in% names(tab_genes)) %>%
  group_by(gene) %>%
  mutate(n = n(), rixes = paste(rix, collapse=",")) %>%
  dplyr::select(-"rix") %>%
  arrange(-n, gene) %>%
  distinct() -> mult_sig_genes
write.csv(mult_sig_genes, file.path(dir, "de_results/mult_VDD_STD_sigif_genes_keepSeg_14nov2019.csv"))
#sig_VDD_STD_rix6_string = read.csv(file.path(dir, "de_results/sig_VDD_STD_rix6_string_24oct2019.csv"))
#sig_VDD_STD_rix6_string = sig_VDD_STD_rix6_string[,-1]
mult_sig_genes %>% filter(n > 3)
plot_genes=mult_sig_genes$gene[1:11]
pdf(paste0("C:/Users/Kathie/Dropbox (ValdarLab)/figures_and_updates/images_meetings/sig_keepSeg_vdd_std_multRIX_countPlots_22nov2019.pdf"))
par(mfrow=c(3,3))
for(i in 1:length(plot_genes)){
  gene = plot_genes[i]
  rixes = mult_sig_genes$rixes[which(mult_sig_genes == gene)]
  rixes = unlist(strsplit(rixes, ","))
  #rix=paste0("rix",c(4,4,6,8,8,9))[i]    #c(4,4,6,8,8,9)
  for(rix in names(keepSeg_dds_diet)){
    if(gene %in% rownames(keepSeg_dds_diet[[rix]])){
      col = ifelse(rix %in% rixes, "blue","black")
      pd <- plotCounts(keepSeg_dds_diet[[rix]], gene=gene, intgroup = "Diet",returnData = T)
      boxplot(count ~ Diet, data=pd, outline=F, ylim=c(min(pd$count), max(pd$count)),
              main=paste(gene, rix), xlab = "Diet")      
      points(x=jitter(as.numeric(pd$Diet)), y=pd$count, col=col)  
    } else {
      plot(0,0)
    }
    
  }
  
  #anova(lm(data=pd, count ~ PO))
}
dev.off()


res = results(dds_cs_sva, contrast=c("Diet", "Vitamin.D.Deficient", "Standard"))
res = res[-which(is.na(res$pvalue)),]
fdr <- fdrtool(res$stat, statistic="normal", plot=F)
res$pval_raw = res$pvalue
res$padj_deseq2 = res$padj
res[,"pvalue"] = fdr$pval
res[,"padj"]  <- p.adjust(res$pvalue, method = "BH")
VDD_STD_rix6_cs = as.data.frame(res)
VDD_STD_rix6_cs$Gene = rownames(res)
sig_VDD_STD_rix6_cs = VDD_STD_rix6_cs %>% filter(padj < 0.05, baseMean > 10) %>% arrange(padj)
#write.csv(sig_VDD_STD_rix6_cs, file.path(dir, "de_results/sig_VDD_STD_rix6_cs_24oct2019.csv"), quote = F)



for(r in names(remSeg_res_diet$res_diet)){
  VDD_STD_rix6_string = as.data.frame(remSeg_res_diet$res_diet$rix6$Standard) 
  VDD_STD_rix6_string$Gene = rownames(remSeg_res_diet$res_diet$rix6$Standard)
  sig_VDD_STD_rix6_string = VDD_STD_rix6_string%>% filter(padj < pval, baseMean > 10) %>% arrange(padj)
  

  
}
VDD_STD_rix6_cs = as.data.frame(remSeg_res_diet$res_diet$rix6$Standard) 
VDD_STD_rix6_cs$Gene = rownames(remSeg_res_diet$res_diet$rix6$Standard)
sig_VDD_STD_rix6_cs = sig_VDD_STD_rix6_cs%>% filter(padj < pval, baseMean > 10) %>% arrange(padj)

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


#############  comparing counts across methods  #################
use_cs <- as_tibble(counts(dds_lst$rix6, normalized=T)) %>% mutate(gene = rownames(counts(dds_cs$rix6))) %>%
  gather(pup, counts, -gene) %>%
  mutate(pup=gsub("^X", "", pup), method = "cs")
use_mm10 <- as_tibble(counts(dds_mm10$rix6, normalized=T)) %>% mutate(gene = rownames(counts(dds_mm10$rix6))) %>%
  gather(pup, counts, -gene) %>%
  mutate(pup=paste0("Pup.ID_", pup), method = "mm10")
use_string <- as_tibble(counts(dds_string$rix6, normalized=T)) %>% mutate(gene = rownames(counts(dds_string$rix6)))  %>% gather(pup, counts, -gene) %>%
  mutate(method = "string")
use_gecco <- gecco_uniq %>% rename("sample_id" = "pup", "TReC" = "counts") %>% dplyr::select(-ensembl_id) %>%
  mutate(method = "gecco") %>% dplyr::select(gene, pup, counts, method) %>%
  filter(!is.na(gene))
compare <- rbind(use_cs, use_mm10, use_string, use_gecco) 

factors <- data.frame(method=c("cs", "mm10", "string", "gecco"),
                      means=c(mean(use_cs$counts),mean(use_mm10$counts),mean(use_string$counts),mean(use_gecco$counts)), 
                      stringsAsFactors = F)
factors$fac <- factors$means/min(factors$means)

n=length(use_gecco$gene)
use1 <- intersect(use_gecco$gene[order(-use_gecco$counts)[1:n]], use_cs$gene[order(-use_cs$counts)])
use2 <- intersect(use_string$gene[order(-use_string$counts)], use1)

compare %>% filter(gene %in% use2) -> plotMat

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

corr <- means %>% spread(method, mean)
cor.test(corr$cs, corr$string)
cor.test(corr$cs, corr$gecco)
cor.test(corr$gecco, corr$string)
ggplot(corr, aes(x=gecco, y=string)) + 
  geom_point()
cor(corr[,-1])

#############  comparing counts across methods  #################
## looking at specific genes in RIX 6 VDD-STD comparison

genes = sig_VDD_STD_rix6_string$Gene

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


################  some old QC code that might be handy to run again  ###############

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

