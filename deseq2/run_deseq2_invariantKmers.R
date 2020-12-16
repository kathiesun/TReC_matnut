setwd("C:/Users/Kathie/rna_seq/kmerSearch")
setwd("~/rna_seq/kmerSearch)/kmerSearch/")

library(tidyverse)
library(DESeq2)
library(sva)

dir <- file.path("C:/Users/Kathie/Dropbox (ValdarLab)/variant_data/")
problemPups <- c(1404, 1716, 1371, 569)  #, 1911, 1951, 1015)
invar_info <- read.csv(file.path(dir,"kmers_based_on_CC/kmers_to_run/GeneProbesInvariant.csv"))

samples <- read.csv(file.path(dir, "2015-10_expression_pups.csv"))
rownames(samples) = samples$Pup.ID
samples$Diet <- gsub(" $", "", samples$Diet)
samples$Diet <- factor(samples$Diet, levels=c("Standard", "Low Protein","Methyl Enriched","Vitamin D Deficient"))
samples$PO <- ifelse(factor(gsub("[0-9]","", samples$RIX)) == "a", 0.5, -0.5)
samples$PO_cat <- factor(gsub("[0-9]","",samples$RIX))

dds_lst <- list()
res <- list()
rix_data = readRDS(file.path(dir, "kmers_based_on_CC/kmer_counts/invar_counts_by_rix_5aug2019.rds"))
rixes = paste0("RIX",c(1:4,6:10))
sva=F
for(rix in rixes[1]){  
  merge_list = rix_data[[rix]]
  comb <- do.call("cbind",lapply(merge_list, function(x) x["sum"]))
  cts <- apply(comb, c(1,2), as.numeric)
  colnames(cts) <- names(merge_list)
  rownames(cts) <- merge_list[[1]]$Gene.Name
  cts = as.matrix(cts)
  
  coldata <- samples %>% 
    dplyr::select(one_of("Pup.ID","RRIX","Diet","PO","PO_cat","Breeding.Batch","Behavior.Batch")) %>%
    mutate(Breeding.Batch  = as.factor(Breeding.Batch), Behavior.Batch = as.factor(Behavior.Batch),
           Pup.ID = paste0("Pup.ID_", Pup.ID)) %>%
    dplyr::filter(Pup.ID %in% colnames(cts)) %>%
    arrange(Pup.ID, colnames(cts)) 
  rownames(coldata) = coldata$Pup.ID
  all(rownames(coldata) == colnames(cts))
  
  dds <- DESeqDataSetFromMatrix(countData = cts,
                                   colData = coldata,
                                   design = ~ Diet + PO_cat)
  dds = DESeq(dds)
  if(sva){
    dat = counts(dds, normalized=T)
    dat = dat[rowMeans(dat) > 1,]
    mod  <- model.matrix(~ Diet + PO, colData(dds))
    mod0 <- model.matrix(~   1, colData(dds))
    svseq <- svaseq(dat, mod, mod0, n.sv = 1)  
    ddssva = dds
    ddssva$SV1 <- svseq$sv[,1]
    
    design(ddssva) <- ~ SV1 + Diet + PO_cat 
    ddssva = DESeq(ddssva)
    ddssva = nbinomWaldTest(ddssva)
    dds = ddssva
  }
  
  dds_lst[[tolower(rix)]] = dds

  print(paste("finished deseq on ",rix))
}

saveRDS(dds_lst, file.path(dataSource, "../de_results/dds_lst_cs_nointeraction-16oct2019.rds"))
#saveRDS(rix_data, file.path(dataSource, "kmer_counts/invar_counts_by_rix_5aug2019.rds"))

##################################################################
## using RIX1
dds_interaction <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ Diet + PO + Diet:PO)
dds_interaction = DESeq(dds_interaction)
dds_nointeraction <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ Diet + PO)
dds_nointeraction = DESeq(dds_nointeraction)


## still using RIX1, sva run above
dds_sva_interaction = dds
dds_sva_nointeraction = dds


res_sva_inter = results(dds_sva_interaction, name="PO")
res_sva_nointer = results(dds_sva_nointeraction, name="PO")
res_inter = results(dds_interaction, name="PO")
res_no_inter = results(dds_nointeraction, name="PO")

pdf(file.path(dir,"../figures_and_updates/kmer_pval_wInterSVA_compare.pdf"))
hist(res_inter$pvalue, main = "Raw p-values with interaction and no SV")
hist(res_no_inter$pvalue, main = "Raw p-values with no interaction and no SV")
hist(res_sva_inter$pvalue, main = "Raw p-values with interaction and SV")
hist(res_sva_nointer$pvalue, main = "Raw p-values with no interaction and SV")
dev.off()


#################  IDK WHAT THIS IS  ###################

sorted_invar_rix1 <- sort_counts(counts = invar_1, snp_info = invar_info)

ratios %>% group_by(Gene.Name) %>% tally() %>% filter(n > 1) %>% dplyr::select("Gene.Name") -> genesKeep
ratios %>%
  filter(Gene.Name %in% t(genesKeep)) %>%
  arrange(Gene.Name, desc(counts)) %>%
  group_by(Gene.Name) %>%
  mutate(max = max(counts), min = min(counts),
         range = max-min) -> ratios

print(ratios %>% arrange(desc(range)) %>%
  filter(row_number() %in% c(1, n())) %>%
  dplyr::select(-one_of("Pup","Gene.ID")), n=50)

metaseq <- read.csv(paste0(dataSource, "/../../../matnut_master/metaseq.csv"))
head(metaseq)
metaseq %>% group_by(pupID) %>% summarize(sum = sum(numreads)) %>% arrange -> sums


######## FUNCTION TO SORT INVARIANT COUNTS ###############

sort_counts_invar <- function(counts, snp_info){
  merged <- counts %>% left_join(snp_info, by="Kmer") %>%
    filter(!is.na(Gene.Name)) 
  merge_list <- list()
  for(i in 1:length(unique(merged$Pup))){
    merge_tmp <- filter(merged, Pup == unique(merged$Pup)[i])
    merge_tmp %>% group_by(Pup, Chromosome, Gene.Name) %>%
      summarize(mean = mean(counts), median = median(counts),
                n = n(), sum=sum(counts), pos = mean(Position)) %>%
      arrange(desc(median)) -> pd
    #pd_2 <- pd %>% filter(Chromosome == 2)
    
    #ggplot(data=pd_2, aes(x=pos, y=log(median+1), alpha = sum)) + 
    #  geom_point(position="jitter") + 
    #  facet_wrap(~ Chromosome, ncol=5, scales = "free") + 
    #  theme_bw()
    merge_list[[paste0("Pup.ID_",unique(merged$Pup)[i])]] <- pd
  }
  return(merge_list)
}

invar_files <- list.files(file.path(dataSource, "kmers_based_on_CC/kmer_counts"), "invar_counts_RIX", full.names = T)
#invar_counts <- lapply(list.files(file.path(dataSource, "jul2019_allRNA_kmer_counts/"), "invar_counts_", full.names = T), read.csv)
