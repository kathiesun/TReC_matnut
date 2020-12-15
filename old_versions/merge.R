#### merge

### on cluster
setwd("~/matnut/src")
library(tidyverse)
library(Biostrings)
library(rtracklayer)
dataSource <- "/nas/depts/006/valdar-lab/users/sunk/"


chrlist <- as.character(c(1:19))
hetPos <- list()

for(i in 1:length(chrlist)){
  chr <- chrlist[i]
  crosses <- readRDS(paste0(dataSource,"CCparental_crosses/chr", chr, ".rds"))
  allcrosses <- do.call("rbind",lapply(crosses, function(x) unique(x %>% 
                                                                     filter(prob>0.8) %>%   
                                                                     select(-one_of("variant_id", "transcript_name","consequence1","consequence2")))))
  allcrosses %>% 
    select(-one_of("gene_name")) %>%
    count(pos) %>% 
    filter(n == length(unique(allcrosses$strain1))) %>% 
    select(pos) %>% distinct() -> completeSetPos
  allcrosses %>% filter(pos %in% completeSetPos$pos) %>% 
    mutate(equal = ifelse(allele1 == allele2, 0, 1)) %>%
    arrange(pos) -> tmp
  sums <- aggregate(equal ~ pos, tmp, sum)
  sums %>% filter(equal > 0) %>%
    select(pos) %>% distinct() -> completeHetPos
  hetPos[[i]] <- allcrosses %>% filter(pos %in% completeHetPos$pos) 
}

cHetPos <- do.call("rbind", hetPos)
cHetPos %>% rename("gene_name"="gene_id") -> cHetPos

genes <- import(paste0(dataSource, "mm10/Mus_musculus.GRCm38.92.gtf"))
getTrans <- as.data.frame(genes[,c("gene_id", "gene_name")])
getTrans <- getTrans %>% select(one_of("gene_id", "gene_name")) %>% distinct()
cHetPos %>% left_join(getTrans, by="gene_id") -> annotHetPos


keep <- which(getTrans$gene_id %in% genelist)
getTrans <- getTrans[keep,]
getTrans %>% 
  group_by(gene_name) %>%
  arrange(desc(width)) %>%
  select(one_of("strand", "gene_id", "gene_name")) %>% 
  dplyr::slice(1) -> uniqGenes
