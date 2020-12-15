## This code can be modified to get diplotypes, inbred lines themselves rather than crosses, and can be easily run in parallel if you need it to (e.g.., to rapidly generate lots of random crosses)
## Let me know if you want to do those things.
setwd("~/matnut/src")

####### ON CLUSTER #######
source("./genomerep/variantdb2/dump_parser.R")
root = "/nas/depts/006/valdar-lab/PUBLIC/2018_1_isvdb/"

##dbdir = fp(root, "exon_1410") ##MGP 1410, Just the exons-- manipulating the exon data is much faster than whole genome.
        ##38.75
##dbdir = fp(root, "full_1504") ##MGP 1504, Whole genome, ensembl build 38.83
dbdir = fp(root, "exon_1504")

db    = db_builder$get.db.lite(dbdir)

allchrs = db$genotype$getChrs()

#strain1s = c("CC001", "CC002", "CC003")
#strain2s = c("CC004", "CC005", "CC006")

strain1s <- c("CC001", "CC041", "CC017", "CC023","CC026","CC014", "CC035", "CC032","CC005")
strain2s <- c("CC011", "CC051", "CC004", "CC047","CC006","CC003", "CC062", "CC042","CC040")

strain1s <- c("A_J", "129S1_SvImJ", "NOD_ShiLtJ","NZO_HlLtJ","CAST_EiJ", "PWK_PhJ", "WSB_EiJ")
strain2s <- c("C57BL6J", "C57BL6J", "C57BL6J", "C57BL6J","C57BL6J","C57BL6J", "C57BL6J")


chrlist = as.character(paste(c(1:19,"X")))
dataSource <- "/nas/depts/006/valdar-lab/users/sunk/"

for(j in 1:length(chrlist)){
  chr <- chrlist[j]
  crosses_dip = list()
  crosses <- list()
  for(i in 1:length(strain1s))
  {
    print(paste0("cross", i))
    strain1 = strain1s[i]
    strain2 = strain2s[i]
    
    dip = db$read(strain1 = strain1, strain2 = strain2, type = "diplotype", phased = T, chr = chr)
    cross = db$read(strain1 = strain1, strain2 = strain2, type = "genotype", phased = T, chr = chr)
    crosses_dip[[i]] = dip
    crosses[[i]] = cross
  }
  
  saveRDS(crosses, paste0(dataSource,"/CCparental_crosses/chr",chr,".rds"))
  saveRDS(crosses_dip, paste0(dataSource,"/CCparental_crosses/dip_chr",chr,".rds"))
  
}


########### PATERNITY CHECK ############
setwd("~/matnut/src")
library(tidyverse)
library(Biostrings)
library(rtracklayer)
dataSource <- "/nas/depts/006/valdar-lab/users/sunk/"


chrlist <- as.character(c(1:19))
patQueries <- list()
genelist <- c("ENSMUSG00000049339","ENSMUSG00000026932","ENSMUSG00000001416","ENSMUSG00000028677",
           "ENSMUSG00000029713","ENSMUSG00000030122","ENSMUSG00000030806",
           "ENSMUSG00000070000","ENSMUSG00000032936","ENSMUSG00000003072",
           "ENSMUSG00000019302","ENSMUSG00000021066","ENSMUSG00000034675",
           "ENSMUSG00000022180","ENSMUSG00000048385","ENSMUSG00000022892",
           "ENSMUSG00000033597","ENSMUSG00000051375","ENSMUSG00000079415",
           "ENSMUSG00000026932")

for(i in 1:length(chrlist)){
  chr <- chrlist[i]
  gene <- genelist[i]
  crosses <- readRDS(paste0(dataSource,"CCparental_crosses/chr", chr, ".rds"))
  allcrosses <- do.call("rbind",lapply(crosses, function(x) unique(x %>% 
                                                                     filter(prob>0.8, nchar(allele1) == 1, gene_name %in% genelist) %>%   
                                                                     select(-one_of("transcript_name","consequence1","consequence2")))))
  allcrosses %>% count(pos) %>% filter(n == 7) %>% select(pos) -> keepPos
  patQueries[[i]] <- allcrosses %>% filter(pos %in% keepPos$pos) 
}

queries <- do.call("rbind", patQueries)
colnames(queries)[6] <- "gene_id"

genes <- import(paste0(dataSource, "mm10/Mus_musculus.GRCm38.92.gtf"))
getTrans <- as.data.frame(genes[,c("transcript_id", "gene_id", "gene_name")])
keep <- which(getTrans$gene_id %in% genelist)
getTrans <- getTrans[keep,]
getTrans %>% 
  group_by(gene_name) %>%
  arrange(desc(width)) %>%
  select(one_of("strand", "gene_id", "gene_name")) %>% 
  dplyr::slice(1) -> uniqGenes
queries %>% left_join(uniqGenes, by="gene_id") -> annotqueries

## remove hets that are too close to each other
klen=35
distBtwHets <- unlist(sapply(1:nrow(annotqueries) - 1, function(x) annotqueries$pos[x+1] - annotqueries$pos[x]))
remove <- which(distBtwHets < ((klen+1)/2))
annotqueries[-remove, ] %>% 
  mutate(start = pos - ((klen-1)/2), end = pos + ((klen-1)/2)) -> kmerSet

 

#write.csv(kmerSet,"../../CCparental_crosses/patTestKmers.csv")
patTestKmers <- read.csv(paste0(dataSource,"/CCparental_crosses/patTestKmers.csv"))
for(i in 1:length(chrlist)){
  chr <- chrlist[i]
  gene <- genelist[i]
  crosses <- readRDS(paste0(dataSource,"CCparental_crosses/chr", chr, ".rds"))
  allcrosses <- do.call("rbind",lapply(crosses, function(x) unique(x %>% 
                                                                     filter(pos %in% patTestKmers$pos) %>%   
                                                                     select(-one_of("transcript_name","consequence1","consequence2")))))
  patQueries[[i]] <- allcrosses
}

queries <- do.call("rbind", patQueries)
colnames(queries)[6] <- "gene_id"


############
setwd("~/matnut/src")
library(refGenome)
library(tidyverse)
library(Biostrings)

dataSource <- file.path("C:/DB Mount","Dropbox\ (ValdarLab)","outputs", "matnut_outputs")
crosses <- readRDS(paste0(dataSource,"/allIsvdbCrosses.rds"))
crosses_dip <- readRDS(paste0(dataSource,"/allIsvdbCrosses_dip.rds"))
xgenes <- readRDS(paste0(dataSource,"/xgenes.rds"))

annot_x <- readRDS(paste0(dataSource,"/annot_x.rds"))
hetsInGenes <- readRDS(paste0(dataSource,"/hetsInGenes.rds"))

annot_x <-  annot_x[,c("gene_id", "gene_name", "gene_biotype", "start", "end")]

## annotated genes with heterozygous eSNPs
hetsInGenes[[1]] %>% dplyr::rename("gene_id"="gene_name") %>%
  left_join(annot_x, by="gene_id") %>%
  mutate(len = end-start) -> annot_het_x

## how many het positions on x chr in each rix?
hetsInGenes <- list()
uniquePosFounders <- list()
uniqueHets <- list()
numhet <- c()

for(i in 1:length(crosses)){
  keep <- intersect(which(crosses[[i]]$allele1 != "None"),
                    which(crosses[[i]]$allele2 != "None"))
  notna <- crosses[[i]][keep,]
  equal <- ifelse(notna$allele1 == notna$allele2, 0, 1)
  het <- as.data.frame(notna[which(equal == 1),])
  
  het %>% filter(prob>0.95) %>% distinct(gene_name) -> confGene
  het %>% filter(prob>0.95) %>% distinct(pos) -> confPos
  
  het %>% as.tibble() %>%
    dplyr::select(-consequence1, -consequence2, -transcript_name) %>%
    filter(gene_name %in% as.vector(t(confGene))) %>% 
    #filter(pos %in% as.vector(t(confPos))) %>%
    group_by(pos) %>%               #
    distinct() %>%          
    top_n(n=1, wt=prob) %>%  
    dplyr:::slice(1) -> uniqueHets[[i]]
  
  uniqueHets[[i]] %>% 
    group_by(gene_name) %>%
    summarize(n_unique = n_distinct(pos)) %>%
    arrange(desc(n_unique)) -> hetsInGenes[[i]]
  hetsInGenes[[i]] %>% top_n(10, n_unique)
    
  crosses_dip[[i]] %>% 
    filter(gene_name %in% as.vector(t(confGene))) %>%
  #crosses_dip[[i]][keep_dip, ]  %>%
    group_by(pos) %>%                 #
    top_n(1, prob) %>%                #
    group_by(founder1, founder2) %>% 
    summarize(n_genes = n_distinct(gene_name),
              n_pos = n_distinct(pos)) -> uniquePosFounders[[i]]
  
  numhet <- c(numhet, length(unique(hetsInGenes[[i]]$gene_name)))

}

uniqueHets[[1]] %>% 
  filter(gene_name == annot_het_x$gene_id[which(annot_het_x$gene_name == "Gdi1")]) %>% 
  filter(prob > 0.5)

topGenes <- c()
for(i in 1:9){
  hetsInGenes[[i]] %>% filter(n_unique>50) -> tops
  topGenes <- c(topGenes, as.character(tops$gene_name))
}

data.frame(table(topGenes)) %>% 
  arrange(desc(Freq)) %>% 
  filter(Freq > 3) -> topGenes_hi
  

## just quantifies how many eSNPs in each gene, precursor to annot_het_x
for(i in 1:9){
  hetsInGenes[[i]] %>% 
    filter(gene_name %in% topGenes_hi$topGenes) %>%
    print()
}

uniqueHets[[1]] %>% filter(gene_name == hetsInGenes[[1]]$gene_name[1]) %>%
  arrange(prob) %>%
  top_n(10, prob)


######### OLD CODE ########

Xce <- matrix(c(4324915, 167000000,
                68073962,	102881732,
                4324915,	102881732,
                91749873,	167000000,
                97885202,	132952034,
                12842817,	109780652,
                4324915,	140042023,
                47905207,	167000000,
                4324915,	141792095,
                34369992,	167000000,
                78336760,	160607985,
                67613211,	127901889,
                17105975,	130598641,
                98769446,	124825747,
                96339885,	148751441,
                80070551,	123299962,
                92265002,	131104444,
                4324915,	136529143), ncol=2, nrow=18, byrow=T)

for(i in 1:length(crosses)){
  #pos1 <- min(Xce[2*i - 1,1], Xce[2*i,1])
  #pos2 <- max(Xce[2*i - 1,2], Xce[2*i,2])
  keep <- intersect(which(crosses[[i]]$allele1 != "None"),
                    which(crosses[[i]]$allele2 != "None"))
  notna <- crosses[[i]][keep,]
  equal <- ifelse(notna$allele1 == notna$allele2, 0, 1)
  het <- notna[which(equal == 1),]
  #het <- ddply(het_temp, .(gene_name), head, n = 1)
  
  #ind_temp <- intersect(which(het$pos > pos1),
  #                      which(het$pos < pos2))
  
  het %>% 
    group_by(gene_name) %>%
    summarize(n_unique = n_distinct(pos)) %>%
    top_n(20, n_unique)
  
  positions <- unique(het$pos)
  ind <- which(het$pos %in% positions)
  fin_het[[i]] <- ddply(het[ind,], .(pos), head, n=1)
  keep_dip <- which(crosses_dip[[i]]$pos %in% fin_het[[i]]$pos)
  fin_dip[[i]] <- ddply(crosses_dip[[i]][keep_dip,], .(pos), head, n=1)
  tally[[i]] <- fin_dip[[i]] %>% group_by(founder1, founder2) %>% tally()
  
  numhet <- c(numhet, length(unique(het$pos)))
  inXist<- c(inXist, length(positions) )
}

rem <- intersect(which(fin_dip[[9]]$founder1 == "C57BL6J"), which(fin_dip[[9]]$founder2 == "C57BL6J"))
final <- fin_het[[9]][-which(fin_het[[9]]$pos %in% fin_dip[[9]][rem,]$pos),]


#### Annotation ####

icr <- read.csv(paste0(dataSource,"/icr_genes.csv"), header=F, stringsAsFactors = F)
colnames(hets)[7] <- "gene_id"

cov_short <- readRDS(paste0(dataSource,"/cov_short_data.rds"))

genes <- colnames(cov_short$df)[grep("ENS", colnames(cov_short$df))]
filter(cov_short$df, RIX == 1) %>% 
  dplyr::select(Pup.ID, dir, Diet, one_of(genes)) -> rix1
         #one_of(annot_het_x$gene_id)) 
  
#keepCol <- c(1:36, which(colnames(rix1) %in% annot_het_x$gene_id))

meanEx <- colMeans(cov_short$df[,genes])
rankedGenes <- meanEx[sort.list(meanEx, decreasing = T)]
## icr genes only
keepCol <- grep("ENS", colnames(rix1))
rix1 <- rix1[,c(4:5, keepCol)]

# create ensemblGenome object for storing Ensembl genomic annotation data
ens <- ensemblGenome()

# read GTF file into ensemblGenome object

### plots for all hets
read.gtf(ens, paste0(dataSource, "/Mus_musculus.GRCm38.91.gtf"))
read.gtf(ens, paste0(dataSource, "/Mus_musculus.GRCm38.dna.chromosome.X.fa"))

genes <- fread(paste0(dataSource, "/Mus_musculus.GRCm38.91.gtf"), sep="\t")
setnames(genes, names(genes), c("chr","source","type","start","end","score","strand","phase","attributes") )

library(rtracklayer)
genes <- import(paste0(dataSource, "/Mus_musculus.GRCm38.91.gtf"))
which(genes$transcript_id %in% kal$transcript_id)
getTrans <- as.data.frame(genes[,c("transcript_id", "gene_id", "gene_name")])
keep <- which(kal$transcript_id %in% getTrans$transcript_id)
shkal <- kal[keep,]
getTrans <- getTrans[which(getTrans$transcript_id %in% shkal$transcript_id),]
getTrans %>% 
  group_by(transcript_id) %>%
  dplyr::slice(1) -> singleTrans
  
unique(getTrans$transcript_id)[1:10]
head(getTrans)
X_genes <- extractSeqids(ens, 'X')
my_gene <- getGenePositions(X_genes)
my_trans <- X_genes@ev$genes[ ,c("start","end","gene_id","gene_name")]
my_trans$transcript_id <- as.character(X_genes@ev$genes$transcript_id)


trans <- kal$transcript_id
annTrans <- sapply(1:5, function(x) extractTranscript(ens, trans[x]))
rbindlist(annTrans)
anno <- inner_join(my_gene, hets, by = "gene_id")
my_gene <- my_gene[,-which(apply(my_gene, 2, function(x) sum(is.na(x))) == nrow(my_gene))]

means <- rix1[,c(1:2,which(colnames(rix1) %in% my_gene$gene_id))] %>% 
  group_by(dir, Diet)

ts <- c()
for(i in 3:ncol(means)){
  a <- as.numeric(t(means[which(means$dir == "a"), i]))
  b <- as.numeric(t(means[which(means$dir == "b"), i]))
  ts <- c(ts, abs(t.test(b, a)$statistic))
}  


ord <- order(ts, decreasing = T)
icr_gene_ids <- colnames(means)[grep("ENS", colnames(means))]
max <- min(10, length(ord))
plotDat <- rix1 %>% select(one_of(c("dir", "Diet", icr_gene_ids[ord[1:max]]))) %>%
  gather(gene_id, mean, icr_gene_ids[ord[1:max]])
plotDat <- left_join(plotDat, my_gene[,c("gene_id","gene_name")], 
                     by = "gene_id")

ggplot(plotDat, aes(dir, mean, col=dir)) + 
  geom_boxplot() + 
  geom_jitter(width=0.1, aes(shape=Diet)) + 
  facet_wrap( ~ gene_name, ncol=max/2, scales = "free") +
  theme_classic()


#### using biostrings ####
library(BSgenome.Mmusculus.UCSC.mm10)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(Biostrings)

#txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
#genome <- BSgenome.Mmusculus.UCSC.mm10
#genome$chrX
#seq(genome$chrX)
#readstring

## needs Biostrings
xch <- readDNAStringSet(paste0(dataSource, "/Mus_musculus.GRCm38.dna.chromosome.X.fa"))
xsep <- paste(xch, collapse = "")
xep <- unlist(strsplit(xsep, ""))

## where the sequences I used to bwt came from
keepGenes <- (annot_het_x %>% filter(n_unique > 5)) %>%
  filter(gene_id %in% names(rankedGenes)[which(rankedGenes > 100)])
uniqueHets[[1]] %>% 
  filter(gene_name %in% keepGenes$gene_id) %>%
  filter(prob > 0.7) %>%
  mutate(start = pos - 17, end = pos + 17) -> rix1_hetPos


## remove hets that are too close to each other
distBtwHets <- unlist(sapply(1:nrow(uniqueHets[[1]]) - 1, function(x) uniqueHets[[1]]$pos[x+1] - uniqueHets[[1]]$pos[x]))
sort(distBtwHets)[1:10]
remove <- which(distBtwHets < 20)
uniqueHets[[1]][-remove, ] %>% 
  filter(gene_name %in% keepGenes$gene_id) %>%
  filter(prob > 0.7) %>%
  mutate(start = pos - 17, end = pos + 17) -> rix1_hetPos
##

seqs <- sapply(1:nrow(rix1_hetPos), function(x) xep[rix1_hetPos$start[x]:rix1_hetPos$end[x]])
seqs1 = seqs2 = seqs
seqs1[18,] <- sapply(1:nrow(rix1_hetPos), function(x) rix1_hetPos$allele1[x])
seqs2[18,] <- sapply(1:nrow(rix1_hetPos), function(x) rix1_hetPos$allele2[x])
seqs1 <- apply(seqs1, 2, function(x) paste(x, collapse = ""))
seqs2 <- apply(seqs2, 2, function(x) paste(x, collapse = ""))

nchar2 <- which(nchar(seqs1) != 35)
for(i in 1:length(nchar2)){
  split <- unlist(strsplit(seqs2[nchar2[i]],""))
  remCharL <- ifelse(length(split) %% 2 != 0, (length(split) - 35) / 2, (length(split) - 35 - 1) / 2)
  remCharR <- ifelse(length(split) %% 2 != 0, (length(split) - 35) / 2, ((length(split) - 35 - 1) / 2) + 1)
  if(remCharL != 0) split <- split[-c(1:remCharL, (length(split)-remCharR + 1):length(split))]
  else split <- split[-c((length(split)-remCharR + 1):length(split))]
  seqs2[nchar2[i]] <- paste(split, collapse="")
}
rix1_hetPos$seq1 <- seqs1
rix1_hetPos$seq2 <- seqs2

allSeqs <- c(seqs1, seqs2)

write.csv(rix1_hetPos, paste0(dataSource,"/rix1_hetPos.csv"))
#write.csv(seqs1, paste0(dataSource,"/rix1_hets_seqs1_v3.csv"))
#write.csv(seqs2, paste0(dataSource,"/rix1_hets_seqs2_v3.csv"))
write.csv(allSeqs, paste0(dataSource,"/rix1_hets_seqsAll_v3.csv"))
rix1_hetPos <- read.csv(paste0(dataSource,"/rix1_hetPos.csv"))

####################
#  icr expression  #
####################
icr_genes <- extractByGeneName(ens, as.vector(t(icr)))
#X_icr_genes <- extractSeqids(icr_genes, "X")
#xist <- extractByGeneName(icr_genes, c("Xist", "Tsix"))
my_gene <- getGenePositions(X_icr_genes)
my_gene <- getGenePositions(icr_genes)

means <- rix1[,c(1:2,which(colnames(rix1) %in% my_gene$gene_id))] %>% 
  group_by(dir, Diet)

ts <- c()
for(i in 3:ncol(means)){
  a <- as.numeric(t(means[which(means$dir == "a"), i]))
  b <- as.numeric(t(means[which(means$dir == "b"), i]))
  ts <- c(ts, abs(t.test(b, a)$statistic))
}  
  

vars <- rix1[,c(1:2,which(colnames(rix1) %in% my_gene$gene_id))] %>% group_by(dir, Diet) %>%
  summarise_all(funs(var))

ord <- order(ts, decreasing = T)
icr_gene_ids <- colnames(means)[grep("ENS", colnames(means))]
max <- min(20, length(ord))
plotDat <- rix1 %>% select(one_of(c("dir", "Diet", icr_gene_ids[ord[1:max]]))) %>%
                 gather(gene_id, mean, icr_gene_ids[ord[1:max]])
plotDat <- left_join(plotDat, my_gene[,c("gene_id","gene_name")], 
                     by = "gene_id")

ggplot(plotDat, aes(dir, mean, col=dir)) + 
  geom_boxplot() + 
  geom_jitter(width=0.1, aes(shape=Diet)) + 
  facet_wrap( ~ gene_name, ncol=max/2, scales = "free") +
  theme_classic()

