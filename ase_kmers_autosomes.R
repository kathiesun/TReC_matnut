options(digits=4)

setwd("C:/Users/Kathie/TReC_matnut/src")
library(tidyverse)
library(fdrtool)

source("prediction_functions.R")
source("summary_functions.R")
source("ase_summary_source.R")
source("stan_pheno_functions.R")


dir <- "C:/Users/Kathie/Dropbox\ (ValdarLab)"
gecco <- read.csv(file.path(dir, "de_results/fullGeccoRnaDump.csv"))

C_diet_4 = read.csv(file.path(dir, "variant_data/C_matrix_4diets.csv"), header=F)
C_diet_2 = read.csv(file.path(dir, "variant_data/C_matrix_2diets.csv"), header=F)
C_diet_3 = matrix(c(0.9659258,-0.2588190,-0.7071068,-0.2588190,0.9659258,-0.7071068),nrow=3,ncol=2)

#gene_count <- read.csv(file.path(dir,'/trec/gene_count_matrix.csv'))
#colnames(gene_count)[1] = "gene_id"
#gene_count$Gene.ID = do.call("rbind",(strsplit(as.character(gene_count$gene_id), "[|]")))[,1]
#gene_count$Gene.Name = do.call("rbind", (strsplit(as.character(gene_count$gene_id), "[|]")))[,2]

gene_count <- read.csv(file.path(dir,'/trec/gene_count_matrix_hetsOnly.csv'))
rownames(gene_count) = gene_count$Gene.Name
samples = colnames(gene_count)[grep("Pup", colnames(gene_count))]

matnut = read.csv(file.path(dir,'matnut_main/AllMice_GeneExpression_SSupdated_11.27.19.csv'))
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

#gene_count = data.frame(gene_count[-which(duplicated(gene_count$Gene.Name)),])
#rownames(gene_count) = gene_count$Gene.Name
counts = gene_count[, grep("Pup.ID", colnames(gene_count))]

colData = matnut[match(colnames(counts),matnut$ID),
                 c("Pup.ID","Breeding.Batch","Behavior.Batch","RIX","Reciprocal","Diet", 
                   "Dam.ID","ID","PO","DietRIX","DietRIXPOq")]
ie_genes = read.csv("C:/Users/Kathie/Dropbox (ValdarLab)/imprinted_genes/2014_05_14_allImprintedGenes.csv")

##########################

regFiles = list.files(pattern = "summary_mnt_50k_29dec2020.rds", file.path(dir, "variant_data/regression_outputs/"),
                      full.names = T)

regs = lapply(regFiles, readRDS)
mu_g = do.call("rbind", lapply(regs, function(x) x[[1]]$summary$mu_g))
mu_g = mu_g[,-intersect(grep("[.]",colnames(mu_g)), grep("^[a-z]", colnames(mu_g)))]
mu_g$Gene.Name = unlist(strsplit(rownames(mu_g),"_"))[c(F,T)]
mu_g$Pup.ID  = as.integer(unlist(strsplit(rownames(mu_g),"_"))[c(T,F)])
mu_g = mu_g %>% left_join(colData, by="Pup.ID")

mu_g$lower_opp = mu_g$lower
mu_g$upper_opp = mu_g$upper
mu_g$lower_opp[which(mu_g$PO < 0)] = 1-mu_g$upper[which(mu_g$PO < 0)]
mu_g$upper_opp[which(mu_g$PO < 0)] = 1-mu_g$lower[which(mu_g$PO < 0)]

imp_perc = mu_g %>% #filter(lower > 0.5 | upper < 0.5) %>%
  group_by(Gene.Name, RIX) %>% 
  summarize(over = sum(lower_opp > 0.5), under = sum(upper_opp < 0.5), total = n(),
            perc = max(over/total, under/total), 
            ratio = max(over/total, under/total) / min(over/total, under/total)) %>%
  filter(perc > 0.25, ratio > 3, total > 3) %>%
  arrange(desc(perc), desc(ratio),desc(total))
  #left_join(mu_g %>% group_by(Gene.Name) %>% tally(), by="Gene.Name") %>%
  #mutate(imp_perc = n.x / n.y) %>%
  #arrange(desc(imp_perc)) %>% 
  #filter(n.x > 3, imp_perc > 0.1)
imp_perc %>% filter(Gene.Name %in% ie_genes$mgi_symbol)

#write.csv(imp_perc, file.path(dir, "trec/priority_ase_genes_11jan2021.csv"))

## look specifically at known imprinted genes
tail((mu_g %>% 
  filter(Gene.Name %in% ie_genes$mgi_symbol) %>%
  group_by(Gene.Name, RIX) %>% 
  summarize(over = sum(lower_opp > 0.5), under = sum(upper_opp < 0.5), total = n(),
            perc = max(over/total, under/total), 
            ratio = max(over/total, under/total) / min(over/total, under/total)) %>%
  filter(perc > 0, !is.nan(ratio), total > 3) %>%
  arrange(desc(perc), desc(ratio),desc(total))), n=50)

mu_g_ie = mu_g %>% filter(Gene.Name %in% imp_perc$Gene.Name)


#imp_genes = imp_perc$Gene.Name[(which(imp_perc$Gene.Name %in% unlist(sig_genes_short)))]
keep_genes = imp_perc$Gene.Name[(which(imp_perc$Gene.Name %in% unlist(sig_genes_list)))]


print(mu_g %>% filter(lower > 0.55 | upper < 0.45) %>%
    group_by(Gene.Name) %>% tally() %>%
    left_join(mu_g %>% group_by(Gene.Name) %>% tally(), by="Gene.Name") %>%
    mutate(imp_perc = n.x / n.y) %>%
    arrange(desc(imp_perc)) %>% 
    filter(n.x > 3, imp_perc > 0.1), n=30)


gregg_genes = read.csv("../gregg_POgenes.csv")
colnames(gregg_genes)[1] = "Gene.Name"

#gregg_genes = gregg_genes %>% filter(sample != "E15", ie_status != "known") %>%
#  select(-ucsc_ID_SNP) %>% distinct()

gregg_ie  = unique(gregg_genes$Gene.Name)[which(unique(gregg_genes$Gene.Name) %in% ie_genes$mgi_symbol)]
gregg_ase = unique(gregg_genes$Gene.Name)[which(unique(gregg_genes$Gene.Name) %in% imp_perc$Gene.Name)]
gregg_des = unique(gregg_genes$Gene.Name)[which(unique(gregg_genes$Gene.Name) %in% unlist(sig_genes_list))]

#de_genes = read.table(file.path(dir, "trec/priority_deseq_genes_26jan2021.txt"))
de_genes = read.csv(file.path(dir, "trec/priority_deseq_genes_30jan2021.csv"))
ase_genes = read.csv(file.path(dir, "trec/priority_ase_genes_11jan2021.csv"))

all_genes = unique(c(gregg_genes$Gene.Name, de_genes$gene, ase_genes$Gene.Name, ie_genes$mgi_symbol))

################################################################


indiv_pups = list.files(file.path(dir, "mini/pup_haplo_blocks_by_CC_parent_jan2021"), pattern="haploBlocks", full.names = T)

phased_CC_haplotype = lapply(indiv_pups, readRDS)
tmp = do.call("rbind", lapply(indiv_pups, function(x) unlist(strsplit(x, "_"))))
names(phased_CC_haplotype) = paste0("Pup.ID_", tmp[,ncol(tmp)-1])

seg_regions = readRDS(file.path(dir,"mini/seg_regions_by_genotyping_perRIX_kmerBased_23nov2019.rds"))
masterSnps_files = list.files(file.path(dir, "variant_data/kmers_to_run"), pattern="masterSnps_chr", full.names = T)
count_files <- list.files(file.path(dir,"variant_data/kmer_counts"), ".csv", full.names = T)

## CC labels
#lab = read.csv("C:/Users/Kathie/Dropbox (ValdarLab)/variant_data/matched_v2_4jun2018.csv")
lab = read.csv(file.path(dir, "matnut_main/matched_v2_4jun2018.csv"))
lab = read.csv(file.path(dir, "matnut_main/CC_labels_paternity.csv"))
#lab %>% dplyr::select(one_of("Pup.ID","RRIX", "CC.1","CC.2")) %>%
#  filter(!is.na(Pup.ID)) -> lab





## Pup demographic info
pupInfo = matnut[,c(1,2,4,5,8,10,12,13,16,17)]
pupInfo %>% mutate(dir = gsub("[0-9]","", Reciprocal)) %>%
  mutate(Diet = gsub(" ","", Diet)) -> pupInfo

data_kmers = ratios_lst = list()

#if(j == 1){
  ## snp info

for(c in 1:19){

  masterSnps = readRDS(masterSnps_files[grep(paste0("chr",c,"_"), masterSnps_files)])
  masterSnps = masterSnps[[1]]
  masterSnps <- do.call("cbind", masterSnps)
  masterSnps$seq.consensus <- paste0(masterSnps$seq.end5, masterSnps$seq.end3)
  #useSnps = masterSnps %>% filter(seq.Gene %in% all_genes)
  
  ## count data 
  chr_files <- count_files[grep(paste0("chr", c,"_"), count_files)]
  ref_file <- chr_files[grep("ref", chr_files)]
  alt_file <- chr_files[grep("alt", chr_files)]
  reforig <- unique(read.csv(ref_file, header = T))
  altorig <- unique(read.csv(alt_file, header = T))
  reforig$X = altorig$X <- NULL
  #reforig = reforig %>% filter(k.mer %in% useSnps$seq.refseq)
  #altorig = altorig %>% filter(k.mer %in% useSnps$seq.altseq)
  reforig$pup = reforig$pup.id 
  altorig$pup = altorig$pup.id 
  
  data_kmers[[c]] = process_and_plot(chr=c, 
                                snp_info=masterSnps, 
                                sample_info=pupInfo, 
                                RIX_info=lab, 
                                ref_counts=reforig, alt_counts=altorig, 
                                phased_CC_haplotype=phased_CC_haplotype, 
                                use_gene=F,
                                problemPups = c(1404, 1716, 1371, 569, 1911, 1951, 1015),
                                seg_regions=seg_regions)
 
  #data_kmers$podietrix = paste(data_kmers$dir, data_kmers$DietRIX, sep="_")
  print(c)

  print(paste(c, "done"))
  
}



pvals = do.call("rbind", lapply(ratios_lst, function(z) 
  do.call("rbind", sapply(1:length(z), function(y) {
    tmpp = data.frame(do.call("rbind", lapply(z[[y]]$freq, function(x) {   
      tmp = data.frame(x$p.value)
      tmp$rix = names(z)[y]
      tmp
    })))
    tmpp$gene = rownames(tmpp)
    tmpp
  }, simplify=F))
))

pval_list = lapply(unique(pvals$gene), function(x) pvals[which(pvals$gene == x),])
#library(metaRNASeq)
#pval_comb = fishercomb(pval_list)

files = list.files(file.path(dir, "trec/data_kmers_from_process_and_plot/"), pattern=".txt",
                   full.names = T)


#data_kmers_list = saveRDS(file.path(dir, "trec/data_kmers_from_process_and_plot/data_kmers_list_out_29jan2021.rds"))
data_kmers_list = readRDS(file.path(dir, "trec/data_kmers_from_process_and_plot/data_kmers_list_out_29jan2021.rds"))


#lab = unique(do.call("rbind", lapply(data_kmers_list, function(x) 
#  x[,c("Pup.ID","CCs","RIX","Reciprocal","DamLine_NewCC_ID","SireLine_NewCC_ID","CC_lab","RRIX","CC.1","CC.2" )]
#  )))



#data_kmers_list = lapply(files, read.table, sep="\n")
#data_kmers_list = lapply(data_kmers_list, function(x) 
#  do.call("rbind", data.frame(apply(x, 1, function(y) strsplit(y," ")))))

#for(i in 1:length(data_kmers_list)){
#  colnames(data_kmers_list[[i]]) = data_kmers_list[[i]][1,]
#  data_kmers_list[[i]] = data_kmers_list[[i]][-1,]
#  rownames(data_kmers_list[[i]]) = NULL
#  data_kmers_list[[i]] = data.frame(data_kmers_list[[i]])
#  data_kmers_list[[i]][,c("CC_1","CC_2","sum","kRat","kLB","kUB")] = 
#    apply(data_kmers_list[[i]][,c("CC_1","CC_2","sum","kRat","kLB","kUB")], 2, as.numeric)
#  data_kmers_list[[i]]$Pup.ID = as.numeric(data_kmers_list[[i]]$Pup.ID)
#  data_kmers_list[[i]] = data_kmers_list[[i]] %>% left_join(lab[, -which(colnames(lab) == "RRIX")], by="Pup.ID")
#  data_kmers_list[[i]]$DamLine_NewCC_ID = unlist(strsplit(data_kmers_list[[i]]$DamLine_NewCC_ID, "/"))[c(T,F)]
#  data_kmers_list[[i]]$SireLine_NewCC_ID = unlist(strsplit(data_kmers_list[[i]]$SireLine_NewCC_ID, "/"))[c(T,F)]
#  data_kmers_list[[i]]$l_count = data_kmers_list[[i]]$CC_1
#  data_kmers_list[[i]]$r_count = data_kmers_list[[i]]$CC_2
#  data_kmers_list[[i]]$mat_count = data_kmers_list[[i]]$CC_1
#  data_kmers_list[[i]]$pat_count = data_kmers_list[[i]]$CC_2
  
#  data_kmers_list[[i]]$mat_count[which(data_kmers_list[[i]]$DamLine_NewCC_ID == data_kmers_list[[i]]$CC.2)] = 
#    data_kmers_list[[i]]$CC_2[which(data_kmers_list[[i]]$DamLine_NewCC_ID == data_kmers_list[[i]]$CC.2)]
#  data_kmers_list[[i]]$pat_count[which(data_kmers_list[[i]]$SireLine_NewCC_ID == data_kmers_list[[i]]$CC.1)] = 
#    data_kmers_list[[i]]$CC_1[which(data_kmers_list[[i]]$SireLine_NewCC_ID == data_kmers_list[[i]]$CC.1)]
#  data_kmers_list[[i]]$CC_1 = data_kmers_list[[i]]$mat_count
#  data_kmers_list[[i]]$CC_2 = data_kmers_list[[i]]$pat_count
#  data_kmers_list[[i]]$mat_count = data_kmers_list[[i]]$pat_count = NULL
#  data_kmers_list[[i]]$CC_lab = paste(data_kmers_list[[i]]$DamLine_NewCC_ID,data_kmers_list[[i]]$SireLine_NewCC_ID, sep="/")
#  data_kmers_list[[i]] = data_kmers_list[[i]] %>% select(one_of("seq.Gene","seq.Chromosome","seq.rsId","seq.Position",
#                            "CC_1","CC_2","sum","CC1_hap","CC2_hap",
#                            "Pup.ID","CCs","logSum","Breeding.Batch","Behavior.Batch","RIX","Reciprocal","Diet",
#                            "DamLine_NewCC_ID","Dam.ID","SireLine_NewCC_ID", "Sire.ID",
#                            "CC_lab","RRIX","pup_gene","DietRIX","CC.1","CC.2","l_count","r_count"))
#}

#######################################################
rixes = levels(pupInfo$RIX)
r=1

test_dat = do.call("rbind", lapply(data_kmers_list, function(x) x %>% filter(RRIX == rixes[r])))
test_dat = test_dat %>% group_by(Pup.ID, seq.Gene) %>% 
  summarize(mat_tot = sum(CC_1), pat_tot = sum(CC_2))
genes = unique(test_dat$seq.Gene)
pups = unique(test_dat$Pup.ID)
mat_tots = test_dat %>% select(c(Pup.ID, seq.Gene, mat_tot)) %>%
  spread(seq.Gene, mat_tot)
mat_tots = t(mat_tots)
colnames(mat_tots) = paste0("mat_",mat_tots[1,])
pat_tots = test_dat %>% select(c(Pup.ID, seq.Gene, pat_tot)) %>%
  spread(seq.Gene, pat_tot)
pat_tots = t(pat_tots)
colnames(pat_tots) = paste0("pat_",pat_tots[1,])
pat_tots = pat_tots[-1,]
mat_tots = mat_tots[-1,]
ase_sums = mat_tots + pat_tots

ase_tots = cbind(mat_tots,pat_tots)
ase_tots_comp = ase_tots[complete.cases(ase_tots),]
ase_sums_comp = ase_sums[complete.cases(ase_tots),]

colData = pupInfo %>% filter(RIX == rixes[r], Pup.ID %in% pups) %>%
  mutate(ID = paste0("Pup.ID_",Pup.ID), 
         RIX = factor(RIX, levels=c(1:4,6:10)),
         Diet = factor(Diet, levels=c("Standard","LowProtein","MethylEnriched","VitaminDDeficient")),
         dir = factor(dir)) 
rownames(colData) = colData$ID

colnames(ase_sums_comp) = gsub("mat_","Pup.ID_", colnames(ase_sums_comp))
reorder_cols = match(colData$ID, colnames(ase_sums_comp))[which(!is.na(match(colData$ID, colnames(ase_sums_comp))))]
ase_sums_comp = ase_sums_comp[,reorder_cols]

all(rownames(colData) == colnames(ase_sums_comp))

dds <- dds <- DESeqDataSetFromMatrix(countData = ase_sums_comp,
                                     colData = colData,
                                     design = ~ dir)
dds <- DESeq(dds)
res <- results(dds, name="dir_b_vs_a")

resLFC <- lfcShrink(dds, coef="dir_b_vs_a", type="apeglm")
resLFC

n <- length(pups)
f <- factor(rep(1:2,each=n))

theta.hat <- 1000 # rough initial estimate of dispersion
x <- model.matrix(~f)
#x = rep(1, n)
niter=5
for (i in 1:niter) {
  param <- cbind(theta.hat, ase_tots_comp)
  #param <- cbind(theta.hat, ase_tots_comp[,1:n])
  
  fit.mle <- apeglm(Y=ase_tots_comp, x=x, log.lik=NULL, param=param,
                    no.shrink=TRUE, log.link=FALSE, method="betabinCR")
  theta.hat <- bbEstDisp(success=ase_tots_comp, size=ase_sums_comp,
                         x=x, beta=fit.mle$map,
                         minDisp=.01, maxDisp=5000)
}

coef <- 2
xlab <- "mean of total counts"
plot(rowMeans(ase_sums_comp), fit.mle$map[,coef], log="x", xlab=xlab, ylab="log odds")


mle <- cbind(fit.mle$map[,coef], fit.mle$sd[,coef])
param <- cbind(theta.hat, ase_sums_comp)
fit2 <- apeglm(Y=ase_tots_comp, x=x, log.lik=NULL, param=param,
               coef=coef, mle=mle, threshold=0.7,
               log.link=FALSE, method="betabinCR")

ylim <- c(-1,1.5)
s.val <- svalue(fit2$thresh) # small-or-false-sign value
cols <- ifelse(s.val < .01, "red", "black")
plot(rowMeans(ase_sums_comp), fit2$map[,coef], main="apeglm",
     log="x", xlab=xlab, ylab="log odds", col=cols, ylim=ylim)
abline(h=0,col=rgb(1,0,0,.5))
fit2$map = cbind(fit2$map, abs(fit2$map[,coef]))
fit2$map = fit2$map[order(fit2$map[,coef+1], decreasing = T),]


ratios_lst = list()
for(c in 1:length(data_kmers_list)){
  ratios_lst[[c]] = run_stan_regress(data_kmers=data_kmers_list[[c]], 
                                     niter=10000, n.thin=5,  
                                     seg_regions=seg_regions,
                                     save_dir=NULL, 
                                     STZ=T, use_gene=F,
                                     no_theta=F, alpha=NULL,
                                     stan=F, stanMod = "ase_mu_g_regr.stan")#ase_mu_g_simple.stan
}

#saveRDS(data_kmers_list, file.path(dir, "trec/data_kmers_from_process_and_plot/data_kmers_list_out_29jan2021.rds"))

ratios_lst = readRDS(file.path(dir, "trec/data_kmers_from_process_and_plot/ratios_lst_out_29jan2021.rds"))


total_counts = lapply(data_kmers_list, function(x) 
  x %>% group_by(seq.Gene, seq.Chromosome, Pup.ID, Reciprocal, Diet) %>%
    summarize(sum_1 = sum(CC_1), sum_2 = sum(CC_2), sum_tot = sum(sum), 
              pos = mean(as.numeric(seq.Position)),
              ratio = sum_1 / sum_tot))

total_counts = do.call("rbind", total_counts)
total_counts$stringtie_cts = NA
total_counts$RRIX = gsub("[a-z]","",total_counts$Reciprocal)
total_counts$PO   = gsub("[0-9]","",total_counts$Reciprocal)



ase_total_counts = matrix(NA, ncol=length(colData$ID), nrow=nrow(gene_count))
colnames(ase_total_counts) = colnames(gene_count)
rownames(ase_total_counts) = rownames(gene_count)

for(p in colnames(ase_total_counts)){
  tmp = total_counts %>% filter(Pup.ID == gsub("Pup.ID_","",p))
  matches = match(tmp$seq.Gene,rownames(ase_total_counts))
  ase_total_counts[matches[!is.na(matches)],which(colnames(ase_total_counts) == p)] = tmp$sum_tot[!is.na(matches)]
}

ase_total_counts = read.csv(file.path(dir,"trec/ase_total_count_matrix.csv"))
rownames(ase_total_counts) = ase_total_counts$X
ase_total_counts$X = NULL

for(i in grep("Pup.ID", colnames(gene_count))){
  coln = grep("Pup.ID", colnames(gene_count))[i]
  pup = unlist(strsplit(colnames(gene_count)[coln], "_"))[c(F,T)]
  tmp = total_counts[which(total_counts$Pup.ID == pup),]
  #gene_count[match(tmp$seq.Gene, rownames(gene_count)), coln]
  total_counts$stringtie_cts[which(total_counts$Pup.ID == pup)] = gene_count[match(tmp$seq.Gene, rownames(gene_count)), coln]
}

plot_tot_counts = total_counts %>% filter(!is.na(stringtie_cts)) %>%
  group_by(Pup.ID) %>% mutate(norm_tot = sum_tot / sd(sum_tot), 
                              norm_str = stringtie_cts / sd(stringtie_cts),
                              log_norm_tot = log(norm_tot+0.5), 
                              log_norm_str = log(norm_str+0.5),
                              ratio = norm_str / (norm_tot + norm_str))
p_list = list()
for(i in unique(plot_tot_counts$seq.Chromosome)){
  pdf = plot_tot_counts %>% filter(seq.Chromosome == i)
  p_list[[i]] = ggplot(pdf, aes(x=pos)) + 
    geom_line(aes(y=log_norm_tot), col="blue") + 
    geom_line(aes(y=log_norm_str), col="red") + 
    #geom_line(aes(y=ratio), col="blue") + 
    #ylim(0, 5) + 
    theme_classic() + 
    facet_wrap( ~ RRIX)
}


counts_per_pup = lapply(unique(total_counts$Pup.ID), function(x)
  total_counts %>% filter(Pup.ID == x, !is.na(stringtie_cts)))

##############  fisher  ##################


binom_test_pvals = lapply(ratios_lst, function(x){   
  adj = do.call("rbind", sapply(1:length(x), function(y) {
    tmp = data.frame(do.call("rbind", lapply(x[[y]]$freq, function(z) 
      c(z$statistic, z$parameter, z$p.value, as.vector(z$conf.int), z$estimate))))
    colnames(tmp) = c("n.success", "n.trial", "p.value", "lower", "upper", "est")
    tmp$Gene.Name = rownames(tmp)
    tmp$RIX = names(x)[[y]]
    tmp
  }, simplify=F)) 
  adj$padj = as.numeric(p.adjust(adj$p.value, method = "BH"))
  adj
})

binom_test_pvals = do.call("rbind", binom_test_pvals)
binom_test_pvals$offset = abs(binom_test_pvals$est - 0.5)

binom_test_pvals_sig = binom_test_pvals %>% filter(padj < 0.05)

pval_list = lapply(unique(binom_test_pvals$Gene.Name), function(x) binom_test_pvals %>% filter(Gene.Name == x))

pmat = data.frame(apply(do.call("rbind", lapply(pval_list, function(g){
  data.frame(-2*sum(log(g$p.value)), nrow(g))
})), 2, as.numeric))
#pmat=do.call("rbind",pmat)
colnames(pmat) = c("fisher_stat", "n")
rownames(pmat) = unique(binom_test_pvals$Gene.Name)


hist(pmat$fisher_stat, breaks=200)
curve(dchisq(x, df = length(unique(binom_test_pvals$RIX))*2)*nrow(pmat), from=0, to=150, col="blue", add=T)

pmat$fisher_p = unlist(lapply(1:nrow(pmat), function(x)
  pchisq(pmat$fisher_stat[x], df=(as.numeric(paste(pmat$n[x]))*2),lower.tail = F)))

fdr = fdrtool(x = pmat$fisher_p, statistic = "pvalue")
pmat$padj_qval = as.numeric(fdr$qval)
pmat$padj_fdr  = as.numeric(fdr$lfdr)
pmat$padj = as.numeric(p.adjust(pmat$fisher_p, method = "BH"))


pmat$ie    = ifelse(rownames(pmat) %in% ie_genes$mgi_symbol, T, F)
pmat$gregg = ifelse(rownames(pmat) %in% gregg_genes$Gene.Name, T, F)
pmat$de    = ifelse(rownames(pmat) %in% de_genes$gene, T, F)
pmat$de    = ifelse(rownames(pmat) %in% genes, T, F)

pmat %>% arrange(fisher_p)

dim(pmat %>% filter(ie, padj_qval<0.05) %>% arrange(desc(fisher_stat)))
pmat %>% filter(!ie, padj_fdr<0.05, de) %>% arrange(desc(fisher_stat))

## Gnas ### Wars, Meg3 (not anymore)
## 100 ie but not de
## 114 de but not ie

genes = rownames(pmat %>% filter(padj_fdr < 1/length(unique(total_counts$seq.Gene)), !ie, n > 2) %>% arrange(padj_fdr))
genes = rownames(pmat %>% filter(padj_fdr < 1e-4, !ie, de) %>% arrange(desc(fisher_stat)))
genes = rownames(pmat %>% filter(!ie, gregg) %>% arrange(desc(fisher_stat)))
genes = rownames(pmat %>% filter(!ie, padj_qval<0.05) %>% arrange(desc(fisher_stat)))
genes = rownames(pmat %>% filter(padj_qval<0.05) %>% arrange(desc(fisher_stat)))
genes = rownames(pmat %>% filter(de) %>% arrange(padj_qval))

hi_rat = total_counts %>% ungroup() %>%
  group_by(seq.Gene) %>%
  mutate(ratio_m = sum_1 / sum_tot) %>%
  summarize(mean_rat = mean(ratio_m), abs_rat = abs(mean_rat - 0.5),
            n=n(), mean_counts = mean(sum_tot), max_count = max(sum_tot)) %>%
  arrange(desc(abs_rat)) %>%
  filter(abs_rat > 0.10, mean_rat < 1, n > 2, max_count > 50)
genes = hi_rat$seq.Gene[which(hi_rat$seq.Gene %in% genes)]
genes = genes[which(genes %in% hi_rat$seq.Gene)]

geneUse = "Eef2"
geneUse = "Inpp5f"
geneUse = "Igf2"
geneUse = "Gatad1"
geneUse = "Inpp5f"

p_list = list()

for(g in genes){
  geneUse = g
  plot_gene = total_counts %>% filter(seq.Gene == geneUse) %>%
    left_join(lab[,-which(colnames(lab) %in% c("RRIX", "Reciprocal", "Diet"))], by="Pup.ID") %>%
    mutate(ratio = sum_1 / sum_tot,
           PO = gsub("[0-9]", "", Reciprocal))
  plot_gene_long = gather(plot_gene, data_type, value, sum_1:sum_tot, factor_key=TRUE)
  jitter <- position_jitter(width = 0.1, height = 0)
  jitter2 <- position_jitter(width = 0.2, height = 0)
  
  colors <- c("Maternal" = "red", "Paternal" = "blue", "Total"="purple")
  colors <- c("firebrick2","royalblue","purple")
  
  p_list[[g]] = ggplot(data=plot_gene_long, aes(x=DamLine_NewCC_ID)) +   
    geom_boxplot(data=plot_gene_long %>% filter(data_type != "sum_tot"), 
                 aes(y=value, fill=data_type)) + 
    geom_point(data=plot_gene_long %>% filter(data_type == "sum_tot"), 
               aes(y=value, col=data_type), size = 3, alpha = 0.3, position=jitter) + 
    facet_grid(~ CCs, scales="free_x") + 
    scale_fill_manual(values=colors, name = "Allele-specific",
                      labels=c("Maternal","Paternal")) + 
    scale_color_manual(name="Total \nexpression", labels = "", values="purple") +
    labs(title = geneUse, y = "Counts", x="Maternal strain") +
    theme_classic()
}



pdf(file.path(dir,'trec/ase_notIE_gene_plots.pdf'))
p_list_non_ie
dev.off()

pdf(file.path(dir,'trec/ase_DE_gene_plots.pdf'))
p_list
dev.off()

pdf(file.path(dir,'trec/ase_IE_gene_plots.pdf'))
p_list
dev.off()

#Fbn1, Lrrc48, Ndnf
#p_list[[g]] = ggplot(plot_gene, aes(x=DamLine_NewCC_ID)) +   
#  geom_point(aes(y=sum_tot,col="Total"), size = 3, alpha = 0.3, position=jitter2) + 
#  geom_point(aes(y=sum_1,col="Maternal"), position=jitter) + 
#  geom_point(aes(y=sum_2,col="Paternal"), position=jitter) + 
#geom_point(aes(y = ratio)) +
#  facet_grid(~ CCs, scales="free_x") + 
#  labs(title = geneUse, color = "Expression",
#       y = "Counts", x="Maternal strain") +
#  scale_color_manual(values = colors) + 
#  theme_classic()


fdr = fdrtool(x = pmat$fisher_p, statistic = "pvalue")
pmat$padj_qval = fdr$qval
pmat$padj_fdr = fdr$pval
fdr = fdrtool(x = pmat$praw, statistic = "pvalue")
pmat$praw_qval = fdr$qval
pmat$praw_fdr = fdr$pval
pmat$padj_bh = p.adjust(pmat$padj_fdr, method = "BH")
pmat$praw_bh = p.adjust(pmat$praw_fdr, method = "BH")


  j=1
  resStan_tmp = ratios_lst[[j]]
  sum_df = summary(resStan_tmp)$summary
  simp_res = do.call("rbind", sapply(1:length(ratios_lst_simp), function(i){
    tmp = data.frame(summary(ratios_lst_simp[[i]])$summary)
    tmp$rix = names(ratios_lst_simp)[i]
    tmp
    }, simplify=F)
  )
  sum_df[-grep("ind|eta|_a|_b|weight", rownames(sum_df)),]
  
  #saveRDS(data_kmers, 
  #        file.path(dir, paste0("/regression_outputs/chr",i,"_data_mnt_29dec2020.rds")))
#} else {
#  data_kmers = readRDS(file.path(dir, paste0("/regression_outputs/chr",i,"_data_mnt_29dec2020.rds")))
#}
