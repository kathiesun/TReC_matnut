options(stringsAsFactors = FALSE)
library(DESeq2)
library(tidyverse)
library(fdrtool)
library(ggrepel)
library(hyper2)
library(BradleyTerry2)
library(BradleyTerryScalable)
#library(plyr)

setwd("C:/Users/Kathie/rna_seq/deseq2")
#setwd("~/rna_seq/deseq/deseq2")
dir <- "C:/Users/Kathie/Dropbox\ (ValdarLab)/"
source("deseq2_functions.R")
source("manhattan_plot.R")


######### read in data ##########
rem_seg_dds = readRDS(file.path(dir, "de_results/dds_lst_string_remSeg_exLowCounts_19feb2020.rds"))
rem_seg_res = compare_diets_PO(dds_lst=rem_seg_dds, fdrtool_adjust=T, 
                               adjust_p=F, ref_diet = "Vitamin.D.Deficient")

problemPups = c(1404, 1716, 1371, 569, 1911, 1951, 1015)
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

regions_list = readRDS(file.path(dir,"mini/seg_regions_by_genotyping_perRIX_kmerBased_23nov2019.rds"))
haplofiles = list.files("C:/Users/Kathie/Dropbox (ValdarLab)/mini/pup_haplo_blocks_by_CC_parent_dec2019", 
                        pattern="haploBlocks.rds", full.names = T)
phased_par_CC_haplotypes = lapply(haplofiles, readRDS)
tmp = do.call("rbind", strsplit(haplofiles, "_"))
names(phased_par_CC_haplotypes) = tmp[,ncol(tmp)-1]
match_damID = read.csv("C:/Users/Kathie/Dropbox (ValdarLab)/matnut_main/MatNut_allRNAseq_sampleinfo.csv")
match_damID = match_damID[,c(5:6,8:9,15:16)]
colnames(match_damID) = c("RIX","dir","Maternal.Strain","Paternal.Strain","Maternal.CC","Paternal.CC")
match_damID %>% filter(Maternal.CC != "") %>% distinct() %>% 
  arrange(RIX, dir) -> match_damID
match_damID$Maternal.CC = do.call("rbind", strsplit(as.character(match_damID$Maternal.CC),"/"))[,1] 
match_damID$Paternal.CC = do.call("rbind", strsplit(as.character(match_damID$Paternal.CC),"/"))[,1] 

CC_lab = read.csv(file.path(dir, "variant_data/matched_v2_4jun2018.csv"))
CC_lab = CC_lab %>% dplyr::select(c("Pup.ID","RIX","CC.1","CC.2")) 
founders = read.csv("C:/Users/Kathie/Dropbox (ValdarLab)/data/cc_founderNameMap1410.txt")

####
ie_genes = read.csv(file.path(dir,"imprinted_genes/use_ie_gene_names.csv"))

annot = read.table("C:/Users/Kathie/Dropbox (ValdarLab)/matnut_main/Mus_musculus.GRCm38.96.gtf", skip=5, fill=T)
annot_genes = annot %>% filter(V15 == "gene_name") %>%
  dplyr::select(one_of(paste0("V",c(1,4,5,7,10,16)))) %>%
  distinct()

colnames(annot_genes) = c("Chr","Start","End","Strand","Gene.ID","Gene.Name")
annot_genes = annot_genes %>% mutate(Start = as.numeric(Start), End = as.numeric(End))
# 
if(length(which(duplicated(annot_genes$Gene.Name))) > 0){
  annot_genes = annot_genes[-which(duplicated(annot_genes$Gene.Name)),]
}

######## signif PO in multiple RIXs #########

##############  fisher  ##################
pcomb = sapply(names(rem_seg_res$res_PO), 
               function(x){
                 tmp = cbind(rix = x,
                             genes = rownames(rem_seg_res$res_PO[[x]]), 
                             praw = as.numeric(rem_seg_res$res_PO[[x]][,c("pval_raw")]),
                             pval = as.numeric(rem_seg_res$res_PO[[x]][,c("pvalue")]))
                 if(ncol(tmp) > 1) tmp
               })

pcomb = data.frame(do.call("rbind",pcomb))
pcomb$praw = as.numeric(paste(pcomb$praw))
pcomb$pval = as.numeric(paste(pcomb$pval))
pcomb$rix = factor(pcomb$rix, levels=c(1:4,6:10))     #paste0("rix",c(1:4,6:10))
all_genes = unique(pcomb$genes)
keep_genes = all_genes
length(which(all_genes %in% keep_genes))
pmat = lapply(keep_genes, function(g){
  test = pcomb %>% filter(genes == g) 
  stat_raw = -2*sum(log(test$praw))
  stat_adj = -2*sum(log(test$pval))
  n=nrow(test)
  data.frame(stat_raw, stat_adj, g, n)
})

pmat=do.call("rbind",pmat)
hist(pmat$stat_raw, breaks=900)
curve(dchisq(x, df = length(unique(pcomb$rix))*2)*nrow(pmat), from=0, to=150, col="blue", add=T)

pmat$padj = unlist(lapply(1:nrow(pmat), function(x)
  pchisq(pmat$stat_adj[x], df=(as.numeric(paste(pmat$n[x]))*2),lower.tail = F)))

pmat$praw = unlist(lapply(1:nrow(pmat), function(x)
  pchisq(pmat$stat_raw[x], df=(as.numeric(paste(pmat$n[x]))*2),lower.tail = F)))
fdr = fdrtool(x = pmat$padj, statistic = "pvalue")
pmat$padj_qval = fdr$qval
pmat$padj_fdr = fdr$pval
fdr = fdrtool(x = pmat$praw, statistic = "pvalue")
pmat$praw_qval = fdr$qval
pmat$praw_fdr = fdr$pval
pmat$padj_bh = p.adjust(pmat$padj_fdr, method = "BH")
pmat$praw_bh = p.adjust(pmat$praw_fdr, method = "BH")

###########################
pval=0.05
#pval=1
convert = data.frame(let = LETTERS[1:8], num = 1:8)
bl = c("Gm11408","Gm11408","Gm13340","Gm13340","Gm4076","Gm13339","Gm47099","Gm13341","Gm11407",
       "Gm11407","Gm11410","Gm28438","Gm28438","Gm28661","Gm26945")


sig_PO_0.05 <- lapply(names(rem_seg_res$res_PO), function(x){
  rem_seg_res$res_PO[[paste(x)]]$gene = rownames(rem_seg_res$res_PO[[paste(x)]])
  rem_seg_res$res_PO[[paste(x)]]$rix = x
  rem_seg_res$res_PO[[paste(x)]] %>% as.data.frame() %>% filter(padj < pval) %>% arrange(padj)
})
sig_PO_df = data.frame(do.call("rbind", sig_PO_0.05))
sig_PO_df %>% arrange(padj, pval_raw, gene, rix) %>% 
  filter(!gene %in% bl) -> sig_PO_df
sig_genes <- table(sig_PO_df$gene)


## 76 genes with signif PO in at least one RIX
## 71 genes with new data
## only one gene found in two RIXs
## no genes found in two RIXs after new data

#plot_genes = new_genes[-grep("Rik|[.]|Gm", names(sig_genes))]
plot_genes = names(sig_genes[which(sig_genes > 1)])
plot_genes = names(sig_genes)[which(names(sig_genes) %in% ie_genes$x)]
plot_genes =  c("Cyp17a1","Auts2","Kdm4b")


pmat_sig = pmat %>% filter(g %in% all_genes, !g %in% bl, padj_qval < pval) %>% arrange(-stat_adj) 
sig_ie = length(pmat_sig$g[(which(pmat_sig$g %in% ie_genes$x))])
plot_genes = unique(c(pmat_sig$g[(which(pmat_sig$g %in% ie_genes$x))],
                      pmat_sig$g[1:11]))
plot_genes = pmat_sig$g
## 107 genes with signif p-val after combining
## 109 genes with new data

plot_genes = unique(c(sig_PO_df$gene, pmat_sig$g))
## 139 with new data

## start plotting


fit=F
plot=F
ase=F
bt_freq=F
clusters_signif = c()
fitlist = p = datalist = bt_post = list()

for(g in plot_genes[19:length(plot_genes)]){   #
  
  ## from meta-analysis below
  fish_p = pmat[pmat$g == g,]     #padj_qval or padj_bh
  #pcomb_g = pcomb %>% filter(genes == g)  #pval
  #colnames(pcomb_g) = toupper(colnames(pcomb_g))
  
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
  masterKmers = readRDS(paste0(dir,"/variant_data/kmers_to_run/masterSnps_chr",find_gene$Chr,"_26aug2019.rds"))
  masterKmers = do.call("cbind", masterKmers[[1]])
  if(ase){
    alt_counts = read.csv(paste0(dir,"/variant_data/kmer_counts/chr",find_gene$Chr,"_alt_counts_19feb2020.csv"))
    ref_counts = read.csv(paste0(dir,"/variant_data/kmer_counts/chr",find_gene$Chr,"_ref_counts_19feb2020.csv"))
  }
  
  if(g %in% masterKmers$seq.Gene){
    minf = data.frame(masterKmers[which(masterKmers$seq.Gene == g),])
    c_this = list()
    if(ase){
      alt_counts %>% filter(k.mer %in% minf$seq.altseq) -> alt_counts
      ref_counts %>% filter(k.mer %in% minf$seq.refseq) -> ref_counts
    }
    
    for(i in 1:nrow(pd)){
      x = pd[i,]
      pup=gsub("[^0-9]", "", x["Pup.ID"])
      if(pup %in% names(phased_par_CC_haplotypes) & pup %in% samples_use$Pup.ID){
        rix = samples_use$RRIX[which(as.character(samples_use$Pup.ID) == pup)]
        po = paste(samples_use$PO_cat[which(as.character(samples_use$Pup.ID) == pup)])
        match = match_damID[which(match_damID$RIX == rix & match_damID$dir == po),]
        cc_order = unlist(lapply(phased_par_CC_haplotypes[[paste(pup)]]$founder_by_parent, function(x) unique(x$cc)))
        mat_cc = phased_par_CC_haplotypes[[pup]]$founder_by_parent[[which(cc_order == match$Maternal.CC)]]
        pat_cc = phased_par_CC_haplotypes[[pup]]$founder_by_parent[[which(cc_order == match$Paternal.CC)]]
  
        mat_cc_found = mat_cc$found[which(mat_cc$chr == find_gene$Chr & 
                                          mat_cc$start < find_gene$Start & mat_cc$end > find_gene$End)]
        pat_cc_found = pat_cc$found[which(pat_cc$chr == find_gene$Chr & 
                                          pat_cc$start < find_gene$Start & pat_cc$end > find_gene$End)]
        mat_cc_found = ifelse(length(mat_cc_found) > 0, mat_cc_found, "N")
        pat_cc_found = ifelse(length(pat_cc_found) > 0, pat_cc_found, "N")
        if(length(mat_cc_found) > 0 & length(pat_cc_found) > 0){
          c_this_tmp = data.frame(mat_cc_found = mat_cc_found,
                                  mat_cc = unique(mat_cc$cc),
                                  pat_cc_found = pat_cc_found,
                                  pat_cc = unique(pat_cc$cc), po=po)
          found_num_1 = convert$num[which(convert$let == c_this_tmp$mat_cc_found)]
          found_num_2 = convert$num[which(convert$let == c_this_tmp$pat_cc_found)]
          mat_snps = ifelse(length(found_num_1) > 0, sum(as.numeric(unlist(minf[paste0("sdp.",found_num_1)]))), NA)  #ifelse(length(found_num_1) > 0, found_num_1, NA)
          pat_snps = ifelse(length(found_num_2) > 0, sum(as.numeric(unlist(minf[paste0("sdp.",found_num_2)]))), NA)
          #sdp_use = convert$num[sort(which(convert$let %in% c(c_this_tmp$mat_cc_found,c_this_tmp$pat_cc_found)))]
          if(length(found_num_1) > 0 & length(found_num_2) > 0 & ase){
            check_counts = cbind(alt_counts %>% filter(pup.id == pup) %>% 
                    rename(c("counts"="alt_counts")) %>%
              left_join(minf[,c("seq.altseq",paste0("sdp.", found_num_1), paste0("sdp.", found_num_2))], 
                        by=c("k.mer" = "seq.altseq")),
              ref_counts %>% filter(pup.id == pup) %>% 
                select("counts") %>% rename(c("counts"="ref_counts"))) %>%
              filter(ref_counts > 3 | alt_counts > 3)
            makes_sense = c()
            mat_counts = pat_counts = c()
            count=1
            if(nrow(check_counts) > 0){
              for(j in 1:nrow(check_counts)){
                y = check_counts[j,]
                col = paste0("sdp.", c(found_num_1, found_num_2))
                if(all(y[col] == 0)){
                  makes_sense[j] = ifelse(y["ref_counts"] > 0, ifelse(y["alt_counts"] < 3, T, F), F)
                  #if(makes_sense[j]){
                  #  mat_counts[count] = pat_counts[count] = y["ref_counts"]
                  #  count=count+1
                  #}                
                } else if(all(y[col] == 1)){
                  makes_sense[j] = ifelse(y["alt_counts"] > 0, ifelse(y["ref_counts"] < 3, T, F), F)
                  #if(makes_sense[j]){
                  #  mat_counts[count] = pat_counts[count] = y["alt_counts"]
                  #  count=count+1
                  #}
                  
                } else {
                  makes_sense[j] = ifelse(all(y["alt_counts"] > 0 & y["ref_counts"] > 0), T, F)
                  if(makes_sense[j]){
                    mat_counts[count] = unlist(ifelse(y[col[1]] == 0, y["ref_counts"], y["alt_counts"]))
                    pat_counts[count] = unlist(ifelse(y[col[2]] == 0, y["ref_counts"], y["alt_counts"]))
                    count=count+1
                  }
                }
              }
              
              c_this_tmp$makes_sense = length(which(!makes_sense))
              c_this_tmp$mat_count = sum(unlist(mat_counts))  #do.call("sum", mat_counts)
              c_this_tmp$pat_count = sum(unlist(pat_counts)) 
              c_this_tmp$cc1_count = ifelse(po == "b", c_this_tmp$pat_count,c_this_tmp$mat_count)
              c_this_tmp$cc2_count = ifelse(po == "a", c_this_tmp$pat_count,c_this_tmp$mat_count)
            } else {
              c_this_tmp$makes_sense = 
                c_this_tmp$mat_count = c_this_tmp$pat_count = 
                c_this_tmp$cc1_count = c_this_tmp$cc2_count = NA
            }
            
          } else {
            c_this_tmp$makes_sense = 
              c_this_tmp$mat_count = c_this_tmp$pat_count = 
              c_this_tmp$cc1_count = c_this_tmp$cc2_count = NA
          }
          
          c_this[[i]] = cbind(Pup.ID = paste0("Pup.ID_",pup),
                              c_this_tmp, 
                              mat_snps = mat_snps, pat_snps = pat_snps)  #kmer_compare[c(found_num_1,found_num_2)]
                                                         #c(paste0("sdp.",found_num_1),paste0("sdp.",found_num_2))
          #names(c_this[[i]])[7:8] = c("sdp_1", "sdp_2")
        } else {
          c_this[[i]] = NULL
        }
        
      }
    }
    c_this <- do.call("rbind", c_this)
    if(ase){
      c_this %>% filter(mat_cc_found != "N", pat_cc_found != "N") -> c_this
      c_this %>% left_join(pd, by="Pup.ID") %>% 
        filter(!is.na(mat_cc)) %>%
        mutate(RIX = factor(RIX)) -> pd
      pd = pd[which(pd$makes_sense < length(makes_sense)/20),]
      pd %>% group_by(mat_cc, RIX, PO) %>%
        summarize(cc2_sum=sum(cc2_count), cc1_sum=sum(cc1_count), pat_sum=sum(pat_count), mat_sum=sum(mat_count), 
                  total = cc1_sum + cc2_sum, mat_rat = mat_sum/total, cc1_rat = cc1_sum/total, 
                  kLB = qbeta(0.05, cc1_sum + 1, cc2_sum + 1),
                  kUB = qbeta(0.95, cc1_sum + 1, cc2_sum + 1),) %>%
        arrange(RIX, PO) -> RIX_pd
      datalist[[g]] = RIX_pd
    } else {
      c_this %>% dplyr::select(-contains("_count"), -c("makes_sense")) -> c_this
      c_this %>% left_join(pd, by="Pup.ID") %>% 
        filter(!is.na(mat_cc)) %>%
        mutate(RIX = factor(RIX)) -> pd
      
      #################################################
      ##                bradley-terry                ##
      #################################################
      ## hyper2 (frequentist)
      pd %>% group_by(RIX) %>% summarize(mean = mean(count), sd = sd(count)) -> pd_means
      pd$norm_counts = apply(pd, 1, function(x) 
        (as.numeric(x["count"]) - pd_means$mean[which(pd_means$RIX == x["RIX"])])/pd_means$sd[which(pd_means$RIX == x["RIX"])]
      )
      pd %>% filter(mat_cc_found != "N", pat_cc_found != "N", mat_cc_found != pat_cc_found) %>%
        mutate(win = ifelse(norm_counts > 0, mat_cc_found, pat_cc_found)) -> pd_use
      if(dim(pd_use)[1] > 1){
        ## BradleyTerry2 (frequentist GLM)
        if(bt_freq){
          pd_use %>% select(RIX, mat_cc_found, pat_cc_found, win) %>% 
            group_by(RIX, mat_cc_found, pat_cc_found, win) %>% tally() %>% 
            arrange(mat_cc_found, pat_cc_found) -> bt_data
          bt_4col = do.call("rbind", lapply(unique(bt_data$RIX), function(x){
            bt_data[which(bt_data$RIX == x)[1],c("mat_cc_found", "pat_cc_found")] -> out 
            out$found1_wins = sum(bt_data$n[which(bt_data$RIX == x & bt_data$win == out$mat_cc_found)])
            out$found2_wins = sum(bt_data$n[which(bt_data$RIX == x & bt_data$win == out$pat_cc_found)])
            out
          }))
          bt_data %>% select(-c("win","n")) %>% distinct() -> fin_bt_data
          
          tmp_winners = apply(fin_bt_data, 1, function(x){
            mat_win = bt_data$n[intersect(which(bt_data$win == x["mat_cc_found"]), 
                                          which(bt_data$pat_cc_found == x["pat_cc_found"]))]
            pat_win = bt_data$n[intersect(which(bt_data$win == x["pat_cc_found"]), 
                                          which(bt_data$mat_cc_found == x["mat_cc_found"]))]
            return(c(mat_win, pat_win))
          }) 
          
          if(class(tmp_winners) == "list"){
            winners = data.frame(do.call("rbind", tmp_winners))
            
          } else {
            winners = data.frame(t(tmp_winners))
          }
          colnames(winners) = c("mat_wins", "pat_wins")
          fin_bt_data$mat_wins = winners$mat_wins
          fin_bt_data$pat_wins = winners$pat_wins
          fin_bt_data$mat_cc_found = factor(fin_bt_data$mat_cc_found, levels=LETTERS[1:8])
          fin_bt_data$pat_cc_found = factor(fin_bt_data$pat_cc_found, levels=LETTERS[1:8])
          
          btModel1 <- BTm(outcome = data.matrix(fin_bt_data[,c("mat_wins","pat_wins")]), 
                          player1 = mat_cc_found, player2 = pat_cc_found, #formula = ~ inheritance, 
                          formula = NULL, 
                          data = fin_bt_data, id = "inheritance")
          #summary(btModel1)
          #BTabilities(btModel1)
          
          #pd_use %>% select(RIX, win) %>% group_by(RIX, win) %>% tally() -> bt_data
          ## BradleyTerryScalable (sparse dtasets, frequentist and Bayesian)
          
          pd_use %>% filter(po == "a") %>%
            select(mat_cc_found, pat_cc_found, RIX) %>% distinct() -> founders
          founders = cbind(founders, t(apply(founders, 1, function(x) sort(c(x["mat_cc_found"], x["pat_cc_found"])))))
          colnames(founders)[(ncol(founders)-1):ncol(founders)] = c("founder1", "founder2")
          
          ## https://cran.r-project.org/web/packages/BradleyTerryScalable/
          toy_btdata <- btdata(bt_4col, return_graph = TRUE) 
          #summary(toy_btdata)
          
          pd_btfit <- btfit(toy_btdata, a=1.1, MAP_by_component = T)
          summary_df = data.frame(summary(pd_btfit, SE=T)$item_summary)
          summary_df$stat = summary_df$estimate / sqrt(summary_df$SE)
          summary_df$pval = dt(summary_df$stat, df=nrow(summary_df))
          
          #coef(pd_btfit)
          #vcov(pd_btfit)
          #btprob(pd_btfit)
          #fitted(pd_btfit, as_df = TRUE)    
        } else {
          post_list = run_bt(df=pd_use)
          ## using Wes's code
          #pd_use = read.table(file.path(dir, "/Inpp5f_norm_counts.txt"), header=T)
          dataset.numerator = table(pd_use$win)
          player.names <- names(dataset.numerator)    #LETTERS[1:8]
          player.names = sort(unique(c(pd_use$mat_cc_found, pd_use$pat_cc_found)))
          N <- length(dataset.numerator)
          #N <- length(player.names)
          comparisons <- t(combn(N, 2))
          rownames(comparisons) <- 1:nrow(comparisons)
          
          #comparisons
          convert = data.frame(let = player.names, num = 1:length(player.names))
          comparisons_let = apply(apply(comparisons, c(1,2), function(x) paste(convert$let[match(x, convert$num)])), 1, 
                                  function(y) paste(y, collapse=""))
          dataset.denominator = rep(0, length(comparisons_let))
          z = table(apply(pd_use, 1, function(x) paste(sort(c(x["mat_cc_found"], x["pat_cc_found"])), collapse="")))
          
          dataset.denominator[match(names(z), comparisons_let)] = -z
          
          #encode the comparisons as a matrix
          comparisons.matrix <- as.matrix(sapply(1:N, function(x){t(matrix(apply(comparisons==x, 1, any)))}))
          if(nrow(comparisons.matrix) != nrow(comparisons)){
            comparisons.matrix = t(comparisons.matrix)
          }
          dataset <- hyper2(c(1:N,lapply(1:nrow(comparisons), function(x){
            comparisons[x,]})), c(dataset.numerator, dataset.denominator), pnames=player.names)
          
          #MLE estimation using 'hyper2'
          #this gives the player ratings p, which are constrained to sum to 1
          mle.p <- maxp(dataset)
          
          ## gibbs sampler
          a <- 2
          b <- 1/2
          #specify number of iterations for the sampler
          iter <- 100000
          #create object to store results
          post.p <- matrix(NA, iter, N)
          colnames(post.p) <- player.names
          #set starting values based on MLE solution
          p <- mle.p
          names(p) <- NULL
          lambda <- a*N*p
          #store posterior hyperparameter a.star
          a.star <- a + dataset.numerator
          #iterate sampler
          for (i in 1:iter){
            #print every 1000 iterations
            if (i%%1000==0){
              print(i)
            }
            #sample latent variable z
            z <- rgamma(nrow(comparisons), -dataset.denominator, sapply(1:nrow(comparisons), function(x){sum(lambda[comparisons[x,]])}))
            #update posterior hyperparameter b.star and unconstrained player ratings lambda
            b.star <- b + sapply(1:N, function(x){sum(z[comparisons.matrix[,x]])})
            lambda <- rgamma(N, a.star, b.star)
            #calculate p given lambda
            p <- lambda/sum(lambda)
            #resample lambda given p
            lambda <- p * rgamma(1, N*a, 1)
            #store posterior sample
            post.p[i,] <- p
          }
          
          #compare marginal posterior means with MLE solutions
          #colMeans(post.p)
          #mle.p
          #calculate posterior probability that Karpov is the worst player
          #unlist(lapply(1:8, function(x) mean(apply(post.p, 1, which.min)==x)))
          #quantiles for each player
          #apply(post.p, 2, quantile, probs=c(0.025, 0.975))
          bt_post[[g]] = post.p
        }
        
        #############################################
        
        
        pd %>% group_by(mat_cc_found) %>% summarize(mean_norm = mean(norm_counts)) -> pd_norm_means
        #pd_norm_means = pd_norm_means %>% filter(!RIX %in% (pd_norm_means %>% group_by(RIX) %>% tally() %>% filter(n==1))$RIX)
        if(nrow(pd_norm_means) > 2){
          pdclus <- kmeans(pd_norm_means$mean_norm, 2, nstart = 20)
          pd_norm_means$clus = pdclus$cluster
          #cat = table(pd_norm_means$mat_cc_found)
          cat = pd_norm_means$mat_cc_found
          lo = table((pd_norm_means %>% filter(clus == 1))$mat_cc_found)
          names_use = sort(unique(cat))
          lo = unlist(lapply(names_use, function(x) ifelse(paste(x) %in% names(lo), lo[paste(x)], 0)))
          names(lo) = names_use
          hi = table((pd_norm_means %>% filter(clus == 2))$mat_cc_found)
          hi = unlist(lapply(names_use, function(x) ifelse(paste(x) %in% names(hi), hi[paste(x)], 0)))
          names(hi) = names_use
          k=length(unique(cat))
          n=length(cat)
          p_k = table(cat)/n    ## prob of getting each founder based on number of observed founders in total pool
          exi = n/2*p_k         ## expected n of founders in each hi or lo group
          
          set.seed(1)
          matlo = mathi = matrix(0, ncol=k, nrow=1000)
          for(i in 1:nrow(matlo)){
            #randlo = rmultinom(n=1000, size=n/2, prob=p_k)
            randlo = table(sample(cat, sum(lo)))
            randlo = unlist(lapply(names_use, function(x) ifelse(paste(x) %in% names(randlo), randlo[paste(x)], 0)))
            names(randlo) = names_use
            #randhi = rmultinom(n=1000, size=n/2, prob=p_k)
            randhi = table(cat) - randlo
            matlo[i,] = randlo
            mathi[i,] = randhi
          }
          
          sim_pr = data.frame(t(sapply(1:nrow(matlo), function(j){
            
            lo_p = dmultinom(matlo[j,], sum(lo), p_k, log = FALSE)
            hi_p = dmultinom(mathi[j,], n-sum(lo), p_k, log = FALSE)
            
            pr = lo_p * hi_p
            return(c(lo=lo_p, hi=hi_p, joint = pr))
          })))
          lo_p = dmultinom(lo, sum(lo), p_k, log = FALSE)
          hi_p = dmultinom(hi, n-sum(lo), p_k, log = FALSE)
          pr = lo_p * hi_p
          
          sim_p = quantile(sim_pr$joint, 0.05)
          signif = ifelse(pr < sim_p, T, F)
          clusters_signif[g] = signif
          #res.ftest <- var.test(mean_norm ~ clus, data = pd_norm_means)
          #res.ftest
          
          
          #ggplot(pd, aes(y=norm_counts, x=RIX, col=mat_cc_found)) + geom_point()
          
          
          #pd %>% group_by(mat_cc, RIX, PO) %>%
          #  summarize(cc2_sum=sum(cc2_count), cc1_sum=sum(cc1_count), pat_sum=sum(pat_count), mat_sum=sum(mat_count), 
          #            total = cc1_sum + cc2_sum, mat_rat = mat_sum/total, cc1_rat = cc1_sum/total, 
          #            kLB = qbeta(0.05, cc1_sum + 1, cc2_sum + 1),
          #            kUB = qbeta(0.95, cc1_sum + 1, cc2_sum + 1),) %>%
          #  arrange(RIX, PO) -> RIX_pd
          #datalist[[g]] = RIX_pd
        }
        
      }
    }
  }
}

      

    
    if(nrow(pd) > 0){
      
      rem_seg_res$res_PO = sapply(1:length(rem_seg_res$res_PO), function(n){
        rem_seg_res$res_PO[[n]]$RIX = names(rem_seg_res$res_PO)[n]
        rem_seg_res$res_PO[[n]]}, simplify=F)
      names(rem_seg_res$res_PO) = names(rem_seg_dds)
      
      stats = do.call("rbind", lapply(rem_seg_res$res_PO, function(x)
        data.frame(x[which(rownames(x) == g),])))
      stats$RIX = rownames(stats)
      which_rix = stats$RIX[which(stats$padj < 0.05)]
      
      #################  regression on PO effects  #######################
      if(fit){
        
        fitlist[[g]] = list()
        
        X <- data.matrix(apply(t(minf[,grep("sdp",colnames(minf))]), c(1,2), as.numeric))
        #X = unique(X[-which(apply(X, 1, function(x) all(x == 0))),-which(duplicated(t(X)))])
        #X.cov <- var(scale(X)) # standardize first
        #X.cov <- var(X)
        #X.mean <- apply(X,2,mean)
        #X.mean = rep(0, dim(X.cov)[1])
        #X.mah <- mahalanobis(X, X.mean, X.cov)
        
        euc_mat = dist(X, method = "euclidian")
        euc_mat=as.matrix(euc_mat, labels=TRUE)
        require(plyr)
        permute_batches = lapply(1:5, function(x) getMatches(df=pd, id="Pup.ID", in_groups="RIX", out_groups="PO"))
        detach("package:plyr", unload=TRUE)
        
        permute_batches = lapply(permute_batches, function(x) {
          x %>% left_join(pd, by=c("Pup.ID_pos"="Pup.ID")) %>%
            left_join(pd, by=c("Pup.ID_neg"="Pup.ID")) %>%
            mutate(diff = count.x-count.y, abs_diff = abs(diff)) %>%
            rename(c("count.x"="count_pos","count.y"="count_neg")) %>%
            select(-contains(".y")) -> x 
          x$euc_d = apply(x, 1, function(y){
            num1 = convert$num[match(y["mat_cc_found.x"], convert$let)]
            num2 = convert$num[match(y["pat_cc_found.x"], convert$let)]
            ifelse(y["mat_cc_found.x"] == y["pat_cc_found.x"], 0, 
                   euc_mat[paste0("sdp.",num1), paste0("sdp.",num2)])
          })
          x
        })
        
        all_batches = do.call("rbind", sapply(1:length(permute_batches), function(x){
          permute_batches[[x]]$batch = x
          permute_batches[[x]]
        }, simplify=F))
        all_batches$RIX.x = factor(all_batches$RIX.x, levels=c(1:4,6:10))
        if(length(unique(all_batches$euc_d)) > 1){
          euc_fit = lm("abs_diff ~ euc_d", all_batches)
          p[[g]] = ggplot(data=all_batches, aes(x=euc_d, y=abs_diff, col=RIX.x, shape=as.factor(batch))) + 
            geom_point() + labs(title=paste("Differences between matched pup pairs in", g),
                                subtitle=paste("R-sq:",round(summary(euc_fit)$adj.r.squared,3), 
                                               ", p=", formatC(summary(euc_fit)$coefficients["euc_d","Pr(>|t|)"], 
                                                               format = "e", digits = 2))) + 
            xlab("Euclidian distance between maternal haplotypes") + 
            ylab(expression(paste("|",Delta,"| expression for different POE"))) + 
            theme_classic() + 
            geom_abline(intercept=euc_fit$coefficients["(Intercept)"], slope=euc_fit$coefficients["euc_d"], col="blue")
          
        }
        
        match_mat = unique(pd[,grep("mat_cc|RIX", colnames(pd))])
        set.seed(20)
        Xclus <- kmeans(X, 2, nstart = 20)
        pc = prcomp(X)
        
        match_clus = cbind(convert, clus = matrix(Xclus$cluster, nrow=length(Xclus$cluster), ncol=1))
        match_clus = cbind(match_clus, pc$x[,which(apply(pc$x, 2, function(x) sum(abs(x))) > 1)])
       
        same_sdp_snp = list()
        same_clus = list()
        pd$mat_group = pd$pat_group = 0
        
        pd = cbind(pd, apply(match_clus[,grep("PC", colnames(match_clus))], 2, function(x){
          x[match(pd$mat_cc_found, match_clus$let)]
        }))
        
        for(r in unique(permute_batches[[1]]$RIX.x)){
          match_mat %>% filter(RIX == r) %>% select(mat_cc_found) -> mats
          sdp1 = convert$num[which(convert$let == mats$mat_cc_found[1])] 
          sdp2 = convert$num[which(convert$let == mats$mat_cc_found[2])]
          same_clus[[r]] = ifelse(Xclus$cluster[sdp1] == Xclus$cluster[sdp2], T, F)
          
          pd$mat_group[which(pd$RIX == r)] = match_clus$clus[match(pd$mat_cc_found[which(pd$RIX == r)], match_clus$let)]
          pd$pat_group[which(pd$RIX == r)] = match_clus$clus[match(pd$pat_cc_found[which(pd$RIX == r)], match_clus$let)]
          
          same_sdp_snp[[r]] = apply(rbind(minf[,paste0("sdp.",sdp1)], minf[,paste0("sdp.",sdp2)]), 2, 
                                    function(x) ifelse(x[1] == x[2], F, T))
        }
        
        if(any(unique(pd$mat_snps) > 0) | any(unique(pd$pat_snps) > 0)){
          snpmat = do.call("rbind", same_sdp_snp)
          test = sapply(2:ncol(snpmat), function(n){
            ifelse(any(sapply(1:(n-1), function(i) ifelse(identical(snpmat[,n-i], snpmat[,n]), T, F))  ), T, F)
          })
          snpmat = snpmat[,-which(c(F, test))]
          sdp_diff = data.frame(RIX=unique(permute_batches[[1]]$RIX.x),snpmat)
          if(!is.null(nrow(snpmat))) {
            snp_load = apply(snpmat, 1, sum)
          } else {
            snp_load = unlist(lapply(snpmat, function(x) ifelse(x, 1, 0)))
          }
          
          snp_load = data.frame(n = snp_load, RIX = names(snp_load))
          
          stats %>% left_join(sdp_diff, by="RIX") %>%
            left_join(snp_load, by="RIX")-> lm_stats
          lm_stats$ab_stat = abs(lm_stats$stat)
          lm_stats$sq_stat = lm_stats$stat * lm_stats$stat
          
          form = paste("pvalue ~", #log2FoldChange
                       paste(colnames(lm_stats)[grep("^X|snpmat", names(lm_stats))], collapse="+"))
          fit_diff = lm(form, data=lm_stats)
          fitlist[[g]]$diff = fit_diff
          #summary(fit_diff)
          
          fit_load = lm("pvalue ~ n", data=lm_stats)
          fitlist[[g]]$load = fit_load
          #summary(fit_load)
          
          
          fit_pat_snp = lm("count ~ pat_group", data=pd)
          fitlist[[g]]$fit_pat_snp = fit_pat_snp
          #summary(fit_pat_snp)
          
          form = paste("count ~", #log2FoldChange
                       paste(colnames(pd)[grep("PC", names(pd))], collapse="+"))
          fit_pc = lm(form, data=pd)
          fitlist[[g]]$pc = fit_pc
          summary(fit_pc)
        }
      }
      
      ##########   plotting   ###########
      if(plot){
        #c_this$het = ifelse(c_this$mat_snps == c_this$pat_snps, F, T)
        pd$xlab = paste0(pd$mat_cc, "(",pd$mat_cc_found,"/",pd$mat_snps,")")
        pd$which_rix = ifelse(pd$RIX %in% which_rix, T, F)
        
        for(fix in which(is.na(pd$mat_cc_found))){
          tmp = pd[fix,]
          replace = which(pd$PO == tmp$PO & pd$RIX == tmp$RIX & !is.na(pd$mat_cc_found))[1]
          pd[fix, 5:ncol(pd)] = pd[replace, 5:ncol(pd)]
        }
        pd = pd %>% filter(!is.na(mat_cc_found)) %>% left_join(stats, by="RIX")
        col_use = add.alpha(c("lightblue"), 0.4)
        pdf(paste0(dir, "/figures_and_updates/images_meetings/multRIX_PO_counts_",g,"_",Sys.Date(),".pdf"),
            width=9, height=7)
        par(oma = c(3, 1, 1, 1), mfrow=c(3,3))
        for(r in c(1:4,6:10)){
          if(r %in% unique(pd$RIX)){
            pd_use = pd %>% filter(RIX == r) %>% arrange(PO)
            pd_use$xlab = factor(pd_use$xlab)
            xlab = paste(pd_use$mat_cc[which(pd_use$po == "a")[1]], "x",  pd_use$pat_cc[which(pd_use$po == "a")[1]])
            numb <- as.numeric(unique(pd_use$pvalue))
            numbf = formatC(numb, format = "e", digits = 3)
            subt = paste0("p=",numbf)
            boxplot(pd_use$count ~ pd_use$xlab , 
                    outline=F, ylim=c(min(pd_use$count), max(pd_use$count)), 
                    main=paste(g, "RIX", r), xlab = xlab, sub=subt, 
                    col=ifelse(pd_use$which_rix, col_use, "white"))  
            points(x=jitter(as.numeric(pd_use$xlab)), y=pd_use$count, 
                   col=ifelse(pd_use$mat_snps > 0 , "seagreen1", "slategray4"), pch=16)  
          } else {
            plot(0,0)
          }
        }
        
        par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
        
        plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n", cex.lab = 6)
        title(paste0("Fisher's combined p-value=", numbf), line=-2)
        
        legend("bottom", c("Signif RIX", "Alt allele"), xpd = TRUE, horiz = TRUE, 
               inset = c(0,0), bty = "n", pch = c(15,20), 
               col = c(col_use, "seagreen1"),cex=1.5, text.col = "black")
        numb <- as.numeric(fish_p$padj_qval)
        numbf = formatC(numb, format = "e", digits = 3)
        dev.off()
      }
    }
    
    
  }
 
}
lapply(datalist, function(x) x[!is.null(x$mat_rat),])


lapply(fitlist, function(x) summary(x$load)$coefficients["n","Pr(>|t|)"])
lapply(fitlist, function(x) summary(x$diff)$coefficients[,"Pr(>|t|)"])    #grep("X", rownames(summary(x$diff)$coefficients))
lapply(fitlist, function(x) x$load$coefficients["n","Pr(>|t|)"])

mat_group_coef = lapply(fitlist, function(x){
  if(!is.null(x[["fit_mat_snp"]])){
    summary(x[["fit_mat_snp"]])$coefficients[,"Pr(>|t|)"]
  }   }) 

table(unlist(lapply(mat_group_coef, function(y)
  ifelse(y["mat_group"]<0.05, T, F))))

sort(unlist(lapply(mat_group_coef, function(y) y["mat_group"])))
###########  meta-analysis ############

sample_sizes = samples_use %>% group_by(RRIX) %>% tally()
sample_sizes$rix = paste0("rix",sample_sizes$RRIX)

rma_comb = sapply(names(rem_seg_res$res_PO), 
                  function(x){
                    #if(length(rownames(dds_string_compare$res_diet[[x]]$Standard)) > 0){
                    cbind(genes = rownames(rem_seg_res$res_PO[[x]]), 
                          rix = x,
                          rem_seg_res$res_PO[[x]][,c("log2FoldChange","lfcSE","stat")])
                    #}
                  })
rma_test = data.frame(do.call("rbind",rma_comb))


for(r in names(rem_seg_res$res_PO)){
  tmp = data.frame(cbind(genes = rownames(rem_seg_res$res_PO[[r]]), 
                         as.numeric(rem_seg_res$res_PO[[r]][,"baseMean"])))
  colnames(tmp) = c("genes", paste0("counts_",r))
  if(nrow(tmp) > 0){
    if(r == names(rem_seg_res$res_diet)[1]){
      baseMean_comb = tmp
    } else {
      baseMean_comb = full_join(baseMean_comb, tmp, by="genes")
    }
  }
}

baseMean_comb[,grep("[0-9]", colnames(baseMean_comb))] = apply(baseMean_comb[,grep("[0-9]", colnames(baseMean_comb))],c(1,2), as.numeric)
baseMean_comb$means = apply(baseMean_comb[,grep("[0-9]", colnames(baseMean_comb))], 1, function(x) mean(x, na.rm=T)) #, na.rm=T
#keep_genes = baseMean_comb$genes[which(baseMean_comb$means > 10)]
#keep_genes = keep_genes[-grep("Gm|Rik|[.]", keep_genes)] 
keep_genes = baseMean_comb$genes

##############  squared z-score method  ##################
rma_test$stat = as.numeric(paste(rma_test$stat))
rma_test$rix = factor(rma_test$rix, levels=unique(rma_test$rix))
all_genes = unique(rma_test$genes)
length(which(keep_genes %in% all_genes))
rma_test = rma_test %>% filter(genes %in% keep_genes)
test_mat = list()
for(r in levels(rma_test$rix)){
  test = rma_test %>% filter(rix == r) 
  zed = (test$stat - mean(test$stat))/sd(test$stat)
  test$zed = zed
  test$zed2 = zed^2
  test$stat2 = test$stat^2
  r = gsub("rix","",r)
  test$n = as.numeric(paste(sample_sizes$n[which(sample_sizes$RRIX == r)]))
  test_mat[[r]] = test
}
test_mat = do.call("rbind", test_mat)
test_mat %>% group_by(genes) %>%
  summarize(zed2_sum = sum(zed2),stat2_sum = sum(stat2)) -> sum_mat
sum_mat$pval_s = pchisq(sum_mat$stat2_sum, df=9, lower.tail = F)
sum_mat$pval_z = pchisq(sum_mat$zed2_sum, df=9, lower.tail = F)

hist(sum_mat$stat2_sum, breaks = 900, xlim=c(1,150))
curve(dchisq(x, df = 9)*nrow(sum_mat), from=0, to=150, col="blue", add=T)

hist(sum_mat$zed2_sum, breaks = 800, xlim=c(1,150))
curve(dchisq(x, df = 9)*nrow(sum_mat), from=0, to=150, col="blue", add=T)

#fdr <- fdrtool(sum_mat$pval, statistic="pvalue", plot=T)
sum_mat$padj_s = p.adjust(sum_mat$pval_s, method = "BH")
sum_mat$padj_z = p.adjust(sum_mat$pval_z, method = "BH")
sum_mat_clean = sum_mat %>% arrange(pval_s)

sum_mat_clean = sum_mat_clean %>% left_join(baseMean_comb[,c("genes","means")], by="genes")
plot_genes = sum_mat_clean$genes[which(sum_mat_clean$padj_s < 0.05)]



## compare with sum_mat_clean with fisher results (pmat)

pmat_sig[unique(c(1:18,which(pmat_sig$g %in% ie_genes$x))),] %>% arrange(padj_qval)
## 6
ie_n = length(which(ie_genes$x %in% pmat$g))
mat = matrix(c(sig_ie,nrow(pmat_sig)-sig_ie,ie_n-sig_ie,nrow(pmat)-ie_n-(nrow(pmat_sig)-sig_ie)),
             nrow=2, ncol=2, byrow = T)
chisq.test(mat)
fisher.test(mat)

####################   simulations   #########################
randmat = matrix(NA, nrow=10000, ncol=8)
for(i in 1:10000){
  randmat[i,] = runif(8, 0,1)
}
randmat = data.frame(randmat)
randmat$stat = unlist(apply(randmat, 1, function(g) -2*sum(log(g))))
randmat$p = pchisq(randmat$stat, df=16, lower.tail = F)
fdr = fdrtool(x = randmat$p, statistic = "pvalue")

randmat$p = pchisq(pmat$stat_raw, df=length(unique(pcomb$rix))*2,lower.tail = F)
pmat$padj_qval = fdr$qval

