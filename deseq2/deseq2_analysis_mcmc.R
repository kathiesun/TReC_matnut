options(stringsAsFactors = FALSE)
library(DESeq2)
library(tidyverse)
library(fdrtool)
library(hyper2)
library(BradleyTerryScalable)
library(coda)

setwd("C:/Users/Kathie/rna_seq/deseq2")
dir <- "C:/Users/Kathie/Dropbox\ (ValdarLab)/"
source("deseq2_functions.R")
source("manhattan_plot.R")
source("../kmerSearch/bradley_terry_mcmc_2.R")


######### read in data ##########
rem_seg_dds = readRDS(file.path(dir, "de_results/dds_lst_string_remSeg_exLowCounts_19feb2020.rds"))
rem_seg_res = compare_diets_PO(dds_lst=rem_seg_dds, fdrtool_adjust=T, 
                               adjust_p=F, ref_diet = "Standard")

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

## length(which(pmat_sig$g %in% sig_PO_df$gene))
## overlap of sets = 41

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
fitlist = p = datalist = list()
post_indiv = bt_post = list()
post_indiv_mcmc = bt_post_mcmc = list()
iter <- 100000

for(g in plot_genes[1:length(plot_genes)]){   #
  
  ## from meta-analysis below
  fish_p = pmat[pmat$g == g,]     #padj_qval or padj_bh
  
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

  #if(g %in% masterKmers$seq.Gene){
    minf = data.frame(masterKmers[which(masterKmers$seq.Gene == g),])
    c_this = list()
    
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
      require(plyr)
      permute_batches = lapply(1:3, function(x) getMatches(df=pd_use, id="Pup.ID", in_groups="RIX", out_groups="PO"))
      detach("package:plyr", unload=TRUE)
      bt_use = lapply(permute_batches, 
                      function(x) x %>% left_join(pd_use[,c("Pup.ID","mat_cc_found", "po", "RIX", "count")], 
                                                  by=c("Pup.ID_pos" = "Pup.ID")) %>% 
                      left_join(pd_use[,c("Pup.ID","mat_cc_found", "po", "RIX", "count")], 
                                                  by=c("Pup.ID_neg" = "Pup.ID")) %>%
                      mutate(diff = count.x - count.y, win = ifelse(diff > 0, mat_cc_found.x, mat_cc_found.y)))
      if(bayes){
        ### each mouse is a match
        post_indiv_mcmc[[g]] = run_bt(pd_use, player1 = "mat_cc_found", player2 = "pat_cc_found", iter=iter)
        
        ### permuted matched pairs
        post_list = lapply(bt_use, function(x) run_bt(x, player1 = "mat_cc_found.x", player2 = "mat_cc_found.y"))
        bt_post_mcmc[[g]] = post_list
      } 
      if(freq){
        ### each mouse is a match
        pd_use %>% group_by(RIX, mat_cc_found, pat_cc_found, win) %>% tally() %>% 
                 arrange(mat_cc_found, pat_cc_found) -> tmp_pd
        pd_4col = lapply(unique(tmp_pd$RIX), function(x){
          tmp_pd[which(tmp_pd$RIX == x)[1],c("mat_cc_found", "pat_cc_found")] -> out 
          out$found1_wins = sum(tmp_pd$n[which(tmp_pd$RIX == x & tmp_pd$win == out$mat_cc_found)])
          out$found2_wins = sum(tmp_pd$n[which(tmp_pd$RIX == x & tmp_pd$win == out$pat_cc_found)])
          out
        })
        pd_4col = do.call("rbind", pd_4col)
        pd_4col$mat_cc_found = factor(pd_4col$mat_cc_found, levels=LETTERS[1:8])
        pd_4col$pat_cc_found = factor(pd_4col$pat_cc_found, levels=LETTERS[1:8])
        single_pd <- btdata(pd_4col, return_graph = TRUE) 
        #summary(single_pd)
        
        single_pdfit = btfit(single_pd, a=1.1, MAP_by_component = T)
        summary_df_sing = data.frame(summary(single_pdfit, SE=T)$item_summary)
        summary_df_sing$stat = summary_df_sing$estimate / sqrt(summary_df_sing$SE)
        summary_df_sing$pval = dt(summary_df_sing$stat, df=nrow(summary_df_sing))

        ### permuted matched pairs
        pd_btfit = lapply(bt_use, function(x){   
          x %>% #select(RIX, mat_cc_found, pat_cc_found, win) %>% 
            group_by(RIX.x, mat_cc_found.x, mat_cc_found.y, win) %>% tally() %>% 
            arrange(mat_cc_found.x, mat_cc_found.y) -> bt_data
          bt_4col = do.call("rbind", lapply(unique(bt_data$RIX.x), function(x){
            bt_data[which(bt_data$RIX.x == x)[1],c("mat_cc_found.x", "mat_cc_found.y")] -> out 
            out$found1_wins = sum(bt_data$n[which(bt_data$RIX.x == x & bt_data$win == out$mat_cc_found.x)])
            out$found2_wins = sum(bt_data$n[which(bt_data$RIX.x == x & bt_data$win == out$mat_cc_found.y)])
            out
          }))
          bt_4col$mat_cc_found.x = factor(bt_4col$mat_cc_found.x, levels=LETTERS[1:8])
          bt_4col$mat_cc_found.y = factor(bt_4col$mat_cc_found.y, levels=LETTERS[1:8])
          
          ## https://cran.r-project.org/web/packages/BradleyTerryScalable/
          toy_btdata <- btdata(bt_4col, return_graph = TRUE) 
          #summary(toy_btdata)
  
          btfit(toy_btdata, a=1.1, MAP_by_component = T)
        })
        summary_df = lapply(pd_btfit, function(x){
          summary_df = data.frame(summary(x, SE=T)$item_summary)
          summary_df$stat = summary_df$estimate / sqrt(summary_df$SE)
          summary_df$pval = dt(summary_df$stat, df=nrow(summary_df))
          summary_df
        })
        post_indiv[[g]] = list(scalable_object = single_pdfit, summary_df = summary_df_sing)
        bt_post[[g]] = list(scalable_object = pd_btfit, summary_df = summary_df)
      }
    }
        
      #############################################
  #}
}


#mean(apply(post.p, 1, which.min)==3)
#mean(apply(post.p, 1, min)==post.p[,3])

#probability that each player is min
lapply(bt_post_mcmc, function(z) 
  lapply(z, function(y) 
  sapply(1:ncol(y$post.p), function(x) mean(apply(y$post.p, 1, min)==y$post.p[,x]))))

#means
lapply(bt_post_mcmc, function(z) 
  lapply(z, function(x) apply(x$post.p, 2, mean)))

#quantiles for each player
lapply(bt_post_mcmc, function(z) 
  lapply(z, function(x) apply(x$post.p, 2, quantile, probs=c(0.025, 0.975))))

#posterior number of groups
lapply(bt_post_mcmc, function(z) 
  lapply(z, function(x) table(apply(x$post.C,1,max))/iter))

#posterior probability of each configuration
lapply(bt_post_mcmc, function(z) 
  lapply(z, function(x) head(-sort(-table(apply(x$post.C,1,paste,collapse=",")))/iter)))
bt_post$Airn$summary_df

######################  extra BradleyTerry2 code  ##########################

#bt_data %>% select(-c("win","n")) %>% distinct() -> fin_bt_data

#tmp_winners = apply(fin_bt_data, 1, function(x){
#  mat_win = bt_data$n[intersect(which(bt_data$win == x["mat_cc_found.x"]), 
#                                which(bt_data$mat_cc_found.y == x["mat_cc_found.y"]))]
#  pat_win = bt_data$n[intersect(which(bt_data$win == x["mat_cc_found.y"]), 
#                                which(bt_data$mat_cc_found.x == x["mat_cc_found.x"]))]
#  return(c(mat_win, pat_win))
#}) 

#if(class(tmp_winners) == "list"){
#  winners = data.frame(do.call("rbind", tmp_winners))
#  
#} else {
#  winners = data.frame(t(tmp_winners))
#}
#colnames(winners) = c("mat_wins", "pat_wins")
#fin_bt_data$mat_wins = winners$mat_wins
#fin_bt_data$pat_wins = winners$pat_wins
#fin_bt_data$mat_cc_found = factor(fin_bt_data$mat_cc_found, levels=LETTERS[1:8])
#fin_bt_data$pat_cc_found = factor(fin_bt_data$pat_cc_found, levels=LETTERS[1:8])

#btModel1 <- BTm(outcome = data.matrix(bt_4col[,3:4]), #data.matrix(fin_bt_data[,c("mat_wins","pat_wins")]), 
#                player1 = mat_cc_found.x, player2 = mat_cc_found.y, #formula = ~ inheritance, 
#                formula = NULL, 
#                data = bt_4col, id = "founder")
#summary(btModel1)
#BTabilities(btModel1)

#pd_use %>% select(RIX, win) %>% group_by(RIX, win) %>% tally() -> bt_data
## BradleyTerryScalable (sparse dtasets, frequentist and Bayesian)

#pd_use %>% 
#  select(mat_cc_found.x, mat_cc_found.y, RIX.x) %>% distinct() -> founders
#founders = cbind(founders, t(apply(founders, 1, function(x) sort(c(x["mat_cc_found.x"], x["mat_cc_found.y"])))))
#colnames(founders)[(ncol(founders)-1):ncol(founders)] = c("founder1", "founder2")

