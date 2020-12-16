options(stringsAsFactors = FALSE)
library(DESeq2)
library(tidyverse)
library(metafor)
library(sva)
library(fdrtool)
library(ggrepel)
library(matrixStats)
library(gridExtra)
library(latticeExtra)


setwd("C:/Users/Kathie/rna_seq/deseq2")
setwd("~/rna_seq/deseq2")
dir="/nas/depts/006/valdar-lab/users/sunk/"
dir <- "C:/Users/Kathie/Dropbox\ (ValdarLab)/"
source("deseq2_functions.R")
######## READ IN DATA #########
countData <- as.matrix(read.csv(file.path(dir,"string_pipe_out/gene_count_matrix.csv"), row.names="X"))

annot = read.table("/nas/depts/006/valdar-lab/users/sunk/mm10/ref_sequence/Mus_musculus.GRCm38.96.gtf", skip=5, fill=T)
  
annot_genes = annot %>% filter(V15 == "gene_name") %>%
	dplyr::select(one_of(paste0("V",c(1,4,5,7,10,16)))) %>%
	distinct()

annot = read.table("C:/Users/Kathie/Dropbox (ValdarLab)/matnut_main/Mus_musculus.GRCm38.96.gtf", skip=5, fill=T)
annot_genes = annot %>% filter(V15 == "gene_name") %>%
  dplyr::select(one_of(paste0("V",c(1,4,5,7,10,16)))) %>%
  distinct()
 
colnames(annot_genes) = c("Chr","Start","End","Strand","Gene.ID","Gene.Name")
annot_genes = annot_genes %>% mutate(Start = as.numeric(Start), End = as.numeric(End))

if(length(which(duplicated(annot_genes$Gene.Name))) > 0){
  annot_genes = annot_genes[-which(duplicated(annot_genes$Gene.Name)),]
}

rownames(countData) = annot_genes$Gene.Name[match(rownames(countData), annot_genes$Gene.ID)]
countData = countData[which(rownames(countData) %in% annot_genes$Gene.Name[which(annot_genes$Chr %in% c(1:19,"X"))]),]

samples <- read.csv(file.path(dir, "variant_data", "2015-10_expression_pups.csv"))
rownames(samples) = samples$Pup.ID
samples$Diet <- gsub(" $", "", samples$Diet)
samples$Diet <- factor(samples$Diet, levels=c("Standard", "Low Protein","Methyl Enriched","Vitamin D Deficient"))
samples$PO <- ifelse(factor(gsub("[0-9]","", samples$RIX)) == "a", 0.5, -0.5)
samples$PO_cat <- factor(gsub("[0-9]","",samples$RIX))

problemPups = c(1404, 1716, 1371, 569, 1911, 1951, 1015)

use_pups <- as.numeric(gsub("Pup.ID_","",colnames(countData)))
use_pups <- use_pups[-which(use_pups %in% problemPups)]

colData <- samples %>% 
	dplyr::filter(Pup.ID %in% use_pups) %>%
	dplyr::select(one_of("Pup.ID","RRIX","Diet","PO","PO_cat","Breeding.Batch","Behavior.Batch")) %>%
	mutate(Breeding.Batch  = as.factor(Breeding.Batch), Behavior.Batch = as.factor(Behavior.Batch))
rownames(colData) = paste0("Pup.ID_", colData$Pup.ID)
all(rownames(colData) %in% colnames(countData))
countData <- countData[, rownames(colData)]
all(rownames(colData) == colnames(countData))

#rs <- rowSums(countData)
#countData = countData[-which(rs < 10), ]
regions_list = readRDS(file.path(dir,"mini/seg_regions_by_genotyping_perRIX_kmerBased_23nov2019.rds"))

########### EDA ###########

hist(log10(rs + 1))
abline(v=1, col="blue", lwd=3)
length(which(rs < 10))

par(mfrow=c(1,2))
plot(countData[,1:2])
logcts <- log10(countData + 1)
plot(logcts[,1:2])

hist(logcts[,1])

rv <- rowVars(logcts)
o <- order(rv, decreasing=TRUE)[1:500]
pc <- prcomp(t(logcts[o,]))
plot(pc$x[,1:2])
dat <- data.frame(pc$x[,1:2], RRIX=factor(colData$RRIX), Diet=factor(colData$Die))
ggplot(dat, aes(PC1, PC2, col=RRIX)) + geom_point() + theme_classic()
ggplot(dat, aes(PC1, PC2, col=Diet)) + geom_point() + theme_classic()


plot(pc$sdev[1:10]^2 / sum(pc$sdev^2), type="b", ylab="% Var")


######## RUN DESEQ2 #########

## new way
#dds_string_sva = run_DESeq2(colData, countData, sva=T, interact=F)
rem_seg_dds  = run_DESeq2(colData, countData, sva=T, interact=F, 
                          contrast="PO", annot=annot_genes, regions_list=regions_list)
keepSeg_dds_diet  = run_DESeq2(colData, countData, sva=T, interact=F, 
                             contrast="Diet", annot=annot_genes, regions_list=NULL)
keep_seg_dds = run_DESeq2(colData, countData, sva=T, interact=F, 
                          contrast="Diet", annot=annot_genes, regions_list=NULL)
rem_seg_res  = compare_diets_PO(dds_lst=rem_seg_dds, fdrtool_adjust=T, 
                                adjust_p=F, ref_diet = "Vitamin.D.Deficient")
remSeg_res_diet = compare_diets_PO(dds_lst=remSeg_dds_diet, fdrtool_adjust=T, 
                                adjust_p=F, ref_diet = "Vitamin.D.Deficient")
keepSeg_res_diet = compare_diets_PO(dds_lst=keepSeg_dds_diet, fdrtool_adjust=T, 
                                   adjust_p=F, ref_diet = "Vitamin.D.Deficient")


remSeg_dds_interDiet  = run_DESeq2(colData, countData, sva=T, interact=T, 
                          contrast="Diet", annot=annot_genes, regions_list=regions_list)
remSeg_res_diet = compare_diets_PO(dds_lst=remSeg_dds_interDiet, fdrtool_adjust=T, 
                                   adjust_p=F, ref_diet = "Vitamin.D.Deficient")
interResLst = lapply(remSeg_dds_interDiet, function(x){
  tmp = results(x)
  dtmp = data.frame(tmp)
  dtmp$gene = rownames(tmp)
  fdr <- fdrtool(dtmp$stat, statistic="normal", plot=F)
  dtmp$pval_raw = dtmp$pvalue
  dtmp$padj_deseq2 = dtmp$padj
  dtmp$pvalue = fdr$pval
  dtmp$padj  <- p.adjust(dtmp$pvalue, method = "BH")
  dtmp %>% arrange(pvalue)
})

interResLst_sig = lapply(names(interResLst), function(x) interResLst[[x]] %>% filter(padj_deseq2 < 0.10) 
                         %>% arrange(pval_raw) %>% mutate(rix = x))
interResLst_sig_df = do.call("rbind", interResLst_sig)
interResLst_sig_df = interResLst_sig_df %>% left_join(annot_genes, by=c("gene" = "Gene.Name"))
write.csv(interResLst_sig_df, file.path(dir, "figures_and_updates/interaction_signif_0.05.csv"))
plotDispEsts(rem_seg_dds$`4`)
par(mfrow=c(3,3))
for(i in names(rem_seg_res$res_PO)){  
  plot(rem_seg_res$res_PO[[paste(i)]]$baseMean+1, -log10(as.numeric(rem_seg_res$res_PO[[paste(i)]]$pval_raw)),
       log="x", xlab="mean of normalized counts",
       ylab=expression(-log[10](pvalue)),
       ylim=c(0,30),
       cex=.4, col=rgb(0,0,0,.3), main=paste("RIX",i))
}


#dds_string_sva = readRDS(file.path(dir, "de_results/allrix_string_deseq_sva_noX_25oct2019.rds"))
#saveRDS(dds_string_sva, file.path(dir, "de_results/dds_lst_string_sva_exLowCounts_6nov2019.rds"))
saveRDS(rem_seg_dds, file.path(dir, "de_results/dds_lst_string_remSeg_exLowCounts_19feb2020.rds"))
#saveRDS(remSeg_dds_diet, file.path(dir, "de_results/dds_lst_string_remSeg_Diet_exLowCounts_19nov2019.rds"))
rem_seg_dds = readRDS(file.path(dir, "de_results/dds_lst_string_remSeg_exLowCounts_19nov2019.rds"))
dds_string_sva = readRDS(file.path(dir, "de_results/dds_lst_string_sva_inclX_30oct2019.rds"))
dds_string_compare = compare_diets_PO(dds_lst=dds_string_sva, fdrtool_adjust=T, 
                                  adjust_p=F, ref_diet = "Vitamin.D.Deficient")

res_string_sva$res_diet$rix6$Standard[which(rownames(res_string_sva$res_diet$rix6$Standard) == "Gprasp1"),]
########## old way ##############
dds_string = res_PO = res_diet = list()
dds_string_sva = list()
countData_list = list()
for(r in unique(colData$RRIX[order(colData$RRIX)])){
  colData_use = data.frame(colData %>% filter(RRIX==r))
  countData_use = countData[,which(colnames(countData) %in% paste0("Pup.ID_", colData_use$Pup.ID))]
  if(!is.null(regions_list)){
    annot_genes$Start = as.numeric(annot_genes$Start)
    annot_genes$End = as.numeric(annot_genes$End)
    
    regions = regions_list[[paste(r)]]
    remove_genes = apply(regions, 1, function(x) 
      intersect(which(annot_genes$Chr %in% x["Chr"]), 
                intersect(which(as.vector(annot_genes$Start) >= as.numeric(x["start"]) - 10000),
                          which(as.vector(annot_genes$End) <= as.numeric(x["end"]) + 10000))))
    remove_genes = unique(unlist(remove_genes))
    remove_names = annot_genes$Gene.Name[remove_genes]
    countData_use[which(rownames(countData_use) %in% remove_names), ] = NA
    countData_list[[paste(r)]] = countData_use
  }
}
countData = do.call("cbind", countData_list)
rem = apply(countData,1, function(x) ifelse(length(which(is.na(x))) > 0, T, F))
length(which(rem))
countData = countData[-which(rem),rownames(colData)]

identical(colnames(countData), rownames(colData))

  dds = DESeqDataSetFromMatrix(countData = countData, colData = colData, 
                               design = ~ RRIX + Diet + PO_cat)    # + Diet:PO)
  countData = countData[rowSums(counts(dds)) >= 10,]
  dds <- estimateSizeFactors(dds)
  dat = counts(dds, normalized=T)
  dat = dat[rowMeans(dat) > 1,]
  mod  <- model.matrix(~   RRIX + Diet + PO_cat, colData(dds))
  mod0 <- model.matrix(~   1            , colData(dds))
  svseq <- svaseq(dat, mod, mod0, n.sv = 1)  
  ddssva = dds
  ddssva$SV1 <- svseq$sv[,1]
  design(ddssva) <- as.formula(~ SV1 + RRIX + Diet + PO_cat)
  dds = ddssva
  dds = DESeq(dds)
  dds = nbinomWaldTest(dds)
  resDiet <- results(dds, contrast=c("Diet","Standard","Vitamin.D.Deficient"))
  if(length(which(is.na(resDiet$pvalue))) > 0) resDiet = resDiet[-which(is.na(resDiet$pvalue)),]
  fdr <- fdrtool(resDiet$stat, statistic="normal", plot=F)#, cutoff.method = "locfdr")
  #fdr <- fdrtool(res$pvalue, statistic="pvalue")#, cutoff.method = "locfdr")
  
  resDiet[,"pval_tool"] = fdr$pval
  resDiet[,"padj_tool"]  <- p.adjust(resDiet$pval_tool, method = "BH")
  res = data.frame(resDiet) %>%
    mutate(gene = rownames(resDiet)) %>%
    #filter(baseMean > 10) %>%     #, padj_tool == min(padj_tool)
    arrange(padj)
  res[1:100,]
  plotGenes = res$gene[1:6]
  
  par(mar=c(8,5,2,2))
  boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)
  plotDispEsts(dds)
  
  plot(resDiet$baseMean+1, -log10(resDiet$pvalue),
       log="x", xlab="mean of normalized counts",
       ylab=expression(-log[10](pvalue)),
       ylim=c(0,30),
       cex=.4, col=rgb(0,0,0,.3))
  
  
  ######## if running LRT ##########
  dds <- estimateDispersions(dds)
  ddsPO <- nbinomLRT(dds, reduced = ~ Diet)
  resPO <- results(ddsPO, name="PO")
  resPO_lim = resPO[-grep("Gm|Rik|[.]", rownames(resPO)),]
  resPO_lim[order(resPO_lim$padj),]
  
  ddsDiet <- nbinomLRT(dds, reduced = ~ PO)
  resDiet <- results(ddsDiet) #, contrast=c("Diet","Standard","Vitamin.D.Deficient"))
  resDiet_lim = resDiet[-grep("Gm|Rik|[.]", rownames(resDiet)),]
  resDiet_lim[order(resDiet_lim$padj),]
  
  
  dds = DESeq(dds)
  
  dat = counts(dds, normalized=T)

  dat = dat[rowMeans(dat) > 1,]
  mod  <- model.matrix(~ Diet + PO_cat, colData(dds))
  mod0 <- model.matrix(~   1, colData(dds))
  svseq <- svaseq(dat, mod, mod0, n.sv = 1)  
  ddssva = dds
  ddssva$SV1 <- svseq$sv[,1]

  design(ddssva) <- ~ SV1 + Diet + PO_cat 
  ddssva = DESeq(ddssva)
  ddssva = nbinomWaldTest(ddssva)
  
  dds = DESeq(dds)
  dds_string_sva[[paste0("rix",r)]] = ddssva
  ddssva = dds_string_sva[[paste0("rix",r)]]
  res <- results(ddssva, contrast = c("PO_cat","a","b"))
  res = res[-which(is.na(res$pvalue)),]
  
  fdr <- fdrtool(res$stat, statistic="normal", plot=F)#, cutoff.method = "locfdr")
  #fdr <- fdrtool(res$pvalue, statistic="pvalue")#, cutoff.method = "locfdr")
  
  res[,"pval_tool"] = fdr$pval
  res[,"padj_tool"]  <- p.adjust(res$pval_tool, method = "BH")
  
  
  
  res_PO[[paste0("rix",r)]] <- res
  res_diet[[paste0("rix",r)]] <- list()
  
  if(length(grep("Standard", colData_use$Diet)) > 0 & length(grep("Vitamin.D.Deficient", colData_use$Diet)) > 0){
     res = results(ddssva, contrast=c("Diet", "Vitamin.D.Deficient", "Standard"))
     res = res[-which(is.na(res$pvalue)),]
     fdr <- fdrtool(res$stat, statistic="normal", plot=F)#, cutoff.method = "locfdr")
     res[,"pval_tool"] = fdr$pval
     res[,"padj_tool"]  <- p.adjust(res$pval_tool, method = "BH")
     res_diet[[paste0("rix",r)]]$Standard = res
  }
  if(length(grep("Low.Protein", colData_use$Diet)) > 0 & length(grep("Vitamin.D.Deficient", colData_use$Diet)) > 0){
    res = results(ddssva, contrast=c("Diet", "Vitamin.D.Deficient", "Low.Protein"))
    res = res[-which(is.na(res$pvalue)),]
    fdr <- fdrtool(res$stat, statistic="normal", plot=F)#, cutoff.method = "locfdr")
    res[,"pval_tool"] = fdr$pval
    res[,"padj_tool"]  <- p.adjust(res$pval_tool, method = "BH")
    res_diet[[paste0("rix",r)]]$Low.Protein = res
    
  }
  if(length(grep("Methyl.Enriched", colData_use$Diet)) > 0 & length(grep("Vitamin.D.Deficient", colData_use$Diet)) > 0){
    res = results(ddssva, contrast=c("Diet", "Vitamin.D.Deficient", "Methyl.Enriched"))
    res = res[-which(is.na(res$pvalue)),]
    fdr <- fdrtool(res$stat, statistic="normal", plot=F)#, cutoff.method = "locfdr")
    res[,"pval_tool"] = fdr$pval
    res[,"padj_tool"]  <- p.adjust(res$pval_tool, method = "BH")
    res_diet[[paste0("rix",r)]]$Methyl.Enriched = res
  }
}

### output objects: dds_string_sva, res_PO, res_diet


######## SOME PLOTS #########
pdf(file.path(dir,"figures_and_updates/pval_fdrtool_diets.pdf"))
par(mfrow=c(2,1))
res_diet = rem_seg_res$res_diet
for(r in names(res_diet)){
  for(d in names(res_diet[[r]])){
    res = res_diet[[r]][[d]]
    hist(res$pvalue, main = paste("P-values no adjustment",r,d))
    hist(res[,"pval_tool"], main = paste("P-values with estimated null model parameters",r,d), xlab="p", 
         sub=paste(paste(colnames(fdr$param), round(fdr$param,3), sep=":"), collapse=", "),
         cex.main=0.7, cex.sub=0.7)
  }
}
dev.off()


sig_PO = lapply(res_PO, function(x) x[which(x$padj_tool < 0.05 & x$baseMean  > 100),])
all_sig_PO_genes = table(unlist(lapply(sig_PO, rownames)))
all_sig_PO_genes[order(-all_sig_PO_genes)]


saveRDS(dds_string_sva, file.path(dir, "de_results/dds_lst_string_sva_11oct2019.rds"))
dds_string = readRDS(file.path(dir, "de_results/dds_lst_string-28sep2019.rds"))
res_diet$rix6$Standard$Gene = rownames(res_diet$rix6$Standard) 
sig_VDD_STD_rix6_string = as.data.frame(res_diet$rix6$Standard) %>% filter(padj < 0.01) %>%
  arrange(padj)

write.csv(sig_VDD_STD_rix6_string, file.path(dir,"figures_and_updates/string_VDD-STD_RIX6_0.01_27sep2019.csv"))

pdf(file.path(dir,"figures_and_updates/string_VDD-STD_RIX6_0.05_27sep2019.pdf"))
par(mfrow=c(2,1))
for(i in 1:10){
  p = plotCounts(dds_string$rix6, gene=sig_VDD_STD_rix6_string$Gene[i],intgroup = "Diet",returnData = F)
}
dev.off()

#################  DE with all RIX at once #######################
alldds = readRDS("C:/Users/Kathie/Dropbox (ValdarLab)/de_results/allrix_string_deseq_sva_noX_25oct2019.rds")
allres = results(alldds, contrast=c("PO_cat","a","b"))
allres = results(alldds, contrast=c("Diet","Vitamin.D.Deficient","Standard"))

allres = allres[-which(is.na(allres$pvalue)),]

fdr <- fdrtool(allres$stat, statistic="normal", plot=F)
allres$pval_raw = allres$pvalue
allres$padj_deseq2 = allres$padj
allres$praw_bh = p.adjust(allres$pval_raw, method = "BH")
allres[,"pvalue"] = fdr$pval
allres[,"padj"]  <- p.adjust(allres$pvalue, method = "BH")

sig_res = allres[which(allres$padj < 0.001),]
sig_res[order(sig_res$padj),]
plot_genes = rownames(allres[which(allres$praw_bh < 0.2 & allres$baseMean > 10),])

p = list()
for(i in 1:length(plot_genes)){
  gene = plot_genes[i]
  pd <- plotCounts(alldds, gene=gene, intgroup = "Diet",returnData = T)
  p[[i]] = ggplot(pd, aes(x=Diet, y=count)) +
    geom_boxplot() +# geom_point() +
    geom_jitter(width = 0.2, col="blue") + 
    ggtitle(paste("Normalized counts across RIXs in",gene)) +   
    theme_bw() + facet_grid()
  #points(x=jitter(as.numeric(as.factor(pd$PO))), y=pd$count)  
  #anova(lm(data=pd, count ~ PO))
  #if(i%%4 == 0){
  #  quad = grid.arrange(p[[i-3]], p[[i-2]], p[[i-1]], p[[i]], ncol=2, nrow=2)
  #  print(quad)
  #}
}
#data.frame(sig_res) %>% 


####################  plot signif genes  ##########################

source("manhattan_plot.R")
ie_genes = read.csv("C:/Users/Kathie/Dropbox (ValdarLab)/imprinted_genes/2014_05_14_allImprintedGenes.csv")
icr = read.table(file.path(dir, "data/imprinted_genes_list/mousebook_imprinted_regions_only_11nov2019.txt"), sep="\n")
icr = lapply(icr$V1, function(x) gsub(" ", "", unlist(strsplit(x, ","))))
icr_list = lapply(icr, function(x) {
  tmp = left_join(data.frame(Gene.Name = x), annot_genes, by="Gene.Name")
  tmp$Start = as.numeric(tmp$Start)
  tmp$End = as.numeric(tmp$End)
  minPos = min(tmp$Start, na.rm=T)
  maxPos = max(tmp$End, na.rm=T)

  if(maxPos - minPos < 5000000 & nrow(tmp) > 1){
    tmp2 = distinct(data.frame(Gene.Name = paste0(tmp$Gene.Name[which(tmp$Start == minPos)], "_",
                                        tmp$Gene.Name[which(tmp$End == maxPos)]),
                      Chr = tmp$Chr,
                      Start = min(tmp$Start,na.rm=T),
                      End = max(tmp$End,na.rm=T),
                      Strand = " ",
                      Gene.ID = " "))
    tmp = tmp2
  } 
  tmp
})

icr = do.call("rbind", icr_list)
icr %>% filter(!is.na(Chr)) %>% distinct() %>%
  arrange(Chr, Start) -> icr
icr[1,] = c("Ddc_Cobl", 11, 11814101, 12464960, "","")
icr[4,] = c("Commd1_Zrsr1", 11, 22896136, 22976496, "","")
icr[14,] = c("Mest_Klf14", 6, 30723547, 30959078, "","")

icr = icr[-c(2,3,5, 15:17),]


haplofiles = list.files("C:/Users/Kathie/Dropbox (ValdarLab)/mini/pup_haplo_blocks_by_CC_parent_dec2019/", 
                        pattern="haploBlocks.rds", full.names = T)
phased_par_CC_haplotypes = lapply(haplofiles, readRDS)
names(phased_par_CC_haplotypes) = do.call("rbind", strsplit(haplofiles, "_"))[,8]
CC_lab = read.csv(file.path(dir, "variant_data/matched_v2_4jun2018.csv"))
CC_lab = CC_lab %>% dplyr::select(c("Pup.ID","RIX","CC.1","CC.2")) %>% filter(!is.na(RIX))

pars = list()
for(r in unique(CC_lab$RIX)[order(unique(CC_lab$RIX))]){
  ppars = list()
  for(p in CC_lab$Pup.ID[which(CC_lab$RIX == r)]){ 
    if(p %in% names(phased_par_CC_haplotypes)){
      par1 = phased_par_CC_haplotypes[[paste(p)]]$founder_by_parent[[1]] %>% filter(nchar(found) > 1)
      par2 = phased_par_CC_haplotypes[[paste(p)]]$founder_by_parent[[2]] %>% filter(nchar(found) > 1)
      allpars = rbind(par1, par2)
      i=1
      repeat {
        if(allpars$chr[i] == allpars$chr[i+1] &
           allpars$group[i+1] == allpars$group[i]+1 & 
           allpars$start[i+1]-allpars$end[i]<50000){
          found = paste(unique(unlist(strsplit(allpars$found[i], ",")),
                               unlist(strsplit(allpars$found[i+1], ","))), collapse=",")
          if(nchar(found) == 3){
            allpars$end[i]=allpars$end[i+1]
            allpars$group[i]=allpars$group[i+1]
            allpars = allpars[-(i+1),]
          } else {
            i=i+1
          }
        } else {
          i=i+1
        }
        if(i == nrow(allpars)-1){
          break
        }
      }
      
      ppars[[p]] = allpars %>% arrange(chr, start) %>% ungroup() %>%
        dplyr::select(-c("group", "cc")) %>% distinct()
    }
  }
  allpars = do.call("rbind", ppars)
  allpars = allpars %>% arrange(chr, start) %>% dplyr::select(-"pup") %>% distinct()
  
  i=1
  repeat {
    if(allpars$chr[i] == allpars$chr[i+1] & allpars$found[i+1] == allpars$found[i] & 
       (allpars$start[i+1] == allpars$start[i] | allpars$end[i+1] == allpars$end[i])){
      allpars$end[i]=max(allpars$end[i], allpars$end[i+1])
      allpars$start[i]=min(allpars$start[i], allpars$start[i+1])
      allpars = allpars[-(i+1),]
    } else {
      i=i+1
    }
    if(i == nrow(allpars)-1){
      break
    }
  }
    
  pars[[paste(r)]] = unique(allpars)
  
}
pars = regions_list
for(i in 1:length(pars)){
  colnames(pars[[i]]) = tolower(colnames(pars[[i]]))
}
#mp = list()
for(r in unique(CC_lab$RIX)[order(unique(CC_lab$RIX))]){
  r=paste(r) 
  mat = data.frame(rem_seg_res$res_PO[[r]])
  mat$Gene.Name = rownames(rem_seg_res$res_PO[[r]])
  mat %>% left_join(annot_genes, by="Gene.Name") -> mat
  mat$Start = as.numeric(mat$Start)
  mat$End = as.numeric(mat$End)
  mat$pos = (mat$Start + mat$End) / 2
  mat$Chr = factor(mat$Chr, levels=c(1:19,"X"))
  ann<-rep(1, length(mat$pvalue))
  for(i in 1:nrow(icr)){
    tmp = icr[i,]
    flag = which(mat$Chr==tmp$chr & mat$pos>=tmp$start & mat$pos<=tmp$end)
    if(length(flag) > 0) ann[flag] = 2
  }
  for(i in 1:nrow(pars[[r]])){
    tmp = pars[[r]][i,]
    flag = which(mat$Chr==tmp$chr & mat$Start>=tmp$start & mat$End<=tmp$end)
    if(length(flag) > 0) ann[flag] = 3
  }
  
  ##compare_ie from deseq2_stringtie.R
  compare_ie = read.csv(file.path(dir,"imprinted_genes/use_ie_gene_names.csv"))
  agree = intersect(compare_ie$x, mat$Gene.Name)
  ie_genes %>% filter(x %in% mat$Gene.Name) %>% dplyr::select("x") -> also_agree
  agree = unique(c(agree, also_agree$mgi_symbol))
  ann[which(mat$Gene.Name %in% agree)]<-2
  ann<-factor(ann, levels=1:3, labels=c("","ICR","SEG"))
  CCs = paste0(CC_lab$CC.1[which(CC_lab$RIX == r)[1]], "/", CC_lab$CC.2[which(CC_lab$RIX == r)[1]])
  png(paste0(dir,"/figures_and_updates/allRIX_manhattan_RIX",r,"_24nov2019.png"),width=1600,height=400,units="px")

  genome_wide_plot(chr=mat$Chr, pos=mat$pos, pvalue=mat$pvalue,
                   annotate=ann, regions=regions_list[[r]], chr_lengths=NULL, 
                   main=paste("Raw p-values in RIX", r))
  dev.off()
    #manhattan.plot(mat$Chr, mat$pos, mat$pval_raw, 
    #                       annotate = list(ann, "ICR"=list(label=list(show=F), col="forestgreen"),
    #                                       "SEG"=list(label=list(show=F), col="dodgerblue")),
    #                       sig.level = 2e-5, 
    #                       should.thin = F, side="bottom", xlab=paste("RIX",r,CCs))
}
#pdf(file.path(dir, "figures_and_updates/allRIX_manhattan_18nov2019.pdf"))
#  apply(mp, function(x) 
#  print(x))


############  meta-analysis  #############

sample_sizes = colData %>% group_by(RRIX) %>% tally()
sample_sizes$rix = paste0("rix",sample_sizes$RRIX)

#dds_string_compare$res_PO
#dds_string_compare$res_diet
rma_comb = sapply(names(rem_seg_res$res_PO), 
                  function(x){
                    #if(length(rownames(dds_string_compare$res_diet[[x]]$Standard)) > 0){
                      cbind(genes = rownames(rem_seg_res$res_PO[[x]]), 
                                    rix = x,
                                    rem_seg_res$res_PO[[x]][,c("log2FoldChange","lfcSE","stat")])
                    #}
                    })
rma_test = data.frame(do.call("rbind",rma_comb))
rma_mesh = rma_test[,-5]
rma_mesh$log2FoldChange = abs(rma_mesh$log2FoldChange)
#write.table(rma_mesh, file.path(dir, "de_results/MeSH_beta_values.txt"), col.names = F, row.names = F, quote = F)
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
  
  #mean(xed)
  #var(xed)
  #hist(zed, prob=T, 30, ylim = c(0,0.5))
  #y = dnorm(seq(-5,5,by=0.1), mean(zed), sd(zed))
  #lines(seq(-5,5,by=0.1), y, type="l", col="blue")
  #fdr <- fdrtool(test$stat, statistic="normal", plot=T)

  #outmat_relm_abs$pval_raw = outmat_relm_abs$pval
  #outmat_relm_abs$pval_b = fdr$pval
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
#if(length(grep("Gm|Rik|[.]", sum_mat$genes)) > 0){
#  sum_mat_clean = sum_mat[-grep("Gm|Rik|[.]", sum_mat$genes),]  %>% arrange(padj_s)
#} else {
  sum_mat_clean = sum_mat %>% arrange(pval_s)
#}

sum_mat_clean = sum_mat_clean %>% left_join(baseMean_comb[,c("genes","means")], by="genes")


plot_genes = sum_mat_clean$genes[which(sum_mat_clean$padj_s < 0.05)]
## 127 genes (47 w/o Gm|Rik|[.])

plot_genes = outmat_relm_abs_sig$genes
plot_genes = sum_mat_clean$genes[1:10]
plot_genes = pmat_sig$g[1:10]
#n=50 #7,8,9,12,40,36,50
diet_leg = data.frame(from=levels(samples_use$Diet), to=c("st","lp","me","vd"))
pdf(file.path(dir,"figures_and_updates/vdd_std_select_genes_30dec2019.pdf"))
forest=F
for(n in plot_genes){
  if(forest){
    f = forest(meta_out_reml_abs[[n]],
               main=paste("Meta-analysis across RIX for",n),
               slab=toupper(levels(rma_test$rix)))
  }
  
  p=list()
  for(r in names(dds_lst)){
    if(n %in% rownames(dds_lst[[r]])){
      pd = plotCounts(dds_lst[[r]], gene=n, intgroup = "Diet",returnData = T)
      pd$Diet = diet_leg$to[match(pd$Diet, diet_leg$from)]
      pval = ifelse(n %in% rownames(rem_seg_res$res_diet[[r]]$Standard), 
                    round(rem_seg_res$res_diet[[r]]$Standard$padj[which(rownames(rem_seg_res$res_diet[[r]]$Standard) == n)],4), 
                    "NA")
      p[[r]] = ggplot(pd, aes(x=Diet, y=count)) +
        geom_boxplot(outlier.shape = NA) +# geom_point() +
        geom_jitter(width = 0.2, col="blue") + 
        ggtitle(label = paste("RIX", toupper(r)), 
                subtitle = paste("st vs vd p-val:", pval)) + 
        theme_classic() + facet_grid() + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    } else {
      p[[r]] = ggplot(data.frame(NA))
    }
  }
  quad = grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],
                      p[[6]],p[[7]],p[[8]],p[[9]], ncol=3, nrow=3,
                      top=paste("Normalized counts of",n))
  
  print(quad)
  #print(p)

}
dev.off()
rma_test %>% filter(genes == "Trim14")
rma_test %>% filter(genes == "mt-Tg")

###########  metafor  ##############

#meta_PO_allgenes 
meta_out_reml_abs <- lapply(all_genes, function(g){
  #})
  #meta_out = list()
  #for(g in all_genes){
  gene = rma_test %>% filter(genes == g) %>% arrange(rix)
  if(length(which(!is.na(gene$log2FoldChange))) > 0){
    rma(yi = abs(gene$log2FoldChange),vi=gene$lfcSE, 
        method="REML", verbose=F, control=list(stepadj=0.5,maxiter=1000))
  } else {
    NA
  }
  #meta_out
})
#names(meta_out) = names(meta_out_reml) = 
names(meta_out_reml_abs) = all_genes
perm_test_all = lapply(meta_out_reml, function(x) permutest(x))

outmat_relm_abs <- do.call("rbind", lapply(meta_out_reml_abs, function(x) x[c("b","se","pval","ci.lb","ci.ub","tau2","se.tau2")]))
outmat_relm_abs = as.data.frame(outmat_relm_abs)
outmat_relm_abs$genes = all_genes
outmat_relm_abs = outmat_relm_abs %>% left_join(baseMean_comb[,c("genes","means")], by="genes") 


fdr <- fdrtool(unlist(outmat_relm_abs$pval), statistic="pvalue", plot=T)
#fdr_b <- fdrtool(unlist(outmat_relm_abs$b), statistic="normal", plot=T)
outmat_relm_abs$pval_raw = outmat_relm_abs$pval
#outmat_relm_abs$pval_b = fdr$pval
outmat_relm_abs$pval = fdr$pval

#outmat_relm_abs$padj_b  <- p.adjust(outmat_relm_abs$pval, method = "BH")
outmat_relm_abs$padj  <- p.adjust(outmat_relm_abs$pval, method = "BH")

outmat_relm_abs %>% filter(genes %in% keep_genes, 
                           !is.na(b), !is.na(se.tau2), !is.na(means),
                           padj < 0.1) %>%
  mutate(pval = as.numeric(pval)) %>% 
  arrange(pval) -> outmat_relm_abs_sig
#if(length(grep("Gm|Rik", outmat_relm_abs_sig$genes)) > 0){
#  outmat_relm_abs_sig = outmat_relm_abs_sig[-grep("Gm|Rik", outmat_relm_abs_sig$genes),]
#}
head(outmat_relm_abs_sig)


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
pcomb$rix = factor(pcomb$rix, levels=paste0("rix",c(1:4,6:10)))
all_genes = unique(pcomb$genes)
length(which(all_genes %in% keep_genes))
pmat = lapply(keep_genes, function(g){
  test = pcomb %>% filter(genes == g) 
  stat_raw = -2*sum(log(test$praw))
  stat_adj = -2*sum(log(test$pval))
  data.frame(stat_raw, stat_adj, g)
})

pmat=do.call("rbind",pmat)
hist(pmat$stat_raw, breaks=96)
#y = dchisq(x=seq(0,200,by=0.1), df = 2*9)
#plot(dchisq(x=seq(0,200,by=0.1), df = 18), type="l", col="blue")
curve(dchisq(x, df = length(pcomb)*2)*nrow(pmat), from=0, to=150, col="blue", add=T)

pmat$padj = pchisq(pmat$stat_adj, df=length(pcomb)*2,lower.tail = F)
pmat$praw = pchisq(pmat$stat_raw, df=length(pcomb)*2,lower.tail = F)
fdr = fdrtool(x = pmat$padj, statistic = "pvalue")
pmat$padj_qval = fdr$qval
pmat$padj_fdr = fdr$pval
fdr = fdrtool(x = pmat$praw, statistic = "pvalue")
pmat$praw_qval = fdr$qval
pmat$praw_fdr = fdr$pval
pmat$padj_bh = p.adjust(pmat$padj_fdr, method = "BH")
pmat$praw_bh = p.adjust(pmat$praw_fdr, method = "BH")

## compare with sum_mat_clean
pmat_sig = pmat %>% filter(g %in% all_genes, praw_bh < 0.05) %>% arrange(padj_qval) 
## 1697, 1140 w/o Gm|Rik|[.]
plot_genes = pmat_sig$g[-grep("Gm|[.]|Rik",pmat_sig$g)][1:20]
which(pmat_sig$g == "Gprasp1")
#write.csv(pmat_sig, file.path(dir, "de_results/combined_VDD_STD_sigif_genes_14nov2019.csv"))


##############  stouffer  #################

zstar = zedst = c()
for(g in keep_genes){
  tg = test_mat %>% filter(genes == g)
  if(nrow(tg) == 9){
    tg$weights = sqrt(tg$n)
    tg$weights = 1 / tg$lfcSE
    zstar[g] = sum(apply(tg,1,function(x) abs(as.numeric(x["weights"])*as.numeric(x["stat"])))) / sqrt(sum(tg$weights^2))
    zedst[g] = sum(apply(tg,1,function(x) abs(as.numeric(x["weights"])*as.numeric(x["zed"])))) / sqrt(sum(tg$weights^2))
  }
}

hist(zstar)
lines(dnorm(x=seq(-5,5,by=0.1),0,1), type="l", col="blue")
hist(zedst)
lines(dnorm(x=seq(-5,5,by=0.1),0,1), type="l", col="blue")
fzstar = folded(zstar)

folded <- function(y) {
  
  ## y is a vector with positive data
  n <- length(y)  ## sample size
  sy2 <- sum(y^2)
  
  sam <- function(para, n, sy2) {
    me <- para[1]   ;   se <- exp( para[2] )
    f <-  - n/2 * log(2/pi/se) + n * me^2 / 2 / se +
      sy2 / 2 / se - sum( log( cosh( me * y/se ) ) )
    f
  }
  
  mod <- optim( c( mean(y), sd(y) ), n = n, sy2 = sy2, sam, control = list(maxit = 2000) )
  mod <- optim( mod$par, sam, n = n, sy2 = sy2, control = list(maxit = 20000) )
  result <- c( -mod$value, mod$par[1], exp(mod$par[2]) )
  names(result) <- c("log-likelihood", "mu", "sigma squared")
  result
  
}



###################  permute  ####################
pups = gsub("Pup.ID_","",colnames(countData))
rixes = unique(samples_use$RRIX)[order(unique(samples_use$RRIX))]
diff_lst <- list()
for(r in rixes){
  use_pups = pups[which(pups %in% samples_use$Pup.ID[which(samples_use$RRIX == r)])]
  pups_a = use_pups[which(use_pups %in% samples_use$Pup.ID[which(samples_use$PO_cat == "a")])]
  pups_b = use_pups[which(use_pups %in% samples_use$Pup.ID[which(samples_use$PO_cat == "b")])]
  pups_a = paste0("Pup.ID_", pups_a)
  pups_b = paste0("Pup.ID_", pups_b)
  N = min(length(pups_a), length(pups_b))
  rand = cbind(sample(1:length(pups_a), N, replace=F), sample(1:length(pups_b), N, replace=F))
  
  diff <- sapply(1:N, function(n){
    countData[,pups_a[rand[n,1]]] - countData[,pups_b[rand[n,2]]]
  })
  colnames(diff) = paste0("rix", r,"_",1:N)
  diff_lst[[paste0("rix",r)]] = diff
}


rix_means <- do.call("cbind", lapply(diff_lst, function(x) apply(x, 1, mean)^2))
rix_vars <- do.call("cbind", lapply(diff_lst, function(x) apply(x, 1, var)^2))
ignore <- apply(rix_vars, 1, function(x) ifelse(length(which(x == 0)) > 0, T, F))

rix_weights <- do.call("cbind", lapply(diff_lst, function(x) apply(x, 1, length)))
rix_weights <- t(apply(rix_weights, 1, function(x) x/sum(x)))
rix_means <- rix_means[-which(ignore),]
rix_vars <- rix_vars[-which(ignore),]
rix_weights <- rix_weights[-which(ignore),]
rix_meta <- sapply(1:nrow(rix_means), function(x){
  data = cbind(means = as.numeric(rix_means[x,]), rix = colnames(rix_means))
  #lm(~ data[,"rix"], data=data.frame(data), weights=rix_weights[x,])
  rma(yi = rix_means[x,], vi=rix_vars[x,], weights = rix_weights[x,],
      method="EB", verbose=F, control=list(stepadj=0.5,maxiter=1000))
}, simplify=F)
outmat_rix_meta <- do.call("rbind", lapply(rix_meta, function(x) x[c("b","se","pval","ci.lb","ci.ub","tau2","se.tau2")]))
outmat_rix_meta = as.data.frame(outmat_rix_meta)
outmat_rix_meta$gene = rownames(rix_means)

outmat_rix_meta %>% filter(pval < 0.05)

diff_set <- t(apply(diff_set, 1, function(x) x + abs(min(x))+1))
colData = cbind(samples=colnames(diff_set), rix=unlist(strsplit(colnames(diff_set),"_"))[c(T,F)])


dds = DESeqDataSetFromMatrix(countData = diff_set, colData = colData,
                             design = ~ rix)
DESeq(dds)


###############   PCA   ################
cleanY = function(y, mod, svs) {
  X = cbind(mod, svs)
  Hat = solve(t(X) %*% X) %*% t(X)
  beta = (Hat %*% t(y))
  rm(Hat)
  gc()
  P = ncol(mod)
  return(y - t(as.matrix(X[,-c(1:P)]) %*% beta[-c(1:P),]))
}

cleaned_RIX_1 <- cleanY(countData_use, mod, svseq$sv)
colData_use$SV1 = svseq$sv[,1]
### nvm this was a bad idea

vsd_string <- lapply(dds_string, function(x) vst(x, blind=FALSE))

pcaData = lapply(vsd_string, function(x) plotPCA(x, intgroup=c("PO","Diet"), returnData=T))
pdf(file.path(dir,"figures_and_updates/string_PCA_by_rix_4outs_rem.pdf"))
for(i in 1:length(pcaData)){
  percentVar <- round(100 * attr(pcaData[[i]], "percentVar"))
  pcaData[[i]]$PO = as.factor(pcaData[[i]]$PO)
  PC_means = pcaData[[i]] %>% group_by(Diet, PO) %>% summarize(med1=median(PC1), med2=median(PC2))
  pcaData[[i]]$far = sapply(1:nrow(pcaData[[i]]), function(d){
    x = ifelse(abs(pcaData[[i]]$PC1[d] - 
                 PC_means$med1[which(PC_means$Diet == pcaData[[i]]$Diet[d] & PC_means$PO == pcaData[[i]]$PO[d])]) >
      sd(PC_means$med1)*5, T, F)
    y = ifelse(abs(pcaData[[i]]$PC2[d] - 
                 PC_means$med2[which(PC_means$Diet == pcaData[[i]]$Diet[d] & PC_means$PO == pcaData[[i]]$PO[d])]) >
                 sd(PC_means$med2)*5, T, F)
    any(x,y)
  })
  p = ggplot(pcaData[[i]], aes(PC1, PC2, color=Diet, shape=PO)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    coord_fixed() + 
    theme_bw() + 
    geom_label_repel(aes(label=ifelse(far,name,'')),
                     #(x>-0.15 & x<0 & y>-0.2 & y<0.2 & (RIX %in% c(2,3)))
                     box.padding   = 0.35, 
                     point.padding = 0.5,
                     segment.color = 'grey50') + 
    ggtitle(paste("PCA of",names(pcaData)[i]))
  print(p)
  
}
dev.off()

###############   SVA   ################

library(sva)

dat_lst = lapply(dds_string, function(x) counts(x, normalized=T))
dat_lst  <- lapply(dat_lst, function(x) x[rowMeans(x) > 1,])
mod  <- lapply(dds_string, function(x) model.matrix(~ Diet + PO_cat, colData(x)))
mod0 <- lapply(dds_string, function(x) model.matrix(~   1, colData(x)))
svseq <- sapply(1:length(dat_lst), function(x) svaseq(dat_lst[[x]], mod[[x]], mod0[[x]], n.sv = 2), simplify=F)              
names(svseq) = names(dat_lst)

#ddssva <- dds_string
ressva_PO=list()
for(r in names(ddssva)){
  ddssva[[r]]$SV1 <- svseq[[r]]$sv[,1]
  ddssva[[r]]$SV2 <- svseq[[r]]$sv[,2]
  
  design(ddssva[[r]]) <- ~ SV1 + Diet + PO_cat 
  ddssva[[r]] = DESeq(ddssva[[r]])
  ressva_PO[[r]] = results(ddssva[[r]], contrast=c("PO_cat","a","b"))
}
pdf(file.path(dir,"figures_and_updates/raw_pvalues_string_SV1_PO.pdf"))
sapply(1:length(rem_seg_res$res_PO), function(x) 
  hist(rem_seg_res$res_PO[[x]]$pval_raw, main=paste("Histogram of p-values in", names(rem_seg_res$res_PO)[x])))
dev.off()

pdf(file.path(dir,"figures_and_updates/pval_fdrtool.pdf"))
for(r in names(res_PO)){
  #[[paste0("rix",r)]] <- res))
  res = res_PO[[r]]
  res = res[-which(is.na(res$pvalue)),]
  fdr <- fdrtool(res$stat, statistic="normal", plot=F)#, cutoff.method = "locfdr")
  #fdr <- fdrtool(res$pvalue, statistic="pvalue")#, cutoff.method = "locfdr")
  
  res[,"pval_tool"] = fdr$pval
  res[,"padj_tool"]  <- p.adjust(res$pval_tool, method = "BH")
  
  par(mfrow=c(2,1))
  hist(res[,"pval_tool"], main = paste("P-values with estimated null model parameters",r), xlab="p", 
       sub=paste(paste(colnames(fdr$param), round(fdr$param,3), sep=":"), collapse=", "),
       cex.main=0.9, cex.sub=0.7)
  hist(res[,"padj_tool"], main = paste("Adj P-values with estimated null model parameters",r), xlab="p", 
       sub=paste(paste(colnames(fdr$param), round(fdr$param,3), sep=":"), collapse=", "),
       cex.main=0.9, cex.sub=0.7)
  
  length(which(res$padj_tool < 0.05))
  res_PO[[r]] = res
}
dev.off()

min(fdr$qval)
min(fdr$lfdr)
lapply(fdr, head)
cut = fndr.cutoff(res$pvalue, statistic="pvalue")
cut=seq(0, 1, by=0.05)
lapply(cut, function(x) censored.fit(res$pvalue, x, statistic="pvalue"))
hist(fdr$pval)
plotMA(res)
