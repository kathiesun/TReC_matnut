
run_DESeq2 <- function(colData, countData, sva=T, interact=T, contrast=c("PO", "Diet"), 
                       rix=NULL, annot=NULL, lowCount = 10, regions_list=NULL){
  if(length(contrast) > 1) contrast = contrast[1]
  dds_lst = list()
  if (is.null(rix)) rix = unique(colData$RRIX[order(colData$RRIX)])
  for(r in rix){
    colData_use = data.frame(colData %>% filter(RRIX==r))
    rownames(colData_use) = paste0("Pup.ID_", colData_use$Pup.ID)
    countData_use = countData[,which(colnames(countData) %in% paste0("Pup.ID_", colData_use$Pup.ID))]
    #print(identical(rownames(colData_use), colnames(countData_use)))

    if(!is.null(regions_list)){
      annot$Start = as.numeric(annot$Start)
      annot$End = as.numeric(annot$End)
      
      regions = data.frame(regions_list[[paste(r)]])
      #remove_genes = apply(regions, 1, function(x) 
      #  intersect(which(annot$Chr %in% x["Chr"]), 
      #            intersect(which(as.vector(annot$Start) >= as.numeric(x["start"])),
      #                  which(as.vector(annot$End) <= as.numeric(x["end"])))))
      remove_genes = c()
      for(i in 1:nrow(regions)){
        remove_genes = c(remove_genes,
                         which(annot$Chr == regions$Chr[i] & 
                                 annot$Start >= regions$start[i] & 
                                 annot$End <= regions$end[i]))
      }
      #remove_genes = apply(regions, 1, function(x) 
      #  intersect(which(annot$Chr %in% x["Chr"]), 
      #            intersect(which(as.vector(annot$Start) >= as.numeric(x["start"])),
      #                      which(as.vector(annot$End) <= as.numeric(x["end"])))))
      
      remove_genes = unique(unlist(remove_genes))
      remove_names = annot_genes$Gene.Name[remove_genes]
      countData_use = countData_use[-which(rownames(countData_use) %in% remove_names),]
    }
    if(contrast == "PO") {
      countData_use = countData_use[-which(rownames(countData_use) %in% annot$Gene.Name[which(annot$Chr == "X")]),]
    }
    if(!all(rownames(colData_use) == colnames(countData_use))) warning("Pups out of order") 
    if(interact){
      base = "Diet + PO + Diet:PO"
    } else {
      base = "Diet + PO_cat"
    }
    dds <- DESeqDataSetFromMatrix(countData = countData_use,
                                  colData = colData_use,
                                  design = as.formula(paste("~", base)))
    keep <- rowSums(counts(dds) >= 10) >= 10
    dds <- dds[keep,]
    dds <- estimateSizeFactors(dds)
    dds <- estimateDispersions(dds)
    #ddsPO <- nbinomLRT(dds, reduced = ~ Diet)
    #dds = DESeq(dds)
    if(sva){
      dat = counts(dds, normalized=T)
      dat = dat[rowMeans(dat) > 1,]
      mod  <- model.matrix(as.formula(paste("~", base)), colData(dds))
      mod0 <- model.matrix(~   1             , colData(dds))
      svseq <- svaseq(dat, mod, mod0, n.sv = 1)  
      ddssva = dds
      ddssva$SV1 <- svseq$sv[,1]
      
      design(ddssva) <- as.formula(paste("~ SV1 +", base))
      dds = ddssva
    }
    #ddssva = DESeq(ddssva)
    #ddssva = nbinomWaldTest(ddssva)
    if(interact){
      reduced = paste(ifelse(sva, "~ SV1 +", "~"), "Diet + PO")  #ifelse(contrast == "PO", "Diet", "PO"))
      dds <- DESeq(dds, test="LRT", reduced=as.formula(reduced))
    } else {
      dds = DESeq(dds)
      dds = nbinomWaldTest(dds)
    }
    
    dds_lst[[tolower(r)]] = dds
    
    print(paste("finished deseq on ",r))
  }
  return(dds_lst)
}



run_DESeq2_old <- function(colData, countData, sva=T, interact=T, rix=NULL){
  dds_lst = list()
  if (is.null(rix)) rix = unique(colData$RRIX[order(colData$RRIX)])
  for(r in rix){
    colData_use = data.frame(colData %>% filter(RRIX==r))
    countData_use = countData[,which(colnames(countData) %in% paste0("Pup.ID_", colData_use$Pup.ID))]
    if(!all(rownames(colData) == colnames(countData))) warning("Pups out of order") 
    if(interact){
      base = "Diet + PO + Diet:PO"
    } else {
      base = "Diet + PO_cat"
    }
    dds <- DESeqDataSetFromMatrix(countData = countData_use,
                                  colData = colData_use,
                                  design = as.formula(paste("~", base)))
    dds = DESeq(dds)
    if(sva){
      dat = counts(dds, normalized=T)
      dat = dat[rowMeans(dat) > 1,]
      mod  <- model.matrix(as.formula(paste("~", base)), colData(dds))
      mod0 <- model.matrix(~   1             , colData(dds))
      svseq <- svaseq(dat, mod, mod0, n.sv = 1)  
      ddssva = dds
      ddssva$SV1 <- svseq$sv[,1]
      
      design(ddssva) <- as.formula(paste("~ SV1 +", base))
      ddssva = DESeq(ddssva)
      ddssva = nbinomWaldTest(ddssva)
      dds = ddssva
    }
    
    dds_lst[[tolower(r)]] = dds
    
    print(paste("finished deseq on ",r))
  }
  return(dds_lst)
}

#####################  Analyze any set of dds_lst result  ###########################

compare_diets_PO <- function(dds_lst, ref_diet="Standard", fdrtool_adjust=T, adjust_p=F){
  res_PO = res_diet = list()
  for(r in names(dds_lst)){
    ddssva = dds_lst[[paste(r)]]
    if(length(grep(":", design(ddssva))) > 0) {
      res = results(ddssva)
      tmp = gsub("\'|\"", "", unlist(strsplit(res@elementMetadata@listData$description[4]," ")))
      contrast = ifelse(tmp[length(tmp)] == "Diet", "PO", "Diet")
      if(fdrtool_adjust){
        if(length(which(is.na(res$pvalue))) > 0) res = res[-which(is.na(res$pvalue)),]
        #if(adjust_p){
          fdr <- fdrtool(res$pvalue, statistic="pvalue", plot=F)
          res$pval_raw = res$pvalue
          res$padj_deseq2 = res$padj
          res$pvalue = fdr$pval
          res$padj  <- fdr$qval
        #} else{
        #  fdr <- fdrtool(res$stat, statistic="normal", plot=F)
        #  res$pval_raw = res$pvalue
        #  res$padj_deseq2 = res$padj
        #  res$pvalue = fdr$pval
        #  res$padj  <- p.adjust(res$pvalue, method = "BH")
        #}
      }
      if(contrast == "PO"){
        res_PO[[paste(r)]] = res
      } else {
        res_diet[[paste(r)]] = NULL
      }
    } else {
      if(length(grep("PO", resultsNames(ddssva)))>0){
        if(length(grep("PO_cat", resultsNames(ddssva)))>0){
          res <- results(ddssva, contrast = c("PO_cat","a","b"))
        } else {
          res <- results(ddssva, name = "PO")
        }
        if(fdrtool_adjust){
          if(length(which(is.na(res$pvalue))) > 0) res = res[-which(is.na(res$pvalue)),]
          if(adjust_p){
            fdr <- fdrtool(res$pvalue, statistic="pvalue", plot=F)
            res$pval_raw = res$pvalue
            res$padj_deseq2 = res$padj
            res$pvalue = fdr$pval
            res$padj  <- fdr$qval
          } else{
            fdr <- fdrtool(res$stat, statistic="normal", plot=F)
            res$pval_raw = res$pvalue
            res$padj_deseq2 = res$padj
            res$pvalue = fdr$pval
            res$padj  <- p.adjust(res$pvalue, method = "BH")
          }
        }
        res_PO[[paste(r)]] <- res
      } else {
        res_PO[[paste(r)]] <- NA
      }
      
      res_diet[[paste(r)]] <- list()
      
      diets <- gsub(" ",".",as.character(unique(ddssva$Diet)))
      other_diets = diets[-which(diets == ref_diet)]
      for(d in other_diets){
        res_lst_diet = list()
        if(length(grep(ref_diet, resultsNames(ddssva)) ) > 0 & length(grep(d, resultsNames(ddssva))) > 0){
          res = results(ddssva, contrast=c("Diet", ref_diet, d))
          if(fdrtool_adjust){
            if(length(which(is.na(res$pvalue))) > 0) res = res[-which(is.na(res$pvalue)),]
            if(adjust_p){
              fdr <- fdrtool(res$pvalue, statistic="pvalue", plot=F)
              res$pval_raw = res$pvalue
              res$padj_deseq2 = res$padj
              res$pvalue = fdr$pval
              res$padj  <- fdr$qval
            } else{
              fdr <- fdrtool(res$stat, statistic="normal", plot=F)
              res$pval_raw = res$pvalue
              res$padj_deseq2 = res$padj
              res$pvalue = fdr$pval
              res$padj  <- p.adjust(res$pvalue, method = "BH")
            }
          }
          res_diet[[paste(r)]][[d]] = res
        }
      }
    }
    
  }
    
  return(list(res_PO = res_PO, res_diet = res_diet))
}


check_sig_PO_SDP <- function(PO_results_lst, 
                             sequence_details, CC_labels, sample_data,
                             phased_CC_haplotypes, alpha = 0.05){
  dict = data.frame(num = 1:8, let=LETTERS[1:8])
  sig_PO <- sapply(1:length(PO_results_lst), function(x){
    if(length(which(PO_results_lst[[x]][,"padj"] < alpha)) > 0)
      cbind(rownames(PO_results_lst[[x]])[which(PO_results_lst[[x]][,"padj"] < alpha)], names(PO_results_lst)[x])
  })
  sig_PO <- do.call("rbind", sig_PO)
  sig_PO = as.data.frame(sig_PO)
  colnames(sig_PO) = c("gene", "RIX")
  sig_PO %>% group_by(gene) %>%
    mutate(count = n()) %>% arrange(-count, gene) %>%
    left_join(sequence_details, by="gene") -> sig_PO
  
  CC_consensus_phases = list()
  rixes = unique(as.numeric(sample_data$RRIX))
  for(r in rixes[order(rixes)]){
    #r=unique(samples_use$RRIX[order(samples_use$RRIX)])[i]
    CCs = unique(c(CC_labels$CC.1[which(CC_labels$RIX == r)],CC_labels$CC.2[which(CC_labels$RIX == r)]))
    CCs = CCs[order(CCs)]
    pups = names(phased_CC_haplotypes)[which(names(phased_CC_haplotypes) %in% CC_labels$Pup.ID[which(CC_labels$RIX == r)])]
    check_ccs = as.vector(unlist(unique(lapply(pups, function(x) lapply(phased_CC_haplotypes[[x]]$founder_by_parent, function(y) unique(y$cc))))))
    identical(CCs, check_ccs)
    all_par1 = do.call("rbind", lapply(pups, function(x) phased_CC_haplotypes[[x]]$founder_by_parent[[1]]))
    all_par2 = do.call("rbind", lapply(pups, function(x) phased_CC_haplotypes[[x]]$founder_by_parent[[2]]))
    all_par1 %>% arrange(chr, start) %>%
      ungroup() %>%
      dplyr::select(-c("pup","group")) %>%
      distinct() %>%
      mutate(rix = r, chr=as.character(chr)) -> all_par1
    all_par2 %>% arrange(chr, start) %>%
      ungroup() %>%
      dplyr::select(-c("pup","group")) %>%
      distinct() %>%
      mutate(rix = r, chr=as.character(chr)) -> all_par2
    CC_consensus_phases[[unique(all_par1$cc)]] = all_par1
    CC_consensus_phases[[unique(all_par2$cc)]] = all_par2
  }
  sig_PO$start = as.numeric(sig_PO$start)
  sig_PO$end = as.numeric(sig_PO$end)
  
  founders = apply(sig_PO, 1, function(x){
    rix = as.numeric(gsub("rix", "", x["RIX"]))
    CC.1 = unique(CC_labels$CC.1[which(CC_labels$RIX == rix)])
    CC.2 = unique(CC_labels$CC.2[which(CC_labels$RIX == rix)])
    CC_consensus_phases[[CC.1]] %>%
      filter(chr %in% x["chr"]) %>% 
      filter(start < as.numeric(x["start"]), end > as.numeric(x["end"])) %>%
      dplyr::select(found) %>%
      distinct() -> found
    found_1 = paste(found$found, collapse=",")
    CC_consensus_phases[[CC.2]] %>%
      filter(chr %in% x["chr"]) %>% 
      filter(start < as.numeric(x["start"]), end > as.numeric(x["end"])) %>%
      dplyr::select(found) %>%
      distinct() -> found
    found_2 = paste(found$found, collapse=",")
    c(found_1, found_2)
  })
  founders = data.frame(t(founders))
  colnames(founders) = c("CC.1","CC.2")
  
  sig_PO_annot <- data.frame(do.call("cbind", c(sig_PO, founders)))
  #sig_PO_annot$SDP = sapply(1:nrow(sig_PO_annot), function(x){
  #  lets = unique(c(unlist(strsplit(sig_PO_annot[x,"CC.1"],",")), unlist(strsplit(sig_PO_annot[x,"CC.2"],","))))
  #  nums = dict$num[which(dict$let %in% lets)]
  #  paste(unique(SDP_data[x,nums]), collapse="")
  #})
  return(sig_PO_annot)
}

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}


####################  matching  ######################

getMatches = function(df, id, in_groups, out_groups){
  levels_in  = apply(unique(df %>% select(in_groups)), 1, function(x) paste(x, collapse="_"))
  levels_out = apply(unique(df %>% select(out_groups)), 1, function(x) paste(x, collapse="_"))
  out = lapply(levels_in, function(x){
    groups = lapply(levels_out, function(y){
      df %>% filter(get(in_groups) == x, get(out_groups) == y) -> tmp
      if(nrow(tmp) > 0){
        tmp$rand = sample(1:nrow(tmp))
        tmp[,c(id, "rand")]
      }
    })
    jd = NULL
    if(!any(unlist(lapply(groups, is.null)))){
      jd = join_all(groups, type='inner', by="rand")
      colnames(jd) = c("Pup.ID_pos", "rand", "Pup.ID_neg")
    }
    jd
  })
  do.call("rbind", out)
}

