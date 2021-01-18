#source("jags_kmers_source.R")  
#source("plotting_ratio_regression.R")

xce <- as.data.frame(rbind(rix1 = c("CC011", "CC001", "e", "a"),
                           rix2 = c("CC041", "CC051", "b", "b"),
                           rix3 = c("CC004", "CC017", "b", "e"),
                           rix4 = c("CC023", "CC047", "b", "b"),
                           rix6 = c("CC026", "CC006", "b", "b"),
                           rix7 = c("CC003", "CC014", "e", "a"),
                           rix8 = c("CC035", "CC062", "a", "a"),
                           rix9 = c("CC032", "CC042", "c", "a"),
                           rix10 = c("CC005", "CC040", "b", "b")))
colnames(xce) = c("num_cc", "dem_cc", "num_let", "dem_let")


xce_df = data.frame(founder = LETTERS[1:8],
                    xce_al = c("a","b","a","n","z","c","e","b"))


## revComp
revComp <- function(str){
  tmp <- unlist(strsplit(str, ""))
  revd <- rev(toupper(tmp))
  comp <- lapply(revd, function(x) ifelse(x == "A", "T", ifelse(x == "T", "A", ifelse(x == "G", "C", ifelse(x == "C", "G", "-")))))
  return(paste(unlist(comp),  collapse=""))
}


process_and_plot <- function(chr, 
                             snp_info, sample_info, RIX_info, 
                             ref_counts, alt_counts, 
                             phased_CC_haplotype=NULL, 
                             use_gene=F,
                             problemPups = c(1404, 1716, 1371, 569, 1911, 1951, 1015),
			                       seg_regions=NULL,#lb = 99e6, ub=104e6
			                       lb = NULL, ub=NULL){
  chr <- ifelse(chr == 20, "X", chr)
  chr_num <- ifelse(chr == "X", 20, as.numeric(chr))
  
  #### XO 1015, 1911 and 1951; low RNA-seq 1404, 1716, 1371, 569
  myCounts <- sort_counts(ref_counts, alt_counts)
  if(length(problemPups) > 0) myCounts <- myCounts %>% filter(!Pup.ID %in% problemPups)
  
  
  if(is.null(phased_CC_haplotype)){
    phased <- phase_counts(snps=snp_info, counts=myCounts, pupInfo=sample_info, CC_labels=RIX_info,
                           haplotypes, cc_haplotypes)
    plotMe <- plot_counts(snps=snp_info, counts=myCounts, pupInfo=sample_info, CC_labels=RIX_info,
                          keepEB=T, keepFS=F, keepFD=F)
  } else {
    plotMe <- plot_counts_based_on_founders(snps=snp_info, counts=myCounts, 
                                            pupInfo=sample_info, 
                                            CC_labels=RIX_info,
                                            keepEB=T, keepFS=F, keepFD=F,
                                            phased_CC_haplotype = phased_CC_haplotype, 
                                            seg_regions = seg_regions,
                                            lb=lb, ub=ub)
  }
  if(use_gene){
    data_kmers = plotMe$plotGenes
  } else {
    data_kmers = plotMe$plotMat
  }

  data_kmers$RRIX <- as.numeric(gsub("[a-z]","", data_kmers$RIX))
  data_kmers$pup_gene <- paste(data_kmers$Pup.ID, data_kmers$seq.Gene, sep="_")
  data_kmers$Diet <- factor(data_kmers$Diet, levels = sort(unique(data_kmers$Diet))[c(3,1,2,4)])
  data_kmers$RRIX <- factor(data_kmers$RRIX, levels = c(1:4, 6:10))
  data_kmers$DietRIX <- factor(paste0(data_kmers$Diet, "_",data_kmers$RRIX), 
                               levels = paste(unique(data_kmers$Diet)[c(3,1,2,4)], rep(c(1:4, 6:10), each=4), sep="_"))
  data_kmers$dir = factor(data_kmers$dir, levels=c("a","b"))
  
  data_kmers %>% arrange(Pup.ID, seq.Position) -> data_kmers
  
  return(data_kmers)
}


sort_counts <- function(ref, alt){
  if(is.null(colnames(ref))){
    colnames(ref) = colnames(alt) <- c("file", "k.mer", "dir", "counts")
    ref$Pup <- unlist(strsplit(as.character(ref$file), "/"))[c(F,F,F,F,T)]
    rem <- grep("old", ref$Pup)
    if(length(rem) > 0) ref <- ref[-rem,]
    ref$Pup.ID <- as.numeric(paste(gsub("[^0-9]", "", ref$Pup)))
    
    alt$Pup <- unlist(strsplit(as.character(alt$file), "/"))[c(F,F,F,F,T)]
    rem <- grep("old", alt$Pup)
    if(length(rem) > 0) alt <- alt[-rem,]
    alt$Pup.ID <- as.numeric(paste(gsub("[^0-9]", "", alt$Pup)))
  }
  alt = alt %>% arrange(pup.id, k.mer)
  ref = ref %>% arrange(pup.id, k.mer)
  tmp <- do.call("rbind", strsplit(as.character(alt$k.mer), ""))
  tmp <- tmp[,-((ncol(tmp)+1)/2)]
  consensus <- apply(tmp, 1, function(x) paste(x, collapse=""))
  alt$consensus = consensus
  alt$pup_consensus = paste0(alt$pup, "_", alt$consensus)
  alt$consensus_dir = sapply(1:nrow(alt), function(x) 
    ifelse(alt$counts[x] > 0, paste(alt$consensus[x], alt$dir[x], sep="_"), 0))
  
  tmp <- do.call("rbind", strsplit(as.character(ref$k.mer), ""))
  tmp <- tmp[,-((ncol(tmp)+1)/2)]
  consensus <- apply(tmp, 1, function(x) paste(x, collapse=""))
  ref$consensus = consensus
  ref$pup_consensus = paste0(ref$pup, "_", ref$consensus)
  ref$consensus_dir = sapply(1:nrow(ref), function(x) 
    ifelse(ref$counts[x] > 0, paste(ref$consensus[x], ref$dir[x], sep="_"), 0))
  
  alt = alt[match(ref$pup_consensus, alt$pup_consensus),]
  all.equal(ref$pup_consensus, alt$pup_consensus)
  
  #matches <- sapply(1:nrow(ref), function(x) 
  #  ifelse(ref$consensus_dir[x] == alt$consensus_dir[x], T,
  #         ifelse(ref$consensus_dir[x] == 0, T, 
  #         ifelse(alt$consensus_dir[x] == 0, T, F))))
  not_match = which(ref$consensus_dir != alt$consensus_dir)
  not_zero = intersect(which(ref$counts[not_match] > 2), which(alt$counts[not_match] > 2))
  get_rid <- not_match[not_zero]
  
  ref$allele = 0
  alt$allele = 1
  if(length(get_rid) > 0) ref = ref[-get_rid,]
  if(length(get_rid) > 0) alt = alt[-get_rid,]
  
  #if(identical(ref$consensus_dir, alt$consensus_dir) && identical(ref$pup.id, alt$pup.id)){
  ref <- ref %>% dplyr::select(-one_of("consensus_dir","allele")) %>%
    dplyr::rename("counts_ref"="counts","k.mer_ref"="k.mer")
  alt <- alt %>% dplyr::select(one_of("counts","k.mer")) %>%
    dplyr::rename("counts_alt"="counts","k.mer_alt"="k.mer")
  
  #all <- left_join(ref, alt, by="consensus") %>% distinct()
  all <- cbind(ref, alt) %>% distinct()
  
  all %>% mutate(counts = counts_ref+counts_alt) %>%
    group_by(consensus, dir) %>% 
    summarise(sum = sum(counts)) %>% 
    arrange(desc(sum)) %>% 
    dplyr::slice(1) %>%
    mutate(filt = paste(consensus, dir, sep="_")) -> keepDir
  
  all %>% 
    dplyr::rename("Pup.ID"="pup.id") %>%
    mutate(filt = paste(consensus, dir, sep="_")) %>%
    filter(filt %in% keepDir$filt) %>%
    dplyr::select(-one_of("filt")) %>%
    mutate(counts = counts_ref+counts_alt) -> counts
  return(counts)
}


#############################################

plot_counts_based_on_founders <- function(snps, counts, pupInfo, CC_labels, lo_bound = 3,
                                          keepEB=T, keepFS=T, keepFD=T,
                                          phased_CC_haplotype=NULL, seg_regions=NULL,
                                          lb = 99e6, ub=104e6){
  dict = data.frame(num = 1:8, let=LETTERS[1:8])
  keep_chr = as.character(unique(snps$seq.Chromosome))
  ref_counts = alt_counts = 
    matrix(NA, nrow=length(unique(counts$k.mer_ref)), ncol=length(unique(counts$Pup)))
  
  pups <- unique(counts$Pup.ID)
  temp_ref <- sapply(1:ncol(ref_counts), function(x){
    counts %>% filter(Pup.ID == pups[x]) -> temp
    temp[match(as.character(snps$seq.refseq), as.character(temp$k.mer_ref)),]$counts_ref
  }, simplify=T)
  
  temp_alt <- sapply(1:ncol(alt_counts), function(x){
    counts %>% filter(Pup.ID == pups[x]) -> temp
    temp[match(snps$seq.refseq, temp$k.mer_ref),]$counts_alt
  }, simplify=T)
  
  consensus <- snps$seq.consensus
  
  colnames(temp_alt) = colnames(temp_ref) = pups
  rownames(temp_alt) = rownames(temp_ref) = consensus
  
  if(length(grep("hetsOrZeros|keep",colnames(snps))) > 0){
    CC_keep <- data.frame(snps[,grep("keep|hetsOrZeros",colnames(snps))])
    rownames(CC_keep) <- snps$seq.consensus
    seq = snps[,grep("seq", colnames(snps))]
    #CCs <- colnames(snps)[grep("CC_ref",colnames(snps))]
    #tmp <- do.call("rbind", strsplit(CCs,"[.]")[c(F,T,T)])[,2:3]
    #colnames(CC_keep) <- apply(tmp, 1, function(x) paste(x, collapse="."))
    xc <- list(ref = temp_ref, alt = temp_alt, 
               seq = seq,
               CC_keep = CC_keep)
  } else{
    xc <- list(ref = temp_ref, alt = temp_alt, 
               seq = snps[,grep("seq", colnames(snps))],
               CC_keep = snps[,grep("sdp", colnames(snps))])
  }
  
  
  remove <- intersect(which(apply(xc$ref, 1, sum) == 0), which(apply(xc$alt, 1, sum) == 0))
  xc <- lapply(xc, function(x) x[-remove,])
  pupn <- colnames(xc$ref)#as.numeric(unlist(strsplit(as.character(colnames(xc$ref)),"_"))[c(F,F,T)])
  
  assign_founders = lapply(colnames(xc$ref), function(x){
    ind = grep(paste0("_",x,"$"), names(phased_CC_haplotype))
    if(ind > 0){     #length(ind) > 0
      if(CC_labels$CC.1[which(CC_labels$Pup.ID == x)] %in% phased_CC_haplotype[[ind]]$founder_by_parent[[1]]$cc){
        par1 = phased_CC_haplotype[[ind]]$founder_by_parent[[1]] %>% filter(chr == keep_chr)
        par2 = phased_CC_haplotype[[ind]]$founder_by_parent[[2]] %>% filter(chr == keep_chr)
      } else {
        par1 = phased_CC_haplotype[[ind]]$founder_by_parent[[2]] %>% filter(chr == keep_chr)
        par2 = phased_CC_haplotype[[ind]]$founder_by_parent[[1]] %>% filter(chr == keep_chr)
      }
      found1 = unlist(lapply(xc$seq$seq.Position, function(y)
        paste(unique(unlist(strsplit(par1$found[intersect(which(par1$start < y), which(par1$end > y))], ","))), collapse=",") ))
      
      found2 = unlist(lapply(xc$seq$seq.Position, function(y)
        paste(unique(unlist(strsplit(par2$found[intersect(which(par2$start < y), which(par2$end > y))], ","))), collapse=",") ))
    } else {
      found1 = found2 = ""
    }
    out = data.frame(cbind(found1, found2))
    colnames(out) = c(CC_labels$CC.1[which(CC_labels$Pup.ID == x)], CC_labels$CC.2[which(CC_labels$Pup.ID == x)])
    CC1_xce_hap = unique(out[which(xc$seq$seq.Position > lb & xc$seq$seq.Position < ub),1])
    CC1_xce_hap = sort(CC1_xce_hap[which(nchar(CC1_xce_hap) > 0)])
    CC2_xce_hap = unique(out[which(xc$seq$seq.Position > lb & xc$seq$seq.Position < ub),2])
    CC2_xce_hap = sort(CC2_xce_hap[which(nchar(CC2_xce_hap) > 0)])
    #out$CC1_xce_hap = unique(out[which(xc$seq$seq.Position > 99e6 & xc$seq$seq.Position < 103e6),1])
    #out$CC2_xce_hap = unique(out[which(xc$seq$seq.Position > 99e6 & xc$seq$seq.Position < 103e6),2])
    CC1_xce_al = factor(xce_df$xce_al[match(CC1_xce_hap, xce_df$founder)], levels = c("n","a","f","e","z","b","c","d"))
                #xce_df$xce_al[which(xce_df$founder %in% CC1_xce_hap)], 
                        
    CC2_xce_al = factor(xce_df$xce_al[match(CC2_xce_hap, xce_df$founder)], levels = c("n","a","f","e","z","b","c","d"))
    out$CC1_xce_hap = paste(CC1_xce_hap, collapse="/")
    out$CC2_xce_hap = paste(CC2_xce_hap, collapse="/")
    out$CC1_xce_al = paste(CC1_xce_al, collapse="/")
    out$CC2_xce_al = paste(CC2_xce_al, collapse="/")
    out$D = ifelse("n" %in% unique(c(paste(CC1_xce_al), paste(CC2_xce_al))), 1, 0)
    out$which_numer = ifelse(any(as.numeric(CC1_xce_al) < as.numeric(CC2_xce_al)), "left", "right")
    return(out)
  })
  
  names(assign_founders) <- colnames(xc$ref)
  assign_founders_num = lapply(assign_founders, function(x){
    cbind(unlist(lapply(x[,1], function(y) ifelse(y == "", "", dict$num[which(dict$let == y)]) )),
          unlist(lapply(x[,2], function(y) ifelse(y == "", "", dict$num[which(dict$let == y)]))))
  })
  
  tossOut <- c("mB", "Z", "H","M")
  if(!keepEB){
    tossOut <- c(tossOut, "eB")
  }
  if(!keepFS){
    tossOut <- c(tossOut, "fS")
  }
  if(!keepFD){
    tossOut <- c(tossOut, "fD")
  }
  
  ratios <- list()
  
  for(i in colnames(xc$ref)){
    alts <- list()
    refs <- list()
    if(!is.null(nrow(assign_founders_num[[i]])) & 
       length(grep(0, assign_founders_num[[i]])) < length(assign_founders_num[[i]])){
      if(!is.null(seg_regions)) rem_reg = seg_regions[[gsub("a|b", "", pupInfo$RIX[pupInfo$Pup.ID == i])]] %>% 
          filter(Chr == keep_chr)
      for(j in 1:2){
        useCol = apply(assign_founders_num[[i]], 1, function(k) as.numeric(k[j]))
        if(any(!is.na(useCol))){
          if(length(grep("keep|posCounts",colnames(snps))) > 0){
            keep <- sapply(1:length(useCol), function(k){
              tmp = xc$CC_keep[k,useCol[k]]
              ifelse(length(tmp) == 0, " ", tmp)
            })
          } else {
            keep <- sapply(1:length(useCol), function(k) ifelse(xc$CC_keep[k,useCol[k]] == 1,"A", 
                                                                ifelse(xc$CC_keep[k,useCol[k]] == 0, "R", xc$CC_keep[k,useCol[k]])))
          }
          problems <- grep(paste(tossOut,collapse="|"), keep)
          if(!is.null(rem_reg)) problems = sort(unique(union(problems, 
                                                             which(xc$seq$seq.Position > rem_reg$start & 
                                                                     xc$seq$seq.Position < rem_reg$end))))
          #problems <- lapply(keep, function(x) grep(paste(tossOut,collapse="|"), x))
          #xc$ref[problems,i] <- NA
          #xc$alt[problems,i] <- NA
          
          if(length(unlist(problems)) > 0){
            keep[which(lapply(problems, length) > 0)] <- NA
          }
          
          #refs[[j]] <- unlist(grep("R", keep))
          #alts[[j]] <- unlist(grep("A", keep))
          refs[[j]] <- which(keep %in% c("0","R"))
          alts[[j]] <- which(keep %in% c("1","A"))
        }
      }
      usable = NULL
      if(!any(unlist(lapply(refs, is.null)))){
        usable <- sort(union(unlist(refs[[1]][which(refs[[1]] %in% alts[[2]])]), 
                             unlist(alts[[1]][which(alts[[1]] %in% refs[[2]])])))
      }
      
      if(length(usable) > 0){
        #useVals <- sapply(1:length(usable), function(y){
        #  ifelse(keep[usable[y]] == "A", xc$ref[usable[y],i]/(xc$ref[usable[y],i]+xc$alt[usable[y],i]),
        #         xc$alt[usable[y],i]/(xc$ref[usable[y],i]+xc$alt[usable[y],i]))
        #})
        sums <- xc$alt[usable,paste(i)] + xc$ref[usable,paste(i)]
        usable <- usable[which(sums > lo_bound*2)]
        CC_numer <- sapply(1:length(usable), function(y){
          ifelse(length(grep("A", keep[usable[y]])) == 1, xc$ref[usable[y],paste(i)], xc$alt[usable[y],paste(i)])})
        CC_denom <- sapply(1:length(usable), function(y){
          ifelse(length(grep("A", keep[usable[y]])) == 1, xc$alt[usable[y],paste(i)], xc$ref[usable[y],paste(i)])})
        ratio <- matrix(NA, nrow=length(keep), ncol=2)  #rep(NA, length(keep))
        ratio[usable,] <- c(CC_numer, CC_denom)                   
        colnames(ratio) = c("CC_1", "CC_2")
        
        ratios[[i]] <- ratio
        
      } else {
        ratios[[i]] <- NULL
      }
    }
  }
  
  plotList <- list()
  
  plotList <- lapply(names(ratios), function(x){ 
    if(!is.null(ratios[[x]])){
      if(!all(is.na(as.matrix(ratios[[x]])))){
        
        s1 <- cbind(xc$seq, as.matrix(ratios[[x]]), 
                    #ref = xc$ref[,x], alt = xc$alt[,x], 
                    sum = (xc$ref[,x] + xc$alt[,x]),
                    CC1_hap = assign_founders[[paste(x)]][,1],
                    CC2_hap = assign_founders[[paste(x)]][,2],
                    assign_founders[[paste(x)]][,-c(1:2)])
        
        s1 <- s1[intersect(which(!is.na(s1$CC_1)), which(s1$sum > lo_bound)),]
        #s1 <- s1[,-c(1:4)]
        if(nrow(s1) > 0){
          s1$Pup.ID <- rep(as.numeric(x))
          s1$CCs <- rep(paste(colnames(assign_founders[[x]])[1:2], collapse = "/"))
          s1
        } else{
          NULL
        }
      } 
    }
  })
 
  plotGenes <- lapply(plotList, function(x){
    if(!is.null(x)){
      x %>% 
        group_by(seq.Gene) %>% 
        dplyr::select(-one_of("seq.TranscriptID","seq.end5","seq.end3","seq.consensus")) %>%
        mutate(minPos = min(seq.Position), maxPos = max(seq.Position)) %>% 
        group_by(seq.Chromosome, Pup.ID, CCs, seq.Gene, minPos, maxPos, CC1_hap, CC2_hap,
                 CC1_xce_hap, CC2_xce_hap, CC1_xce_al, CC2_xce_al, D, which_numer) %>% 
        summarize(sumTot = sum(sum), sum1 = sum(CC_1), sum2 = sum(CC_2)) 
      
      
    } else {
      NULL
    }
  })
  plotGenes <- unique(do.call("rbind", plotGenes))
  plotGenes$logSumTot <- log(plotGenes$sumTot)
  
  ### add pup demographic info
  pupInfo$Diet <- gsub(" $","", pupInfo$Diet, perl=T)
  
  
  plotGenes <- left_join(plotGenes, pupInfo, by="Pup.ID") %>%
    mutate(Xist = if_else(seq.Gene == "Xist", T, F))
  
  plotMat <- unique(do.call("rbind", plotList))
  plotMat$logSum <- log(plotMat$sum)
  # plotMat$Pup.ID <- as.numeric(unlist(strsplit(as.character(plotMat$Pup),"_"))[c(F,F,T)])
  #plotMat$Pup.ID <- as.numeric(plotMat$Pup)
  
  plotMat <- left_join(plotMat, pupInfo, by="Pup.ID") %>%
    mutate(Xist = if_else(seq.Gene == "Xist", T, F))
  
  #####################################
  plotMat$fin1 <- plotMat$CC_1
  plotMat$fin2 <- plotMat$CC_2
  plotMat[which(plotMat$Xist),"fin1"] <- plotMat[which(plotMat$Xist),"CC_2"]
  plotMat[which(plotMat$Xist),"fin2"] <- plotMat[which(plotMat$Xist),"CC_1"]
  
  tmp <- plotMat[which(plotMat$which_numer == "right"),"fin1"]
  plotMat[which(plotMat$which_numer == "right"),"fin1"] <- plotMat[which(plotMat$which_numer == "right"),"fin2"]
  plotMat[which(plotMat$which_numer == "right"),"fin2"] <- tmp
  plotMat$kRat <- plotMat$fin1 / plotMat$sum
  plotMat$kLB = qbeta(0.05,plotMat$fin1 + 1, plotMat$fin2 + 1)
  plotMat$kUB = qbeta(0.95,plotMat$fin1 + 1, plotMat$fin2 + 1)
  plotMat$kLB[which(is.nan(plotMat$kRat))] = plotMat$kUB[which(is.nan(plotMat$kRat))] <- 0
  
  plotMat$CC_lab = plotMat$CCs
  plotMat[which(plotMat$which_numer == "right"),"CC_lab"] = apply(do.call("rbind", strsplit(plotMat$CC_lab[which(plotMat$which_numer == "right")],"/"))[,2:1], 1,
                                                                        function(x) paste(x,collapse="/"))
  
  plotGenes$fin1 <- plotGenes$sum1
  plotGenes$fin2 <- plotGenes$sum2
  plotGenes$sum = plotGenes$fin1 + plotGenes$fin2
  plotGenes[which(plotGenes$Xist),"fin1"] <- plotGenes[which(plotGenes$Xist),"sum2"]
  plotGenes[which(plotGenes$Xist),"fin2"] <- plotGenes[which(plotGenes$Xist),"sum1"]
  
  tmp = plotGenes[which(plotGenes$which_numer == "right"),"fin1"]
  plotGenes[which(plotGenes$which_numer == "right"),"fin1"] <- plotGenes[which(plotGenes$which_numer == "right"),"fin2"]
  plotGenes[which(plotGenes$which_numer == "right"),"fin2"] <- tmp
  plotGenes$seq.Position = apply(plotGenes, 1, function(x) (as.numeric(x["minPos"]) + as.numeric(x["maxPos"]))/2)
  plotGenes = plotGenes %>% as.data.frame %>% 
    mutate(gRat = fin1/sumTot, 
           gLB = qbeta(0.05, sum1 + 1, sum2 + 1),
           gUB = qbeta(0.95, sum1 + 1, sum2 + 1))
  
  plotGenes$CC_lab = plotGenes$CCs
  plotGenes[which(plotGenes$which_numer == "right"),"CC_lab"] = apply(do.call("rbind", strsplit(plotGenes$CC_lab[which(plotGenes$which_numer == "right")],"/"))[,2:1], 1,
                                                                  function(x) paste(x,collapse="/"))

  return(list(plotMat=plotMat, plotGenes=plotGenes, metrics=NULL))
}


