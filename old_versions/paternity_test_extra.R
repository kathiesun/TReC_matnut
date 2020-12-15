######## Paternity test
setwd("/nas/depts/006/valdar-lab/users/sunk")
#setwd("C:/DB Mount/Dropbox\ (ValdarLab)/outputs/")
library(tidyverse)

# exon boundaries
exons <- read.csv("CCVariants/SNPsInExons.csv", colClasses=c("FounderSDP"="character"))

# CC mice in matnut
myCC <- read.csv("CCVariants/MatNut_RIX-Xceinfo_020315.csv")
myCC <- as.character(myCC$Alt.Name)[which(myCC$Alt.Name != "")]

# founder haplotypes for CC
filenames <- list.files("CCVariants/foundHaps", pattern="*.csv", full.names=TRUE)
ldf <- lapply(filenames, read.csv)
tmp <- do.call("rbind", strsplit(filenames, "/"))
names(ldf) <- unlist(strsplit(tmp[,ncol(tmp)], "_"))[c(T,F)]
founderHaps <- ldf
founders <- c("AJ", "C57BL6J", "129", "NOD","NZO","CAST","PWK","WSB")

curateKmers <- function(exons, myCC, founders=c("AJ", "C57BL6J", "129", "NOD","NZO","CAST","PWK","WSB"), 
			variants, founderHaps){
  chr <- unique(variants$Chromosome)
  ### read in data
	exonsUse <- exons[which(exons$Chromosome == chr),]
	founderHapsUse <- lapply(founderHaps, function(x) na.omit(x[which(x$chromosome == chr),]))
	
	### define haplotype blocks
	aa <- grep("AA", colnames(founderHapsUse[[1]]))
	hh <- grep("HH", colnames(founderHapsUse[[1]]))
	hapBlocks <- lapply(founderHapsUse, function(x) unlist(sapply(1:nrow(x), function(y) {
	  tmp <- which(x[y,c(aa:hh)] > 0.99)
	  ifelse(length(tmp) > 0, tmp, NA)}, simplify=T)))
	for(i in 1:length(founderHapsUse)){
	  tmp<- na.omit(cbind(founderHapsUse[[i]], hapBlocks = as.matrix(hapBlocks[[i]])))
	  founderHapsUse[[i]] <- tmp[order(tmp$position.B38.),]
	}
	for(i in 1:length(founderHapsUse)){
	  founderHapsUse[[i]]$place <- 0
	  founderHapsUse[[i]]$group_num <- 1
	    
	  for(j in 1:length(founderHapsUse[[i]]$position.B38.)){
	    if(j==1){
	      curBlock=founderHapsUse[[i]]$hapBlocks[j]
	      group_num=1
	      founderHapsUse[[i]]$place[j]=1
	    } else if (founderHapsUse[[i]]$hapBlocks[j] != curBlock){
	      curBlock=founderHapsUse[[i]]$hapBlocks[j]
	      founderHapsUse[[i]]$place[j]=1
	      founderHapsUse[[i]]$place[j-1]=2
	      group_num=group_num+1
	    }
	    founderHapsUse[[i]]$group_num[j]=group_num
	  }
	  founderHapsUse[[i]]$place[j]=2
	}
	
	hapSimp <- lapply(founderHapsUse, function(x)
	                  x %>% arrange(position.B38.) %>% group_by(hapBlocks) %>% 
	                    filter(place %in% c(1,2)) %>%
	                    select(one_of("chromosome","position.B38.","hapBlocks","place","group_num")) %>%
	                    arrange(position.B38.)) 
	                  
	blockLim <- lapply(hapSimp, function(x) arrange(spread(x, key="place", value="position.B38."), group_num))

	### Get rid of variants within 25 bp to exon start/stop
	avoid <- rbind(cbind((unique(exonsUse$start) - 25), (unique(exonsUse$start) + 25)), cbind((unique(exonsUse$end) - 25), (unique(exonsUse$end) + 25)))
	tooClose <- sapply(1:length(variants$PositionB38), function(x) 
			   ifelse(length(intersect(which(variants$PositionB38[x] < avoid[,2]), which(variants$PositionB38[x] > avoid[,1]))) > 0, F, T))
	variants <- variants[tooClose,]
	
	### Separate out ref and alt alleles, CC and founder counts
	ccMiceNum <- grep("CC", colnames(variants)) 
	ccMice <- unlist(strsplit(colnames(variants)[ccMiceNum], "[.]"))[c(T,F)]
	keepccMiceNum <- ccMiceNum[which(ccMice %in% myCC)]
	ccMiceCounts <- variants[,keepccMiceNum]
	ccMiceCounts_ref <- apply(ccMiceCounts,2,function(x) as.numeric(paste(unlist(strsplit(x, "/"))[c(T,F)])))
	ccMiceCounts_alt <- apply(ccMiceCounts,2,function(x) as.numeric(paste(unlist(strsplit(x, "/"))[c(F,T)])))

	founderMiceNum <- unlist(lapply(founders, function(x) grep(x, colnames(variants))[1]))
	founderMice <- colnames(variants)[founderMiceNum]
	founderCounts <- variants[,founderMiceNum]
	founderCounts_ref <- apply(founderCounts,2,function(x) as.numeric(paste(unlist(strsplit(x, "/"))[c(T,F)])))
	founderCounts_alt <- apply(founderCounts,2,function(x) as.numeric(paste(unlist(strsplit(x, "/"))[c(F,T)])))

	seq <- as.data.frame(do.call("rbind",strsplit(as.character(variants$ProbeSeq), "\\[|\\]|/")))
	colnames(seq) <- c("end5","ref","alt","end3")
	seq$refseq <- paste0(seq$end5, seq$ref, seq$end3)
	seq$altseq <- paste0(seq$end5, seq$alt, seq$end3)

	sdp <- do.call("rbind", strsplit(variants$FounderSDP,""))
	isAlt <- sapply(1:nrow(sdp), function(x) intersect(which(sdp[x,] == 1), which(founderCounts_ref[x,] == 0)))
	isRef <- sapply(1:nrow(sdp), function(x) intersect(which(sdp[x,] == 0), which(founderCounts_ref[x,] > 0)))
	isAlt2 <- sapply(1:nrow(sdp), function(x) intersect(which(sdp[x,] == 1), which(founderCounts_alt[x,] > 0)))
	isRef2 <- sapply(1:nrow(sdp), function(x) intersect(which(sdp[x,] == 0), which(founderCounts_alt[x,] == 0)))
	altTrue <- sapply(1:nrow(sdp), function(x) ifelse(identical(isAlt[[x]], isAlt2[[x]]), 1, 0), simplify = T)
	refTrue <- sapply(1:nrow(sdp), function(x) ifelse(identical(isRef[[x]], isRef2[[x]]), 1, 0), simplify = T)
	bothTrue <- intersect(which(altTrue == 1), which(refTrue == 1))
	
	### trim out variants that disagree with founder sdp
	var_trim <- list()
	var_trim[["snp_id"]] <- variants[bothTrue, 1:4]
	var_trim[["seq"]] <- seq[bothTrue,]
	var_trim[["sdp"]] <- sdp[bothTrue,]
	var_trim[["founder_ref"]] <- founderCounts_ref[bothTrue,]
	var_trim[["founder_alt"]] <- founderCounts_alt[bothTrue,]
	var_trim[["CC_ref"]] <- ccMiceCounts_ref[bothTrue,]
	var_trim[["CC_alt"]] <- ccMiceCounts_alt[bothTrue,]
	var_trim[["founder_map"]] <- matrix(NA, nrow=dim(var_trim[["CC_alt"]])[1], ncol=dim(var_trim[["CC_alt"]])[2])
	var_trim[["keep"]] <- matrix(NA, nrow=dim(var_trim[["CC_alt"]])[1], ncol=dim(var_trim[["CC_alt"]])[2])
	
	### make sure counts match with expected founders
  for(i in 1:ncol(var_trim[["CC_ref"]])){
    CC=unlist(strsplit(colnames(var_trim[["CC_ref"]])[i],"[.]"))[c(T,F)]
    var_trim[["founder_map"]][,i] <- sapply(1:length(var_trim[["CC_ref"]][,i]), function(y){
      lo=which(var_trim[["snp_id"]][y,]$PositionB38 > blockLim[[CC]]$`1`)
      hi=which(var_trim[["snp_id"]][y,]$PositionB38 < blockLim[[CC]]$`2`)
      groupInd=ifelse(length(lo)==0, 1, ifelse(length(hi)==0, nrow(blockLim[[CC]]),
                                               intersect(which(var_trim[["snp_id"]][y,]$PositionB38 > blockLim[[CC]]$`1`), 
                                                         which(var_trim[["snp_id"]][y,]$PositionB38 < blockLim[[CC]]$`2`))))
      as.numeric(blockLim[[CC]][groupInd,"hapBlocks"])
    })
    
    var_trim[["keep"]][,i] <- sapply(1:length(var_trim[["CC_ref"]][,i]), function(y){
      ifelse(var_trim[["sdp"]][y,var_trim[["founder_map"]][y,i]] == "0" && var_trim[["CC_ref"]][y,i] > 0, T,
           ifelse(var_trim[["sdp"]][y,var_trim[["founder_map"]][y,i]] == "1" && var_trim[["CC_ref"]][y,i] == 0, T, F))
    })
  }
    
	keepFoReal <- which(apply(var_trim[["keep"]], 1, sum) == ncol(var_trim[["keep"]]))
	
	### make sure counts are ~ 30
	realCounts <- sapply(1:length(keepFoReal), function(x){
                  	   ref <- which(var_trim[["sdp"]][keepFoReal[x],] == "0")
                  	   alt <- which(var_trim[["sdp"]][keepFoReal[x],] == "1")
                  	   refInd <- which(var_trim[["founder_map"]][keepFoReal[x],] %in% ref)
                  	   altInd <- which(var_trim[["founder_map"]][keepFoReal[x],] %in% alt)
                  	   
                  	   c(var_trim[["CC_ref"]][keepFoReal[x],refInd], var_trim[["CC_alt"]][keepFoReal[x],altInd])
                  	  })
  kmer_means <- colMeans(realCounts)	
  kmer_sd <- apply(realCounts, 2, sd)
  kmer_range <- apply(realCounts, 2, range)
  z <- (30 - kmer_means) / kmer_sd
  remove <- which(pnorm(30, kmer_means, kmer_sd) < 0.05)
  if(length(remove) > 0) keepFoReal <- keepFoReal[-remove]
  
  ### final list of variants 
  #var_trim <- list()
  var_trim[["snp_id"]] <- var_trim[["snp_id"]][keepFoReal,]
  var_trim[["seq"]] <- var_trim[["seq"]][keepFoReal,]
  var_trim[["sdp"]] <- var_trim[["sdp"]][keepFoReal,]
  var_trim[["founder_ref"]] <- var_trim[["founder_ref"]][keepFoReal,]
  var_trim[["founder_alt"]] <- var_trim[["founder_alt"]][keepFoReal,]
  var_trim[["CC_ref"]] <- var_trim[["CC_ref"]][keepFoReal,]
  var_trim[["CC_alt"]] <- var_trim[["CC_alt"]][keepFoReal,]
  var_trim[["founder_map"]] <- var_trim[["founder_map"]][keepFoReal,]
  var_trim[["keep"]] <- var_trim[["keep"]][keepFoReal,]
  if(sum(var_trim[["keep"]] == F) > 0){
    print("something went wrong")
  } else {
    var_trim[["keep"]] <- NULL
  }
  return(var_trim)
}

# vcf files

perfectKmers <- list()
kmersToRun_ref <- list()
kmersToRun_alt <- list()

for(i in 1:20){
  buffer = ""
  if(nchar(i) == 1) buffer = "0"
  varUse <- read.csv(paste0("CCVariants/CCVariants",buffer,i,".csv"), colClasses=c("FounderSDP"="character"))
  perfectKmers[[i]] <- curateKmers(exons, myCC, founders, variants=varUse, founderHaps)
  kmersToRun_ref[[i]] <- perfectKmers[[i]]$seq$refseq
  kmersToRun_alt[[i]] <- perfectKmers[[i]]$seq$altseq
  saveRDS(kmersToRun_ref, "CCVariants/kmersToRun_ref.rds")
  saveRDS(kmersToRun_alt, "CCVariants/kmersToRun_alt.rds")
  saveRDS(perfectKmers, "CCVariants/perfectKmers.rds")
}


