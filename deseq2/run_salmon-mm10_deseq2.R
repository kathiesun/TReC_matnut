#BiocManager::install("Biostrings")
library(tidyverse)
library(tximport)
library(tximportData)
library(DESeq2)
library(GenomicFeatures)
library(tximeta)
library(Biostrings)
library(seqinr)
library(ensembldb)


setwd("C:/Users/Kathie/rna_seq/kmerSearch")
setwd("~/rna_seq/kmerSearch)/kmerSearch/")
#dataSource <- "C:/Users/Kathie/Dropbox\ (ValdarLab)/variant_data"
dir <- "/nas/depts/006/valdar-lab/users/sunk/"

#source("masterKmers_source.R")
#source("~/matnut/src/matnut/plot.hpd_ks.R")
#source("plotting_ratio_regression.R")
#source("ase_summary_single.R")
#library(cmdline)
#sim <- cmdline.integer("sim")

samples <- read.csv(file.path(dir, "variant_data", "2015-10_expression_pups.csv"))
rownames(samples) = samples$Pup.ID
samples$Diet <- gsub(" $", "", samples$Diet)
samples$Diet <- factor(samples$Diet, levels=c("Standard", "Low Protein","Methyl Enriched","Vitamin D Deficient"))
samples$PO <- ifelse(factor(gsub("[0-9]","", samples$RIX)) == "a", 0.5, -0.5)
samples$PO_cat <- factor(gsub("[0-9]","",samples$RIX))
#use_pups <- list.files(file.path(dir,"mm10/salmon_ref_bt"), "Pup.ID")#, paste0("Pup.ID_",rownames(samples))))
#use_pups <- list.files(file.path(dir,"matnut_sea_outputs/salmon"), "Pup.ID")#, paste0("Pup.ID_",rownames(samples))))
use_pups <- list.files(file.path(dir,"mm10_transcriptome_quantification"), "Pup.ID")#, paste0("Pup.ID_",rownames(samples))))
use_pups <- as.numeric(gsub("Pup.ID_","",use_pups))
allfiles <- file.path(dir,"mm10_transcriptome_quantification/", paste0("Pup.ID_",rownames(samples)), "quant.sf")
#allfiles <- file.path(dir,"matnut_sea_outputs/salmon", paste0("Pup.ID_",rownames(samples)),"transcripts_quant", "quant.sf")
allfiles <- data.frame(Pup.ID=rownames(samples), files=allfiles)
files <- allfiles[which(allfiles$Pup.ID %in% use_pups),]
files$Pup.ID = as.numeric(paste(files$Pup.ID))
#pups <- lapply(use_files, function(x){ tmp=strsplit(x, "/")[[1]]
#	       tmp[grep("Pup.ID", tmp)]
#	})
#pups <- unlist(strsplit(unlist(pups),"_"))[c(F,T)]
		      
samples_use <- samples %>% 
	dplyr::filter(Pup.ID %in% as.numeric(files$Pup.ID)) %>%
	dplyr::select(one_of("Pup.ID","RRIX","Diet","PO","PO_cat","Breeding.Batch","Behavior.Batch")) %>%
	mutate(Breeding.Batch  = as.factor(Breeding.Batch), Behavior.Batch = as.factor(Behavior.Batch))
coldata <- left_join(samples_use, files, by="Pup.ID")

ref_fa <- readDNAStringSet(file.path(dir, "mm10/ref_sequence/Mus_musculus.GRCm38.cdna.all.fa"))

seq_name = names(ref_fa)  
seq_det = lapply(seq_name, function(x) {
	unlist(strsplit(x, " "))[c(1,3,4,7)]
	#keep = grep("chr|gene:|gene_symbol|ENS", tmp)
	#tmp[keep]
})
seq_det <- do.call("rbind", seq_det)

#### tx2gene for ref mm10 ####
tx2gene <- data.frame(TXNAME=seq_det[,1], GENEID=gsub("gene_symbol:","",seq_det[,4]))

transcripts <- do.call("rbind", strsplit(seq_det[,1], "[.]"))
chr <- do.call("rbind", strsplit(seq_det[,2], "[:]"))
gene <- do.call("rbind", strsplit(seq_det[,3], "[:|.]"))
gene_sym <- do.call("rbind", strsplit(seq_det[,4], "[:]"))
seq_det <- cbind(transcripts, chr, gene, gene_sym)
colnames(seq_det) = c("transcript","tr_num", "del","build","chr","start","end","dir","del2","gene","gene_num","del3","gene_name")
seq_det <- seq_det[,-grep("del", colnames(seq_det))]
seq_det <- data.frame(seq_det)
seq_det %>% mutate(start = as.numeric(paste(start)), 
		   end = as.numeric(paste(end))) -> seq_det
seq_det$sequence = paste(ref_fa)


#### for combined or novel (?) transcriptomes using tximeta ####
txome <- makeLinkedTxome(file.path(dir, "mm10/salmon/salmon.index"), source="scallop,ensembl", organism="Mus musculus", release="96", genome="GRCm38",
			 fasta="/nas/depts/006/valdar-lab/users/sunk/mm10/ref_sequence/Mus_musculus.GRCm38.cdna.all.fa", 
			 gtf="/nas/depts/006/valdar-lab/users/sunk/mm10/ref_sequence/Mus_musculus.GRCm38.96.gtf",
	      write = TRUE) #, jsonFile=paste0(dir, "/mm10/deseq"))

txome <- makeLinkedTxome(file.path(dir, "mm10/salmon_ref_bt/mm10_index"), source="Ensembl", organism="Mus musculus", release="96", genome="GRCm38",
			 fasta="/nas/depts/006/valdar-lab/users/sunk/mm10/ref_sequence/Mus_musculus.GRCm38.cdna.all.fa", 
			 gtf="/nas/depts/006/valdar-lab/users/sunk/mm10/ref_sequence/Mus_musculus.GRCm38.96.gtf",
	      write = TRUE) #, jsonFile=paste0(dir, "/mm10/deseq"))

#### IGNORE: tx2gene for reference mm10 transcriptome ####
txi.tx <- tximport(as.character(files$files), type = "salmon", txOut = TRUE)

###### for mm10 ref #######
txdb <-makeTxDbFromGFF(file.path(dir, "mm10/ref_sequence/Mus_musculus.GRCm38.96.gtf"), format="gtf", 
		       dataSource="Ensembl", organism="Mus musculus")
k <- keys(txdb, keytype = "TXNAME")
tx2gene_gtf <- select(txdb, k, "GENEID", "TXNAME")

#txdb <- makeTxDbFromEnsembl(organism="Mus musculus")

#sal_dir <- system.file(paste0(dir,"/mm10/salmon_ref_bt/Pup.ID_", pup,"/transcripts_quant/quant.sf"))
file <- file.path(sal_dir, "quant.sf")

###### for 188 scallop'd merged thing ######
txdb <-makeTxDbFromGFF(file.path(dir, "/mm10/scallop/gtf_merge/gffall.annotated.gtf"), organism="Mus musculus", dataSource ="Ensembl 3/14/19")
#se <- tximeta(coldata, type="salmon")

txdump <- as.list(txdb)
k <- keys(txdb, keytype = "GENEID")
df <- select(txdb, keys = k,  columns = "TXNAME", keytype = "GENEID")


######################
dds_lst <- list()
res <- list()
keepCounts_list <- list()
res0.1 <- list()
for(i in unique(coldata$RRIX)){
	coldata_use <- coldata %>% dplyr::filter(RRIX == i)    
	pups <- as.character(coldata_use$Pup.ID)
	rixFiles <- as.character(coldata_use$files)
	names(rixFiles) <- pups
#	txdb <-makeTxDbFromGFF(paste0(dir, "/mm10/scallop/bt_output/RIX_",i,"/gffall.annotated.gtf"), format="gtf", organism="Mus musculus", dataSource =paste0("scallop_RIX ",i))
#se <- tximeta(coldata, type="salmon")

#txdump <- as.list(txdb)
#k <- keys(txdb, keytype = "GENEID")
#df <- select(txdb, keys = k,  columns = "TXNAME", keytype = "GENEID")


#x <- read.table(paste0(dir, "/mm10/scallop/bt_output/RIX_",i,"/gffall.annotated.gtf"), header=F, sep="\t")
#x <- apply(x, 1, function(x){
#	tmp <- unlist(strsplit(x, ";| "))
#	tmp[which(tmp != "")]
#})
#len <- unlist(lapply(x, length)) 
#tx_14 <- t(sapply(which(len == 14), function(n) x[[n]]))
#tx_18 <- t(sapply(which(len == 18), function(n) x[[n]]))
#tx_22 <- t(sapply(which(len == 22), function(n) x[[n]])) 
#tx_18_new <- cbind(tx_18[,1:12], matrix(NA, nrow=nrow(tx_18), ncol=2), tx_18[,13:14], 
#		 matrix(NA, nrow=nrow(tx_18), ncol=2), tx_18[,15:18])
#tx <- data.frame(rbind(tx_18_new, tx_22))
#which(!tx_14[,"V95"] %in% c(tx_18[,"V95"], tx_22[,"V95"]))
#tx2gene_gtf <- data.frame(TXNAME=tx[,"V92"], GENEID=tx[,"V14"])
	#x <- do.call("rbind", x)
#	txmap <- read.table(paste0(dir, "/mm10/scallop/bt_output/RIX_",i,"/gffall.scallop.gtf.tmap"), header=T)
#	txmap <- read.table(paste0(dir, "matnut_sea_outputs/scallop/output/RIX_",i,"/gffall.scallop.gtf.tmap"), header=T)
	#tx2gene <- data.frame(TXNAME=txmap$qry_id, GENEID=txmap$ref_gene_id)  
	#head(tx2gene)
	#tx2gene <- apply(tx2gene, c(1,2), as.character)
	#tx2gene[which(tx2gene[,"GENEID"] == "-"), "GENEID"] = as.character(tx2gene[,"TXNAME"][which(tx2gene[,"GENEID"] == "-")])
	#tx2gene <- data.frame(tx2gene)
	
	##### getting RIX-specific transcript #####
	#assem_fa <- readDNAStringSet(paste0(dir, "/mm10/scallop/bt_output/RIX_",i,"/union.fa"))

	#assem_name = names(assem_fa)  
	#assem = lapply(assem_name, function(x) {
	#	tmp = unlist(strsplit(x, " "))
	#	if(length(tmp) > 1) {
	#		tmp[c(1,3,4,7)]
	#	} else {
	#		c(tmp, rep(NA, 3))
	#	}
	#keep = grep("chr|gene:|gene_symbol|ENS", tmp)
	#tmp[keep]
	#})
	#assem <- do.call("rbind", assem)

	#unique <- grep("gene", assem[,1])
	#add_uniq <- lapply(assem_name[unique], function(x)
	#		   c(x, rep(NA,12)))
	#add_uniq <- do.call("rbind", add_uniq)
	#transcripts <- do.call("rbind", strsplit(assem[-unique,1], "[.]"))
	#chr <- do.call("rbind", strsplit(assem[-unique,2], "[:]"))
	#gene <- do.call("rbind", strsplit(assem[-unique,3], "[:|.]"))
	#gene_sym <- do.call("rbind", strsplit(assem[-unique,4], "[:]"))
	#assem <- cbind(transcripts, chr, gene, gene_sym)
	#assem <- rbind(add_uniq, assem)
	#colnames(assem) = c("transcript","tr_num", "del","build","chr","start","end","dir","del2","gene","gene_num","del3","gene_name")
	#assem <- assem[,-grep("del", colnames(assem))]
	#assem <- data.frame(assem)
	#assem %>% mutate(start = as.numeric(paste(start)), 
	#	end = as.numeric(paste(end))) -> assem
	#assem$sequence = paste(assem_fa)
	

	txi <- tximport(rixFiles, type="salmon", tx2gene=tx2gene)

	use <- which(colnames(txi$counts) %in% coldata_use$Pup.ID)
	txi_use <- list(abundance=txi$abundance[,use], 
			counts=txi$counts[,use], 
			length=txi$length[,use], 
			countsFromAbundance=txi$countsFromAbundance)
	dedata <- DESeqDataSetFromTximport(txi_use,
                                   colData = coldata_use,
                                   design = ~ Diet + PO + Diet:PO)
	dds <- DESeq(dedata) 
	dds <- nbinomWaldTest(dds)
	#res <- results(dds, c("PO_cat","b","a"))
	res <- results(dds, name="PO")
	res[[paste0("rix",i)]] <- res
	res0.1[[paste0("rix",i)]] <- res[which(res$padj < 0.1),]
	keepCounts <- rowSums(counts(dds)) >= 10
	keepCounts <- counts(dds)[which(keepCounts==T),]
	keepCounts <- keepCounts[order(-rowSums(keepCounts)),]
	keepCounts_list[[paste0("rix",i)]] <- keepCounts
	dds_lst[[paste0("rix",i)]] <- dds
}
#dds_ref <- dds_lst
saveRDS(dds_lst, file.path(dir, "mm10_transcriptome_quantification/deseq/dds_lst-15sep2019.rds"))

########## scallop RIX-specific transcripts ###############

ref_rix <- readDNAStringSet(file.path(dir, "mm10/scallop/bt_output/RIX_1/unique.fa"))

which(names(ref_rix) == "gene.79697.0.0")
rix_name = names(ref_rix) 
sequences <- data.frame(ref_rix) 
rix_det <- data.frame(names = rix_name, 
	width = unlist(lapply(sequences, nchar)),
	seq = sequences$ref_rix)	
rix_det[12156,]
