seq1 = read.table(paste0(dir, "/xci_paper_data/phylo_seqs/beast1/sequence1_align.txt"), sep="\n")
seq1 = as.character(seq1$V1[-grep("\\*", seq1$V1)])
seq1 = gsub(" +", " ", seq1)
seq1 = data.frame(do.call("rbind", strsplit(as.character(seq1), " |\t", seq1)))
colnames(seq1) = c("CC", "seq","n")
if(length(which(seq1$CC == "")) > 0) seq1 = seq1[-which(seq1$CC == ""),]
CCs = do.call("rbind", strsplit(as.character(unique(seq1$CC)), "-|_"))[,1]
if(length(which(toupper(CCs) == "X")) > 0){
  CCs[which(toupper(CCs) == "X")] = "B6"
  seq1$CC = gsub("X", "B6", toupper(seq1$CC))
}

strings = data.frame(CC = CCs, string = "")
strings$CC = as.character(strings$CC)
strings$string = apply(strings, 1, function(x){
  paste(seq1$seq[grep(as.character(paste(x["CC"])), as.character(seq1$CC))], collapse="")
})
strings = strings[order(strings$CC), ]
header = paste0("#NEXUS

BEGIN DATA;
\tDIMENSIONS NTAX=", length(unique(strings$CC)), " NCHAR=", unique(unlist(lapply(strings$string, nchar))), ";
\tFORMAT MISSING=N GAP=- DATATYPE=DNA;
\tMATRIX")

ender = "\t;
END;"

seq1_out = paste(header, paste0("\t", apply(strings, 1, function(x) paste0(x, collapse="\t")), collapse="\n"), ender, sep="\n")
write.table(seq1_out, col.names = F, row.names = F, quote = F, 
            file.path(dir, "xci_paper_data/phylo_seqs/sequence1_102709036-102711871.nex"))


#####################################################################
