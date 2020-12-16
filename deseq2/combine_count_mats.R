options(stringsAsFactors = FALSE)
library(DESeq2)
library(tidyverse)
library(sva)

setwd("~/rna_seq/deseq2")
dir="/nas/depts/006/valdar-lab/users/sunk/"
source("deseq2_functions.R")


countData<- as.matrix(read.csv(file.path(dir,"string_pipe_out/gene_count_matrix_old.csv"), row.names="gene_id"))

countData_2 <- as.matrix(read.csv(file.path(dir,"string_pipe_out/gene_count_matrix_redo.csv"), row.names="gene_id"))
redo_pups = colnames(countData_2)

countData_2 = countData_2[match(rownames(countData), rownames(countData_2)), ]
identical(rownames(countData), rownames(countData_2))
identical(colnames(countData)[match(redo_pups, colnames(countData))], colnames(countData_2))
countData[,match(redo_pups, colnames(countData))] = countData_2
identical(countData[,match(redo_pups, colnames(countData))], countData_2)

write.csv(countData, file.path(dir, "string_pipe_out/gene_count_matrix.csv"))

######################

transcriptData<- as.matrix(read.csv(file.path(dir,"string_pipe_out/transcript_count_matrix_old.csv"), row.names="transcript_id"))

transcriptData_2 <- as.matrix(read.csv(file.path(dir,"string_pipe_out/transcript_count_matrix_redo.csv"), row.names="transcript_id"))
redo_pups = colnames(transcriptData_2)

transcriptData_2 = transcriptData_2[match(rownames(transcriptData), rownames(transcriptData_2)), ]
identical(rownames(transcriptData), rownames(transcriptData_2))
identical(colnames(transcriptData)[match(redo_pups, colnames(transcriptData))], colnames(transcriptData_2))
transcriptData[,match(redo_pups, colnames(transcriptData))] = transcriptData_2
identical(transcriptData[,match(redo_pups, colnames(transcriptData))], transcriptData_2)

write.csv(transcriptData, file.path(dir, "string_pipe_out/transcript_count_matrix.csv"))
