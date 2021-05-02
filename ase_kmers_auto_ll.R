options(digits=4)
args <- commandArgs(trailingOnly = TRUE) 
r = as.numeric(args[1])

setwd("/nas/longleaf/home/kys6/TReC_matnut/src")
library(tidyverse)
library(fdrtool)
library(DESeq2)
library(rjags)
library(apeglm)
library(VGAM)


source("prediction_functions.R")
source("summary_functions.R")
source("ase_summary_source.R")
#source("stan_pheno_functions.R")


dir <- "/nas/depts/006/valdar-lab/users/sunk"
data_kmers_list = readRDS(file.path(dir, "trec/data_kmers_from_process_and_plot/data_kmers_list_out_29jan2021.rds"))


#for (r in c(2:4,6:10)){
	test_dat = do.call("rbind", lapply(data_kmers_list, function(x) x %>% filter(RRIX == r)))
	test_dat = test_dat %>% group_by(Pup.ID, seq.Gene, Diet, Reciprocal, pup_gene) %>% 
	      		summarize(mat_tot = sum(CC_1), pat_tot = sum(CC_2)) %>%
	          	mutate(Pup.ID = factor(Pup.ID), 
				seq.Gene = factor(seq.Gene), 
				Diet = factor(Diet), 
				Reciprocal = factor(Reciprocal))
    	test_dat$pup_gene = factor(test_dat$pup_gene)
      
      	test_dat$tot = test_dat$mat_tot + test_dat$pat_tot
      	encoded = getEncoding(test_dat, terms=c("Pup.ID","seq.Gene", "Diet", "Reciprocal", "pup_gene") )

        fit_list = list()
        for (g in encoded$Level[which(encoded$Variable == "seq.Gene")]){
		temp_dat = test_dat %>% filter(seq.Gene == g)
		try(fit_list[[g]] <- vglm(cbind(mat_tot, pat_tot) ~ 1, betabinomial (zero = 2),
			data = temp_dat, trace = TRUE),
		        silent=T)
	        
	        #coef(fit2)
	}
	
	saveRDS(fit_list, paste0(dir, "/ase_vgam/ase_intOnly_rix",r,"_1may2021.rds"))  
#}
