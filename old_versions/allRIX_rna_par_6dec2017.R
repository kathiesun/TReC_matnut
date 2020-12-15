setwd("/nas02/home/k/y/kys6/matnut/src")
sumFunc <- "~/matnut/src/matnut/summary_functions.R"
source(sumFunc) 

cov_short <- readRDS("~/Data/matnut_outputs/cov_short_data.rds")
library(car)
library(cmdline)
sim <- cmdline.integer("sim")

#bsub -J [1-2] R CMD BATCH --vanilla --args --sim=\$LSB_JOBINDEX allRIX_killdevil_6dec2017.R


######## Define function to run in accumulator

output_dir <- "~/Data/matnut_outputs/rna"

runAll <- function(df,
                   output_dir,
		               tryLam=1, 
                   perJob=500,
                   sim) {  
  system.time ({
    # remove missing
    data <- df$df
    encoded <- df$encoded
    origPhen <- grep("ENS", colnames(data))
    missing <- which(colSums(data[, origPhen]) == 0)
    if(length(missing)>0) data <- data[, -missing]
    
    # rank all other genes
    origPhen <- grep("ENS", colnames(data))
    noPhen <- setdiff(seq(1, ncol(data), by=1), origPhen)
    allrank <- cbind( data[,noPhen],
                      t(apply(data[,origPhen], 1, function(x) rank(x, ties.method = "min")/length(!is.na(x)))) )
    total <- length(origPhen)
 
    # determine genes for this job
    a = (perJob * (sim-1)) + length(noPhen)+1
    b = a+perJob - 1
    if(sim==(floor(total/perJob)+1)) b = total + length(noPhen)
    phenInd <- seq(a,b,by=1)
    datause <- allrank[,c(noPhen, phenInd)]
    
    allSummary <- makeSummary(datalist=datause, phenotype=colnames(datause)[-noPhen], tryLam=c(1), normd=T, plot=F,
                        randvar=c("RIX", "Diet:RIX"), fixvar="Diet", POvar=c("RIX", "Diet:RIX"), encoded=encoded)
   allSummary$mcmcObject <- NULL
 
    return(allSummary)
  })
}

test <- runAll(df=cov_short, 
       output_dir=output_dir,
       sim=sim)

saveRDS(test, paste(output_dir, paste0("rnaseq_allRIX_",sim,".rds"), sep="/"))
