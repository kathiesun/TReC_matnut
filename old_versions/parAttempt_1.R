setwd("/nas02/home/k/y/kys6/matnut/src")
source("~/matnut/src/matnut/parArgs.R")
source("~/matnut/src/matnut/boxcox_functions.R")
lmerFunc <- "~/matnut/src/matnut/lmer_functions_rna.R"
source(lmerFunc) 
cov_short <- readRDS("~/Data/matnut_outputs/cov_short_data.rds")
library(car)

runLMER <- function(cov_short, rix, tryLam=c(0, 0.5)){
  system.time ({
    datause <- cov_short$df[which(paste(cov_short$df$RIX) == rix), ]
    phenInd <- grep("ENS",colnames(datause))
    fitLMER <- BC.model(y.mat=datause[,phenInd], data=datause[,-phenInd], 
                        indvariable="~ PO + Diet", print=T, 
                        transformParams=getMatnutTransformParams(tryLam = c(0, 0.5), normd = T))
    lmerSum <- lmer.getDecodedSummary(lmerOb=fitLMER, formula="~ PO + Diet", phenotypes=fitLMER$phenotypes)
    rnaseq[[j]] <- return(list(summary=lmerSum, lmerObj=fitLMER))
  })
}

tmpdir = "~/job.tmp/"

system.type = "killdevil"

accum = parallel$get.cluster.accum(system.type       = system.type,
                                   func = runLMER,
                                   sharedVariables   = list(cov_short=cov_short),
                                   
                                   filesToSource     = lmerFunc,
                                   batchSize         = 1,
                                   timeLimit.hours   = .15,
                                   cpuMemLimit.GB    = .5,
                                   coresPerJob       = 1,
                                   maxSimulJobs      = 100,
                                   systemOpts        = c(),
                                   outdir            = tmpdir,
                                   retryFailing      = T,
                                   saveProp          = T)
 
rixes <- unique(cov_short$df$RIX)
seq <- seq(1:9)

for (i in seq){
	accum$addCall(list(rix = rixes[i]))
}

outputs = accum$runAll()

iter = accum$getOutputIterator(outputs)

output_dir <- "~/Data/matnut_outputs/rna"
counter=1

while(iter$hasNext())
{
  item = iter$nextItem()
  saveRDS(item, paste(output_dir, paste0("rnaseq_out2_",rixes[seq[counter]],".rds"), sep="/"))
  counter=counter+1
}
