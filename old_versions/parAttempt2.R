setwd("/nas02/home/k/y/kys6/matnut/src")
source("~/matnut/src/matnut/parArgs.R")

sumFunc <- "~/matnut/src/matnut/summary_functions.R"
source(sumFunc) 

cov_short <- readRDS("~/Data/matnut_outputs/cov_short_data.rds")
library(car)
library(cmdline)
sim <- cmdline.integer("sim")


######## Define function to run in accumulator

# remove missing
origPhen <- grep("ENS", colnames(cov_short$df))
missing <- which(colSums(cov_short$df[, origPhen]) == 0)
if(length(missing)>0) cov_short$df <- cov_short$df[, -missing]

# rank all other genes
origPhen <- grep("ENS", colnames(cov_short$df))
noPhen <- setdiff(seq(1, ncol(cov_short$df), by=1), origPhen)
allrank <- cbind( cov_short$df[,noPhen],
                   t(apply(cov_short$df[,origPhen], 1, function(x) rank(x, ties.method = "min")/length(!is.na(x)))) )

indvariable="~ -1 + Diet + (PO + 0 | RIX) + (PO + 0|Diet:RIX)"

runAll <- function(rankdara=allrank, noPhen=noPhen, phenInd, encoded=cov_short$encoded, 
                   indvariable=indvariable, tryLam=1){
  system.time ({
    datause <- rankdata[,c(noPhen, phenInd)]
    fitLMER <- BC.model(y.mat=datause[,-noPhen], data=datause,
                         indvariable=indvariable, 
                         transformParams=getMatnutTransformParams(tryLam = tryLam, normd = T))
    
    lmerSum <- lmer.getDecodedSummary(lmerOb=fitLMER2, formula=indvariable, phenotypes=fitLMER2$phenotypes)
    
    fitLMER$JAGSformula <- lmer.to.jags.form(indvariable)
    jagsFit2 <- runJagsModels_cov(datalist=datause, 
                                 testLMER = fitLMER2, encoded=encoded, phenotype=fitLMER2$phenotypes, n.iter=200000)

    
    #jagsSum <- list()
    #all.y <- matrix(NA, nrow=nrow(datause), ncol=length(fitLMER2$y.transform))
    #for(i in 1:length(fitLMER2$y.transform)){
    #  all.y[,i] <- fitLMER2$y.transform[[i]]
    #  all.mcmc <- mcmc.stack(jagsFit2[[i]]$fit)
    #  jagsSum[[i]] <- jags.getDecodedSummary(all.mcmc, encoded)
    #}
    
    #allSummary <- list(lmer=lmerSum, jags=jagsSum)
    
    #allSummary$compare <- sapply(1:length(lmerSum), function(x) data.frame(compareSummary(allSummary$lmer[[x]],allSummary$jags[[x]])), 
    #                             simplify = F)
    return(lmerSum)   #allSummary)
  })
}


######### Define accumulator


tmpdir = "~/job.tmp/"

system.type = "killdevil"

accum = parallel$get.cluster.accum(system.type       = system.type,
                                   func = runAll,
                                   sharedVariables   = list(rankdata=allrank, 
                                                            noPhen=noPhen, 
                                                            encoded=cov_short$encoded,
                                                            indvariable=indvariable),
                                   filesToSource     = sumFunc,
                                   batchSize         = 1,
                                   timeLimit.hours   = .15,
                                   cpuMemLimit.GB    = .5,
                                   coresPerJob       = 1,
                                   maxSimulJobs      = 100,
                                   systemOpts        = c(),
                                   outdir            = tmpdir,
                                   retryFailing      = T,
                                   saveProp          = T)

perJob = 500
total = length(origPhen)

for (i in 1:2){    #seq(1:(floor(total/perJob)+1))){
  a = (perJob * (i-1)) + length(noPhen)+1
  b = a+perJob - 1
  if(i==(floor(total/perJob)+1)) b = total + length(noPhen)
  phenInd <- seq(a,b,by=1)
  #test <- runAll(rankdata=allrank, noPhen=noPhen, phenInd=phenInd, encoded=cov_short$encoded, 
  #               indvariable=indvariable)
  #saveRDS(test, paste(output_dir, paste0("rnaseq_allRIX_",i,".rds"), sep="/"))
  
  accum$addCall(list(phenInd=phenInd))
}

######### Store outputs
outputs = accum$runAll()
harvested = accum$getAllOutputs(outputs, removeFailing = T)


#outputs = accum$runAll()

#iter = accum$getOutputIterator(outputs)

output_dir <- "~/Data/matnut_outputs/rna"
saveRDS(harvested, paste(output_dir, "rnaseq_allRIX.rds", sep="/"))
#counter=1

#while(iter$hasNext())
#{
#  item = iter$nextItem()
#  saveRDS(item, paste(output_dir, paste0("rnaseq_allRIX_",counter,".rds"), sep="/"))
#  counter=counter+1
#}
