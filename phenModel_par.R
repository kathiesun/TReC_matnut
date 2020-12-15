library(rjags)
library(coda)
library(ggplot2)
library(data.table)
library(lmerTest)

setwd("/nas02/home/k/y/kys6/matnut/src")

source("./matnut/summary_functions.R")

library(cmdline)
sim <- cmdline.integer("sim")

###### Read in data
dataSource <- file.path("~/Data/matnut_outputs")
matnut <- readRDS(file.path(dataSource,'matnut_data.rds'))

####################
#   Run function   #
####################
runAll <- function(df,
                   output_dir,
                   tryLam=1, 
                   perJob=500,
                   sim) { 

  myPhen <- makeSummary(datalist=df, phenotype=df$ptypes[sim], 
                        randvar=c("DamID", "RIX", "Diet:RIX"), fixvar="Diet", POvar=c("RIX", "Diet:RIX"),
                        tryLam=c(-1, 0, .25, .33, .5, 1, 2, 3), contrasts=T, 
                        chains=2, n.adapt=20000, n.iter=50000)
  
  saveRDS(myPhen, paste0(output_dir,"/modelFit_phen",sim,".rds"))				
}

runAll(matnut, output_dir=paste(dataSource, "par_out",sep="/"), 
	tryLam=c(-1, 0, .25, .33, .5, 1, 2, 3), sim=sim)


testPhen <- list()
for(j in 1:20){
  myPhen <- readRDS(paste0(dataSource,"/par_out/modelFit_phen",j,".rds"))
  phenName <- names(myPhen[["allSummary"]])
  for (i in 1:length(names(myPhen))){
    testName <- names(myPhen)[i]
    #if(testName == "mcmcObject") testPhen[[testName]][[phenName]]<- myPhen[[testName]][[phenName]][[phenName]]
    if(testName == "lmer_obj") testPhen[[testName]][[phenName]]<- myPhen[[testName]] 
    else testPhen[[testName]][[phenName]]<- myPhen[[testName]][[phenName]]
  }
}


