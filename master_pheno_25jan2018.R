library(rjags)
library(coda)
library(ggplot2)
library(tidyverse)
library(data.table)
library(lmerTest)

setwd("C:/Users/Kathie/matnut/src")

source("./matnut/summary_functions.R")
source("./matnut/tables_n_plots.R")

## Read in data
dir <- "C:/Users/Kathie/Dropbox\ (ValdarLab)"
data  <- read.csv(file.path(dir, "matnut_main/AllPhenotypes_Matnut4.csv"))
ptypes <- colnames(data)[44:ncol(data)]
terms <- c("RIX","PORIX","Diet","DietRIX", "PODietRIX","DamID","BreedingBatch","BehaviorBatch")
data$DamID <- factor(data$DamID)
data$Diet <- factor(gsub(" ", "",data$Diet), levels=unique(gsub(" ", "",data$Diet))[c(3,1,2,4)])
POlevels <- unlist(lapply(paste(c(1:4, 6:10)), function(x) paste0(as.character(paste(x)), c("a","b"))))
data$Reciprocal <- factor(data$Reciprocal, levels=POlevels)
data$RIX <- factor(data$RIX)
data$PO <- ifelse(gsub("[0-9]","",data$Reciprocal) == "a", 0.5, -0.5)
data$PORIX <- factor(paste0("PO", data$RIX), levels=paste0("PO", levels(data$RIX)))
levels <- unlist(lapply(levels(data$RIX), function(x) paste0(levels(data$Diet), x)))[-c(33,36)]
data$DietRIX <- factor(paste0(data$Diet, data$RIX), levels=levels)
levels <- unlist(lapply(levels(data$PORIX), function(x) paste(x, levels(data$Diet))))[-c(65,68,69,72)]
data$PODietRIX <- factor(paste0("PO",data$DietRIX), levels=paste0("PO", levels(data$DietRIX)))

encoded <- getEncoding(data, terms)
matnut <- list(df=data, 
               ptypes=ptypes,
               encoded=encoded)
#location of data on killdevil: dataSource <- "~/Data/matnut_outputs/"
#matnut <- readRDS(file.path(dir,'/phenotype_analysis/matnut_data.rds'))


## Run linear model
myPhen <- makeSummary(datalist=matnut, phenotype=matnut$ptypes[1:2], 
                        randvar=c("DamID", "RIX", "Diet:RIX"), fixvar="Diet", POvar=c("RIX", "Diet:RIX"),
                      tryLam=c(-1, 0, .25, .33, .5, 1, 2, 3), contrasts=F, 
                        chains=2, n.adapt=200, n.iter=500)
#saveRDS(myPhen, paste0(dataSource, "weight_phen_redo.rds"))

## Using previously generated data
#myPhen <- readRDS(file.path(dir, "phenotype_analysis/regTest_all_14dec2018.rds")
myPhen <- readRDS(file.path(dir, "phenotype_analysis/complete_phen_model.rds"))

#myPhen <- lapply(myOrig, function(x) x[[1]])

## Generate all tables
tables <- writeTables(results = myPhen, df=matnut, output=paste0(dir,"/phenotype_analysis/out/"), print=T)

## Generate all ribbon plots

ribPlots <- generateRibbonPlots(results=myPhen, df=matnut, ptypes=names(myPhen$allSummary)[1], 
                                output=paste0(dir,"/phenotype_analysis/out/"), print=F)
dietFlags <- flagPlots(myPhen, matnut, ptypes=names(myPhen$allSummary), print=T, byVar = "Diet", 
                       output=paste0(dir,"/phenotype_analysis/out/"))
rixFlags  <- flagPlots(myPhen, matnut, ptypes=names(myPhen$allSummary), print=T, byVar = "RIX", 
                       output=paste0(dir,"/phenotype_analysis/out/"))

## Generate catepillar plots
graphics.off()
pdf(file.path(paste0(dir,"allPheno_JLComparison_9may2018.pdf")), 
    onefile=TRUE, width = 8.5, height = 11)
for(i in 1:length(myPhen$plot)){  
  print(myPhen$plot[[i]])
}
graphics.off() 

