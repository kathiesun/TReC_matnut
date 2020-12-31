library(tidyverse)

setwd("C:/Users/Kathie/TReC_matnut//src")

source("./matnut/prediction_functions.R")
source("./matnut/summary_functions.R")

dir <- "C:/Users/Kathie/Dropbox\ (ValdarLab)"
gecco <- read.csv(file.path(dir, "de_results/fullGeccoRnaDump.csv"))
regFiles = list.files(pattern = "summary_mnt_50k_29dec2020.rds", file.path(dir, "variant_data/regression_outputs/"))

readRDS(grep("chr17", regfiles)