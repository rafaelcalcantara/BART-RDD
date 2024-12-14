setwd("~/Git/BART-RDD")
library(doParallel)
library(XBART)
### Parallelization
no_cores <- 125
## Generate data
source("R/simulation_data.R")
rm(list=ls())
gc()
no_cores <- 125
## Create results folder
if (!dir.exists("Results")) dir.create("Results")
## Get simulation script names
scripts <- list.files("R")[grep("simulation",list.files("R"))]
scripts <- scripts[-grep("master|data|results",scripts)]
scripts <- paste0("R/",scripts)
### Parameters
p_categorical <- 0
ntrees        <- 20
num_sweeps    <- 150
burnin        <- 50
Ow            <- c(0.03,0.04,0.007,0.02)
s0 <- 1
s1 <- 125
for (i in scripts)
{
  print(paste0("Script: ",i))
  time <- system.time({
    source(i)
  })
  print(paste0("Elapsed time: ",time))
}
