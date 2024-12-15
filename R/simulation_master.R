setwd("~/Git/BART-RDD")
library(doParallel)
library(XBART)
## Take command line arguments
args <- commandArgs(trailingOnly = T)
if (length(args) > 0)
{
  s0 <- args[1]
  s1 <- args[2]
} else
{
  s0 <- 1
  s1 <- 1000
}
### Parallelization
no_cores <- 125
# ## Generate data
# source("R/simulation_data.R")
# rm(list=ls())
# gc()
# no_cores <- 125
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
Ow            <- c(0.066,0.0425,0.029,0.01)
for (j in scripts)
{
  print(paste0("Script: ",j))
  time <- system.time({
    source(j)
  })
  print(paste0("Elapsed time: ",time))
}
