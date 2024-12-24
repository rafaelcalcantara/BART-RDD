setwd("~/Git/BART-RDD")
library(doParallel)
library(XBART)
### Parallelization
no_cores <- 11
## Create results folder
if (!dir.exists("Results")) dir.create("Results")
## Get simulation script names
scripts <- list.files("R")[grep("simulation",list.files("R"))]
scripts <- scripts[-grep("master|data|results|polynomial_1|polynomial_2|polynomial_3",scripts)]
scripts <- paste0("R/",scripts)
## Samples to estimate
s0 <- 1
# s1 <- 1000
## Take command line arguments
args <- commandArgs(trailingOnly = T)
if (length(args) > 0)
{
  files <- args
} else
{
  files <- 1:length(list.files("Data"))
}
### Parameters
p_categorical <- 0
ntrees        <- 10
num_sweeps    <- 150
burnin        <- 50
# Ow            <- c(0.015,0.007,0.003,0.002)
for (j in scripts)
{
  print(paste0("Script: ",j))
  time <- system.time({
    source(j)
  })
  print(paste0("Elapsed time: ",time))
}
