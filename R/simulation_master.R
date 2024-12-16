setwd("~/Git/BART-RDD")
library(doParallel)
library(XBART)
### Parallelization
no_cores <- 125
## Create results folder
if (!dir.exists("Results")) dir.create("Results")
## Get simulation script names
scripts <- list.files("R")[grep("simulation",list.files("R"))]
scripts <- scripts[-grep("master|data|results",scripts)]
scripts <- paste0("R/",scripts)
scripts <- "R/simulation_oracle.R"
## Samples to estimate
s0 <- 1
s1 <- 1000
## Take command line arguments
args <- commandArgs(trailingOnly = T)
if (length(args) > 0)
{
  files <- args
} else
{
  files <- length(list.files("Data"))
}
### Parameters
p_categorical <- 0
ntrees        <- 10
num_sweeps    <- 150
burnin        <- 50
Ow            <- c(0.029,0.018,0.029,0.01)
for (j in scripts)
{
  print(paste0("Script: ",j))
  time <- system.time({
    source(j)
  })
  print(paste0("Elapsed time: ",time))
}
