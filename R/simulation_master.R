setwd("~/Git/BART-RDD")
library(doParallel)
### Parallelization
no_cores <- 10
## Create results folder
if (!dir.exists("Results")) dir.create("Results")
## Get simulation script names
scripts <- list.files("R")[grep("simulation",list.files("R"))]
scripts <- scripts[grep("barddt_v2|tbart|sbart|polynomial",scripts)]
scripts <- paste0("R/",scripts)
## Samples to estimate
s0 <- 1
## Take command line arguments
args <- commandArgs(trailingOnly = T)
if (length(args) > 0)
{
  files <- args[1]
  s1 <- args[2]
} else
{
  files <- 1:length(list.files("Data"))
  s1 <- 10
}
## Stochtree paramlist
global.parmlist <- list(standardize=T,sample_sigma_global=TRUE,sigma2_global_init=0.01)
mean.parmlist <- list(num_trees=50, min_samples_leaf=20, alpha=0.95, beta=10, max_depth=10, sample_sigma2_leaf=FALSE)
## Run estimation scripts
for (j in scripts)
{
  print(paste0("Script: ",j))
  time <- system.time({
    source(j)
  })
  print(paste0("Elapsed time: ",time))
}
