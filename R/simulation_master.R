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
## Run estimation scripts
for (j in scripts)
{
  print(paste0("Script: ",j))
  time <- system.time({
    source(j)
  })
  print(paste0("Elapsed time: ",time))
}
