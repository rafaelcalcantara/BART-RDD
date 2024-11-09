setwd("~/Git/BART-RDD/")
source("R/simulation_data.R")
library(doParallel)
library(XBART)
### Parallelization
no_cores <- detectCores()-1
## Get simulation script names
scripts <- list.files("R")[grep("win",list.files("R"))]
# scripts <- scripts[-grep("master|data|llr0|llr2|cgs0|results|bcf",scripts)]
scripts <- scripts[-grep("master|data|llr0|llr1|llr2|cgs0|cgs1|results|bcf",scripts)]
scripts <- paste0("R/",scripts)
p_categorical <- 1
Owidth        <- 0.2
for (i in scripts)
{
  print(paste0("Script: ",i))
  time <- system.time({
    source(i)
  })
  print(paste0("Elapsed time: ",time))
}