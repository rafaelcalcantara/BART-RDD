setwd("~/../Git/BART-RDD/")
source("R/simulation_data.R")
library(doParallel)
### Parallelization
no_cores <- detectCores()-1
## Get simulation script names
scripts <- list.files("R")[grep("win",list.files("R"))]
scripts <- scripts[-grep("master|data|llr0|llr2|cgs0|results|bcf",scripts)]
scripts <- paste0("R/",scripts)
p_categorical <- 2
Owidth        <- 0.1
for (i in scripts)
{
  print(paste0("Script: ",i))
  time <- system.time({
    source(i)
  })
  print(paste0("Elapsed time: ",time))
}