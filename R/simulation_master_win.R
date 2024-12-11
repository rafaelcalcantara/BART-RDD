# if (substr(sessioninfo::os_name(),1,3) == "Win")
# {
#   setwd("~/../Git/BART-RDD/")
# } else
# {
#   setwd("~/Git/BART-RDD/")
# }
setwd("../")
source("R/simulation_data.R")
library(doParallel)
library(XBART)
### Parallelization
no_cores <- detectCores()-1
## Create results folder
if (!dir.exists("Results")) dir.create("Results")
## Get simulation script names
scripts <- list.files("R")[grep("win",list.files("R"))]
scripts <- scripts[-grep("master|data|llr0|llr1|llr2|cgs0|cgs1|results|bcf",scripts)]
scripts <- paste0("R/",scripts)
### Parameters
p_categorical <- 0
ntrees        <- 20
num_sweeps    <- 150
burnin        <- 50
Ow            <- c(0.066,0.034,0.03,0.014)
for (i in scripts)
{
  print(paste0("Script: ",i))
  time <- system.time({
    source(i)
  })
  print(paste0("Elapsed time: ",time))
}