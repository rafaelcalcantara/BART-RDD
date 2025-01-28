cluster <- TRUE
n <- 4000
s <- 200
source("generate_master_scripts.R")
for (i in 1:14) source(paste0("simulation_master",i,".R"))