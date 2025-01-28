cluster <- TRUE
n <- 4000
s <- 100
source("generate_master_scripts.R")
for (i in 1:1) source(paste0("simulation_master",i,".R"))