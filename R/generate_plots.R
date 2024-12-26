setwd("~/Git/BART-RDD")
rm(list = ls())
gc()
N <- c(500,1000,1500)
rho <- c(0.7,0.94)
pts_in_window <- c(75,75)
source("R/simulation_data.R")
rm(list=ls())
gc()
Omin <- 2
Opct <- 0.4
s1 <- 50
source("R/simulation_master.R")
rm(list=ls())
gc()
samples <- 1:50
source("R/simulation_results.R")
rm(list=ls())
gc()
rho <- c(0.7,0.94)
source("R/sim_plots_tables.R")