setwd("~/Git/BART-RDD")
rm(list = ls())
gc()
N <- c(500,1000,2500)
sig_error <- 3
pts_in_window <- 75
source("R/simulation_data.R")
rm(list=ls())
gc()
Omin <- 5
Opct <- 0.3
source("R/simulation_master.R")
rm(list=ls())
gc()
source("R/simulation_results.R")
rm(list=ls())
gc()
source("R/sim_plots_tables.R")