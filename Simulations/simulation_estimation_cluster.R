args <- commandArgs(trailingOnly = TRUE)
s0 <- as.numeric(args[1])
s1 <- as.numeric(args[2])
dgp <- args[3]
n <- as.numeric(args[4])
c <- as.numeric(args[5])
Owidth <- as.numeric(args[6])
source("simulation_estimator_functions.R")
source("simulation_process_results")
## Create relevant folders for current DGP
if (dgp %in% list.files("Results") == FALSE) dir.create(paste0("Results/",dgp))
if (dgp %in% list.files("Time") == FALSE) dir.create(paste0("Time/",dgp))
if (dgp %in% list.files("Results/RMSE") == FALSE) dir.create(paste0("Results/RMSE/",dgp))
if (dgp %in% list.files("Results/Screenshots") == FALSE) dir.create(paste0("Results/Screenshots/",dgp))
## Run simulations
time <- system.time({
  for (i in s0:s1)
  {
    fit_general(i)
  }
  screenshots(s0,s1)
})
write.table(time[3],paste0("Time/",dgp,"/total_cluster.csv"), append=TRUE, row.names = FALSE, col.names = FALSE, sep = ",")