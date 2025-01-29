source("simulation_estimator_functions.R")
source("simulation_process_results.R")
## Create relevant folders for current DGP
if (dgp %in% list.files("Results") == FALSE) dir.create(paste0("Results/",dgp))
if (dgp %in% list.files("Time") == FALSE) dir.create(paste0("Time/",dgp))
if (dgp %in% list.files("Results/RMSE") == FALSE) dir.create(paste0("Results/RMSE/",dgp))
if (dgp %in% list.files("Results/Screenshots") == FALSE) dir.create(paste0("Results/Screenshots/",dgp))
if (dgp %in% list.files("Results/Fits") == FALSE) dir.create(paste0("Results/Fits/",dgp))
## Run simulations
cl <- makeCluster(no_cores,type="SOCK")
registerDoParallel(cl)
clusterExport(cl,varlist=ls())
time <- system.time({
  out <- parLapply(cl,1:s,fit_general)
})
stopCluster(cl)
write.table(time[3],paste0("Time/",dgp,"/total_parallel.csv"), append=TRUE, row.names = FALSE, col.names = FALSE, sep = ",")