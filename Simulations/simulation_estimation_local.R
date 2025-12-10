source("simulation_estimator_functions.R")
source("simulation_process_results.R")
## Create relevant folders for current DGP
if (dgp %in% list.files("Time") == FALSE) dir.create(paste0("Time/",dgp))
if (dgp %in% list.files("Results/RMSE") == FALSE) dir.create(paste0("Results/RMSE/",dgp))
# if (dgp %in% list.files("Results/Screenshots") == FALSE) dir.create(paste0("Results/Screenshots/",dgp))
if (dgp %in% list.files("Results/Fits") == FALSE) dir.create(paste0("Results/Fits/",dgp))
if ("cate" %in% list.files(paste0("Results/Fits/",dgp)) == FALSE) dir.create(paste0("Results/Fits/",dgp,"/cate"))
## Run simulations
cl <- makeCluster(no_cores,type="SOCK")
registerDoParallel(cl)
clusterExport(cl,varlist=ls())
time <- system.time({
  out <- parLapply(cl,1:s,fit_general)
})
stopCluster(cl)
# time <- system.time({
#   for (j in 1:s)
#   {
#     fit_general(j)
#   }
#   # screenshot(s0,s1)
# })
write.table(time[3],paste0("Time/",dgp,"/total_parallel.csv"), append=TRUE, row.names = FALSE, col.names = FALSE, sep = ",")
