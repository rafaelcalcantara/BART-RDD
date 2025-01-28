rm(list=ls()[ls() %in% c("dgp","s","no_cores","n","c","Owidth") == FALSE])
source("simulation_estimator_functions.R")
## Create Results and Time folders to store outputs and runtimes
if ("Results" %in% list.files() == FALSE) dir.create("Results")
if ("Time" %in% list.files() == FALSE) dir.create("Time")
## Create Results and Time subfolders for current DGP
if (dgp %in% list.files("Results") == FALSE) dir.create(paste0("Results/",dgp))
if (dgp %in% list.files("Time") == FALSE) dir.create(paste0("Time/",dgp))
## Run simulations
cl <- makeCluster(no_cores,type="SOCK")
registerDoParallel(cl)
clusterExport(cl,varlist=ls())
time <- system.time({
  out <- parLapply(cl,1:s,fit_general)
})
stopCluster(cl)
write.table(time[3],paste0("Time/",dgp,"/total_parallel.csv"), append=TRUE, row.names = FALSE, col.names = FALSE, sep = ",")