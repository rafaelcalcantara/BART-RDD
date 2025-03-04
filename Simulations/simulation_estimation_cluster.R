args <- commandArgs(trailingOnly = TRUE)
s0 <- as.numeric(args[1])
s1 <- as.numeric(args[2])
dgp <- args[3]
n <- as.numeric(args[4])
c <- as.numeric(args[5])
Owidth <- as.numeric(args[6])
k1 <- as.numeric(args[7])
k2 <- as.numeric(args[8])
k3 <- as.numeric(args[9])
k4 <- as.numeric(args[10])
k5 <- as.numeric(args[11])
p <- as.numeric(args[12])
rho <- as.numeric(args[13])
source("simulation_estimator_functions.R")
source("simulation_process_results.R")
## Create relevant folders for current DGP
if (dgp %in% list.files("Results") == FALSE) dir.create(paste0("Results/",dgp))
if (dgp %in% list.files("Time") == FALSE) dir.create(paste0("Time/",dgp))
if (dgp %in% list.files("Results/RMSE") == FALSE) dir.create(paste0("Results/RMSE/",dgp))
if (dgp %in% list.files("Results/Screenshots") == FALSE) dir.create(paste0("Results/Screenshots/",dgp))
if (dgp %in% list.files("Results/Fits") == FALSE) dir.create(paste0("Results/Fits/",dgp))
if ("cate" %in% list.files(paste0("Results/",dgp)) == FALSE) dir.create(paste0("Results/",dgp,"/cate"))
if ("cate" %in% list.files(paste0("Results/Fits/",dgp)) == FALSE) dir.create(paste0("Results/Fits/",dgp,"/cate"))
if ("yhat" %in% list.files(paste0("Results/",dgp)) == FALSE) dir.create(paste0("Results/",dgp,"/yhat"))
if ("yhat" %in% list.files(paste0("Results/Fits/",dgp)) == FALSE) dir.create(paste0("Results/Fits/",dgp,"/yhat"))
## Run simulations
time <- system.time({
  for (i in s0:s1)
  {
    fit_general(i)
  }
  # screenshot(s0,s1)
})
write.table(time[3],paste0("Time/",dgp,"/total_cluster.csv"), append=TRUE, row.names = FALSE, col.names = FALSE, sep = ",")