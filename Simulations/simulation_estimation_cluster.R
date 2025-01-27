source("simulation_estimator_functions.R")
## Fit from cluster
if ("leaf.rdd" %in% models)
{
  barddt.global.parmlist <- list(standardize=T,sample_sigma_global=TRUE,sigma2_global_init=0.01)
  barddt.mean.parmlist <- list(num_trees=150, min_samples_leaf=20, alpha=0.95, beta=2,
                               max_depth=20, sample_sigma2_leaf=FALSE)
  barddt.var.parmlist <- list(num_trees = 2,min_samples_leaf = 10)
  print("BARDDT fit")
  print(paste(c("k1","k2","k3","k4","k5","p","rho"),args[1:7],sep=": "))
  cl <- makeCluster(no_cores,type="SOCK")
  registerDoParallel(cl)
  clusterExport(cl,varlist=ls())
  time <- system.time({
    out <- parLapply(cl,1:s,fit_cluster_barddt)
  })
  stopCluster(cl)
  print("BARDDT fit done!")
  print(time)
  # saveRDS(out,paste0("Results/barddt_",dgp,".rds"))
}
if ("tbart" %in% models)
{
  tbart.global.parmlist <- list(standardize=T,sample_sigma_global=TRUE,sigma2_global_init=0.01)
  tbart.mean.parmlist <- list(num_trees=150, min_samples_leaf=20, alpha=0.95, beta=2,
                              max_depth=20, sample_sigma2_leaf=FALSE)
  tbart.var.parmlist <- list(num_trees = 2,min_samples_leaf = 10)
  print("T-BART fit")
  print(paste(c("k1","k2","k3","k4","k5","p","rho"),args[1:7],sep=": "))
  cl <- makeCluster(no_cores,type="SOCK")
  registerDoParallel(cl)
  clusterExport(cl,varlist=ls())
  time <- system.time({
    out <- parLapply(cl,1:s,fit_cluster_tbart)
  })
  stopCluster(cl)
  print("T-BART fit done!")
  print(time)
  # saveRDS(out,paste0("Results/tbart_",dgp,".rds"))
}
if ("sbart" %in% models)
{
  sbart.global.parmlist <- list(standardize=T,sample_sigma_global=TRUE,sigma2_global_init=0.01)
  sbart.mean.parmlist <- list(num_trees=150, min_samples_leaf=20, alpha=0.95, beta=2,
                              max_depth=20, sample_sigma2_leaf=FALSE)
  sbart.var.parmlist <- list(num_trees = 2,min_samples_leaf = 10)
  print("S-BART fit")
  print(paste(c("k1","k2","k3","k4","k5","p","rho"),args[1:7],sep=": "))
  cl <- makeCluster(no_cores,type="SOCK")
  registerDoParallel(cl)
  clusterExport(cl,varlist=ls())
  time <- system.time({
    out <- parLapply(cl,1:s,fit_cluster_sbart)
  })
  stopCluster(cl)
  print("S-BART fit done!")
  print(time)
  # saveRDS(out,paste0("Results/sbart_",dgp,".rds"))
}
if ("polynomial" %in% models)
{
  print("Polynomial fit")
  print(paste(c("k1","k2","k3","k4","k5","p","rho"),args[1:7],sep=": "))
  cl <- makeCluster(no_cores,type="SOCK")
  registerDoParallel(cl)
  clusterExport(cl,varlist=ls())
  time <- system.time({
    out <- parLapply(cl,1:s,fit_cluster_polynomial)
  })
  stopCluster(cl)
  print("Polynomial fit done!")
  print(time)
  # saveRDS(out,paste0("Results/polynomial_",dgp,".rds"))
}