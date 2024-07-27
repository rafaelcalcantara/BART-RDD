## Setup
set.seed(0)
## devtools::install_github("JingyuHe/XBART@XBCF-RDD")
library(XBART)
setwd("~/../Git/BART-RDD")
if (!dir.exists("Results")) dir.create("Results") ## Create results folder
if (length(list.files("Results")[grep("bcf_",list.files("Results"))])!=0) ## Clean up folder
{
  files <- paste0("Results/",list.files("Results")[grep("bcf_",list.files("Results"))])
  for (i in files) file.remove(i)
}
### Parameters
ntrees        <- 5
Nmin          <- 5
num_sweeps    <- 120
burnin        <- 20
### Functions
fit <- function(i)
{
  print(paste0("Sample: ",i))
  ys <- data$y[,i]
  ifelse(is.list(data$w),ws <- data$w[[i]],ws <- subset(data$w,select=i))
  xs <- data$x[,i]
  zs <- data$z[,i]
  fit <- XBART::XBCF.discrete(y=ys, Z=zs, X_con = cbind(xs,ws), X_mod = cbind(xs,ws),
                              pihat=zs,
                              num_trees_mod = ntrees,
                              num_trees_con = ntrees,
                              num_cutpoints = n,
                              num_sweeps = num_sweeps,
                              burnin = burnin, Nmin = Nmin,
                              p_categorical_con = p_categorical,
                              p_categorical_mod = p_categorical,
                              tau_con = 2*var(ys)/ntrees,
                              tau_mod = 0.5*var(ys)/ntrees)
  test <- -Owidth+c<=xs & xs<=Owidth+c
  test.sample <- cbind(c,ws)[test,]
  pred <- XBART::predict.XBCFdiscrete(fit,X_con = test.sample, X_mod = test.sample,Z=zs[test],pihat=zs[test],burnin=burnin)
  pred$tau.adj[,(burnin+1):num_sweeps]
}
##
### BEGIN LOOP
files <- length(list.files("Data"))
for (i in 1:files)
{
  print(paste0("DGP: ",i))
  data <- readRDS(paste0("Data/dgp_",i,".rds"))
  n <- nrow(data$y)
  s <- ncol(data$y)
  c <- data$c
  # test <- readRDS(paste0("Data/test_dgp_",i,".rds"))
  # test.w <- test$w
  # test.sample <- cbind(c,test.w)
  cl <- makeCluster(no_cores,type="SOCK")
  registerDoParallel(cl)
  clusterExport(cl,varlist=ls())
  time <- system.time({
    out <- parLapply(cl,1:s,fit)
  })
  stopCluster(cl)
  print(time)
  saveRDS(list(results=out,time=time),paste0("Results/bcf_",i,".rds"))
  rm(out)
  gc()
}
##
print("Done!")