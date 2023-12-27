## Setup
set.seed(0)
library(doParallel)
## devtools::install_github("JingyuHe/XBART@XBCF-RDD")
library(XBART)
setwd("~/../Git/BART-RDD")
if (!dir.exists("Results")) dir.create("Results") ## Create results folder
if (length(list.files("Results")[grep("sbart_",list.files("Results"))])!=0) ## Clean up folder
{
  files <- paste0("Results/",list.files("Results")[grep("sbart_",list.files("Results"))])
  for (i in files) file.remove(i)
}
### Parallelization
no_cores <- detectCores()-1
### Parameters
ntrees        <- 10
Nmin          <- 5
num_sweeps    <- 120
burnin        <- 20
p_categorical <- 2
### Functions
fit <- function(i)
{
  print(paste0("Sample: ",i))
  ys <- data$y[,i]
  ws <- data$w[[i]]
  xs <- data$x[,i]
  zs <- data$z[,i]
  fit <- XBART::XBART(ys, cbind(xs,ws,zs), num_trees = ntrees,
                      num_cutpoints = n, num_sweeps = num_sweeps,
                      burnin = burnin, Nmin = Nmin,
                      p_categorical = p_categorical,
                      tau = var(ys)/ntrees, parallel=F)
  pred1 <- XBART::predict.XBART(fit,cbind(rep(0,n),ws,rep(1,n)))[,(burnin+1):num_sweeps]
  pred0 <- XBART::predict.XBART(fit,cbind(rep(0,n),ws,rep(0,n)))[,(burnin+1):num_sweeps]
  pred1-pred0
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
  h <- apply(data$x,2,sd)*0.5
  cl <- makeCluster(no_cores,type="SOCK")
  registerDoParallel(cl)
  clusterExport(cl,varlist=ls())
  time <- system.time({
    out <- parLapply(cl,1:s,fit)
  })
  stopCluster(cl)
  print(time)
  saveRDS(list(results=out,time=time),paste0("Results/sbart_",i,".rds"))
  rm(out)
  gc()
}
##
print("Done!")