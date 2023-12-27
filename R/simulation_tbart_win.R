## Setup
set.seed(0)
library(doParallel)
## devtools::install_github("JingyuHe/XBART@XBCF-RDD")
library(XBART)
setwd("~/../Git/BART-RDD")
if (!dir.exists("Results")) dir.create("Results") ## Create results folder
if (length(list.files("Results")[grep("tbart_",list.files("Results"))])!=0) ## Clean up folder
{
  files <- paste0("Results/",list.files("Results")[grep("tbart_",list.files("Results"))])
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
  fit1 <- XBART::XBART(ys[zs==1], cbind(xs,ws)[zs==1,], num_trees = ntrees,
                num_cutpoints = sum(zs==1), num_sweeps = num_sweeps,
                burnin = burnin, Nmin = Nmin,
                p_categorical = p_categorical,
                tau = var(ys[zs==1])/ntrees, parallel=F)
  fit0 <- XBART::XBART(ys[zs==0], cbind(xs,ws)[zs==0,], num_trees = ntrees,
                num_cutpoints = sum(zs==0), num_sweeps = num_sweeps,
                burnin = burnin, Nmin = Nmin,
                p_categorical = p_categorical,
                tau = var(ys[zs==0])/ntrees, parallel=F)
  pred1 <- XBART::predict.XBART(fit1,cbind(rep(0,sum(zs==1)),ws[zs==1,]))[,(burnin+1):num_sweeps]
  pred0 <- XBART::predict.XBART(fit0,cbind(rep(0,sum(zs==0)),ws[zs==0,]))[,(burnin+1):num_sweeps]
  colMeans(pred1)-colMeans(pred0)
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
  saveRDS(list(results=out,time=time),paste0("Results/tbart_",i,".rds"))
  rm(out)
  gc()
}
##
print("Done!")