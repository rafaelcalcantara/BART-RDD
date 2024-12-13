## Setup
set.seed(0)
library(doParallel)
setwd("~/../Git/BART-RDD")
if (!dir.exists("Results")) dir.create("Results") ## Create results folder
if (length(list.files("Results")[grep("llr0_",list.files("Results"))])!=0) ## Clean up folder
{
  files <- paste0("Results/",list.files("Results")[grep("llr0_",list.files("Results"))])
  for (i in files) file.remove(i)
}
### Parallelization
no_cores <- detectCores()-1
### Functions
fit <- function(i)
{
  print(paste0("Sample: ",i))
  rdrobust::rdrobust(data$y[,i],data$x[,i],c=c)
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
  cl <- makeCluster(no_cores,type="SOCK")
  registerDoParallel(cl)
  clusterExport(cl,varlist=ls())
  time <- system.time({
    out <- parLapply(cl,1:s,fit)
  })
  stopCluster(cl)
  print(time)
  saveRDS(list(results=out,time=time),paste0("Results/llr0_",i,".rds"))
  rm(out)
  gc()
}
##
print("Done!")