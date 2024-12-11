## Setup
set.seed(0)
setwd("../")
if (!dir.exists("Results")) dir.create("Results") ## Create results folder
if (length(list.files("Results")[grep("llr1_",list.files("Results"))])!=0) ## Clean up folder
{
  files <- paste0("Results/",list.files("Results")[grep("llr1_",list.files("Results"))])
  for (i in files) file.remove(i)
}
### Functions
fit <- function(i)
{
  print(paste0("Sample: ",i))
  ifelse(is.list(data$w),ws <- data$w[[i]],ws <- subset(data$w,select=i))
  rdrobust::rdrobust(data$y[,i],data$x[,i],c=c,covs=ws)
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
  saveRDS(list(results=out,time=time),paste0("Results/llr1_",i,".rds"))
  rm(out)
  gc()
}
##
print("Done!")