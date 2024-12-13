## Setup
set.seed(0)
### Functions
fit <- function(i)
{
  print(paste0("Sample: ",i))
  ys <- data$y[,i]
  ws <- as.matrix(data$w[,i])
  xs <- data$x[,i]
  zs <- data$z[,i]
  fit <- XBART::XBART(ys, cbind(xs,ws,zs), num_trees = ntrees,
                      num_cutpoints = n, num_sweeps = num_sweeps,
                      burnin = burnin, p_categorical = p_categorical,
                      tau = var(ys)/ntrees, parallel=F)
  test <- -Owidth+c<=xs & xs<=Owidth+c
  test.sample.0 <- cbind(c,ws,0)[test,]
  test.sample.1 <- cbind(c,ws,1)[test,]
  pred1 <- XBART::predict.XBART(fit,test.sample.1)[,(burnin+1):num_sweeps]
  pred0 <- XBART::predict.XBART(fit,test.sample.0)[,(burnin+1):num_sweeps]
  pred1-pred0
}
##
### BEGIN LOOP
files <- length(list.files("Data"))
for (i in 1:files)
{
  print(paste0("DGP: ",i))
  data <- readRDS(paste0("Data/dgp_",i,".rds"))
  if (paste0("Results/sbart_",i,".rds") %in% list.files("Results/") == T)
  {
    res <- readRDS(paste0("Results/sbart_",i,".rds"))
  } else
  {
    res <- list(results=vector("list",s))
  }
  n <- data$n
  s <- ncol(data$y)
  if (n==500)
  {
    Owidth <- Ow[1]
  } else if (n==1000)
  {
    Owidth <- Ow[2]
  } else if (n==2500)
  {
    Owidth <- Ow[3]
  } else
  {
    Owidth <- Ow[4]
  }
  c <- data$c
  cl <- makeCluster(no_cores,type="SOCK")
  registerDoParallel(cl)
  clusterExport(cl,varlist=ls())
  time <- system.time({
    out <- parLapply(cl,s0:s1,fit)
  })
  stopCluster(cl)
  print(time)
  res$results[s0:s1] <- out
  saveRDS(list(results=res$results,time=time/(s1-s0+1)),paste0("Results/sbart_",i,".rds"))
  rm(out)
  gc()
}
##
print("Done!")