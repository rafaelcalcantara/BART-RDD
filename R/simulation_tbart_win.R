## Setup
set.seed(0)
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
  ws <- as.matrix(data$w[,i])
  xs <- data$x[,i]
  zs <- data$z[,i]
  fit1 <- XBART::XBART(ys[zs==1], cbind(xs,ws)[zs==1,], num_trees = ntrees,
                num_cutpoints = sum(zs==1), num_sweeps = num_sweeps,
                burnin = burnin, Nmin = Nmin,
                p_categorical = p_categorical, max_depth = max_depth,
                tau = var(ys[zs==1])/ntrees, parallel=F)
  fit0 <- XBART::XBART(ys[zs==0], cbind(xs,ws)[zs==0,], num_trees = ntrees,
                num_cutpoints = sum(zs==0), num_sweeps = num_sweeps,
                burnin = burnin, Nmin = Nmin,
                p_categorical = p_categorical, max_depth = max_depth,
                tau = var(ys[zs==0])/ntrees, parallel=F)
  test <- -Owidth+c<=xs & xs<=Owidth+c
  test.sample <- cbind(c,ws)[test,]
  pred1 <- XBART::predict.XBART(fit1,test.sample)[,(burnin+1):num_sweeps]
  pred0 <- XBART::predict.XBART(fit0,test.sample)[,(burnin+1):num_sweeps]
  pred1-pred0
}
##
### BEGIN LOOP
files <- length(list.files("Data"))
s0 <- 1
for (i in 1:files)
{
  print(paste0("DGP: ",i))
  data <- readRDS(paste0("Data/dgp_",i,".rds"))
  if (paste0("Results/tbart_",i,".rds") %in% list.files("Results/") == T)
  {
    res <- readRDS(paste0("Results/tbart_",i,".rds"))
  } else
  {
    res <- list(results=vector("list",s))
  }
  n <- data$n
  lvl <- data$level
  Owidth <- Ow[lvl]
  s <- ncol(data$y)
  s1 <- s
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
  saveRDS(list(results=res$results,time=time/(s1-s0+1)),paste0("Results/tbart_",i,".rds"))
  rm(out)
  gc()
}
##
print("Done!")