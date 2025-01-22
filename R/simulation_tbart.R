## Setup
set.seed(0)
### Functions
fit <- function(i)
{
  print(paste0("Sample: ",i))
  Owidth <- data$h[i]
  ys <- data$y[,i]
  ws <- as.matrix(data$w[[i]])
  xs <- data$x[,i]
  zs <- data$z[,i]
  tbart.fit.0 = stochtree::bart(X_train= as.matrix(cbind(xs,ws)[zs==0,]), y_train=ys[zs==0],
                                mean_forest_params=mean.parmlist,
                                general_params=global.parmlist,num_mcmc=1000,num_gfr=30)
  tbart.fit.1 = stochtree::bart(X_train= as.matrix(cbind(xs,ws)[zs==1,]), y_train=ys[zs==1],
                                mean_forest_params=mean.parmlist,
                                general_params=global.parmlist,num_mcmc=1000,num_gfr=30)
  test <- -Owidth+c<=xs & xs<=Owidth+c
  xmat_test <- cbind(c,ws)[test,]
  pred1 <- predict(tbart.fit.1,xmat_test)$y_hat
  pred0 <- predict(tbart.fit.0,xmat_test)$y_hat
  pred1-pred0
}
##
### BEGIN LOOP
for (i in files)
{
  print(paste0("DGP: ",i))
  data <- readRDS(paste0("Data/dgp_",i,".rds"))
  s <- ncol(data$y)
  if (paste0("tbart_",i,".rds") %in% list.files("Results/") == T)
  {
    res <- readRDS(paste0("Results/tbart_",i,".rds"))
  } else
  {
    res <- list(results=vector("list",s))
  }
  n <- data$n
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
  saveRDS(list(results=res$results,time=time),paste0("Results/tbart_",i,".rds"))
  rm(out)
  gc()
}
##
print("Done!")