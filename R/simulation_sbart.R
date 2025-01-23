## Setup
set.seed(0)
### Functions
fit <- function(i)
{
  print(paste0("Sample: ",i))
  ys <- data$y[,i]
  ws <- as.matrix(data$w[[i]])
  xs <- data$x[,i]
  zs <- data$z[,i]
  Owidth <- rdrobust::rdbwselect(ys,xs,c)$bws[2]
  sbart.fit = stochtree::bart(X_train= as.matrix(cbind(xs,ws,zs)), y_train=ys,
                                mean_forest_params=mean.parmlist,
                                general_params=global.parmlist,num_mcmc=1000,num_gfr=30)
  test <- -Owidth+c<=xs & xs<=Owidth+c
  xmat_test_1 <- cbind(c,ws,1)[test,]
  xmat_test_0 <- cbind(c,ws,0)[test,]
  pred1 <- predict(sbart.fit,xmat_test_1)$y_hat
  pred0 <- predict(sbart.fit,xmat_test_0)$y_hat
  pred1-pred0
}
### BEGIN LOOP
for (i in files)
{
  print(paste0("DGP: ",i))
  data <- readRDS(paste0("Data/dgp_",i,".rds"))
  s <- ncol(data$y)
  if (paste0("sbart_",i,".rds") %in% list.files("Results/") == T)
  {
    res <- readRDS(paste0("Results/sbart_",i,".rds"))
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
  saveRDS(list(results=res$results,time=time),paste0("Results/sbart_",i,".rds"))
  rm(out)
  gc()
}
##
print("Done!")