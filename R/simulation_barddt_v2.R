## Setup
set.seed(0)
### Fit function
fit <- function(i)
{
  print(paste0("Sample: ",i))
  ys <- data$y[,i]
  ws <- as.matrix(data$w[[i]])
  xs <- data$x[,i]
  zs <- data$z[,i]
  Owidth <- rdrobust::rdbwselect(ys,xs,c)$bws[2]
  B <- cbind(zs*xs, xs, zs,rep(1,n))
  s <- 1/sqrt(sum(apply(B,2,var)))
  B <- s*B
  test <- -Owidth+c<=xs & xs<=Owidth+c
  B1 <- s*cbind(rep(c,n), rep(c,n), rep(1,n), rep(1,n))
  B0 <- s*cbind(rep(0,n), rep(c,n), rep(0,n), rep(1,n))
  barddt.fit = stochtree::bart(X_train= as.matrix(cbind(xs,ws)), y_train=ys,
                               W_train = B, mean_forest_params=mean.parmlist,
                               general_params=global.parmlist,
                               num_mcmc=1000,num_gfr=30)
  B1 <- B1[test,]
  B0 <- B0[test,]
  xmat_test <- cbind(rep(c,n),ws)[test,]
  pred1 <- predict(barddt.fit,xmat_test,B1)$y_hat
  pred0 <- predict(barddt.fit,xmat_test,B0)$y_hat
  pred1 - pred0
}
### BEGIN LOOP
for (i in files)
{
  print(paste0("DGP: ",i))
  data <- readRDS(paste0("Data/dgp_",i,".rds"))
  s <- ncol(data$y)
  if (paste0("bart_rdd_",i,".rds") %in% list.files("Results/") == T)
  {
    res <- readRDS(paste0("Results/bart_rdd_",i,".rds"))
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
  saveRDS(list(results=res$results,time=time),paste0("Results/bart_rdd_",i,".rds"))
  rm(out)
  gc()
}
##
print("Done!")
