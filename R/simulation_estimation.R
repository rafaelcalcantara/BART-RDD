# Functions that run the models
## BARDDT
fit.barddt <- function(y,x,w,z)
{
  barddt.global.parmlist <- list(standardize=T,sample_sigma_global=TRUE,sigma2_global_init=0.01)
  barddt.mean.parmlist <- list(num_trees=150, min_samples_leaf=20, alpha=0.95, beta=2,
                               max_depth=20, sample_sigma2_leaf=FALSE)
  barddt.var.parmlist <- list(num_trees = 2,min_samples_leaf = 10)
  B <- cbind(z*x,(1-z)*x, z,rep(1,n))
  test <- -Owidth+c<=x & x<=Owidth+c
  B1 <- cbind(rep(c,n), rep(0,n), rep(1,n), rep(1,n))
  B0 <- cbind(rep(0,n), rep(c,n), rep(0,n), rep(1,n))
  s <- 1/sqrt(sum(apply(B,2,var)))
  B <- s*B
  barddt.fit = stochtree::bart(X_train= as.matrix(cbind(x,w)), y_train=y,
                               W_train = B, mean_forest_params=barddt.mean.parmlist,
                               general_params=barddt.global.parmlist,
                               variance_forest_params=barddt.var.parmlist,
                               num_mcmc=1000,num_gfr=30)
  B1 <- B1[test,]
  B0 <- B0[test,]
  xmat_test <- cbind(rep(0,n),w)[test,]
  pred1 <- predict(barddt.fit,xmat_test,B1)$y_hat
  pred0 <- predict(barddt.fit,xmat_test,B0)$y_hat
  return(list(pred1=pred1,pred0=pred0))
}
## T-BART
fit.tbart <- function(y,x,w,z)
{
  tbart.global.parmlist <- list(standardize=T,sample_sigma_global=TRUE,sigma2_global_init=0.01)
  tbart.mean.parmlist <- list(num_trees=150, min_samples_leaf=20, alpha=0.95, beta=2,
                              max_depth=20, sample_sigma2_leaf=FALSE)
  tbart.var.parmlist <- list(num_trees = 2,min_samples_leaf = 10)
  test <- -Owidth+c<=x & x<=Owidth+c
  tbart.fit.0 = stochtree::bart(X_train= as.matrix(cbind(x,w)[z==0,]), y_train=y[z==0],
                                mean_forest_params=tbart.mean.parmlist,
                                general_params=tbart.global.parmlist,
                                variance_forest_params=tbart.var.parmlist,
                                num_mcmc=1000,num_gfr=30)
  tbart.fit.1 = stochtree::bart(X_train= as.matrix(cbind(x,w)[z==1,]), y_train=y[z==1],
                                mean_forest_params=tbart.mean.parmlist,
                                general_params=tbart.global.parmlist,
                                variance_forest_params=tbart.var.parmlist,
                                num_mcmc=1000,num_gfr=30)
  xmat_test <- cbind(c,w)[test,]
  pred1 <- predict(tbart.fit.1,xmat_test)$y_hat
  pred0 <- predict(tbart.fit.0,xmat_test)$y_hat
  return(list(pred1=pred1,pred0=pred0))
}
## S-BART
fit.sbart <- function(y,x,w,z)
{
  sbart.global.parmlist <- list(standardize=T,sample_sigma_global=TRUE,sigma2_global_init=0.01)
  sbart.mean.parmlist <- list(num_trees=150, min_samples_leaf=20, alpha=0.95, beta=2,
                              max_depth=20, sample_sigma2_leaf=FALSE)
  sbart.var.parmlist <- list(num_trees = 2,min_samples_leaf = 10)
  test <- -Owidth+c<=x & x<=Owidth+c
  sbart.fit = stochtree::bart(X_train= as.matrix(cbind(x,z,w)), y_train=y,
                              mean_forest_params=sbart.mean.parmlist,
                              general_params=sbart.global.parmlist,
                              variance_forest_params=sbart.var.parmlist,,
                              num_mcmc=1000,num_gfr=30)
  xmat_test.1 <- cbind(c,1,w)[test,]
  xmat_test.0 <- cbind(c,0,w)[test,]
  pred1 <- predict(sbart.fit,xmat_test.1)$y_hat
  pred0 <- predict(sbart.fit,xmat_test.0)$y_hat
  return(list(pred1=pred1,pred0=pred0))
}
## Polynomial
fit.polynomial <- function(y,x,w,z)
{
  test <- -Owidth+c<=x & x<=Owidth+c
  dfw <- data.frame(w=w)
  fmla <- as.formula(paste('y~(',paste(paste('poly(', names(dfw), ', 4)', sep=''),collapse="+"),')*poly(x,1)*z + poly(x,3)'))
  df <- data.frame(x=x,w=w,y=y,z=z)
  df$z <- as.factor(df$z)
  h <- rdrobust::rdbwselect(y,x,c)$bws[1]
  df.train <- subset(df,c-h<=x & x<c+h)
  poly.fit <- lm(fmla,data = df.train)
  df.test <- df[test,]
  xmat_test.1 <- xmat_test.0 <- df.test
  xmat_test.1$x <- c
  xmat_test.1$z <- "1"
  xmat_test.0$x <- c
  xmat_test.0$z <- "0"
  pred1 <- predict(poly.fit,xmat_test.1)
  pred0 <- predict(poly.fit,xmat_test.0)
  return(list(pred1=pred1,pred0=pred0))
}
# Parallelized functions
## Run regressions from R
fit_r <- function(model)
{
  if (model=="leaf.rdd")
  {
    xs <- x[,1]
    ys <- y[,1]
    zs <- z[,1]
    out <- fit.barddt(ys,xs,w,zs)
  } else if (model=="tbart")
  {
    xs <- x[,1]
    ys <- y[,1]
    zs <- z[,1]
    out <- fit.tbart(ys,xs,w,zs)
  } else if (model=="sbart")
  {
    xs <- x[,1]
    ys <- y[,1]
    zs <- z[,1]
    out <- fit.sbart(ys,xs,w,zs)
  } else if (model=="polynomial")
  {
    xs <- x[,1]
    ys <- y[,1]
    zs <- z[,1]
    out <- fit.polynomial(ys,xs,w,zs)
  }
  out$pred1 - out$pred0
}
## Run regressions from cluster
fit_cluster_barddt <- function(sample)
{
  xs <- x[,sample]
  ys <- y[,sample]
  zs <- z[,sample]
  out <- fit.barddt(ys,xs,w,zs)
  out$pred1 - out$pred0
}
fit_cluster_tbart <- function(sample)
{
  xs <- x[,sample]
  ys <- y[,sample]
  zs <- z[,sample]
  out <- fit.tbart(ys,xs,w,zs)
  out$pred1 - out$pred0
}
fit_cluster_sbart <- function(sample)
{
  xs <- x[,sample]
  ys <- y[,sample]
  zs <- z[,sample]
  out <- fit.sbart(ys,xs,w,zs)
  out$pred1 - out$pred0
}
fit_cluster_polynomial <- function(sample)
{
  xs <- x[,sample]
  ys <- y[,sample]
  zs <- z[,sample]
  out <- fit.polynomial(ys,xs,w,zs)
  out$pred1 - out$pred0
}
## Run the models
if (isTRUE(cmd.line))
{
  ## Create/clean results folder
  if (!dir.exists("Results")) dir.create("Results") ## Create results folder, if non-existent
  if (length(list.files("Results"))!=0) ## Clean up folder
  {
    files <- paste0("Results/",list.files("Results"))
    for (i in files) file.remove(i)
  }
  ## Fit from cluster
  if ("leaf.rdd" %in% models)
  {
    print("BARDDT fit")
    print(paste(c("k1","k2","k3","k4","k5","p","rho"),args[1:7],sep=": "))
    cl <- makeCluster(no_cores,type="SOCK")
    registerDoParallel(cl)
    clusterExport(cl,varlist=ls())
    time <- system.time({
      out <- parLapply(cl,1:s,fit_cluster_barddt)
    })
    stopCluster(cl)
    print("BARDDT fit done!")
    print(time)
    saveRDS(out,paste0("Results/barddt_",dgp,".rds"))
  }
  if ("tbart" %in% models)
  {
    print("T-BART fit")
    print(paste(c("k1","k2","k3","k4","k5","p","rho"),args[1:7],sep=": "))
    cl <- makeCluster(no_cores,type="SOCK")
    registerDoParallel(cl)
    clusterExport(cl,varlist=ls())
    time <- system.time({
      out <- parLapply(cl,1:s,fit_cluster_tbart)
    })
    stopCluster(cl)
    print("T-BART fit done!")
    print(time)
    saveRDS(out,paste0("Results/tbart_",dgp,".rds"))
  }
  if ("sbart" %in% models)
  {
    print("T-BART fit")
    print(paste(c("k1","k2","k3","k4","k5","p","rho"),args[1:7],sep=": "))
    cl <- makeCluster(no_cores,type="SOCK")
    registerDoParallel(cl)
    clusterExport(cl,varlist=ls())
    time <- system.time({
      out <- parLapply(cl,1:s,fit_cluster_sbart)
    })
    stopCluster(cl)
    print("T-BART fit done!")
    print(time)
    saveRDS(out,paste0("Results/sbart_",dgp,".rds"))
  }
  if ("polynomial" %in% models)
  {
    print("Polynomial fit")
    print(paste(c("k1","k2","k3","k4","k5","p","rho"),args[1:7],sep=": "))
    cl <- makeCluster(no_cores,type="SOCK")
    registerDoParallel(cl)
    clusterExport(cl,varlist=ls())
    time <- system.time({
      out <- parLapply(cl,1:s,fit_cluster_polynomial)
    })
    stopCluster(cl)
    print("Polynomial fit done!")
    print(time)
    saveRDS(out,paste0("Results/polynomial_",dgp,".rds"))
  }
} else
{
  ## Fit from R
  cl <- makeCluster(no_cores,type="SOCK")
  registerDoParallel(cl)
  clusterExport(cl,varlist=ls())
  time <- system.time({
    out <- parLapply(cl,models,fit_r)
  })
  stopCluster(cl)
  print(time)
  names(out) <- models
  ### Visualizing
  test <- c-Owidth <= x[,1] & x[,1] <= c+Owidth
  cate <- cate[test,]
  matplot(w[test,1],cbind(cate,rowMeans(out$leaf.rdd),rowMeans(out$tbart)),
          col=c("black","orange","maroon"),pch=19,cex=0.8)
  matplot(w[test,1],cbind(cate,rowMeans(out$leaf.rdd),rowMeans(out$sbart)),
          col=c("black","orange","red"),pch=19,cex=0.8)
  matplot(w[test,1],cbind(cate,rowMeans(out$leaf.rdd),out$polynomial),
          col=c("black","orange","pink"),pch=19,cex=0.8)
  
  rmse.ate <- sqrt(mean((ate-cate)^2))
  rmse.barddt <- sqrt(mean((rowMeans(out$leaf.rdd)-cate)^2))
  rmse.tbart <- sqrt(mean((rowMeans(out$tbart)-cate)^2))
  rmse.sbart <- sqrt(mean((rowMeans(out$sbart)-cate)^2))
  rmse.polynomial <- sqrt(mean((out$polynomial-cate)^2))
  print("RMSE BARDDT:")
  print(rmse.barddt/rmse.ate)
  print("RMSE TBART:")
  print(rmse.tbart/rmse.ate)
  print("RMSE SBART:")
  print(rmse.sbart/rmse.ate)
  print("RMSE Polynomial:")
  print(rmse.polynomial/rmse.ate)
}