# Functions that run the models
## BARDDT
fit.barddt <- function(y,x,w,z,test,c)
{
  barddt.global.parmlist <- list(standardize=T,sample_sigma_global=TRUE,sigma2_global_init=0.1)
  barddt.mean.parmlist <- list(num_trees=150, min_samples_leaf=20, alpha=0.95, beta=2,
                               max_depth=20, sample_sigma2_leaf=FALSE, sigma2_leaf_init = diag(rep(0.1/150,4)))
  # barddt.var.parmlist <- list(num_trees = 2,min_samples_leaf = 10)
  B <- cbind(z*x,(1-z)*x, z,rep(1,n))
  B1 <- cbind(rep(c,n), rep(0,n), rep(1,n), rep(1,n))
  B0 <- cbind(rep(0,n), rep(c,n), rep(0,n), rep(1,n))
  barddt.fit = stochtree::bart(X_train= as.matrix(cbind(x,w)), y_train=y,
                               W_train = B, mean_forest_params=barddt.mean.parmlist,
                               general_params=barddt.global.parmlist,
                               # variance_forest_params=barddt.var.parmlist,
                               num_mcmc=100,num_gfr=30)
  B1 <- B1[test,]
  B0 <- B0[test,]
  xmat_test <- as.matrix(cbind(rep(0,n),w)[test,])
  pred1 <- predict(barddt.fit,xmat_test,B1)$y_hat
  pred0 <- predict(barddt.fit,xmat_test,B0)$y_hat
  return(pred1-pred0)
}
## T-BART
fit.tbart <- function(y,x,w,z,test,c)
{
  tbart.global.parmlist <- list(standardize=T,sample_sigma_global=TRUE,sigma2_global_init=0.01)
  tbart.mean.parmlist <- list(num_trees=150, min_samples_leaf=20, alpha=0.95, beta=2,
                              max_depth=20, sample_sigma2_leaf=TRUE)
  # tbart.var.parmlist <- list(num_trees = 2,min_samples_leaf = 10)
  tbart.fit.0 = stochtree::bart(X_train= as.matrix(cbind(x,w)[z==0,]), y_train=y[z==0],
                                mean_forest_params=tbart.mean.parmlist,
                                general_params=tbart.global.parmlist,
                                # variance_forest_params=tbart.var.parmlist,
                                num_mcmc=100,num_gfr=30)
  tbart.fit.1 = stochtree::bart(X_train= as.matrix(cbind(x,w)[z==1,]), y_train=y[z==1],
                                mean_forest_params=tbart.mean.parmlist,
                                general_params=tbart.global.parmlist,
                                # variance_forest_params=tbart.var.parmlist,
                                num_mcmc=100,num_gfr=30)
  xmat_test <- as.matrix(cbind(c,w)[test,])
  pred1 <- predict(tbart.fit.1,xmat_test)$y_hat
  pred0 <- predict(tbart.fit.0,xmat_test)$y_hat
  return(pred1-pred0)
}
## S-BART
fit.sbart <- function(y,x,w,z,test,c)
{
  sbart.global.parmlist <- list(standardize=T,sample_sigma_global=TRUE,sigma2_global_init=0.01)
  sbart.mean.parmlist <- list(num_trees=150, min_samples_leaf=20, alpha=0.95, beta=2,
                              max_depth=20, sample_sigma2_leaf=TRUE)
  # sbart.var.parmlist <- list(num_trees = 2,min_samples_leaf = 10)
  sbart.fit = stochtree::bart(X_train= as.matrix(cbind(x,z,w)), y_train=y,
                              mean_forest_params=sbart.mean.parmlist,
                              general_params=sbart.global.parmlist,
                              # variance_forest_params=sbart.var.parmlist,,
                              num_mcmc=100,num_gfr=30)
  xmat_test.1 <- as.matrix(cbind(c,1,w)[test,])
  xmat_test.0 <- as.matrix(cbind(c,0,w)[test,])
  pred1 <- predict(sbart.fit,xmat_test.1)$y_hat
  pred0 <- predict(sbart.fit,xmat_test.0)$y_hat
  return(pred1-pred0)
}
## Polynomial
fit.polynomial <- function(y,x,w,z,h,test,c)
{
  dfw <- data.frame(w=w)
  fmla <- as.formula(paste('y~(',paste(paste('poly(', names(dfw), ', 4)', sep=''),collapse="+"),')*poly(x,1)*z + poly(x,3)'))
  df <- data.frame(x=x,w=w,y=y,z=z)
  df$z <- as.factor(df$z)
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
  return(pred1-pred0)
}
fit.ate <- function(y,x)
{
  return(rdrobust::rdrobust(y,x,c))
}
# Function to run all models for 1 sample
fit_general <- function(sample)
{
  set.seed(sample)
  source("simulation_data.R")
  ## Saving data
  saveRDS(list(y=y,x=x,z=z,w=w,cate=cate),paste0("Data/dgp_",dgp,"_sample_",sample,".rds"))
  ###
  w <- read.table(paste0("Data/w_",p,".csv"), sep = ",")
  ate <- fit.ate(y,x)
  h <- ate$bws[2,2]
  ate <- ate$coef[3]
  test <- -Owidth+c<=x & x<=Owidth+c
  write.table(ate,paste0("Results/",dgp,"/ate_sample_",sample,".csv"), row.names = FALSE, col.names = FALSE, sep = ",")
  time.barddt <- system.time({
    pred.barddt <- fit.barddt(y,x,w,z,test,c)
    write.table(pred.barddt,paste0("Results/",dgp,"/barddt_sample_",sample,".csv"), row.names = FALSE, col.names = FALSE, sep = ",")
  })
  time.tbart <- system.time({
    pred.tbart <- fit.tbart(y,x,w,z,test,c)
    write.table(pred.tbart,paste0("Results/",dgp,"/tbart_sample_",sample,".csv"), row.names = FALSE, col.names = FALSE, sep = ",")
  })
  time.sbart <- system.time({
    pred.sbart <- fit.sbart(y,x,w,z,test,c)
    write.table(pred.sbart,paste0("Results/",dgp,"/sbart_sample_",sample,".csv"), row.names = FALSE, col.names = FALSE, sep = ",")
  })
  time.polynomial <- system.time({
    pred.polynomial <- fit.polynomial(y,x,w,z,h,test,c)
    write.table(pred.polynomial,paste0("Results/",dgp,"/polynomial_sample_",sample,".csv"), row.names = FALSE, col.names = FALSE, sep = ",")
  })
  writeLines(c(dim(pred.barddt),dim(pred.tbart),dim(pred.sbart),length(pred.polynomial)),paste0("Logs/",dgp,"/sizes_",sample,".txt"))
  if (nrow(pred.barddt)!=nrow(pred.tbart) | nrow(pred.barddt)!=nrow(pred.sbart) | nrow(pred.barddt)!=length(pred.polynomial)
      | nrow(pred.tbart)!=nrow(pred.sbart) | nrow(pred.tbart)!=length(pred.polynomial) | nrow(pred.sbart)!=length(pred.polynomial))
  {
    stop("Dim mismatch")
  }
  # write.table(c(time.barddt[3],sample),paste0("Time/",dgp,"/barddt.csv"), append=TRUE, row.names = FALSE, col.names = FALSE, sep = ",")
  # write.table(c(time.tbart[3],sample),paste0("Time/",dgp,"/tbart.csv"), append=TRUE, row.names = FALSE, col.names = FALSE, sep = ",")
  # write.table(c(time.sbart[3],sample),paste0("Time/",dgp,"/sbart.csv"), append=TRUE, row.names = FALSE, col.names = FALSE, sep = ",")
  # write.table(c(time.polynomial[3],sample),paste0("Time/",dgp,"/polynomial.csv"), append=TRUE, row.names = FALSE, col.names = FALSE, sep = ",")
  # write.table(c(time.barddt[3]+time.tbart[3]+time.sbart[3]+time.polynomial[3],sample),paste0("Time/",dgp,"/total_per_sample.csv"), append=TRUE, row.names = FALSE, col.names = FALSE, sep = ",")
  # Processing results
  # calc.rmse(sample,ate,test)
  # point.est(sample,test)
}