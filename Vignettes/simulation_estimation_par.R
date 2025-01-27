library(doParallel)
no_cores <- 10
## Common parameters
p_categorical <- 0
burnin <- 50
num_sweeps <- 150
# BARDDT parameters
Omin <- 2
Opct <- 0.1
ntrees_con <- 150
ntrees_mod <- 150
# T-BART parameters
ntrees <- 200
## Set Owidth and test set
h <- apply(x,2,function(i) h.grid(i,c,pts_in_window))
test <- sapply(1:s, function(i) c-h[i] <= x[,i] & x[,i] <= c+h[i])
cate <- sapply(1:s, function(i) cate[test[,i],i])

barddt <- tbart <- list(s)

## BARDDT fit------------------------------------------------------------------
fit.barddt <- function(i)
{
  ys <- y[,i]
  ws <- as.matrix(w[[i]])
  xs <- x[,i]
  Owidth <- h[i]
  sample <- 1000*Owidth
  train <- c-sample < xs & xs < c+sample
  fit <- XBART::XBCF.rd(ys[train], ws[train,], xs[train], c,
                        Owidth = Owidth, Omin = Omin, Opct = Opct,
                        num_cutpoints = sum(train), max_depth = 30,
                        num_trees_con = ntrees_con, num_trees_mod = ntrees_mod,
                        num_sweeps = num_sweeps, update_tau = TRUE,
                        alpha_mod = 0.95, alpha_con = 0.95, beta_mod = 1.25, beta_con = 1.25,
                        burnin = burnin, a_scaling = FALSE, b_scaling = FALSE,
                        Nmin = 1,parallel = FALSE,
                        p_categorical_con = p_categorical,
                        p_categorical_mod = p_categorical,
                        tau_con =0.5/ntrees_con, tau_mod = 0.5/ntrees_mod)
  test <- -Owidth+c<=xs & xs<=Owidth+c
  pred <- XBART::predict.XBCFrd(fit,ws[test,],rep(c,sum(test)))
  pred$tau.adj[,(burnin+1):num_sweeps]
}  
cl <- makeCluster(no_cores,type="SOCK")
registerDoParallel(cl)
clusterExport(cl,varlist=ls())
time.barddt <- system.time({
  barddt <- parLapply(cl,1:s,fit.barddt)
})
stopCluster(cl)
## T-BART----------------------------------------------------------------------
fit.tbart <- function(i)
{
  
  print(paste0("Sample: ",i))
  ys <- y[,i]
  ws <- as.matrix(w[[i]])
  xs <- x[,i]
  zs <- z[,i]
  Owidth <- h[i]
  sample <- max(abs(xs))
  train <- c-sample < xs & xs < c+sample
  ys <- ys[train]
  xs <- xs[train]
  ws <- ws[train,]
  zs <- zs[train]
  test <- -Owidth+c<=xs & xs<=Owidth+c
  test.sample <- cbind(c,ws)[test,]
  # fit1 <- XBART::XBART(ys[zs==1], cbind(xs,ws)[zs==1,], num_trees = ntrees,
  #                      num_cutpoints = sum(zs==1), num_sweeps = num_sweeps,
  #                      burnin = burnin, p_categorical = p_categorical,
  #                      tau = var(ys[zs==1])/ntrees, parallel=F)
  # fit0 <- XBART::XBART(ys[zs==0], cbind(xs,ws)[zs==0,], num_trees = ntrees,
  #                      num_cutpoints = sum(zs==0), num_sweeps = num_sweeps,
  #                      burnin = burnin, p_categorical = p_categorical,
  #                      tau = var(ys[zs==0])/ntrees, parallel=F)
  # pred1 <- XBART::predict.XBART(fit1,test.sample)[,(burnin+1):num_sweeps]
  # pred0 <- XBART::predict.XBART(fit0,test.sample)[,(burnin+1):num_sweeps]
  fit1 <- stochtree::bart(cbind(xs,ws)[zs==1,], ys[zs==1], X_test = cbind(c,ws)[test,],
                          num_gfr = 0, num_burnin = burnin, num_mcmc = num_sweeps)
  fit0 <- stochtree::bart(cbind(xs,ws)[zs==0,], ys[zs==0], X_test = cbind(c,ws)[test,],
                          num_gfr = 0, num_burnin = burnin, num_mcmc = num_sweeps)
  pred1 <- fit1$y_hat_test[,(burnin+1):num_sweeps]
  pred0 <- fit0$y_hat_test[,(burnin+1):num_sweeps]
  pred1-pred0
}

cl <- makeCluster(no_cores,type="SOCK")
registerDoParallel(cl)
clusterExport(cl,varlist=ls())
time.tbart <- system.time({
  tbart <- parLapply(cl,1:s,fit.tbart)
})
stopCluster(cl)
## BCF-------------------------------------------------------------------------
# fit.bcf <- function(i)
# {
#   
#   print(paste0("Sample: ",i))
#   ys <- y[,i]
#   ws <- as.matrix(w[[i]])
#   xs <- x[,i]
#   zs <- z[,i]
#   Owidth <- h[i]
#   sample <- 10*Owidth
#   train <- c-sample < xs & xs < c+sample
#   ys <- ys[train]
#   xs <- xs[train]
#   ws <- ws[train,]
#   zs <- zs[train]
#   test <- -Owidth+c<=xs & xs<=Owidth+c
#   test.sample <- cbind(c,ws)[test,]
#   # fit <- XBART::XBCF.discrete(ys, zs, cbind(xs,ws), cbind(xs,ws),
#   #                             pihat = zs,
#   #                             num_trees_con = ntrees_con,
#   #                             num_trees_mod = ntrees_mod,
#   #                             num_cutpoints = sum(train), num_sweeps = num_sweeps,
#   #                             burnin = burnin, p_categorical_con = p_categorical,
#   #                             p_categorical_mod = p_categorical,
#   #                             tau_con = var(ys)/ntrees, tau_mod = var(ys)/ntrees,
#   #                             parallel=F)
#   # pred <- XBART::predict.XBCFdiscrete(fit,test.sample,test.sample,zs[test],pihat = zs[test])
#   # pred$tau.adj[,(burnin+1):num_sweeps]
#   fit <- stochtree::bcf(cbind(xs,ws), zs, ys, pi_train = zs,
#                         X_test = test.sample, Z_test = zs[test], pi_test = zs[test],
#                         num_gfr = 0, num_burnin = burnin, num_mcmc = num_sweeps)
#   fit$tau_hat_test[,(burnin+1):num_sweeps]
# }
# 
# cl <- makeCluster(no_cores,type="SOCK")
# registerDoParallel(cl)
# clusterExport(cl,varlist=ls())
# time.tbart <- system.time({
#   bcf <- parLapply(cl,1:s,fit.bcf)
# })
# stopCluster(cl)
## Calculate RMSE--------------------------------------------------------------
barddt <- sapply(barddt,rowMeans)
tbart <- sapply(tbart,rowMeans)
# bcf <- sapply(bcf,rowMeans)
rmse.barddt <- mapply(function(fit,tau) sqrt(mean((fit-tau)^2)),barddt,cate)
rmse.tbart <- mapply(function(fit,tau) sqrt(mean((fit-tau)^2)),tbart,cate)
# rmse.bcf <- mapply(function(fit,tau) sqrt(mean((fit-tau)^2)),bcf,cate)
rmse.ate <- sapply(cate, function(tau) sqrt(mean((tau-ate)^2)))
print(c("RMSE for BARDDT: ",mean(rmse.barddt/rmse.ate)))
print(c("RMSE for T-BART: ",mean(rmse.tbart/rmse.ate)))
# print(c("RMSE for BCF: ",mean(rmse.bcf/rmse.ate)))
# plot(rmse.barddt/rmse.ate,rmse.tbart/rmse.ate,bty="n",pch=19,col="darkblue")
# abline(a=0,b=1,col="red")
##
sample <- 1
# 
# plot(cate[[sample]],barddt[[sample]],bty="n",pch=19,col="darkblue")
# abline(a=0,b=1,col="red")
# plot(cate[[sample]],tbart[[sample]],bty="n",pch=19,col="darkblue")
# abline(a=0,b=1,col="red")
matplot(sapply(w,rowMeans)[test[,sample],sample],cbind(cate[[sample]],tbart[[sample]],barddt[[sample]]),pch=20,col=c("black","maroon","orange"))