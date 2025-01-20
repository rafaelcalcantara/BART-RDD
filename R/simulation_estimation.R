library(doParallel)
no_cores <- 10
## Common parameters
p_categorical <- 0
burnin <- 50
num_sweeps <- 150
# BARDDT parameters
Omin <- 1
Opct <- 0.25
ntrees_con <- 30
ntrees_mod <- 30
# T-BART parameters
ntrees <- 30
## Set Owidth and test set
h <- apply(x,2,function(i) h.grid(i,c,pts_in_window))
test <- sapply(1:s, function(i) c-h[i] <= x[,i] & x[,i] <= c+h[i])
cate <- sapply(1:s, function(i) cate[test[,i],i])
## BARDDT fit------------------------------------------------------------------
fit.barddt <- function(i)
{
  print(paste0("Sample: ",i))
  ys <- y[,i]
  ws <- as.matrix(w[[i]])
  xs <- x[,i]
  Owidth <- h[i]
  sample <- 5*Owidth
  train <- c-sample < xs & xs < c+sample
  fit <- XBART::XBCF.rd(ys[train], ws[train,], xs[train], c,
                        Owidth = Owidth, Omin = Omin, Opct = Opct,
                        num_cutpoints = n,
                        num_trees_con = ntrees_con, num_trees_mod = ntrees_mod,
                        num_sweeps = num_sweeps,
                        burnin = burnin,
                        Nmin = 1,
                        p_categorical_con = p_categorical,
                        p_categorical_mod = p_categorical,
                        tau_con = var(ys[train])/ntrees_con, tau_mod = var(ys[train])/ntrees_mod,
                        alpha_mod = 0.95, beta_mod = 1.25, update_tau = TRUE,
                        parallel=F)
  test <- -Owidth+c<=xs & xs<=Owidth+c
  pred <- XBART::predict.XBCFrd(fit,ws[test,],rep(c,sum(test)))
  list(pred=pred$tau.adj[,(burnin+1):num_sweeps],
       count1=fit$count_fail_1[(burnin+1):num_sweeps],
       count2=fit$count_fail_2[(burnin+1):num_sweeps],
       cutoff_nodes_con=fit$cutoff_nodes_con[(burnin+1):num_sweeps,],
       invalid_nodes_1_con=fit$invalid_nodes_1_con[(burnin+1):num_sweeps,],
       invalid_nodes_2_con=fit$invalid_nodes_2_con[(burnin+1):num_sweeps,],
       cutoff_nodes_mod=fit$cutoff_nodes_mod[(burnin+1):num_sweeps,],
       invalid_nodes_1_mod=fit$invalid_nodes_1_mod[(burnin+1):num_sweeps,],
       invalid_nodes_2_mod=fit$invalid_nodes_2_mod[(burnin+1):num_sweeps,])
}
cl <- makeCluster(no_cores,type="SOCK",outfile="")
registerDoParallel(cl)
clusterExport(cl,varlist=ls())
time.barddt <- system.time({
  barddt <- parLapply(cl,1:s,fit.barddt)
})
stopCluster(cl)
### Save count of splits rejected per constraint
constraint.fails.1 <- lapply(barddt, function(i) i$count1)
constraint.fails.2 <- lapply(barddt, function(i) i$count2)
cutoff.nodes.con <- lapply(barddt, function(i) i$cutoff_nodes_con)
invalid.nodes.1.con <- lapply(barddt, function(i) i$invalid_nodes_1_con)
invalid.nodes.2.con <- lapply(barddt, function(i) i$invalid_nodes_2_con)
cutoff.nodes.mod <- lapply(barddt, function(i) i$cutoff_nodes_mod)
invalid.nodes.1.mod <- lapply(barddt, function(i) i$invalid_nodes_1_mod)
invalid.nodes.2.mod <- lapply(barddt, function(i) i$invalid_nodes_2_mod)
## T-BART fit------------------------------------------------------------------
fit.tbart <- function(i)
{
  print(paste0("Sample: ",i))
  ys <- y[,i]
  ws <- as.matrix(w[[i]])
  xs <- x[,i]
  zs <- z[,i]
  Owidth <- h[i]
  sample <- 5*Owidth
  train <- c-sample < xs & xs < c+sample
  ys <- ys[train]
  xs <- xs[train]
  ws <- ws[train,]
  zs <- zs[train]
  fit1 <- XBART::XBART(ys[zs==1], cbind(xs,ws)[zs==1,], num_trees = ntrees,
                       num_cutpoints = sum(zs==1), num_sweeps = num_sweeps,
                       burnin = burnin, p_categorical = p_categorical,
                       tau = var(ys[zs==1])/ntrees, parallel=F)
  fit0 <- XBART::XBART(ys[zs==0], cbind(xs,ws)[zs==0,], num_trees = ntrees,
                       num_cutpoints = sum(zs==0), num_sweeps = num_sweeps,
                       burnin = burnin, p_categorical = p_categorical,
                       tau = var(ys[zs==0])/ntrees, parallel=F)
  test <- -Owidth+c<=xs & xs<=Owidth+c
  test.sample <- cbind(c,ws)[test,]
  pred1 <- XBART::predict.XBART(fit1,test.sample)[,(burnin+1):num_sweeps]
  pred0 <- XBART::predict.XBART(fit0,test.sample)[,(burnin+1):num_sweeps]
  pred1-pred0
}
cl <- makeCluster(no_cores,type="SOCK",outfile="")
registerDoParallel(cl)
clusterExport(cl,varlist=ls())
time.tbart <- system.time({
  tbart <- parLapply(cl,1:s,fit.tbart)
})
stopCluster(cl)
## Calculate RMSE--------------------------------------------------------------
barddt <- sapply(barddt,function(i) rowMeans(i$pred))
tbart <- sapply(tbart,rowMeans)
rmse.barddt <- mapply(function(fit,tau) sqrt(mean((fit-tau)^2)),barddt,cate)
rmse.tbart <- mapply(function(fit,tau) sqrt(mean((fit-tau)^2)),tbart,cate)
rmse.ate <- sapply(cate, function(tau) sqrt(mean((tau-ate)^2)))
print(c("RMSE for BARDDT: ",mean(rmse.barddt/rmse.ate)))
print(c("RMSE for T-BART: ",mean(rmse.tbart/rmse.ate)))
plot(rmse.barddt/rmse.ate,rmse.tbart/rmse.ate,bty="n",pch=19,col="darkblue")
abline(a=0,b=1,col="red")
##
plot(cate[[sample]],barddt[[sample]],bty="n",pch=19,col="darkblue")
abline(a=0,b=1,col="red")
plot(cate[[sample]],tbart[[sample]],bty="n",pch=19,col="darkblue")
abline(a=0,b=1,col="red")
matplot(sapply(w,rowMeans)[test[,sample],sample],cbind(cate[[sample]],tbart[[sample]],barddt[[sample]]),pch=20,col=c("black","maroon","orange"))