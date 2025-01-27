source("simulation_estimator_functions.R")
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