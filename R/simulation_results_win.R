setwd("~/../Git/BART-RDD/")
files <- length(list.files("Data"))
## RMSE
rmse.bart <- rep(0,files)
sb.bart <- rep(0,files)
v.bart <- rep(0,files)
ci.bart <- rep(0,files)
cov.bart <- rep(0,files)
zero.bart <- rep(0,files)
rmse.llr <- matrix(0,files,3)
sb.llr <- matrix(0,files,3)
v.llr <- matrix(0,files,3)
ci.llr <- matrix(0,files,3)
cov.llr <- matrix(0,files,3)
zero.llr <- matrix(0,files,3)
for (i in 1:files)
{
  data <- readRDS(paste0("Data/dgp_",i,".rds"))
  ate <- data$ate
  bart <- readRDS(paste0("Results/bart_rdd_",i,".rds"))
  llr <- readRDS(paste0("Results/llr0_",i,".rds"))
  rmse.bart[i] <- sqrt(mean((sapply(bart$results,mean)-ate)^2))
  sb.bart[i] <- mean(sapply(bart$results,mean)-ate)^2
  v.bart[i] <- var(sapply(bart$results,mean))
  ci.bart[i] <- mean(sapply(bart$results,function(i) quantile(colMeans(i),0.975)-quantile(colMeans(i),0.025)))
  cov.bart[i] <- mean(sapply(bart$results,function(i) quantile(colMeans(i),0.975))-ate>0 & sapply(bart$results,function(i) quantile(colMeans(i),0.025))-ate<0)
  zero.bart[i] <- mean(sapply(bart$results,function(i) quantile(colMeans(i),0.975))>0 & sapply(bart$results,function(i) quantile(colMeans(i),0.025))<0)
  rmse.llr[i,] <- rowMeans((sapply(llr$results,function(i) i$coef)-ate)^2)
  sb.llr[i,] <- rowMeans(sapply(llr$results,function(i) i$coef)-ate)^2
  v.llr[i,] <- apply(sapply(llr$results,function(i) i$coef),1,var)
  ci.llr[i,] <- rowMeans(sapply(llr$results,function(i) i$ci[,2]-i$ci[,1]))
  cov.llr[i,] <- rowMeans(sapply(llr$results,function(i) i$ci[,2])-ate>0 & sapply(llr$results,function(i) i$ci[,1])-ate<0)
  zero.llr[i,] <- rowMeans(sapply(llr$results,function(i) i$ci[,2])>0 & sapply(llr$results,function(i) i$ci[,1])<0)
}