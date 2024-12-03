setwd("../")
files <- length(list.files("Data"))
Ow <- c(0.021,0.066)
params <- data.frame(delta.mu=rep(0,files),delta.tau=rep(0,files),
                     level=rep(0,files),n=rep(0,files),sig_error=rep(0,files))
##
bart.rdd.cate <- vector("list",files)
bart.rdd.cate.rmse <- rep(0,files)
bart.rdd.cate.bias <- rep(0,files)
bart.rdd.cate.var <- rep(0,files)
# bart.rdd.cate.cov <- rep(0,files)
# bart.rdd.cate.ci <- rep(0,files)
##
sbart.cate <- vector("list",files)
sbart.cate.rmse <- rep(0,files)
sbart.cate.bias <- rep(0,files)
sbart.cate.var <- rep(0,files)
# sbart.cate.cov <- rep(0,files)
# sbart.cate.ci <- rep(0,files)
##
tbart.cate <- vector("list",files)
tbart.cate.rmse <- rep(0,files)
tbart.cate.bias <- rep(0,files)
tbart.cate.var <- rep(0,files)
# tbart.cate.cov <- rep(0,files)
# tbart.cate.ci <- rep(0,files)
##
oracle.cate <- vector("list",files)
oracle.cate.rmse <- rep(0,files)
oracle.cate.bias <- rep(0,files)
oracle.cate.var <- rep(0,files)
# oracle.cate.cov <- rep(0,files)
# oracle.cate.ci <- rep(0,files)
##
polynomial.cate <- vector("list",files)
polynomial.cate.rmse <- rep(0,files)
polynomial.cate.bias <- rep(0,files)
polynomial.cate.var <- rep(0,files)
polynomial.cate.cov <- rep(0,files)
polynomial.cate.ci <- rep(0,files)
##
for (i in 1:files)
{
  print(paste0("DGP: ",i))
  data <- readRDS(paste0("Data/dgp_",i,".rds"))
  n <- data$n
  Owidth <- ifelse(n==500,Ow[2],Ow[1])
  test <- apply(data$x,2,function(i) i >= -Owidth & Owidth >= i)
  cate <- sapply(1:ncol(data$tau.x), function(i) data$tau.x[test[,i],i])
  params[i,] <- c(data$delta_mu,data$delta_tau,data$level,n,data$sig_error)
  bart.rdd <- readRDS(paste0("Results/bart_rdd_",i,".rds"))
  sbart <- readRDS(paste0("Results/sbart_",i,".rds"))
  tbart <- readRDS(paste0("Results/tbart_",i,".rds"))
  oracle <- readRDS(paste0("Results/oracle_",i,".rds"))
  polynomial <- readRDS(paste0("Results/polynomial_",i,".rds"))
  ## BART- RDD
  bart.rdd.cate[[i]] <- sapply(bart.rdd$results,rowMeans)
  bart.rdd.cate.rmse[i] <- mean(sqrt(mapply(function(a,b) mean((a-b)^2),bart.rdd.cate[[i]],cate)))
  bart.rdd.cate.bias[i] <- mean(unlist(mapply(function(a,b) a-b,bart.rdd.cate[[i]],cate)))
  bart.rdd.cate.var[i] <- mean(sapply(bart.rdd.cate[[i]],var))
  ## S-BART
  sbart.cate[[i]] <- sapply(sbart$results,rowMeans)
  sbart.cate.rmse[i] <- mean(sqrt(mapply(function(a,b) mean((a-b)^2),sbart.cate[[i]],cate)))
  sbart.cate.bias[i] <- mean(unlist(mapply(function(a,b) a-b,sbart.cate[[i]],cate)))
  sbart.cate.var[i] <- mean(sapply(sbart.cate[[i]],var))
  ## T-BART
  tbart.cate[[i]] <- sapply(tbart$results,rowMeans)
  tbart.cate.rmse[i] <- mean(sqrt(mapply(function(a,b) mean((a-b)^2),tbart.cate[[i]],cate)))
  tbart.cate.bias[i] <- mean(unlist(mapply(function(a,b) a-b,tbart.cate[[i]],cate)))
  tbart.cate.var[i] <- mean(sapply(tbart.cate[[i]],var))
  ## Oracle
  oracle.cate[[i]] <- oracle$results
  oracle.cate.rmse[i] <- mean(sqrt(mapply(function(a,b) mean((a-b)^2),oracle.cate[[i]],cate)))
  oracle.cate.bias[i] <- mean(unlist(mapply(function(a,b) a-b,oracle.cate[[i]],cate)))
  oracle.cate.var[i] <- mean(sapply(oracle.cate[[i]],var))
  ## Polynomial
  polynomial.cate[[i]] <- polynomial$results
  polynomial.cate.rmse[i] <- mean(sqrt(mapply(function(i,j) mean((i-j)^2),polynomial.cate[[i]],cate)))
  polynomial.cate.bias[i] <- mean(unlist(mapply(function(i,j) i-j,polynomial.cate[[i]],cate)))
  polynomial.cate.var[i] <- mean(sapply(polynomial.cate[[i]],var))
}
# Aggregating results
## RMSE
rmse.cate <- cbind(`BART-RDD`=bart.rdd.cate.rmse,
                   `S-BART`=sbart.cate.rmse,
                   `T-BART`=tbart.cate.rmse,
                   `Oracle`=oracle.cate.rmse,
                   `Polynomial`=polynomial.cate.rmse)
## Bias
bias.cate <- cbind(`BART-RDD`=bart.rdd.cate.bias,
                   `S-BART`=sbart.cate.bias,
                   `T-BART`=tbart.cate.bias,
                   `Oracle`=oracle.cate.bias,
                   `Polynomial`=polynomial.cate.bias)
## Variance
var.cate <- cbind(`BART-RDD`=bart.rdd.cate.var,
                   `S-BART`=sbart.cate.var,
                   `T-BART`=tbart.cate.var,
                   `Oracle`=oracle.cate.var,
                   `Polynomial`=polynomial.cate.var)
save.image("Results/sims.RData")