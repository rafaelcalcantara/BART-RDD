setwd("~/Git/BART-RDD")
files <- length(list.files("Data"))
params <- data.frame(n=rep(0,files),rho=rep(0,files))
# samples <- 1:1000
##
bart.rdd.cate <- vector("list",files)
bart.rdd.cate.rmse <- rep(0,files)
bart.rdd.cate.bias <- rep(0,files)
bart.rdd.cate.var <- rep(0,files)
##
sbart.cate <- vector("list",files)
sbart.cate.rmse <- rep(0,files)
sbart.cate.bias <- rep(0,files)
sbart.cate.var <- rep(0,files)
##
tbart.cate <- vector("list",files)
tbart.cate.rmse <- rep(0,files)
tbart.cate.bias <- rep(0,files)
tbart.cate.var <- rep(0,files)
##
oracle.cate <- vector("list",files)
oracle.cate.rmse <- rep(0,files)
oracle.cate.bias <- rep(0,files)
oracle.cate.var <- rep(0,files)
##
# polynomial.1.cate <- vector("list",files)
# polynomial.1.cate.rmse <- rep(0,files)
# polynomial.1.cate.bias <- rep(0,files)
# polynomial.1.cate.var <- rep(0,files)
# polynomial.1.cate.dt <- rep(0,files)
# ##
# polynomial.2.cate <- vector("list",files)
# polynomial.2.cate.rmse <- rep(0,files)
# polynomial.2.cate.bias <- rep(0,files)
# polynomial.2.cate.var <- rep(0,files)
# polynomial.2.cate.dt <- rep(0,files)
# ##
# polynomial.3.cate <- vector("list",files)
# polynomial.3.cate.rmse <- rep(0,files)
# polynomial.3.cate.bias <- rep(0,files)
# polynomial.3.cate.var <- rep(0,files)
# polynomial.3.cate.dt <- rep(0,files)
##
polynomial.full.cate <- vector("list",files)
polynomial.full.cate.rmse <- rep(0,files)
polynomial.full.cate.bias <- rep(0,files)
polynomial.full.cate.var <- rep(0,files)
##
ate.cate <- rep(0,files)
##
for (i in 1:files)
{
  print(paste0("DGP: ",i))
  data <- readRDS(paste0("Data/dgp_",i,".rds"))
  n <- data$n
  Owidth <- data$h
  test <- sapply(samples,function(i) -Owidth[i] <= data$x[,i] & data$x[,i] <= Owidth[i])
  cate <- sapply(samples, function(i) data$tau.x[test[,i],i])
  params[i,] <- c(n,data$rho)
  ate.cate[i] <- mean(sapply(cate,function(i) sqrt(mean((i-data$tau)^2))))
  bart.rdd <- readRDS(paste0("Results/bart_rdd_",i,".rds"))
  sbart <- readRDS(paste0("Results/sbart_",i,".rds"))
  tbart <- readRDS(paste0("Results/tbart_",i,".rds"))
  oracle <- readRDS(paste0("Results/oracle_",i,".rds"))
  # polynomial.1 <- readRDS(paste0("Results/polynomial_1_",i,".rds"))
  # polynomial.2 <- readRDS(paste0("Results/polynomial_2_",i,".rds"))
  # polynomial.3 <- readRDS(paste0("Results/polynomial_3_",i,".rds"))
  polynomial.full <- readRDS(paste0("Results/polynomial_full_",i,".rds"))
  ## BART- RDD
  bart.rdd.cate[[i]] <- sapply(bart.rdd$results[samples],rowMeans)
  bart.rdd.cate.rmse[i] <- mean(sqrt(mapply(function(a,b) mean((a-b)^2),bart.rdd.cate[[i]],cate)))
  bart.rdd.cate.bias[i] <- mean(unlist(mapply(function(a,b) a-b,bart.rdd.cate[[i]],cate)))
  bart.rdd.cate.var[i] <- mean(sapply(bart.rdd.cate[[i]],var))
  ## S-BART
  sbart.cate[[i]] <- sapply(sbart$results[samples],rowMeans)
  sbart.cate.rmse[i] <- mean(sqrt(mapply(function(a,b) mean((a-b)^2),sbart.cate[[i]],cate)))
  sbart.cate.bias[i] <- mean(unlist(mapply(function(a,b) a-b,sbart.cate[[i]],cate)))
  sbart.cate.var[i] <- mean(sapply(sbart.cate[[i]],var))
  ## T-BART
  tbart.cate[[i]] <- sapply(tbart$results[samples],rowMeans)
  tbart.cate.rmse[i] <- mean(sqrt(mapply(function(a,b) mean((a-b)^2),tbart.cate[[i]],cate)))
  tbart.cate.bias[i] <- mean(unlist(mapply(function(a,b) a-b,tbart.cate[[i]],cate)))
  tbart.cate.var[i] <- mean(sapply(tbart.cate[[i]],var))
  ## Oracle
  oracle.cate[[i]] <- oracle$results[samples]
  oracle.cate.rmse[i] <- mean(sqrt(mapply(function(a,b) mean((a-b)^2),oracle.cate[[i]],cate)))
  oracle.cate.bias[i] <- mean(unlist(mapply(function(a,b) a-b,oracle.cate[[i]],cate)))
  oracle.cate.var[i] <- mean(sapply(oracle.cate[[i]],var))
  # ## Polynomial 1
  # polynomial.1.cate[[i]] <- polynomial.1$results[samples]
  # polynomial.1.cate.rmse[i] <- mean(sqrt(mapply(function(i,j) mean((i-j)^2),polynomial.1.cate[[i]],cate)))
  # polynomial.1.cate.bias[i] <- mean(unlist(mapply(function(i,j) i-j,polynomial.1.cate[[i]],cate)))
  # polynomial.1.cate.var[i] <- mean(sapply(polynomial.1.cate[[i]],var))
  # polynomial.1.cate.dt[i] <- mean(sapply(polynomial.1.cate[[i]],function(tau.hat) sqrt(mean((tau.hat-data$tau)^2))))
  # ## Polynomial 2
  # polynomial.2.cate[[i]] <- polynomial.2$results[samples]
  # polynomial.2.cate.rmse[i] <- mean(sqrt(mapply(function(i,j) mean((i-j)^2),polynomial.2.cate[[i]],cate)))
  # polynomial.2.cate.bias[i] <- mean(unlist(mapply(function(i,j) i-j,polynomial.2.cate[[i]],cate)))
  # polynomial.2.cate.var[i] <- mean(sapply(polynomial.2.cate[[i]],var))
  # polynomial.2.cate.dt[i] <- mean(sapply(polynomial.2.cate[[i]],function(tau.hat) sqrt(mean((tau.hat-data$tau)^2))))
  # ## Polynomial 3
  # polynomial.3.cate[[i]] <- polynomial.3$results[samples]
  # polynomial.3.cate.rmse[i] <- mean(sqrt(mapply(function(i,j) mean((i-j)^2),polynomial.3.cate[[i]],cate)))
  # polynomial.3.cate.bias[i] <- mean(unlist(mapply(function(i,j) i-j,polynomial.3.cate[[i]],cate)))
  # polynomial.3.cate.var[i] <- mean(sapply(polynomial.3.cate[[i]],var))
  # polynomial.3.cate.dt[i] <- mean(sapply(polynomial.3.cate[[i]],function(tau.hat) sqrt(mean((tau.hat-data$tau)^2))))
  ## Polynomial full
  polynomial.full.cate[[i]] <- polynomial.full$results[samples]
  polynomial.full.cate.rmse[i] <- mean(sqrt(mapply(function(i,j) mean((i-j)^2),polynomial.full.cate[[i]],cate)))
  polynomial.full.cate.bias[i] <- mean(unlist(mapply(function(i,j) i-j,polynomial.full.cate[[i]],cate)))
  polynomial.full.cate.var[i] <- mean(sapply(polynomial.full.cate[[i]],var))
}
# Aggregating results
## RMSE
rmse.cate <- cbind(`Oracle`=oracle.cate.rmse,
                   `BART-RDD`=bart.rdd.cate.rmse,
                   `Poly_full`=polynomial.full.cate.rmse,
                   `T-BART`=tbart.cate.rmse,
                   `S-BART`=sbart.cate.rmse
                   # `Poly_1`=polynomial.1.cate.rmse,
                   # `Poly_2`=polynomial.2.cate.rmse,
                   # `Poly_3`=polynomial.3.cate.rmse
                   )
# ## Bias
bias.cate <- cbind(`Oracle`=oracle.cate.bias,
                   `BART-RDD`=bart.rdd.cate.bias,
                   `Poly_full`=polynomial.full.cate.bias,
                   `T-BART`=tbart.cate.bias,
                   `S-BART`=sbart.cate.bias
                   # `Poly_1`=polynomial.1.cate.bias,
                   # `Poly_2`=polynomial.2.cate.bias,
                   # `Poly_3`=polynomial.3.cate.bias
                   )
## Variance
var.cate <- cbind(`Oracle`=oracle.cate.var,
                   `BART-RDD`=bart.rdd.cate.var,
                  `Poly_full`=polynomial.full.cate.var,
                   `T-BART`=tbart.cate.var,
                   `S-BART`=sbart.cate.var
                   # `Poly_1`=polynomial.1.cate.var,
                   # `Poly_2`=polynomial.2.cate.var,
                   # `Poly_3`=polynomial.3.cate.var,
                   )
# ## Delta_tau
# dt.cate <- cbind(`Oracle`=oracle.cate.dt,
#                    `BART-RDD`=bart.rdd.cate.dt,
#                    `T-BART`=tbart.cate.dt,
#                    `S-BART`=sbart.cate.dt,
#                    `Poly_1`=polynomial.1.cate.dt,
#                    `Poly_2`=polynomial.2.cate.dt,
#                    `Poly_3`=polynomial.3.cate.dt,
#                    `Poly_full`=polynomial.full.cate.dt)
save.image("Results/sims.RData")