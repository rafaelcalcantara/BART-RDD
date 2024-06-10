setwd("~/../Git/BART-RDD/")
files <- length(list.files("Data"))
Owidth <- 0.1
##
params <- data.frame(tau=rep(0,files),delta.mu=rep(0,files),delta.tau=rep(0,files),kappa=rep(0,files))
##
bart.rdd.ate <- vector("list",files)
bcf.ate <- vector("list",files)
sbart.ate <- vector("list",files)
tbart.ate <- vector("list",files)
cgs.ate <- vector("list",files)
llr.ate <- vector("list",files)
##
rmse.bart <- rep(0,files)
sb.bart <- rep(0,files)
v.bart <- rep(0,files)
ci.bart <- rep(0,files)
cov.bart <- rep(0,files)
zero.bart <- rep(0,files)
bart.rdd.cate.bias <- vector("list",files)
bart.rdd.cate.var <- vector("list",files)
bart.rdd.cate <- rep(0,files)
bart.rdd.cate.cov <- rep(0,files)
bart.rdd.cate.ci <- rep(0,files)
##
rmse.bcf <- rep(0,files)
sb.bcf <- rep(0,files)
v.bcf <- rep(0,files)
ci.bcf <- rep(0,files)
cov.bcf <- rep(0,files)
zero.bcf <- rep(0,files)
bcf.cate.bias <- vector("list",files)
bcf.cate.var <- vector("list",files)
bcf.cate <- rep(0,files)
bcf.cate.cov <- rep(0,files)
bcf.cate.ci <- rep(0,files)
##
rmse.sbart <- rep(0,files)
sb.sbart <- rep(0,files)
v.sbart <- rep(0,files)
ci.sbart <- rep(0,files)
cov.sbart <- rep(0,files)
zero.sbart <- rep(0,files)
sbart.cate.bias <- vector("list",files)
sbart.cate.var <- vector("list",files)
sbart.cate <- rep(0,files)
sbart.cate.cov <- rep(0,files)
sbart.cate.ci <- rep(0,files)
##
rmse.tbart <- rep(0,files)
sb.tbart <- rep(0,files)
v.tbart <- rep(0,files)
ci.tbart <- rep(0,files)
cov.tbart <- rep(0,files)
zero.tbart <- rep(0,files)
tbart.cate.bias <- vector("list",files)
tbart.cate.var <- vector("list",files)
tbart.cate <- rep(0,files)
tbart.cate.cov <- rep(0,files)
tbart.cate.ci <- rep(0,files)
##
# rmse.llr0 <- matrix(0,files,3)
# sb.llr0 <- matrix(0,files,3)
# v.llr0 <- matrix(0,files,3)
# ci.llr0 <- matrix(0,files,3)
# cov.llr0 <- matrix(0,files,3)
# zero.llr0 <- matrix(0,files,3)
##
rmse.llr1 <- matrix(0,files,3)
sb.llr1 <- matrix(0,files,3)
v.llr1 <- matrix(0,files,3)
ci.llr1 <- matrix(0,files,3)
cov.llr1 <- matrix(0,files,3)
zero.llr1 <- matrix(0,files,3)
##
# rmse.llr2 <- matrix(0,files,3)
# sb.llr2 <- matrix(0,files,3)
# v.llr2 <- matrix(0,files,3)
# ci.llr2 <- matrix(0,files,3)
# cov.llr2 <- matrix(0,files,3)
# zero.llr2 <- matrix(0,files,3)
##
# rmse.cgs0 <- rep(0,files)
# sb.cgs0 <- rep(0,files)
# v.cgs0 <- rep(0,files)
# ci.cgs0 <- rep(0,files)
# cov.cgs0 <- rep(0,files)
# zero.cgs0 <- rep(0,files)
##
rmse.cgs1 <- rep(0,files)
sb.cgs1 <- rep(0,files)
v.cgs1 <- rep(0,files)
ci.cgs1 <- rep(0,files)
cov.cgs1 <- rep(0,files)
zero.cgs1 <- rep(0,files)
##
for (i in 1:files)
{
  print(paste0("DGP: ",i))
  data <- readRDS(paste0("Data/dgp_",i,".rds"))
  ## test <- readRDS(paste0("Data/test_dgp_",i,".rds"))
  ## cate <- test$cate
  test <- apply(data$x,2,function(i) i >= -Owidth & Owidth >= i)
  ate <- data$tau
  cate <- sapply(1:ncol(data$tau.x), function(i) data$tau.x[test[,i],i])
  ate <- sapply(cate,mean)
  params[i,] <- c(data$tau,data$delta_mu,data$delta_tau,data$kappa)
  bart <- readRDS(paste0("Results/bart_rdd_",i,".rds"))
  bcf <- readRDS(paste0("Results/bcf_",i,".rds"))
  sbart <- readRDS(paste0("Results/sbart_",i,".rds"))
  tbart <- readRDS(paste0("Results/tbart_",i,".rds"))
  # llr0 <- readRDS(paste0("Results/llr0_",i,".rds"))
  llr1 <- readRDS(paste0("Results/llr1_",i,".rds"))
  ## llr2 <- readRDS(paste0("Results/llr2_",i,".rds"))
  # cgs0 <- readRDS(paste0("Results/cgs0_",i,".rds"))
  cgs1 <- readRDS(paste0("Results/cgs1_",i,".rds"))
  ##
  bart.rdd.ate[[i]] <- sapply(bart$results,colMeans)
  rmse.bart[i] <- sqrt(mean((colMeans(bart.rdd.ate[[i]])-ate)^2))
  sb.bart[i] <- mean(colMeans(bart.rdd.ate[[i]])-ate)
  v.bart[i] <- var(colMeans(bart.rdd.ate[[i]]))
  ci.bart[i] <- mean(sapply(bart$results,function(i) quantile(colMeans(i),0.975)-quantile(colMeans(i),0.025)))
  cov.bart[i] <- mean(sapply(bart$results,function(i) quantile(colMeans(i),0.975))-ate>0 & sapply(bart$results,function(i) quantile(colMeans(i),0.025))-ate<0)
  zero.bart[i] <- mean(sapply(bart$results,function(i) quantile(colMeans(i),0.975))>0 & sapply(bart$results,function(i) quantile(colMeans(i),0.025))<0)
  ## CATE
  bart.cate <- sapply(bart$results,rowMeans)
  bart.rdd.cate.bias[[i]] <- mapply(function(i,j) i-j,bart.cate,cate)
  bart.rdd.cate.var[[i]] <- sapply(bart.cate,var)
  bart.rdd.cate[i] <- mean(sqrt(mapply(function(i,j) mean((i-j)^2),bart.cate,cate)))
  bart.rdd.cate.int <- lapply(bart$results,function(i) t(apply(i,1,quantile,c(0.025,0.975))))
  bart.rdd.cate.cov[i] <- mean(sapply(1:ncol(data$x), function(i) mean(bart.rdd.cate.int[[i]][,1]<cate[[i]] & cate[[i]] < bart.rdd.cate.int[[i]][,2])))
  bart.rdd.cate.ci[i] <- mean(sapply(1:ncol(data$x), function(i) mean(-bart.rdd.cate.int[[i]][,1] + bart.rdd.cate.int[[i]][,2])))
  ##
  bcf.ate[[i]] <- sapply(bcf$results,colMeans)
  rmse.bcf[i] <- sqrt(mean((colMeans(bcf.ate[[i]])-ate)^2))
  sb.bcf[i] <- mean(colMeans(bcf.ate[[i]])-ate)
  v.bcf[i] <- var(colMeans(bcf.ate[[i]]))
  ci.bcf[i] <- mean(sapply(bcf$results,function(i) quantile(colMeans(i),0.975)-quantile(colMeans(i),0.025)))
  cov.bcf[i] <- mean(sapply(bcf$results,function(i) quantile(colMeans(i),0.975))-ate>0 & sapply(bcf$results,function(i) quantile(colMeans(i),0.025))-ate<0)
  zero.bcf[i] <- mean(sapply(bcf$results,function(i) quantile(colMeans(i),0.975))>0 & sapply(bcf$results,function(i) quantile(colMeans(i),0.025))<0)
  ## CATE
  bcf.cate.temp <- sapply(bcf$results,rowMeans)
  bcf.cate.bias[[i]] <- mapply(function(i,j) i-j,bcf.cate.temp,cate)
  bcf.cate.var[[i]] <- sapply(bcf.cate.temp,var)
  bcf.cate[i] <- mean(sqrt(mapply(function(i,j) mean((i-j)^2),bcf.cate.temp,cate)))
  bcf.cate.int <- lapply(bcf$results,function(i) t(apply(i,1,quantile,c(0.025,0.975))))
  bcf.cate.cov[i] <- mean(sapply(1:ncol(data$x), function(i) mean(bcf.cate.int[[i]][,1]<cate[[i]] & cate[[i]]<bcf.cate.int[[i]][,2])))
  bcf.cate.ci[i] <- mean(sapply(1:ncol(data$x), function(i) mean(-bcf.cate.int[[i]][,1] + bcf.cate.int[[i]][,2])))
  ##
  sbart.ate[[i]] <- sapply(sbart$results,colMeans)
  rmse.sbart[i] <- sqrt(mean((colMeans(sbart.ate[[i]])-ate)^2))
  sb.sbart[i] <- mean(colMeans(sbart.ate[[i]])-ate)
  v.sbart[i] <- var(colMeans(sbart.ate[[i]]))
  ci.sbart[i] <- mean(sapply(sbart$results,function(i) quantile(colMeans(i),0.975)-quantile(colMeans(i),0.025)))
  cov.sbart[i] <- mean(sapply(sbart$results,function(i) quantile(colMeans(i),0.975))-ate>0 & sapply(sbart$results,function(i) quantile(colMeans(i),0.025))-ate<0)
  zero.sbart[i] <- mean(sapply(sbart$results,function(i) quantile(colMeans(i),0.975))>0 & sapply(sbart$results,function(i) quantile(colMeans(i),0.025))<0)
  ## CATE
  sbart.cate.temp <- sapply(sbart$results,rowMeans)
  sbart.cate.bias[[i]] <- mapply(function(i,j) i-j,sbart.cate.temp,cate)
  sbart.cate.var[[i]] <- sapply(sbart.cate.temp,var)
  sbart.cate[i] <- mean(sqrt(mapply(function(i,j) mean((i-j)^2),sbart.cate.temp,cate)))
  sbart.cate.int <- lapply(sbart$results,function(i) t(apply(i,1,quantile,c(0.025,0.975))))
  sbart.cate.cov[i] <- mean(sapply(1:ncol(data$x), function(i) mean(sbart.cate.int[[i]][,1]<cate[[i]] & cate[[i]]<sbart.cate.int[[i]][,2])))
  sbart.cate.ci[i] <- mean(sapply(1:ncol(data$x), function(i) mean(-sbart.cate.int[[i]][,1] + sbart.cate.int[[i]][,2])))
  ##
  tbart.ate[[i]] <- sapply(tbart$results,colMeans)
  rmse.tbart[i] <- sqrt(mean((colMeans(tbart.ate[[i]])-ate)^2))
  sb.tbart[i] <- mean(colMeans(tbart.ate[[i]])-ate)
  v.tbart[i] <- var(colMeans(tbart.ate[[i]]))
  ci.tbart[i] <- mean(sapply(tbart$results,function(i) quantile(colMeans(i),0.975)-quantile(colMeans(i),0.025)))
  cov.tbart[i] <- mean(sapply(tbart$results,function(i) quantile(colMeans(i),0.975))-ate>0 & sapply(tbart$results,function(i) quantile(colMeans(i),0.025))-ate<0)
  zero.tbart[i] <- mean(sapply(tbart$results,function(i) quantile(colMeans(i),0.975))>0 & sapply(tbart$results,function(i) quantile(colMeans(i),0.025))<0)
  ## CATE
  tbart.cate.temp <- sapply(tbart$results,rowMeans)
  tbart.cate.bias[[i]] <- mapply(function(i,j) i-j,tbart.cate.temp,cate)
  tbart.cate.var[[i]] <- sapply(tbart.cate.temp,var)
  tbart.cate[i] <- mean(sqrt(mapply(function(i,j) mean((i-j)^2),tbart.cate.temp,cate)))
  tbart.cate.int <- lapply(tbart$results,function(i) t(apply(i,1,quantile,c(0.025,0.975))))
  tbart.cate.cov[i] <- mean(sapply(1:ncol(data$x), function(i) mean(tbart.cate.int[[i]][,1]<cate[[i]] & cate[[i]]<tbart.cate.int[[i]][,2])))
  tbart.cate.ci[i] <- mean(sapply(1:ncol(data$x), function(i) mean(-tbart.cate.int[[i]][,1] + tbart.cate.int[[i]][,2])))
  ##
  # rmse.llr0[i,] <- sqrt(rowMeans((sapply(llr0$results,function(i) i$coef)-ate)^2))
  # sb.llr0[i,] <- rowMeans(sapply(llr0$results,function(i) i$coef)-ate)^2
  # v.llr0[i,] <- apply(sapply(llr0$results,function(i) i$coef),1,var)
  # ci.llr0[i,] <- rowMeans(sapply(llr0$results,function(i) i$ci[,2]-i$ci[,1]))
  # cov.llr0[i,] <- rowMeans(sapply(llr0$results,function(i) i$ci[,2])-ate>0 & sapply(llr0$results,function(i) i$ci[,1])-ate<0)
  # zero.llr0[i,] <- rowMeans(sapply(llr0$results,function(i) i$ci[,2])>0 & sapply(llr0$results,function(i) i$ci[,1])<0)
  ##
  llr.ate[[i]] <- t(sapply(llr1$results,function(i) c(i$coef[1],i$ci[1,])))
  rmse.llr1[i,] <- sqrt(mean((llr.ate[[i]][,1]-ate)^2))
  sb.llr1[i,] <- mean(llr.ate[[i]][,1]-ate)^2
  v.llr1[i,] <- var(llr.ate[[i]][,1])
  ci.llr1[i,] <- mean(llr.ate[[i]][,3]-llr.ate[[i]][,2])
  cov.llr1[i,] <- mean(llr.ate[[i]][,3]-ate>0 & llr.ate[[i]][,2]-ate<0)
  zero.llr1[i,] <- mean(llr.ate[[i]][,3]>0 & llr.ate[[i]][,2]<0)
  ##
  # rmse.llr2[i,] <- rowMeans((sapply(llr2$results,function(i) i$coef)-ate)^2)
  # sb.llr2[i,] <- rowMeans(sapply(llr2$results,function(i) i$coef)-ate)^2
  # v.llr2[i,] <- apply(sapply(llr2$results,function(i) i$coef),1,var)
  # ci.llr2[i,] <- rowMeans(sapply(llr2$results,function(i) i$ci[,2]-i$ci[,1]))
  # cov.llr2[i,] <- rowMeans(sapply(llr2$results,function(i) i$ci[,2])-ate>0 & sapply(llr2$results,function(i) i$ci[,1])-ate<0)
  # zero.llr2[i,] <- rowMeans(sapply(llr2$results,function(i) i$ci[,2])>0 & sapply(llr2$results,function(i) i$ci[,1])<0)
  ##
  # cgs.res <- t(apply(sapply(cgs0$results,function(i) i$atem),2,function(i) c(mean(i),quantile(i,c(0.025,0.975)))))
  # rmse.cgs0[i] <- sqrt(mean((cgs.res[,1]-ate)^2))
  # sb.cgs0[i] <- mean(cgs.res[,1]-ate)^2
  # v.cgs0[i] <- var(cgs.res[,1])
  # ci.cgs0[i] <- mean(cgs.res[,3]-cgs.res[,2])
  # cov.cgs0[i] <- mean(cgs.res[,3]-ate>0 & cgs.res[,2]-ate<0)
  # zero.cgs0[i] <- mean(cgs.res[,3]>0 & cgs.res[,2]<0)
  ##
  cgs.ate[[i]] <- t(sapply(cgs1$results, function(i) c(mean(i$atem),quantile(i$atem,c(0.025,0.975)))))
  cgs.res <- t(apply(sapply(cgs1$results,function(i) i$atem),2,function(i) c(mean(i),quantile(i,c(0.025,0.975)))))
  rmse.cgs1[i] <- sqrt(mean((cgs.res[,1]-ate)^2))
  sb.cgs1[i] <- mean(cgs.res[,1]-ate)^2
  v.cgs1[i] <- var(cgs.res[,1])
  ci.cgs1[i] <- mean(cgs.res[,3]-cgs.res[,2])
  cov.cgs1[i] <- mean(cgs.res[,3]-ate>0 & cgs.res[,2]-ate<0)
  zero.cgs1[i] <- mean(cgs.res[,3]>0 & cgs.res[,2]<0)
}
##
names <- c("BART-RDD","BCF","S-BART","T-BART","CGS","LLR")
###
rmse <- cbind(rmse.bart,rmse.bcf,rmse.sbart,rmse.tbart,rmse.cgs1,rmse.llr1[,1])
colnames(rmse) <- names
table(apply(rmse,1,function(i) which(i==min(i))))
###
sb <- cbind(sb.bart,sb.bcf,sb.sbart,sb.tbart,sb.cgs1,sb.llr1[,1])
colnames(sb) <- names
table(apply(sb,1,function(i) which(i==min(i))))
###
v <- cbind(v.bart,v.bcf,v.sbart,v.tbart,v.cgs1,v.llr1[,1])
colnames(v) <- names
table(apply(v,1,function(i) which(i==min(i))))
###
ci <- cbind(ci.bart,ci.bcf,ci.sbart,ci.tbart,ci.cgs1,ci.llr1[,1])
colnames(ci) <- names
table(apply(ci,1,function(i) which(i==min(i))))
###
cov <- cbind(cov.bart,cov.bcf,cov.sbart,cov.tbart,cov.cgs1,cov.llr1[,1])
colnames(cov) <- names
zero <- cbind(zero.bart,zero.bcf,zero.sbart,zero.tbart,zero.cgs1,zero.llr1[,1])
colnames(zero) <- names
###
bart.rdd.int <- sapply(lapply(bart.rdd.ate,function(i) apply(i,2,function(j) c(0.05-quantile(j,0.025),quantile(j,0.975)-0.05))),function(i) quantile(apply(i,2,function(j) ifelse(min(j)>0,0,abs(min(j)))),0.95))
bcf.int <- sapply(lapply(bcf.ate,function(i) apply(i,2,function(j) c(0.05-quantile(j,0.025),quantile(j,0.975)-0.05))),function(i) quantile(apply(i,2,function(j) ifelse(min(j)>0,0,abs(min(j)))),0.95))
sbart.int <- sapply(lapply(sbart.ate,function(i) apply(i,2,function(j) c(0.05-quantile(j,0.025),quantile(j,0.975)-0.05))),function(i) quantile(apply(i,2,function(j) ifelse(min(j)>0,0,abs(min(j)))),0.95))
tbart.int <- sapply(lapply(tbart.ate,function(i) apply(i,2,function(j) c(0.05-quantile(j,0.025),quantile(j,0.975)-0.05))),function(i) quantile(apply(i,2,function(j) ifelse(min(j)>0,0,abs(min(j)))),0.95))
llr.int <- sapply(llr.ate, function(i) quantile(apply(i,1,function(j) ifelse(0.05-j[2]>0 & j[3]-0.05>0,0,abs(min(0.05-j[2],j[3]-0.05)))),0.95))
cgs.int <- sapply(cgs.ate, function(i) quantile(apply(i,1,function(j) ifelse(0.05-j[2]>0 & j[3]-0.05>0,0,abs(min(0.05-j[2],j[3]-0.05)))),0.95))
int <- cbind(bart.rdd.int,bcf.int,sbart.int,tbart.int,cgs.int,llr.int)
colnames(int) <- names
###
rmse.cate <- cbind(bart.rdd.cate,bcf.cate,sbart.cate,tbart.cate)
colnames(rmse.cate) <- names[1:4]
table(apply(rmse.cate,1,function(i) which(i==min(i))))
###
cov.cate <- cbind(bart.rdd.cate.cov,bcf.cate.cov,sbart.cate.cov,tbart.cate.cov)
colnames(cov.cate) <- names[1:4]
###
ci.cate <- cbind(bart.rdd.cate.ci,bcf.cate.ci,sbart.cate.ci,tbart.cate.ci)
colnames(ci.cate) <- names[1:4]
###
save.image("Results/sims.RData")