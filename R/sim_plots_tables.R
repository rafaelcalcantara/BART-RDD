setwd("~/Git/BART-RDD")
load("Results/sims.RData")
###
rmse.plot <- cbind(params,rmse.cate/ate.cate)
par(mfrow=c(2,2))
for (i in c(500,1000,1500))
{
  rmse <- subset(rmse.plot,n==i)
  matplot(t(rmse[,-c(1:2,7)]),bty="n",pch=19,col=1:2,
          lty=2,type="b",ylab="RMSE/RMSE (ATE)",
          xaxt="n",main=bquote(N==.(i)))
  axis(1,at=1:(ncol(rmse)-3),labels=names(rmse)[-c(1:2,7)],cex.axis=0.8)
  grid(nx = NA,
       ny = NULL,
       lty = 2, col = "gray", lwd = 1)
  # matplot(t(rmse[,-c(1:2)]),bty="n",pch=19,col=1:2,
  #         lty=2,type="b",ylab="RMSE/RMSE (ATE)",
  #         xaxt="n",main=bquote(N==.(i)))
  # axis(1,at=1:(ncol(rmse)-2),labels=names(rmse)[-c(1:2)],cex.axis=0.8)
  # grid(nx = NA,
  #      ny = NULL,
  #      lty = 2, col = "gray", lwd = 1)
}
plot.new()
legend("center",col=1:3,legend=c("0.75","0.95"),
       title=bquote(rho),ncol=2,lty=1,pch=19,cex=0.75,lwd=2)
###
# sample <- 1
dgp <- 6
# bart.rdd.cate.est <- bart.rdd.cate[[dgp]][[sample]]
# sbart.cate.est <- sbart.cate[[dgp]][[sample]]
# tbart.cate.est <- tbart.cate[[dgp]][[sample]]
# oracle.cate.est <- oracle.cate[[dgp]][[sample]]
# polynomial.1.cate.est <- polynomial.1.cate[[dgp]][[sample]]
# polynomial.2.cate.est <- polynomial.2.cate[[dgp]][[sample]]
# polynomial.3.cate.est <- polynomial.3.cate[[dgp]][[sample]]
# polynomial.full.cate.est <- polynomial.full.cate[[dgp]][[sample]]
bart.rdd.cate.est <- unlist(bart.rdd.cate[[dgp]])
sbart.cate.est <- unlist(sbart.cate[[dgp]])
tbart.cate.est <- unlist(tbart.cate[[dgp]])
oracle.cate.est <- unlist(oracle.cate[[dgp]])
polynomial.full.cate.est <- unlist(polynomial.full.cate[[dgp]])
cate.estimates <- cbind(`BART-RDD`=bart.rdd.cate.est,`S-BART`=sbart.cate.est,
                        `T-BART`=tbart.cate.est,
                        # Oracle = oracle.cate.est,
                        # Polynomial.1 = polynomial.1.cate.est,
                        # Polynomial.2 = polynomial.2.cate.est,Polynomial.3 = polynomial.3.cate.est,
                        Polynomial.full = polynomial.full.cate.est)
cate.estimates <- as.data.frame(cate.estimates)
data <- readRDS(paste0("Data/dgp_",dgp,".rds"))
Owidth <- data$h
test <- sapply(samples,function(i) -Owidth[i] <= data$x[,i] & data$x[,i] <= Owidth[i])
w <- unlist(sapply(samples,function(i) data$w[test[,i],i]))
cate <- data$tau.x
cate <- unlist(sapply(samples,function(i) cate[test[,i],i]))
par(mfrow=c(2,2))
for (i in names(cate.estimates))
{
  matplot(w,cbind(cate.estimates[[i]],cate),pch=19,ylab=i,col=2:1)
}