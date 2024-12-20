setwd("~/Git/BART-RDD")
load("Results/sims.RData")
###
rmse.plot <- cbind(params,rmse.cate/ate.cate)
par(mfrow=c(2,2))
for (i in c(500,1000,2500,5000))
{
  rmse <- subset(rmse.plot,n==i)
  matplot(t(rmse[,-c(1:4)]),bty="n",pch=16,col=1:3,
          lty=2,type="b",ylab="RMSE/RMSE(ATE)",
          xaxt="n",main=bquote(N==.(i)))
  axis(1,at=1:(ncol(rmse)-4),labels=names(rmse)[-c(1:4)],cex.axis=0.8)
  grid(nx = NA,
       ny = NULL,
       lty = 2, col = "gray", lwd = 1)
  if (i==500)
  {
    legend("topleft",col=1:3,legend=unique(params$delta.tau),
           title=expression(tau~"'/"~mu~"'"),ncol=2,lty=2,cex=0.75,lwd=2,
           text.width = 0.4*c(strwidth("0.3"),strwidth("0.5"),strwidth("0.7")),
           y.intersp = 0.5)
  }
}
###
# sample <- 1
dgp <- 10
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
                        `T-BART`=tbart.cate.est, Oracle = oracle.cate.est,
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