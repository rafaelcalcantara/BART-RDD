setwd("../")
load("Results/sims.RData")
###
rmse.plot <- cbind(params,rmse.cate)
dt <- 1
par(mfrow=c(2,2))
for (i in c(500,1000,2500,5000))
{
  rmse <- subset(rmse.plot,n==i & delta.tau==dt)
  matplot(t(rmse[,-c(1:5)])/dt,bty="n",pch=16,col=rmse$level+1,
          lty=2,type="b",ylab=bquote("RMSE/"~delta[tau]),
          xaxt="n",main=bquote(N==.(i)))
  axis(1,at=1:(ncol(rmse)-5),labels=names(rmse)[-c(1:5)],cex.axis=0.8)
  grid(nx = NA,
       ny = NULL,
       lty = 2, col = "gray", lwd = 1)
  if (i==500)
  {
    legend(x=par("usr")[1],y=par("usr")[4],col=2:3,legend=c(1,2),
           title=bquote(tau),ncol=1,lty=2,cex=0.75,lwd=2,
           text.width = 0.2*c(strwidth("1"),strwidth("2")),y.intersp = 0.75)
  }
}
###
sample <- 1000
dgp <- 2
Owidth <- Ow[1]
bart.rdd.cate.est <- bart.rdd.cate[[dgp]][[sample]]
sbart.cate.est <- sbart.cate[[dgp]][[sample]]
tbart.cate.est <- tbart.cate[[dgp]][[sample]]
oracle.cate.est <- oracle.cate[[dgp]][[sample]]
polynomial.1.cate.est <- polynomial.1.cate[[dgp]][[sample]]
polynomial.2.cate.est <- polynomial.2.cate[[dgp]][[sample]]
polynomial.3.cate.est <- polynomial.3.cate[[dgp]][[sample]]
polynomial.full.cate.est <- polynomial.full.cate[[dgp]][[sample]]
cate.estimates <- cbind(`BART-RDD`=bart.rdd.cate.est,`S-BART`=sbart.cate.est,
                        `T-BART`=tbart.cate.est, Oracle = oracle.cate.est,
                        Polynomial.1 = polynomial.1.cate.est,
                        Polynomial.2 = polynomial.2.cate.est,Polynomial.3 = polynomial.3.cate.est,
                        Polynomial.full = polynomial.full.cate.est)
cate.estimates <- as.data.frame(cate.estimates)
data <- readRDS(paste0("Data/dgp_",dgp,".rds"))
test <- -Owidth <= data$x[,sample] & data$x[,sample] <= Owidth
w <- data$w[,sample][test]
cate <- data$tau.x
par(mfrow=c(2,2))
for (i in names(cate.estimates))
{
  matplot(w,cbind(cate[test,sample],cate.estimates[[i]]),pch=19,ylab=i)
}