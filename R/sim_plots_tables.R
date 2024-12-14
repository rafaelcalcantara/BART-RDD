setwd("../")
load("Results/sims.RData")
###
# par(mfrow=c(2,2),bty="n")
# matplot(unlist(sapply(1:11, function(i) data$w[test[,i],i])),cbind(unlist(cate),unlist(bart.rdd.cate[[9]])),pch=19)
# matplot(unlist(sapply(1:11, function(i) data$w[test[,i],i])),cbind(unlist(cate),unlist(sbart.cate[[9]])),pch=19)
# matplot(unlist(sapply(1:11, function(i) data$w[test[,i],i])),cbind(unlist(cate),unlist(tbart.cate[[9]])),pch=19)
# matplot(unlist(sapply(1:11, function(i) data$w[test[,i],i])),cbind(unlist(cate),unlist(polynomial.cate[[9]])),pch=19)
###
# rmse.plot <- cbind(params,rmse.cate)
# par(mfrow=c(2,2),bty="n")
# for (i in c(6,7,8,10))
# {
#   rmse <- subset(rmse.plot,n==5000)
#   plot(rmse$delta.tau,rmse[,i],type="p",col=rmse$level+1,pch=rmse$sig_error+15,
#        ylab=colnames(rmse)[i],xlab=expression(delta[tau]),ylim=c(0,0.9))
#   # lines(unique(rmse$delta.tau),rmse[,i][rmse$sig_error==0.25 & rmse$level==1],col=2)
#   # lines(unique(rmse$delta.tau),rmse[,i][rmse$sig_error==0.25 & rmse$level==2],col=3)
#   # lines(unique(rmse$delta.tau),rmse[,i][rmse$sig_error==0.5 & rmse$level==1],col=2)
#   # lines(unique(rmse$delta.tau),rmse[,i][rmse$sig_error==0.5 & rmse$level==2],col=3)
#   lines(unique(rmse$delta.tau),rmse[,i][rmse$sig_error==1 & rmse$level==1],col=2)
#   lines(unique(rmse$delta.tau),rmse[,i][rmse$sig_error==1 & rmse$level==2],col=3)
#   lines(unique(rmse$delta.tau),rmse[,i][rmse$sig_error==2 & rmse$level==1],col=2)
#   lines(unique(rmse$delta.tau),rmse[,i][rmse$sig_error==2 & rmse$level==2],col=3)
#   grid(nx = NA,
#        ny = NULL,
#        lty = 2, col = "gray", lwd = 1)
#   points(rmse$delta.tau,rmse$Oracle,type="p",col="gray",pch=rmse$sig_error+15,
#          ylab=colnames(rmse)[i],xlab=expression(delta[tau]))
#   # lines(unique(rmse$delta.tau),rmse$Oracle[rmse$sig_error==0.25 & rmse$level==1],col="gray")
#   # lines(unique(rmse$delta.tau),rmse$Oracle[rmse$sig_error==0.25 & rmse$level==2],col="gray")
#   # lines(unique(rmse$delta.tau),rmse$Oracle[rmse$sig_error==0.5 & rmse$level==1],col="gray")
#   # lines(unique(rmse$delta.tau),rmse$Oracle[rmse$sig_error==0.5 & rmse$level==2],col="gray")
#   lines(unique(rmse$delta.tau),rmse$Oracle[rmse$sig_error==1 & rmse$level==1],col="gray")
#   lines(unique(rmse$delta.tau),rmse$Oracle[rmse$sig_error==1 & rmse$level==2],col="gray")
#   lines(unique(rmse$delta.tau),rmse$Oracle[rmse$sig_error==2 & rmse$level==1],col="gray")
#   lines(unique(rmse$delta.tau),rmse$Oracle[rmse$sig_error==2 & rmse$level==2],col="gray")
# }
# mtext(expression(tau(w) == sin(10*pi*w)), side = 3, line = - 2, outer = TRUE)
###
rmse.plot <- cbind(params,rmse.cate)
par(mfrow=c(2,2))
for (i in c(500,1000,2500,5000))
{
  rmse <- subset(rmse.plot,n==i & delta.tau==1)
  matplot(t(rmse[,-c(1:5)]),bty="n",pch=rmse$sig_error+14,col=rmse$level+1,
          lty=2,type="b",ylab="RMSE",
          xaxt="n",main=bquote(N==.(i)))
  axis(1,at=1:(ncol(rmse)-5),labels=names(rmse)[-c(1:5)],cex.axis=0.8)
  grid(nx = NA,
       ny = NULL,
       lty = 2, col = "gray", lwd = 1)
  if (i==500)
  {
    legend("topleft",col=2:3,legend=c(1,2),title=bquote(tau),ncol=2,lty=2,cex=0.75,lwd=2)
    legend("top",pch=15:16,legend=c(1,2),title=bquote(sigma),ncol=2,cex=0.75)
  }
}
###
sample <- 3
dgp <- 2
Owidth <- Ow[2]
bart.rdd.cate.est <- bart.rdd.cate[[dgp]][[sample]]
sbart.cate.est <- sbart.cate[[dgp]][[sample]]
tbart.cate.est <- tbart.cate[[dgp]][[sample]]
oracle.cate.est <- oracle.cate[[dgp]][[sample]]
polynomial.cate.est <- polynomial.cate[[dgp]][[sample]]
cate.estimates <- cbind(`BART-RDD`=bart.rdd.cate.est,`S-BART`=sbart.cate.est,
                        `T-BART`=tbart.cate.est, Polynomial = polynomial.cate.est)
cate.estimates <- as.data.frame(cate.estimates)
data <- readRDS(paste0("Data/dgp_",dgp,".rds"))
test <- -Owidth <= data$x[,sample] & data$x[,sample] <= Owidth
w <- data$w[,sample][test]
cate <- data$tau.x
for (i in names(cate.estimates))
{
  matplot(w,cbind(cate[test,sample],cate.estimates[[i]]),pch=19,ylab=i)
}