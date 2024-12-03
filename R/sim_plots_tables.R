setwd("../")
load("Results/sims.RData")
###
# par(mfrow=c(2,2),bty="n")
# matplot(unlist(sapply(1:11, function(i) data$w[test[,i],i])),cbind(unlist(cate),unlist(bart.rdd.cate[[9]])),pch=19)
# matplot(unlist(sapply(1:11, function(i) data$w[test[,i],i])),cbind(unlist(cate),unlist(sbart.cate[[9]])),pch=19)
# matplot(unlist(sapply(1:11, function(i) data$w[test[,i],i])),cbind(unlist(cate),unlist(tbart.cate[[9]])),pch=19)
# matplot(unlist(sapply(1:11, function(i) data$w[test[,i],i])),cbind(unlist(cate),unlist(polynomial.cate[[9]])),pch=19)
###
rmse.plot <- cbind(params,rmse.cate)
for (i in c(6,7,8,10))
{
  rmse <- subset(rmse.plot,n==1000 & sig_error == 0.5)
  plot(rmse$delta.tau,rmse[,i],type="p",col=rmse$level+1,pch=rmse$delta.mu+15,
       ylab=colnames(rmse)[i],xlab=expression(delta[tau]),ylim=c(0,1))
  lines(unique(rmse$delta.tau),rmse[,i][rmse$delta.mu==0.5 & rmse$level==1],col=2)
  lines(unique(rmse$delta.tau),rmse[,i][rmse$delta.mu==0.5 & rmse$level==2],col=3)
  lines(unique(rmse$delta.tau),rmse[,i][rmse$delta.mu==2 & rmse$level==1],col=2)
  lines(unique(rmse$delta.tau),rmse[,i][rmse$delta.mu==2 & rmse$level==2],col=3)
  grid(nx = NA,
       ny = NULL,
       lty = 2, col = "gray", lwd = 1)
  points(rmse$delta.tau,rmse$Oracle,type="p",col="gray",pch=rmse$delta.mu+15,
         ylab=colnames(rmse)[i],xlab=expression(delta[tau]))
  lines(unique(rmse$delta.tau),rmse$Oracle[rmse$delta.mu==0.5 & rmse$level==1],col="gray")
  lines(unique(rmse$delta.tau),rmse$Oracle[rmse$delta.mu==0.5 & rmse$level==2],col="gray")
  lines(unique(rmse$delta.tau),rmse$Oracle[rmse$delta.mu==2 & rmse$level==1],col="gray")
  lines(unique(rmse$delta.tau),rmse$Oracle[rmse$delta.mu==2 & rmse$level==2],col="gray")
}
mtext(expression(tau(w) == sin(10*pi*w)), side = 3, line = - 2, outer = TRUE)
###
sample <- 20
dgp <- 3
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