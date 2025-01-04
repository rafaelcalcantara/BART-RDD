setwd("~/Git/BART-RDD")
load("Results/sims.RData")
###
for (i in 1:nrow(params))
{
  rmse <- do.call("rbind",lapply(mget(ls(pattern = ".cate.rmse")),"[",i,))
  rmse <- rmse[c(2,1,5,3,4),]
  axis.names <- rownames(rmse)
  axis.names <- gsub("\\.cate.rmse","",axis.names)
  axis.names <- gsub("\\.rdd","-rdd",axis.names)
  axis.names <- gsub("\\.full"," full",axis.names)
  axis.names <- toupper(axis.names)
  matplot(rmse,bty="n",pch=19,col=1,lty=c(2,1,2),type=c("l","b","l"),
          ylab="RMSE/RMSE (ATE)",xaxt="n",
          main=paste(names(params),params[i,],sep=": "))
  axis(1,at=1:nrow(rmse),labels=axis.names,cex.axis=0.8)
  grid(nx = NA,
       ny = NULL,
       lty = 2, col = "gray", lwd = 1)
}
###
dgp <- 15
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