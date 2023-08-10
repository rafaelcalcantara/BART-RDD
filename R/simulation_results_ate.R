## Setup
library(lattice)
s      <- 500
sample <- c(500,1000)
model  <- 3:6
xi <- nu <- kappa <- c(0.25,2)
i <- 1
## BART-RDD
res.mat <- as.data.frame(matrix(0,32,5,dimnames=list(1:32,c("BART-RDD","Model","Xi","Nu","Kappa"))))
rmse.bart.rdd <- res.mat
cov.bart.rdd <- res.mat
length.bart.rdd <- res.mat
alpha.bart.rdd <- res.mat
interval.bart.rdd <- cbind(LI=0,UI=0,res.mat[,-1])
distance.bart.rdd <- matrix(0,s,32)
index <- 0
for (j in 1:length(model))
{
    for (k in 1:length(xi))
    {
        for (l in 1:length(nu))
        {
            for (m in 1:length(kappa))
            {
                index <- index+1
                dgp <- c(model[j],xi[k],nu[l],kappa[m])
                file <- paste0("Results/bart_rdd_",sample[i],"_",model[j],"_",xi[k],"_",nu[l],"_",kappa[m])
                bart.rdd <- readRDS(paste0(file,".rds"))
                bart.rdd <- lapply(bart.rdd, function(x) x$pred)
                a <- sapply(bart.rdd, function(x) mean(x<=xi[k]))
                alpha.bart.rdd[index,] <- c(mean(a),dgp)
                bart.rdd <- t(sapply(bart.rdd,function(x) c(mean(colMeans(x)),quantile(colMeans(x),c(0.025,0.975)))))
                rmse.bart.rdd[index,] <- c(sqrt(mean((bart.rdd[,1]-xi[k])^2)),dgp)
                cov.bart.rdd[index,] <- c(mean(bart.rdd[,2] <= xi[k] & xi[k] <= bart.rdd[,3]),dgp)
                length.bart.rdd[index,] <- c(mean(bart.rdd[,3] - bart.rdd[,2]),dgp)
                interval.bart.rdd[index,] <- c(mean(bart.rdd[,2]),mean(bart.rdd[,3]),dgp)
                distance.bart.rdd[,index] <- apply(bart.rdd,1,function(x) ifelse(x[2]<=xi[k] & xi[k]<=x[3],0,min(abs(x[2]-xi[k]),abs(x[3]-xi[k]))))
            }
        }
    }
}
interval.bart.rdd <- reshape(interval.bart.rdd,direction="long",varying=c("LI","UI"),v.names="BART-RDD",idvar=c("Model","Xi","Nu","Kappa"),timevar="Interval",times=c("LI","UI"))
## CKT
res.mat <- as.data.frame(matrix(0,32,5,dimnames=list(1:32,c("CKT","Model","Xi","Nu","Kappa"))))
rmse.ckt <- res.mat
cov.ckt <- res.mat
length.ckt <- res.mat
interval.ckt <- cbind(LI=0,UI=0,res.mat[,-1])
distance.ckt <- matrix(0,s,32)
index <- 0
for (j in 1:length(model))
{
    for (k in 1:length(xi))
    {
        for (l in 1:length(nu))
        {
            for (m in 1:length(kappa))
            {
                index <- index+1
                dgp <- c(model[j],xi[k],nu[l],kappa[m])
                file <- paste0("Results/cct_",sample[i],"_",model[j],"_",xi[k],"_",nu[l],"_",kappa[m])
                ckt <- readRDS(paste0(file,".rds"))
                ckt <- lapply(ckt, function(x) x$pred)
                ckt <- do.call("rbind",ckt)[1:s,]
                rmse.ckt[index,] <- c(sqrt(mean((ckt[,1]-xi[k])^2)),dgp)
                cov.ckt[index,] <- c(mean(ckt[,2] <= xi[k] & xi[k] <= ckt[,3]),dgp)
                length.ckt[index,] <- c(mean(ckt[,3] - ckt[,2]),dgp)
                interval.ckt[index,] <- c(mean(ckt[,2]),mean(ckt[,3]),dgp)
                distance.ckt[,index] <- apply(ckt,1,function(x) ifelse(x[2]<=xi[k] & xi[k]<=x[3],0,min(abs(x[2]-xi[k]),abs(x[3]-xi[k]))))
            }
        }
    }
}
interval.ckt <- reshape(interval.ckt,direction="long",varying=c("LI","UI"),v.names="CKT",idvar=c("Model","Xi","Nu","Kappa"),timevar="Interval",times=c("LI","UI"))
## KR
res.mat <- as.data.frame(matrix(0,32,5,dimnames=list(1:32,c("KR","Model","Xi","Nu","Kappa"))))
rmse.kr <- res.mat
cov.kr <- res.mat
length.kr <- res.mat
interval.kr <- cbind(LI=0,UI=0,res.mat[,-1])
distance.kr <- matrix(0,s,32)
index <- 0
for (j in 1:length(model))
{
    for (k in 1:length(xi))
    {
        for (l in 1:length(nu))
        {
            for (m in 1:length(kappa))
            {
                index <- index+1
                dgp <- c(model[j],xi[k],nu[l],kappa[m])
                file <- paste0("Results/kr_",sample[i],"_",model[j],"_",xi[k],"_",nu[l],"_",kappa[m])
                kr <- readRDS(paste0(file,".rds"))
                dt <- lapply(kr, function(x) x$dt)
                kr <- lapply(kr, function(x) x$pred)
                kr <- do.call("rbind",kr)[1:s,]
                rmse.kr[index,] <- c(sqrt(mean((kr[,1]-xi[k])^2)),dgp)
                cov.kr[index,] <- c(mean(kr[,2] <= xi[k] & xi[k] <= kr[,3]),dgp)
                length.kr[index,] <- c(mean(kr[,3] - kr[,2]),dgp)
                interval.kr[index,] <- c(mean(kr[,2]),mean(kr[,3]),dgp)
                distance.kr[,index] <- apply(kr,1,function(x) ifelse(x[2]<=xi[k] & xi[k]<=x[3],0,min(abs(x[2]-xi[k]),abs(x[3]-xi[k]))))
            }
        }
    }
}
interval.kr <- reshape(interval.kr,direction="long",varying=c("LI","UI"),v.names="KR",idvar=c("Model","Xi","Nu","Kappa"),timevar="Interval",times=c("LI","UI"))
## BART1
res.mat <- as.data.frame(matrix(0,32,5,dimnames=list(1:32,c("BART1","Model","Xi","Nu","Kappa"))))
rmse.bart.1 <- res.mat
cov.bart.1 <- res.mat
length.bart.1 <- res.mat
alpha.bart.1 <- res.mat
interval.bart.1 <- cbind(LI=0,UI=0,res.mat[,-1])
distance.bart.1 <- matrix(0,s,32)
index <- 0
for (j in 1:length(model))
{
    for (k in 1:length(xi))
    {
        for (l in 1:length(nu))
        {
            for (m in 1:length(kappa))
            {
                index <- index+1
                dgp <- c(model[j],xi[k],nu[l],kappa[m])
                file <- paste0("Results/bart1_",sample[i],"_",model[j],"_",xi[k],"_",nu[l],"_",kappa[m])
                bart.1 <- readRDS(paste0(file,".rds"))
                bart.1 <- lapply(bart.1, function(x) x$pred)
                a <- sapply(bart.1, function(x) mean(x<=xi[k]))
                alpha.bart.1[index,] <- c(mean(a),dgp)
                bart.1 <- t(sapply(bart.1,function(x) c(mean(colMeans(x)),quantile(colMeans(x),c(0.025,0.975)))))[1:s,]
                rmse.bart.1[index,] <- c(sqrt(mean((bart.1[,1]-xi[k])^2)),dgp)
                cov.bart.1[index,] <- c(mean(bart.1[,2] <= xi[k] & xi[k] <= bart.1[,3]),dgp)
                length.bart.1[index,] <- c(mean(bart.1[,3] - bart.1[,2]),dgp)
                interval.bart.1[index,] <- c(mean(bart.1[,2]),mean(bart.1[,3]),dgp)
                distance.bart.1[,index] <- apply(bart.1,1,function(x) ifelse(x[2]<=xi[k] & xi[k]<=x[3],0,min(abs(x[2]-xi[k]),abs(x[3]-xi[k]))))
            }
        }
    }
}
interval.bart.1 <- reshape(interval.bart.1,direction="long",varying=c("LI","UI"),v.names="BART1",idvar=c("Model","Xi","Nu","Kappa"),timevar="Interval",times=c("LI","UI"))
## BART2
res.mat <- as.data.frame(matrix(0,32,5,dimnames=list(1:32,c("BART2","Model","Xi","Nu","Kappa"))))
rmse.bart.2 <- res.mat
cov.bart.2 <- res.mat
length.bart.2 <- res.mat
alpha.bart.2 <- res.mat
interval.bart.2 <- cbind(LI=0,UI=0,res.mat[,-1])
distance.bart.2 <- matrix(0,s,32)
index <- 0
for (j in 1:length(model))
{
    for (k in 1:length(xi))
    {
        for (l in 1:length(nu))
        {
            for (m in 1:length(kappa))
            {
                index <- index+1
                dgp <- c(model[j],xi[k],nu[l],kappa[m])
                file <- paste0("Results/bart2_",sample[i],"_",model[j],"_",xi[k],"_",nu[l],"_",kappa[m])
                bart.2 <- readRDS(paste0(file,".rds"))
                bart.2 <- lapply(bart.2, function(x) x$pred)
                a <- sapply(bart.2, function(x) mean(x<=xi[k]))
                alpha.bart.2[index,] <- c(mean(a),dgp)
                bart.2 <- t(sapply(bart.2,function(x) c(mean(x),quantile(x,c(0.025,0.975)))))[1:s,]
                rmse.bart.2[index,] <- c(sqrt(mean((bart.2[,1]-xi[k])^2)),dgp)
                cov.bart.2[index,] <- c(mean(bart.2[,2] <= xi[k] & xi[k] <= bart.2[,3]),dgp)
                length.bart.2[index,] <- c(mean(bart.2[,3] - bart.2[,2]),dgp)
                interval.bart.2[index,] <- c(mean(bart.2[,2]),mean(bart.2[,3]),dgp)
                distance.bart.2[,index] <- apply(bart.2,1,function(x) ifelse(x[2]<=xi[k] & xi[k]<=x[3],0,min(abs(x[2]-xi[k]),abs(x[3]-xi[k]))))
            }
        }
    }
}
interval.bart.2 <- reshape(interval.bart.2,direction="long",varying=c("LI","UI"),v.names="BART2",idvar=c("Model","Xi","Nu","Kappa"),timevar="Interval",times=c("LI","UI"))
## BCF
bcf <- vector("list",32)
res.mat <- as.data.frame(matrix(0,32,5,dimnames=list(1:32,c("BCF","Model","Xi","Nu","Kappa"))))
rmse.bcf <- res.mat
cov.bcf <- res.mat
length.bcf <- res.mat
alpha.bcf <- res.mat
interval.bcf <- cbind(LI=0,UI=0,res.mat[,-1])
distance.bcf <- matrix(0,s,32)
index <- 0
for (j in 1:length(model))
{
    for (k in 1:length(xi))
    {
        for (l in 1:length(nu))
        {
            for (m in 1:length(kappa))
            {
                index <- index+1
                dgp <- c(model[j],xi[k],nu[l],kappa[m])
                file <- paste0("Results/bcf_",sample[i],"_",model[j],"_",xi[k],"_",nu[l],"_",kappa[m])
                bcf <- readRDS(paste0(file,".rds"))
                bcf <- lapply(bcf, function(x) x$pred)
                a <- sapply(bcf, function(x) mean(x<=xi[k]))
                alpha.bcf[index,] <- c(mean(a),dgp)
                bcf <- t(sapply(bcf,function(x) c(mean(colMeans(x)),quantile(colMeans(x),c(0.025,0.975)))))[1:s,]
                rmse.bcf[index,] <- c(sqrt(mean((bcf[,1]-xi[k])^2)),dgp)
                cov.bcf[index,] <- c(mean(bcf[,2] <= xi[k] & xi[k] <= bcf[,3]),dgp)
                length.bcf[index,] <- c(mean(bcf[,3] - bcf[,2]),dgp)
                interval.bcf[index,] <- c(mean(bcf[,2]),mean(bcf[,3]),dgp)
                distance.bcf[,index] <- apply(bcf,1,function(x) ifelse(x[2]<=xi[k] & xi[k]<=x[3],0,min(abs(x[2]-xi[k]),abs(x[3]-xi[k]))))
            }
        }
    }
}
interval.bcf <- reshape(interval.bcf,direction="long",varying=c("LI","UI"),v.names="BCF",idvar=c("Model","Xi","Nu","Kappa"),timevar="Interval",times=c("LI","UI"))
## Merge results
load("Tables/cgs.RData")
rmse <- merge(rmse.bart.rdd,rmse.bart.1)
rmse <- merge(rmse,rmse.bart.2)
rmse <- merge(rmse,rmse.bcf)
rmse <- merge(rmse,rmse.ckt)
rmse <- merge(rmse,rmse.kr)
rmse <- merge(rmse,rmse.cgs)
###
cov <- merge(cov.bart.rdd,cov.bart.1)
cov <- merge(cov,cov.bart.2)
cov <- merge(cov,cov.bcf)
cov <- merge(cov,cov.ckt)
cov <- merge(cov,cov.kr)
cov <- merge(cov,cov.cgs)
###
length <- merge(length.bart.rdd,length.bart.1)
length <- merge(length,length.bart.2)
length <- merge(length,length.bcf)
length <- merge(length,length.ckt)
length <- merge(length,length.kr)
length <- merge(length,length.cgs)
###
alpha <- merge(alpha.bart.rdd,alpha.bart.1)
alpha <- merge(alpha,alpha.bart.2)
alpha <- merge(alpha,alpha.bcf)
alpha <- merge(alpha,alpha.cgs)
###
interval <- merge(interval.bart.rdd,interval.bart.1)
interval <- merge(interval,interval.bart.2)
interval <- merge(interval,interval.bcf)
interval <- merge(interval,interval.ckt)
interval <- merge(interval,interval.kr)
## interval <- merge(interval,interval.cgs)
## Plots
plot.rmse <- reshape(rmse,direction="long",varying=c("BART-RDD","BART1","BART2","BCF","CKT","KR","CGS"),v.names="RMSE",idvar=c("Model","Xi","Nu","Kappa"),timevar="Method",times=c("BART-RDD","BART1","BART2","BCF","CKT","KR","CGS"))
plot.cov <- reshape(cov,direction="long",varying=c("BART-RDD","BART1","BART2","BCF","CKT","KR","CGS"),v.names="Coverage",idvar=c("Model","Xi","Nu","Kappa"),timevar="Method",times=c("BART-RDD","BART1","BART2","BCF","CKT","KR","CGS"))
plot.length <- reshape(length,direction="long",varying=c("BART-RDD","BART1","BART2","BCF","CKT","KR","CGS"),v.names="Length",idvar=c("Model","Xi","Nu","Kappa"),timevar="Method",times=c("BART-RDD","BART1","BART2","BCF","CKT","KR","CGS"))
plot.interval <- reshape(interval,direction="long",varying=c("BART-RDD","BART1","BART2","BCF","CKT","KR"),v.names="Value",idvar=c("Model","Xi","Nu","Kappa","Interval"),timevar="Method",times=c("BART-RDD","BART1","BART2","BCF","CKT","KR"))
plot.interval$Method <- as.factor(plot.interval$Method)
### RMSE
plot.rmse$Method <- as.factor(plot.rmse$Method)
png("Figures/rmse1.png")
xyplot(RMSE~Method|factor(Nu,labels=paste("Nu=",c(0.25,2),sep=""))+factor(Kappa,labels=paste("Kappa=",c(0.25,2),sep="")),data=subset(plot.rmse,Model %in% 3:4),groups=factor(Xi,labels=paste("Xi=",c(0.25,2),sep="")),auto.key=list(column=2,cex=0.75),scales=list(cex=0.5),pch=21)
dev.off()
png("Figures/rmse2.png")
xyplot(RMSE~Method|factor(Nu,labels=paste("Nu=",c(0.25,2),sep=""))+factor(Kappa,labels=paste("Kappa=",c(0.25,2),sep="")),data=subset(plot.rmse,Model %in% 5:6),groups=factor(Xi,labels=paste("Xi=",c(0.25,2),sep="")),auto.key=list(column=2,cex=0.75),scales=list(cex=0.5),pch=21)
dev.off()
### Cov x Length
p <- merge(plot.cov,plot.length)
p$Length <- log(p$Length)
png("Figures/coverage1.png")
xyplot(Coverage~Length|factor(Nu,labels=paste("Nu=",c(0.25,2),sep=""))+factor(Kappa,labels=paste("Kappa=",c(0.25,2),sep="")),data=subset(p,Model %in% 3:4),groups=factor(Method),auto.key=list(column=7,cex=0.5),scales=list(cex=0.5),pch=21)
dev.off()
png("Figures/coverage2.png")
xyplot(Coverage~Length|factor(Nu,labels=paste("Nu=",c(0.25,2),sep=""))+factor(Kappa,labels=paste("Kappa=",c(0.25,2),sep="")),data=subset(p,Model %in% 5:6),groups=factor(Method),auto.key=list(column=7,cex=0.5),scales=list(cex=0.5),pch=21)
dev.off()
### Interval inflation
d1 <- data.frame(Method="BART-RDD",delta=apply(distance.bart.rdd,2,quantile,0.95))
d2 <- data.frame(Method="BART1",delta=apply(distance.bart.1,2,quantile,0.95))
d3 <- data.frame(Method="BART2",delta=apply(distance.bart.2,2,quantile,0.95))
d4 <- data.frame(Method="BCF",delta=apply(distance.bcf,2,quantile,0.95))
d5 <- data.frame(Method="CKT",delta=apply(distance.ckt,2,quantile,0.95))
d6 <- data.frame(Method="KR",delta=apply(distance.kr,2,quantile,0.95))
d7 <- data.frame(Method="CGS",delta=apply(distance.cgs,2,quantile,0.95))
d <- rbind(d1,d2,d3,d4,d5,d6,d7)
d <- cbind(d,subset(plot.length,select=c("Method","Model","Xi","Nu","Kappa","Length")))
d$`Corrected Interval Length` <- log(2*d$delta+d$Length)
d$Method <- as.factor(d$Method)
png("Figures/length1.png")
xyplot(`Corrected Interval Length`~Method|factor(Nu,labels=paste("Nu=",c(0.25,2),sep=""))+factor(Kappa,labels=paste("Kappa=",c(0.25,2),sep="")),data=subset(d,Model %in% 3:4),groups=factor(Xi,labels=paste("Xi=",c(0.25,2),sep="")),auto.key=list(column=2,cex=0.75),scales=list(cex=0.5),pch=21)
dev.off()
png("Figures/length2.png")
xyplot(`Corrected Interval Length`~Method|factor(Nu,labels=paste("Nu=",c(0.25,2),sep=""))+factor(Kappa,labels=paste("Kappa=",c(0.25,2),sep="")),data=subset(d,Model %in% 5:6),groups=factor(Xi,labels=paste("Xi=",c(0.25,2),sep="")),auto.key=list(column=2,cex=0.75),scales=list(cex=0.5),pch=21)
dev.off()
