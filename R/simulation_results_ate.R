## Setup
s      <- 500
sample <- c(500,1000)
model  <- 3:6
xi <- nu <- kappa <- c(0.25,2)
## BART-RDD
bart.rdd <- vector("list",32)
res.mat <- as.data.frame(matrix(0,32,5,dimnames=list(1:32,c("BART-RDD","Model","Xi","Nu","Kappa"))))
rmse.bart.rdd <- res.mat
cov.bart.rdd <- res.mat
length.bart.rdd <- res.mat
for (i in 1:length(500))
{
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
                    bart.rdd[[index]] <- readRDS(paste0(file,".rds"))
                    bart.rdd[[index]] <- t(sapply(bart.rdd[[index]],function(x) c(mean(colMeans(x)),quantile(colMeans(x),c(0.025,0.975)))))
                    rmse.bart.rdd[index,] <- c(sqrt(mean((bart.rdd[[index]][,1]-xi[k])^2)),dgp)
                    cov.bart.rdd[index,] <- c(mean(bart.rdd[[index]][,2] <= xi[k] & xi[k] <= bart.rdd[[index]][,3]),dgp)
                    length.bart.rdd[index,] <- c(mean(bart.rdd[[index]][,3] - bart.rdd[[index]][,2]),dgp)
                }
            }
        }
    }
}
## CKT
ckt <- vector("list",32)
res.mat <- as.data.frame(matrix(0,32,5,dimnames=list(1:32,c("CKT","Model","Xi","Nu","Kappa"))))
rmse.ckt <- res.mat
cov.ckt <- res.mat
length.ckt <- res.mat
for (i in 1:length(500))
{
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
                    ckt[[index]] <- do.call("rbind",ckt[[index]])
                    rmse.ckt[index,] <- c(sqrt(mean((ckt[[index]][,1]-xi[k])^2)),dgp)
                    cov.ckt[index,] <- c(mean(ckt[[index]][,2] <= xi[k] & xi[k] <= ckt[[index]][,3]),dgp)
                    length.ckt[index,] <- c(mean(ckt[[index]][,3] - ckt[[index]][,2]),dgp)
                }
            }
        }
    }
}
## KR
kr <- vector("list",32)
res.mat <- as.data.frame(matrix(0,32,5,dimnames=list(1:32,c("KR","Model","Xi","Nu","Kappa"))))
rmse.kr <- res.mat
cov.kr <- res.mat
length.kr <- res.mat
for (i in 1:length(500))
{
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
                    kr[[index]] <- readRDS(paste0(file,".rds"))
                    kr[[index]] <- do.call("rbind",kr[[index]])
                    rmse.kr[index,] <- c(sqrt(mean((kr[[index]][,1]-xi[k])^2)),dgp)
                    cov.kr[index,] <- c(mean(kr[[index]][,2] <= xi[k] & xi[k] <= kr[[index]][,3]),dgp)
                    length.kr[index,] <- c(mean(kr[[index]][,3] - kr[[index]][,2]),dgp)
                }
            }
        }
    }
}
## BART1
bart.1 <- vector("list",32)
res.mat <- as.data.frame(matrix(0,32,5,dimnames=list(1:32,c("BART1","Model","Xi","Nu","Kappa"))))
rmse.bart.1 <- res.mat
cov.bart.1 <- res.mat
length.bart.1 <- res.mat
for (i in 1:length(500))
{
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
                    dt <- sapply(bart.1,function(x) x$dt)
                    bart.1 <- lapply(bart.1, function(x) x$pred)
                    bart.1 <- t(sapply(bart.1,function(x) c(mean(colMeans(x)),quantile(colMeans(x),c(0.025,0.975)))))
                    rmse.bart.1[index,] <- c(sqrt(mean((bart.1[,1]-xi[k])^2)),dgp)
                    cov.bart.1[index,] <- c(mean(bart.1[,2] <= xi[k] & xi[k] <= bart.1[,3]),dgp)
                    length.bart.1[index,] <- c(mean(bart.1[,3] - bart.1[,2]),dgp)
                }
            }
        }
    }
}
## BART2
bart.2 <- vector("list",32)
res.mat <- as.data.frame(matrix(0,32,5,dimnames=list(1:32,c("BART2","Model","Xi","Nu","Kappa"))))
rmse.bart.2 <- res.mat
cov.bart.2 <- res.mat
length.bart.2 <- res.mat
for (i in 1:length(500))
{
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
                    dt <- sapply(bart.2,function(x) x$dt)
                    bart.2 <- lapply(bart.2, function(x) x$pred)
                    bart.2 <- t(sapply(bart.2,function(x) c(mean(x),quantile(x,c(0.025,0.975)))))
                    rmse.bart.2[index,] <- c(sqrt(mean((bart.2[,1]-xi[k])^2)),dgp)
                    cov.bart.2[index,] <- c(mean(bart.2[,2] <= xi[k] & xi[k] <= bart.2[,3]),dgp)
                    length.bart.2[index,] <- c(mean(bart.2[,3] - bart.2[,2]),dgp)
                }
            }
        }
    }
}
## BCF
bcf <- vector("list",32)
res.mat <- as.data.frame(matrix(0,32,5,dimnames=list(1:32,c("BCF","Model","Xi","Nu","Kappa"))))
rmse.bcf <- res.mat
cov.bcf <- res.mat
length.bcf <- res.mat
for (i in 1:length(500))
{
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
                    dt <- sapply(bcf,function(x) x$dt)
                    bcf <- lapply(bcf, function(x) x$pred)
                    bcf <- t(sapply(bcf,function(x) c(mean(colMeans(x)),quantile(colMeans(x),c(0.025,0.975)))))
                    rmse.bcf[index,] <- c(sqrt(mean((bcf[,1]-xi[k])^2)),dgp)
                    cov.bcf[index,] <- c(mean(bcf[,2] <= xi[k] & xi[k] <= bcf[,3]),dgp)
                    length.bcf[index,] <- c(mean(bcf[,3] - bcf[,2]),dgp)
                }
            }
        }
    }
}
## Merge results
rmse <- merge(rmse.bart.rdd,rmse.ckt)
## rmse <- merge(rmse,rmse.kr)
rmse <- merge(rmse,rmse.bart.1)
rmse <- merge(rmse,rmse.bart.2)
rmse <- merge(rmse,rmse.bcf)
cov <- merge(cov.bart.rdd,cov.ckt)
cov <- merge(cov,cov.kr)
cov <- merge(cov,cov.bart.1)
cov <- merge(cov,cov.bart.2)
cov <- merge(cov,cov.bcf)
length <- merge(length.bart.rdd,length.ckt)
length <- merge(length,length.kr)
length <- merge(length,length.bart.1)
length <- merge(length,length.bart.2)
length <- merge(length,length.bcf)
## Plots
plot.rmse <- reshape(rmse,direction="long",varying=c("BART-RDD","CKT","KR","BART1","BART2","BCF"),v.names="RMSE",idvar=c("Model","Xi","Nu","Kappa"),timevar="Method",times=c("BART-RDD","CKT","KR","BART1","BART2","BCF"))
plot.cov <- reshape(cov,direction="long",varying=c("BART-RDD","CKT","KR","BART1","BART2","BCF"),v.names="Coverage",idvar=c("Model","Xi","Nu","Kappa"),timevar="Method",times=c("BART-RDD","CKT","KR","BART1","BART2","BCF"))
plot.length <- reshape(length,direction="long",varying=c("BART-RDD","CKT","KR","BART1","BART2","BCF"),v.names="Length",idvar=c("Model","Xi","Nu","Kappa"),timevar="Method",times=c("BART-RDD","CKT","KR","BART1","BART2","BCF"))
### RMSE
plot.rmse$Method <- as.factor(plot.rmse$Method)
coplot(RMSE~Method|as.factor(Nu)*as.factor(Kappa),data=plot.rmse,panel = function(x, y, ...){boxplot(y ~ x, add=TRUE,names=F)},xlim=c(0,7),number=1,overlap=0.1)
### Cov x Length
p <- merge(plot.cov,plot.length)
cols <- rainbow(7)
coplot(Coverage~Length|as.factor(Nu)*as.factor(Kappa),data=p,
       bg=cols,pch=21,show.given=T)
legend("topright",legend=c("BART-RDD","CKT","KR","BART1","BART2","BCF"),
       pt.bg=cols,pch=21,cex=0.75)
#####
bart.test <- bart.rdd[rmse.bart.rdd$Xi==0.25 & rmse.bart.rdd$Nu==2 & rmse.bart.rdd$Kappa==0.25]
bart.test <- do.call("rbind",bart.test)
bart.test <- bart.test[bart.test[,2]>=0.25|0.25>=bart.test[,3],]
mean(apply(bart.test,1,function(x) max(x[2]-0.25,0.25-x[3])))
bart.test <- bart.rdd[rmse.bart.rdd$Xi==2 & rmse.bart.rdd$Nu==2 & rmse.bart.rdd$Kappa==0.25]
bart.test <- do.call("rbind",bart.test)
bart.test <- bart.test[bart.test[,2]>=2|2>=bart.test[,3],]
mean(apply(bart.test,1,function(x) max(x[2]-2,2-x[3])))
