## Setup
library(parallel)
library(foreach)
library(doParallel)
library(lattice)
s      <- 500
sample <- c(500,1000)
model  <- 3:6
xi <- nu <- kappa <- c(0.25,2)
calc.stats <- function(data,index,results,s)
{
    foreach (n=1:s,.multicombine=T) %dopar%
        {
            rmse.temp <- sqrt(mean((data[[index]][[n]] - rowMeans(results[[n]]))^2))
            cov.temp <- mean(apply(results[[n]],1,function(x) quantile(x,0.025))<=data[[index]][[n]] & data[[index]][[n]]<=apply(results[[n]],1, function(x) quantile(x,0.975)))
            length.temp <- mean(-apply(results[[n]],1,function(x) quantile(x,0.025)) + apply(results[[n]],1, function(x) quantile(x,0.975)))
            return(list(rmse=rmse.temp,cov=cov.temp,length=length.temp))
        }
}
## Parallelization
no_cores <- detectCores() - 1
registerDoParallel(no_cores)
## Load data
data <- vector("list",32)
index <- 0
for (i in 1:length(model))
{
    for (j in 1:length(xi))
    {
        for (k in 1:length(nu))
        {
            for (l in 1:length(kappa))
            {
                index <- index+1
                for (m in 1:s)
                {
                    dgp <- paste(model[i],xi[j],nu[k],kappa[l],m,sep="_")
                    dgp <- readRDS(paste0("Data/DGP_500_",dgp))
                    test <- -0.1<=dgp$X & dgp$X<=0.1
                    data[[index]][[m]] <- dgp$cate[test,]
                }
            }
        }
    }
}
## BART-RDD
bart.rdd <- vector("list",32)
res.mat <- as.data.frame(matrix(0,32,5,dimnames=list(1:32,c("BART-RDD","Model","Xi","Nu","Kappa"))))
rmse.bart.rdd <- res.mat
cov.bart.rdd <- res.mat
length.bart.rdd <- res.mat
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
                print(paste(c("Model","Xi","Nu","Kappa"),dgp,sep=": "))
                file <- paste(model[j],xi[k],nu[l],kappa[m],sep="_")
                file <- paste0("Results/bart_rdd_500_",file)
                bart.rdd <- readRDS(paste0(file,".rds"))
                bart.rdd <- lapply(bart.rdd, function(x) x$pred)
                temp <- calc.stats(data,index,bart.rdd,s)
                temp <- do.call("rbind",temp)
                temp <- apply(temp,2,unlist)
                temp <- colMeans(temp)
                rmse.bart.rdd[index,] <- c(temp[1],dgp)
                cov.bart.rdd[index,] <- c(temp[2],dgp)
                length.bart.rdd[index,] <- c(temp[3],dgp)
            }
        }
    }
}
## BART1
res.mat <- as.data.frame(matrix(0,32,5,dimnames=list(1:32,c("BART1","Model","Xi","Nu","Kappa"))))
rmse.bart.1 <- res.mat
cov.bart.1 <- res.mat
length.bart.1 <- res.mat
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
                print(paste(c("Model","Xi","Nu","Kappa"),dgp,sep=": "))
                file <- paste(model[j],xi[k],nu[l],kappa[m],sep="_")
                file <- paste0("Results/bart1_500_",file)
                bart.1 <- readRDS(paste0(file,".rds"))
                bart.1 <- lapply(bart.1, function(x) x$pred)[1:s]
                temp <- calc.stats(data,index,bart.1,500)
                temp <- do.call("rbind",temp)
                temp <- apply(temp,2,unlist)
                temp <- colMeans(temp)
                rmse.bart.1[index,] <- c(temp[1],dgp)
                cov.bart.1[index,] <- c(temp[2],dgp)
                length.bart.1[index,] <- c(temp[3],dgp)
            }
        }
    }
}
## BCF
res.mat <- as.data.frame(matrix(0,32,5,dimnames=list(1:32,c("BCF","Model","Xi","Nu","Kappa"))))
rmse.bcf <- res.mat
cov.bcf <- res.mat
length.bcf <- res.mat
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
                print(paste(c("Model","Xi","Nu","Kappa"),dgp,sep=": "))
                file <- paste(model[j],xi[k],nu[l],kappa[m],sep="_")
                file <- paste0("Results/bcf_500_",file)
                bcf <- readRDS(paste0(file,".rds"))
                bcf <- lapply(bcf, function(x) x$pred)[1:s]
                temp <- calc.stats(data,index,bcf,500)
                temp <- do.call("rbind",temp)
                temp <- apply(temp,2,unlist)
                temp <- colMeans(temp)
                rmse.bcf[index,] <- c(temp[1],dgp)
                cov.bcf[index,] <- c(temp[2],dgp)
                length.bcf[index,] <- c(temp[3],dgp)
            }
        }
    }
}
## Merge results
rmse <- merge(rmse.bart.rdd,rmse.bart.1)
rmse <- merge(rmse,rmse.bcf)
cov <- merge(cov.bart.rdd,cov.bart.1)
cov <- merge(cov,cov.bcf)
length <- merge(length.bart.rdd,length.bart.1)
length <- merge(length,length.bcf)
## Plots
plot.rmse <- reshape(rmse,direction="long",varying=c("BART-RDD","BART1","BCF"),v.names="RMSE",idvar=c("Model","Xi","Nu","Kappa"),timevar="Method",times=c("BART-RDD","BART1","BCF"))
plot.cov <- reshape(cov,direction="long",varying=c("BART-RDD","BART1","BCF"),v.names="Coverage",idvar=c("Model","Xi","Nu","Kappa"),timevar="Method",times=c("BART-RDD","BART1","BCF"))
plot.length <- reshape(length,direction="long",varying=c("BART-RDD","BART1","BCF"),v.names="Length",idvar=c("Model","Xi","Nu","Kappa"),timevar="Method",times=c("BART-RDD","BART1","BCF"))
### RMSE
plot.rmse$Method <- as.factor(plot.rmse$Method)
png("Figures/rmse_cate.png")
xyplot(RMSE~Method|factor(Nu,labels=paste("Nu=",c(0.25,2),sep=""))+factor(Kappa,labels=paste("Kappa=",c(0.25,2),sep="")),data=plot.rmse,groups=factor(Xi,labels=paste("Xi=",c(0.25,2),sep="")),auto.key=list(column=2,cex=0.75),scales=list(cex=0.5),pch=21)
dev.off()
### Cov x Length
p <- merge(plot.cov,plot.length)
png("Figures/coverage_cate.png")
xyplot(Coverage~Length|factor(Nu,labels=paste("Nu=",c(0.25,2),sep=""))+factor(Kappa,labels=paste("Kappa=",c(0.25,2),sep="")),data=p,groups=Method,pch=21,auto.key=list(column=3,cex=0.75))
dev.off()
