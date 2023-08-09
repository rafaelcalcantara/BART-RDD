## Setup
library(lattice)
s      <- 1000
sample <- c(500,1000)
model  <- 3:6
xi <- nu <- kappa <- c(0.25,2)
i <- 1
## BART-RDD
res.mat <- as.data.frame(matrix(0,32,5,dimnames=list(1:32,c("BART-RDD","Model","Xi","Nu","Kappa"))))
bart.rdd.dt <- res.mat
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
                bart.rdd <- sapply(bart.rdd, function(x) as.numeric(x$dt,units="secs"))
                bart.rdd.dt[index,] <- c(mean(bart.rdd),dgp)
            }
        }
    }
}
## CKT
res.mat <- as.data.frame(matrix(0,32,5,dimnames=list(1:32,c("CKT","Model","Xi","Nu","Kappa"))))
ckt.dt <- res.mat
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
                ckt <- sapply(ckt, function(x) as.numeric(x$dt,units="secs"))
                ckt.dt[index,] <- c(mean(ckt),dgp)
            }
        }
    }
}
## KR
res.mat <- as.data.frame(matrix(0,32,5,dimnames=list(1:32,c("KR","Model","Xi","Nu","Kappa"))))
kr.dt <- res.mat
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
                kr <- sapply(kr, function(x) as.numeric(x$dt,units="secs"))
                kr.dt[index,] <- c(mean(kr),dgp)
            }
        }
    }
}
## BART1
res.mat <- as.data.frame(matrix(0,32,5,dimnames=list(1:32,c("BART1","Model","Xi","Nu","Kappa"))))
bart1.dt <- res.mat
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
                bart1 <- readRDS(paste0(file,".rds"))
                bart1 <- sapply(bart1, function(x) as.numeric(x$dt,units="secs"))
                bart1.dt[index,] <- c(mean(bart1),dgp)
            }
        }
    }
}
## BART2
res.mat <- as.data.frame(matrix(0,32,5,dimnames=list(1:32,c("BART2","Model","Xi","Nu","Kappa"))))
bart2.dt <- res.mat
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
                bart2 <- readRDS(paste0(file,".rds"))
                bart2 <- sapply(bart2, function(x) as.numeric(x$dt,units="secs"))
                bart2.dt[index,] <- c(mean(bart2),dgp)
            }
        }
    }
}
## BCF
res.mat <- as.data.frame(matrix(0,32,5,dimnames=list(1:32,c("BCF","Model","Xi","Nu","Kappa"))))
bcf.dt <- res.mat
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
                bcf <- sapply(bcf, function(x) as.numeric(x$dt,units="secs"))
                bcf.dt[index,] <- c(mean(bcf),dgp)
            }
        }
    }
}
## Merge results
load("Tables/cgs.RData")
dt <- merge(bart.rdd.dt,bart1.dt)
dt <- merge(dt,bart2.dt)
dt <- merge(dt,bcf.dt)
dt <- merge(dt,ckt.dt)
dt <- merge(dt,kr.dt)
## dt <- merge(dt,cgs.dt)
## Plots
plot.dt <- reshape(dt,direction="long",varying=c("BART-RDD","BART1","BART2","BCF","CKT","KR"),v.names="Runtime",timevar="Method",times=c("BART-RDD","BART1","BART2","BCF","CKT","KR"))
### RMSE
plot.dt$Method <- as.factor(plot.dt$Method)
## png("Figures/runtime.png")
boxplot(Runtime~Method,data=plot.dt,axes=F)
axis(1,at=1:6,labels=c("BART-RDD","BART1","BART2","BCF","CKT","KR"),cex.axis=0.75)
axis(2,cex.axis=0.75)
## dev.off()
