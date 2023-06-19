## Setup
library(rdrobust)
library(parallel)
library(foreach)
library(doParallel)
n <- 500
s <- 5000
### Parallelization
no_cores <- detectCores() - 1
registerDoParallel(no_cores)
###
fit1 <- function(s)
{
    foreach(i=1:s,.multicombine=T,.export="n") %dopar%
        {
            print(paste0("Iteration ",i," for fit 1"))
            x <- 2*rbeta(n,2,4)-1
            z <- x>=0
            mu <- function(x) 1/(1+exp(-5*x))
            tau <- function(x) 2 + x
            yl <- mu(x) + tau(x)*z + rnorm(n)
            rdrobust(yl,x,c=0)
        }
}
fit2 <- function(s)
{
    foreach(i=1:s,.multicombine=T,.export="n") %dopar%
        {
            print(paste0("Iteration ",i," for fit 2"))
            x <- 2*rbeta(n,2,4)-1
            z <- x>=0
            mu <- function(x) 1/(1+exp(-5*x))
            tau <- function(x) 2 + x
            y <- mu(x) + tau(x)*z + rnorm(n)
            x1 <- round(x,digits=3)
            rdrobust(y,x1,c=0)
        }
}
fit3 <- function(s)
{
    foreach(i=1:s,.multicombine=T,.export="n") %dopar%
        {
            print(paste0("Iteration ",i," for fit 3"))
            x <- 2*rbeta(n,2,4)-1
            z <- x>=0
            mu <- function(x) 1/(1+exp(-5*x))
            tau <- function(x) 2 + x
            y <- mu(x) + tau(x)*z + rnorm(n)
            x1 <- round(x,digits=2)
            rdrobust(y,x1,c=0)
        }
}
fit4 <- function(s)
{
    foreach(i=1:s,.multicombine=T,.export="n") %dopar%
        {
            print(paste0("Iteration ",i," for fit 4"))
            x <- 2*rbeta(n,2,4)-1
            z <- x>=0
            mu <- function(x) 1/(1+exp(-5*x))
            tau <- function(x) 2 + x
            y <- mu(x) + tau(x)*z + rnorm(n)
            x1 <- round(x,digits=1)
            rdrobust(y,x1,c=0)
        }
}
fit5 <- function(s)
{
    set.seed(0)
    foreach(i=1:s,.multicombine=T,.export="n") %dopar%
        {
            print(paste0("Iteration ",i," for fit 5"))
            x <- 2*rbeta(n,2,4)-1
            z <- x>=0
            mu <- function(x) 1/(1+exp(-5*x))
            tau <- function(x) 2 + x
            yl <- mu(x) + tau(x)*z + rnorm(n)
            y <- yl>=0
            rdrobust(as.numeric(y),x,c=0)
        }
}
cct1 <- fit1(s)
cct2 <- fit2(s)
cct3 <- fit3(s)
cct4 <- fit4(s)
cct5 <- fit5(s)
stopImplicitCluster()
###
cct1 <- sapply(cct1, function(i) c(i$Estimate[1],i$Estimate[1]-1.96*i$Estimate[3],i$Estimate[1]+1.96*i$Estimate[3]))
cct2 <- sapply(cct2, function(i) c(i$Estimate[1],i$Estimate[1]-1.96*i$Estimate[3],i$Estimate[1]+1.96*i$Estimate[3]))
cct3 <- sapply(cct3, function(i) c(i$Estimate[1],i$Estimate[1]-1.96*i$Estimate[3],i$Estimate[1]+1.96*i$Estimate[3]))
cct4 <- sapply(cct4, function(i) c(i$Estimate[1],i$Estimate[1]-1.96*i$Estimate[3],i$Estimate[1]+1.96*i$Estimate[3]))
cct5 <- sapply(cct5, function(i) c(i$Estimate[1],i$Estimate[1]-1.96*i$Estimate[3],i$Estimate[1]+1.96*i$Estimate[3]))
cct1 <- t(cct1)
cct2 <- t(cct2)
cct3 <- t(cct3)
cct4 <- t(cct4)
cct5 <- t(cct5)
###
cols <- rainbow(5,s=0.5)
p <- data.frame(Fit=c(1,1,2,2,3,3,4,4,5,5),
                V=c("Squared Bias","Variance"),
                Value=c(mean(cct1[,1]-2)^2,
                        var(cct1[,1]),
                        mean(cct2[,1]-2)^2,
                        var(cct2[,1]),
                        mean(cct3[,1]-2)^2,
                        var(cct3[,1]),
                        mean(cct4[,1]-2)^2,
                        var(cct4[,1]),
                        mean(cct5[,1]-2)^2,
                        var(cct5[,1])))
barplot(p$Value~p$Fit+p$V,beside=T,xlab="",ylab="",col=cols)
legend(x="topright",ncol=3,legend=1:5,fill=cols)
