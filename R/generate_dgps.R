## Setup
if (Sys.info()[1]=="Linux")
{
  setwd("~/Documents/Git/XBCF-RDD")
} else
{
  setwd("~/Git/XBCF-RDD/")
}
set.seed(0)
s <- 1000
n <- 500
## DGP functions
### 1) CGS
dgp.cgs <- function(n,p)
{
    ## Additional covariates are unnecessary
    w  <- matrix(rnorm(n*p),n,p)
    x  <- 2*rbeta(n,2,4)-1
    z  <- x>=0
    mu <- function(x)
    {
        0.48 + 1.27*x + 7.18*x^2 + 20.21*x^3 + 21.54*x^4 + 7.33*x^5
    }
    tau <- function(x)
    {
        0.04 - 0.43*x - 10.18*x^2 - 12.22*x^3 - 30.55*x^4 - 3.77*x^5
    }
    y  <- mu(x) + tau(x)*z + 0.1295*rt(n,3)
    ate <- 0.04
    tau.x <- tau(x)
    return(list(y=y,x=x,z=z,w=w,ate=ate,tau.x=tau.x))
}
### 2) FH
dgp.fh <- function(n,p=2)
{
    ## We consider only the case where the covariates affect Y but don't change across treatment states
    x <- 2*rbeta(n,2,4)-1
    z <- x >= 0
    w <- matrix(rnorm(n*p,sd=0.5),n,p)
    a <- 2/(1:p)^2
    mu <- function(x,w)
    {
        0.5*x + 0.25*x^2 + 0.4*w%*%a + 0.2*w^2%*%a
    }
    tau <- function(x)
    {
        1 - 0.25*x
    }
    y <- mu(x,w) + tau(x)*z + rnorm(n,sd=0.25)
    ate <- 1
    tau.x <- tau(x)
    return(list(y=y,x=x,z=z,w=w,ate=ate,tau.x=tau.x))
}
### 3) KR
dgp.kr <- function(n,p=200)
{
    sig.e <- 0.1295^2
    sig.z <- 0.1353^2
    v     <- (0.8*sqrt(6)*sig.e/pi)/1:p
    var   <- matrix(0, p+1, p+1)
    var[1,1] <- sig.e
    var[1,2:(p+1)] <- v
    var[2:(p+1),1] <- v
    var[2:(p+1),2:(p+1)] <- sig.z*diag(p)
    eZ <- mnormt::rmnorm(n,varcov=var)
    e <- eZ[,1]
    w <- eZ[,-1]
    a <- 2/(1:p)^2
    x  <- 2*rbeta(n,2,4)-1
    z  <- x>=0
    mu <- function(x,w,a)
    {
        0.36 + 0.96*x + 5.47*x^2 + 15.28*x^3 + 5.14*x^5 + as.vector(0.22*w%*%a)
    }
    tau <- function(x,w,a)
    {
        2 - 0.34*x - 8.31*x^2 - 6.86*x^3 - 0.83*x^5 + as.vector(0.06*w%*%a)
    }
    y  <- mu(x,w,a) + tau(x,w,a)*z + e
    ate <- mean(tau(0,w,a))
    tau.x <- tau(x,w,a)
    return(list(y=y,x=x,z=z,w=w,ate=ate,tau.x=tau.x))
}
### 4) W and X independent, interact in tau
dgp4.fun <- function(n,p)
{
    x <- 2*rbeta(n,2,4)-1
    z <- x >= 0
    w <- matrix(rnorm(n*p,0,0.5),n,p)
    mu <- function(x,w)
    {
        1/(1+exp(-5*x)) + rowMeans(w)
    }
    tau <- function(x,w)
    {
        1 + 0.5*x*rowMeans(cos(w)) + sin(rowSums(w))
    }
    y <- mu(x,w) + z*tau(x,w) + rnorm(n)
    ate <- mean(tau(0,w))
    tau.x <- tau(x,w)
    return(list(y=y,x=x,z=z,w=w,ate=ate,tau.x=tau.x))
}
### 5) W and X dependent, don't interact in tau
dgp5.fun <- function(n,p)
{
    x <- 2*rbeta(n,2,4)-1
    z <- x >= 0
    w2 <- 0.5*x + rnorm(n,0,0.5)
    w3 <- x^2 - 0.1*x + rnorm(n,0,0.5)
    w4 <- sin(x) + rnorm(n,0,0.5)
    w5 <- x^4 + 0.2*x^2 - x + rnorm(n,0,0.5)
    if (p!=0)
    {
        w1 <- matrix(rnorm(n*p,0,0.5),n,p)
        w <- cbind(w1,w2)
    } else {
        w <- w2
    }
    w <- cbind(w,w3)
    w <- cbind(w,w4)
    w <- cbind(w,w5)
    mu <- function(x,w)
    {
        1/(1+exp(-5*x)) + rowMeans(w)
    }
    tau <- function(x,w)
    {
        1 + x + sin(rowSums(w))
    }
    y <- mu(x,w) + z*tau(x,w) + rnorm(n)
    ate <- mean(tau(0,w))
    tau.x <- tau(x,w)
    return(list(y=y,x=x,z=z,w=w,ate=ate,tau.x=tau.x))
}
### 6) W and X dependent, interact in tau
dgp6.fun <- function(n,p)
{
    x <- 2*rbeta(n,2,4)-1
    z <- x >= 0
    w2 <- 0.5*x + rnorm(n,0,0.5)
    w3 <- x^2 - 0.1*x + rnorm(n,0,0.5)
    w4 <- sin(x) + rnorm(n,0,0.5)
    w5 <- x^4 + 0.2*x^2 - x + rnorm(n,0,0.5)
    if (p!=0)
    {
        w1 <- matrix(rnorm(n*p,0,0.5),n,p)
        w <- cbind(w1,w2)
    } else {
        w <- w2
    }
    w <- cbind(w,w3)
    w <- cbind(w,w4)
    w <- cbind(w,w5)
    mu <- function(x,w)
    {
        1/(1+exp(-5*x)) + rowMeans(w)
    }
    tau <- function(x,w)
    {
        1 + 0.5*x*rowMeans(cos(w)) + sin(rowSums(w)) + 0.7*(w2>=0) - 0.5*(w3<0 & w4>=0.1)
    }
    y <- mu(x,w) + z*tau(x,w) + rnorm(n)
    ate <- mean(tau(0,w))
    tau.x <- tau(x,w)
    return(list(y=y,x=x,z=z,w=w,ate=ate,tau.x=tau.x))
}
## DGP1
for (P in c(4,6,10))
{
    dgp <- vector("list",s)
    for (i in 1:s) dgp[[i]] <- dgp.cgs(n,P)
    saveRDS(dgp,paste0("Data/DGP1_",P,".rds"))
}
## DGP2
for (P in c(4,6,10))
{
    dgp <- vector("list",s)
    for (i in 1:s) dgp[[i]] <- dgp.fh(n,P)
    saveRDS(dgp,paste0("Data/DGP2_",P,".rds"))
}
## DGP3
for (P in c(4,6,10))
{
    dgp <- vector("list",s)
    for (i in 1:s) dgp[[i]] <- dgp.kr(n,P)
    saveRDS(dgp,paste0("Data/DGP3_",P,".rds"))
}
## DGP4
for (P in c(4,6,10))
{
    dgp <- vector("list",s)
    for (i in 1:s) dgp[[i]] <- dgp4.fun(n,P)
    saveRDS(dgp,paste0("Data/DGP4_",P,".rds"))
}
## DGP5
for (P in c(0,2,6))
{
    dgp <- vector("list",s)
    for (i in 1:s) dgp[[i]] <- dgp5.fun(n,P)
    saveRDS(dgp,paste0("Data/DGP5_",P+4,".rds"))
}
## DGP6
for (P in c(0,2,6))
{
    dgp <- vector("list",s)
    for (i in 1:s) dgp[[i]] <- dgp6.fun(n,P)
    saveRDS(dgp,paste0("Data/DGP6_",P+4,".rds"))
}
