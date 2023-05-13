## Setup
s <- 1000
n <- 500
p <- 10
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
        0.36 + 0.96*x + 5.47*x^2 + 15.28*x^3 + 5.14*x^5 + as.vector(w%*%a)
    }
    tau <- function(x,w,a)
    {
        0.02 - 0.34*x - 8.31*x^2 - 6.86*x^3 - 0.83*x^5 + as.vector(0.06*w%*%a)
    }
    y  <- mu(x,w,a) + tau(x,w,a)*z + e
    ate <- mean(tau(0,w,a))
    tau.x <- tau(x,w,a)
    return(list(y=y,x=x,z=z,w=w,ate=ate,tau.x=tau.x))
}
## DGP1a
dgp1a <- vector("list",s)
for (i in 1:s) dgp1a[[i]] <- dgp.cgs(n,p)
saveRDS(dgp1a,"Data/DGP1a.rds")
## DGP1b
dgp1b <- vector("list",s)
for (i in 1:s) dgp1b[[i]] <- dgp1a[[i]]
saveRDS(dgp1b,"Data/DGP1b.rds")
## DGP2
dgp2 <- vector("list",s)
for (i in 1:s) dgp2[[i]] <- dgp.fh(n,p)
saveRDS(dgp2,"Data/DGP2.rds")
## DGP3
dgp3 <- vector("list",s)
for (i in 1:s) dgp3[[i]] <- dgp.kr(n,p)
saveRDS(dgp3,"Data/DGP3.rds")
