## setwd("~/Git/XBCF-RDD/")
set.seed(0)
## General settings
s <- 500
n <- c(500,2000)
p <- 7
m <- 1:4
xi <- c(0.05,0.2)
nu <- kappa <- c(0.25,2)
x.dist <- 1:4
vars <- function(n,p,x.dist)
{
    if (x.dist==1) X <- 2*rbeta(n,2,4)-1
    if (x.dist==2) X <- 2*rbeta(n,2,4)-0.5
    if (x.dist==3) X <- rnorm(n,0,0.5)
    if (x.dist==4) X <- runif(n,-1,1)
    W1 <- matrix(rnorm(n*2,0,0.25),n,2)
    Wp <- lapply(3:(p-3),function(i) rnorm(n,(i-1)*X/i,1))
    Wp <- do.call("cbind",Wp)
    W2 <- cbind(rbinom(n,1,0.7),rbinom(n,1,0.3),as.numeric(X>=0.3))
    W <- cbind(W1,Wp,W2)
    Z <- X>=0
    return(list(X=X,W=W,Z=Z))
}
mu <- function(X,W,n,p,m)
{
    a <- 2/1:p
    if (m==1)
    {
        mu <- 0.1*X - 0.2*X^2 + 0.5*X^3 + W%*%a
    } else if (m==2)
    {
        mu <- 0.1*X - 0.2*X^2 + 0.5*X^3 + W%*%a + W[,p]*X - 0.5*W[,p-1]*X + 0.3*W[,p-2]
    } else if (m==3)
    {
        mu <- exp(X) + sqrt(abs(W))%*%a
    } else
    {
        mu <- exp(X) + W[,p-2]*(sqrt(abs(W))%*%a) + W[,p]*X - 0.5*W[,p-1]*X + 0.3*W[,p-2]
    }
    return(mu)
}
tau <- function(X,W,n,p,m)
{
    b <- 1/1:p
    if (m==1)
    {
        tau <- 0.7*X + 0.4*X^2 - 0.1*X^3 + W%*%b
    } else if (m==2)
    {
        tau <- 0.7*X + 0.4*X^2 - 0.1*X^3 + W%*%b + 0.5*W[,p]*X - W[,p-1]*X
    } else if (m==3)
    {
        tau <- sin(X) + sqrt(abs(W))%*%b
    } else
    {
        tau <- sin(X) + W[,p-2]*(sqrt(abs(W))%*%b) + 0.5*W[,p]*X - W[,p-1]*X
    }
    return(tau)
}
sig.mu <- function(W,n,p,m) sd(mu(0,W,n,p,m))
sig.tau <- function(W,n,p,m) sd(tau(0,W,n,p,m))
tau.bar <- function(W,n,p,m) mean(tau(0,W,n,p,m))
mu.std <- function(X,W,n,p,m) mu(X,W,n,p,m)/sig.mu(W,n,p,m)
tau.std <- function(X,W,n,p,m,xi,nu)
{
    return(xi + nu*(tau(X,W,n,p,m)-tau.bar(W,n,p,m))/sig.tau(W,n,p,m))
}
out <- function(n,Z,mu,tau,kappa) mu + tau*Z + kappa*rnorm(n)
## Function to generate data
dgp <- function(n,p,m,xi,nu,kappa,x.dist)
{
    v <- vars(n,p,x.dist)
    X <- v$X
    W <- v$W
    Z <- v$Z
    M <- mu.std(X,W,n,p,m)
    T <- tau.std(X,W,n,p,m,xi,nu)
    Y <- out(n,Z,M,T,kappa)
    cate <- tau.std(0,W,n,p,m,xi,nu)
    ate <- mean(cate)
    return(list(Y=Y,X=X,W=W,Z=Z,cate=cate,ate=ate))
}
## Saving data
### name: DGP_n_m_xi_nu_kappa_s
for (i in 1:s)
{
    for (j in n)
    {
        for (k in m)
        {
            for (l in xi)
            {
                for (o in nu)
                {
                    for (q in kappa)
                    {
                        for (r in x.dist)
                            {
                                data <- dgp(j,p,k,l,o,q,r)
                                saveRDS(data,paste0("Data/DGP_",j,"_",k,"_",l,"_",o,"_",q,"_",r,"_",i,".rds"))
                            }
                    }
                }
            }
        }
    }
}
