## Setup
set.seed(0)
library(parallel)
library(foreach)
library(doParallel)
devtools::install_github("JingyuHe/XBART@XBCF-RDD")
library(XBART)
library(rdrobust)
## Parallelization
no_cores <- detectCores() - 1
registerDoParallel(no_cores)
##
s             <- 500
n             <- 500
c             <- 0
### DGP
sig1 <- 0.04
sig2 <- 0.03
sig3 <- 0.02
u1 <- u2 <- u3 <- runif(n,0,1)
w1 <- w2 <- w3 <- vector("list",s)
for (i in 1:s)
{
    W1 <- runif(n,u1,u1+1)
    W2 <- runif(n,0,0.5)
    W3 <- rbinom(n,2,u1)+1
    W4 <- rbinom(n,1,0.6)+1
    w1[[i]] <- cbind(W1,W2,W3,W4)
}
for (i in 1:s)
{
    W1 <- runif(n,u2,u2+1)
    W2 <- runif(n,0,0.5)
    W3 <- rbinom(n,2,u2)+1
    W4 <- rbinom(n,1,0.6)+1
    w2[[i]] <- cbind(W1,W2,W3,W4)
}
for (i in 1:s)
{
    W1 <- runif(n,u3,u3+1)
    W2 <- runif(n,0,0.5)
    W3 <- rbinom(n,2,u3)+1
    W4 <- rbinom(n,1,0.6)+1
    w3[[i]] <- cbind(W1,W2,W3,W4)
}
x1 <- matrix(2*rbeta(n*s,2,4)-u1-0.2,n,s)
x2 <- matrix(2*rbeta(n*s,2,4)-u2-0.2,n,s)
x3 <- matrix(2*rbeta(n*s,2,4)-u3-0.2,n,s)
z1 <- apply(x1,2,function(i) as.numeric(i>=c))
z2 <- apply(x2,2,function(i) as.numeric(i>=c))
z3 <- apply(x3,2,function(i) as.numeric(i>=c))
e1 <- matrix(rnorm(n*s,0,sig1),n,s)
e2 <- matrix(rnorm(n*s,0,sig2),n,s)
e3 <- matrix(rnorm(n*s,0,sig3),n,s)
mu <- function(x,w,m)
{
    w34 <- (-1)^(w[,4])*w[,3]
    w34 <- 0.1 + (w34-min(w34))/(max(w34)-min(w34))
    w34 <- round(w34,digits=1)
    if(m==1) out <- 0.01*sin((w[,1]+w[,2])*pi) + 1 + 0.1*x - 0.1*x^2 + 0.1*x^3
    if(m==2) out <- 0.01*w[,4]*sin((w[,1]+w[,2])*pi) + 1 + 0.1*w[,4] + 0.1*w[,4]*x - 0.1*w[,4]*x^2 + 0.1*w[,4]*x^3
    if(m==3) out <- 0.01*w34*sin((w[,1]+w[,2])*pi) + 1 + 0.1*w34 + 0.1*w34*x - 0.1*w34*x^2 + 0.1*w34*x^3
    return(out)
}
tau <- function(x,w,m)
{
    w34 <- (-1)^(w[,4])*w[,3]
    w34 <- 0.1 + (w34-min(w34))/(max(w34)-min(w34))
    w34 <- round(w34,digits=1)
    if(m==1) out <- 0.01*cos((w[,1]+w[,2])*pi) + 0.03 - 0.2*x + 0.1*x^2 - 0.01*x^3
    if(m==2) out <- 0.01*w[,4]*cos((w[,1]+w[,2])*pi) + 0.03*w[,4] - 0.2*w[,4]*x + 0.1*w[,4]*x^2 - 0.01*w[,4]*x^3
    if(m==3) out <- 0.01*w34*cos((w[,1]+w[,2])*pi) + 0.03*w34 - 0.2*w34*x + 0.1*w34*x^2 - 0.01*w34*x^3
    return(out)
}
y1 <- matrix(0,n,s)
y2 <- matrix(0,n,s)
y3 <- matrix(0,n,s)    
for (i in 1:s)
{
    y1[,i] <- mu(x1[,i],w1[[i]],1) + tau(x1[,i],w1[[i]],1)*z1[,i] + e1[,i]
    y2[,i] <- mu(x2[,i],w2[[i]],2) + tau(x2[,i],w2[[i]],2)*z2[,i] + e2[,i]
    y3[,i] <- mu(x3[,i],w3[[i]],3) + tau(x3[,i],w3[[i]],3)*z3[,i] + e3[,i]
}
## Plotting data and functions
X1 <- x1[,1]
X2 <- x2[,1]
X3 <- x3[,1]
Y1 <- y1[,1]
Y2 <- y2[,1]
Y3 <- y3[,1]
Z1 <- z1[,1]
Z2 <- z2[,1]
Z3 <- z3[,1]
m1 <- mu(X1,w1[[1]],1)
m2 <- mu(X2,w2[[1]],2)
m3 <- mu(X3,w3[[1]],3)
t1 <- tau(X1,w1[[1]],1)
t2 <- tau(X2,w2[[1]],2)
t3 <- tau(X3,w3[[1]],3)
w34 <- (-1)^(w3[[1]][,4])*w3[[1]][,3]
w34 <- 0.1 + (w34-min(w34))/(max(w34)-min(w34))
w34 <- as.factor(w34)
###
hist(X1,main="",freq=F,xlab="X")
###
par(mfrow=c(2,2))
hist(w1[[1]][,1],main="",freq=F,xlab=expression(W[1]))
hist(w1[[1]][,2],main="",freq=F,xlab=expression(W[2]))
barplot(table(w1[[1]][,3]),main="",xlab=expression(W[3]))
barplot(table(w1[[1]][,4]),main="",xlab=expression(W[4]))
###
par(mfrow=c(2,3))
plot(X1,t1,xlab="X",ylab=expression(tau),bty="n",pch=20)
abline(v=c,lty=2)
plot(X2,t2,xlab="X",ylab=expression(tau),bty="n",pch=20,col=w2[[1]][,4])
abline(v=c,lty=2)
plot(X3,t3,xlab="X",ylab=expression(tau),bty="n",pch=20,col=w34)
abline(v=c,lty=2)
###
plot(X1,m1+t1*Z1,xlab="X",ylab="E(Y)",bty="n",pch=20)
abline(v=c,lty=2)
plot(X2,m2+t2*Z2,xlab="X",ylab="E(Y)",bty="n",pch=20,col=w2[[1]][,4])
abline(v=c,lty=2)
plot(X3,m3+t3*Z3,xlab="X",ylab="E(Y)",bty="n",pch=20,col=w34)
abline(v=c,lty=2)
###
par(mfrow=c(1,3))
plot(X1,Y1,xlab="X",ylab="Y",bty="n",pch=20,col=Z1+1)
abline(v=c,lty=2)
plot(X2,Y2,xlab="X",ylab="Y",bty="n",pch=20,col=Z2+1)
abline(v=c,lty=2)
plot(X3,Y3,xlab="X",ylab="Y",bty="n",pch=20,col=Z3+1)
abline(v=c,lty=2)
###
y <- y + e
plot(x[,1],y[,1],bty="n",pch=20,col=z[,1]+1,xlab="X",ylab="Y")
#############################
## ## setwd("~/Git/XBCF-RDD/")
## set.seed(0)
## ## General settings
## s <- 500
## n <- c(500,2000)
## p <- 7
## m <- 1:4
## xi <- c(0.05,0.2)
## nu <- kappa <- c(0.5,2)
## x.dist <- 1:4
## vars <- function(n,p,x.dist)
## {
##     if (x.dist==1) X <- 2*rbeta(n,2,4)-1
##     if (x.dist==2) X <- 2*rbeta(n,2,4)-0.5
##     if (x.dist==3) X <- rnorm(n,0,0.5)
##     if (x.dist==4) X <- runif(n,-1,1)
##     W1 <- matrix(rnorm(n*2,0,0.25),n,2)
##     Wp <- lapply(3:(p-3),function(i) rnorm(n,(i-1)*X/i,1))
##     Wp <- do.call("cbind",Wp)
##     W2 <- cbind(rbinom(n,1,0.7),rbinom(n,1,0.3),as.numeric(X>=0.3))
##     W <- cbind(W1,Wp,W2)
##     Z <- X>=0
##     return(list(X=X,W=W,Z=Z))
## }
## mu <- function(X,W,n,p,m)
## {
##     a <- 2/1:p
##     if (m==1)
##     {
##         mu <- 0.1*X - 0.2*X^2 + 0.5*X^3 + W%*%a
##     } else if (m==2)
##     {
##         mu <- 0.1*X - 0.2*X^2 + 0.5*X^3 + W%*%a + W[,p]*X - 0.5*W[,p-1]*X + 0.3*W[,p-2]
##     } else if (m==3)
##     {
##         mu <- exp(X) + sqrt(abs(W))%*%a
##     } else
##     {
##         mu <- exp(X) + W[,p-2]*(sqrt(abs(W))%*%a) + W[,p]*X - 0.5*W[,p-1]*X + 0.3*W[,p-2]
##     }
##     return(mu)
## }
## tau <- function(X,W,n,p,m)
## {
##     b <- 1/1:p
##     if (m==1)
##     {
##         tau <- 0.7*X + 0.4*X^2 - 0.1*X^3 + W%*%b
##     } else if (m==2)
##     {
##         tau <- 0.7*X + 0.4*X^2 - 0.1*X^3 + W%*%b + 0.5*W[,p]*X - W[,p-1]*X
##     } else if (m==3)
##     {
##         tau <- sin(X) + sqrt(abs(W))%*%b
##     } else
##     {
##         tau <- sin(X) + W[,p-2]*(sqrt(abs(W))%*%b) + 0.5*W[,p]*X - W[,p-1]*X
##     }
##     return(tau)
## }
## sig.mu <- function(X,W,n,p,m) sd(mu(X,W,n,p,m))
## sig.tau <- function(X,W,n,p,m) sd(tau(X,W,n,p,m))
## tau.bar <- function(X,W,n,p,m) mean(tau(X,W,n,p,m))
## mu.std <- function(X,W,n,p,m) mu(X,W,n,p,m)/sig.mu(X,W,n,p,m)
## tau.std <- function(X,W,n,p,m,xi,nu)
## {
##     return(xi + sqrt(nu)*tau.bar(X,W,n,p,m)*(tau(X,W,n,p,m)-tau.bar(0,W,n,p,m))/sig.tau(X,W,n,p,m))
## }
## out <- function(n,Z,mu,tau,kappa) mu + tau*Z + kappa*rnorm(n)
## ## Function to generate data
## dgp <- function(n,p,m,xi,nu,kappa,x.dist)
## {
##     v <- vars(n,p,x.dist)
##     X <- v$X
##     W <- v$W
##     Z <- v$Z
##     M <- mu.std(X,W,n,p,m)
##     T <- tau.std(X,W,n,p,m,xi,nu)
##     Y <- out(n,Z,M,T,kappa)
##     cate <- tau.std(0,W,n,p,m,xi,nu)
##     ate <- mean(cate)
##     return(list(Y=Y,X=X,W=W,Z=Z,cate=cate,ate=ate))
## }
## ## Saving data
## ### name: DGP_n_m_xi_nu_kappa_s
## for (i in 1:s)
## {
##     for (j in n)
##     {
##         for (k in m)
##         {
##             for (l in xi)
##             {
##                 for (o in nu)
##                 {
##                     for (q in kappa)
##                     {
##                         for (r in x.dist)
##                             {
##                                 data <- dgp(j,p,k,l,o,q,r)
##                                 saveRDS(data,paste0("Data/DGP_",j,"_",k,"_",l,"_",o,"_",q,"_",r,"_",i,".rds"))
##                             }
##                     }
##                 }
##             }
##         }
##     }
## }

