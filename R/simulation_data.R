#### NOTES
##### TEST DIFFERENT DIST OF X
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
s <- 1000
n <- 500
c <- 0
### DGP
sig1 <- 0.05
sig2 <- 0.05
sig3 <- 0.05
u1 <- u2 <- u3 <- matrix(runif(n*s,0,1),n,s)
w1 <- w2 <- w3 <- vector("list",s)
for (i in 1:s)
{
    W1 <- runif(n,u1[,i],u1[,i]+1)
    W2 <- runif(n,0,0.5)
    W3 <- rbinom(n,2,u1[,i])+1
    W4 <- rbinom(n,1,0.6)+1
    w1[[i]] <- cbind(W1,W2,W3,W4)
}
for (i in 1:s)
{
    W1 <- runif(n,u2[,i],u2[,i]+1)
    W2 <- runif(n,0,0.5)
    W3 <- rbinom(n,2,u2[,i])+1
    W4 <- rbinom(n,1,0.6)+1
    w2[[i]] <- cbind(W1,W2,W3,W4)
}
for (i in 1:s)
{
    W1 <- runif(n,u3[,i],u3[,i]+1)
    W2 <- runif(n,0,0.5)
    W3 <- rbinom(n,2,u3[,i])+1
    W4 <- rbinom(n,1,0.6)+1
    w3[[i]] <- cbind(W1,W2,W3,W4)
}
x1 <- matrix(2*rbeta(n*s,2,4)-u1[,i]-0.2,n,s)
x2 <- matrix(2*rbeta(n*s,2,4)-u2[,i]-0.2,n,s)
x3 <- matrix(2*rbeta(n*s,2,4)-u3[,i]-0.2,n,s)
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
    if(m==1) out <- 0.1875*sin((w[,1]+w[,2])*pi) + 1 + 1.875*x - 1.25*x^2 + 1.75*x^3
    if(m==2) out <- 0.1*w[,4]*sin((w[,1]+w[,2])*pi) + 1 + 0.9*w[,4] + w[,4]*x - 0.9*w[,4]*x^2 + w[,4]*x^3
    if(m==3) out <- 0.2*w34*sin((w[,1]+w[,2])*pi) + 1 + 2*w34 + 2.27*w34*x - 1.13*w34*x^2 + 2*w34*x^3
    return(out)
}
tau <- function(x,w,m)
{
    w34 <- (-1)^(w[,4])*w[,3]
    w34 <- 0.1 + (w34-min(w34))/(max(w34)-min(w34))
    w34 <- round(w34,digits=1)
    if(m==1) out <- 0.025*cos((w[,1]+w[,2])*pi) + 0.05 - 2.8*x + 1.4*x^2 - 0.14*x^3
    if(m==2) out <- 0.0125*w[,4]*cos((w[,1]+w[,2])*pi) + 0.05 - 1.8*w[,4]*x + 0.9*w[,4]*x^2 - 0.09*w[,4]*x^3
    if(m==3) out <- 0.05*w34*cos((w[,1]+w[,2])*pi) + 0.05 - 3.4*w34*x + 1.7*w34*x^2 - 0.17*w34*x^3
    return(out)
}
y1 <- y2 <- y3 <- cate1 <- cate2 <- cate3 <- matrix(0,n,s)
for (i in 1:s)
{
    y1[,i] <- mu(x1[,i],w1[[i]],1) + tau(x1[,i],w1[[i]],1)*z1[,i] + 0.5*u1[,i] + e1[,i]
    cate1[,i] <- tau(0,w1[[i]],1)
    y2[,i] <- mu(x2[,i],w2[[i]],2) + tau(x2[,i],w2[[i]],2)*z2[,i] + 0.5*u2[,i] + e2[,i]
    cate2[,i] <- tau(0,w2[[i]],2)
    y3[,i] <- mu(x3[,i],w3[[i]],3) + tau(x3[,i],w3[[i]],3)*z3[,i] + 0.5*u3[,i] + e3[,i]
    cate3[,i] <- tau(0,w3[[i]],3)
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
## Saving data
saveRDS(list(y=y1,x=x1,w=w1,z=z1,c=c,cate=cate1), "Data/dgp1.rds")
saveRDS(list(y=y2,x=x2,w=w2,z=z2,c=c,cate=cate2), "Data/dgp2.rds")
saveRDS(list(y=y3,x=x3,w=w3,z=z3,c=c,cate=cate3), "Data/dgp3.rds")
