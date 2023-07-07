## Setup
library(parallel)
library(foreach)
library(doParallel)
devtools::install_github("rafaelcalcantara/XBART@XBCF-RDD")
library(XBART)
library(rpart)
library(rpart.plot)
library(xtable)
## Data
set.seed(0)
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
attach(dgp6.fun(500,0))
###
n             <- 500
c             <- 0
Omin          <- as.integer(0.03*n)
h             <- quantile(abs(x),0.2)
Opct          <- 0.9
m             <- 10
Nmin          <- 10
num_sweeps    <- 1020
burnin        <- 20
p_categorical <- 0
num_cutpoints <- n
###
fit <- XBCF.rd(y, w, x, c, Owidth = h, Omin = Omin, Opct = Opct,
               num_trees_mod = m, num_trees_con = m,
               num_cutpoints = n, num_sweeps = num_sweeps,
               burnin = burnin, Nmin = Nmin,
               p_categorical_con = p_categorical,
               p_categorical_mod = p_categorical,
               tau_con = 2*var(y)/m,
               tau_mod = 0.5*var(y)/m)
pred <- predict.XBCFrd(fit,w,rep(0,n))
pred <- pred$tau.adj[,(burnin+1):num_sweeps]
ind.pred <- rowMeans(pred)
## Rpart
cart <- rpart(y~.,data=data.frame(y=ind.pred,w1=w[,1],
                                  w2=w[,2],w3=w[,3],w4=w[,4]))
png("Figures/cart_example.png")
rpart.plot(cart)
dev.off()
###
by(tau.x,w[,1] < -0.048 & w[,4] < -0.094 & w[,2] < -0.12,mean)
by(tau.x,w[,1] < -0.048 & w[,4] < -0.094 & w[,2] >= -0.12,mean)
