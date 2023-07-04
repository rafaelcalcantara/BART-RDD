## Setup
library(parallel)
library(foreach)
library(doParallel)
devtools::install_github("jingyuhe/XBART@XBCF-RDD")
library(XBART)
library(rpart)
library(rpart.plot)
library(xtable)
## Data
set.seed(0)
n <- 500
x <- 2*rbeta(n,2,4)-1
w1 <- matrix(rnorm(n*2,0,0.5),n,2)
w2 <- rbinom(n,1,0.7)
## w2 <- rnorm(n,0,0.5)
w <- cbind(w1,w2)
z <- x>=0
mu <- function(x,w1,w2) x/(1+x) + rowSums(w1)
tau <- function(x,w1,w2) 1/(2+exp(-0.5*x)) + .1*rowSums(w1) + 0.5*w2
y <- mu(x,w1,w2) + tau(x,w1,w2)*z + rnorm(n)
ate <- mean(tau(0,w1,w2))
cate0 <- mean(tau(0,w1,w2)[w2==0])
cate1 <- mean(tau(0,w1,w2)[w2==1])
###
c             <- 0
Omin          <- as.integer(0.03*n)
h             <- quantile(abs(x),0.2)
Opct          <- 0.9
m             <- 10
Nmin          <- 10
num_sweeps    <- 120
burnin        <- 20
p_categorical <- 1
num_cutpoints <- n
###
### Parallelization
no_cores <- detectCores() - 1
registerDoParallel(no_cores)
###
fit <- XBCF.rd(y, w, x, c, Owidth = h, Omin = Omin, Opct = Opct,
               num_trees_mod = m, num_trees_con = m,
               num_cutpoints = n, num_sweeps = num_sweeps,
               burnin = burnin, Nmin = Nmin,
               p_categorical_con = p_categorical,
               p_categorical_mod = p_categorical,
               tau_con = 2*var(y)/m,
               tau_mod = 0.5*var(y)/m, parallel=T,
               nthread = no_cores)
pred <- predict.XBCFrd(fit,w,rep(0,n))
###
head(pred$tau.adj[,(burnin+1):num_sweeps])
summary(colMeans(pred$tau.adj[,(burnin+1):num_sweeps]))
ate
summary(colMeans(pred$tau.adj[w2==0,(burnin+1):num_sweeps]))
cate0
summary(colMeans(pred$tau.adj[w2==1,(burnin+1):num_sweeps]))
cate1
fit$importance_treatment
## table(-h<=x & x<=0,w2)
## table(0<=x & x<=h,w2)
Omin
h
###
stopImplicitCluster()
## Rpart
cart <- rpart(y~.,data=data.frame(y=pred$tau.adj.mean,w=w))
rpart.plot(cart)
