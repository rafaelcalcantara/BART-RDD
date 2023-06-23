## Setup
library(rdrobust)
library(parallel)
library(foreach)
library(doParallel)
devtools::install_github("rafaelcalcantara/XBART@XBCF-RDD")
library(XBART)
library(rdrobust)
###
setwd("~/Documents/Git/XBCF-RDD")
data <- read.csv("Application/headstart.csv")
data <- data[complete.cases(data),]
attach(data)
y <- mort_age59_related_postHS
x <- povrate60
w <- cbind(census1960_pop, census1960_pctsch1417, census1960_pctsch534,
  census1960_pctsch25plus, census1960_pop1417, census1960_pop534,
  census1960_pop25plus, census1960_pcturban, census1960_pctblack)
out <- matrix(NA,9,3)
###
## CCT
cct1 <- rdrobust(y,x,c=59.1968)
cct2 <- rdrobust(y,x,c=59.1968,covs=w)
## XBCF estimation
fit.general <- function(h,y,w,x)
{
    foreach(i=1:length(h),.multicombine=T,.export=c("p","c","Omin","Opct","m","n","num_sweeps","burnin","Nmin","p_categorical")) %dopar%
        {
            fit <- XBCF.rd(y, w, x, c, Owidth = h[i], Omin = Omin, Opct = Opct,
                           num_trees_mod = m, num_trees_con = m,
                           num_cutpoints = n, num_sweeps = num_sweeps,
                           burnin = burnin, Nmin = Nmin,
                           p_categorical_con = p_categorical, p_categorical_mod = p_categorical,
                           tau_con = 2*var(y)/m,
                           tau_mod = 0.5*var(y)/m, parallel=F,random_seed=0)
            predict.XBCFrd(fit,w,rep(0,n))
        }
}
fit.xbcf <- function(y,w,x,p_cat)
{
    t0 <- Sys.time()
    h <- quantile(abs(x),0.125)
    fit <- XBCF.rd(y, w, x, c, Owidth = h, Omin = Omin, Opct = Opct,
                           num_trees_mod = m, num_trees_con = m,
                           num_cutpoints = n, num_sweeps = num_sweeps,
                           burnin = burnin, Nmin = Nmin,
                           p_categorical_con = p_cat, p_categorical_mod = p_cat,
                           tau_con = 2*var(y)/m,
                   tau_mod = 0.5*var(y)/m, parallel=T,
                   nthread = no_cores)
    pred <- predict.XBCFrd(fit,w,rep(0,n))
    post <- colMeans(pred$tau.adj,na.rm=T)
    t1 <- Sys.time()
    dt <- difftime(t1,t0)
    print(paste0("Elapsed time: ",round(dt,2)," seconds"))
    return(list(ate.post=post,pred=fit,Owidth=h,time=dt))
}
###
c             <- 59.1968
n             <- length(y)
Omin          <- as.integer(0.03*n)
Opct          <- 0.9
m             <- 10
Nmin          <- 10
num_sweeps    <- 120
burnin        <- 20
p_categorical <- 0
num_cutpoints <- n
###
### Parallelization
no_cores <- detectCores() - 1
registerDoParallel(no_cores)
###
xbcf1 <- fit.xbcf(y,NULL,x,0)
ate.post1 <- xbcf1$ate.post
xbcf2 <- fit.general(quantile(abs(x),seq(0.05,0.2,0.05)),y,w,x)
ate.post2 <- xbcf2$ate.post
