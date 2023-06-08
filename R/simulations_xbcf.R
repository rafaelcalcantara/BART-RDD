## Setup
library(parallel)
library(foreach)
library(doParallel)
devtools::install_github("rafaelcalcantara/XBART@XBCF-RDD")
library(XBART)
### Functions to find good Owidth values
params <- function(x)
{
    if (is.null(x$theta)) c(params(x$left),params(x$right))
    else x$theta
}
nbot <- function(t) length(params(t))
nbot.trees <- function(t,m) sapply(0:(m-1),function(y) nbot(t[[as.character(y)]]))
nbot.sweeps <- function(t,num_sweeps,m) sapply(1:num_sweeps, function(y) nbot.trees(t[[y]],m))
avg.nbot <- function(Owidth,y,w,x)
{
    foreach(i=1:length(Owidth),.multicombine=T,.export=c("p","c","Omin","Opct","m","n","num_sweeps","burnin","Nmin","p_categorical")) %dopar%
        {
            fit <- XBCF.rd(y, w, x, c, Owidth = Owidth[i], Omin = Omin, Opct = Opct,
                           num_trees_mod = 1, num_trees_con = 1,
                           num_cutpoints = n, num_sweeps = 100,
                           burnin = burnin, Nmin = Nmin,
                           p_categorical_con = p_categorical, p_categorical_mod = p_categorical,
                           tau_con = 2*var(y),
                           tau_mod = 0.5*var(y), parallel=F)
            trees <- jsonlite::parse_json(fit$tree_json_mod)$trees
            median(nbot.sweeps(trees,100,1))
        }
}
findOwidth <- function(Owidth,y,w,x,t)
{
    nbots <- sapply(1:5,function(a) unlist(avg.nbot(Owidth,y,w,x)))
    nbots <- t(nbots)
    nbots <- apply(nbots,2,function(x) sum(x>t))
    nbots <- nbots/5
    ## nbots <- unlist(avg.nbot(Owidth,y,w,x))
    nbots <- rbinom(length(Owidth),1,nbots)
    return(Owidth[as.logical(nbots)])
}
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
                           tau_mod = 0.5*var(y)/m, parallel=F)
            predict.XBCFrd(fit,w,rep(0,n))
        }
}
fit.xbcf <- function(y,w,x)
{
    t0 <- Sys.time()
    h <- quantile(abs(x),0.11)
    fit <- XBCF.rd(y, w, x, c, Owidth = h, Omin = Omin, Opct = Opct,
                           num_trees_mod = m, num_trees_con = m,
                           num_cutpoints = n, num_sweeps = num_sweeps,
                           burnin = burnin, Nmin = Nmin,
                           p_categorical_con = p_categorical, p_categorical_mod = p_categorical,
                           tau_con = 2*var(y)/m,
                           tau_mod = 0.5*var(y)/m, parallel=F)
    pred <- predict.XBCFrd(fit,w,rep(0,n))
    post <- colMeans(pred$tau.adj,na.rm=T)
    t1 <- Sys.time()
    dt <- difftime(t1,t0)
    print(paste0("Elapsed time: ",round(dt,2)," seconds"))
    return(list(ate.post=post,pred=fit,Owidth=h,time=dt))
}
fit.1a <- function(s,data)
{
    foreach(i=1:s,.multicombine=T,.export=c("p","c","Omin","Opct","m","n","num_sweeps","burnin","Nmin","p_categorical")) %dopar%
        {
            print(paste0("Simulation ",i," for DGP 1a"))
            fit <- fit.xbcf(data[[i]]$y,NULL,data[[i]]$x)
            saveRDS(fit,paste0("Results/xbcf_dgp1a_",i,".rds"))
        }
}
fit.1b <- function(s,data)
{
    foreach(i=1:s,.multicombine=T,.export=c("p","c","Omin","Opct","m","n","num_sweeps","burnin","Nmin","p_categorical")) %dopar%
        {
            print(paste0("Simulation ",i," for DGP 1b"))
            fit <- fit.xbcf(data[[i]]$y,data[[i]]$w,data[[i]]$x)
            saveRDS(fit,paste0("Results/xbcf_dgp1b_",i,".rds"))
        }
}
fit.2 <- function(s,data)
{
    foreach(i=1:s,.multicombine=T,.export=c("p","c","Omin","Opct","m","n","num_sweeps","burnin","Nmin","p_categorical")) %dopar%
        {
            print(paste0("Simulation ",i," for DGP 2"))
            fit <- fit.xbcf(data[[i]]$y,data[[i]]$w,data[[i]]$x)
            saveRDS(fit,paste0("Results/xbcf_dgp2_",i,".rds"))
        }
}
fit.3 <- function(s,data)
{
    foreach(i=1:s,.multicombine=T,.export=c("p","c","Omin","Opct","m","n","num_sweeps","burnin","Nmin","p_categorical")) %dopar%
        {
            print(paste0("Simulation ",i," for DGP 3"))
            fit <- fit.xbcf(data[[i]]$y,data[[i]]$w,data[[i]]$x)
            saveRDS(fit,paste0("Results/xbcf_dgp3_",i,".rds"))
        }
}
####
c             <- 0
## Owidth        <- function(x) quantile(x,seq(0.05,0.5,0.05))
Omin          <- as.integer(0.03*n)
Opct          <- 0.9
m             <- 10
Nmin          <- as.integer(0.02*n)
num_sweeps    <- 50
burnin        <- 20
p_categorical <- 0
num_cutpoints <- n
### Parallelization
no_cores <- detectCores() - 1
registerDoParallel(no_cores)
## DGP1a
dgp <- readRDS("Data/DGP1.rds")
fit.1a(s,dgp)
## DGP1b
fit.1b(s,dgp)
## DGP2
dgp <- readRDS("Data/DGP2.rds")
fit.2(s,dgp)
## DGP3
dgp <- readRDS("Data/DGP3.rds")
fit.3(s,dgp)
####
stopImplicitCluster()
