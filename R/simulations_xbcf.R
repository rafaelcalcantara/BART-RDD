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
    foreach(i=1:length(Owidth),.multicombine=T,.export=c("c","Omin","Opct","m","n","num_sweeps","burnin","Nmin","p_categorical")) %dopar%
        {
            fit <- XBCF.rd(y, w, x, c, Owidth = Owidth[i], Omin = Omin, Opct = Opct,
                           num_trees_mod = 1, num_trees_con = 1,
                           num_cutpoints = n, num_sweeps = num_sweeps,
                           burnin = burnin, Nmin = Nmin,
                           p_categorical_con = p_categorical, p_categorical_mod = p_categorical,
                           tau_con = 1,
                           tau_mod = 1, parallel=F)
            trees <- jsonlite::parse_json(fit$tree_json_mod)$trees
            min(nbot.sweeps(trees,num_sweeps,1))
        }
}
findOwidth <- function(Owidth,y,w,x)
{
    nbots <- unlist(avg.nbot(Owidth,y,w,x))
    return(Owidth[nbots>2])
}
fit.general <- function(h,y,w,x)
{
    foreach(i=1:length(h),.multicombine=T,.export=c("c","Omin","Opct","m","n","num_sweeps","burnin","Nmin","p_categorical")) %dopar%
        {
            fit <- XBCF.rd(y, w, x, c, Owidth = h[i], Omin = Omin, Opct = Opct,
                           num_trees_mod = m, num_trees_con = m,
                           num_cutpoints = n, num_sweeps = num_sweeps,
                           burnin = burnin, Nmin = Nmin,
                           p_categorical_con = p_categorical, p_categorical_mod = p_categorical,
                           tau_con = 1/m,
                           tau_mod = 1/m, parallel=F)
            predict.XBCFrd(fit,w,rep(0,n))
        }
}
fit.xbcf <- function(Owidth,y,w,x)
{
    t0 <- Sys.time()
    h <- findOwidth(Owidth,y,w,x)
    while(length(h)==0) h <- findOwidth(Owidth,y,w,x)
    fit <- fit.general(h,y,w,x)
    pred <- lapply(1:length(h), function(x) fit[[x]]$tau.adj)
    pred <- do.call("cbind",pred)
    post <- colMeans(pred,na.rm=T)
    t1 <- Sys.time()
    dt <- difftime(t1,t0)
    print(paste0("Selected ",length(h)," out of ",length(Owidth)," values provided for Owidth"))
    print(paste0("Elapsed time: ",round(dt,2)," seconds"))
    return(list(ate.post=post,pred=fit,Owidth=h,time=dt))
}
####
c             <- 0
Owidth        <- seq(0.05,1,0.05)
Omin          <- 2
Opct          <- 0.9
m             <- 10
Nmin          <- 10
num_sweeps    <- 20
burnin        <- 10
p_categorical <- 0
num_cutpoints <- n
### Parallelization
no_cores <- detectCores() - 1
registerDoParallel(no_cores)
## DGP1a
dgp <- readRDS("Data/DGP1.rds")
for (i in 1:s)
{
    print(paste0("Simulation ",i," for DGP 1a"))
    data <- dgp[[i]]
    list2env(data,globalenv())
    fit <- fit.xbcf(Owidth,y,NULL,x)
    saveRDS(fit,paste0("Results/xbcf_dgp1a_",i,".rds"))
}
## DGP1b
for (i in 1:s)
{
    print(paste0("Simulation ",i," for DGP 1b"))
    data <- dgp[[i]]
    list2env(data,globalenv())
    fit <- fit.xbcf(Owidth,y,w,x)
    saveRDS(fit,paste0("Results/xbcf_dgp1b_",i,".rds"))
}
## DGP2
dgp <- readRDS("Data/DGP2.rds")
for (i in 1:s)
{
    print(paste0("Simulation ",i," for DGP 2"))
    data <- dgp[[i]]
    list2env(data,globalenv())
    fit <- fit.xbcf(Owidth,y,w,x)
    saveRDS(fit,paste0("Results/xbcf_dgp2_",i,".rds"))
}
## DGP3
dgp <- readRDS("Data/DGP3.rds")
for (i in 1:s)
{
    print(paste0("Simulation ",i," for DGP 3"))
    data <- dgp[[i]]
    list2env(data,globalenv())
    fit <- fit.xbcf(Owidth,y,w,x)
    saveRDS(fit,paste0("Results/xbcf_dgp3_",i,".rds"))
}
####
stopImplicitCluster()
