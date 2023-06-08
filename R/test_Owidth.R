## Setup
library(parallel)
library(foreach)
library(doParallel)
devtools::install_github("rafaelcalcantara/XBART@XBCF-RDD")
library(XBART)
s <- 1
n <- 500
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
                           num_cutpoints = n, num_sweeps = 500,
                           burnin = 0, Nmin = Nmin,
                           p_categorical_con = p_categorical, p_categorical_mod = p_categorical,
                           tau_con = 2*var(y),
                           tau_mod = 0.5*var(y), parallel=F)
            trees <- jsonlite::parse_json(fit$tree_json_mod)$trees
            median(nbot.sweeps(trees,500,1))
        }
}
findOwidth <- function(Owidth,y,w,x)
{
    ## nbots <- sapply(1:5,function(a) unlist(avg.nbot(Owidth,y,w,x)))
    ## nbots <- t(nbots)
    ## nbots <- apply(nbots,2,function(x) sum(x>t))
    ## nbots <- nbots/5
    nbots <- unlist(avg.nbot(Owidth,y,w,x))
    ## nbots <- rbinom(length(Owidth),1,nbots)
    return(Owidth[nbots>4])
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
fit.xbcf <- function(Owidth,y,w,x)
{
    t0 <- Sys.time()
    ## h <- findOwidth(Owidth,y,w,x)
    ## while(length(h)==0) h <- findOwidth(Owidth,y,w,x)
    nbots <- avg.nbot(Owidth,y,w,x)
    nbots <- unlist(nbots)
    h <- Owidth[min(which(nbots==max(nbots)))]
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
Owidth        <- function(x) quantile(x,seq(0.05,0.5,0.05))
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
for (i in s)
{
    print(paste0("Simulation ",i," for DGP 1a"))
    data <- dgp[[i]]
    list2env(data,globalenv())
    fit <- fit.general(Owidth(abs(x)),y,NULL,x)
    ## saveRDS(fit,paste0("Results/xbcf_dgp1a_",i,".rds"))
}
h <- Owidth(abs(x))
pred <- sapply(fit, function(i) mean(colMeans(i$tau.adj)))
plot(h,pred,"b")
nbots <- avg.nbot(h,y,NULL,x)
nbots <- unlist(nbots)
Error <- sqrt((sapply(fit,function(i) mean(colMeans(i$tau.adj)))-data$ate)^2)
plot(nbots,Error,xlab="median nbots")
p1a <- recordPlot()
## DGP1b
for (i in s)
{
    print(paste0("Simulation ",i," for DGP 1b"))
    data <- dgp[[i]]
    list2env(data,globalenv())
    fit <- fit.general(Owidth(abs(x)),y,w,x)
    ## saveRDS(fit,paste0("Results/xbcf_dgp1b_",i,".rds"))
}
h <- Owidth(abs(x))
pred <- sapply(fit, function(i) mean(colMeans(i$tau.adj)))
plot(h,pred,"b")
nbots <- avg.nbot(h,y,w,x)
nbots <- unlist(nbots)
Error <- sqrt((sapply(fit,function(i) mean(colMeans(i$tau.adj)))-data$ate)^2)
plot(apply(nbots,2,median),Error,xlab="median nbots")
p1b <- recordPlot()
## DGP2
dgp <- readRDS("Data/DGP2.rds")
for (i in s)
{
    print(paste0("Simulation ",i," for DGP 2"))
    data <- dgp[[i]]
    list2env(data,globalenv())
    fit <- fit.general(Owidth(abs(x)),y,w,x)
    ## saveRDS(fit,paste0("Results/xbcf_dgp2_",i,".rds"))
}
h <- Owidth(abs(x))
pred <- sapply(fit, function(i) mean(colMeans(i$tau.adj)))
plot(h,pred,"b")
nbots <- avg.nbot(h,y,w,x)
nbots <- unlist(nbots)
Error <- sqrt((sapply(fit,function(i) mean(colMeans(i$tau.adj)))-data$ate)^2)
plot(nbots,Error,xlab="median nbots")
p2 <- recordPlot()
## DGP3
dgp <- readRDS("Data/DGP3.rds")
for (i in s)
{
    print(paste0("Simulation ",i," for DGP 3"))
    data <- dgp[[i]]
    list2env(data,globalenv())
    fit <- fit.general(Owidth(abs(x)),y,w,x)
    ## saveRDS(fit,paste0("Results/xbcf_dgp3_",i,".rds"))
}
h <- Owidth(abs(x))
pred <- sapply(fit, function(i) mean(colMeans(i$tau.adj)))
plot(h,pred,"b")
nbots <- avg.nbot(h,y,w,x)
nbots <- unlist(nbots)
Error <- sqrt((sapply(fit,function(i) mean(colMeans(i$tau.adj)))-data$ate)^2)
plot(apply(nbots,2,median),Error,xlab="median nbots")
p3 <- recordPlot()
####
cowplot::plot_grid(plotlist=list(p1a,p1b,p2,p3),nrow=2,ncol=2)
####
stopImplicitCluster()
####
## nbots <- do.call("rbind",avg.nbot(Owidth,y,w,x))
## nbots <- t(nbots)
## colnames(nbots) <- Owidth
## par(mfrow=c(1,2))
## boxplot(nbots)
## plot(Owidth,colMeans(nbots),type="b")
## par(mfrow=c(1,1))
####
## a <- fit$pred
## a <- sapply(a,function(x) colMeans(x$tau.adj))
## colnames(a) <- fit$Owidth
## boxplot(a)
## abline(h=data$ate,lty=1)
## mean(fit$ate.post)
## quantile(fit$ate.post,c(0.025,0.975))
