set.seed(0)
## Setup
library(parallel)
library(foreach)
library(doParallel)
devtools::install_github("rafaelcalcantara/XBART@XBCF-RDD")
library(XBART)
setwd("~/Documents/Git/XBCF-RDD")
## Functions
temp <- function(tree,w,x,x0,y0,x1,y1)
{
    if (length(tree)==0) return()
    if (tree$var == 1)
    {
        segments(x0,w[tree$val],x1,w[tree$val])
        y0.l <- y0
        y1.l <- w[tree$val]
        y0.r <- w[tree$val]
        y1.r <- y1
        x0.l <- x0
        x0.r <- x0
        x1.l <- x1
        x1.r <- x1
    } else
    {
        segments(x[tree$val],y0,x[tree$val],y1)
        x0.l <- x0
        x1.l <- x[tree$val]
        x0.r <- x[tree$val]
        x1.r <- x1
        y0.l <- y0
        y0.r <- y0
        y1.l <- y1
        y1.r <- y1
    }
    temp(tree$left,w,x,x0.l,y0.l,x1.l,y1.l)
    temp(tree$right,w,x,x0.r,y0.r,x1.r,y1.r)
}
plot.trees <- function(w,x,tree,Owidth)
{
    plot(x,w,"n")
    points(x[x>=0],w[x>=0],pch=24,cex=0.6,bg="gray50")
    points(x[x<0],w[x<0],pch=21,cex=0.6,bg="gray50")
    temp(tree,w,x,par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4])
    abline(v=-Owidth,lty=2)
    abline(v=Owidth,lty=2)
    abline(v=0,lty=3)
}   
###
split <- function(tree,var,val)
{
    if (length(tree) != 0) stop("Not a leaf node")
    tree$var <- var
    tree$val <- val
    tree$left <- list()
    tree$right <- list()
    return(tree)
}
nbot <- function(tree)
{
    if (length(tree)==0) return(1)
    nbot(tree$left)+nbot(tree$right)
}
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
nbot.par <- function(Owidth,y,w,x)
{
    foreach(i=1:length(Owidth),.multicombine=T,.export=c("p","c","Omin","Opct","m","n","num_sweeps","burnin","Nmin","p_categorical")) %dopar%
        {
            fit <- XBCF.rd(y, w, x, c, Owidth = Owidth[i], Omin = Omin, Opct = Opct,
                           num_trees_mod = 1, num_trees_con = 1,
                           num_cutpoints = n, num_sweeps = 1000,
                           burnin = 0, Nmin = Nmin,
                           p_categorical_con = p_categorical, p_categorical_mod = p_categorical,
                           tau_con = 2*var(y),
                           tau_mod = 0.5*var(y), parallel=F)
            trees <- jsonlite::parse_json(fit$tree_json_mod)$trees
            nbot.sweeps(trees,1000,1)
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
fit.xbcf <- function(y,w,x)
{
    t0 <- Sys.time()
    h <- quantile(abs(x),0.125)
    fit <- XBCF.rd(y, w, x, c, Owidth = h, Omin = Omin, Opct = Opct,
                   num_trees_mod = m, num_trees_con = m,
                   num_cutpoints = n, num_sweeps = num_sweeps,
                   burnin = burnin, Nmin = Nmin,
                   p_categorical_con = p_categorical, p_categorical_mod = p_categorical,
                   tau_con = 0.5*var(y)/m,
                   tau_mod = 2*var(y)/m, parallel=F)
    pred <- predict.XBCFrd(fit,w,rep(0,n))
    post <- colMeans(pred$tau.adj,na.rm=T)
    t1 <- Sys.time()
    dt <- difftime(t1,t0)
    print(paste0("Elapsed time: ",round(dt,2)," seconds"))
    return(list(ate.post=post,pred=fit,Owidth=h,time=dt))
}
## Example data: w independent of x
h <- 0.25
n <- 100
x <- runif(n,-1,1)
w <- runif(n,-1,1)
### Plots
#### This plot makes it clear that Owidth=0.25 is too narrow for
#### this dataset
## plot(x,w)
## abline(v=-h,lty=2)
## abline(v=h,lty=2)
## abline(v=0,lty=3)
#### This tree is invalid because it partitions the window
root1 <- list()
root1 <- split(root1,1,order(w)[50])
root1$left <- split(root1$left,2,order(x)[25])
root1$right <- split(root1$right,1,order(w)[75])
root1$left$left <- split(root1$left$left,1,order(w)[20])
root1$left$right <- split(root1$left$right,2,order(x)[65])
root1$right$left <- split(root1$right$left,2,order(x)[40])
root1$right$right <- split(root1$right$right,2,order(x)[60])
###
#### This tree is invalid because, although the partitions don't
#### cut the window, the middle partition has no points inside it
#### to the right of the cutoff. Besides, even if there were points
#### there, it's unlikely that this tree would be valid because of
#### the Opct condition. This condition ensures the partitions are
#### tighter around the window
root2 <- list()
root2 <- split(root2,1,order(w)[50])
root2$left <- split(root2$left,2,order(x)[25])
root2$right <- split(root2$right,1,order(w)[75])
root2$left$left <- split(root2$left$left,1,order(w)[20])
root2$left$right <- split(root2$left$right,2,order(x)[65])
root2$right$left <- split(root2$right$left,2,order(x)[20])
root2$right$right <- split(root2$right$right,2,order(x)[60])
###
#### This tree is invalid because, although the partitions don't
#### cut the window and the middle partition has points inside it
#### to the right of the cutoff, there are many more points
#### in the middle partition outside the window, violating the
#### Opct condition
root3 <- list()
root3 <- split(root3,1,order(w)[50])
root3$left <- split(root3$left,2,order(x)[20])
root3$right <- split(root3$right,2,order(x)[10])
root3$left$left <- split(root3$left$left,1,order(w)[20])
root3$left$right <- split(root3$left$right,2,order(x)[70])
root3$right$left <- split(root3$right$left,2,order(x)[25])
root3$right$right <- split(root3$right$right,2,order(x)[65])
###
root4 <- list()
root4 <- split(root4,1,order(w)[50])
root4$left <- split(root4$left,2,order(x)[20])
root4$right <- split(root4$right,2,order(x)[10])
root4$left$left <- split(root4$left$left,1,order(w)[20])
root4$left$right <- split(root4$left$right,2,order(x)[71])
root4$right$left <- split(root4$right$left,2,order(x)[24])
root4$right$right <- split(root4$right$right,2,order(x)[68])
###
png("Figures/bad_trees_1a.png")
plot.trees(w,x,root1,h)
dev.off()
png("Figures/bad_trees_1b.png")
plot.trees(w,x,root2,h)
dev.off()
png("Figures/bad_trees_1c.png")
plot.trees(w,x,root3,h)
dev.off()
###
png("Figures/good_tree_1.png")
plot.trees(w,x,root4,0.36)
dev.off()
## Fit
n <- 5000
x <- rnorm(n,0,0.25)
w <- rnorm(n,0,0.25)
z <- x>=0
mu.fun <- function(W, X){return(0.1 * w + 1/(1+exp(-5*X)))} 
tau.fun <- function(W, X) return( sin(mu.fun(W, X)) +1) # make sure the treatment effect is non-zero
y <- mu.fun(w, x) + tau.fun(w, x)*z + rnorm(n, 0, 0.2)
####
c             <- 0
Owidth        <- quantile(abs(x),seq(0.01,1,0.01))
Omin          <- as.integer(0.03*n)
Opct          <- 0.9
m             <- 10
Nmin          <- 10
num_sweeps    <- 70
burnin        <- 20
p_categorical <- 0
num_cutpoints <- n
## Error
## ### Parallelization
no_cores <- detectCores() - 1
registerDoParallel(no_cores)
###
fit <- fit.general(Owidth,y,w,x)
stopImplicitCluster()
###
pred <- sapply(fit,function(i) colMeans(i$tau.adj)[(burnin+1):num_sweeps])
tau.var <- apply(pred,2,var)
tau.hat <- apply(pred,2,function(i) c(mean(i),quantile(i,c(0.025,0.975))))
tau.hat <- t(tau.hat)
Error <- (tau.hat[,1]-mean(tau.fun(w,0)))^2
###
png("Figures/example_error.png")
plot(Owidth[1:75],Error[1:75],type="b",lty=1,pch=19,col="blue",cex=0.75,xlab="h",ylab="Squared error")
dev.off()
png("Figures/example_variance.png")
plot(Owidth,tau.var,type="l",ylab="Posterior variance",xlab="h")
dev.off()
###
fit <- fit.xbcf(y,w,x)
quantile(fit$ate.post,0.975)-quantile(fit$ate.post,0.025)
(mean(fit$ate.post)-mean(tau.fun(w,0)))^2
pred <- predict.XBCFrd(fit$pred,w,x)
png("Figures/good_owidth.png")
matplot(sort(x),cbind(rowMeans(pred$tau.adj[,(burnin+1):num_sweeps]),tau.fun(w,x))[order(x),],type=c("p","l"),pch=1,lty=1,col=c("blue","black"),ylab=expression(tau(w,x)),xlab="X")
abline(v=fit$Owidth,lty=2)
abline(v=-fit$Owidth,lty=2)
abline(v=0,lty=3)
dev.off()
###
## Error <- (sapply(fit, function(i) mean(colMeans(i$tau.adj)))-mean(tau.fun(w,0)))^2
## Error <- sqrt(Error)
## nbots <- avg.nbot(Owidth,y,w,x)
## nbots <- unlist(nbots)
## png("Figures/error.png")
## plot(Owidth,Error,"b",xlab="h")
## dev.off()
## stopImplicitCluster()
## ## Good Owidth
## fit <- XBCF.rd(y, w, x, c=0, Owidth = Owidth[9], Omin = 10, Opct = 0.95,
##                num_trees_mod = 10, num_trees_con = 10,
##                num_cutpoints = n, num_sweeps = 100,
##                burnin = 10, Nmin = 20,
##                p_categorical_con = 0, p_categorical_mod = 0,
##                tau_con = 2*var(y)/10,
##                tau_mod = 0.5*var(y)/10, parallel=F,
##                random_seed=0)
## tau <- predict.XBCFrd(fit,w,rep(0,n))
## mean((colMeans(tau$tau.adj) - mean(tau.fun(w,0)))^2)
## tau <- predict.XBCFrd(fit,w,x)
## ###
## png("Figures/good_owidth.png")
## plot(sort(x),tau$tau.adj.mean[order(x)],col="blue",
##      xlab = "x", ylab=expression(tau(w,x)),type="p",
##      ylim = range(tau$tau.adj.mean,tau.fun(w,x)))
## lines(sort(x),tau.fun(w,x)[order(x)])
## abline(v=-0.1,lty=2)
## abline(v=0.1,lty=2)
## dev.off()
## ###
## nbots <- nbot.par(Owidth,y,w,x)
## nbots <- do.call("rbind",nbots)
## plot(apply(nbots,1,median),Error,xlab="Median nbot")
#### Maybe counting depths is better than IQR?
###
## ## Example data: w and x related
## ### In this setting, the partitions must reflect the relationship
## ### between w and x otherwise we end up with empty nodes.
## ### For example, in this case high w equals high x, meaning
## ### partitions at high values of w are less likely to feature
## ### x near the cutoff and thus inside the window. Similarly for
## ### low w and low x.
## h <- 0.25
## n <- 100
## x <- runif(n,-1,1)
## w <- sapply(x,function(i) rnorm(1,i,0.5))
## y <- w + x + (x>=0)*(3 + w + x) + rnorm(n)
## ### Plots
## #### This plot makes it clear that Owidth=0.25 is too narrow for
## #### this dataset
## plot(x,w)
## abline(v=-h,lty=2)
## abline(v=h,lty=2)
## abline(v=0,lty=3)
## #### This tree is invalid because it creates empty nodes in window
## root1 <- list()
## root1 <- split(root1,1,order(w)[49])
## root1$left <- split(root1$left,2,order(x)[30])
## root1$right <- split(root1$right,2,order(x)[32])
## root1$left$left <- split(root1$left$left,1,order(w)[15])
## root1$left$right <- split(root1$left$right,2,order(x)[65])
## root1$right$left <- split(root1$right$left,2,order(x)[25])
## root1$right$right <- split(root1$right$right,1,order(w)[85])
## ###
## root2 <- list()
## root2 <- split(root2,1,order(w)[49])
## root2$left <- split(root2$left,2,order(x)[30])
## root2$right <- split(root2$right,2,order(x)[32])
## root2$left$left <- split(root2$left$left,1,order(w)[15])
## root2$left$right <- split(root2$left$right,2,order(x)[65])
## root2$right$left <- split(root2$right$left,2,order(x)[25])
## root2$right$right <- split(root2$right$right,2,order(x)[67])
## ###
## png("Figures/bad_tree_2.png")
## plot.trees(w,x,root1,h)
## dev.off()
## png("Figures/good_tree_2.png")
## plot.trees(w,x,root2,h)
## dev.off()
## ###
