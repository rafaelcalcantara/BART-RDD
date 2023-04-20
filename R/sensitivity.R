set.seed(000)
###
library(devtools)
## Install latest version of XBART package if available
install_github("Rafael-C-Alcantara/XBART@XBCF-RDD")
library(XBART)
## Install RATest package (if needed) for Lee (2008) data
## install.packages("RATest")
## library(RATest)
## data(lee2008)
## ### Plot Lee data to see if it matches the paper
## leePlot <- cbind(lee2008$difdemshare,lee2008$demsharenext)
## leePlot <- leePlot[order(leePlot[,1]),]
## leePlot <- cbind(cut(leePlot[,1],seq(-0.25,0.25,length.out=101)),leePlot)
## leePlot <- do.call("rbind",by(leePlot[,2:3],leePlot[,1],colMeans))
## plot(leePlot,pch=19,col="black")
## abline(v=0,lty=2)
## rm(leePlot)
### Plot matches fig 4.1 (a) in Lee (2008)
## Helper functions
### Create synthetic data based on Lee (2008) data
## leeSynth <- function(data,n)
## {
##     obs <- sample(1:nrow(data),n)
##     data <- data[obs,]
##     out <- list(y=data[,"demsharenext"],x=data[,"difdemshare"],
##                 z=data[,"difdemshare"]>=0,
##                 w=data[,colnames(data) %in% c("difdemshare","demsharenext")==F])
##     out$w <- apply(out$w,2,function(x) (x-min(x))/(max(x)-min(x)))
##     out$y <- out$x*rowMeans(out$w)+out$z*rowMeans(out$w)+rnorm(length(out$y),sd=0.5*sd(out$y))
##     out$tau <- mean(rowMeans(out$w))
##     return(out)
## }
    
### Extract 95% CI and posterior mean from pred.XBCFrd
tau.pred <- function(x,burnin,num_sweeps)
{
    ate.posterior <- colMeans(x$tau.adj[,(burnin+1):num_sweeps])
    return(mean(ate.posterior))
}
### Minimize MSE for XBCF-RDD for a given grid of Owidth
ate.xbcf <- function(g,n,burnin,num_sweeps)
{
    pred.vec <- vector("numeric",length(g))
    for (i in 1:length(g))
    {
        ## XBCF
        fit.XBCFrd <- XBCF.rd(y, w, x, c, Owidth = Owidth[i], Omin = Omin, Opct = Opct,
                              pcat_con = 0, pcat_mod = 0,
                              num_trees_mod = num_trees_mod, num_trees_con = num_trees_con,
                              num_cutpoints = num_cutpoints, num_sweeps = num_sweeps,
                              burnin = burnin, Nmin = Nmin)

### Given the original (w,x) used for fitting the model, we take (w,x=0) as the test set
        pred.vec[i] <- tau.pred(predict.XBCFrd(fit.XBCFrd, W = w, X = rep(0,n)),burnin,num_sweeps)
    }
    return(pred.vec)
}
### Calculate P(-h<X<h)
## probX <- function(x,Owidth)
## {
##     out <- vector("numeric",length(Owidth))
##     for (i in 1:length(Owidth))
##     {
##         out[i] <- sum(x >= -Owidth[i] & x <= Owidth[i])/length(x)
##     }
##     return(out)
## }
### Function to recover top N combinations of Owidth and Omin
## topN <- function(x,N)
## {
##     ind <- match(head(sort(x),N),x)
##     col <- sapply(ind/nrow(x),function(x) ifelse(is.integer(x),x,as.integer(x)+1))
##     row <- ind %% nrow(x)
##     return(cbind("Owidth"=rownames(x)[row],"Omin"=colnames(x)[col]))
## }
## ### Get ratio Omin/P(-h<x<h)
## OminRatio <- function(x)
## {
##     data <- names(x)
##     num <- sapply(data,function(y) as.numeric(x[[y]][,2]))
##     for (i in data)
##     {
##         if (i=="Gaussian")
##         {
##             den.g <- pnorm(as.numeric(x[[i]][,1]),sd=.5) - pnorm(-as.numeric(x[[i]][,1]),sd=.5)
##         }
##         if (i=="Uniform")
##         {
##             den.u <- punif(as.numeric(x[[i]][,1]),min=-2,max=2) - punif(-as.numeric(x[[i]][,1]),min=-2,max=2)
##         }
##     }
##     den <- cbind(den.g,den.u)
##     return(num/den/100)
## }
### mu and tau functions
mu <- function(W, X) return(0.1 * rowSums(W) + 1/(1+exp(-5*X)))
tau <- function(W, X) return( sin(mu(W, X)) + 1)
### Plotting function
plot.res <- function(Owidth,y,dist)
{
    matplot(x=Owidth,y=y,type="b",lty=1,col="blue",main=dist,pch=16)
}
## Estimation settings
### Owidth: overlap bandiwdth
### Omin: minimum number of observations inside overlap region for each leaf node
### Opct: If nb is the number of obs in node b, Opct*nb of them have to be inside overlap region
## Owidth        <- seq(0.06,0.6,0.01)
Owidth <- seq(0.01,0.2,0.01)
Omin          <- 2
Opct          <- 0.9
num_trees_mod <- 10
num_trees_con <- 10
Nmin          <- 10
num_sweeps    <- 20
burnin        <- 10
## Example data
n <- 100
p <- 2
c <- 0 # Cutoff
num_cutpoints <- n
s <- 50

results <- matrix(0,length(Owidth),2)
## pX <- matrix(0,length(Owidth),s,dimnames=list(Owidth,1:s))
for (h in 1:2)
{
    if (h==1)
    {
        dist.x <- function(x) rnorm(x,sd=.5)
        d <- "Gaussian"
    } else if (h==2)
    {
        dist.x <- function(x) runif(x,-2,2)
        d <- "Uniform"
    }
    temp.results <- vector("numeric",length(Owidth))
    for (i in 1:s)
    {
        print(paste0("Distibution of x: ",d))
        print(paste0("Simulation ",i))
        w <- matrix(rexp(n*p), n, p)
        x <- dist.x(n)
        z <- x >= c
        y <- mu(w, x) + tau(w, x)*z + rnorm(n, 0, 0.1)
        true.ate <- mean(tau(w,0))
        ate <- ate.xbcf(Owidth,n,burnin,num_sweeps)
        temp.results <- temp.results+(ate-true.ate)^2
    }
    results[,h] <- temp.results/s
}
###
## temp.results <- matrix(0,length(Owidth),length(Omin))
## for (i in 1:s)
## {
##     print("Lee (2008) data")
##     print(paste0("Simulation ",i))
##     data <- leeSynth(lee2008,n)
##     w <- data$w
##     x <- data$x
##     z <- data$z
##     y <- data$y
##     true.ate <- data$tau
##     ate <- ate.xbcf(Owidth,Omin,n,burnin,num_sweeps)
##     temp.results <- temp.results+(ate-true.ate)^2
##     pX[,i] <- probX(x,Owidth)
## }
## results[["Lee (2008)"]] <- temp.results/s
###
cutpoints <- vector("numeric",length(Owidth))
for (i in 1:length(Owidth))
{
    ## XBCF
    fit.XBCFrd <- XBCF.rd(y, w, x, c, Owidth = Owidth[i],
                          Omin = Omin, Opct = Opct,
                          pcat_con = 0, pcat_mod = 0,
                          num_trees_mod = 1, num_trees_con = 1,
                          num_cutpoints = num_cutpoints,
                          num_sweeps = 1,
                          burnin = 0, Nmin = Nmin)
    trees_json <- jsonlite::fromJSON(fit.XBCFrd$tree_json_mod,simplifyVector=F)
    cutpoints[i] <- trees_json$trees[["0"]][["0"]][["avail.cutpoints"]]
}
results <- cbind(results,cutpoints)
rownames(results) <- Owidth
colnames(results) <- c("Gaussian","Uniform","Cutpoints")
###
results[order(results[,1]),c(1,3)]
results[order(results[,2]),c(2,3)]
###
saveRDS(results,"R/sensitivity.rds")
## saveRDS(pX,"R/prob_X_Lee.rds")
###
## png("MSE_opct.png")
## par(mfrow=c(1,2))
## plot.res(Owidth,results$Gaussian,"Gaussian")
## plot.res(Owidth,results$Uniform,"Uniform")
## par(mfrow=c(1,1))
## dev.off()
###
## top10 <- lapply(results,topN,10)
## top10
## round(OminRatio(top10),2)
## ###
## par(mfrow=c(1,2))
## plot(OminRatio(lapply(results,topN,150))[,1],type="b",ylab=expression(O[min]/N[h]),main="Gaussian")
## plot(OminRatio(lapply(results,topN,150))[,2],type="b",ylab=expression(O[min]/N[h]),main="Uniform")
## par(mfrow=c(1,1))
###
## plot(lee2008$difdemshare,type="n",ylim=c(0,1),xlim=c(min(lee2008$difdemshare),max(lee2008$difdemshare)),xlab="",ylab="")
## ####
## lines(density(lee2008$difdemshare))
## abline(v=-as.numeric(top10[["Lee (2008)"]][1,1]))
## abline(v=as.numeric(top10[["Lee (2008)"]][1,1]))
## ####
## lines(density(rnorm(nrow(lee2008),sd=.5)),col="blue")
## abline(v=-as.numeric(top10[["Gaussian"]][1,1]),col="blue")
## abline(v=as.numeric(top10[["Gaussian"]][1,1]),col="blue")
## ####
## lines(density(runif(nrow(lee2008),min=-2,max=2)),col="red")
## abline(v=-as.numeric(top10[["Uniform"]][1,1]),col="red")
## abline(v=as.numeric(top10[["Uniform"]][1,1]),col="red")
## ####
## legend("topright",legend=c("Lee (2008)","Gaussian","Uniform"),col=c("black","blue","red"),lty=1)
###
## minMSE <- do.call("rbind",lapply(results,function(x) apply(x,2,min)))
## colnames(minMSE) <- n
## minOmin <- do.call("rbind",lapply(results,function(x) apply(x,2,function(x) Omin[which(x==min(x))])))
## colnames(minOmin) <- n
