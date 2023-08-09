## Setup
library(parallel)
library(foreach)
library(doParallel)
devtools::install_github("JingyuHe/XBART@XBCF-RDD")
library(XBART)
model  <- 3:6
xi <- nu <- kappa <- c(0.25,2)
## Parallelization
no_cores <- detectCores() - 1
registerDoParallel(no_cores)
c             <- 0
Omin          <- 10
h             <- 0.1
Opct          <- 0.9
ntrees        <- 10
Nmin          <- 10
num_sweeps    <- 120
burnin        <- 20
p_categorical <- 1
## ## h
## fit <- function(s,n,m,xi,nu,kappa)
## {
##     foreach(i=0:s,.multicombine=T,.export=c("n","c","Omin","h","Opct","m","n","num_sweeps","burnin","Nmin","p_categorical")) %dopar%
##         {
##             h <- 0.05 + i/20
##             print(paste0("N=",n,"; Model=",m,"; xi=",xi,"; nu=",nu,"; kappa=",kappa,"; h=",h))
##             mse <- rep(0,5)
##             for (sample in 1:5)
##             {
##                 data <- readRDS(paste0("Data/DGP_",n,"_",m,"_",xi,"_",nu,"_",kappa,"_",sample))
##                 fit <- XBCF.rd(data$Y, data$W, data$X, c,
##                                Owidth = h, Omin = Omin, Opct = Opct,
##                                num_trees_mod = ntrees,
##                                num_trees_con = ntrees,
##                                num_cutpoints = n,
##                                num_sweeps = num_sweeps,
##                                burnin = burnin, Nmin = Nmin,
##                                p_categorical_con = p_categorical,
##                                p_categorical_mod = p_categorical,
##                                tau_con = 2*var(data$Y)/ntrees,
##                                tau_mod = 0.5*var(data$Y)/ntrees, parallel=F)
##                 test <- -h<=data$X & data$X<=h
##                 pred <- predict.XBCFrd(fit,data$W[test,],rep(0,sum(test)))
##                 pred <- pred$tau.adj[,(burnin+1):num_sweeps]
##                 pred <- colMeans(pred)
##                 mse[sample] <- (mean(pred)-xi)^2
##             }
##             return(sqrt(mean(mse)))
##         }
## }
## ###
## results.h <- vector("list",32)
## index <- 0
## for (j in model)
## {
##     for (k in xi)
##     {
##         for (l in nu)
##         {
##             for (m in kappa)
##             {
##                 index <- index+1
##                 results.h[[index]] <- fit(4,500,j,k,l,m)
##             }
##         }
##     }
## }
## ###
## saveRDS(results.h,"Tables/rmse_h.rds")
## rmse <- apply(do.call("rbind",lapply(results.h,unlist)),2,function(x) mean(sqrt(x)))
## rmse <- data.frame(h=c(0.05,0.1,0.15,0.2,0.25),RMSE=rmse)
## ###
## yticks <- seq(min(rmse$RMSE),max(rmse$RMSE),(-min(rmse$RMSE)+max(rmse$RMSE))/5)
## png("Figures/rmse_h.png")
## plot(rmse,type="b",col="blue",pch=19,cex=0.75,axes=F)
## axis(1,at=rmse$h,labels=F)
## text(x=rmse$h,par("usr")[3]-0.002,labels=rmse$h,pos=1,xpd=T,cex=0.75)
## axis(2,at=yticks,labels=F)
## text(y=yticks,par("usr")[1]-0.002,labels=round(yticks,2),pos=2,xpd=T,cex=0.75)
## dev.off()
## Omin
fit <- function(s,n,m,xi,nu,kappa)
{
    foreach(i=1:s,.multicombine=T,.export=c("n","c","Omin","h","Opct","m","n","num_sweeps","burnin","Nmin","p_categorical")) %dopar%
        {
            Omin <- i*(i==1) + 5*(i-1)*(i>1)
            print(paste0("N=",n,"; Model=",m,"; xi=",xi,"; nu=",nu,"; kappa=",kappa,"; Omin=",Omin))
            mse <- rep(0,5)
            for (sample in 1:5)
            {
                data <- readRDS(paste0("Data/DGP_",n,"_",m,"_",xi,"_",nu,"_",kappa,"_",sample))
                fit <- XBCF.rd(data$Y, data$W, data$X, c,
                               Owidth = h, Omin = Omin, Opct = Opct,
                               num_trees_mod = ntrees,
                               num_trees_con = ntrees,
                               num_cutpoints = n,
                               num_sweeps = num_sweeps,
                               burnin = burnin, Nmin = Nmin,
                               p_categorical_con = p_categorical,
                               p_categorical_mod = p_categorical,
                               tau_con = 2*var(data$Y)/ntrees,
                               tau_mod = 0.5*var(data$Y)/ntrees,
                               parallel=F)
                test <- -h<=data$X & data$X<=h
                pred <- predict.XBCFrd(fit,data$W[test,],rep(0,sum(test)))
                pred <- pred$tau.adj[,(burnin+1):num_sweeps]
                pred <- colMeans(pred)
                mse[sample] <- (mean(pred)-xi)^2
            }
            return(sqrt(mean(mse)))
        }
}
###
results.omin <- vector("list",32)
index <- 0
for (j in model)
{
    for (k in xi)
    {
        for (l in nu)
        {
            for (m in kappa)
            {
                index <- index+1
                results.omin[[index]] <- fit(5,500,j,k,l,m)
            }
        }
    }
}
###
saveRDS(results.omin,"Tables/rmse_omin.rds")
rmse <- apply(do.call("rbind",lapply(results.omin,unlist)),2,function(x) mean(sqrt(x)))
rmse <- data.frame(N=c(1,5,10,15,20),RMSE=rmse)
###
yticks <- seq(min(rmse$RMSE),max(rmse$RMSE),(-min(rmse$RMSE)+max(rmse$RMSE))/5)
png("Figures/rmse_omin.png")
plot(rmse,type="b",col="blue",pch=19,cex=0.75,axes=F)
axis(1,at=rmse$N,labels=F)
text(x=rmse$N,par("usr")[3]-0.002,labels=rmse$N,pos=1,xpd=T,cex=0.75)
axis(2,at=yticks,labels=F)
text(y=yticks,par("usr")[1]-0.002,labels=round(yticks,2),pos=2,xpd=T,cex=0.75)
dev.off()
