## Setup
library(parallel)
library(foreach)
library(doParallel)
devtools::install_github("Rafael-C-Alcantara/XBART@XBCF-RDD")
library(XBART)
## Helper functions
tau.pred <- function(x,burnin,num_sweeps)
{
    ate.posterior <- colMeans(x$tau.adj[,(burnin+1):num_sweeps])
    return(c(mean(ate.posterior),quantile(ate.posterior,c(.025,.975))))
}
### Return ATE prediction for XBCF-RDD
pred.ate.xbcf <- function()
{
    ## XBCF
    fit.XBCFrd <- XBCF.rd(y, w, x, c, Owidth = Owidth, Omin = Omin, Opct = Opct,
                          pcat_con = 0, pcat_mod = 0,
                          num_trees_mod = num_trees_mod, num_trees_con = num_trees_con,
                          num_cutpoints = num_cutpoints, num_sweeps = num_sweeps,
                          burnin = burnin, Nmin = Nmin)

### Given the original (w,x) used for fitting the model, we take (w,x=0) as the test set
    pred <- predict.XBCFrd(fit.XBCFrd, W = w, X = rep(0,n))

    ## ATE summary for pred.list
    return(tau.pred(pred,burnin,num_sweeps))
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
                           tau_mod = 0.5*var(y)/m, parallel=F,random_seed=0)
            predict.XBCFrd(fit,w,rep(0,n))
        }
}
fit.xbcf <- function(y,w,x)
{
    if (is.null(w)) p_categorical=0
    t0 <- Sys.time()
    h <- quantile(abs(x),0.125)
    fit <- XBCF.rd(y, w, x, c, Owidth = h, Omin = Omin, Opct = Opct,
                   num_trees_mod = m, num_trees_con = m,
                   num_cutpoints = n, num_sweeps = num_sweeps,
                   burnin = burnin, Nmin = Nmin,
                   p_categorical_con = p_categorical, p_categorical_mod = p_categorical,
                   tau_con = 0.5*var(y)/m,
                   tau_mod = 2*var(y)/m, parallel=F, random_seed=0)
    pred <- predict.XBCFrd(fit,w,rep(0,n))
    post <- colMeans(pred$tau.adj,na.rm=T)
    t1 <- Sys.time()
    dt <- difftime(t1,t0)
    print(paste0("Elapsed time: ",round(dt,2)," seconds"))
    return(list(ate.post=post,pred=fit,Owidth=h,time=dt))
}
##
data <- readstata13::read.dta13("Data/Merged_1853_1864_data.dta")
data$dist_netw <- data$dist_netw/100
data$dist_2 <- ifelse(data$broad==0,-data$dist_netw,data$dist_netw)
vars <- c("log_rentals_1864","dist_2","dist_cent","dist_square",
          "dist_fire","dist_thea","dist_police","dist_urinal",
          "dist_pub","dist_church","dist_bank",
          "dist_vent","dist_pump","dist_pit_fake",
          "no_sewer","old_sewer")
covs <- c("dist_cent","dist_square",
          "dist_fire","dist_thea","dist_police","dist_urinal",
          "dist_pub","dist_church","dist_bank",
          "dist_vent","dist_pump","dist_pit_fake",
          "no_sewer","old_sewer")
data <- data[,vars]
data <- data[complete.cases(data),]
y <- data$log_rentals_1864
x <- data$dist_2
w <- data[,covs]
n <- nrow(data)
c <- 0
##
l0 <- loess(y~x,data=data.frame(y=y,x=x)[x<0,])
l1 <- loess(y~x,data=data.frame(y=y,x=x)[x>=0,])
l0 <- predict(l0,se=T)
l1 <- predict(l1,se=T)
y0.hat <- cbind(l0$fit,l0$fit-1.96*l0$se.fit)
y0.hat <- cbind(y0.hat,l0$fit+1.96*l0$se.fit)
y1.hat <- cbind(l1$fit,l1$fit-1.96*l1$se.fit)
y1.hat <- cbind(y1.hat,l1$fit+1.96*l1$se.fit)
png("Figures/cholera.png")
plot(x,y,type="p",pch=21,lwd=0.5,col="azure4",bg="azure",axes=F,
     xlab="Distance to Broad St. Pump",
     ylab="(Log) Rental Prices in 1864")
axis(1)
axis(2)
matlines(sort(x[x<0]),y0.hat[order(x[x<0]),],
         lwd=1.5,col="gray10",lty=c(1,2,2))
matlines(sort(x[x>=0]),y1.hat[order(x[x>=0]),],
         lwd=1.5,col="gray10",lty=c(1,2,2))
abline(v=0,lty=2)
dev.off()
##
cct1 <- rdrobust::rdrobust(y,x)
cct2 <- rdrobust::rdrobust(y,x,covs=w)
##
Omin          <- as.integer(0.03*n)
Opct          <- 0.9
## Owidth        <- quantile(abs(x),0.05)
m             <- 10
num_cutpoints <- n
Nmin          <- 10
num_sweeps    <- 150
burnin        <- 50
p_categorical <- 2
h <- quantile(abs(x),seq(0.01,1,0.01))
### Parallelization
no_cores <- detectCores() - 1
registerDoParallel(no_cores)
fit1 <- fit.general(h,y,w,x)
stopImplicitCluster()
### Plot estimates for different Owidth
pred <- sapply(fit1,function(i) colMeans(i$tau.adj)[(burnin+1):num_sweeps])
tau.hat <- apply(pred,2,function(i) c(mean(i),quantile(i,c(0.025,0.975))))
tau.hat <- t(tau.hat)
###
png("Figures/cholera_owidth.png")
matplot(h[1:70],tau.hat[1:70,],lty=c(1,2,2),pch=19,type=c("b","l","l"),col="blue",xlab="h",ylab=expression(tau))
dev.off()
png("Figures/cholera_good_owidth.png")
matplot(h[10:20],tau.hat[10:20,],lty=c(1,2,2),pch=19,type=c("b","l","l"),col="blue",xlab="h",ylab=expression(tau))
dev.off()
png("Figures/cholera_post_variance.png")
plot(h,apply(pred,2,var),type="l",ylab="Posterior variance")
abline(v=h[10],lty=2)
abline(v=h[20],lty=2)
dev.off()
### Obtain estimate for Owidth 12.5%
fit2 <- fit.xbcf(y,NULL,x)
ate.post1 <- fit2$ate.post[(burnin+1):num_sweeps]
fit3 <- fit.xbcf(y,w,x)
ate.post2 <- fit3$ate.post[(burnin+1):num_sweeps]
### Plot results
p <- data.frame(Method = c("CCT1","CCT2","XBCF1","XBCF2"),
                Est = c(cct1$coef[1],cct2$coef[1],
                        mean(ate.post1),mean(ate.post2)),
                LI = c(cct1$ci[1,1],cct2$ci[1,1],
                       quantile(ate.post1,0.025),
                       quantile(ate.post2,0.025)),
                UI = c(cct1$ci[1,2],cct2$ci[1,2],
                       quantile(ate.post1,0.975),
                       quantile(ate.post2,0.975)),
                Position = c(1,3,5,7))
p1 <- reshape(p,direction="long",idvar=c("Method","Position"),varying=c("Est","LI","UI"),times=c("Est","LI","UI"),v.name="value")
cols <- c("firebrick","dodgerblue")
plot(p1$Position,p1$value,"n",xlab="",ylab="",axes=F)
axis(2)
axis(1,at=c(1,2,6,8),labels=c("","CCT","XBCF",""))
segments(x0=p$Position,y0=p$LI,y1=p$UI,col=cols,lwd=1.5)
segments(x0=p$Position-0.05,x1=p$Position+0.05,y0=p$LI,col=cols,lwd=1.5)
segments(x0=p$Position-0.05,x1=p$Position+0.05,y0=p$UI,col=cols,lwd=1.5)
points(x=p$Position,y=p$Est,col=cols,pch=19)
legend("topleft",legend=c("Yes","No"),cex=0.7,
       col=c(cols[2],cols[1]),lty=1,lwd=1.5,pch=19,
       title="Includes covariates",horiz=T)
