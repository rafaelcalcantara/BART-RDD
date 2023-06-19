## Setup
library(parallel)
library(foreach)
library(doParallel)
devtools::install_github("rafaelcalcantara/XBART@XBCF-RDD")
library(XBART)
library(rdrobust)
###
setwd("~/Documents/Git/XBCF-RDD/")
data <- foreign::read.dta("Data/CKT_2023_SIM--CostSharing.dta")
c <- 0
y <- data$op_n_per
x <- data$week
w <- data[c("sh_male", "h_inc_per", "sh_tpe", "b_year")]
n <- length(y)
## Plot data
l0 <- loess(y~x,data=data.frame(y=y[x<c],x=x[x<c]))
l0 <- predict(l0,se=T)
l0 <- data.frame(Est=l0$fit, LI=l0$fit-1.96*l0$se.fit,
                 UI=l0$fit+1.96*l0$se.fit)
l0 <- l0[order(x[x<c]),]
l1 <- loess(y~x,data=data.frame(y=y[x>=c],x=x[x>=c]))
l1 <- predict(l1,se=T)
l1 <- data.frame(Est=l1$fit, LI=l1$fit-1.96*l1$se.fit,
                 UI=l1$fit+1.96*l1$se.fit)
l1 <- l1[order(x[x>=c]),]
###
plot(x,y,pch=21,bg="azure",cex=0.5,
     xlab="Weeks from third birthday",
     ylab="Number of visits per 10.000")
matlines(x[x<c],l0,col="black",lty=1,lwd=c(1.5,1,1))
matlines(x[x>=c],l1,col="black",lty=1,lwd=c(1.5,1,1))
abline(v=c,lty=2)
## RDRobust estimation
cct1 <- rdrobust(y,x,c)
cct2 <- rdrobust(y,x,c,covs=w)
## XBCF estimation
fit.xbcf <- function(y,w,x)
{
    t0 <- Sys.time()
    h <- quantile(abs(x),0.125)
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
###
Omin          <- as.integer(0.03*n)
Opct          <- 0.9
m             <- 10
Nmin          <- 10
num_sweeps    <- 120
burnin        <- 20
p_categorical <- 1
num_cutpoints <- n
###
xbcf1 <- fit.xbcf(y,NULL,x)
ate.post1 <- xbcf1$ate.post
xbcf2 <- fit.xbcf(y,w,x)
ate.post2 <- xbcf2$ate.post
## Comparing estimates
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
