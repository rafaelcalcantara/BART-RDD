## Implement product kernel local regression of FH
### lambda (order of covariate specific kernels) must be adapted for different values of p
library(np)  
library(rdd)
library(matrixStats)
library(xtable)
library(boot)
##############
# Procedures #
# rdd.x: Sharp RDD procedure for the application
# input:
# y: outcome
# z: running variable
# x: covariates (control variables)
# boot: number of bootstrap replications
# bw0, bw1: bandwidths for the estimation of E(Y|X,Z) below and above the cut-off
# bwz: bandwidth for the kernel function on Z
# if no bandwidths are provided (NULL), cross-validation is used
# output:
# in rows: estimate (est), standard error (se), and p-value (p-val)
# in columns: estimates when using bandiwdths  bw0, bw1 and half of bandwidth bwz (bw/2),
# estimates when using bandiwdths  bw0, bw1, and bwz (bw),
# estimates when using bandiwdths  bw0, bw1 and twice bandwidth bwz (bw*2)

rdd.x=function(y,z,x, boot=200, bw0=NULL, bw1=NULL, bwz=NULL){
    d=1*(z>=0)
    xz=data.frame(x,z)
    xzcutoff=data.frame(x,rep(0,length(d)))
    xz0=xz[d==0,]; xz1=xz[d==1,]; d1=d[d==1]; d0=d[d==0]; y1=y[d==1]; y0=y[d==0];
    if (is.null(bw0)==1) regbw0 <- npregbw(ydat=y0, xdat=xz0, bwmethod="cv.ls", ckertype="epanechnikov", regtype="ll")$bw
    if (is.null(bw0)==0) regbw0 <- bw0
    if (is.null(bw1)==1) regbw1 <- npregbw(ydat=y1, xdat=xz1, bwmethod="cv.ls", ckertype="epanechnikov", regtype="ll")$bw
    if (is.null(bw1)==0) regbw1 <- bw1
    if (is.null(bwz)==1) regbw <- npregbw(ydat=y, xdat=z, bwmethod="cv.ls", ckertype="epanechnikov", regtype="lc")$bw
    if (is.null(bwz)==0) regbw <- bwz
    est=rdd.x.est(y,z,x, regbw0, regbw1, regbw)
    se=sqrt(colVars(rdd.x.boot(y=y,z=z,x=x, bw0=regbw0, bw1=regbw1, bwz=regbw, boot=boot)))
    return(c(est,est-1.96*se,est+1.96*se))
}

rdd.x.est=function(y,z,x, bw0, bw1, bwz){
    d=1*(z>=0)
    xz=data.frame(x,z)
    xzcutoff=data.frame(x,rep(0,length(d)))
    xz0=xz[d==0,]; xz1=xz[d==1,]; d1=d[d==1]; d0=d[d==0]; y1=y[d==1]; y0=y[d==0]; 
    reg0<-npreg(bws=bw0, tydat=y0, txdat=xz0, exdat=xzcutoff, bwmethod="cv.ls", ckertype="epanechnikov", regtype="ll")$mean 
    ## reg0u<-reg0; reg0o<-reg0
    reg1<-npreg(bws=bw1, tydat=y1, txdat=xz1, exdat=xzcutoff, bwmethod="cv.ls", ckertype="epanechnikov", regtype="ll")$mean
    ## reg1u<-reg1; reg1o<-reg1
    kernwgt=npksum(bws=bwz, tydat=y, txdat=z, exdat=0, ckertype="epanechnikov", regtype="lc", return.kernel.weights=TRUE )$kw
    ## kernwgtu=npksum(bws=bwz/2, tydat=y, txdat=z, exdat=0, ckertype="epanechnikov", regtype="lc", return.kernel.weights=TRUE )$kw
    ## kernwgto=npksum(bws=bwz*2, tydat=y, txdat=z, exdat=0, ckertype="epanechnikov", regtype="lc", return.kernel.weights=TRUE )$kw
    mu2=0.1; mu1=3/16
    kernwgt=(mu2-mu1*z)*kernwgt; ## kernwgtu=(mu2-mu1*z)*kernwgtu
    ## kernwgto=(mu2-mu1*z)*kernwgto
    effect=(sum((reg1-reg0)*kernwgt))/(sum(kernwgt)); ## effectu=(sum((reg1u-reg0u)*kernwgtu))/(sum(kernwgtu))
    ## effecto=(sum((reg1o-reg0o)*kernwgto))/(sum(kernwgto))
    ## output=c(effectu, effect, effecto) 
    ## output
    return(effect)
}



                                        # further functions (called by the main functions above)
rdd.x.boot<-function(y,z,x, bw0, bw1, bwz, boot=1999){
    obs<-length(y)
    mc=c(); i=1
    while(i<=boot){
        sboot<-sample(1:obs,obs,TRUE)
        yb=y[sboot]
        zb<-z[sboot]
        if (length(x)==length(y)) xb<-x[sboot]
        if (length(x)!=length(y)) xb<-x[sboot,]
        
        est<-c(rdd.x.est(yb,zb,xb, bw0, bw1, bwz))
        if (sum(is.na(est))==0) mc<-rbind(mc, est)
        i=i+1
    }
    mc
}
