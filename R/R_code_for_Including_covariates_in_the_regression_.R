library(np)  
library(rdd)
library(matrixStats)
library(xtable)
library(boot)
##############
                                        # Procedures #
##############

                                        # rdd.x.sim: RDD procedure for the simulation
                                        # input:
                                        # y: outcome
                                        # z: running variable
                                        # x: covariates (control variables)
                                        # t: treatment - only needs to be specified in the fuzzy RDD (if t is set to NULL, the sharp RDD is assumed)
                                        # bw0, bw1: bandwidths for the estimation of E(Y|X,Z) below and above the cut-off
                                        # bwz: bandwidth for the kernel function on Z
                                        # bw0t, bw1t: bandwidths for the estimation of E(T|X,Z) below and above the cut-off

                                        # output:
                                        # effect: effect based on CV bandwidths
                                        # effectu: effect based on undersmoothing (CV/2) 
                                        # effecto: effect based on oversmoothing (CV*2)
                                        # neffect, neffectu, neffecto: "naive effects" without boundary kernel

rdd.x.sim=function(y,z,x, t=NULL, bw0=NULL, bw1=NULL, bwz=NULL, bw0t=NULL, bw1t=NULL){
    ## NOTE: we don't consider under- or over-smoothing (regXo and regXu, i.e. CV/2 and 2*CV)
    ## NOTE: we don't consider the naive estimator
    ## Setting
    d=1*(z>=0)
    xz=data.frame(x,z)
    xzcutoff=data.frame(x,rep(0,length(d)))
    xz0=xz[d==0,]; xz1=xz[d==1,]; d1=d[d==1]; d0=d[d==0]; y1=y[d==1]; y0=y[d==0];
    ## bw0
    if (is.null(bw0))
    {
        regbw0 <- npregbw(ydat=y0, xdat=xz0, bwmethod="cv.ls", ckertype="epanechnikov", regtype="ll")
    } else {
        regbw0 <- bw0
    }
    ## reg0
    reg0<-npreg(bws=regbw0, tydat=y0, txdat=xz0, exdat=xzcutoff, bwmethod="cv.ls", ckertype="epanechnikov", regtype="ll")$mean
    ### reg0u<-npreg(bws=regbw0$bw/2, tydat=y0, txdat=xz0, exdat=xzcutoff, bwmethod="cv.ls", ckertype="epanechnikov", regtype="ll")$mean
    ### reg0o<-npreg(bws=regbw0$bw*2, tydat=y0, txdat=xz0, exdat=xzcutoff, bwmethod="cv.ls", ckertype="epanechnikov", regtype="ll")$mean
    ## bw1
    if (is.null(bw1))
    {
        regbw1 <- npregbw(ydat=y1, xdat=xz1, bwmethod="cv.ls", ckertype="epanechnikov", regtype="ll")
    } else {
        regbw1 <- bw1
    }
    ## reg1
    reg1<-npreg(bws=regbw1, tydat=y1, txdat=xz1, exdat=xzcutoff, bwmethod="cv.ls", ckertype="epanechnikov", regtype="ll")$mean
    ### reg1u<-npreg(bws=regbw1$bw/2, tydat=y1, txdat=xz1, exdat=xzcutoff, bwmethod="cv.ls", ckertype="epanechnikov", regtype="ll")$mean
### reg1o<-npreg(bws=regbw1$bw*2, tydat=y1, txdat=xz1, exdat=xzcutoff, bwmethod="cv.ls", ckertype="epanechnikov", regtype="ll")$mean
    ## bwz
    if (is.null(bwz))
    {
        regbw <- npregbw(ydat=y, xdat=z, bwmethod="cv.ls", ckertype="epanechnikov", regtype="lc")
    } else {
        regbw <- bwz
    }
    ## kernwgt
    kernwgt=npksum(bws=regbw, tydat=y, txdat=z, exdat=0, ckertype="epanechnikov", regtype="lc", return.kernel.weights=TRUE )$kw
    ### kernwgtu=npksum(bws=regbw$bw/2, tydat=y, txdat=z, exdat=0, ckertype="epanechnikov", regtype="lc", return.kernel.weights=TRUE )$kw
    ### kernwgto=npksum(bws=regbw$bw*2, tydat=y, txdat=z, exdat=0, ckertype="epanechnikov", regtype="lc", return.kernel.weights=TRUE )$kw
    if (is.null(t)==1){
        ## neffect=(sum((reg1-reg0)*kernwgt))/(sum(kernwgt))
        ## neffectu=(sum((reg1u-reg0u)*kernwgtu))/(sum(kernwgtu))
        ## neffecto=(sum((reg1o-reg0o)*kernwgto))/(sum(kernwgto))
        mu2=0.1; mu1=3/16
        kernwgt=(mu2-mu1*z)*kernwgt
        kernwgtu=(mu2-mu1*z)*kernwgtu
        kernwgto=(mu2-mu1*z)*kernwgto
        effect=(sum((reg1-reg0)*kernwgt))/(sum(kernwgt))
        ## effectu=(sum((reg1u-reg0u)*kernwgtu))/(sum(kernwgtu))
        ## effecto=(sum((reg1o-reg0o)*kernwgto))/(sum(kernwgto))
    }
    ## Comment out the region below because it's for fuzzy RDD only
    ## if (is.null(t)!=1){
    ##     y1=t[d==1]; y0=t[d==0];
    ##     if (is.null(bw0t)==1) regbw0t <- npregbw(ydat=y0, xdat=xz0, bwmethod="cv.ls", ckertype="epanechnikov", regtype="ll")
    ##     if (is.null(bw0t)==0) regbw0t <- bw0t
    ##     reg0t<-npreg(bws=regbw0t, tydat=y0, txdat=xz0, exdat=xzcutoff, bwmethod="cv.ls", ckertype="epanechnikov", regtype="ll")$mean 
    ##     reg0ut<-npreg(bws=regbw0t$bw/2, tydat=y0, txdat=xz0, exdat=xzcutoff, bwmethod="cv.ls", ckertype="epanechnikov", regtype="ll")$mean 
    ##     reg0ot<-npreg(bws=regbw0t$bw*2, tydat=y0, txdat=xz0, exdat=xzcutoff, bwmethod="cv.ls", ckertype="epanechnikov", regtype="ll")$mean
    ##     if (is.null(bw1t)==1) regbw1t <- npregbw(ydat=y1, xdat=xz1, bwmethod="cv.ls", ckertype="epanechnikov", regtype="ll")
    ##     if (is.null(bw1t)==0) regbw1t <- bw1t
    ##     reg1t<-npreg(bws=regbw1t, tydat=y1, txdat=xz1, exdat=xzcutoff, bwmethod="cv.ls", ckertype="epanechnikov", regtype="ll")$mean
    ##     reg1ut<-npreg(bws=regbw1t$bw/2, tydat=y1, txdat=xz1, exdat=xzcutoff, bwmethod="cv.ls", ckertype="epanechnikov", regtype="ll")$mean
    ##     reg1ot<-npreg(bws=regbw1t$bw*2, tydat=y1, txdat=xz1, exdat=xzcutoff, bwmethod="cv.ls", ckertype="epanechnikov", regtype="ll")$mean
    ##     neffect=(sum((reg1-reg0)*kernwgt))/(sum((reg1t-reg0t)*kernwgt)); neffectu=(sum((reg1u-reg0u)*kernwgtu))/(sum((reg1ut-reg0ut)*kernwgtu))
    ##     neffecto=(sum((reg1o-reg0o)*kernwgto))/(sum((reg1ot-reg0ot)*kernwgto))
    ##     mu2=0.1; mu1=3/16
    ##     kernwgt=(mu2-mu1*z)*kernwgt; kernwgtu=(mu2-mu1*z)*kernwgtu
    ##     kernwgto=(mu2-mu1*z)*kernwgto
    ##     effect=(sum((reg1-reg0)*kernwgt))/(sum((reg1t-reg0t)*kernwgt)); effectu=(sum((reg1u-reg0u)*kernwgtu))/(sum((reg1ut-reg0ut)*kernwgtu))
    ##     effecto=(sum((reg1o-reg0o)*kernwgto))/(sum((reg1ot-reg0ot)*kernwgto))
    ## }
    output=c(effect, effectu, effecto, neffect, neffectu, neffecto)
    output
}


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

rdd.x=function(y,z,x, boot=1999, bw0=NULL, bw1=NULL, bwz=NULL){
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
    est=rdd.x.est(y,z,x, bw0, bw1, bwz)
    se=sqrt(colVars(rdd.x.boot(y=y,z=z,x=x, bw0=bw0, bw1=bw1, bwz=bwz, boot=boot)))
    results=rbind(est,se, 2*pnorm(-abs(est/se)))
    colnames(results)=c("bw/2", "bw", "bw*2")
    rownames(results)=c("est", "se", "p-val")
    results
}

rdd.x.est=function(y,z,x, bw0, bw1, bwz){
    d=1*(z>=0)
    xz=data.frame(x,z)
    xzcutoff=data.frame(x,rep(0,length(d)))
    xz0=xz[d==0,]; xz1=xz[d==1,]; d1=d[d==1]; d0=d[d==0]; y1=y[d==1]; y0=y[d==0]; 
    reg0<-npreg(bws=bw0, tydat=y0, txdat=xz0, exdat=xzcutoff, bwmethod="cv.ls", ckertype="epanechnikov", regtype="ll")$mean 
    reg0u<-reg0; reg0o<-reg0
    reg1<-npreg(bws=bw1, tydat=y1, txdat=xz1, exdat=xzcutoff, bwmethod="cv.ls", ckertype="epanechnikov", regtype="ll")$mean
    reg1u<-reg1; reg1o<-reg1 
    kernwgt=npksum(bws=bwz, tydat=y, txdat=z, exdat=0, ckertype="epanechnikov", regtype="lc", return.kernel.weights=TRUE )$kw
    kernwgtu=npksum(bws=bwz/2, tydat=y, txdat=z, exdat=0, ckertype="epanechnikov", regtype="lc", return.kernel.weights=TRUE )$kw
    kernwgto=npksum(bws=bwz*2, tydat=y, txdat=z, exdat=0, ckertype="epanechnikov", regtype="lc", return.kernel.weights=TRUE )$kw
    mu2=0.1; mu1=3/16
    kernwgt=(mu2-mu1*z)*kernwgt; kernwgtu=(mu2-mu1*z)*kernwgtu
    kernwgto=(mu2-mu1*z)*kernwgto
    effect=(sum((reg1-reg0)*kernwgt))/(sum(kernwgt)); effectu=(sum((reg1u-reg0u)*kernwgtu))/(sum(kernwgtu))
    effecto=(sum((reg1o-reg0o)*kernwgto))/(sum(kernwgto))
    output=c(effectu, effect, effecto) 
    output
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

bs.rdd2 = function(data,indices) {
    dat = data[indices,] 
    effect=RDestimate(dat[,1]~dat[,2] , kernel = "epanechnikov", bw=bwz)$est
    return(effect)
} 


## ##############
##                                         # Simulation #
## ##############

##                                         #simulation parameters
## n=1000             # sample size
## repetitions=1000   # number of simulations
## alpha=0            # association between X and treatment (assignment)
## beta=0.4           # effect of X on Y 
## fuzzy=0            # fuzzy RDD if equal to 1 


##                                         #run simulation
## i=1
## results=c(); compliers=c()
## while(i<=repetitions){
##     set.seed(i)
##     z=rnorm(n)
##     u=rnorm(n)
##     d=z>=0
##     x1<-0.5*rnorm(n)+alpha*d
##     x2<-0.5*rnorm(n)+alpha*d
##     x=data.frame(c(x1), c(x2))
##     if (fuzzy!=1) y=d+0.5*z-0.25*d*z+0.25*z^2+beta*(x1+x2)+0.5*beta*(x1^2+x2^2)+u
##     if (fuzzy==1){
##         t=(2*d+0.5*u+rnorm(n)>1)*1
##         comp=((2+0.5*u+rnorm(n)>1) & (0.5*u+rnorm(n)<=1) ) *1
##         y=t+0.5*z-0.25*t*z+0.25*z^2+beta*(x1+x2)+0.5*beta*(x1^2+x2^2)+u
##     } 
##     xz=data.frame(x,z)
##     xzcutoff=data.frame(x,rep(0,length(d)))
##     xz0=xz[d==0,]; xz1=xz[d==1,]; d1=d[d==1]; d0=d[d==0]; y1=y[d==1]; y0=y[d==0];
##     bw0 <- npregbw(ydat=y0, xdat=xz0, bwmethod="cv.ls", ckertype="epanechnikov", regtype="ll")
##     bw1 <- npregbw(ydat=y1, xdat=xz1, bwmethod="cv.ls", ckertype="epanechnikov", regtype="ll")
##     bwz <- npregbw(ydat=y, xdat=z, bwmethod="cv.ls", ckertype="epanechnikov", regtype="lc")
##     if (fuzzy!=1) temp=c(rdd.x.sim(y,z,x, bw0=bw0, bw1=bw1, bwz=bwz), RDestimate(y~z, kernel = "epanechnikov", bw = bwz$bw)$est, RDestimate(y~z, kernel = "epanechnikov")$est)
##     if (fuzzy==1) temp=c(rdd.x.sim(y,z,x, t=t, bw0=bw0, bw1=bw1, bwz=bwz), RDestimate(y~z+t, kernel = "epanechnikov", bw = bwz$bw)$est, RDestimate(y~z+t, kernel = "epanechnikov")$est)
##     results=rbind(results,temp)
##     if (fuzzy==1) compliers=c(compliers,comp)
##     i=i+1
## }    

##                                         #compute bias, standard deviation, and RMSE
## bias=colMeans(results)-1
## sdev=sqrt(colVars(results))
## rmse=sqrt(bias^2+sdev^2)

##                                         # show simulation results
## results=rbind(bias, sdev, rmse)
## xtable(results[,c(1:3,7:12)], round=3)


## ###############
##                                         # Application #
## ###############

##                                         #load data and define variables
## load("C:/Lalive2008.RData")
## attach(data)
## y=unemployment_duration[female == 1 & period == 1 & tr==1]
## z=dage_1[female == 1 & period == 1 & tr==1]
## marrstatus=1*marr+2*single 
## educ=1*educ_med+2*educ_hi 
## sector=1*landw+2*versorg+3*nahrung+4*textil+5*holzind+6*elmasch+7*andfabr+8*bau+9*gasthand+10*verkehr+11*dienstl
## x=data.frame(factor(marrstatus), ordered(educ), factor(foreign), c(rr), c(lwage_ljob), c(previous_experience), factor(white_collar), factor(sector))[female == 1 & period == 1 & tr==1,]
## xpara=cbind(marr, single, educ_med, educ_hi, foreign, rr, lwage_ljob, previous_experience, white_collar, landw, versorg, nahrung, textil, holzind, elmasch, andfabr, bau, gasthand, verkehr, dienstl)[female == 1 & period == 1 & tr==1,]
## d=1*(z>=0)
## xz=data.frame(x,z)
## xzcutoff=data.frame(x,rep(0,length(d)))
## xz0=xz[d==0,]; xz1=xz[d==1,]; d1=d[d==1]; d0=d[d==0]; y1=y[d==1]; y0=y[d==0];

##                                         # test covariate balance at cut-off
## testout=c()
## for (i in 1:ncol(xpara)){
##     out=RDestimate(xpara[,i]~z , kernel = "epanechnikov")
##     temp=c(mean(xpara[,i]),out$est[1], out$p[1], out$est[2], out$p[2])
##     testout=rbind(testout,temp)
## }
## rownames(testout)=c("married (binary)", "single (binary)", "education: medium (binary)", "education: high (binary)", "foreign (binary)", "replacement rate", "log wage in last job", "actual to potential work experience", "white collar (binary)", "sector: agriculture (binary)", "sector: utilities (binary)", "sector: food (binary)", "sector: textiles (binary)", "sector: wood (binary)", "sector: machines (binary)", "sector: other manufactoring (binary)", "sector: construction (binary)", "sector: tourism (binary)", "sector: traffic (binary)", "sector: services (binary)")
## xtable(testout, round=3)

##                                         # define bandwidth for regressions below and above the cutoff, as well as Z
## regbw0 <- npregbw(ydat=y0, xdat=xz0, bwmethod="cv.ls", ckertype="epanechnikov", regtype="ll")$bw
## regbw1 <- npregbw(ydat=y1, xdat=xz1, bwmethod="cv.ls", ckertype="epanechnikov", regtype="ll")$bw
## bwz=0.2

##                                         # run RDD with X
## results=rdd.x(y,z,x, bw0=regbw0, bw1=regbw1, bwz=bwz, boot=999)
## xtable(results, round=3)

##                                         # run RDD without X
## temp = boot(data=data.frame(y,z), statistic=bs.rdd2, R=999)
## results2=rbind(round(results3$t0,3),round(colSds(results3$t),3), round(2*pnorm(-abs(results3$t0/colSds(results3$t))),3))
## rownames(results2)=c("est","se","p-val")
## xtable(results2)
