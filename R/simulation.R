## Setup
library(parallel)
library(foreach)
library(doParallel)
s      <- 200
sample <- c(500,2000)
model  <- c(2,4)
xi <- c(0.05,0.2)
nu <- kappa <- c(0.5,2)
x.dist <- 1
## Parallelization
no_cores <- detectCores() - 1
registerDoParallel(no_cores)
if (Sys.info()["sysname"] == "Linux")
{
    ## Setup
    devtools::install_github("JingyuHe/XBART@XBCF-RDD")
    library(XBART)
    library(rdrobust)
    devtools::install_github("kolesarm/RDHonest")
    library(RDHonest)
    devtools::install_github("akreiss/HighDimRD")
    library(HighDimRD)
    ## BART Prior settings
    c             <- 0
    Omin          <- 1
    Opct          <- 0.9
    ntrees        <- 10
    Nmin          <- 10
    num_sweeps    <- 150
    burnin        <- 50
    p_categorical <- 3
    ## BART-RDD
    fit <- function(s,ncut,Model,xi,nu,kappa,x.dist)
    {
        foreach(i=1:s,.multicombine=T,.export=c("c","Omin","Opct","m","num_sweeps","burnin","Nmin","p_categorical")) %dopar%
            {
                print("BART-RDD simulations")
                print(paste0("N=",ncut,"; Model=",Model,"; xi=",xi,"; nu=",nu,"; kappa=",kappa,"; X distribution=",x.dist))
                print(paste0("Simulation ",i))
                data <- readRDS(paste0("Data/DGP_",ncut,"_",Model,"_",xi,"_",nu,"_",kappa,"_",x.dist,"_",i,".rds"))
                h <- quantile(data$X,seq(0.1,0.9,0.05))
                h <- abs(h)
                h <- sort(h)[2:4]
                dt <- 1:3
                pred <- vector("list",3)
                for (j in 1:3)
                {
                    test <- -h[j]<=data$X & data$X<=h[j]
                    print(paste0("Test set size: ",sum(test)))
                    if (sum(test)>10)
                    {
                        t0 <- Sys.time()
                        fit.obj <- XBCF.rd(data$Y, data$W, data$X, c,
                                           Owidth = h[j],
                                           Omin = Omin, Opct = Opct,
                                           num_trees_mod = ntrees,
                                           num_trees_con = ntrees,
                                           num_cutpoints = ncut,
                                           num_sweeps = num_sweeps,
                                           burnin = burnin, Nmin = Nmin,
                                           p_categorical_con = p_categorical,
                                           p_categorical_mod = p_categorical,
                                           tau_con = 2*var(data$Y)/ntrees,
                                           tau_mod = 0.5*var(data$Y)/ntrees,
                                           parallel=F)
                        pred[[j]] <- predict.XBCFrd(fit.obj,data$W[test,],rep(0,sum(test)))
                        pred[[j]] <- pred[[j]]$tau.adj[,(burnin+1):num_sweeps]
                        t1 <- Sys.time()
                        dt[j] <- difftime(t1,t0)
                    } else
                    {
                        print("Identification strip has less than 10 points; dropping current h")
                        h[j] <- NA
                    }
                }
                pred <- pred[complete.cases(h)]
                print("Candidate estimates:")
                print(sapply(pred,function(a) mean(colMeans(a))))
                post.var <- sapply(pred, function(a) var(colMeans(a)))
                pred <- pred[[which(post.var==max(post.var))]]
                dt <- dt[which(post.var==max(post.var))]
                print(paste0("Chosen h: ",h[which(post.var==max(post.var))]))
                print(paste0("True tau: ",xi))
                print(paste0("Tau estimate: ", mean(colMeans(pred))))
                return(list(pred=pred,dt=dt))
            }
    }
###
    for (i in sample)
    {
        for (j in model)
        {
            for (k in xi)
            {
                for (l in nu)
                {
                    for (m in kappa)
                    {
                        for (n in x.dist)
                        {
                            bart.rdd <- fit(s,i,j,k,l,m,n)
                            saveRDS(bart.rdd,paste0("Results/bart_rdd_",i,"_",j,"_",k,"_",l,"_",m,"_",n,".rds"))
                        }
                    }
                }
            }
        }
    }
        ## CCT
        fit <- function(s,n,m,xi,nu,kappa,x.dist)
        {
            foreach(i=1:s,.multicombine=T,.export=c("n","c","Omin","h","Opct","m","n","num_sweeps","burnin","Nmin","p_categorical")) %dopar%
                {
                    print("CCT simulations")
                    print(paste0("N=",n,"; Model=",m,"; xi=",xi,"; nu=",nu,"; kappa=",kappa,"; X distribution=",x.dist))
                    print(paste0("Simulation ",i))
                    data <- readRDS(paste0("Data/DGP_",n,"_",m,"_",xi,"_",nu,"_",kappa,"_",x.dist,"_",i,".rds"))
                    t0 <- Sys.time()
                    fit <- rdrobust(data$Y,data$X,0,covs=data$W)
                    t1 <- Sys.time()
                    dt <- difftime(t1,t0)
                    pred <- c(fit$coef[1],fit$ci[3,])
                    return(list(pred=pred,dt=dt))
                }
        }
    ###
        for (i in sample)
        {
            for (j in model)
            {
                for (k in xi)
                {
                    for (l in nu)
                    {
                        for (m in kappa)
                        {
                            for (n in x.dist)
                                {
                                    cct <- fit(s,i,j,k,l,m,n)
                                    saveRDS(cct,paste0("Results/cct_",i,"_",j,"_",k,"_",l,"_",m,"_",n,".rds"))
                                }
                        }
                    }
                }
            }
        }
        ## KR
        fit <- function(s,n,m,xi,nu,kappa,x.dist)
        {
            foreach(i=1:s,.multicombine=T,.export=c("n","c","Omin","h","Opct","m","n","num_sweeps","burnin","Nmin","p_categorical")) %dopar%
                {
                    print("KR simulations")
                    print(paste0("N=",n,"; Model=",m,"; xi=",xi,"; nu=",nu,"; kappa=",kappa,"; X distribution=",x.dist))
                    print(paste0("Simulation ",i))
                    data <- readRDS(paste0("Data/DGP_",n,"_",m,"_",xi,"_",nu,"_",kappa,"_",x.dist,"_",i,".rds"))
                    t0 <- Sys.time()
                    w1 <- fourier_basis(matrix(data$X),4)
                    w_HighDim <- cbind(data$W,interaction_terms(data$W),w1,interaction_terms(w1))
                    ## Estimation
                    fit  <- HighDim_rd(data$Y,data$X,w_HighDim,
                                       tpc="CV",rd="robust")
                    t1 <- Sys.time()
                    dt <- difftime(t1,t0)
                    pred <- c(fit$rd[["Estimate"]][,"tau.bc"],fit$rd[["ci"]]["Robust",])
                    return(list(pred=pred,dt=dt))
                }
        }
    ###
        for (i in sample)
        {
            for (j in model)
            {
                for (k in xi)
                {
                    for (l in nu)
                    {
                        for (m in kappa)
                        {
                            for (n in x.dist)
                                {
                                    kr <- fit(s,i,j,k,l,m,n)
                                    saveRDS(kr,paste0("Results/kr_",i,"_",j,"_",k,"_",l,"_",m,"_",n,".rds"))
                                }
                        }
                    }
                }
            }
        }
    ## BART1
    fit <- function(s,n,Model,xi,nu,kappa,x.dist)
    {
        foreach(i=1:s,.multicombine=T,.export=c("n","c","Omin","h","Opct","m","n","num_sweeps","burnin","Nmin","p_categorical")) %dopar%
            {
                print("BART1 simulations")
                print(paste0("N=",n,"; Model=",Model,"; xi=",xi,"; nu=",nu,"; kappa=",kappa,"; X distribution=",x.dist))
                print(paste0("Simulation ",i))
                data <- readRDS(paste0("Data/DGP_",n,"_",Model,"_",xi,"_",nu,"_",kappa,"_",x.dist,"_",i,".rds"))
                X <- cbind(data$X,data$W,data$Z)
                t0 <- Sys.time()
                fit <- XBART(data$Y, X, num_trees = m,
                             num_cutpoints = n, num_sweeps = num_sweeps,
                             burnin = burnin, Nmin = Nmin,
                             p_categorical_con = p_categorical+1,
                             ## add 1 for Z
                             p_categorical_mod = p_categorical+1,
                             tau_con = 2*var(y)/m,
                             tau_mod = 0.5*var(y)/m, parallel=F)
                test <- -h<=data$X & data$X<=h
                pred1 <- predict.XBART(fit,cbind(rep(0,sum(test)),data$W[test,],rep(1,sum(test))))[,(burnin+1):num_sweeps]
                pred0 <- predict.XBART(fit,cbind(rep(0,sum(test)),data$W[test,],rep(0,sum(test))))[,(burnin+1):num_sweeps]
                pred <- pred1-pred0
                t1 <- Sys.time()
                dt <- difftime(t1,t0)
                return(list(pred=pred,dt=dt))
            }
    }
###
    for (i in sample)
    {
        for (j in model)
        {
            for (k in xi)
            {
                for (l in nu)
                {
                    for (m in kappa)
                    {
                        for (n in x.dist)
                        {
                            bart1 <- fit(s,i,j,k,l,m,n)
                            saveRDS(bart1,paste0("Results/bart1_",i,"_",j,"_",k,"_",l,"_",m,"_",n,".rds"))
                        }
                    }
                }
            }
        }
    }
    ## BART2
    fit <- function(s,n,Model,xi,nu,kappa,x.dist)
    {
        foreach(i=1:s,.multicombine=T,.export=c("n","c","Omin","h","Opct","m","n","num_sweeps","burnin","Nmin","p_categorical")) %dopar%
            {
                print("BART2 simulations")
                print(paste0("N=",n,"; Model=",Model,"; xi=",xi,"; nu=",nu,"; kappa=",kappa,"; X distribution=",x.dist))
                print(paste0("Simulation ",i))
                data <- readRDS(paste0("Data/DGP_",n,"_",Model,"_",xi,"_",nu,"_",kappa,"_",x.dist,"_",i,".rds"))
                X <- cbind(data$X,data$W)
                t0 <- Sys.time()
                fit0 <- XBART(data$Y[!data$Z], X[!data$Z,],
                              num_trees = m,
                              num_cutpoints = n,
                              num_sweeps = num_sweeps,
                              burnin = burnin, Nmin = Nmin,
                              p_categorical_con = p_categorical,
                              p_categorical_mod = p_categorical,
                              tau_con = 2*var(data$Y[!data$Z])/m,
                              tau_mod = 0.5*var(data$Y[!data$Z])/m,
                              parallel=F)
                fit1 <- XBART(data$Y[data$Z], X[data$Z,],
                              num_trees = m,
                              num_cutpoints = n,
                              num_sweeps = num_sweeps,
                              burnin = burnin, Nmin = Nmin,
                              p_categorical_con = p_categorical,
                              p_categorical_mod = p_categorical,
                              tau_con = 2*var(data$Y[data$Z])/m,
                              tau_mod = 0.5*var(data$Y[data$Z])/m,
                              parallel=F)
                test0 <- -h<=data$X
                test1 <- data$X<=h
                pred1 <- predict.XBART(fit1,cbind(rep(0,sum(data$Z & test1)),data$W[test1,]))
                pred0 <- predict.XBART(fit0,cbind(rep(0,sum(!data$Z & test0)),data$W[test0,]))
                pred <- colMeans(pred1[,(burnin+1):num_sweeps])-colMeans(pred0[,(burnin+1):num_sweeps])
                t1 <- Sys.time()
                dt <- difftime(t1,t0)
                return(list(pred=pred,dt=dt))
            }
    }
###
    for (i in sample)
    {
        for (j in model)
        {
            for (k in xi)
            {
                for (l in nu)
                {
                    for (m in kappa)
                    {
                        for (n in x.dist)
                            {
                                bart2 <- fit(s,i,j,k,l,m,n)
                                saveRDS(bart2,paste0("Results/bart2_",i,"_",j,"_",k,"_",l,"_",m,"_",n,".rds"))
                            }
                    }
                }
            }
        }
    }
    ## BCF
    fit <- function(s,n,Model,xi,nu,kappa,x.dist)
    {
        foreach(i=1:s,.multicombine=T,.export=c("n","c","Omin","h","Opct","m","n","num_sweeps","burnin","Nmin","p_categorical")) %dopar%
            {
                print("BCF simulations")
                print(paste0("N=",n,"; Model=",Model,"; xi=",xi,"; nu=",nu,"; kappa=",kappa,"; X distribution=",x.dist))
                print(paste0("Simulation ",i))
                data <- readRDS(paste0("Data/DGP_",n,"_",Model,"_",xi,"_",nu,"_",kappa,"_",x.dist,"_",i,".rds"))
                X <- cbind(data$X,data$W)
                t0 <- Sys.time()
                fit <- XBCF.discrete(data$Y, data$Z, X_con=X, X_mod=X,
                                     num_trees_con = m,
                                     num_trees_mod = m,
                                     num_cutpoints = n,
                                     num_sweeps = num_sweeps,
                                     burnin = burnin, Nmin = Nmin,
                                     p_categorical_con = p_categorical,
                                     p_categorical_mod = p_categorical,
                                     tau_con = 2*var(data$Y)/m,
                                     tau_mod = 0.5*var(data$Y)/m,
                                     parallel=F)
                test <- -h<=data$X & data$X<=h
                pred <- predict.XBCFdiscrete(fit,cbind(rep(0,sum(test)),data$W[test,]),
                                             cbind(rep(0,sum(test)),data$W[test,]),
                                             data$Z[test],pihat=data$Z[test],
                                             burnin=burnin)
                pred <- pred$tau.adj[,(burnin+1):num_sweeps]
                t1 <- Sys.time()
                dt <- difftime(t1,t0)
                return(list(pred=pred,dt=dt))
            }
    }
###
    for (i in sample)
    {
        for (j in model)
        {
            for (k in xi)
            {
                for (l in nu)
                {
                    for (m in kappa)
                    {
                        for (n in x.dist)
                            {
                                bcf <- fit(s,i,j,k,l,m,n)
                                saveRDS(bcf,paste0("Results/bcf_",i,"_",j,"_",k,"_",l,"_",m,"_",n,".rds"))
                            }
                    }
                }
            }
        }
    }
} else
{
    setwd("~/Git/XBCF-RDD/")
    ## CGS
    install.packages("http://apps.olin.wustl.edu/faculty/chib/rpackages/rdd/bayesrdd_1.0.zip",repo=NULL,source=T)
    library(bayesrdd)
    c <- 0
    burn <- 50
    nsamples <- 200
    fit <- function(s,p,n,m,xi,nu,kappa,x.dist)
    {
        foreach(i=1:s,.multicombine=T,.packages = "bayesrdd",.export=c("burn","nsamples")) %do%
            {
                print("CGS simulations")
                print(paste0("N=",n,"; Model=",m,"; xi=",xi,"; nu=",nu,"; kappa=",kappa,"; X distribution=",x.dist))
                print(paste0("Simulation ",i))
                data <- readRDS(paste0("Data/DGP_",n,"_",m,"_",xi,"_",nu,"_",kappa,"_",x.dist,"_",i,".rds"))
                P = c(.70,.30)
                mz = c(5,5)
                mztau = c(5,5)
                s2mean0_ = .3
                s2sd0_ = 1
                lamstmean0_ = rep(1,p+2)
                lamstsd0_ = 5*rep(1,p+2)
                t0 <- Sys.time()
                fit  <- bayesrddest(y = data$Y,
                                    z = data$X,
                                    V = as.matrix(data$W[,5]),
                                    W = data$W[,1:4],
                                    mw = rep(5,p-1), ## as per help file
                                    tau = c,
                                    p = P,
                                    mz = mz,
                                    mztau = mztau,
                                    beta0_ = rep(0,5+2*(p-1)),
                                    lamstmean0_ = lamstmean0_,
                                    lamstsd0_ = lamstsd0_,
                                    d = rep(1,p+2),
                                    s2mean0_ = s2mean0_,
                                    s2sd0_ = s2sd0_,
                                    distribution = "gaussian",
                                    nuval=3,
                                    hetero = FALSE,
                                    n0=burn,
                                    m=nsamples)
                t1 <- Sys.time()
                pred <- fit$atem
                dt <- difftime(t1,t0)
                return(list(pred=pred,dt=dt))
            }
    }
###
    for (i in sample)
    {
        for (j in model)
        {
            for (k in xi)
            {
                for (l in nu)
                {
                    for (m in kappa)
                    {
                        for (n in x.dist)
                            {
                                cgs <- fit(s,7,i,j,k,l,m,n)
                                saveRDS(cgs,paste0("Results/cgs_",i,"_",j,"_",k,"_",l,"_",m,"_",n,".rds"))
                            }
                    }
                }
            }
        }
    }
}
##
print("Simulation ended!")
