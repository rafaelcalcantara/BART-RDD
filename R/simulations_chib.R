## Install and load CGS package
install.packages("http://apps.olin.wustl.edu/faculty/chib/rpackages/rdd/bayesrdd_1.0.zip",
                 repo=NULL,source=T)
library(bayesrdd)
set.seed(000)
## Setup
mu <- function(W, X) return(0.1 * rowSums(W) + 1/(1+exp(-5*X)))
tau <- function(W, X) return( sin(mu(W, X)) + 1)
n <- 100
p <- 2
c <- 0 # Cutoff
P = c(.70,.30);
mz = c(3, 3);
mztau = c(3, 2);
lamstmean0_ = c(1,1);
lamstsd0_ = 5*c(1,1);
s2mean0_ = .3;
s2sd0_ = 1;
nudgp = 3;
nuval = 100;
burn <- 1000
nsamples <- 100
## Simulations
s <- 100 ## samples
###
results <- list("ATE"=vector("numeric",s),"XBCF-RDD (1)"=NA,"XBCF-RDD (2)"=NA,
                "XBCF-RDD (3)"=NA,"CGS"=matrix(0,s,3),"KR"=NA,"FH"=NA)
###
for (i in 1:s)
{
    print(paste0("Iteration ",i))
    w <- matrix(rnorm(n*p), n, p)
    x <- rnorm(n,sd=.5)
    z <- x >= c
    y <- mu(w, x) + tau(w, x)*z + rnorm(n, 0, 0.1)
    results[["ATE"]][i] <- mean(tau(w,0))
    ## Estimation
    ate.cgs <- invisible(bayesrddest(y = y,
                                     z = x,
                                     tau = c,
                                     p = P,
                                     mz = mz,
                                     mztau = mztau,
                                     beta0_ = rep(0,4),
                                     lamstmean0_ = lamstmean0_,
                                     lamstsd0_ = lamstsd0_,
                                     d = c(1,1),
                                     s2mean0_ = s2mean0_,
                                     s2sd0_ = s2sd0_,
                                     distribution = "gaussian",
                                     nuval = nuval,
                                     hetero = FALSE,
                                     n0=burn,
                                     m=nsamples))
    ate.cgs <- c(mean(ate.cgs$atem),quantile(ate.cgs$atem,c(.025,.975)))
    ## Store results
    results$CGS[i,] <- ate.cgs
}
##
saveRDS(results,"~/Git/XBCF-RDD/R/results_cgs.rds")
