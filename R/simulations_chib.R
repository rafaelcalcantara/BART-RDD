## Install and load CGS package
install.packages("http://apps.olin.wustl.edu/faculty/chib/rpackages/rdd/bayesrdd_1.0.zip",
                 repo=NULL,source=T)
library(bayesrdd)
set.seed(000)
## Setup
mu <- function(W, X) return(0.1 * rowSums(W) + 1/(1+exp(-5*X)))
tau <- function(W, X) return( sin(mu(W, X)) + 1)
n <- 2000
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
results <- list(mse=NA,cont.tau=NA,cont.zero=NA,int.length=NA,pe=NA)
dims <- list(NULL,c("XBCF-RDD (1)","XBCF-RDD (2)","XBCF-RDD (3)",
                    "CGS","KR","FH"))
results$mse <- matrix(0,nrow=s,ncol=6,dimnames=dims)
results$cont.tau <- matrix(0,nrow=s,ncol=6,dimnames=dims)
results$cont.zero <- matrix(0,nrow=s,ncol=6,dimnames=dims)
results$int.length <- matrix(0,nrow=s,ncol=6,dimnames=dims)
results$pe <- matrix(0,nrow=s,ncol=6,dimnames=dims)
###
for (i in 1:s)
{
  w <- matrix(rnorm(n*p), n, p)
  x <- rnorm(n,sd=.5)
  z <- x >= c
  y <- mu(w, x) + tau(w, x)*z + rnorm(n, 0, 0.1)
  true.ate <- mean(tau(w,0))
  ## Estimation
  ate.cgs <- bayesrddest(y = y,
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
                         m=nsamples)
  ate.cgs <- c(quantile(ate.cgs$atem,c(.025,.975)),mean(ate.cgs$atem))
  ## Store results
  ### mse
  results$mse[i,"CGS"] <- (ate.cgs[3]-true.ate)^2
  ### cont.tau
  results$cont.tau[i,"CGS"] <- true.ate >= ate.cgs[1] & true.ate <= ate.cgs[2]
  ### cont.zero
  results$cont.zero[i,"CGS"] <- 0 >= ate.cgs[1] & 0 <= ate.cgs[2]
  ### int.length
  results$int.length[i,"CGS"] <- ate.cgs[2] - ate.cgs[1]
### pe (point estimate)
  results$pe[i,"CGS"] <- ate.cgs[3]
}
##
saveRDS(results,"~/Git/XBCF-RDD/R/results_cgs.rds")
