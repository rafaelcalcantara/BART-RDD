set.seed(7)
setwd("~/Git/BART-RDD/")
library(XBART)
### Setup--------------------------------------------------
par(bty="L")
Owidth        <- 0.2
Omin          <- 5
Opct          <- 0.75
ntrees_con    <- 15
ntrees_mod    <- 5
ntrees        <- 5
Nmin          <- 5
max_depth     <- 250
num_sweeps    <- 250
burnin        <- 150
p_categorical <- 0
###
mu0.x <- function(x) 3*x^5 - 2.5*x^4 - 1.5*x^3 + 2*x^2 + 3*x + 2
mu0.w <- function(w,mu.x) w*sd(mu.x)/sd(w)
tau0.x <- function(x,c) (-exp(x)+exp(c))/(1+exp(2*x))
tau0.w <- function(w,type,ate) {
  if (type==1) out <- w < mean(w)
  if (type==2) out <- as.integer(cut(w,quantile(w,seq(0,1,length.out=4)),include.lowest = T))
  out <- out %% 2 == 0
  out <- (out-mean(out))/sd(out)
  return(out)
}
mu <- function(x,w,delta) {
  mu.x <- mu0.x(x)
  mu.w <- mu0.w(w,mu.x)
  (mu.x + mu.w)/sd(mu.x + mu.w)*delta
}
tau <- function(x,w,ate,c,type,delta) {
  tau.x <- tau0.x(x,c)
  tau.w <- tau0.w(w,type,ate)
  ate + (tau.x + tau.w)/sd(tau.x + tau.w)*delta
  # ate + tau.x
}
#### Parameters
ate <- 0.2
delta_mu <- 0.5
delta_tau <- 2
c <- 0
type <- 2
n <- 1000
#### Data
x <- 2*rbeta(n,2,4)-0.75
z <- as.numeric(x>=c)
w <- rnorm(n) + z + 10*x^3 - 3*x^2 - x
y <- mu(x,w,delta_mu) + tau(x,w,ate,c,type,delta_tau)*z
Ey <- y
y <- y + rnorm(n)
#####
par(mfrow=c(2,2))
plot(x,y,col=z+1,pch=19)
abline(v=c,lty=2)
plot(x,w)
### BART-RDD
fit <- XBART::XBCF.rd(y, w, x, c,
                      Owidth = Owidth, Omin = Omin, Opct = Opct,
                      num_trees_mod = ntrees_mod,
                      num_trees_con = ntrees_con,
                      num_cutpoints = n,
                      num_sweeps = num_sweeps,
                      burnin = burnin, Nmin = Nmin,
                      p_categorical_con = p_categorical,
                      p_categorical_mod = p_categorical, max_depth = max_depth,
                      tau_con = 2*var(y)/ntrees,
                      tau_mod = 0.5*var(y)/ntrees,parallel=T,nthread=10)
#### Plot
cate <- tau(c,w,ate,c,type,delta_tau)
test <- -Owidth+c <= x & x <= c+Owidth
cate.test <- cate[test]
# cate.test <- cate
pred <- XBART::predict.XBCFrd(fit,w[test],rep(c,sum(test)))
pred.bart.rdd <- pred$tau.adj[,(burnin+1):num_sweeps]
matplot(w[test],cbind(cate.test,rowMeans(pred.bart.rdd)),
        col=c("black","blue"),ylab=bquote(tau(X==c,W)),xlab="W",pch=19,
        cex.axis=1,cex.lab=1,cex=1,bty="L",main="BART-RDD")
##
cate <- tau(x,w,ate,c,type,delta_tau)
test <- -Owidth+c <= x & x <= c+Owidth
cate.test <- cate[test]
pred <- XBART::predict.XBCFrd(fit,w[test],x[test])
pred.bart.rdd <- pred$tau.adj[,(burnin+1):num_sweeps]
matplot(x[test],cbind(cate.test,rowMeans(pred.bart.rdd)),
        col=c("black","blue"),ylab=bquote(tau(X,W)),xlab="X",pch=19,
        cex.axis=1,cex.lab=1,cex=1,bty="L",main="BART-RDD")
# ### S-BART
# fit <- XBART::XBART(y, cbind(x,w,z), num_trees = ntrees,
#                     num_cutpoints = n, num_sweeps = num_sweeps,
#                     burnin = burnin, Nmin = Nmin,
#                     p_categorical = p_categorical+1, max_depth = max_depth,
#                     tau = var(y)/ntrees,parallel=T,nthread=10)
# test.sample <- cbind(c,w,z)
# #### Plot
# test.sample.1 <- cbind(c,w,1)[test,]
# test.sample.0 <- cbind(c,w,0)[test,]
# pred1 <- XBART::predict.XBART(fit,test.sample.1)[,(burnin+1):num_sweeps]
# pred0 <- XBART::predict.XBART(fit,test.sample.0)[,(burnin+1):num_sweeps]
# pred.sbart <- pred1-pred0
# matplot(w[test],cbind(cate.test,rowMeans(pred.sbart)),
#         col=c("black","blue"),ylab=bquote(tau(X==c,w)),xlab="W",pch=19,
#         cex.axis=1,cex.lab=1,cex=1,bty="L",main="S-BART")
# ### T-BART
# fit1 <- XBART::XBART(y[z==1], cbind(x,w)[z==1,], num_trees = ntrees,
#                      num_cutpoints = sum(z==1), num_sweeps = num_sweeps,
#                      burnin = burnin, Nmin = Nmin,
#                      p_categorical = p_categorical, max_depth = max_depth,
#                      tau = var(y[z==1])/ntrees,parallel=T,nthread=10)
# fit0 <- XBART::XBART(y[z==0], cbind(x,w)[z==0,], num_trees = ntrees,
#                      num_cutpoints = sum(z==0), num_sweeps = num_sweeps,
#                      burnin = burnin, Nmin = Nmin,
#                      p_categorical = p_categorical, max_depth = max_depth,
#                      tau = var(y[z==0])/ntrees,parallel=T,nthread=10)
# test.sample <- cbind(c,w)[test,]
# #### Plot
# pred1 <- XBART::predict.XBART(fit1,test.sample)[,(burnin+1):num_sweeps]
# pred0 <- XBART::predict.XBART(fit0,test.sample)[,(burnin+1):num_sweeps]
# pred.tbart <- pred1-pred0
# matplot(w[test],cbind(cate.test,rowMeans(pred.tbart)),
#         col=c("black","blue"),ylab=bquote(tau(X==c,W)),xlab="W",pch=19,
#         cex.axis=1,cex.lab=1,cex=1,bty="L",main="T-BART")