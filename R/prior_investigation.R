set.seed(7)
setwd("~/Git/BART-RDD/")
devtools::install_github("JingyuHe/XBART@XBCF-RDD")
library(XBART)
### Setup--------------------------------------------------
par(bty="L")
Omin          <- 1
Opct          <- 0.9
ntrees_con    <- 30
ntrees_mod    <- 30
Nmin          <- 5
num_sweeps    <- 250
burnin        <- 150
p_categorical <- 0
###
mu <- function(x,w)
{
  out <- 10*w + 6.5*x^5 - 2.6*x^3 + 1.25*x + 0.5
  # 32.5*x^4 - 7.8*x^2 + 1.25
  # out <- out/sd(out)*2
  return(out)
}
tau <- function(x,w,c)
{
  out <- log(x-c+1) + 2*sin(3*pi*w)-mean(2*sin(3*pi*w))
  # out <- out/sd(out)*5
  return(out)
}
h.grid <- function(x,c,grid)
{
  abs.x <- sort(abs(x-c))
  out <- rep(0,length(grid))
  names(out) <- grid
  x.right <- sum(c < x)
  x.left <- sum(x < c)
  x.tot <- length(x)
  for(total in grid)
  {
    i <- 1
    sum.right <- sum.left <- 0
    while(sum.right < total/2 | sum.left < total/2)
    {
      sum.left <- sum(c-abs.x[i] <= x & x < c)
      sum.right <- sum(c < x & x <= c+abs.x[i])
      if (sum.left == sum(x<c) & sum.right == sum(c<x)) break
      i <- i+1
    }
    out[as.character(total)] <- abs.x[i]
  }
  return(out)
}
ate <- 2
c <- 0
n <- 10000
#### Data
x <- 2*rbeta(n,2,4)-0.75
Owidth <- h.grid(x,c,50)
z <- as.numeric(x>=c)
w <- runif(n,0,1)
Ey <- mu(x,w) + (ate+tau(x,w,c))*z
y <- Ey + rnorm(n,0,sqrt(1))
#####
par(mfrow=c(1,1))
plot(x,y,col=z+1,pch=19)
abline(v=c,lty=2)
#####
train <- c-h.grid(x,c,500) < x & x < c+h.grid(x,c,500)
x <- x[train]
y <- y[train]
w <- w[train]
### BART-RDD
fit <- XBART::XBCF.rd(y, w, x, c,
                      Owidth = Owidth, Omin = Omin, Opct = Opct,
                      num_trees_mod = ntrees_mod,
                      num_trees_con = ntrees_con,
                      num_cutpoints = n,
                      num_sweeps = num_sweeps,
                      burnin = burnin, Nmin = Nmin,
                      p_categorical_con = p_categorical,
                      p_categorical_mod = p_categorical,
                      tau_con = 2*var(y)/ntrees_con,
                      tau_mod = 0.5*var(y)/ntrees_mod,parallel=T,nthread=10)
#### Plot
cate <- tau(c,w,c)+ate
test <- -Owidth+c <= x & x <= c+Owidth
cate.test <- cate[test]
# cate.test <- cate
pred <- XBART::predict.XBCFrd(fit,w[test],rep(c,sum(test)))
pred.bart.rdd <- pred$tau.adj[,(burnin+1):num_sweeps]
matplot(w[test],cbind(cate.test,rowMeans(pred.bart.rdd)),
        col=c("black","blue"),ylab=bquote(tau(X==c,W)),xlab="W",pch=19,
        cex.axis=1,cex.lab=1,cex=1,bty="L",main="BART-RDD")