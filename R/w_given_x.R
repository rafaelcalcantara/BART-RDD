set.seed(7)
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
    while(sum.right < total | sum.left < total) 
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
c <- 0
n <- 1000
x <- 2*rbeta(n,2,4)-0.75
mtemp <- (x+0.75)/2
stemp <- 50
w <- rbeta(n,mtemp*stemp,(1-mtemp)*stemp)
test <- c-h.grid(x,c,75) <= x & x <= c+h.grid(x,c,75)
par(mfrow=c(1,2))
plot(x,w)
plot(x[test],w[test])
cor(x,w)
summary(w)
var(w)
####
rho <- 0.75
u1 <- rnorm(n)
u2 <- rnorm(n,rho*u1,sqrt(1-rho^2))
u <- pnorm(cbind(u1,u2))
# x <- 2*qbeta(u[,1],2,4)-0.75
# w <- qbeta(u[,2],2.2,5.13)
x <- qunif(u[,1],-0.2,1)
w <- qunif(u[,2],0.1,1)
test <- c-h.grid(x,c,75) <= x & x <= c+h.grid(x,c,75)
par(mfrow=c(1,2))
plot(x,w)
plot(x[test],w[test])
cor(x,w,meth="spearman")
summary(w)