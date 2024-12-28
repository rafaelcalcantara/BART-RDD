set.seed(7)
### Functions
mu0.x <- function(x) 1.5*x^5 - 0.6*x^3 + 0.25*x + 0.5
mu0.w <- function(w) -15*sin(w)
tau0.x <- function(x,c) log(x-c+1)
tau0.w <- function(w) sin(pi*w)
mu <- function(x,w) {
  mu.w <- mu0.w(w)
  mu0.x(x) + sqrt(3)*mu.w/sd(mu.w)
}
tau <- function(x,c,w,ate) {
  tau.w <- tau0.w(w)
  tau0.x(x,c) + (tau.w-mean(tau.w))/sd(tau.w) + ate
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
n <- 500
ate <- 1
####
# rho <- 0.4
# u1 <- rnorm(n)
# u2 <- rnorm(n,rho*u1,sqrt(1-rho^2))
# u <- pnorm(cbind(u1,u2))
# x <- 2*qbeta(u[,1],5,1)-1
# w <- 2*qbeta(u[,2],1,5)-1
k <- 1
x <- 2*rbeta(n,2,4)-1
w <- rnorm(n,k*x,1)
c <- 0
z <- as.numeric(x>=c)
y <- mu(x,w) + tau(x,c,w,ate)*z + rnorm(n)
test <- c-h.grid(x,c,75) <= x & x <= c+h.grid(x,c,75)
par(mfrow=c(2,2))
plot(x,y,col=z+1,pch=19)
abline(v=c,lty=2)
plot(x,w)
abline(v=0)
plot(x[test],w[test])
abline(v=0)
plot(w,tau(c,c,w,ate))
# plot(w,mu0.w(w))
cor(x,w,meth="spearman")
summary(tau(c,c,w,ate))
sd(w[test])/mean(w[test])
mean(w[test])
sd(w[test])