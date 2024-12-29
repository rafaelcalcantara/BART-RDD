set.seed(7)
### Functions
#### Steep mu.x
# mu0.x <- function(x) 0.05*x^5 - 0.6*x^3 + 1.5*x^2 + 0.25*x + 0.5
#### Not steep mu.x
mu0.x <- function(x,k) 1.5*x^5 - 0.6*x^3 + 0.25*x + 0.5
mu0.w <- function(w) sin(w)
#### Steep tau.x
# tau0.x <- function(x,c) sin(0.5*pi*x)
#### Not steep tau.x
tau0.x <- function(x,c) 0.2*x^2 + 1.2*x
tau0.w <- function(w) sin(w)
mu <- function(x,w,k) mu0.x(x,k) + mu0.w(w)
tau <- function(x,c,w,ate) tau0.x(x,c) + tau0.w(w) + ate
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
ate <- 1
n <- 500
# x <- 2*rbeta(n,2,4)-1
# mtemp <- (x+1)/2
# stemp <- 20
# w <- rbeta(n,mtemp*stemp,(1-mtemp)*stemp)-0.5
rho <- 0.9
k <- 0.2
u1 <- rnorm(n)
u2 <- rnorm(n,rho*u1,sqrt(1-rho^2))
u <- pnorm(cbind(u1,u2))
x <- qnorm(u[,1])
w <- qnorm(u[,2],0,sqrt(1))
# rho <- 0.75
# k <- 0.1
# x <- rnorm(n)
# w <- rnorm(n,rho*x,sqrt(1-rho^2))
test <- c-h.grid(x,c,75) <= x & x <= c+h.grid(x,c,75)
par(mfrow=c(1,2))
plot(x,mu(x,w,k))
plot(x,tau(x,c,w,ate))
cor(x,w)
summary(w[test])
var(w)
sd(mu(c,w,k))
sd(tau(c,c,w,ate))
sd(mu(c,w,k))/sd(tau(c,c,w,ate))
### Checking functions near x=c
plot(x[test],mu(x,w,k)[test])
plot(x[test],tau(x,c,w,ate)[test])
### Checking functions at x=c
plot(w[test],mu(c,w,k)[test])
plot(w[test],tau(c,c,w,ate)[test])
## Checking y
z <- as.numeric(x>=c)
y <- mu(x,w,k) + tau(x,c,w,ate)*z + rnorm(n,0,sqrt(1))
par(mfrow=c(1,1))
plot(x,y,col=z+1,pch=19,bty="l")
abline(v=c,lty=2)