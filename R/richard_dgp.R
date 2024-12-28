set.seed(7)
### Functions
mu0.x <- function(x) 1.5*x^5 - 0.6*x^3 + 0.25*x + 0.5
mu0.w <- function(w) -15*sin(w)
tau0.x <- function(x,c) x
tau0.w <- function(w) sin(0.5*pi*w)
mu <- function(x,w) {
  mu.w <- mu0.w(w)
  mu0.x(x) + mu.w
}
tau <- function(x,c,w,ate) {
  tau0.x(x,c) + tau0.w(w) + ate
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
ate <- 1
n <- 1500
# x <- 2*rbeta(n,2,4)-1
# mtemp <- (x+1)/2
# stemp <- 20
# w <- rbeta(n,mtemp*stemp,(1-mtemp)*stemp)-0.5
# rho <- 0
# u1 <- rnorm(n)
# u2 <- rnorm(n,rho*u1,sqrt(1-rho^2))
# u <- pnorm(cbind(u1,u2))
# x <- 2*qbeta(u[,1],2,4)-1
# w <- qbeta(u[,2],2.2,5.13)-0.5
rho <- 1
x <- rnorm(n)
w <- rnorm(n,rho*x,sqrt(1-rho^2))
test <- c-h.grid(x,c,75) <= x & x <= c+h.grid(x,c,75)
par(mfrow=c(1,2))
plot(x,w)
plot(x,tau(x,c,w,ate))
cor(x,w)
summary(w)
var(w)
sd(mu(c,w))
sd(tau(c,c,w,ate))
sd(mu(c,w))/sd(tau(c,c,w,ate))