m <- 1 # marginal mean of x
beta <- rep(1,p)/sqrt(p)
K <- 2*toeplitz(seq(1,0,length.out = p))

beta <- rho*beta/sqrt(as.numeric(beta%*%K%*%beta))

wstar <- function(wval) rowSums(wval)/sqrt(ncol(wval)) # should be same as wval%*%beta


mu0.x <- function(x) (x+1)^3


mu0.w <- function(wval)
{
  return((wstar(wval)+2)^2)
}


mu <- function(xval,wval,k1,k3,sf) sf*(k1*mu0.x(xval)+mu0.w(wval)*(sign(xval+1)*sqrt(abs(xval+1)))^k3)

tau0 <- function(w) {
  pnorm(2*w[,1]+3,0,1)/2 + dnorm(w[,1],0,1)
}


tau <- function(wval,minval,k2,tau.bar) k2*(tau0(wval) - tau.bar) + minval

## Demeaning tau and adjusting parameters
Kcond <- K - K%*%beta%*%t(beta)%*%K
wcond <- MASS::mvrnorm(10000,-beta*m, Kcond)

sf <- 1/sd(mu(c,wcond,k1,k3,1))

print(sd(mu(c,wcond,k1,k3,sf)))
k2.new <- k2/sd(tau0(wcond))

tau.bar <- mean(tau0(wcond))



sigma_y <- k4

mintau <- -k2.new*min(tau0(wcond) - tau.bar) + k5


## Sampling
# w <- MASS::mvrnorm(n,rep(0,p),K)
x <- rnorm(n, m + w%*%beta,sqrt(1-rho^2))
z <- as.numeric(x>=c)
y <- mu(x,w,k1,k3,sf) + tau(w,mintau,k2.new,tau.bar)*z + rnorm(n,0,sigma_y)
cate <- tau(w,mintau,k2.new,tau.bar)


# print("sd(tau)")
# print(sd(tau(wcond,ate,k2.new,tau.bar)))
# print("sd(mu)")
# print(sd(mu(c,wcond,k1,k3)))
# print("sigma")
# print(sigma_y)