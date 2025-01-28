set.seed(007)
library(pracma)
library(MASS)
# Generate data
m <- -1 # marginal mean of x
beta <- rep(1,p)/sqrt(p)
K <- 2*toeplitz(seq(1,0,length.out = p))

beta <- rho*beta/sqrt(as.numeric(beta%*%K%*%beta))

wstar <- function(wval) rowSums(wval)/sqrt(ncol(wval)) # should be same as wval%*%beta


mu0.x <- function(x) (x+1)^3


mu0.w <- function(wval)
{
  return((wstar(wval)+2)^2)
}


mu <- function(xval,wval,k1,k3) k1*mu0.x(xval)+mu0.w(wval)*(sign(xval+1)*abs(xval+1))^k3

tau0 <- function(w) {
  pnorm(2*w[,1]+3,0,1)/2 + dnorm(w[,1],0,1)
}


tau <- function(wval,ate,k2,tau.bar) k2*(tau0(wval) - tau.bar) + ate

## Demeaning tau and adjusting parameters
Kcond <- K - K%*%beta%*%t(beta)%*%K
wcond <- mvrnorm(10000,-beta*m, Kcond)

k2.new <- k2*sd(mu(c,wcond,k1,k3))/sd(tau0(wcond))

tau.bar <- mean(tau0(wcond))

sigma_y <- k4*sd(mu(c,wcond,k1,k3))

ate <- k5*sigma_y

## Sampling
w <- mvrnorm(n,rep(0,p),K)
x <- matrix(rnorm(n*s, m + w%*%beta,1-rho^2),n,s)
z <- apply(x,2,function(xs) as.numeric(xs>=c))
y <- sapply(1:s, function(i) mu(x[,i],w,k1,k3) + tau(w,ate,k2.new,tau.bar)*z[,i] + rnorm(n,0,sigma_y))
cate <- sapply(1:s, function(i) tau(w,ate,k2.new,tau.bar))

## Saving data
### Create data folder if it doesn't exist
if ("Data" %in% list.files() == FALSE) dir.create("Data")
for (i in 1:s)
{
  ### Create folder for current DGP if it doesn't exist
  if (dgp %in% list.files("Data") == FALSE) dir.create(paste0("Data/",dgp))
  filename <- paste0("Data/",dgp,"/sample_",i,".rds")
  saveRDS(list(y=y[,i],x=x[,i],z=z[,i],w=w,cate=cate[,i],ate=ate),filename)
}


# print("sd(tau)")
# print(sd(tau(wcond,ate,k2.new,tau.bar)))
# print("sd(mu)")
# print(sd(mu(c,wcond,k1,k3)))
# print("sigma")
# print(sigma_y)