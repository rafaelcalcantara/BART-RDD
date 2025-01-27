## Set parameters and demean data
x0 <- rnorm(n,x.center,1)
w0 <- matrix(rnorm(n*p,rho*x.center,sqrt(1-rho^2)),n,p)
k2 <- k2*ate
k3 <- k3*sd(tau0(c,c,w0,ate,k2))/sd(mu0(c,w0,k1,1,ate))
mu.bar <- mean(mu0(c,w0,k1,k3,ate))
tau.bar <- mean(tau0(c,c,w0,ate,k2))
sig_error <- sig_error*max(abs(mean(tau(c,c,w0,ate,k2,tau.bar))),2*sd(tau(c,c,w0,ate,k2,tau.bar)))
## Samples
x <- matrix(rnorm(n*s,x.center,1),n,s)
z <- apply(x,2,function(i) as.numeric(i>=c))
w <- lapply(1:s, function(i) matrix(rnorm(n*p,rep(x[,i],p)*rho,sqrt(1-rho^2)),n,p))
prog <- sapply(1:s, function(i) mu(x[,i],w[[i]],k1,k3,ate,mu.bar))
cate <- sapply(1:s, function(i) tau(x[,i],c,w[[i]],ate,k2,tau.bar))
y <- sapply(1:s, function(i) prog[,i] + cate[,i]*z[,i] + rnorm(n,0,sig_error))
## Plots to check data
j <- 1 # sample
# plot(x[,j],y[,j],col=z[,j]+1,pch=19)
### Looking at functions at x=c
prog <- sapply(1:s, function(i) mu(c,w[[i]],k1,k3,ate,mu.bar))
cate <- sapply(1:s, function(i) tau(c,c,w[[i]],ate,k2,tau.bar))
# plot(rowMeans(w[[j]]),prog[,j])
# plot(rowMeans(w[[j]]),cate[,j])
# plot(rowMeans(w[[j]]), prog[,j] + cate[,j]*z[,j],col=z[,j]+1)
if (print.params)
{
  print(paste(c("K1: ","K2: ","K3: ","sigma: ","ATE: ","sd(mu): ","sd(tau): "),
              c(k1,k2,k3,sig_error,mean(cate[,j]),
                sd(prog[,j]),sd(cate[,j]))))
}