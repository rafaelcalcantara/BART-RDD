set.seed(0)
if (!dir.exists("Data")) dir.create("Data") ## Create data folder
if (length(list.files("Data"))!=0) ## Clean up folder
{
  files <- paste0("Data/",list.files("Data"))
  for (i in files) file.remove(i)
}
### Functions
mu0.x <- function(x) 3*x^5 - 2.5*x^4 - 1.5*x^3 + 2*x^2 + 3*x + 2
mu0.w <- function(w) w
tau0.x <- function(x,c) (exp(x)-exp(c))/(1+exp(x))
tau0.w <- function(w) {
  out <- sin(10*pi*w)
  out <- out-mean(out)
  return(out)
}
mu <- function(x,w,delta) {
  mu.x <- mu0.x(x)
  mu.w <- mu0.w(w)
  (mu.x + mu.w)/sd(mu.x + mu.w)*delta
}
tau <- function(x,c,w,delta) {
  tau.x <- tau0.x(x,c)
  tau.w <- tau0.w(w)
  ate + (tau.x + tau.w)/sd(tau.x + tau.w)*delta
}
## Parameters
s <- 200
n <- 5000
c <- 0
ate <- 0.5
delta_mu <- c(0.5,1,2)
delta_tau <- c(0.1,0.5,0.75)
ind <- 0
out <- vector("list",length(delta_mu)*length(delta_tau))
for (dm in 1:length(delta_mu))
{
  for (dt in 1:length(delta_tau))
  {
    ind <- ind+1
    print(paste0("DGP ",ind))
    ## Generate data
    x <- matrix(2*rbeta(n*s,2,4)-0.75,n,s)
    z <- apply(x,2,function(i) as.numeric(i>=c))
    w <- matrix(runif(n*s),n,s)
    cate <- apply(w, 2, function(i) tau(c,c,i,delta_tau[dt]))
    y <- sapply(1:s, function(i) mu(x[,i],w[,i],delta_mu[dm]) + tau(x[,i],c,w[,i],delta_tau[dt])*z[,i] + rnorm(n,0,sqrt(0.5)))
    ## Save data
    out[[ind]] <- list(y=y,x=x,z=z,w=w,c=c,tau.x=cate,tau=ate,delta_mu=delta_mu[dm],delta_tau=delta_tau[dt])
    saveRDS(out[[ind]],paste0("Data/dgp_",ind,".rds"))
  }
}