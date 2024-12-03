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
tau0.w <- function(w,level) {
  if (level==1) out <- sin(3*pi*w) ## 6th degree polynomial on W
  if (level==2) out <- sin(7*pi*w)*(w-0.5)/(5+exp(-2*w)) ## 14th degree polynomial on W
  out <- out-mean(out)
  return(out)
}
mu <- function(x,w,delta) {
  mu.x <- mu0.x(x)
  mu.w <- mu0.w(w)
  (mu.x + mu.w)/sd(mu.x + mu.w)*delta
}
tau <- function(x,c,w,delta,level) {
  tau.x <- tau0.x(x,c)
  tau.w <- tau0.w(w,level)
  ate + (tau.x + tau.w)/sd(tau.x + tau.w)*delta
}
## Parameters
s <- 200
c <- 0
ate <- 0.5
delta_mu <- c(0.5,2)
delta_tau <- c(0.25,0.75)
level <- 1:2
N <- c(500,1000)
sig_error <- c(0.25,0.5)
out <- vector("list",length(delta_mu)*length(delta_tau)*length(level)*length(N)*length(sig_error))
ind <- 0
for (dm in delta_mu)
{
  for (dt in delta_tau)
  {
    for (lvl in level)
    {
      for (n in N)
      {
        for (sig in sig_error)
        {
          ind <- ind+1
          print(paste0("DGP ",ind))
          ## Generate data
          x <- matrix(2*rbeta(n*s,2,4)-0.75,n,s)
          z <- apply(x,2,function(i) as.numeric(i>=c))
          w <- matrix(runif(n*s),n,s)
          cate <- apply(w, 2, function(i) tau(c,c,i,dt,lvl))
          y <- sapply(1:s, function(i) mu(x[,i],w[,i],dm) + tau(x[,i],c,w[,i],dt,lvl)*z[,i] + rnorm(n,0,sqrt(sig)))
          ## Save data
          out[[ind]] <- list(y=y,x=x,z=z,w=w,c=c,tau.x=cate,tau=ate,delta_mu=dm,delta_tau=dt,level=lvl,n=n,sig_error=sig)
          saveRDS(out[[ind]],paste0("Data/dgp_",ind,".rds"))
        }
      }
    }
  }
}