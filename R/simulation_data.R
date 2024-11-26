set.seed(0)
if (!dir.exists("Data")) dir.create("Data") ## Create data folder
if (length(list.files("Data"))!=0) ## Clean up folder
{
  files <- paste0("Data/",list.files("Data"))
  for (i in files) file.remove(i)
}
### Functions
mu0.x <- function(x) 3*x^5 - 2.5*x^4 - 1.5*x^3 + 2*x^2 + 3*x + 2
mu0.w <- function(w,mu.x) w*sd(mu.x)/sd(w)
tau0.x <- function(x,c) (-exp(x)+exp(c))/(1+exp(2*x))
tau0.w <- function(w,type,ate) {
  if (type==1) out <- w < mean(w)
  # if (type==2) out <- as.integer(cut(w,quantile(w,seq(0,1,length.out=5)),include.lowest = T))
  if (type == 2) out <- as.integer(w*10)
  out <- out %% 2 == 0
  out <- (out-mean(out))/sd(out)
  return(out)
}
mu <- function(x,w,delta) {
  mu.x <- mu0.x(x)
  mu.w <- mu0.w(w,mu.x)
  (mu.x + mu.w)/sd(mu.x + mu.w)*delta
}
tau <- function(x,w,ate,c,type,delta) {
  tau.x <- tau0.x(x,c)
  tau.w <- tau0.w(w,type,ate)
  ate + (tau.x + (x*tau.w-mean(x*tau.w)) + tau.w)/sd(tau.x + (x*tau.w-mean(x*tau.w)) + tau.w)*delta
}
## Parameters
s <- 100
n <- 500
c <- 0
ate <- 0.5
delta_mu <- c(0.5,1,2)
delta_tau <- c(0.1,0.2,0.3,0.4,0.5)
ind <- 0
type <- 1:2
out <- vector("list",length(type)*length(delta_mu)*length(delta_tau))
for (tp in 1:length(type))
{
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
      # w <- matrix(as.integer(runif(n*s,1,classes[k])),n,s)
      cate <- apply(w, 2, function(i) tau(c,i,ate,c,type[tp],delta_tau[dt]))
      y <- sapply(1:s, function(i) mu(x[,i],w[,i],delta_mu[dm]) + tau(x[,i],w[,i],ate,c,type[tp],delta_tau[dt])*z[,i] + rnorm(n))
      if (tp == 1) ws <- apply(w,2,function(i) as.integer(i < mean(i)))
      # if (tp == 2) ws <- apply(w,2,function(i) as.integer(cut(i,quantile(i,seq(0,1,length.out=5)),include.lowest=T)) %% 2 == 0)
      if (tp == 2) ws <- apply(w,2,function(i) as.integer(as.integer(i*10) %% 2 == 0))
      ## Save data
      out[[ind]] <- list(y=y,x=x,z=z,w=w,ws=ws,c=c,tau.x=cate,tau=ate,type=type[tp],delta_mu=delta_mu[dm],delta_tau=delta_tau[dt])
      # out[[ind]] <- list(y=y,x=x,z=z,w=w,c=c,tau.x=tau(0,w,kappa[k],ate,c,classes[m]),tau=ate,delta_mu=delta_mu[i],delta_tau=delta_tau[j],kappa=kappa[k])
      saveRDS(out[[ind]],paste0("Data/dgp_",ind,".rds"))
    }
  }
}