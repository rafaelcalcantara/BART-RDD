set.seed(7)
setwd("~/../Git/BART-RDD/")
if (!dir.exists("Data")) dir.create("Data") ## Create data folder
if (length(list.files("Data"))!=0) ## Clean up folder
{
  files <- paste0("Data/",list.files("Data"))
  for (i in files) file.remove(i)
}
## Common features
s <- 1000
c <- 0
n <- c(750) ## Sample size
ate <- c(0.05,0.5)
het <- 1:3 ## Heterogeneity level
str <- c(0.5,1.5) ## Signal-to-noise ratio
nonlin <- c(0,1)
## Functions
mu <- function(x=NULL,w=NULL,het=NULL,nonlin=NULL)
{
  if (het==1) out <- 0.1875*sin((w[,1]+w[,2])*pi) + 1 + 1.875*x - 1.25*x^2 + 1.75*x^3 + nonlin*dnorm(x,0,2)*sin(2*pi*dnorm(x,0,0.05))
  if (het==2) out <- 0.1*w[,4]*sin((w[,1]+w[,2])*pi) + 1 + 0.9*w[,4] + w[,4]*x - 0.9*w[,4]*x^2 + w[,4]*x^3 + nonlin*w[,4]*dnorm(x,0,2)*sin(2*pi*dnorm(x,0,0.05))
  if (het==3)
  {
    f <- rep(0,length(x))
    for (i in 1:length(f))
    {
      if (w[i,3]==1 & w[i,4]==1) f[i] <- 0.43
      if (w[i,3]==2 & w[i,4]==1) f[i] <- 0.27
      if (w[i,3]==3 & w[i,4]==1) f[i] <- 0.1
      if (w[i,3]==1 & w[i,4]==2) f[i] <- 0.77
      if (w[i,3]==2 & w[i,4]==2) f[i] <- 0.93
      if (w[i,3]==3 & w[i,4]==2) f[i] <- 1.1
    }
    out <- 0.2*f*sin((w[,1]+w[,2])*pi) + 1 + 2*f + 2.27*f*x - 1.13*f*x^2 + 2*f*x^3 + nonlin*f*dnorm(x,0,2)*sin(2*pi*dnorm(x,0,0.05))
  }
  return(out)
}
tau <- function(x=NULL,w=NULL,het=NULL,ate=NULL)
{
  if (het==1) out <- 0.025*cos((w[,1]+w[,2])*pi) + ate - 2.8*x + 1.4*x^2 - 0.14*x^3
  if (het==2) out <- 0.0125*w[,4]*cos((w[,1]+w[,2])*pi) + ate + 0.05*(w[,4]-mean(w[,4])) - 1.8*w[,4]*x + 0.9*w[,4]*x^2 - 0.09*w[,4]*x^3
  if (het==3)
  {
    f <- rep(0,length(x))
    for (i in 1:length(f))
    {
      if (w[i,3]==1 & w[i,4]==1) f[i] <- 0.43
      if (w[i,3]==2 & w[i,4]==1) f[i] <- 0.27
      if (w[i,3]==3 & w[i,4]==1) f[i] <- 0.1
      if (w[i,3]==1 & w[i,4]==2) f[i] <- 0.77
      if (w[i,3]==2 & w[i,4]==2) f[i] <- 0.93
      if (w[i,3]==3 & w[i,4]==2) f[i] <- 1.1
    }
    out <- 0.05*f*cos((w[,1]+w[,2])*pi) + ate + 0.05*(f-mean(f)) - 3.4*f*x + 1.7*f*x^2 - 0.17*f*x^3
  }
  return(out)
}
dgp.fun <- function(n=NULL,ate=NULL,het=NULL,str=NULL,nonlin=NULL,s=NULL,c=NULL)
{
  u <- matrix(runif(n*s),n,s)
  w1 <- matrix(runif(n*s,as.numeric(u),as.numeric(u+1)),n,s)
  w2 <- matrix(runif(n*s,0,0.5),n,s)
  w3 <- matrix(rbinom(n*s,2,as.numeric(u))+1,n,s)
  w4 <- matrix(rbinom(n*s,1,0.6)+1,n,s)
  x <- matrix(2*rbeta(n*s,2,4)-as.numeric(u)-0.2+c,n,s)
  z <- apply(x,2,function(i) as.numeric(i>c))
  y <- matrix(0,n,s)
  sate <- rep(0,s)
  v <- rep(0,s)
  W <- vector("list",s)
  STR <- rep(0,s)
  for (i in 1:s)
  {
    w <- cbind(w1[,i],w2[,i],w3[,i],w4[,i])
    temp <- mu(x[,i],w,het,nonlin) + tau(x[,i],w,het,ate)*z[,i]
    sate[i] <- mean(tau(c,w,het,ate))
    v[i] <- var(tau(c,w,het,ate))
    y[,i] <- temp + rnorm(n,0,1/str)
    W[[i]] <- w
    STR[i] <- sd(temp)/sd(y[,i]-temp)
  }
  out <- list(y=y,x=x,z=z,w=W,c=c,ate=sate,str=STR,var=v)
  return(out)
}
## Generating DGPs
ind <- 0
out <- vector("list",length(n)*length(ate)*length(het)*length(str))
for (i in n)
{
  for (j in ate)
  {
    for (k in het)
    {
      for (l in str)
      {
        for (m in nonlin)
        {
          ind <- ind+1
          print(paste0("DGP ",ind))
          out[[ind]] <- dgp.fun(i,j,k,l,m,s,c)
          saveRDS(out[[ind]],paste0("Data/dgp_",ind,".rds"))
        }
      }
    }
  }
}
## Diagnostics
sapply(out, function(i) mean(i$ate))
sapply(out, function(i) mean(i$str))
sapply(out, function(i) mean(i$v))/sapply(out, function(i) mean(i$ate))
par(mfrow=c(1,2))
plot(out[[23]]$x[,1],out[[23]]$y[,1],bty="n",col=out[[23]]$z[,1]+1,pch=20,xlab="",ylab="")
plot(out[[24]]$x[,1],out[[24]]$y[,1],bty="n",col=out[[24]]$z[,1]+1,pch=20,xlab="",ylab="")