set.seed(7)
setwd("~/../Git/BART-RDD/")
if (!dir.exists("Data")) dir.create("Data") ## Create data folder
if (length(list.files("Data"))!=0) ## Clean up folder
{
  files <- paste0("Data/",list.files("Data"))
  for (i in files) file.remove(i)
}
## Common features
s <- 500
c <- 0
n <- 1000 ## Sample size
ate <- 0.05
str <- c(0.5,1,2) ## Signal-to-noise ratio
## xw: X dependent on W?
## het.mu.1/het.tau.1: dependence on W but no interaction with X
## het.mu.2/het.tau.2: dependence on W through interaction with X
xw <- 0
het <- c(0,1,2)
## Functions
mu <- function(x=NULL,w=NULL,het=NULL,ate=NULL,c)
{
  w34 <- rep(0,length(x))
  for (i in 1:length(x))
  {
    if (w[i,3]==1 & w[i,4]==1)
    {
      w34[i] <- 1
    } else if (w[i,3]==1 & w[i,4]==2)
    {
      w34[i] <- 2
    } else if (w[i,3]==2 & w[i,4]==1)
    {
      w34[i] <- 3
    } else if (w[i,3]==2 & w[i,4]==2)
    {
      w34[i] <- 4
    } else if (w[i,3]==3 & w[i,4]==1)
    {
      w34[i] <- 5
    } else if (w[i,3]==3 & w[i,4]==2)
    {
      w34[i] <- 6
    }
  }
  fw0 <- w[,1] + w[,2]
  fw <- (fw0+w34)/max(abs(fw0+w34))
  # out <- dnorm(x,c,0.1*sdx)*fw1*(1+het)*sqrt(ate)/5 + 1 + 1.875*x - 1.25*x^2 + 2*x^3
  out <- het*fw + pnorm(x,c,0.5)
  return(out)
}
tau <- function(x=NULL,w=NULL,het=NULL,ate=NULL,c=NULL)
{
  w34 <- rep(0,length(x))
  for (i in 1:length(x))
  {
    if (w[i,3]==1 & w[i,4]==1)
    {
      w34[i] <- 1
    } else if (w[i,3]==1 & w[i,4]==2)
    {
      w34[i] <- 2
    } else if (w[i,3]==2 & w[i,4]==1)
    {
      w34[i] <- 3
    } else if (w[i,3]==2 & w[i,4]==2)
    {
      w34[i] <- 4
    } else if (w[i,3]==3 & w[i,4]==1)
    {
      w34[i] <- 5
    } else if (w[i,3]==3 & w[i,4]==2)
    {
      w34[i] <- 6
    }
  }
  # fw <- 5*w[,1] + 2*w[,1]^2 + 10/(5-exp(-10*w[,2])) + 0.75*(w[,3]+w[,4])
  fw0 <- w[,1] + w[,1]^2 + log(w[,2])
  fw0 <- fw0 - mean(fw0)
  fw1 <- w34-mean(w34)
  fw <- (fw0+fw1)/max(abs(fw0+fw1))
  out <- het*fw*ate*0.5 + 2*ate/(2+x-c)
  return(out)
}
dgp.fun <- function(n=NULL,ate=NULL,het=NULL,str=NULL,xw=NULL,s=NULL,c=NULL)
{
  w1 <- matrix(runif(n*s,-2,2),n,s)
  w2 <- matrix(runif(n*s,0,4),n,s)
  w3 <- matrix(rbinom(n*s,2,0.7)+1,n,s)
  w4 <- matrix(rbinom(n*s,1,0.6)+1,n,s)
  if (xw==0) x <- matrix(2*rbeta(n*s,2,4)-0.75,n,s)
  if (xw==1)
  {
    alpha <- 2
    beta <- 21/(2*(sqrt((w1+w2)^2)))
    x <- seq(-2,1,length=n)
    p <- dunif(w1,-2,2)*dunif(w2,0,4)*dbeta(x,alpha,beta)
    ind <- apply(p,2,function(i) sample(1:n,n,T,i))
    x <- apply(ind,2,function(i) x[i]-0.25) + rnorm(n,0,0.1)
    w1 <- sapply(1:s,function(i) w1[ind[,i],i]) + rnorm(n,0,0.1)
    w2 <- sapply(1:s,function(i) w2[ind[,i],i]) + rnorm(n,0,0.1)
  }
  z <- apply(x,2,function(i) as.numeric(i>c))
  y <- matrix(0,n,s)
  W <- vector("list",s)
  prog <- matrix(0,n,s)
  cate <- matrix(0,n,s)
  for (i in 1:s)
  {
    w <- cbind(w1[,i],w2[,i],w3[,i],w4[,i])
    W[[i]] <- w
    temp <- mu(x[,i],w,het,ate,c) + tau(x[,i],w,het,ate,c)*z[,i]
    prog[,i] <- mu(rep(c,n),w,het,ate,c)
    cate[,i] <- tau(rep(c,n),w,het,ate,c)
    y[,i] <- temp + rnorm(n,0,sd(temp)/str)
  }
  out <- list(y=y,x=x,z=z,w=W,c=c,mu=prog,tau=cate,ate=colMeans(cate),het=het,str=str,xw=xw)
  return(out)
}
## Generating DGPs
ind <- 0
out <- vector("list",length(n)*length(ate)*length(het)*length(str)*length(xw))
for (i1 in n)
{
  for (i2 in ate)
  {
    for (i3 in het)
    {
      for (i4 in str)
      {
        for (i5 in xw)
        {
          ind <- ind+1
          print(paste0("DGP ",ind))
          out[[ind]] <- dgp.fun(i1,i2,i3,i4,i5,s,c)
          saveRDS(out[[ind]],paste0("Data/dgp_",ind,".rds"))
        }
      }
    }
  }
}
## Diagnostics
sapply(out,function(i) colMeans(i$tau))
sapply(out, function(i) colMeans(i$z))
cor.tab <- t(sapply(out, function(i) c(cor(i$y[,1],i$x[,1]),cor(i$y[,1],i$w[[1]]))))
cor.tab <- cbind(c("",1,"","",2,"","",3,""),round(cor.tab,2))
cor.tab <- rbind(c("k","X","W1","W2","W3","W4"),cor.tab)
stargazer::stargazer(cor.tab, summary = F, rownames = F, colnames = F)
##
pdf("Figures/sim_dgp_1.pdf")
par(mfrow=c(2,3))
for (i in c(3,6,9))
{
  plot(out[[i]]$x[,1],mu(out[[i]]$x[,1],out[[i]]$w[[1]],out[[i]]$het,ate,c),bty="n",col=out[[i]]$z[,1]+1,pch=20,xlab="X",ylab=expression(mu))
}
for (i in c(3,6,9))
{
  plot(out[[i]]$x[,1],tau(out[[i]]$x[,1],out[[i]]$w[[1]],out[[i]]$het,ate,c),bty="n",col=out[[i]]$z[,1]+1,pch=20,xlab="X",ylab=expression(tau))
}
dev.off()
###
pdf("Figures/sim_dgp_2.pdf")
par(mfrow=c(1,3))
for (i in c(3,6,9))
{
  plot(out[[i]]$x[,1],out[[i]]$y[,1],bty="n",col=out[[i]]$z[,1]+1,pch=20,xlab="",ylab="Y")
}
dev.off()
par(mfrow=c(1,1))