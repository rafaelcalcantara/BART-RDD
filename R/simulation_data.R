set.seed(0)
setwd("~/../Git/BART-RDD/")
library(XBART)
###
if (!dir.exists("Data")) dir.create("Data") ## Create data folder
if (length(list.files("Data"))!=0) ## Clean up folder
{
  files <- paste0("Data/",list.files("Data"))
  for (i in files) file.remove(i)
}
### Common DGP features
s <- 500
n <- 1000
c <- 0
mu0.x <- function(x) 3*x^5 - 2.5*x^4 - 1.5*x^3 + 2*x^2 + 3*x + 2
mu0.w <- function(w) {
  if(is.list(w)) sapply(w, function(i) rowMeans(i[,1:2])+rowMeans(i[,3:4]))
  else rowMeans(w[,1:2])+rowMeans(w[,3:4])
}
tau0.x <- function(x) -0.1*x
tau0.w <- function(w) {
  if(is.list(w)) sapply(w, function(i) rowMeans(i[,3:4]))
  else rowMeans(w[,3:4])
}
mu <- function(x,w,kappa,delta) (mu0.x(x) + kappa*mu0.w(w))/sd(mu0.x(x) + kappa*mu0.w(w))*delta
tau <- function(x,w,kappa,delta,tau.bar) tau.bar + (tau0.x(x) + kappa*tau0.w(w))/sd(tau0.x(x) + kappa*tau0.w(w))*delta
#### Parameters
ate <- 0.75
delta_mu <- 1
delta_tau <- 0.3
kappa <- 1
#### Testing the settings
x <- 2*rbeta(n,2,4)-0.75
z <- as.numeric(x>=c)
test <- -0.07 <= x & x <= 0.07
w <- cbind(runif(n,-0.1,0.1),rnorm(n,0,0.2),rbinom(n,1,0.4)-0.4,rbinom(n,1,dnorm(x,c,0.5)))
w[,4] <- w[,4]-mean(w[,4])
#### X
mu.test <- mu0.x(x)
tau.test <- tau0.x(x) + ate
plot(x,mu.test)
plot(x,tau.test)
plot(x,mu.test + z*tau.test,col=z+1,pch=19,ylab=expression(E[X](Y[z])),xlab="X",bty="n")
abline(v=c,lty=2)
points(x[x>c],mu.test[x>c],col="lightgray",pch=21)
points(x[x<c],(mu.test+tau.test)[x<c],col="lightgray",pch=21)
#### W
mu.test <- mu0.x(x) + mu0.w(w)
tau.test <- tau0.x(x) + tau0.w(w) + ate
plot(x,mu.test)
plot(x,tau.test)
plot(x,mu.test + z*tau.test,col=z+1,pch=19,ylab=expression(E[X](Y[z])),xlab="X",bty="n")
abline(v=c,lty=2)
points(x[x>c],mu.test[x>c],col="lightgray",pch=21)
points(x[x<c],(mu.test+tau.test)[x<c],col="lightgray",pch=21)
#### All together now
mu.test <- mu(x,w,kappa,delta_mu)
tau.test <- tau(x,w,kappa,delta_tau,ate)
plot(x,mu.test)
plot(x,tau.test)
plot(x,mu.test + z*tau.test,col=z+1,pch=19,ylab=expression(E[X](Y[z])),xlab="X",bty="n")
abline(v=c,lty=2)
points(x[x>c],mu.test[x>c],col="lightgray",pch=21)
points(x[x<c],(mu.test+tau.test)[x<c],col="lightgray",pch=21)
#### Add error
plot(x,mu.test + z*tau.test + rnorm(n),col=z+1,pch=19,ylab="Y",xlab="X",bty="n")
abline(v=c,lty=2)
### Creating and saving the data
ate <- c(0.2,0.5)
delta_mu <- c(0.5,1.25)
delta_tau <- c(0.1,0.3)
kappa <- 1
ind <- 0
out <- vector("list",length(ate)*length(delta_mu)*length(delta_tau)*length(kappa))
for (m in 1:length(ate))
{
  for (i in 1:length(delta_mu))
  {
    for (j in 1:length(delta_tau))
    {
      for (k in 1:length(kappa))
      {
        ind <- ind+1
        print(paste0("DGP ",ind))
        ## Generate data
        x <- matrix(2*rbeta(n*s,2,4)-0.75,n,s)
        z <- apply(x,2,function(i) as.numeric(i>=c))
        w <- lapply(1:s,function(i) cbind(runif(n,-0.1,0.1),rnorm(n,0,0.2),rbinom(n,1,0.4)-0.4,rbinom(n,1,dnorm(x,c,0.5))))
        w <- lapply(w, function(i) cbind(i[,1:3],i[,4]-mean(i[,4])))
        y <- mu(x,w,kappa[k],delta_mu[i]) + tau(x,w,kappa[k],delta_tau[j],ate[m])*z
        y <- y + matrix(rnorm(n*s,0,1),n,s)
        ## Save data
        out[[ind]] <- list(y=y,x=x,z=z,w=w,c=c,tau.x=tau(0,w,kappa[k],delta_tau[j],ate[m]),tau=ate[m],delta_mu=delta_mu[i],delta_tau=delta_tau[j],kappa=kappa[k])
        saveRDS(out[[ind]],paste0("Data/dgp_",ind,".rds"))
      }
    }
  }
}
## lapply(out, function(i) apply(i$tau.x,2,table))
# ### Creating test samples
# n <- 1000
# s <- 1
# ind <- 0
# out <- vector("list",length(ate)*length(delta_mu)*length(delta_tau)*length(kappa))
# for (m in 1:length(ate))
# {
#   for (i in 1:length(delta_mu))
#   {
#     for (j in 1:length(delta_tau))
#     {
#       for (k in 1:length(kappa))
#       {
#         ind <- ind+1
#         print(paste0("DGP ",ind))
#         ## Generate data
#         x <- 2*rbeta(n*s,2,4)-0.75
#         w <- cbind(runif(n,-0.1,0.1),rnorm(n,0,0.2),rbinom(n,1,0.4)-0.4,rbinom(n,2,0.7)-2*0.7)
#         cate <- tau(c,w,kappa[k],delta_tau[j],ate[m])
#         ## Save data
#         out[[ind]] <- list(x=x,w=w,cate=cate,delta_mu=delta_mu[i],delta_tau=delta_tau[j],kappa=kappa[k])
#         saveRDS(out[[ind]],paste0("Data/test_dgp_",ind,".rds"))
#       }
#     }
#   }
# }