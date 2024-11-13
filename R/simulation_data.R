set.seed(0)
## Needed to add this to ~/.R/Makevars file to install XBART from github,
### after installing gcc, boost and gsl
# LDFLAGS+=-L/opt/homebrew/lib
# CPPFLAGS+=-I/opt/homebrew/include
# FC = /opt/homebrew/Cellar/gcc/14.2.0_1/bin/gfortran
# F77 = /opt/homebrew/Cellar/gcc/14.2.0_1/bin/gfortran
# FLIBS = -L/opt/homebrew/Cellar/gcc/14.2.0_1/lib/gcc/14
# install_github("JingyuHe/XBART@XBCF-RDD")
library(XBART)
###
if (!dir.exists("Data")) dir.create("Data") ## Create data folder
if (length(list.files("Data"))!=0) ## Clean up folder
{
  files <- paste0("Data/",list.files("Data"))
  for (i in files) file.remove(i)
}
### Common DGP features
s <- 1000
n <- 500
c <- 0
mu0.x <- function(x) 3*x^5 - 2.5*x^4 - 1.5*x^3 + 2*x^2 + 3*x + 2
mu0.w <- function(w) w*sd(mu0.x(x))/sd(w)
tau0.x <- function(x,c) (-exp(x)+exp(c))*0.1/(1+exp(2*x))
tau0.w <- function(w,tau.bar,level) {
  w <- as.integer(w)
  if (level == 1) het <- w < mean(w)
  if (level == 2) het <- w %% 2 == 0
  if (level == 3) het <- w %% 2 == 0 | w %% 3 == 0
  if (level == 4) het <- (w %% 2 == 0 | w %% 3 == 0) & w %% 5 == 0
  if (level == 5) het <- (w %% 2 == 0) + (w %% 3 == 0) + (w %% 5 == 0)
  # het <- w
  het <- het-mean(het)
  # k <- 2*tau.bar/classes
  k <- 0.7*tau.bar/abs(min(het))
  return(k*het)
  # het <- cut(w,quantile(w,probs=seq(0,1,length.out=classes+1)),include.lowest=T)
  # het <- as.numeric(het)
  # het <- het-mean(het)
  # return(het/classes)
}
mu <- function(x,w,kappa,delta,tau) (mu0.x(x) + kappa*mu0.w(w))/sd(mu0.x(x) + kappa*mu0.w(w))*delta*sd(tau)
tau <- function(x,w,kappa,tau.bar,c,classes,level) tau.bar + tau0.x(x,c) + kappa*tau0.w(w,tau.bar,level)
#### Parameters
ate <- 1
delta_mu <- 0.5
delta_sig <- 0.5
kappa <- 1
classes <- 10
level <- 3
p <- 0.4
#### Testing the settings
x <- 2*rbeta(n,2,4)-0.75
z <- as.numeric(x>=c)
test <- -0.07 <= x & x <= 0.07
# w <- rbinom(n,classes,p) + 1
# w <- rnorm(n)
w <- as.integer(runif(n,1,classes+1))
#### X
mu.test <- mu0.x(x)
tau.test <- tau0.x(x,c) + ate
plot(x,mu.test)
plot(x,tau.test)
plot(x,mu.test + z*tau.test,col=z+1,pch=19,ylab=expression(E[X](Y[z])),xlab="X",bty="n")
abline(v=c,lty=2)
points(x[x>c],mu.test[x>c],col="lightgray",pch=21)
points(x[x<c],(mu.test+tau.test)[x<c],col="lightgray",pch=21)
#### W
mu.test <- mu0.x(x) + mu0.w(w)
tau.test <- tau0.x(x,c) + tau0.w(w,ate,level) + ate
plot(x,mu.test)
plot(x,tau.test)
plot(x,mu.test + z*tau.test,col=z+1,pch=19,ylab=expression(E[X](Y[z])),xlab="X",bty="n")
abline(v=c,lty=2)
points(x[x>c],mu.test[x>c],col="lightgray",pch=21)
points(x[x<c],(mu.test+tau.test)[x<c],col="lightgray",pch=21)
#### All together now
tau.test <- tau(x,w,kappa,ate,c,classes,level)
mu.test <- mu(x,w,kappa,delta_mu,tau.test)
plot(x,mu.test)
plot(x,tau.test)
plot(x,mu.test + z*tau.test,col=z+1,pch=19,ylab=expression(E[X](Y[z])),xlab="X",bty="n")
abline(v=c,lty=2)
points(x[x>c],mu.test[x>c],col="lightgray",pch=21)
points(x[x<c],(mu.test+tau.test)[x<c],col="lightgray",pch=21)
#### Add error
plot(x,mu.test + z*tau.test + rnorm(n,0,delta_sig*sd(tau.test)),col=z+1,pch=19,ylab="Y",xlab="X",bty="n")
abline(v=c,lty=2)
## Creating and saving the data
# ate <- seq(0.1,1,length.out=10)
# delta_mu <- c(0.5,1.25)
# delta_tau <- c(0.1,0.3)
ate <- 0.2
delta_mu <- 0.5
delta_tau <- 0.5
delta_sig <- 1
kappa <- 1
ind <- 0
classes <- c(4,8,12,16)
levels <- 1:2
out <- vector("list",length(classes)*length(delta_mu)*length(delta_tau)*length(kappa))
for (m in 1:length(levels))
{
  for (i in 1:length(delta_mu))
  {
    for (j in 1:length(delta_tau))
    {
      for (k in 1:length(classes))
      {
        ind <- ind+1
        print(paste0("DGP ",ind))
        ## Generate data
        x <- matrix(2*rbeta(n*s,2,4)-0.75,n,s)
        z <- apply(x,2,function(i) as.numeric(i>=c))
        # w <- matrix(rbinom(n*s,classes,p),n,s)
        # w <- matrix(rnorm(n*s,2,1),n,s)
        w <- matrix(runif(n*s,0,classes[k]+1),n,s)
        # w <- matrix(as.integer(runif(n*s,1,classes+1)),n,s)
        t <- tau(x,w,kappa,ate,c,classes[k],levels[m])
        y <- mu(x,w,kappa,delta_mu[i],t) + t*z
        y <- y + matrix(rnorm(n*s,0,delta_sig*sd(t)),n,s)
        ## Save data
        out[[ind]] <- list(y=y,x=x,z=z,w=w,c=c,tau.x=apply(w,2,function(i) tau(0,i,kappa,ate,c,classes[k],levels[m])),tau=ate,delta_mu=delta_mu[i],delta_tau=delta_tau[j],kappa=kappa)
        # out[[ind]] <- list(y=y,x=x,z=z,w=w,c=c,tau.x=tau(0,w,kappa[k],ate,c,classes[m]),tau=ate,delta_mu=delta_mu[i],delta_tau=delta_tau[j],kappa=kappa[k])
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