set.seed(000)
###
library(jsonlite)
library(devtools)
## Install latest version of XBART package if available
install_github("Rafael-C-Alcantara/XBART@XBCF-RDD")
library(XBART)
###
mu <- function(W, X) return(0.1 * rowSums(W) + 1/(1+exp(-5*X)))
tau <- function(W, X) return( sin(mu(W, X)) + 1)
n <- 100
p <- 2
c <- 0 # Cutoff
w <- matrix(rnorm(n*p), n, p)
x <- rnorm(n,sd=.5)
z <- x >= c
y <- mu(w, x) + tau(w, x)*z + rnorm(n, 0, 0.1)
true.ate <- mean(tau(w,0))
## Estimation settings
### Owidth: overlap bandiwdth
### Omin: minimum number of observations inside overlap region for each leaf node
### Opct: If nb is the number of obs in node b, Opct*nb of them have to be inside overlap region
Owidth        <- seq(0.01,0.99,0.01)
Omin          <- 2
Opct          <- 0.9
num_trees_mod <- 1
num_trees_con <- 1
num_cutpoints <- n
Nmin          <- 10
num_sweeps    <- 1
burnin        <- 0
## XBCF-RDD
cutpoints <- vector("numeric",length(Owidth))
for (i in 1:length(Owidth))
{
    fit.XBCFrd <- XBCF.rd(y, w, x, c,
                          Owidth = Owidth[i], Omin = Omin,
                          Opct = Opct,
                          pcat_con = 0, pcat_mod = 0,
                          num_trees_mod = num_trees_mod,
                          num_trees_con = num_trees_con,
                          num_cutpoints = num_cutpoints,
                          num_sweeps = num_sweeps,
                          burnin = burnin, Nmin = Nmin)
    ## JSON
    trees_json <- fromJSON(fit.XBCFrd$tree_json_mod,simplifyVector=F)
### trees_json[[i]][[j]]: sweep i of tree j
    cutpoints[i] <- trees_json$trees[["0"]][["0"]][["avail.cutpoints"]]
}
## Plot available cutpoints at root node
plot(x=Owidth,y=cutpoints,type="b",col="dodgerblue4",pch=16)
