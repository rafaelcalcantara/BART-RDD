## Setup
devtools::install_github("rafaelcalcantara/XBART@XBCF-RDD")
library(XBART)
### Return available cutpoints at root node for a given Owidth
cutpoints <- function(Owidth)
{
    fit <- XBCF.rd(y, w, x, c, Owidth = Owidth, Omin = Omin, Opct = Opct,
                   num_trees_mod = 1, num_trees_con = 1,
                   num_cutpoints = num_cutpoints, num_sweeps = 1,
                   p_categorical_con = p_categorical, p_categorical_mod = p_categorical,
                   tau_con = 1/num_trees_con,
                   tau_mod = 1/num_trees_mod,
                   burnin = 0, Nmin = Nmin)
    trees_json <- jsonlite::fromJSON(fit$tree_json_mod,simplifyVector=F)
    return(trees_json$trees[["0"]][["0"]][["avail.cutpoints"]])
}
### Find minimum Owidth that leads to more than 1 available cutppoint at root node
findOwidth <- function(seq)
{
    i <- 0.01
    ## . > 1: no obs inside bandwidth (numerical error)
    ## . < 0.75: sufficient num of cutpoints available (not too restrictive)
    while(cutpoints(i)/cutpoints(sqrt(max(x^2))) > 1|
          cutpoints(i)/cutpoints(sqrt(max(x^2)))<0.75) i <- i+seq
    return(i)
}
###
c             <- 0
Omin          <- 2
Opct          <- 0.9
num_trees_mod <- 10
num_trees_con <- 10
Nmin          <- 10
num_sweeps    <- 50
burnin        <- 20
p_categorical <- 0
## DGP1a
dgp <- readRDS("Data/DGP1a.rds")
num_cutpoints <- length(dgp[[1]]$y) ## Same for all DGPs assuming n doesn't change
s <- length(dgp) ## number of samples
for (i in 1:s)
{
    data <- dgp[[i]]
    list2env(data,globalenv())
    Owidth <- findOwidth(0.01)
    t0 <- Sys.time()
    fit <- XBCF.rd(y, NULL, x, c, Owidth = Owidth, Omin = Omin, Opct = Opct,
                   num_trees_mod = num_trees_mod, num_trees_con = num_trees_con,
                   num_cutpoints = num_cutpoints, num_sweeps = num_sweeps,
                   burnin = burnin, Nmin = Nmin,
                   p_categorical_con = p_categorical, p_categorical_mod = p_categorical,
                   tau_con = 1/num_trees_con,
                   tau_mod = 1/num_trees_mod)
    t1 <- Sys.time()
    write.table(Owidth,"Results/Owidth_dgp1a.csv",append=T)
    write.table(t1-t0,"Results/time_xbcf_dgp1a.csv",append=T)
    saveRDS(fit,paste0("Results/xbcf_dgp1a_",i,".rds"))
}
## DGP1b
dgp <- readRDS("Data/DGP1b.rds")
s <- length(dgp) ## number of samples
for (i in 1:s)
{
    data <- dgp[[i]]
    list2env(data,globalenv())
    Owidth <- findOwidth(0.01)
    t0 <- Sys.time()
    fit <- XBCF.rd(y, w, x, c, Owidth = Owidth, Omin = Omin, Opct = Opct,
                   num_trees_mod = num_trees_mod, num_trees_con = num_trees_con,
                   num_cutpoints = num_cutpoints, num_sweeps = num_sweeps,
                   burnin = burnin, Nmin = Nmin,
                   p_categorical_con = p_categorical, p_categorical_mod = p_categorical,
                   tau_con = 1/num_trees_con,
                   tau_mod = 1/num_trees_mod)
    t1 <- Sys.time()
    write.table(Owidth,"Results/Owidth_dgp1b.csv",append=T)
    write.table(t1-t0,"Results/time_xbcf_dgp1b.csv",append=T)
    saveRDS(fit,paste0("Results/xbcf_dgp1b_",i,".rds"))
}
## DGP2
dgp <- readRDS("Data/DGP2.rds")
s <- length(dgp) ## number of samples
for (i in 1:s)
{
    data <- dgp[[i]]
    list2env(data,globalenv())
    Owidth <- findOwidth(0.01)
    t0 <- Sys.time()
    fit <- XBCF.rd(y, w, x, c, Owidth = Owidth, Omin = Omin, Opct = Opct,
                   num_trees_mod = num_trees_mod, num_trees_con = num_trees_con,
                   num_cutpoints = num_cutpoints, num_sweeps = num_sweeps,
                   burnin = burnin, Nmin = Nmin,
                   p_categorical_con = p_categorical, p_categorical_mod = p_categorical,
                   tau_con = 1/num_trees_con,
                   tau_mod = 1/num_trees_mod)
    t1 <- Sys.time()
    ## write.table(Owidth,"Results/Owidth_dgp2.csv",append=T)
    ## write.table(t1-t0,"Results/time_xbcf_dgp2.csv",append=T)
    saveRDS(list(fit=fit,wts=(1/x^2)/sum(1/x^2),
                 Owidth = Owidth, time=t1-t0),
            paste0("Results/xbcf_dgp2_",i,".rds"))
}
## DGP3
dgp <- readRDS("Data/DGP3.rds")
s <- length(dgp) ## number of samples
for (i in 1:s)
{
    data <- dgp[[i]]
    list2env(data,globalenv())
    Owidth <- findOwidth(0.01)
    t0 <- Sys.time()
    fit <- XBCF.rd(y, w, x, c, Owidth = Owidth, Omin = Omin, Opct = Opct,
                   num_trees_mod = num_trees_mod, num_trees_con = num_trees_con,
                   num_cutpoints = num_cutpoints, num_sweeps = num_sweeps,
                   burnin = burnin, Nmin = Nmin,
                   p_categorical_con = p_categorical, p_categorical_mod = p_categorical,
                   tau_con = 1/num_trees_con,
                   tau_mod = 1/num_trees_mod)
    t1 <- Sys.time()
    write.table(Owidth,"Results/Owidth_dgp3.csv",append=T)
    write.table(t1-t0,"Results/time_xbcf_dgp3.csv",append=T)
    saveRDS(fit,paste0("Results/xbcf_dgp3_",i,".rds"))
}
