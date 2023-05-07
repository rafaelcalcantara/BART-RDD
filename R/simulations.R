## Master script for simulations with DGPs from the other papers
## The order in the time vectors generated is DGP1, DGP3, DGP2
rm(list=ls())
n <- 500
p <- 10
c <- 0
s <- 1000 # samples
if (Sys.info()["sysname"] == "Linux")
{
    setwd("~/Documents/Git/XBCF-RDD/")
    ## XBCF-RDD simulations
    source("R/simulations_xbcf.R")
    gc()
    ## KR simulations
    source("R/simulations_kr.R")
    gc()
    ## FH simulations
    source("R/simulations_fh.R")
    gc()
} else
{
    setwd("~/Git/XBCF-RDD/")
    ## CGS simulations
    source("R/simulations_cgs.R")
    gc()
}
