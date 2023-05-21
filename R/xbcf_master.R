## Master file to run everything necessary for XBCF (fit, weights and predictions)
setwd("~/Documents/Git/XBCF-RDD")
s <- 1000
n <- 500
p <- 10
if (Sys.info()[1]=="Linux")
    {
        source("R/simulations_xbcf.R")
        source("R/simulations_kr.R")
    } else
    {
        source("R/simulations_cgs.R")
    }
## source("R/results_xbcf.R")
