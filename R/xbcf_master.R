## Master file to run everything necessary for XBCF (fit, weights and predictions)
s <- 1000
n <- 500
p <- 10
if (Sys.info()[1]=="Linux")
{
  setwd("~/Documents/Git/XBCF-RDD")
  source("R/simulations_xbcf.R")
  source("R/simulations_kr.R")
} else
{
  setwd("~/Git/XBCF-RDD/")
  source("R/simulations_cgs.R")
}
## source("R/results_xbcf.R")
